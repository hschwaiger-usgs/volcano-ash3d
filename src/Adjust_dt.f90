!##############################################################################
!
!     Adjust_DT
!
!     This subroutine is used for calculating the maximal time step allowed
!     given the velocities, cell sizes, CFL number, and the numerical 
!     scheme.
!     This subroutine is called twice from Ash3d.F90; once prior to the time
!     loop to anticipate the total number of number of steps the simulation
!     might need, and once within the time loop after new velocities are 
!     determinined.
!
!##############################################################################

      subroutine Adjust_DT(mesostep)
!      subroutine Adjust_DT

      use precis_param

      use io_units

      use global_param,  only : &
         DEG2RAD,PI,EPS_SMALL,CFL,DT_MAX,DT_MIN,MPS_2_KMPHR, &
         useDiffusion,useVarDiffH,useVarDiffV,useCN,VERB

      use Tephra,        only : &
         n_gs_aloft
         
      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,dx,dy,dz_vec_pd,IsLatLon,&
         kappa_pd,sigma_nx_pd,sigma_ny_pd,sigma_nz_pd

      use solution,      only : &
         vx_pd,vy_pd,vz_pd,vf_pd

      use wind_grid,     only : &
         vx_meso_next_step_sp,vy_meso_next_step_sp,vz_meso_next_step_sp, &
         vf_meso_next_step_sp
         
      use time_data,     only : &
         time,dt,dt_meso_next,                                    &
         dtodx,dtodxdx,dtody,dtodydy,dtodz,dtodzdz,      &
         Simtime_in_hours
         
      use io_data,       only : &
         NextWriteTime

      use Diffusion,     only : &
         kx,ky,kz,Imp_DT_fac

      implicit none

      logical, intent(in), optional :: mesostep

      integer       :: i,j,k
      real(kind=ip) :: tmp1,tmp2,tmp3
      real(kind=ip),save :: time_diffuse
      real(kind=ip) :: time_advect
      real(kind=ip) :: dx2,dy2,dz2,tmp_sum,dt_tmp
      real(kind=ip) :: vxmax,vxmax_dx
      real(kind=ip) :: vymax,vymax_dy
      !real(kind=ip) :: khmax,kvmax
      integer       :: fac
      real(kind=ip) :: vzmax_dz
      real(kind=ip) :: diffx_dx,diffy_dy,diffz_dz
      real(kind=ip) :: minsig
      logical       :: CheckMesoVel
      logical,save  :: have_DT_diffus = .false.

      !-------------------------------------------------------
      !  ADVECTION
      !-------------------------------------------------------
      ! Advection time restriction follows Eq. 4.16 of LeVeque
      !   nu = abs( u dt/dx) < 1  where nu is the CFL number
      ! We impose the CFL number (~0.8 or so, but < 1.0) so we have
      !   dt = CFL / abs(u/dx)
      ! We need to find the largest u/dx quotient across the grid and
      ! similarly for v/dy and w/dz to get the most restrictive dt

      ! Get the constraining velocities relative to the grid size
      ! The bigger v[x,y,z]max_d[x,y,z] are, the smaller the time step required
      if(present(mesostep))then
        CheckMesoVel = mesostep
      else
        CheckMesoVel = .false.
      endif

      if(CheckMesoVel)then
        ! In this block, we find the conditions based on velocities at the next
        ! meso time step
        if(.not.IsLatLon) then
          vxmax = real(maxval(abs(vx_meso_next_step_sp(1:nxmax,1:nymax,1:nzmax))),kind=ip)
          vxmax = vxmax*MPS_2_KMPHR
          vxmax_dx = vxmax/dx
          vymax = real(maxval(abs(vy_meso_next_step_sp(1:nxmax,1:nymax,1:nzmax))),kind=ip)
          vymax = vymax*MPS_2_KMPHR
          vymax_dy = vymax/dy
        else
          vxmax_dx = 0.0_ip
          vymax_dy = 0.0_ip
        endif
        vzmax_dz = 0.0_ip
        ! This branch looks at conditions cell-by-cell
        ! Use vx_meso_1_sp and vx_meso_2_sp
        do i=1,nxmax
          do j=1,nymax
            do k=1,nzmax

              if(IsLatLon)then
                ! Advect in x
                minsig = minval(sigma_nx_pd(i:i+1,j,k))
                tmp1 = real(abs(vx_meso_next_step_sp(i,j,k)),kind=ip)*MPS_2_KMPHR*minsig/kappa_pd(i,j,k)
                if(tmp1.gt.vxmax_dx)vxmax_dx=tmp1
                ! Advect in y
                minsig = minval(sigma_ny_pd(i,j:j+1,k))
                tmp2 = real(abs(vy_meso_next_step_sp(i,j,k)),kind=ip)*MPS_2_KMPHR*minsig/kappa_pd(i,j,k)
                if(tmp2.gt.vymax_dy)vymax_dy=tmp2

              ! Advect in z
              ! Note: for this to work, you really need to set vf_pd=0 for all
              !       species that are flushed out of the system, otherwise, this
              !       will always be dominated by the large grain sizes with the
              !       highest fall velocities
              tmp3 =   (real(abs(vz_meso_next_step_sp(i,j,k)) + &
                     maxval(abs(vf_meso_next_step_sp(i,j,k,1:nsmax))),kind=ip)) &
                       *MPS_2_KMPHR / dz_vec_pd(k)
              if(tmp3.gt.vzmax_dz) vzmax_dz = tmp3

              endif
            enddo
          enddo
        enddo
      else
        ! This is the normal block which looks over all cells using the
        ! velocities for this time step
        if(.not.IsLatLon) then
          vxmax = maxval(abs(vx_pd(1:nxmax,1:nymax,1:nzmax)))
          vxmax_dx = vxmax/dx
          vymax = maxval(abs(vy_pd(1:nxmax,1:nymax,1:nzmax)))
          vymax_dy = vymax/dy
        else
          vxmax_dx = 0.0_ip
          vymax_dy = 0.0_ip
        endif
        ! This branch looks at conditions cell-by-cell
        ! Use vx_meso_1_sp and vx_meso_2_sp
        ! First initialize vzmax_dz to 0.0 so we can test for the largest
        vzmax_dz = 0.0_ip
        do i=1,nxmax
          do j=1,nymax
            do k=1,nzmax

              if(IsLatLon)then
                ! Advect in x
                minsig = minval(sigma_nx_pd(i:i+1,j,k))
                tmp1 = abs(vx_pd(i,j,k))*minsig/kappa_pd(i,j,k)
                if(tmp1.gt.vxmax_dx) vxmax_dx=tmp1
                ! Advect in y
                minsig = minval(sigma_ny_pd(i,j:j+1,k))
                tmp2 = abs(vy_pd(i,j,k))*minsig/kappa_pd(i,j,k)
                if(tmp2.gt.vymax_dy) vymax_dy=tmp2
              endif
              ! Advect in y
              ! Note: for this to work, you really need to set vf_pd=0 for all
              !       species that are flushed out of the system, otherwise, this
              !       will always be dominated by the large grain sizes with the
              !       highest fall velocities
              tmp3 =       (abs(vz_pd(i,j,k)) + &
                    maxval(abs(vf_pd(i,j,k,1:nsmax))))/dz_vec_pd(k)
              if(tmp3.gt.vzmax_dz) vzmax_dz = tmp3

            enddo
          enddo
        enddo
      endif
      if(n_gs_aloft.eq.0)vzmax_dz = 0.0_ip
      !time_advect = max(vxmax/dx,vymax/dy,vzmax_dz)
      time_advect = 1.0_ip/max(vxmax_dx,vymax_dy,vzmax_dz)
      !write(*,*)"time_advect = ",vxmax_dx,vymax_dy,vzmax_dz,time_advect

      !-------------------------------------------------------
      !  DIFFUSION
      !-------------------------------------------------------

      if (useDiffusion)then
        if (useVarDiffH.or.useVarDiffV)then
          ! if we are using variable diffusivity based on current wind velocities,
          ! then we always need to calculate DT (at least on each met step
          have_DT_diffus = .false.
        endif
  
        if (.not.have_DT_diffus)then
          ! Recall that dx and dy are set in Calc_Mesh to be the most restrictive
          dx2 = dx*dx
          dy2 = dy*dy
          dz2 = minval(dz_vec_pd(1:nzmax))**2.0_ip

          diffx_dx = DT_MAX
          diffy_dy = DT_MAX
          diffz_dz = DT_MAX
          do i=1,nxmax
            do j=1,nymax
              do k=1,nzmax

                !if(IsLatLon)then
                  ! Diffusion in x
                  minsig = minval(sigma_nx_pd(i:i+1,j,k))
                  tmp1 = ((kappa_pd(i,j,k)/minsig)**2.0_ip)/kx(i,j,k)
                  if(tmp1.lt.diffx_dx) diffx_dx=tmp1
                  ! Diffusion in y
                  minsig = minval(sigma_ny_pd(i,j:j+1,k))
                  tmp2 = ((kappa_pd(i,j,k)/minsig)**2.0_ip)/ky(i,j,k)
                  if(tmp2.lt.diffy_dy) diffy_dy=tmp2
                  ! Diffusion in z
                  minsig = minval(sigma_nz_pd(i,j,k:k+1))
                  tmp3 = ((kappa_pd(i,j,k)/minsig)**2.0_ip)/kz(i,j,k)
                  if(tmp3.lt.diffz_dz) diffz_dz=tmp3
                !else
              enddo
            enddo
          enddo
          !if (useVarDiffH)then
          !  ! Here is where we should actually loop over all cells
          !  tmp1 = dx2/maxval(kx)
          !  tmp2 = dy2/maxval(ky)
          !else
          !  tmp1 = dx2/maxval(kx)
          !  tmp2 = dy2/maxval(ky)
          !endif
  
          !if(useVarDiffV)then
          !  ! Here is where we should actually loop over all cells
          !  tmp3 = dz2/maxval(kz)
          !else
          !  tmp3 = dz2/maxval(kz)
          !endif
          time_diffuse = min(tmp1,tmp2,tmp3)

          ! Now set the time step restriction
          if (useCN) then
            ! Crank-Nicolson places no stability restriction on the time step,
            ! but accuracy considerations require dt ~ fac * dx2/k where fac is 1 to 4 or so
            ! This also applies if Backward Euler is used instead of C-N
            time_diffuse = time_diffuse * Imp_DT_fac
            !write(*,*)"Appling time factor of ",Imp_DT_fac
          else
            ! Forward Euler requires dt < 0.5 * dx2/k
            !write(*,*)"Appling time factor of 0.5"
            time_diffuse = time_diffuse * 0.5_ip
          endif

          if(VERB.gt.-1)then
            if(time_diffuse.gt.time_advect)then
              ! Diffusion conditions allow a larger time step than advection
              write(global_info,*)"Time step limited by advection"
              if(vxmax_dx.gt.vymax_dy.and.vxmax_dx.gt.vzmax_dz)then
                write(global_info,*)"   Restriction set by dx/vx"
              elseif(vymax_dy.gt.vxmax_dx.and.vymax_dy.gt.vzmax_dz)then
                write(global_info,*)"   Restriction set by dy/vy"
              elseif(vzmax_dz.gt.vxmax_dx.and.vzmax_dz.gt.vymax_dy)then
                write(global_info,*)"   Restriction set by dz/vz"
              endif
            else
              ! Advection is the dominant restriction on time steps
              write(global_info,*)"Time step limited by diffusion"
              if(tmp1.lt.tmp2.and.tmp1.lt.tmp3)then
                write(global_info,*)"   Restriction set by dx2/kx"
              elseif(tmp2.lt.tmp1.and.tmp2.lt.tmp3)then
                write(global_info,*)"   Restriction set by dy2/ky"
              elseif(tmp3.lt.tmp1.and.tmp3.lt.tmp2)then
                write(global_info,*)"   Restriction set by dz2/kz"
              endif
            endif
          endif

            ! Now set logical flag so we don't need to calculate this again
          have_DT_diffus = .true.
        endif ! .not.have_DT_diffus
      else
        ! if we are not using diffusion, set the diffusion restriction to the
        ! largest time step allowed.
        time_diffuse   = DT_MAX
      endif ! useDiffusion
      !-------------------------------------------------------

      !tmp_sum =  time_advect + time_diffuse

      dt_tmp = min(time_advect,time_diffuse)
!      if(tmp_sum.gt.1.0_ip/DT_MIN)then
      if(dt_tmp.lt.DT_MIN)then
        write(global_info,*)"WARNING: Calculated DT is too low"
        write(global_info,*)"         Setting DT to dt_min = ",DT_MIN
        write(global_info,*)"         CFL condition probably violated."
        write(global_info,*)"         Check for high or NaN velocities."
        write(global_info,*)"     time_advect = ",time_advect
        write(global_info,*)"    time_diffuse = ",time_diffuse
        write(global_info,*)"    max vel.s/dx = ",vxmax_dx,vymax_dy,vzmax_dz
        write(global_info,*)"             CFL = ",CFL
        write(global_info,*)"   calculated dt = ",min(1.0_ip,CFL)*tmp_sum
        dt = DT_MIN
!      elseif(tmp_sum.lt.1.0_ip/DT_MAX)then
      elseif(dt_tmp.gt.DT_MAX)then
        write(global_info,*)"WARNING: Calculated DT is too high"
        write(global_info,*)"         Setting DT to dt_max = ",DT_MAX
        write(global_info,*)"         CFL condition probably violated."
        write(global_info,*)"         Check for zero or NaN velocities."
        write(global_info,*)"     time_advect = ",time_advect
        write(global_info,*)"    time_diffuse = ",time_diffuse
        write(global_info,*)"    max vel.s/dx = ",vxmax_dx,vymax_dy,vzmax_dz
        write(global_info,*)"             CFL = ",CFL
        write(global_info,*)"   calculated dt = ",min(1.0_ip,CFL)*tmp_sum
        dt = DT_MAX
      else
!        dt = min(1.0_ip,CFL)/tmp_sum
        dt = min(1.0_ip,CFL) * dt_tmp
      endif
      
      ! Hardwire these dt values for test cases
      !dt = 1.0e-5_ip

      ! Reset dt to be an integer multiple of DT_MIN
      fac = int(dt/DT_MIN)
      dt = DT_MIN*fac

      if(.not.IsLatLon) then
          ! Precalculate some terms used for diffusion
          !   The LatLon branch requires cell-specific values
        dtodx = dt/dx
        dtody = dt/dy
        dtodz = dt/minval(dz_vec_pd(:))

        dtodxdx = dtodx/dx
        dtodydy = dtody/dy
        dtodzdz = dtodz/minval(dz_vec_pd(:))
      endif !IsLatLon

      if (((NextWriteTime-time).gt.EPS_SMALL).and.(NextWriteTime-time.lt.dt)) then
          dt = NextWritetime-time
      endif

      if(time+dt.gt.Simtime_in_hours)then
        ! Don't let time advance past the requested stop time
        dt = Simtime_in_hours-time
      endif

      if(CheckMesoVel)then
        dt_meso_next = dt
      endif

      end subroutine Adjust_DT
