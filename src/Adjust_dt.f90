!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Adjust_DT(mesostep)
!
!  Called from: MesoInterpolater
!  Arguments:
!    mesostep (optional) : logical value that indicates if a new mesostep is read.
!
!     This subroutine is used for calculating the maximal time step allowed given
!     the velocities, cell sizes, CFL number, diffusivity and the numerical scheme.
!     This subroutine is called at the end of MesoInterpolator.  The first
!     call to MesoInterpolator (and therefore to Adjust_DT) is before the time
!     loop to anticipate the total number of steps the simulation
!     might need, and then once within the time loop after new velocities are 
!     determined
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Adjust_DT(mesostep)

      use precis_param

      use io_units

      use global_param,  only : &
         EPS_SMALL,CFL,DT_MAX,DT_MIN,MPS_2_KMPHR, &
         useDiffusion,useVarDiffH,useVarDiffV,useCN

      use io_data,       only : &
         NextWriteTime

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,dx,dy,dz_vec_pd,IsLatLon,&
         kappa_pd,sigma_nx_pd,sigma_ny_pd,sigma_nz_pd

      use solution,      only : &
         vx_pd,vy_pd,vz_pd,vf_pd

      use time_data,     only : &
         time,dt,dt_ip,dt_meso_next,Simtime_in_hours

      use wind_grid,     only : &
         vx_meso_next_step_sp,vy_meso_next_step_sp,vz_meso_next_step_sp, &
         vf_meso_next_step_sp
         
      use Tephra,        only : &
         n_gs_aloft

      use Diffusion,     only : &
         kx,ky,kz,Imp_DT_fac

      implicit none

      logical, intent(in), optional :: mesostep

      integer       :: i,j,k
      real(kind=ip) :: tmp1,tmp2,tmp3,tmp4
      real(kind=dp) :: time_advect
      real(kind=dp) :: dt_tmp
      real(kind=ip) :: dx2,dy2,dz2
      real(kind=ip) :: vxmax,vxmax_dx
      real(kind=ip) :: vymax,vymax_dy
      integer       :: fac
      real(kind=ip) :: vzmax_dz
      real(kind=ip) :: minsig
      real(kind=ip) :: maxdiffus
      logical       :: CheckMesoVel

      real(kind=dp),save :: time_diffuse
      logical      ,save :: have_DT_diffus = .false.

      !-------------------------------------------------------
      !  ADVECTION
      !-------------------------------------------------------
      ! Advection time restriction follows Eq. 4.16 of LeVeque
      !   nu = abs( u dt/dx) < 1  where nu is the CFL number
      ! We impose the CFL number (default is 0.8 or so, but must be < 1.0) so we have
      !   dt = CFL / abs(u/dx)
      ! We need to find the largest u/dx quotient across the grid and
      ! similarly for v/dy and w/dz to get the most restrictive dt

      ! Get the constraining velocities relative to the grid size.
      ! The bigger v[x,y,z]max_d[x,y,z] are, the smaller the time step required
      if(present(mesostep))then
        CheckMesoVel = mesostep
      else
        CheckMesoVel = .false.
      endif
      vxmax_dx = 0.0_ip
      vymax_dy = 0.0_ip
      vzmax_dz = 0.0_ip

      ! First, determine if we are running a FAST_DT case where we are only evaluating
      ! the dt based on the bracketing Met steps, or if we are calculating dt every time
      ! the solution advances
      if(CheckMesoVel)then
        ! Met step only
        if(.not.IsLatLon) then
          ! regular projected grids have a uniform x and y geometry so it is faster
          ! to calculate these terms here
          vxmax = real(maxval(abs(vx_meso_next_step_sp(1:nxmax,1:nymax,1:nzmax))),kind=ip)
          vxmax = vxmax*MPS_2_KMPHR
          vxmax_dx = vxmax/dx
          vymax = real(maxval(abs(vy_meso_next_step_sp(1:nxmax,1:nymax,1:nzmax))),kind=ip)
          vymax = vymax*MPS_2_KMPHR
          vymax_dy = vymax/dy
        endif
        ! This branch looks at conditions cell-by-cell which is needed for Lon/Lat girds
        ! and also all z calculations (because we might have a logarithmic or variable dz)
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
              endif
              ! Advect in z
              ! Note: for this to work, we really need to set vf_meso_next_step_sp=0 for all
              !       species that are flushed out of the system.  Otherwise, this
              !       will always be dominated by the large grain sizes with the
              !       highest fall velocities
              tmp3 =   (real(abs(vz_meso_next_step_sp(i,j,k)) + &
                     maxval(abs(vf_meso_next_step_sp(i,j,k,1:nsmax))),kind=ip)) &
                       *MPS_2_KMPHR / dz_vec_pd(k)
              if(tmp3.gt.vzmax_dz) vzmax_dz = tmp3
            enddo
          enddo
        enddo
      else
        ! This is the default block which looks over all cells using the
        ! velocities for this time step
        if(.not.IsLatLon) then
          ! regular projected grids have a uniform x and y geometry so it is faster
          ! to calculate these terms here
          vxmax = maxval(abs(vx_pd(1:nxmax,1:nymax,1:nzmax)))
          vxmax_dx = vxmax/dx
          vymax = maxval(abs(vy_pd(1:nxmax,1:nymax,1:nzmax)))
          vymax_dy = vymax/dy
        endif
        ! This branch looks at conditions cell-by-cell which is needed for Lon/Lat girds
        ! and also all z calculations (because we might have a logarithmic or variable dz)
        ! Using vx_pd and vy_pd
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
              ! Advect in z
              ! Note: for this to work, we really need to set vf_pd=0 for all
              !       species that are flushed out of the system. Otherwise, this
              !       will always be dominated by the large grain sizes with the
              !       highest fall velocities.  This is done in Set_Vf_Meso.
              tmp3 =       (abs(vz_pd(i,j,k)) + &
                    maxval(abs(vf_pd(i,j,k,1:nsmax))))/dz_vec_pd(k)
              if(tmp3.gt.vzmax_dz) vzmax_dz = tmp3
            enddo
          enddo
        enddo
      endif
      ! We might have a special case where we need to continue advection even when no ash
      ! is present (between eruptive pulses, resuspension cases where source hasn't been activated)
      ! Reset vzmaz_dz if this case is special.
      if(n_gs_aloft.eq.0)vzmax_dz = 0.0_ip
      time_advect = 1.0_ip/max(vxmax_dx,vymax_dy,vzmax_dz)

      !-------------------------------------------------------
      !  DIFFUSION
      !-------------------------------------------------------
      if (useDiffusion)then
        if (useVarDiffH.or.useVarDiffV)then
          ! If we are using variable diffusivity based on current wind velocities,
          ! then we always need to calculate DT (or each met step if using FAST_DT).
          ! If neither the grid nor the diffusivities change, then we only need to
          ! do this once.
          have_DT_diffus = .false.
        endif
  
        if (.not.have_DT_diffus)then  ! If we don't have dt precalculated, enter this branch
          ! Recall that dx and dy are set in Calc_Mesh to be the most restrictive
          !   -- it is better to test cell-by-cell with the same ds measure as
          !      the diffusion routines
          dx2 = dx*dx
          dy2 = dy*dy
          dz2 = minval(dz_vec_pd(1:nzmax))**2.0_ip

          ! Note: we might want to separate these into different diffusion scenarios depending
          !       on if horizontal and vertical diffusion are independently turned on or off.
          !         1-d: only vert
          !         2-d: only horz
          !         3-d: both
          !       Current implementation assume 3-d diffusion with non-zero diffusivities
          time_diffuse = DT_MAX
          do i=1,nxmax
            do j=1,nymax
              do k=1,nzmax
                ! Diffusion in x
                minsig = minval(sigma_nx_pd(i:i+1,j,k))
                tmp1 = ((kappa_pd(i,j,k)/minsig)**2.0_ip)/kx(i,j,k)
                ! Diffusion in y
                minsig = minval(sigma_ny_pd(i,j:j+1,k))
                tmp2 = ((kappa_pd(i,j,k)/minsig)**2.0_ip)/ky(i,j,k)
                ! Diffusion in z
                minsig = minval(sigma_nz_pd(i,j,k:k+1))
                tmp3 = ((kappa_pd(i,j,k)/minsig)**2.0_ip)/kz(i,j,k)
                tmp4 = tmp1+tmp2+tmp3
                if(time_diffuse.gt.tmp4)time_diffuse=tmp4
              enddo
            enddo
          enddo

          ! Now set the time step restriction
          if (useCN) then
            ! Crank-Nicolson places no stability restriction on the time step,
            ! but accuracy considerations require dt ~ fac * dx2/k where fac is 1 to 4 or so
            ! This also applies if Backward Euler is used instead of C-N
            time_diffuse = time_diffuse * Imp_DT_fac
          else
            ! Forward Euler requires dt < 0.5 * dx2/k
            time_diffuse = time_diffuse * 0.5_dp
          endif
          maxdiffus    = min(dx2,dy2,dz2)/real(time_advect,kind=ip)

          ! Write to output streams a comment about which of advection or
          ! diffusion is the dominant constraint on the time-step.
          ! Note: we don't want to do this for the variable diffusivity cases
          !       since this will write out this comment every time-step.
          if(.not.useVarDiffH.and.&
             .not.useVarDiffV)then
            if (useCN) then
              do io=1,2;if(VB(io).le.verbosity_info)then
                write(outlog(io),'(a24,f8.2)')"Applying time factor of ",Imp_DT_fac
              endif;enddo
            else
              do io=1,2;if(VB(io).le.verbosity_info)then
                write(outlog(io),*)"Applying time factor of 0.5"
              endif;enddo
            endif

            do io=1,2;if(VB(io).le.verbosity_info)then
              if(time_diffuse.gt.time_advect)then
                ! Diffusion conditions allow a larger time step than advection
                write(outlog(io),*)"Time step limited by advection"
                if(vxmax_dx.gt.vymax_dy.and.vxmax_dx.gt.vzmax_dz)then
                  write(outlog(io),*)"   Restriction set by dx/vx"
                elseif(vymax_dy.gt.vxmax_dx.and.vymax_dy.gt.vzmax_dz)then
                  write(outlog(io),*)"   Restriction set by dy/vy"
                elseif(vzmax_dz.gt.vxmax_dx.and.vzmax_dz.gt.vymax_dy)then
                  write(outlog(io),*)"   Restriction set by dz/vz"
                endif
              else
                ! Advection is the dominant restriction on time steps
                write(outlog(io),*)"Time step limited by diffusion"
                if(tmp1.lt.tmp2.and.tmp1.lt.tmp3)then
                  write(outlog(io),*)"   Restriction set by dx2/kx"
                elseif(tmp2.lt.tmp1.and.tmp2.lt.tmp3)then
                  write(outlog(io),*)"   Restriction set by dy2/ky"
                elseif(tmp3.lt.tmp1.and.tmp3.lt.tmp2)then
                  write(outlog(io),*)"   Restriction set by dz2/kz"
                endif
              endif
              write(outlog(io),'(a42,f10.3)')"   Velocities allow up a diffusivity up to",maxdiffus
            endif;enddo
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

      dt_tmp = min(time_advect,time_diffuse)
      if(dt_tmp.lt.DT_MIN)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"WARNING: Calculated DT is too low"
          write(outlog(io),*)"         Setting DT to dt_min = ",DT_MIN
          write(outlog(io),*)"         CFL condition probably violated."
          write(outlog(io),*)"         Check for high or NaN velocities."
          write(outlog(io),*)"     time_advect = ",time_advect
          write(outlog(io),*)"    time_diffuse = ",time_diffuse
          write(outlog(io),*)"    max vel.s/dx = ",vxmax_dx,vymax_dy,vzmax_dz
          write(outlog(io),*)"             CFL = ",CFL
          write(outlog(io),*)"   calculated dt = ",min(1.0_ip,CFL)*dt_tmp
        endif;enddo
        dt = DT_MIN
      elseif(dt_tmp.gt.DT_MAX)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"WARNING: Calculated DT is too high"
          write(outlog(io),*)"         Setting DT to dt_max = ",DT_MAX
          write(outlog(io),*)"         CFL condition probably violated."
          write(outlog(io),*)"         Check for zero or NaN velocities."
          write(outlog(io),*)"     time_advect = ",time_advect
          write(outlog(io),*)"    time_diffuse = ",time_diffuse
          write(outlog(io),*)"    max vel.s/dx = ",vxmax_dx,vymax_dy,vzmax_dz
          write(outlog(io),*)"             CFL = ",CFL
          write(outlog(io),*)"   calculated dt = ",min(1.0_ip,CFL)*dt_tmp
        endif;enddo
        dt = DT_MAX
      else
        dt = real(min(1.0_ip,CFL),kind=dp) * dt_tmp
      endif
      
      ! Reset dt to be an integer multiple of DT_MIN
      fac = int(dt/DT_MIN)
      dt = DT_MIN*fac

      if (((NextWriteTime-time).gt.EPS_SMALL).and.(NextWriteTime-time.lt.dt)) then
        dt = NextWritetime-time
      endif

      if(time+dt.gt.Simtime_in_hours)then
        ! Don't let time advance past the requested stop time
        dt = Simtime_in_hours-time
      endif
      dt_ip = real(dt,kind=ip)

      if(CheckMesoVel)then
        dt_meso_next = dt
      endif

      end subroutine Adjust_DT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
