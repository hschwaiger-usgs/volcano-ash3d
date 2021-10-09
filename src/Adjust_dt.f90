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
         useDiffusion,useCN

      use Tephra,        only : &
         n_gs_aloft
         
      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,dx,dy,dz_vec_pd,IsLatLon,&
         kappa_pd,sigma_nx_pd,sigma_ny_pd

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
         kx,ky,kz

      implicit none

      logical, intent(in), optional :: mesostep

      integer       :: i,j,k
      real(kind=ip) :: tmp
      real(kind=ip) :: time_diffuse, time_advect
      real(kind=ip) :: dx2,dy2,dz2,tmp_sum
      real(kind=ip) :: vxmax,vxmax_dx
      real(kind=ip) :: vymax,vymax_dy
      !real(kind=ip) :: khmax,kvmax
      integer       :: fac
      real(kind=ip) :: vzmax_dz
      real(kind=ip) :: minsig
      logical       :: CheckMesoVel

      time_diffuse   = 0.0_ip

      dx2 = 2.0_ip/(dx*dx)
      dy2 = 2.0_ip/(dy*dy)
      dz2 = 2.0_ip/(minval(dz_vec_pd(1:nzmax)**2.0_ip))

      ! Get the constraining velocities relative to the grid size
      ! The bigger v[x,y,z]max_d[x,y,z] are, the smaller the time step required
      if(present(mesostep))then
        CheckMesoVel = mesostep
      else
        CheckMesoVel = .false.
      endif

      vy_meso_next_step_sp = 0.0
      vz_pd = 0.0
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
              ! Note: for this to work, you really need to set vf_pd=0 for all
              !       species that are flushed out of the system, otherwise, this
              !       will always be dominated by the large grain sizes with the
              !       highest fall velocities
              tmp =   (real(abs(vz_meso_next_step_sp(i,j,k)) + &
                     maxval(abs(vf_meso_next_step_sp(i,j,k,1:nsmax))),kind=ip)) &
                       *MPS_2_KMPHR / dz_vec_pd(k)
              if(tmp.gt.vzmax_dz) vzmax_dz = tmp
              if(IsLatLon)then
                minsig = minval(sigma_nx_pd(i:i+1,j,k))
                tmp = real(abs(vx_meso_next_step_sp(i,j,k)),kind=ip)*MPS_2_KMPHR*minsig/kappa_pd(i,j,k)
                if(tmp.gt.vxmax_dx)vxmax_dx=tmp
                minsig = minval(sigma_ny_pd(i,j:j+1,k))
                tmp = real(abs(vy_meso_next_step_sp(i,j,k)),kind=ip)*MPS_2_KMPHR*minsig/kappa_pd(i,j,k)
                if(tmp.gt.vymax_dy)vymax_dy=tmp
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
        vzmax_dz = 0.0_ip
        ! This branch looks at conditions cell-by-cell
        ! Use vx_meso_1_sp and vx_meso_2_sp
        do i=1,nxmax
          do j=1,nymax
            do k=1,nzmax
              ! Note: for this to work, you really need to set vf_pd=0 for all
              !       species that are flushed out of the system, otherwise, this
              !       will always be dominated by the large grain sizes with the
              !       highest fall velocities
              tmp =        abs(vz_pd(i,j,k)) + &
                    maxval(abs(vf_pd(i,j,k,1:nsmax)))/dz_vec_pd(k)
              if(tmp.gt.vzmax_dz)vzmax_dz = tmp
              if(IsLatLon)then
                minsig = minval(sigma_nx_pd(i:i+1,j,k))
                tmp = abs(vx_pd(i,j,k))*minsig/kappa_pd(i,j,k)
                if(tmp.gt.vxmax_dx)vxmax_dx=tmp
                minsig = minval(sigma_ny_pd(i,j:j+1,k))
                tmp = abs(vy_pd(i,j,k))*minsig/kappa_pd(i,j,k)
                if(tmp.gt.vymax_dy)vymax_dy=tmp
              endif
            enddo
          enddo
        enddo
      endif
      if(n_gs_aloft.eq.0)vzmax_dz = 0.0_ip

      ! Get the constraining diffusivities
      !khmax = max( maxval(abs(kx)), maxval(abs(ky)) )
      !kvmax = maxval(abs(kz))

      if (useDiffusion.and..not.useCN) then
          !  This is the case where diffusion limits as kx/dx*dx
          !    Note: Crank-Nicolson doesn't have this limitation
        !time_diffuse = khmax*(dx2+dy2) + kvmax*dz2
        time_diffuse = maxval(abs(kx))*dx2 + &
                       maxval(abs(ky))*dy2 + &
                       maxval(abs(kz))*dz2
      else
          ! Running case where advection limits CFL
          !  i.e. no diffusion or CN diffusion
        time_diffuse = 0.0_ip
      endif
      !time_advect = max(vxmax/dx,vymax/dy,vzmax_dz)
      time_advect = max(vxmax_dx,vymax_dy,vzmax_dz)

      tmp_sum =  time_advect + time_diffuse

      if(tmp_sum.gt.1.0_ip/DT_MIN)then
        write(global_info,*)"WARNING: Calculated DT is too low"
        write(global_info,*)"         Setting DT to dt_min = ",DT_MIN
        write(global_info,*)"         CFL condition probably violated."
        write(global_info,*)"         Check for high or NaN velocities."
        write(global_info,*)"     time_advect = ",time_advect
        write(global_info,*)"    time_diffuse = ",time_diffuse
        write(global_info,*)"    max vel.s/dx = ",vxmax_dx,vymax_dy,vzmax_dz
        write(global_info,*)"             CFL = ",CFL
        write(global_info,*)"   calculated dt = ",min(1.0_ip,CFL)/tmp_sum
        dt = DT_MIN
      elseif(tmp_sum.lt.1.0_ip/DT_MAX)then
        write(global_info,*)"WARNING: Calculated DT is too high"
        write(global_info,*)"         Setting DT to dt_max = ",DT_MAX
        write(global_info,*)"         CFL condition probably violated."
        write(global_info,*)"         Check for zero or NaN velocities."
        write(global_info,*)"     time_advect = ",time_advect
        write(global_info,*)"    time_diffuse = ",time_diffuse
        write(global_info,*)"    max vel.s/dx = ",vxmax_dx,vymax_dy,vzmax_dz
        write(global_info,*)"             CFL = ",CFL
        write(global_info,*)"   calculated dt = ",min(1.0_ip,CFL)/tmp_sum
        dt = DT_MAX
      else
        dt = min(1.0_ip,CFL)/tmp_sum
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
