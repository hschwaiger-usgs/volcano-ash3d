      module AdvectionHorz_DCU

      use precis_param
      use global_param,  only : &
         EPS_THRESH

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,dx,dy,dz_vec_pd,ts0,ts1,&
         IsLatLon,sigma_ny_pd,sigma_nx_pd,kappa_pd,IsPeriodic

      use solution,      only : &
         concen_pd,vx_pd,vy_pd, &
         outflow_yz1_pd,outflow_yz2_pd,outflow_xz1_pd,outflow_xz2_pd,&
         SpeciesID,IsAloft

      use time_data,     only : &
         dt

      integer, parameter :: fluc_l = 1
      integer, parameter :: fluc_r = 2

      contains

!******************************************************************************

      subroutine advect_x

      ! Explicit advection routine, 2nd order upwind with limiter.
      ! RP Denlinger and HF Schwaiger

      !!!$ USE omp_lib

      implicit none

      integer       :: j,k,n  ! These are the indeces mapping to the global arrays
      integer       :: l        ! This is the index along the particular advection direction
      integer       :: ncells
      integer       :: idx_dum
      real(kind=ip) :: ubar     ! u at interface
      !real(kind=ip) :: au       ! absolute value of u at interface
      real(kind=ip) :: aus      ! absolute value of usig at interface
      real(kind=ip) :: rm2,rm1,rp1,rp2
      real(kind=ip) :: ldq,dqu,theta
      real(kind=ip) :: update
      real(kind=ip) :: divu_p, divu_m
      real(kind=ip) :: LFluct_Rbound,RFluct_Lbound
      real(kind=ip) :: LimFlux_Rbound,LimFlux_Lbound

       ! arrays that live on cell-centers: Note that we have 2 ghost cells
      real(kind=ip),dimension(-1:nxmax+2)               :: q_cc      ! concen
      real(kind=ip),dimension(-1:nxmax+2)               :: vel_cc    ! vel
      real(kind=ip),dimension(-1:nxmax+2)               :: dt_vol_cc ! dt on local cell volume
       ! arrays that live on cell interfaces
       !  Note: interface I for cell i is at (i-1/2); i.e. the left or negative side of i
       !        We only need the interfaces up to the boundary of the domain (not the ghost cells)
      real(kind=ip),dimension(-1:nxmax+2)               :: usig_I   ! vel*(interface area)
      real(kind=ip),dimension( 0:nxmax+2)               :: fss_I    ! fluctuations (A delQ)+-
      real(kind=ip),dimension( 0:nxmax+2)               :: dq_I     ! delta q
      real(kind=ip),dimension( 0:nxmax+2,fluc_l:fluc_r) :: fs_I     ! second-order term of Taylor S.(~F)

      !integer OMP_GET_MAX_THREADS
      !integer OMP_GET_NUM_THREADS
      !integer OMP_GET_THREAD_NUM
      !integer :: nthreads,thread_num
      !logical :: OMP_get_nested

      INTERFACE
        subroutine Set_BC
        end subroutine
      END INTERFACE

      call Set_BC

      ! We are advecting in x so set the length of the cell list accordingly
      ncells = nxmax

      concen_pd(:,:,:,:,ts1) = 0.0_ip

      ! Here we bundle the n and k loops into one do loop to more
      ! efficiently enable parallelization 
      !$OMP PARALLEL do &
      !$OMP  DEFAULT(NONE) &
      !$OMP  SHARED(ncells,nzmax,nymax,IsLatLon,dx,dy,dz,dt,concen_pd,kappa_pd, &
      !$OMP         vx_pd,sigma_nx_pd,outflow_yz1_pd,outflow_yz2_pd,IsPeriodic) &
      !$OMP  PRIVATE(k,n,fs_I,fss_I,usig_I, &
      !$OMP          rm2,rm1,rp1,rp2,dq,divu_p,divu_m, &
      !$OMP          dqu,theta,ldq,au,aus,update) &
      !$OMP  FIRSTPRIVATE(vol,dt_vol)
      do idx_dum=1,nsmax*nzmax
      !  do k=1,nzmax
      ! Now recover n and k
        k = mod(idx_dum-1,nzmax) +1
        n = floor(real((idx_dum-1)/nzmax,kind=ip))+1
        if(.not.IsAloft(n)) cycle

        do j=1,nymax
            ! Initialize cell-centered values for this x-row
            ! Note: ghost cells should contain q_cc=0 and vel_cc=edge
          q_cc(  -1:ncells+2) = concen_pd(-1:ncells+2,j,k,n,ts0)
          vel_cc(-1:ncells+2) =     vx_pd(-1:ncells+2,j,k)

          ! Ghost cells were set in Set_BC.f90, but could be reset here
          ! if desired or for testing.  Tests showed that velocities
          ! should either have a constant or a linear extrapolation, but
          ! concentrations should be set to zero.
          !if(IsPeriodic)then
          !  q_cc(        -1:0       ) =   q_cc(ncells-1:ncells)
          !  q_cc(  ncells+1:ncells+2) =   q_cc(       1:2)
          !  vel_cc(      -1:0       ) = vel_cc(ncells-1:ncells)
          !  vel_cc(ncells+1:ncells+2) = vel_cc(       1:2)
          !else
              ! No extrapolation, just zeros
            !q_cc(        -1:0       ) = 0.0_ip
            !q_cc(  ncells+1:ncells+2) = 0.0_ip
            !vel_cc(      -1:0       ) = 0.0_ip
            !vel_cc(ncells+1:ncells+2) = 0.0_ip
              ! zero order extrapolation (constant from edge)
            !q_cc(        -1) =   q_cc(1     )
            !q_cc(         0) =   q_cc(1     )
            !q_cc(  ncells+1) =   q_cc(ncells)
            !q_cc(  ncells+2) =   q_cc(ncells)
            !vel_cc(      -1) = vel_cc(1     )
            !vel_cc(       0) = vel_cc(1     )
            !vel_cc(ncells+1) = vel_cc(ncells)
            !vel_cc(ncells+2) = vel_cc(ncells)
              ! first order extrapolation (linear from edge)
            !q_cc(       0) = 2.0_ip*q_cc(       1)-q_cc(2)
            !q_cc(       1) = 2.0_ip*q_cc(       0)-q_cc(1)
            !q_cc(ncells+1) = 2.0_ip*q_cc(ncells  )-q_cc(ncells-1)
            !q_cc(ncells+2) = 2.0_ip*q_cc(ncells+1)-q_cc(ncells  )
            !vel_cc(       0) = 2.0_ip*vel_cc(       1)-vel_cc(2)
            !vel_cc(       1) = 2.0_ip*vel_cc(       0)-vel_cc(1)
            !vel_cc(ncells+1) = 2.0_ip*vel_cc(ncells  )-vel_cc(ncells-1)
            !vel_cc(ncells+2) = 2.0_ip*vel_cc(ncells+1)-vel_cc(ncells  )
          !endif

           ! Calculate \Delta t / \kappa
              ! using kappa of cell
          if(IsLatLon)then
            dt_vol_cc(-1:ncells+2) = dt/kappa_pd(-1:ncells+2,j,k)
          else
            dt_vol_cc(-1:ncells+2) = dt/(dx*dy*dz_vec_pd(k))
          endif
            ! Next, initialize interface values
          usig_I(-1:ncells+2)     = 0.0_ip
          dq_I(   0:ncells+2)     = 0.0_ip
          fss_I(  0:ncells+2)     = 0.0_ip
          fs_I(   0:ncells+2,1:2) = 0.0_ip

          ! First loop over all interfaces and set usig_I and fss_I
          ! Interface l is on the left (negative) side of cell l
          do l=1,ncells+1
            ubar = 0.5_ip*(vel_cc(l-1)+vel_cc(l))
            if(IsLatLon)then
              usig_I(l) = ubar*sigma_nx_pd(l,j,k)
            else
              usig_I(l) = ubar*dz_vec_pd(k)*dy
            endif

            ! Now find the limited delta-q at this interface
            rp1 = max(0.0_ip,q_cc(l  ))  ! cell on the +side of interface I
            rp2 = max(0.0_ip,q_cc(l+1))  ! second cell to right of I
            rm1 = max(0.0_ip,q_cc(l-1))  ! cell on the -side of interface I
            rm2 = max(0.0_ip,q_cc(l-2))  ! second cell to the left of I

              ! This cycling is good for production runs, but causes
              ! problems with convergence tests
#ifndef NOCYCLE
            if (rp1.le.EPS_THRESH.and.rm1.le.EPS_THRESH) cycle
#endif
            dq_I(l) = rp1-rm1

             ! Apply a limiter to the flux at interface l (Eq 6.39b)
#ifdef LIM_NONE
              ! No high-res limiter, just upwind (linear)
            ldq = 0.0_ip
#elif LIM_LAXWEN
              ! Lax-Wendrof (linear)
            ldq = dq_I(l)
#else
              ! Only calculate dqu and theta if the limiter is being used
              ! requires it.
              ! Get delta Q at upwind interface relative to interface l
            if(usig_I(l).gt.0.0_ip)then
              dqu = rm1-rm2 ! upwind is interface l-1
            else
              dqu = rp2-rp1 ! upwind is interface l+1
            endif

              ! Make sure that theta is not singular
            if(abs(dqu).gt.3.0_ip*abs(dq_I(l)).or.abs(dq_I(l)).lt.EPS_THRESH)then
              theta = sign(3.0_ip,dqu*dq_I(l))
            else
              theta = dqu/dq_I(l)
            endif
#endif

#ifdef LIM_BW
              ! Beam-Warming (linear)
            ldq = dqu
#endif
#ifdef LIM_FROMM
              ! Fromm (linear)
            ldq = dq_I(l)*(0.5_ip*(1.0_ip+theta))
#endif
#ifdef LIM_MINMOD
              ! minmod (non-linear)
            ldq = dq_I(l)*max(0.0_ip,min(1.0_ip,theta))
#endif
#ifdef LIM_SUPERBEE
              ! superbee (non-linear)
            ldq = dq_I(l)*max(0.0_ip,min(1.0_ip,2.0_ip*theta),min(2.0_ip,theta))
#endif
#ifdef LIM_MC
              ! MC (non-linear)
            ldq = dq_I(l)*max(0.0_ip,min((1.0_ip+theta)/2.0_ip,2.0_ip,2.0_ip*theta))
#endif

#ifdef LIM_NONE
            ! With no limiter, ldq=0 and this whole term evaporates.
            fss_I(l) = 0.0_ip
#else
            aus = abs(usig_I(l))
              ! Using cell-centered volume (average over interface performs more poorly)
            fss_I(l) = 0.5_ip*aus*(1.0_ip-dt_vol_cc(l)*aus)*ldq
#endif
          enddo  ! loop over l (interfaces)

          ! Next, loop over interfaces of the main grid and set the
          ! left and right fs_I.  Fluctuations between ghost cells remains initialized at 0
          ! Interface l is on the left (negative) side of cell l
          do l=0,ncells+1
              ! Set flux based on upwind velocity for color equation
              !  (equals conservative form if div.v=0)
            fs_I(l,fluc_r) = max(0.0_ip,usig_I(l))*dq_I(l) ! flux OUT OF l-1 cell to l cell
            fs_I(l,fluc_l) = min(0.0_ip,usig_I(l))*dq_I(l) ! flux OUT OF l cell to l-1 cell
              ! Modification for conservative form in
              ! divergent/convergent velocities
            !if (l.eq.1)then
            !  if(IsPeriodic)then
            !    divu_p = max(0.0_ip,usig_I(l  )) - max(0.0_ip,usig_I(ncells+1))
            !  else
            !    divu_p = max(0.0_ip,usig_I(l  ))
            !  endif
            !else
              divu_p = max(0.0_ip,usig_I(l  )) - max(0.0_ip,usig_I(l-1))
            !endif
            !if (l.eq.ncells+1)then
            !  if(IsPeriodic)then
            !    divu_m =  min(0.0_ip,usig_I(1)) - min(0.0_ip,usig_I(l  ))
            !  else
                divu_m = -min(0.0_ip,usig_I(l  ))
            !  endif
            !else
              divu_m =  min(0.0_ip,usig_I(l+1)) - min(0.0_ip,usig_I(l  ))
            !endif
            fs_I(l,fluc_r) = fs_I(l,fluc_r) + q_cc(l  ) * divu_m
            fs_I(l,fluc_l) = fs_I(l,fluc_l) + q_cc(l-1) * divu_p
          enddo  ! loop over l (interfaces)

          !--------------------------------------------------------
          ! Now loop over the cell indicies and apply these fluxes to volume l
          do l=0,ncells+1  ! Note: we loop on 1 ghost too to get outflow boundary fluxes

            !if (l.eq.0)then
            !  ! Interface fluctuation/limited-q at left boundary
            !  !  i.e. nothing is coming in from the left
            !  RFluct_Lbound  = 0.0_ip
            !  LimFlux_Lbound = 0.0_ip
            !else
            !  ! Interface fluctuation/limited-q at left cell interface
              RFluct_Lbound  =  fs_I(l  ,fluc_r)
              LimFlux_Lbound = fss_I(l)
            !endif
            !if (l.eq.ncells+1)then
            !  ! Interface fluctuation/limited-q at right boundary
            !  !  i.e. nothing is coming in from the right
            !  LFluct_Rbound  = 0.0_ip
            !  LimFlux_Rbound = 0.0_ip
            !else
              ! Interface fluctuation/limited-q at right cell interface
              LFluct_Rbound  =  fs_I(l+1,fluc_l)
              LimFlux_Rbound = fss_I(l+1)
            !endif
            ! Building Eq 6.59 of LeVeque
            !   Apply the first- and second-order term of Taylor series
            update = -dt_vol_cc(l)*( &
                        LFluct_Rbound + & ! leftward fluctuation at right boundary
                        RFluct_Lbound + & ! rightward fluctuation at left boundary
                        LimFlux_Rbound- & ! limited flux function at right boundary
                        LimFlux_Lbound)   ! limited flux function at left boundary
            concen_pd(l,j,k,n,ts1) = concen_pd(l,j,k,n,ts0) + update
            if(.not.IsPeriodic)then
              if(l.eq.0)then
                ! Flux out the - side of advection row  (W)
                outflow_yz1_pd(j,k,n) = outflow_yz1_pd(j,k,n) + update
              elseif(l.eq.ncells+1)then
                ! Flux out the + side of advection row  (E)
                outflow_yz2_pd(j,k,n) = outflow_yz2_pd(j,k,n) + update
              endif
            endif ! IsPeriodic
          enddo ! loop over l (cell centers)

        enddo ! loop over j=1,nymax
      enddo ! loop over idx_dum
      !$OMP END PARALLEL do

      concen_pd(  1:nxmax,1:nymax,1:nzmax,1:nsmax,ts0) = &
        concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts1)

      end subroutine advect_x

!******************************************************************************

      subroutine advect_y

      ! Explicit advection routine, 2nd order upwind with limiter.
      ! RP Denlinger and HF Schwaiger

      !!!$ USE omp_lib

      implicit none

      integer       :: i,k,n  ! These are the indeces mapping to the global arrays
      integer       :: l        ! This is the index along the particular advection direction
      integer       :: ncells
      integer       :: idx_dum
      real(kind=ip) :: ubar     ! u at interface
      !real(kind=ip) :: au       ! absolute value of u at interface
      real(kind=ip) :: aus      ! absolute value of usig at interface
      real(kind=ip) :: rm2,rm1,rp1,rp2
      real(kind=ip) :: ldq,dqu,theta
      real(kind=ip) :: update
      real(kind=ip) :: divu_p, divu_m
      real(kind=ip) :: LFluct_Rbound,RFluct_Lbound
      real(kind=ip) :: LimFlux_Rbound,LimFlux_Lbound

       ! arrays that live on cell-centers: Note that we have 2 ghost cells
      real(kind=ip),dimension(-1:nymax+2)               :: q_cc      ! concen
      real(kind=ip),dimension(-1:nymax+2)               :: vel_cc    ! vel
      real(kind=ip),dimension(-1:nymax+2)               :: dt_vol_cc ! dt on local cell volume
       ! arrays that live on cell interfaces
       !  Note: interface I for cell i is at (i-1/2); i.e. the left or negative side of i
       !        We only need the interfaces up to the boundary of the domain (not the ghost cells)
      real(kind=ip),dimension(-1:nymax+2)               :: usig_I   ! vel*(interface area)
      real(kind=ip),dimension( 0:nymax+2)               :: fss_I    ! fluctuations (A delQ)+-
      real(kind=ip),dimension( 0:nymax+2)               :: dq_I     ! delta q
      real(kind=ip),dimension( 0:nymax+2,fluc_l:fluc_r) :: fs_I     ! second-order term of Taylor S.(~F)
      !integer OMP_GET_MAX_THREADS
      !integer OMP_GET_NUM_THREADS
      !integer OMP_GET_THREAD_NUM
      !integer :: nthreads,thread_num
      !logical :: OMP_get_nested

      INTERFACE
        subroutine Set_BC
        end subroutine
      END INTERFACE

      call Set_BC

      ! We are advecting in y so set the length of the cell list accordingly
      ncells = nymax

      concen_pd(:,:,:,:,ts1) = 0.0_ip

      ! Here we bundle the n and k loops into one do loop to more
      ! efficiently enable parallelization 
      !$OMP PARALLEL do &
      !$OMP  DEFAULT(NONE) &
      !$OMP  SHARED(ncells,nzmax,nxmax,IsLatLon,dx,dy,dz,dt,concen_pd,kappa_pd, &
      !$OMP         vy_pd,sigma_ny_pd,outflow_xz1_pd,outflow_xz2_pd) &
      !$OMP  PRIVATE(k,n,fs_I,fss_I,usig_I, &
      !$OMP          rm2,rm1,rp1,rp2,dq_I,divu_p,divu_m, &
      !$OMP          dqu,theta,ldq,au,aus,update) &
      !$OMP  FIRSTPRIVATE(vol,dt_vol)
      do idx_dum=1,nsmax*nzmax
      !  do k=1,nzmax
      ! Now recover n and k
        k = mod(idx_dum-1,nzmax) +1
        n = floor(real((idx_dum-1)/nzmax,kind=ip))+1
        if(.not.IsAloft(n)) cycle

        do i=1,nxmax
            ! Initialize cell-centered values for this y-row
            ! Note: ghost cells should contain q_cc=0 and vel_cc=edge
          q_cc(  -1:ncells+2) = concen_pd(i,-1:ncells+2,k,n,ts0)
          vel_cc(-1:ncells+2) =     vy_pd(i,-1:ncells+2,k)

          ! Ghost cells were set in Set_BC.f90, but could be reset here
          ! if desired or for testing.  Tests showed that velocities
          ! should either have a constant or a linear extrapolation, but
          ! concentrations should be set to zero.
          !if(IsPeriodic)then
          !  q_cc(        -1:0       ) =   q_cc(ncells-1:ncells)
          !  q_cc(  ncells+1:ncells+2) =   q_cc(       1:2)
          !  vel_cc(      -1:0       ) = vel_cc(ncells-1:ncells)
          !  vel_cc(ncells+1:ncells+2) = vel_cc(       1:2)
          !else
              ! No extrapolation, just zeros
            !q_cc(        -1:0       ) = 0.0_ip
            !q_cc(  ncells+1:ncells+2) = 0.0_ip
            !vel_cc(      -1:0       ) = 0.0_ip
            !vel_cc(ncells+1:ncells+2) = 0.0_ip
              ! zero order extrapolation (constant from edge)
            !q_cc(        -1) =   q_cc(1     )
            !q_cc(         0) =   q_cc(1     )
            !q_cc(  ncells+1) =   q_cc(ncells)
            !q_cc(  ncells+2) =   q_cc(ncells)
            !vel_cc(      -1) = vel_cc(1     )
            !vel_cc(       0) = vel_cc(1     )
            !vel_cc(ncells+1) = vel_cc(ncells)
            !vel_cc(ncells+2) = vel_cc(ncells)
              ! first order extrapolation (linear from edge)
            !q_cc(       0) = 2.0_ip*q_cc(       1)-q_cc(2)
            !q_cc(       1) = 2.0_ip*q_cc(       0)-q_cc(1)
            !q_cc(ncells+1) = 2.0_ip*q_cc(ncells  )-q_cc(ncells-1)
            !q_cc(ncells+2) = 2.0_ip*q_cc(ncells+1)-q_cc(ncells  )
            !vel_cc(       0) = 2.0_ip*vel_cc(       1)-vel_cc(2)
            !vel_cc(       1) = 2.0_ip*vel_cc(       0)-vel_cc(1)
            !vel_cc(ncells+1) = 2.0_ip*vel_cc(ncells  )-vel_cc(ncells-1)
            !vel_cc(ncells+2) = 2.0_ip*vel_cc(ncells+1)-vel_cc(ncells  )

            !write(*,*)q_cc(-1:0),q_cc(ncells+1:ncells+2)
            !write(*,*)vel_cc(-1:1),vel_cc(ncells:ncells+2)

          !endif

          if(IsLatLon)then
            ! Calculate \Delta t / \kappa
              ! using kappa of cell
            dt_vol_cc(-1:ncells+2) = dt/kappa_pd(i,-1:ncells+2,k)
          else
            dt_vol_cc(-1:ncells+2) = dt/(dx*dy*dz_vec_pd(k))
          endif
            ! Next, initialize interface values
          usig_I(-1:ncells+2)     = 0.0_ip
          dq_I(   0:ncells+2)     = 0.0_ip
          fss_I(  0:ncells+2)     = 0.0_ip
          fs_I(   0:ncells+2,1:2) = 0.0_ip

          ! First loop over all interfaces and set usig_I and fss_I
          ! Interface l is on the left (negative) side of cell l
          do l=1,ncells+1
            ubar = 0.5_ip*(vel_cc(l-1)+vel_cc(l))
            if(IsLatLon)then
              usig_I(l) = ubar*sigma_ny_pd(i,l,k)
            else
              usig_I(l) = ubar*dz_vec_pd(k)*dx
            endif

            ! Now find the limited delta-q at this interface
            rp1 = max(0.0_ip,q_cc(l  ))  ! cell on the +side of interface I
            rp2 = max(0.0_ip,q_cc(l+1))  ! second cell to right of I
            rm1 = max(0.0_ip,q_cc(l-1))  ! cell on the -side of interface I
            rm2 = max(0.0_ip,q_cc(l-2))  ! second cell to the left of I

              ! This cycling is good for production runs, but causes
              ! problems with convergence tests
#ifndef NOCYCLE
            if (rp1.le.EPS_THRESH.and.rm1.le.EPS_THRESH) cycle
#endif
            dq_I(l) = rp1-rm1

             ! Apply a limiter to the flux at interface l (Eq 6.39b)
#ifdef LIM_NONE
              ! No high-res limiter, just upwind (linear)
            ldq = 0.0_ip
#elif LIM_LAXWEN
              ! Lax-Wendrof (linear)
            ldq = dq_I(l)
#else
              ! Only calculate dqu and theta if the limiter is being used
              ! requires it.
              ! Get delta Q at upwind interface relative to interface l
            if(usig_I(l).gt.0.0_ip)then
              dqu = rm1-rm2 ! upwind is interface l-1
            else
              dqu = rp2-rp1 ! upwind is interface l+1
            endif

              ! Make sure that theta is not singular
            if(abs(dqu).gt.3.0_ip*abs(dq_I(l)).or.abs(dq_I(l)).lt.EPS_THRESH)then
              theta = sign(3.0_ip,dqu*dq_I(l))
            else
              theta = dqu/dq_I(l)
            endif
#endif

#ifdef LIM_BW
              ! Beam-Warming (linear)
            ldq = dqu
#endif
#ifdef LIM_FROMM
              ! Fromm (linear)
            ldq = dq_I(l)*(0.5_ip*(1.0_ip+theta))
#endif
#ifdef LIM_MINMOD
              ! minmod (non-linear)
            ldq = dq_I(l)*max(0.0_ip,min(1.0_ip,theta))
#endif
#ifdef LIM_SUPERBEE
              ! superbee (non-linear)
            ldq = dq_I(l)*max(0.0_ip,min(1.0_ip,2.0_ip*theta),min(2.0_ip,theta))
#endif
#ifdef LIM_MC
              ! MC (non-linear)
            ldq = dq_I(l)*max(0.0_ip,min((1.0_ip+theta)/2.0_ip,2.0_ip,2.0_ip*theta))
#endif

#ifdef LIM_NONE
            ! With no limiter, ldq=0 and this whole term evaporates.
            fss_I(l) = 0.0_ip
#else
            aus = abs(usig_I(l))
              ! Using cell volume
            fss_I(l) = 0.5_ip*aus*(1.0_ip-dt_vol_cc(l)*aus)*ldq
#endif
          enddo  ! loop over l (interfaces)

          ! Next, loop over interfaces of the main grid and set the
          ! left and right fs_I.  Fluctuations between ghost cells remains initialized at 0
          ! Interface l is on the left (negative) side of cell l
          do l=0,ncells+1
              ! Set flux based on upwind velocity for color equation
              !  (equals conservative form if div.v=0)
            fs_I(l,fluc_r) = max(0.0_ip,usig_I(l))*dq_I(l) ! flux OUT OF l-1 cell to l cell
            fs_I(l,fluc_l) = min(0.0_ip,usig_I(l))*dq_I(l) ! flux OUT OF l cell to l-1 cell
              ! Modification for conservative form in
              ! divergent/convergent velocities
            divu_p   = max(0.0_ip,usig_I(l  ))  - max(0.0_ip,usig_I(l-1))
            divu_m   =  min(0.0_ip,usig_I(l+1)) - min(0.0_ip,usig_I(l))

            fs_I(l,fluc_r) = fs_I(l,fluc_r) + q_cc(l  ) * divu_m
            fs_I(l,fluc_l) = fs_I(l,fluc_l) + q_cc(l-1) * divu_p
          enddo  ! loop over l (interfaces)

          !--------------------------------------------------------
          ! Now loop over the cell indicies and apply these fluxes to volume l
          do l=0,ncells+1  ! Note: we loop on 1 ghost too to get outflow boundary fluxes

            ! Interface fluctuation/limited-q at left cell interface
            RFluct_Lbound  =  fs_I(l  ,fluc_r)
            LimFlux_Lbound = fss_I(l)
            ! Interface fluctuation/limited-q at right cell interface
            LFluct_Rbound  =  fs_I(l+1,fluc_l)
            LimFlux_Rbound = fss_I(l+1)

            ! Building Eq 6.59 of LeVeque
            !   Apply the first- and second-order term of Taylor series
            update = -dt_vol_cc(l)*( &
                        LFluct_Rbound + & ! leftward fluctuation at right boundary
                        RFluct_Lbound + & ! rightward fluctuation at left boundary
                        LimFlux_Rbound- & ! limited flux function at right boundary
                        LimFlux_Lbound)   ! limited flux function at left boundary
            concen_pd(i,l,k,n,ts1) = concen_pd(i,l,k,n,ts0) + update
            if(l.eq.0)then
              ! Flux out the - side of advection row  (N)
              outflow_xz1_pd(i,k,n) = outflow_xz1_pd(i,k,n) + update
            elseif(l.eq.ncells+1)then
              ! Flux out the + side of advection row  (S)
              outflow_xz2_pd(i,k,n) = outflow_xz2_pd(i,k,n) + update
            endif
          enddo ! loop over l (cell centers)

        enddo ! loop over i=1,nxmax
      enddo ! loop over idx_dum
      !$OMP END PARALLEL do

      concen_pd(  1:nxmax,1:nymax,1:nzmax,1:nsmax,ts0) = &
        concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts1)

      end subroutine advect_y

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module AdvectionHorz_DCU
