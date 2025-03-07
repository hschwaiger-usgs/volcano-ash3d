!##############################################################################
!
! AdvectionHorz_DCU
!
! This module contains the subroutines used for advecting material horizontally
! via the donor-cell-upwind
!
!      subroutine advect_x
!      subroutine advect_y
!
!##############################################################################

      module AdvectionHorz_DCU

      use precis_param

      use io_units

      use global_param,  only : &
         EPS_THRESH

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,dx,dy,dz_vec_pd,ts0,ts1,&
         sigma_ny_pd,sigma_nx_pd,kappa_pd,IsPeriodic,ZScaling_ID

      use solution,      only : &
         concen_pd,vx_pd,vy_pd, &
         outflow_yz1_pd,outflow_yz2_pd,outflow_xz1_pd,outflow_xz2_pd,&
         SpeciesID,IsAloft,imin,imax,jmin,jmax,kmin,kmax

      use time_data,     only : &
         dt

      use Topography,    only : &
         DelDxonD_cc,DelDyonD_cc

      implicit none

        ! Set everything to private by default
      private

        ! Publicly available subroutines/functions
      public advect_x,                  &
             advect_y

      integer, parameter :: fluc_l = 1
      integer, parameter :: fluc_r = 2

      contains

      !------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  advect_x()
!
!  Called from: AdvectHorz
!  Arguments:
!    none
!
!  This subroutine applies a 2nd order upwind advection routine with a limiter
!  determined via preprocessor flags.  The concentration array is updated in
!  concen_pd(:,:,:,:,t=2) then copied back to concen_pd(:,:,:,:,t=1).  This
!  subroutine has the equivalent structure as advect_y and advect_z
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine advect_x

      !!!$ use omp_lib

      integer       :: j,k,n  ! These are the indices mapping to the global arrays
      integer       :: i_I    ! This is the index along interfaces in the particular advection direction
      integer       :: i_cc   ! This is the index along cell-centers in the particular advection direction
      integer       :: ncells

       ! arrays that live on cell-centers: Note that we have 2 ghost cells
      real(kind=ip),dimension(-1:nxmax+2)               :: update_cc
      real(kind=ip),dimension(-1:nxmax+2)               :: q_cc      ! concen
      real(kind=ip),dimension(-1:nxmax+2)               :: vel_cc    ! vel
      real(kind=ip),dimension(-1:nxmax+2)               :: sig_I     ! cell area
      real(kind=ip),dimension(-1:nxmax+2)               :: kap_cc    ! cell volume
      real(kind=ip),dimension(-1:nxmax+2)               :: dt_vol_cc ! dt on local cell volume
      real(kind=ip),dimension(-1:nxmax+2)               :: DelDonD_cc

       ! arrays that live on cell interfaces
       !  Note: interface I for cell i is at (i-1/2); i.e. the left or negative side of i
       !        We only need the interfaces up to the boundary of the domain (not the ghost cells)
      real(kind=ip),dimension(-1:nxmax+2)               :: usig_I   ! vel*(interface area)
      ! This block is only needed for the in-line 1-d advection code as
      ! opposed to the function call
      real(kind=ip),dimension( 0:nxmax+2)     :: dq_I
      real(kind=ip),dimension( 0:nxmax+2)     :: fss_I    ! fluctuations (A delQ)+-
      real(kind=ip),dimension( 0:nxmax+2,1:2) :: fs_I     ! second-order term of Taylor S.(~F)
      real(kind=ip) :: ldq_I ! limited Delta Q
      real(kind=ip) :: dqu_I ! Delta Q at upwind interface

      real(kind=ip) :: aus      ! absolute value of usig at interface
      real(kind=ip) :: theta
      real(kind=ip) :: divu_p, divu_m
      real(kind=ip) :: LFluct_Rbound,RFluct_Lbound
      real(kind=ip) :: LimFlux_Rbound,LimFlux_Lbound
      integer :: rmin, rmax     ! min and max indicies of the row

      INTERFACE
        subroutine Set_BC(bc_code)
          integer,intent(in) :: bc_code ! 1 for advection, 2 for diffusion
        end subroutine Set_BC
      END INTERFACE

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine advect_x"
      endif;enddo

      call Set_BC(1)

      ! We are advecting in x so set the length of the cell list accordingly
      rmin = imin
      rmax = imax
      ncells = rmax - rmin + 1

      concen_pd(:,:,:,:,ts1) = 0.0_ip

      do n=1,nsmax
        if(.not.IsAloft(n)) cycle
        !$OMP PARALLEL DO                                             &
        !$OMP SCHEDULE (static)                                       &
        !$OMP DEFAULT(NONE)                                           &
        !$OMP SHARED(n,kmin,kmax,jmin,jmax,rmin,rmax,nxmax,ncells,    &
        !$OMP        dt,concen_pd,kappa_pd,IsPeriodic,                &
        !$OMP        vx_pd,sigma_nx_pd,outflow_yz1_pd,outflow_yz2_pd, &
        !$OMP        DelDxonD_cc,ZScaling_ID) &
        !$OMP PRIVATE(j,k,q_cc,vel_cc,dt_vol_cc,usig_I,update_cc,     &
        !$OMP         dq_I,fs_I,fss_I,ldq_I,dqu_I,i_I,i_cc,kap_cc,    &
        !$OMP         aus,theta,divu_p,divu_m,sig_I,                  &
        !$OMP         LFluct_Rbound,RFluct_Lbound,                    &
        !$OMP         LimFlux_Rbound,LimFlux_Lbound,DelDonD_cc)       &
        !$OMP COLLAPSE(2)
        do k=kmin,kmax
          do j=jmin,jmax
            ! Initialize cell-centered values for this x-row
            ! Note: ghost cells should contain q_cc=0 and vel_cc=edge
            q_cc(    rmin-2:rmin-1+ncells+2) = concen_pd(rmin-2:rmin-1+ncells+2,j,k,n,ts0)
            vel_cc(  rmin-2:rmin-1+ncells+2) =     vx_pd(rmin-2:rmin-1+ncells+2,j,k)
            sig_I(rmin-2:rmin-1+ncells+2)    = sigma_nx_pd(rmin-2:rmin-1+ncells+2,j,k)
            kap_cc(rmin-2:rmin-1+ncells+2)   = kappa_pd(rmin-2:rmin-1+ncells+2,j,k)

            ! Ghost cells were set in Set_BC.f90, but could be reset here
            ! if desired or for testing.  Tests showed that velocities
            ! should either have a constant or a linear extrapolation, but
            ! concentrations should be set to zero.

            ! Calculate \Delta t / \kappa
              ! using kappa of cell
            dt_vol_cc(rmin-2:rmin-1+ncells+2) = real(dt,kind=ip) / &
                                                kap_cc(rmin-2:rmin-1+ncells+2)

            ! Make sure to initialize this since we are only setting it where is matters
            usig_I      = 0.0_ip
            usig_I(rmin-1:rmin-1+ncells+1) = 0.5_ip*(vel_cc(rmin-2:rmin-1+ncells  ) + &
                                                     vel_cc(rmin-1:rmin-1+ncells+1))* &
                                                      sig_I(rmin-1:rmin-1+ncells+1)
            DelDonD_cc  = 0.0_ip
            if (ZScaling_ID.eq.2) then
              DelDonD_cc(rmin-1:rmin-1+ncells+1)=DelDxonD_cc(rmin-1:rmin-1+ncells+1,j)
            endif

            dq_I(rmin-1:rmin-1+ncells+2) = q_cc(rmin-1:rmin-1+ncells+2) - &
                                           q_cc(rmin-2:rmin-1+ncells+1)

            ! First get the limited Delta Q, in we are using high-order
            ! methods
            ldq_I = 0.0_ip
            fs_I( rmin-1:rmin-1+ncells+2,fluc_l) = 0.0_ip
            fs_I( rmin-1:rmin-1+ncells+2,fluc_r) = 0.0_ip
            fss_I(rmin-1:rmin-1+ncells+2)        = 0.0_ip

#ifndef LIM_NONE
            do i_I = rmin,rmin-1+ncells+1
              ! This cycling is good for production runs, but causes
              ! problems with convergence tests
#ifndef NOCYCLE
              if (abs(dq_I(i_I)).le.EPS_THRESH) cycle
#endif

#ifdef LIM_LAXWEN
              ! Lax-Wendroff (linear)
              ldq_I = dq_I(i_I)
#else
              ! Only calculate dqu_I and theta if the limiter is being used
              ! requires it.
              ! Get delta Q at upwind interface relative to interface I
              if(usig_I(i_I).gt.0.0_ip)then
                dqu_I = dq_I(i_I-1) ! upwind is interface I-1
              else
                dqu_I = dq_I(i_I+1) ! upwind is interface I+1
              endif

              ! Make sure that theta is not singular
              if(abs(dqu_I).gt.3.0_ip*abs(dq_I(i_I)).or. &
                 abs(dq_I(i_I)).lt.EPS_THRESH)then
                theta = sign(3.0_ip,dqu_I*dq_I(i_I))
              else
                theta = dqu_I/dq_I(i_I)
              endif
#endif
#ifdef LIM_BW
              ! Beam-Warming (linear)
              ldq_I = dqu_I
#endif
#ifdef LIM_FROMM
              ! Fromm (linear)
              ldq_I = dq_I(i_I)*(0.5_ip*(1.0_ip+theta))
#endif
#ifdef LIM_MINMOD
              ! minmod (non-linear)
              ldq_I = dq_I(i_I)*max(0.0_ip,min(1.0_ip,theta))
#endif
#ifdef LIM_SUPERBEE
              ! superbee (non-linear)
              ldq_I = dq_I(i_I)*max(0.0_ip,min(1.0_ip,2.0_ip*theta),min(2.0_ip,theta))
#endif
#ifdef LIM_MC
              ! MC (non-linear)
              ldq_I = dq_I(i_I)*max(0.0_ip,min((1.0_ip+theta)/2.0_ip,2.0_ip,2.0_ip*theta))
#endif
              aus = abs(usig_I(i_I))
                ! Using cell-centered volume (average over interface performs
                ! more poorly)
              fss_I(i_I) = 0.5_ip*aus*(1.0_ip-dt_vol_cc(i_I)*aus)*ldq_I
            enddo
#endif

            ! Next, loop over interfaces of the main grid and set the
            ! left and right fs_I.  Fluctuations between ghost cells
            ! remains initialized at 0
            ! Interface i_I is on the left (negative) side of cell i_cc
            do i_I=rmin-1,rmin-1+ncells+1
                ! Set flux based on upwind velocity for color equation
                !  (equals conservative form if div.v=0)
              fs_I(i_I,fluc_r) = max(0.0_ip,usig_I(i_I))*dq_I(i_I) ! flux OUT OF l-1 cell to l cell
              fs_I(i_I,fluc_l) = min(0.0_ip,usig_I(i_I))*dq_I(i_I) ! flux OUT OF l cell to l-1 cell
                ! Modification for conservative form in
                ! divergent/convergent velocities
              divu_p = (max(0.0_ip,usig_I(i_I  )) - max(0.0_ip,usig_I(i_I-1)))
              divu_m = (min(0.0_ip,usig_I(i_I+1)) - min(0.0_ip,usig_I(i_I  )))

              fs_I(i_I,fluc_r) = fs_I(i_I,fluc_r) + q_cc(i_I  ) * divu_m
              fs_I(i_I,fluc_l) = fs_I(i_I,fluc_l) + q_cc(i_I-1) * divu_p

              ! Topo bit
              fs_I(i_I,fluc_r) = fs_I(i_I,fluc_r) + &
                                   0.5_ip*vel_cc(i_I  )*(sig_I(i_I  )+sig_I(i_I+1)) * &
                                          q_cc(i_I  )*DelDonD_cc(i_I)
              fs_I(i_I,fluc_l) = fs_I(i_I,fluc_l) - &
                                   0.5_ip*vel_cc(i_I-1)*(sig_I(i_I-1)+sig_I(i_I  )) * &
                                          q_cc(i_I-1)*DelDonD_cc(i_I-1)

            enddo  ! loop over i_I (interfaces)

            !--------------------------------------------------------
            ! Now loop over the cell indicies and apply these fluxes to
            ! volume i_cc
            do i_cc=rmin-1,rmin-1+ncells+1  ! Note: we additionally loop on 1 ghost to get
                                            !       outflow boundary fluxes

               ! Interface fluctuation/limited-q at left cell interface
              RFluct_Lbound  =  fs_I(i_cc,fluc_r)
              LimFlux_Lbound = fss_I(i_cc)
               ! Interface fluctuation/limited-q at right cell interface
              LFluct_Rbound  =  fs_I(i_cc+1,fluc_l)
              LimFlux_Rbound = fss_I(i_cc+1)
              ! Building Eq 6.59 of LeVeque
              !   Apply the first- and second-order term of Taylor series
              update_cc(i_cc) = -dt_vol_cc(i_cc)*( &
                          LFluct_Rbound + & ! leftward fluctuation at right boundary
                          RFluct_Lbound + & ! rightward fluctuation at left boundary
                          LimFlux_Rbound- & ! limited flux function at right boundary
                          LimFlux_Lbound)   ! limited flux function at left boundary

            enddo ! loop over l (cell centers)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! Now update concentration for all interior cells
            concen_pd(rmin:rmin-1+ncells,j,k,n,ts1) = &
               concen_pd(rmin:rmin-1+ncells,j,k,n,ts0) + &
               update_cc(rmin:rmin-1+ncells)

            if(rmin.eq.1.and..not.IsPeriodic) &
              ! Flux out the - side of advection row  (W)
              outflow_yz1_pd(j,k,n) = outflow_yz1_pd(j,k,n) + update_cc(0)
            if(rmax.eq.nxmax.and..not.IsPeriodic) &
              ! Flux out the + side of advection row  (E)
              outflow_yz2_pd(j,k,n) = outflow_yz2_pd(j,k,n) + update_cc(nxmax+1)

          enddo ! loop over j=jmin,jmax
        enddo ! loop over k=kmin,kmax
      !$OMP END PARALLEL DO

      enddo ! loop over n

      concen_pd(  1:nxmax,1:nymax,1:nzmax,1:nsmax,ts0) = &
        concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts1)

      end subroutine advect_x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  advect_y()
!
!  Called from: AdvectHorz
!  Arguments:
!    none
!
!  This subroutine applies a 2nd order upwind advection routine with a limiter
!  determined via preprocessor flags.  The concentration array is updated in
!  concen_pd(:,:,:,:,t=2) then copied back to concen_pd(:,:,:,:,t=1).  This
!  subroutine has the equivalent structure as advect_x and advect_z
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine advect_y

      integer       :: i,k,n    ! These are the indices mapping to the global arrays
      integer       :: i_I    ! This is the index along interfaces in the particular advection direction
      integer       :: i_cc   ! This is the index along cell-centers in the particular advection direction
      integer       :: ncells

       ! arrays that live on cell-centers: Note that we have 2 ghost cells
      real(kind=ip),dimension(-1:nymax+2)               :: update_cc
      real(kind=ip),dimension(-1:nymax+2)               :: q_cc      ! concen
      real(kind=ip),dimension(-1:nymax+2)               :: vel_cc    ! vel
      real(kind=ip),dimension(-1:nymax+2)               :: sig_I     ! cell area
      real(kind=ip),dimension(-1:nymax+2)               :: kap_cc    ! cell volume
      real(kind=ip),dimension(-1:nymax+2)               :: dt_vol_cc ! dt on local cell volume
      real(kind=ip),dimension(-1:nymax+2)               :: DelDonD_cc
       ! arrays that live on cell interfaces
       !  Note: interface I for cell i is at (i-1/2); i.e. the left or negative side of i
       !        We only need the interfaces up to the boundary of the domain (not the ghost cells)
      real(kind=ip),dimension(-1:nymax+2)               :: usig_I   ! vel*(interface area)
      ! This block is only needed for the in-line 1-d advection code as
      ! opposed to the function call
      real(kind=ip),dimension( 0:nymax+2)     :: dq_I
      real(kind=ip),dimension( 0:nymax+2)     :: fss_I    ! fluctuations (A delQ)+-
      real(kind=ip),dimension( 0:nymax+2,1:2) :: fs_I     ! second-order term of Taylor S.(~F)
      real(kind=ip) :: ldq_I ! limited Delta Q
      real(kind=ip) :: dqu_I ! Delta Q at upwind interface

      real(kind=ip) :: aus      ! absolute value of usig at interface
      real(kind=ip) :: theta
      real(kind=ip) :: divu_p, divu_m
      real(kind=ip) :: LFluct_Rbound,RFluct_Lbound
      real(kind=ip) :: LimFlux_Rbound,LimFlux_Lbound
      integer :: rmin, rmax     ! min and max indicies of the row

      INTERFACE
        subroutine Set_BC(bc_code)
          integer,intent(in) :: bc_code ! 1 for advection, 2 for diffusion
        end subroutine Set_BC
      END INTERFACE

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine advect_y"
      endif;enddo

      call Set_BC(1)

      ! We are advecting in y so set the length of the cell list accordingly
      rmin = jmin
      rmax = jmax
      ncells = rmax - rmin + 1

      concen_pd(:,:,:,:,ts1) = 0.0_ip

      do n=1,nsmax
        if(.not.IsAloft(n)) cycle
        !$OMP PARALLEL DO                                             &
        !$OMP SCHEDULE (static)                                       &
        !$OMP DEFAULT(NONE)                                           &
        !$OMP SHARED(n,kmin,kmax,imin,imax,rmin,rmax,nymax,ncells,    &
        !$OMP        dt,concen_pd,kappa_pd,DelDyonD_cc,ZScaling_ID,   &
        !$OMP        vy_pd,sigma_ny_pd,outflow_xz1_pd,outflow_xz2_pd) &
        !$OMP PRIVATE(i,k,q_cc,vel_cc,dt_vol_cc,usig_I,update_cc,   &
        !$OMP         dq_I,fs_I,fss_I,ldq_I,dqu_I,i_I,i_cc,kap_cc,    &
        !$OMP         aus,theta,divu_p,divu_m,sig_I,                  &
        !$OMP         LFluct_Rbound,RFluct_Lbound,                    &
        !$OMP         LimFlux_Rbound,LimFlux_Lbound,DelDonD_cc)       &
        !$OMP COLLAPSE(2)

        do k=kmin,kmax
          do i=imin,imax
            ! Initialize cell-centered values for this y-row
            ! Note: ghost cells should contain q_cc=0 and vel_cc=edge
            q_cc(     rmin-2:rmin-1+ncells+2) = concen_pd(i,rmin-2:rmin-1+ncells+2,k,n,ts0)
            vel_cc(   rmin-2:rmin-1+ncells+2) =     vy_pd(i,rmin-2:rmin-1+ncells+2,k)
            sig_I(rmin-2:rmin-1+ncells+2)     = sigma_ny_pd(i,rmin-2:rmin-1+ncells+2,k)
            kap_cc(rmin-2:rmin-1+ncells+2)    = kappa_pd(i,rmin-2:rmin-1+ncells+2,k)

            ! Ghost cells were set in Set_BC.f90, but could be reset here
            ! if desired or for testing.  Tests showed that velocities
            ! should either have a constant or a linear extrapolation, but
            ! concentrations should be set to zero.

            ! Calculate \Delta t / \kappa
              ! using kappa of cell
            dt_vol_cc(rmin-2:rmin-1+ncells+2) = real(dt,kind=ip) / &
                                                kap_cc(rmin-2:rmin-1+ncells+2)

            ! Make sure to initialize this since we are only setting it where is matters
            usig_I = 0.0_ip
            usig_I(rmin-1:rmin-1+ncells+1) = 0.5_ip*(vel_cc(rmin-2:rmin-1+ncells  ) + &
                                                     vel_cc(rmin-1:rmin-1+ncells+1))* &
                                                      sig_I(rmin-1:rmin-1+ncells+1)
            DelDonD_cc  = 0.0_ip
            if (ZScaling_ID.eq.2) then
              DelDonD_cc(rmin-1:rmin-1+ncells+1)=DelDyonD_cc(i,rmin-1:rmin-1+ncells+1) 
            endif

            dq_I(rmin-1:rmin-1+ncells+2) = q_cc(rmin-1:rmin-1+ncells+2) - &
                                           q_cc(rmin-2:rmin-1+ncells+1)

            ! First get the limited Delta Q, in we are using high-order
            ! methods
            ldq_I = 0.0_ip
            fs_I( rmin-1:rmin-1+ncells+2,fluc_l) = 0.0_ip
            fs_I( rmin-1:rmin-1+ncells+2,fluc_r) = 0.0_ip
            fss_I(rmin-1:rmin-1+ncells+2) = 0.0_ip

#ifndef LIM_NONE
            do i_I = rmin,rmin-1+ncells+1
              ! This cycling is good for production runs, but causes
              ! problems with convergence tests
#ifndef NOCYCLE
              if (abs(dq_I(i_I)).le.EPS_THRESH) cycle
#endif

#ifdef LIM_LAXWEN
              ! Lax-Wendroff (linear)
              ldq_I = dq_I(i_I)
#else
              ! Only calculate dqu_I and theta if the limiter is being used
              ! requires it.
              ! Get delta Q at upwind interface relative to interface I
              if(usig_I(i_I).gt.0.0_ip)then
                dqu_I = dq_I(i_I-1) ! upwind is interface I-1
              else
                dqu_I = dq_I(i_I+1) ! upwind is interface I+1
              endif
              ! Make sure that theta is not singular
              if(abs(dqu_I).gt.3.0_ip*abs(dq_I(i_I)).or. &
                 abs(dq_I(i_I)).lt.EPS_THRESH)then
                theta = sign(3.0_ip,dqu_I*dq_I(i_I))
              else
                theta = dqu_I/dq_I(i_I)
              endif
#endif
#ifdef LIM_BW
              ! Beam-Warming (linear)
              ldq_I = dqu_I
#endif
#ifdef LIM_FROMM
              ! Fromm (linear)
              ldq_I = dq_I(i_I)*(0.5_ip*(1.0_ip+theta))
#endif
#ifdef LIM_MINMOD
              ! minmod (non-linear)
              ldq_I = dq_I(i_I)*max(0.0_ip,min(1.0_ip,theta))
#endif
#ifdef LIM_SUPERBEE
              ! superbee (non-linear)
              ldq_I = dq_I(i_I)*max(0.0_ip,min(1.0_ip,2.0_ip*theta),min(2.0_ip,theta))
#endif
#ifdef LIM_MC
              ! MC (non-linear)
              ldq_I = dq_I(i_I)*max(0.0_ip,min((1.0_ip+theta)/2.0_ip,2.0_ip,2.0_ip*theta))
#endif
              aus = abs(usig_I(i_I))
                ! Using cell-centered volume (average over interface performs
                ! more poorly)
              fss_I(i_I) = 0.5_ip*aus*(1.0_ip-dt_vol_cc(i_I)*aus)*ldq_I
            enddo
#endif

            ! Next, loop over interfaces of the main grid and set the
            ! left and right fs_I.  Fluctuations between ghost cells
            ! remains initialized at 0
            ! Interface i_I is on the left (negative) side of cell i_cc
            do i_I=rmin-1,rmin-1+ncells+1
                ! Set flux based on upwind velocity for color equation
                !  (equals conservative form if div.v=0)
              fs_I(i_I,fluc_r) = max(0.0_ip,usig_I(i_I))*dq_I(i_I) ! flux OUT OF l-1 cell to l cell
              fs_I(i_I,fluc_l) = min(0.0_ip,usig_I(i_I))*dq_I(i_I) ! flux OUT OF l cell to l-1 cell
                ! Modification for conservative form in
                ! divergent/convergent velocities
              divu_p = max(0.0_ip,usig_I(i_I  )) - max(0.0_ip,usig_I(i_I-1))
              divu_m = min(0.0_ip,usig_I(i_I+1)) - min(0.0_ip,usig_I(i_I))

              fs_I(i_I,fluc_r) = fs_I(i_I,fluc_r) + q_cc(i_I  ) * divu_m
              fs_I(i_I,fluc_l) = fs_I(i_I,fluc_l) + q_cc(i_I-1) * divu_p

              ! Topo bit
              fs_I(i_I,fluc_r) = fs_I(i_I,fluc_r) + &
                                   0.5_ip*vel_cc(i_I  )*(sig_I(i_I  )+sig_I(i_I+1)) * &
                                          q_cc(i_I  )*DelDonD_cc(i_I)
              fs_I(i_I,fluc_l) = fs_I(i_I,fluc_l) - &
                                   0.5_ip*vel_cc(i_I-1)*(sig_I(i_I-1)+sig_I(i_I  )) * &
                                          q_cc(i_I-1)*DelDonD_cc(i_I-1)

            enddo  ! loop over i_I (interfaces)

            !--------------------------------------------------------
            ! Now loop over the cell indicies and apply these fluxes to
            ! volume i_cc
            do i_cc=rmin-1,rmin-1+ncells+1  ! Note: we additionally loop on 1 ghost to get
                                            !       outflow boundary fluxes

               ! Interface fluctuation/limited-q at left cell interface
              RFluct_Lbound  =  fs_I(i_cc,fluc_r)
              LimFlux_Lbound = fss_I(i_cc)
               ! Interface fluctuation/limited-q at right cell interface
              LFluct_Rbound  =  fs_I(i_cc+1,fluc_l)
              LimFlux_Rbound = fss_I(i_cc+1)
              ! Building Eq 6.59 of LeVeque
              !   Apply the first- and second-order term of Taylor series
              update_cc(i_cc) = -dt_vol_cc(i_cc)*( &
                          LFluct_Rbound + & ! leftward fluctuation at right boundary
                          RFluct_Lbound + & ! rightward fluctuation at left boundary
                          LimFlux_Rbound- & ! limited flux function at right boundary
                          LimFlux_Lbound)   ! limited flux function at left boundary

            enddo ! loop over l (cell centers)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Now scale the update back to z

            ! Now update concentration for all interior cells
            concen_pd(i,rmin:rmin-1+ncells,k,n,ts1) = &
               concen_pd(i,rmin:rmin-1+ncells,k,n,ts0) + &
               update_cc(rmin:rmin-1+ncells)

            if(rmin.eq.1) &
              ! Flux out the - side of advection row  (N)
              outflow_xz1_pd(i,k,n) = outflow_xz1_pd(i,k,n) + update_cc(0)
            if(rmax.eq.nymax) &
              ! Flux out the + side of advection row  (S)
              outflow_xz2_pd(i,k,n) = outflow_xz2_pd(i,k,n) + update_cc(nymax+1)

          enddo ! loop over i=imin,imax
        enddo ! loop over k=kmin,kmax
      !$OMP END PARALLEL DO

      enddo ! loop over n

      concen_pd(  1:nxmax,1:nymax,1:nzmax,1:nsmax,ts0) = &
        concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts1)

      end subroutine advect_y

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module AdvectionHorz_DCU

!##############################################################################
