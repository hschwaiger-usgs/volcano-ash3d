!##############################################################################
!
! AdvectionVert_DCU
!
! This module manages the advection routine used for vertical advection.
! Only one subroutine is in this module, advect_z, which is the vertical
! advection routine using the 1-d donor-cell-upwind (DCU) method.  Alternate
! schemes could be invoked, but we have only implemented DCU.
!
!##############################################################################

      module AdvectionVert_DCU

      use precis_param

      use global_param,  only : &
         EPS_THRESH

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,dx,dy,dz_vec_pd,ts0,ts1,&
         sigma_nz_pd,kappa_pd 

      use solution,      only : &
         concen_pd,vz_pd,vf_pd, &
         outflow_xy1_pd,outflow_xy2_pd,DepositGranularity,&
         SpeciesID,IsAloft,imin,imax,jmin,jmax,kmin,kmax

      use time_data,     only : &
         dt

      implicit none

        ! Set everything to private by default
      private

        ! Publicly available subroutines/functions
      public advect_z

      integer, parameter :: fluc_l = 1
      integer, parameter :: fluc_r = 2

      contains

      !------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  advect_z()
!
!  Called from: Ash3d.F90
!  Arguments:
!    none
!
!  This subroutine applies a 2nd order upwind advection routine with a limiter
!  determined via preprocessor flags.  The concentration array is updated in
!  concen_pd(:,:,:,:,t=2) then copied back to concen_pd(:,:,:,:,t=1).  This
!  subroutine has the equivalent structure as advect_x and advect_y
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine advect_z

      !!!$ use omp_lib

      integer       :: i,j,n    ! These are the indeces mapping to the global arrays
      integer       :: l        ! This is the index along the particular advection direction
      integer       :: ncells

       ! arrays that live on cell-centers: Note that we have 2 ghost cells
      real(kind=ip),dimension(-1:nzmax+2)               :: update_cc
      real(kind=ip),dimension(-1:nzmax+2)               :: q_cc      ! concen
      real(kind=ip),dimension(-1:nzmax+2)               :: vel_cc    ! vel
      real(kind=ip),dimension(-1:nzmax+2)               :: dt_vol_cc ! dt on local cell volume
       ! arrays that live on cell interfaces
       !  Note: interface I for cell i is at (i-1/2); i.e. the left or negative side of i
       !        We only need the interfaces up to the boundary of the domain (not the ghost cells)
      real(kind=ip),dimension(-1:nzmax+2)               :: usig_I   ! vel*(interface area)

      ! This block is only needed for the in-line 1-d advection code as
      ! opposed to the function call
      real(kind=ip),dimension( 0:nzmax+2)     :: dq_I
      real(kind=ip),dimension( 0:nzmax+2)     :: fss_I    ! fluctuations (A delQ)+-
      real(kind=ip),dimension( 0:nzmax+2,1:2) :: fs_I     ! second-order term of Taylor S.(~F)
      real(kind=ip) :: ldq_I ! limited Delta Q
      real(kind=ip) :: dqu_I ! Delta Q at upwind interface
      integer :: i_I
      integer :: i_cc
      real(kind=ip) :: aus      ! absolute value of usig at interface
      real(kind=ip) :: theta
      real(kind=ip) :: divu_p, divu_m
      real(kind=ip) :: LFluct_Rbound,RFluct_Lbound
      real(kind=ip) :: LimFlux_Rbound,LimFlux_Lbound
      integer :: rmin, rmax     ! min and max indecies of the row

      !integer OMP_GET_MAX_THREADS
      !integer OMP_GET_NUM_THREADS
      !integer OMP_GET_THREAD_NUM
      !integer :: nthreads,thread_num
      !logical :: OMP_get_nested

      INTERFACE
        subroutine Set_BC(bc_code)
          integer,intent(in) :: bc_code ! 1 for advection, 2 for diffusion
        end subroutine Set_BC
        !function AdvectUpdate_1d(ncells,q_cc,dt_vol_cc,usig_I)
        !  integer :: ncells
        !  real(kind=8),dimension(-1:ncells+2) :: AdvectUpdate_1d
        !  real(kind=8),dimension(-1:ncells+2) :: q_cc
        !  real(kind=8),dimension(-1:ncells+2) :: dt_vol_cc
        !  real(kind=8),dimension(-1:ncells+2) :: usig_I
        !end function AdvectUpdate_1d
      END INTERFACE

      call Set_BC(1)

      ! We are advecting in z so set the length of the cell list accordingly
      rmin = kmin
      rmax = kmax
      ncells = rmax - rmin + 1

      concen_pd(:,:,:,:,ts1) = 0.0_ip

! There are a few strategies here: loop over a combined index for n and
! j, loop over a combined index for n, j and i, or loop first over n,
! then launch opm threads over nested i,j with a collapse flag
!
!      ! Here we bundle the n and j loops into one do loop to more
!      ! efficiently enable parallelization 
!      do idx_dum=1,nsmax*nymax
!      ! Now recover n and j
!        j = mod(idx_dum-1,nymax) +1
!        n = floor(real((idx_dum-1)/nymax,kind=ip))+1
!        if(.not.IsAloft(n)) cycle
!        do i=1,nxmax

!  OR

!      do idx_dum=1,nsmax*nymax*nxmax
!         ! Now recover i,j and n
!         ! Note: we can recover these indecies in any order we choose,
!         ! as long as all are accounted for.  So recover in optimal
!         ! (column-major) order
!         i = (idx_dum-1)/(nymax*nsmax)+1
!         j = (idx_dum - (i-1)*nymax*nsmax-1)/nsmax+1
!         n = idx_dum - (j-1) * nsmax - (i-1)*nymax*nsmax

      do n=1,nsmax
        if(.not.IsAloft(n)) cycle
        !$OMP PARALLEL DO &
        !$OMP DEFAULT(NONE) &
        !$OMP SHARED(n,nxmax,nymax,ncells,nsmax,dt,concen_pd,kappa_pd,&
        !$OMP vz_pd,sigma_nz_pd,outflow_xy1_pd,outflow_xy2_pd,vf_pd),&
        !$OMP PRIVATE(l,i,j,q_cc,vel_cc,dt_vol_cc,usig_I,update_cc,&
        !$OMP dq_I,fs_I,fss_I,ldq_I,dqu_I,i_I,i_cc,&
        !$OMP aus,theta,divu_p,divu_m,&
        !$OMP LFluct_Rbound,RFluct_Lbound,&
        !$OMP LimFlux_Rbound,LimFlux_Lbound),&
        !$OMP collapse(2)
        do j=jmin,jmax
          do i=imin,imax
            ! Initialize cell-centered values for this z-column
            ! Note: ghost cells should contain q_cc=0 and vel_cc=edge
            q_cc(  rmin-2:rmin-1+ncells+2) = concen_pd(i,j,rmin-2:rmin-1+ncells+2,n,ts0)
            vel_cc(rmin-2:rmin-1+ncells+2) =     vz_pd(i,j,rmin-2:rmin-1+ncells+2) + &
                                                 vf_pd(i,j,rmin-2:rmin-1+ncells+2,n)

            ! Ghost cells were set in Set_BC.f90, but could be reset here
            ! if desired or for testing.  Tests showed that velocities
            ! should either have a constant or a linear extrapolation, but
            ! concentrations should be set to zero.

            ! Calculate \Delta t / \kappa
              ! using kappa of cell
            dt_vol_cc(rmin-2:rmin-1+ncells+2) = real(dt,kind=ip) / &
                                                kappa_pd(i,j,rmin-2:rmin-1+ncells+2)

            ! Make sure to initialize this since we are only setting it where is matters
            usig_I = 0.0_ip
            do l=rmin,rmin-1+ncells+1
              usig_I(l) = 0.5_ip*(vel_cc(l-1)+vel_cc(l))*sigma_nz_pd(i,j,l)
            enddo

            ! This calculates the update in a row in one function call
            !update_cc(-1:ncells+2) = AdvectUpdate_1d(ncells,q_cc,dt_vol_cc,usig_I)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! This is the branch that includes the update in-line
            dq_I(rmin-1:rmin-1+ncells+2) = q_cc(rmin-1:rmin-1+ncells+2) - &
                                           q_cc(rmin-2:rmin-1+ncells+1)

            ! First get the limited Delta Q, in we are using high-order
            ! methods
            ldq_I = 0.0_ip
            fs_I( rmin-1:rmin-1+ncells+2,1:2) = 0.0_ip
            fss_I(rmin-1:rmin-1+ncells+2) = 0.0_ip

#ifndef LIM_NONE
            do i_I = rmin,rmin-1+ncells+1
              ! This cycling is good for production runs, but causes
              ! problems with convergence tests
#ifndef NOCYCLE
              if (abs(dq_I(i_I)).le.EPS_THRESH) cycle
#endif

#if LIM_LAXWEN
              ! Lax-Wendrof (linear)
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
              fs_I(i_I,2) = max(0.0_ip,usig_I(i_I))*dq_I(i_I) ! flux OUT OF l-1 cell to l cell
              fs_I(i_I,1) = min(0.0_ip,usig_I(i_I))*dq_I(i_I) ! flux OUT OF l cell to l-1 cell
                ! Modification for conservative form in
                ! divergent/convergent velocities
              divu_p = max(0.0_ip,usig_I(i_I  )) - max(0.0_ip,usig_I(i_I-1))
              divu_m = min(0.0_ip,usig_I(i_I+1)) - min(0.0_ip,usig_I(i_I))

              fs_I(i_I,2) = fs_I(i_I,2) + q_cc(i_I  ) * divu_m
              fs_I(i_I,1) = fs_I(i_I,1) + q_cc(i_I-1) * divu_p
            enddo  ! loop over i_I (interfaces)

            !--------------------------------------------------------
            ! Now loop over the cell indicies and apply these fluxes to
            ! volume i_cc
            do i_cc=rmin-1,rmin-1+ncells+1  ! Note: we additionally loop on 1 ghost to get
                                            !       outflow boundary fluxes

               ! Interface fluctuation/limited-q at left cell interface
              RFluct_Lbound  =  fs_I(i_cc,2)
              LimFlux_Lbound = fss_I(i_cc)
               ! Interface fluctuation/limited-q at right cell interface
              LFluct_Rbound  =  fs_I(i_cc+1,1)
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
            concen_pd(i,j,rmin:rmin-1+ncells,n,ts1) = &
               concen_pd(i,j,rmin:rmin-1+ncells,n,ts0) + &
               update_cc(rmin:rmin-1+ncells)

              if(rmin.eq.1) &
                ! Flux out the - side of advection row  (copied to deposit)
                outflow_xy1_pd(i,j,n) = outflow_xy1_pd(i,j,n) + update_cc(0)
              if(rmax.eq.nzmax) &
                ! Flux out the + side of advection row  (top of domain)
                outflow_xy2_pd(i,j,n) = outflow_xy2_pd(i,j,n) + update_cc(ncells+1)

          enddo
        enddo
        !$OMP END PARALLEL DO
      enddo ! loop over idx_dum

      concen_pd(  1:nxmax,1:nymax,1:nzmax,1:nsmax,ts0) = &
        concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts1)

      DepositGranularity(1:nxmax,1:nymax,1:nsmax)=outflow_xy1_pd(1:nxmax,1:nymax,1:nsmax)

      end subroutine advect_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module AdvectionVert_DCU

!##############################################################################
