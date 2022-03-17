      module Diffusion

      use precis_param

      use io_units

      use global_param,  only : &
         useVarDiffH,useVarDiffV

      implicit none

      real(kind=ip) :: diffusivity_horz    ! horizontal diffusion coefficient (m2/s)
      real(kind=ip) :: diffusivity_vert    ! vertical diffusion coefficient (m2/s)

         ! Factor that controls how much of the
         ! n+1 time step is used in the solution
         !  0.0 = Forward-Euler
         !  0.5 = Crank-Nicolson
         !  1.0 = Backward-Euler
         ! Note: if 0.0 is used, the method is numerically equivalant to EXPLDIFF,
         !       but solved with the lapack solvers with A being the identity
         !       matrix, not as efficient, but useful for checking.  In this
         !       case, Imp_DT_fac would need to be 0.5 as required by the
         !       explicit solver.  If either Imp_fac = 0.5 or 1.0, then the
         !       method is unconditionally stable, but accuracy requires a
         !       Imp_DT_fac to be around 2.0
      real(kind=ip) :: Imp_fac = 0.5_ip    
      real(kind=ip) :: Imp_DT_fac  = 4.0_ip

#ifdef USEPOINTERS
      real(kind=ip),dimension(:,:,:),pointer :: kx
      real(kind=ip),dimension(:,:,:),pointer :: ky
      real(kind=ip),dimension(:,:,:),pointer :: kz
#else
      real(kind=ip),dimension(:,:,:),allocatable:: kx
      real(kind=ip),dimension(:,:,:),allocatable:: ky
      real(kind=ip),dimension(:,:,:),allocatable:: kz
#endif      

      contains

!******************************************************************************

      subroutine Allocate_Diff(nx,ny,nz)

      implicit none

      integer :: nx,ny,nz
      !integer :: ngridnode

      !ngridnode = (nx+2)*(ny+2)*(nz+2)

      ! Initialize diffusivity arrays with diffusivity_horz and diffusivity_vert
      ! These will change if useVarDiff = .true.
      allocate(kx(0:nx+1,0:ny+1,0:nz+1)); kx = diffusivity_horz
      allocate(ky(0:nx+1,0:ny+1,0:nz+1)); ky = diffusivity_horz
      allocate(kz(0:nx+1,0:ny+1,0:nz+1)); kz = diffusivity_vert

      end subroutine Allocate_Diff

!******************************************************************************

      subroutine Deallocate_Diff

      implicit none

      deallocate(kx)
      deallocate(ky)
      deallocate(kz)

      end subroutine Deallocate_Diff

!******************************************************************************

      subroutine DiffuseHorz(i)

      use global_param,  only : &
         useCN,VERB

      implicit none

      integer :: i

      if(useCN)then
        if(mod(i,2).eq.0) then
          if(VERB.gt.1)write(global_info,*)"Ash3d: Calling diffCN_zxy"
          call diffCN_x
          call diffCN_y
        else
          if(VERB.gt.1)write(global_info,*)"Ash3d: Calling diffCN_zyx"
          call diffCN_y
          call diffCN_x
        endif
      else
        if(mod(i,2).eq.0) then
          if(VERB.gt.1)write(global_info,*)"Ash3d: Calling diff_zxy"
          call diff_x
          call diff_y
        else
          if(VERB.gt.1)write(global_info,*)"Ash3d: Calling diff_zyx"
          call diff_y
          call diff_x
        endif
      endif

      end subroutine DiffuseHorz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine DiffuseVert

      use global_param,  only : &
         useCN

      implicit none

      if(useCN)then
        call diffCN_z
      else
        call diff_z
      endif

      end subroutine DiffuseVert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine diff_x

      ! Explicit diffusion routine.
      ! RP Denlinger and HF Schwaiger

      !!!$ USE omp_lib

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,ts0,ts1,kappa_pd,sigma_nx_pd,IsPeriodic

      use solution,      only : &
         concen_pd,IsAloft,imin,imax,jmin,jmax,kmin,kmax

      use time_data,     only : &
         dt

      implicit none

      integer :: j,k,n  ! These are the indices mapping to the global arrays
      integer :: l_I    ! This is the interface index along the particular diffusion direction
      integer :: l_cc   ! This is the cell-centered index along the particular diffusion direction
      integer :: ncells

       ! arrays that live on cell-centers: Note that we only have 1 ghost cell
      real(kind=ip),dimension(0:nxmax+1)     :: update_cc
      real(kind=ip),dimension(0:nxmax+1)     :: q_cc      ! concen
      real(kind=ip),dimension(0:nxmax+1)     :: vol_cc    ! volume of cell

       ! arrays that live on cell interfaces
       !  Note: interface I for cell i is at (i-1/2); i.e. the left or negative side of i
       !        We only need the interfaces up to the boundary of the domain
       !        (not the ghost cells)
           ! These should be only from 1 to ncells + 1
      real(kind=ip),dimension( 1:nxmax+1)     :: sig_I   ! area of interface
      real(kind=ip),dimension( 1:nxmax+1)     :: ds_I    ! length measure used at interface
      real(kind=ip),dimension( 1:nxmax+1)     :: dq_I    ! step in concentration
      real(kind=ip),dimension( 1:nxmax+1)     :: k_ds_I  ! k/ds

      real(kind=ip) :: LFluct_Rbound,RFluct_Lbound
      integer :: rmin, rmax     ! min and max indices of the row

      !integer OMP_GET_MAX_THREADS
      !integer OMP_GET_NUM_THREADS
      !integer OMP_GET_THREAD_NUM
      !integer :: nthreads,thread_num
      !logical :: OMP_get_nested

      ! We are diffusing in x so set the length of the cell list accordingly
      !rmin = imin
      !rmax = imax
      ! Since sub-grid for diffusion has not been tested, use the full length
      rmin = 1
      rmax = nxmax
      ncells = rmax - rmin + 1

      concen_pd(:,:,:,:,ts1) = 0.0_ip
      ! Neuman boundary conditions in x
        !***  Left/Right (X)
      concen_pd(     -1,:,:,:,ts0) = concen_pd(    1,:,:,:,ts0)
      concen_pd(      0,:,:,:,ts0) = concen_pd(    1,:,:,:,ts0)
      concen_pd(nxmax+1,:,:,:,ts0) = concen_pd(nxmax,:,:,:,ts0)
      concen_pd(nxmax+2,:,:,:,ts0) = concen_pd(nxmax,:,:,:,ts0)
      ! Apply periodicity
      if(IsPeriodic)then
        concen_pd(-1     ,:,:,:,ts0) = concen_pd(nxmax-1,:,:,:,ts0)
        concen_pd( 0     ,:,:,:,ts0) = concen_pd(nxmax  ,:,:,:,ts0)
        concen_pd(nxmax+1,:,:,:,ts0) = concen_pd(1      ,:,:,:,ts0)
        concen_pd(nxmax+2,:,:,:,ts0) = concen_pd(2      ,:,:,:,ts0)
      endif

      do n=1,nsmax
        if(.not.IsAloft(n)) cycle
      !!!$OMP PARALLEL DO &
      !!!$OMP DEFAULT(NONE) &
      !!!$OMP SHARED(n,nymax,nzmax,ncells,nsmax,dt,concen_pd,kappa_pd,&
      !!!$OMP sigma_nx_pd,kx,&
      !!!$OMP IsPeriodic),&
      !!!$OMP PRIVATE(l,j,k,q_cc,update_cc,ds,k_ds2_I,&
      !!!$OMP dq_I,l_cc,k_ds2_I)&
      !!!$OMP collapse(2)
        do k=kmin,kmax
          do j=jmin,jmax
            ! Initialize cell-centered values for this x-row
            ! Note: ghost cells should contain q_cc values at edge (Neumann)
            update_cc(0:ncells+1) = 0.0_ip
            vol_cc(rmin-1:rmin-1+ncells+1) = kappa_pd(   rmin-1:rmin-1+ncells+1,j,k)
            q_cc(  rmin-1:rmin-1+ncells+1) = concen_pd(  rmin-1:rmin-1+ncells+1,j,k,n,ts0)
            sig_I( rmin  :rmin-1+ncells+1) = sigma_nx_pd(rmin  :rmin-1+ncells+1,j,k)
            dq_I(  rmin  :rmin-1+ncells+1) = q_cc(       rmin  :rmin-1+ncells+1) - &
                                             q_cc(       rmin-1:rmin-1+ncells)

              ! Loop over interfaces and get geometry term
            do l_I = rmin,rmin-1+ncells+1
              l_cc = l_I
                ! ds is the 1/dx for this cell along x
                ! For approximating the dx with respect to the gradient of k, we
                ! use a symmetric ds with the area of the interface and the
                ! average of the volumes across the interface
              ds_I(l_I) = 2.0_ip*sig_I(l_I)/(vol_cc(l_cc-1)+vol_cc(l_cc))
                ! get an average diffusivity using the arithmetic average
              k_ds_I(l_I) = 0.5_ip*(kx(l_cc-1,j,k)+kx(l_cc,j,k))*ds_I(l_I)
            enddo
              ! Loop over cells and update
            do l_cc=rmin,rmin-1+ncells
              l_I = l_cc
                ! Eq 4.11 LeVeque02
                ! Note that the ds used for the diffusive flux uses the
                ! interface of the flux, but the kappa of the updated cell
              LFluct_Rbound = dt*k_ds_I(l_I+1)*dq_I(l_I+1)* &
                               (sig_I(l_I+1)/vol_cc(l_cc))
              RFluct_Lbound = -dt*k_ds_I(l_I  )*dq_I(l_I  )* &
                               (sig_I(l_I)/vol_cc(l_cc))

              update_cc(l_cc) = LFluct_Rbound + RFluct_Lbound

            enddo ! loop over l (cell centers)

            concen_pd(   rmin:rmin-1+ncells,j,k,n,ts1) = &
               concen_pd(rmin:rmin-1+ncells,j,k,n,ts0) + &
               update_cc(rmin:rmin-1+ncells)

          enddo ! 
        enddo ! loop over j=1,nymax
     !!!! !$OMP END PARALLEL do

      enddo ! loop over idx_dum

      ! Note: need to fix diffusion for periodic grids

      concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts0) = &
        concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts1)

      end subroutine diff_x

!******************************************************************************

      subroutine diff_y

      ! Explicit diffusion routine.
      ! RP Denlinger and HF Schwaiger

      !!!$ USE omp_lib

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,ts0,ts1,kappa_pd,sigma_ny_pd

      use solution,      only : &
         concen_pd,IsAloft,imin,imax,jmin,jmax,kmin,kmax

      use time_data,     only : &
         dt

      implicit none

      integer :: i,k,n  ! These are the indices mapping to the global arrays
      integer :: l_I    ! This is the interface index along the particular diffusion direction
      integer :: l_cc   ! This is the cell-centered index along the particular diffusion direction
      integer :: ncells

       ! arrays that live on cell-centers: Note that we only have 1 ghost cell
      real(kind=ip),dimension(0:nymax+1)     :: update_cc
      real(kind=ip),dimension(0:nymax+1)     :: q_cc      ! concen
      real(kind=ip),dimension(0:nymax+1)     :: vol_cc    ! volume of cell

       ! arrays that live on cell interfaces
       !  Note: interface I for cell i is at (i-1/2); i.e. the left or negative side of i
       !        We only need the interfaces up to the boundary of the domain
       !        (not the ghost cells)
           ! These should be only from 1 to ncells + 1
      real(kind=ip),dimension( 1:nymax+1)     :: sig_I   ! area of interface
      real(kind=ip),dimension( 1:nymax+1)     :: ds_I    ! length measure usedbat interface
      real(kind=ip),dimension( 1:nymax+1)     :: dq_I    ! step in concentration
      real(kind=ip),dimension( 1:nymax+1)     :: k_ds_I  ! k/ds

      real(kind=ip) :: LFluct_Rbound,RFluct_Lbound
      integer :: rmin, rmax     ! min and max indices of the row

      !integer OMP_GET_MAX_THREADS
      !integer OMP_GET_NUM_THREADS
      !integer OMP_GET_THREAD_NUM
      !integer :: nthreads,thread_num
      !logical :: OMP_get_nested

      ! We are diffusing in y so set the length of the cell list accordingly
      !rmin = jmin
      !rmax = jmax
      ! Since sub-grid for diffusion has not been tested, use the full length
      rmin = 1
      rmax = nymax
      ncells = rmax - rmin + 1

      concen_pd(:,:,:,:,ts1) = 0.0_ip
      ! Neuman boundary conditions in y
        !***  Up/Down (Y)
      concen_pd(:,     -1,:,:,ts0) = concen_pd(:,    1,:,:,ts0)
      concen_pd(:,      0,:,:,ts0) = concen_pd(:,    1,:,:,ts0)
      concen_pd(:,nymax+1,:,:,ts0) = concen_pd(:,nymax,:,:,ts0)
      concen_pd(:,nymax+2,:,:,ts0) = concen_pd(:,nymax,:,:,ts0)

      do n=1,nsmax
        if(.not.IsAloft(n)) cycle
      !!!$OMP PARALLEL DO &
      !!!$OMP DEFAULT(NONE) &
      !!!$OMP SHARED(n,nxmax,nzmax,ncells,nsmax,dt,concen_pd,kappa_pd,&
      !!!$OMP sigma_ny_pd,ky,&
      !!!$OMP IsPeriodic),&
      !!!$OMP PRIVATE(l,j,k,q_cc,update_cc,ds,k_ds2_I,&
      !!!$OMP dq_I,l_cc,k_ds2_I)&
      !!!$OMP collapse(2)
        do k=kmin,kmax
          do i=imin,imax
            ! Initialize cell-centered values for this y-row
            ! Note: ghost cells should contain q_cc values at edge (Neumann)
            update_cc(0:ncells+1) = 0.0_ip
            vol_cc(rmin-1:rmin-1+ncells+1) = kappa_pd(   i,rmin-1:rmin-1+ncells+1,k)
            q_cc(  rmin-1:rmin-1+ncells+1) = concen_pd(  i,rmin-1:rmin-1+ncells+1,k,n,ts0)
            sig_I(rmin:rmin-1+ncells+1)    = sigma_ny_pd(i,rmin  :rmin-1+ncells+1,k)
            dq_I(rmin:rmin-1+ncells+1)     = q_cc(         rmin  :rmin-1+ncells+1) - &
                                             q_cc(         rmin-1:rmin-1+ncells)

              ! Loop over interfaces and get geometry term
            do l_I = rmin,rmin-1+ncells+1
              l_cc = l_I
                ! ds is the 1/dy for this cell along y
                ! For approximating the dy with respect to the gradient of k, we
                ! use a symmetric ds with the area of the interface and the
                ! average of the volumes across the interface
              ds_I(l_I) = 2.0_ip*sig_I(l_I)/(vol_cc(l_cc-1)+vol_cc(l_cc))
                ! get an average diffusivity using the arithmetic average
              k_ds_I(l_I) = 0.5_ip*(ky(i,l_cc-1,k)+ky(i,l_cc,k))*ds_I(l_I)
            enddo
              ! Loop over cells and update
            do l_cc=rmin,rmin-1+ncells
              l_I = l_cc
                ! Eq 4.11 LeVeque02
                ! Note that the ds used for the diffusive flux uses the
                ! interface of the flux, but the kappa of the updated cell
              LFluct_Rbound = dt*k_ds_I(l_I+1)*dq_I(l_I+1)* &
                               (sig_I(l_I+1)/vol_cc(l_cc))
              RFluct_Lbound = -dt*k_ds_I(l_I  )*dq_I(l_I  )* &
                               (sig_I(l_I)/vol_cc(l_cc))

              update_cc(l_cc) = LFluct_Rbound + RFluct_Lbound
            enddo ! loop over l (cell centers)

            concen_pd(   i,rmin:rmin-1+ncells,k,n,ts1) = &
               concen_pd(i,rmin:rmin-1+ncells,k,n,ts0) + &
               update_cc(  rmin:rmin-1+ncells)

          enddo ! 
        enddo ! loop over j=1,nymax
     !!!! !$OMP END PARALLEL do

      enddo ! loop over idx_dum

      concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts0) = &
        concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts1)

      end subroutine diff_y

!******************************************************************************

      subroutine diff_z

      ! Explicit diffusion routine.
      ! RP Denlinger and HF Schwaiger

      !!!$ USE omp_lib

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,ts0,ts1,kappa_pd,sigma_nz_pd

      use solution,      only : &
         concen_pd,IsAloft,imin,imax,jmin,jmax,kmin,kmax

      use time_data,     only : &
         dt

      implicit none

      integer :: i,j,n  ! These are the indices mapping to the global arrays
      integer :: l_I    ! This is the interface index along the particular diffusion direction
      integer :: l_cc   ! This is the cell-centered index along the particular diffusion direction
      integer :: ncells

       ! arrays that live on cell-centers: Note that we only have 1 ghost cell
      real(kind=ip),dimension(0:nzmax+1)     :: update_cc
      real(kind=ip),dimension(0:nzmax+1)     :: q_cc      ! concen
      real(kind=ip),dimension(0:nzmax+1)     :: vol_cc    ! volume of cell

       ! arrays that live on cell interfaces
       !  Note: interface I for cell i is at (i-1/2); i.e. the left or negative side of i
       !        We only need the interfaces up to the boundary of the domain
       !        (not the ghost cells)
           ! These should be only from 1 to ncells + 1
      real(kind=ip),dimension( 1:nzmax+1)     :: sig_I   ! area of interface
      real(kind=ip),dimension( 1:nzmax+1)     :: ds_I    ! length measure used at interface
      real(kind=ip),dimension( 1:nzmax+1)     :: dq_I    ! step in concentration
      real(kind=ip),dimension( 1:nzmax+1)     :: k_ds_I  ! k/ds

      real(kind=ip) :: LFluct_Rbound,RFluct_Lbound
      integer :: rmin, rmax     ! min and max indices of the row

      !integer OMP_GET_MAX_THREADS
      !integer OMP_GET_NUM_THREADS
      !integer OMP_GET_THREAD_NUM
      !integer :: nthreads,thread_num
      !logical :: OMP_get_nested

      ! We are diffusing in z so set the length of the cell list accordingly
      !rmin = kmin
      !rmax = kmax
      ! Since sub-grid for diffusion has not been tested, use the full length
      rmin = 1
      rmax = nzmax
      ncells = rmax - rmin + 1

      concen_pd(:,:,:,:,ts1) = 0.0_ip
      ! Neuman boundary conditions in z
        !***  Bottom/Top (Z)
      concen_pd(:,:,     -1,:,ts0) = concen_pd(:,:,    1,:,ts0)
      concen_pd(:,:,      0,:,ts0) = concen_pd(:,:,    1,:,ts0)
      concen_pd(:,:,nzmax+1,:,ts0) = concen_pd(:,:,nzmax,:,ts0)
      concen_pd(:,:,nzmax+2,:,ts0) = concen_pd(:,:,nzmax,:,ts0)

      do n=1,nsmax
        if(.not.IsAloft(n)) cycle
      !!!$OMP PARALLEL DO &
      !!!$OMP DEFAULT(NONE) &
      !!!$OMP SHARED(n,nymax,nzmax,ncells,nsmax,dt,concen_pd,kappa_pd,&
      !!!$OMP sigma_nz_pd,kz,&
      !!!$OMP IsPeriodic),&
      !!!$OMP PRIVATE(l,j,k,q_cc,update_cc,ds,k_ds2_I,&
      !!!$OMP dq_I,l_cc,k_ds2_I)&
      !!!$OMP collapse(2)
        do j=jmin,jmax
          do i=imin,imax
            ! Initialize cell-centered values for this z-row
            ! Note: ghost cells should contain q_cc values at edge (Neumann)
            update_cc(0:ncells+1) = 0.0_ip
            vol_cc(rmin-1:rmin-1+ncells+1) = kappa_pd(   i,j,rmin-1:rmin-1+ncells+1)
            q_cc(  rmin-1:rmin-1+ncells+1) = concen_pd(  i,j,rmin-1:rmin-1+ncells+1,n,ts0)
            sig_I( rmin  :rmin-1+ncells+1) = sigma_nz_pd(i,j,rmin  :rmin-1+ncells+1)
            dq_I(  rmin  :rmin-1+ncells+1) = q_cc(           rmin  :rmin-1+ncells+1) - &
                                             q_cc(           rmin-1:rmin-1+ncells)

              ! Loop over interfaces and get geometry term
            do l_I = rmin,rmin-1+ncells+1
              l_cc = l_I
                ! ds is the 1/dz for this cell along z
                ! For approximating the dz with respect to the gradient of k, we
                ! use a symmetric ds with the area of the interface and the
                ! average of the volumes across the interface
              ds_I(l_I) = 2.0_ip*sig_I(l_I)/(vol_cc(l_cc-1)+vol_cc(l_cc))
                ! get an average diffusivity using the arithmetic average
              k_ds_I(l_I) = 0.5_ip*(kz(i,j,l_cc-1)+kz(i,j,l_cc))*ds_I(l_I)
            enddo

              ! Loop over cells and update
            do l_cc=rmin,rmin-1+ncells
              l_I = l_cc
                ! Eq 4.11 LeVeque02
                ! Note that the ds used for the diffusive flux uses the
                ! interface of the flux, but the kappa of the updated cell
              LFluct_Rbound = dt*k_ds_I(l_I+1)*dq_I(l_I+1)* &
                               (sig_I(l_I+1)/vol_cc(l_cc))
              RFluct_Lbound = -dt*k_ds_I(l_I  )*dq_I(l_I  )* &
                               (sig_I(l_I)/vol_cc(l_cc))

              update_cc(l_cc) = LFluct_Rbound + RFluct_Lbound

            enddo ! loop over l (cell centers)

            concen_pd(i,j,rmin:rmin-1+ncells,n,ts1) = &
               concen_pd(i,j,rmin:rmin-1+ncells,n,ts0) + &
               update_cc(rmin:rmin-1+ncells)

          enddo ! 
        enddo ! loop over j=1,nymax
     !!!! !$OMP END PARALLEL do

      enddo ! loop over idx_dum

      concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts0) = &
        concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts1)

      end subroutine diff_z

!******************************************************************************

      subroutine diffCN_x

      ! Implicit Crank-Nicolson diffusion routine
      ! Implements Eq 4.13 of LeVeque02

      use mesh,          only : &
         IsLatLon,nxmax,nymax,nzmax,nsmax,ts0,ts1,kappa_pd,sigma_nx_pd,IsPeriodic

      use solution,      only : &
         concen_pd,IsAloft,imin,imax,jmin,jmax,kmin,kmax

      use time_data,     only : &
         dt

      implicit none

      integer :: j,k,n  ! These are the indices mapping to the global arrays
      integer :: l_I    ! This is the interface index along the particular diffusion direction
      integer :: l_cc   ! This is the cell-centered index along the particular diffusion direction
      integer :: ncells

      real(kind=ip) :: BC_left_t0,BC_left_t1
      real(kind=ip) :: BC_right_t0,BC_right_t1
      real(kind=sp),allocatable,dimension(:) :: DL_s,D_s,DU_s,B_s
      real(kind=dp),allocatable,dimension(:) :: DL_d,D_d,DU_d,B_d
      integer :: nlineq,nrhs,ldb,info
      real(kind=ip) :: LeftFac,CenterFac,RightFac

       ! arrays that live on cell-centers: Note that we only have 1 ghost cell
      real(kind=ip),dimension(0:nxmax+1)     :: q_cc      ! concen
      real(kind=ip),dimension(0:nxmax+1)     :: vol_cc    ! volume of cell

       ! arrays that live on cell interfaces
       !  Note: interface I for cell i is at (i-1/2); i.e. the left or negative
       !  side of i
       !        We only need the interfaces up to the boundary of the domain
       !        (not the ghost cells)
           ! These should be only from 1 to ncells + 1
      real(kind=ip),dimension( 1:nxmax+1)     :: sig_I   ! area of interface
      real(kind=ip),dimension( 1:nxmax+1)     :: ds_I    ! length measure used at interface
      real(kind=ip),dimension( 1:nxmax+1)     :: vavg_I
      real(kind=ip),dimension( 1:nxmax+1)     :: kavg_I
      real(kind=ip),dimension( 1:nxmax+1)     :: k_ds_I  ! k/ds
!      real(kind=ip),dimension( 1:nxmax+1)     :: ksig2_vol_I  ! k*sig*sig/volavg 

      integer :: rmin, rmax     ! min and max indices of the row

#ifdef CRANKNIC
      ! Note: The only reason not to use Crank-Nicolson is if you
      !       don't have blas and lapack installed.  This pre-proc.
      !       directive allows this section to be turned off.
      INTERFACE
        subroutine sgtsv(N,NRHS,DL,D,DU,B,LDB,INFO)
          integer                         ,intent(in)    :: N
          integer                         ,intent(in)    :: NRHS
          real(kind=4),dimension(N-1)     ,intent(inout) :: DL
          real(kind=4),dimension(N)       ,intent(inout) :: D
          real(kind=4),dimension(N-1)     ,intent(inout) :: DU
          integer                         ,intent(in)    :: LDB
          real(kind=4),dimension(ldb,NRHS),intent(inout) :: B
          integer                         ,intent(out)   :: INFO
        end subroutine
        subroutine dgtsv(N,NRHS,DL,D,DU,B,LDB,INFO)
          integer                         ,intent(in)    :: N
          integer                         ,intent(in)    :: NRHS
          real(kind=8),dimension(N-1)     ,intent(inout) :: DL
          real(kind=8),dimension(N)       ,intent(inout) :: D
          real(kind=8),dimension(N-1)     ,intent(inout) :: DU
          integer                         ,intent(in)    :: LDB
          real(kind=8),dimension(ldb,NRHS),intent(inout) :: B
          integer                         ,intent(out)   :: INFO
        end subroutine
        subroutine sptsv(N,NRHS,D,E,B,LDB,INFO)
          integer                         ,intent(in)    :: N
          integer                         ,intent(in)    :: NRHS
          real(kind=4),dimension(N)       ,intent(inout) :: D
          real(kind=4),dimension(N-1)     ,intent(inout) :: E
          integer                         ,intent(in)    :: LDB
          real(kind=4),dimension(ldb,NRHS),intent(inout) :: B
          integer                         ,intent(out)   :: INFO
        end subroutine
        subroutine dptsv(N,NRHS,D,E,B,LDB,INFO)
          integer                         ,intent(in)    :: N
          integer                         ,intent(in)    :: NRHS
          real(kind=8),dimension(N)       ,intent(inout) :: D
          real(kind=8),dimension(N-1)     ,intent(inout) :: E
          integer                         ,intent(in)    :: LDB
          real(kind=8),dimension(ldb,NRHS),intent(inout) :: B
          integer                         ,intent(out)   :: INFO
        end subroutine
      END INTERFACE
#endif

      ! We are diffusing in x so set the length of the cell list accordingly
      !rmin = imin
      !rmax = imax
      ! Since sub-grid for diffusion has not been tested, use the full length
      rmin = 1
      rmax = nxmax
      ncells = rmax - rmin + 1

      concen_pd(:,:,:,:,ts1) = 0.0_ip
      ! Neuman boundary conditions in x
        !***  Left/Right (X)
      concen_pd(     -1,:,:,:,ts0) = concen_pd(    1,:,:,:,ts0)
      concen_pd(      0,:,:,:,ts0) = concen_pd(    1,:,:,:,ts0)
      concen_pd(nxmax+1,:,:,:,ts0) = concen_pd(nxmax,:,:,:,ts0)
      concen_pd(nxmax+2,:,:,:,ts0) = concen_pd(nxmax,:,:,:,ts0)
      ! Apply periodicity
      if(IsPeriodic)then
        concen_pd(-1     ,:,:,:,ts0) = concen_pd(nxmax-1,:,:,:,ts0)
        concen_pd( 0     ,:,:,:,ts0) = concen_pd(nxmax  ,:,:,:,ts0)
        concen_pd(nxmax+1,:,:,:,ts0) = concen_pd(1      ,:,:,:,ts0)
        concen_pd(nxmax+2,:,:,:,ts0) = concen_pd(2      ,:,:,:,ts0)
      endif

      if(ncells.gt.1)then

      nlineq = ncells
      ldb = nlineq  ! leading dimension of b is num of equations
      if (ip.eq.4)then
        allocate(DL_s(nlineq-1));
        allocate(D_s(nlineq));
        allocate(DU_s(nlineq-1));
        allocate(B_s(nlineq));
      endif
      allocate(DL_d(nlineq-1));
      allocate(D_d(nlineq));
      allocate(DU_d(nlineq-1));
      allocate(B_d(nlineq));

      do n=1,nsmax
        if(.not.IsAloft(n)) cycle

        do k=kmin,kmax
          do j=jmin,jmax
            ! solve the problem in x for each y
            vol_cc(rmin-1:rmin-1+ncells+1) = kappa_pd(      rmin-1:rmin-1+ncells+1,j,k)
            q_cc(  rmin-1:rmin-1+ncells+1) = concen_pd(     rmin-1:rmin-1+ncells+1,j,k,n,ts0)
            sig_I( rmin  :rmin-1+ncells+1) = sigma_nx_pd(   rmin  :rmin-1+ncells+1,j,k)
            vavg_I(rmin  :rmin-1+ncells+1) = 0.5_ip*(vol_cc(rmin-1:rmin-1+ncells ) + &
                                                     vol_cc(rmin  :rmin-1+ncells+1))
            kavg_I(rmin  :rmin-1+ncells+1) = 0.5_ip*(kx(    rmin-1:rmin-1+ncells  ,j,k) + &
                                                     kx(    rmin  :rmin-1+ncells+1,j,k))
            ds_I(  rmin  :rmin-1+ncells+1) = sig_I(         rmin  :rmin-1+ncells+1) / &
                                            vavg_I(         rmin  :rmin-1+ncells+1)
            k_ds_I(rmin  :rmin-1+ncells+1) = kavg_I(        rmin  :rmin-1+ncells+1)* &
                                                    ds_I(   rmin  :rmin-1+ncells+1)

            ! Note:  nrhs will be 1 in this case, but we can apply this
            ! to all rows at once using nrhs = j*k*n
            nrhs = 1
            ! Loop over all cells in this x-row
            do l_cc=rmin,rmin-1+ncells
              l_I = l_cc  ! Interface ID refers to the Left side of the cell (k-1/2)
              LeftFac    = dt*k_ds_I(l_I  )*sig_I(l_I  )/vol_cc(l_cc)
              RightFac   = dt*k_ds_I(l_I+1)*sig_I(l_I+1)/vol_cc(l_cc)
              CenterFac = LeftFac + RightFac

              if(l_cc.eq.1) then
                !DL_d = nothing     :: No lower diagonal for first row
                D_d(l_cc)  = 1.0_ip + (Imp_fac)*CenterFac & ! This is part of the normal stencil
                                     -(Imp_fac)*LeftFac     ! This line is the BC that ensures
                                                            ! that q(ncell)=q(ncell+1) at t1
                                                            ! i.e. no outward flux
                DU_d(l_cc) =        - (Imp_fac)*RightFac

                  ! RHS contains left boundary term
                BC_left_t0 = q_cc(1)   ! This sets a Neuman condition at t0
                !BC_left_t1 ~= q_cc(1) ! This condition is baked into the matrix
                                       ! stencil so that BC_left_t1=q_cc(1) at t1
                B_d(l_cc)  =        (1.0_ip-Imp_fac)*LeftFac    * BC_left_t0 + &
                          (1.0_ip - (1.0_ip-Imp_fac)*CenterFac) * q_cc(l_cc) + &
                                    (1.0_ip-Imp_fac)*RightFac   * q_cc(l_cc+1)
              elseif(l_cc.lt.ncells)then
                DL_d(l_cc-1) =         - (Imp_fac)*LeftFac
                D_d(l_cc)  = 1.0_ip + (Imp_fac)*CenterFac
                DU_d(l_cc)   =         - (Imp_fac)*RightFac

                B_d(l_cc)    =           (1.0_ip-Imp_fac)*LeftFac    * q_cc(l_cc-1) + &
                               (1.0_ip - (1.0_ip-Imp_fac)*CenterFac) * q_cc(l_cc  ) + &
                                         (1.0_ip-Imp_fac)*RightFac   * q_cc(l_cc+1)
              elseif(l_cc.eq.ncells)then
                DL_d(l_cc-1) =         - (Imp_fac)*LeftFac
                D_d(l_cc)  = 1.0_ip + (Imp_fac)*CenterFac  & ! This is part of the normal stencil
                                     -(Imp_fac)*RightFac     ! This line is the BC that ensures
                                                             ! that q(ncell)=q(ncell+1) at t1
                                                             ! i.e. no outward flux
                !DU_d = nothing   :: No upper diagonal for last row

                  ! RHS contains right boundary term
                BC_right_t0 = q_cc(ncells)    ! This sets a Neuman condition at t0
                !BC_right_t1 ~= q_cc(ncells)  ! This condition is baked into the matrix
                                              ! stencil so that BC_right_t1=q_cc(ncells) at t1
                B_d(l_cc)    =           (1.0_ip-Imp_fac)*LeftFac    * q_cc(l_cc-1) + &
                               (1.0_ip - (1.0_ip-Imp_fac)*CenterFac) * q_cc(l_cc ) + &
                                         (1.0_ip-Imp_fac)*RightFac   * BC_right_t0
              endif
 
            enddo
#ifdef CRANKNIC
      ! Note: The only reason not to use Crank-Nicolson is if you
      !       don't have blas and lapack installed.  This pre-proc.
      !       directive allows this section to be turned off.
            if(useVarDiffH.or.IsLatLon)then
              ! This is the call for solving single or double
              ! precision general tridiagonal Ax=b
              ! Note: This is the function to call if kx or vol is spatially
              ! variable
              if(ip.eq.4)then
                DL_s = real(DL_d,kind=4)
                D_s  = real(D_d ,kind=4)
                DU_s = real(DU_d,kind=4)
                B_s  = real(B_d ,kind=4)
                call sgtsv(     &
                      nlineq, &  !i The order of the matrix A.  N >= 0.
                      nrhs,   &  !i The number of right hand sides
                      DL_s,   &  !b array, dimension (N-1) sub-diagonal
                      D_s,    &  !b dimension (N) diagonal
                      DU_s,   &  !b array, dimension (N-1) super-diagonal
                      B_s,    &  !b dimension (LDB,NRHS) On entry, the N by NRHS matrix of right hand side matrix B.
                      ldb,    &  !i The leading dimension of the array B. LDB >= max(1,N)
                      info)      !o
              elseif(ip.eq.8)then
                call dgtsv(     &
                      nlineq, &  !i The order of the matrix A.  N >= 0.
                      nrhs,   &  !i The number of right hand sides
                      DL_d,   &  !b array, dimension (N-1) sub-diagonal
                      D_d,    &  !b dimension (N) diagonal
                      DU_d,   &  !b array, dimension (N-1) super-diagonal
                      B_d,    &  !b dimension (LDB,NRHS) On entry, the N by NRHS matrix of right hand side matrix B.
                      ldb,    &  !i The leading dimension of the array B. LDB >= max(1,N)
                      info)      !o
              endif
            else
              ! This is the call for solving single or double
              ! precision symmetric positive definite tridiagonal Ax=b
              ! Note: A will only be symmetric if kx and vol are homogeneous
              !       This is really not much faster than dgtsv
              if(ip.eq.4)then
                D_s  = real(D_d ,kind=4)
                DU_s = real(DU_d,kind=4)
                B_s  = real(B_d ,kind=4)
                call sptsv(     &
                      nlineq, &  !i The order of the matrix A.  N >= 0.
                      nrhs,   &  !i The number of right hand sides
                      D_s,    &  !b dimension (N) diagonal
                      DU_s,   &  !b array, dimension (N-1) sub or super-diagonal
                      B_s,    &  !b dimension (LDB,NRHS) On entry, the N by NRHS matrix of right hand side matrix B.
                      ldb,    &  !i The leading dimension of the array B. LDB >= max(1,N)
                      info)      !o
              elseif(ip.eq.8)then
                call dptsv(     &
                      nlineq, &  !i The order of the matrix A.  N >= 0.
                      nrhs,   &  !i The number of right hand sides
                      D_d,    &  !b dimension (N) diagonal
                      DU_d,   &  !b array, dimension (N-1) sub or super-diagonal
                      B_d,    &  !b dimension (LDB,NRHS) On entry, the N by NRHS matrix of right hand side matrix B.
                      ldb,    &  !i The leading dimension of the array B. LDB >= max(1,N)
                      info)      !o
              endif
            endif
#endif
            concen_pd(1:nxmax,j,k,n,ts1) = B_d

          enddo ! loop over j
        enddo ! loop over k
      enddo ! loop over n

      if (ip.eq.4)then
        deallocate(DL_s,D_s,DU_s,B_s)
      endif
      deallocate(DL_d,D_d,DU_d,B_d)

      concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts0) = &
        concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts1)

      endif

      end subroutine diffCN_x

!******************************************************************************

      subroutine diffCN_y

      ! Implicit Crank-Nicolson diffusion routine
      ! Implements Eq 4.13 of LeVeque02

      use mesh,          only : &
         IsLatLon,nxmax,nymax,nzmax,nsmax,ts0,ts1,kappa_pd,sigma_ny_pd

      use solution,      only : &
         concen_pd,IsAloft,imin,imax,jmin,jmax,kmin,kmax

      use time_data,     only : &
         dt

      implicit none

      integer :: i,k,n  ! These are the indices mapping to the global arrays
      integer :: l_I    ! This is the interface index along the particular diffusion direction
      integer :: l_cc   ! This is the cell-centered index along the particular diffusion direction
      integer :: ncells

      real(kind=ip) :: BC_left_t0,BC_left_t1
      real(kind=ip) :: BC_right_t0,BC_right_t1
      real(kind=sp),allocatable,dimension(:) :: DL_s,D_s,DU_s,B_s
      real(kind=dp),allocatable,dimension(:) :: DL_d,D_d,DU_d,B_d
      integer :: nlineq,nrhs,ldb,info
      real(kind=ip) :: LeftFac,CenterFac,RightFac

       ! arrays that live on cell-centers: Note that we only have 1 ghost cell
      real(kind=ip),dimension(0:nymax+1)     :: q_cc      ! concen
      real(kind=ip),dimension(0:nymax+1)     :: vol_cc    ! volume of cell

       ! arrays that live on cell interfaces
       !  Note: interface I for cell i is at (i-1/2); i.e. the left or negative
       !  side of i
       !        We only need the interfaces up to the boundary of the domain
       !        (not the ghost cells)
           ! These should be only from 1 to ncells + 1
      real(kind=ip),dimension( 1:nymax+1)     :: sig_I   ! area of interface
      real(kind=ip),dimension( 1:nymax+1)     :: ds_I    ! length measure used at interface
      real(kind=ip),dimension( 1:nymax+1)     :: vavg_I
      real(kind=ip),dimension( 1:nymax+1)     :: kavg_I
      real(kind=ip),dimension( 1:nymax+1)     :: k_ds_I  ! k/ds

      integer :: rmin, rmax     ! min and max indices of the row

#ifdef CRANKNIC
      ! Note: The only reason not to use Crank-Nicolson is if you
      !       don't have blas and lapack installed.  This pre-proc.
      !       directive allows this section to be turned off.
      INTERFACE
        subroutine sgtsv(N,NRHS,DL,D,DU,B,LDB,INFO)
          integer                         ,intent(in)    :: N
          integer                         ,intent(in)    :: NRHS
          real(kind=4),dimension(N-1)     ,intent(inout) :: DL
          real(kind=4),dimension(N)       ,intent(inout) :: D
          real(kind=4),dimension(N-1)     ,intent(inout) :: DU
          integer                         ,intent(in)    :: LDB
          real(kind=4),dimension(ldb,NRHS),intent(inout) :: B
          integer                         ,intent(out)   :: INFO
        end subroutine
        subroutine dgtsv(N,NRHS,DL,D,DU,B,LDB,INFO)
          integer                         ,intent(in)    :: N
          integer                         ,intent(in)    :: NRHS
          real(kind=8),dimension(N-1)     ,intent(inout) :: DL
          real(kind=8),dimension(N)       ,intent(inout) :: D
          real(kind=8),dimension(N-1)     ,intent(inout) :: DU
          integer                         ,intent(in)    :: LDB
          real(kind=8),dimension(ldb,NRHS),intent(inout) :: B
          integer                         ,intent(out)   :: INFO
        end subroutine
        subroutine sptsv(N,NRHS,D,E,B,LDB,INFO)
          integer                         ,intent(in)    :: N
          integer                         ,intent(in)    :: NRHS
          real(kind=4),dimension(N)       ,intent(inout) :: D
          real(kind=4),dimension(N-1)     ,intent(inout) :: E
          integer                         ,intent(in)    :: LDB
          real(kind=4),dimension(ldb,NRHS),intent(inout) :: B
          integer                         ,intent(out)   :: INFO
        end subroutine
        subroutine dptsv(N,NRHS,D,E,B,LDB,INFO)
          integer                         ,intent(in)    :: N
          integer                         ,intent(in)    :: NRHS
          real(kind=8),dimension(N)       ,intent(inout) :: D
          real(kind=8),dimension(N-1)     ,intent(inout) :: E
          integer                         ,intent(in)    :: LDB
          real(kind=8),dimension(ldb,NRHS),intent(inout) :: B
          integer                         ,intent(out)   :: INFO
        end subroutine
      END INTERFACE
#endif

      ! We are diffusing in y so set the length of the cell list accordingly
      !rmin = jmin
      !rmax = jmax
      ! Since sub-grid for diffusion has not been tested, use the full length
      rmin = 1
      rmax = nymax
      ncells = rmax - rmin + 1

      concen_pd(:,:,:,:,ts1) = 0.0_ip
      ! Neuman boundary conditions in y
        !***  Left/Right (Y)
      concen_pd(:,     -1,:,:,ts0) = concen_pd(:,    1,:,:,ts0)
      concen_pd(:,      0,:,:,ts0) = concen_pd(:,    1,:,:,ts0)
      concen_pd(:,nymax+1,:,:,ts0) = concen_pd(:,nymax,:,:,ts0)
      concen_pd(:,nymax+2,:,:,ts0) = concen_pd(:,nymax,:,:,ts0)

      if(ncells.gt.1)then

      nlineq = nymax
      ldb = nlineq  ! leading dimension of b is num of equations
      if (ip.eq.4)then
        allocate(DL_s(nlineq-1));
        allocate(D_s(nlineq));
        allocate(DU_s(nlineq-1));
        allocate(B_s(nlineq));
      endif
      allocate(DL_d(nlineq-1));
      allocate(D_d(nlineq));
      allocate(DU_d(nlineq-1));
      allocate(B_d(nlineq));

      do n=1,nsmax
        if(.not.IsAloft(n)) cycle

        do k=kmin,kmax
          do i=imin,imax
            ! solve the problem in y for each x
            vol_cc(rmin-1:rmin-1+ncells+1) = kappa_pd(    i,rmin-1:rmin-1+ncells+1,k)
            q_cc(  rmin-1:rmin-1+ncells+1) = concen_pd(   i,rmin-1:rmin-1+ncells+1,k,n,ts0)
            sig_I( rmin  :rmin-1+ncells+1) = sigma_ny_pd( i,rmin  :rmin-1+ncells+1,k)
            vavg_I(rmin  :rmin-1+ncells+1) = 0.5_ip*(vol_cc(rmin-1:rmin-1+ncells  ) + &
                                                     vol_cc(rmin  :rmin-1+ncells+1))
            kavg_I(rmin  :rmin-1+ncells+1) = 0.5_ip*(ky(  i,rmin-1:rmin-1+ncells  ,k) + &
                                                     ky(  i,rmin  :rmin-1+ncells+1,k))
            ds_I(  rmin  :rmin-1+ncells+1) = sig_I(         rmin  :rmin-1+ncells+1) / &
                                            vavg_I(         rmin  :rmin-1+ncells+1)
            k_ds_I(rmin  :rmin-1+ncells+1) = kavg_I(        rmin  :rmin-1+ncells+1)* &
                                                    ds_I(   rmin  :rmin-1+ncells+1)

            ! Note:  nrhs will be 1 in this case, but we can apply this
            ! to all rows at once using nrhs = i*k*n
            nrhs = 1
            ! Loop over all cells in this y-row
            do l_cc=rmin,rmin-1+ncells
              l_I = l_cc  ! Interface ID refers to the Left side of the cell (k-1/2)
              LeftFac    = dt*k_ds_I(l_I  )*sig_I(l_I  )/vol_cc(l_cc)
              RightFac   = dt*k_ds_I(l_I+1)*sig_I(l_I+1)/vol_cc(l_cc)
              CenterFac = LeftFac + RightFac

              if(l_cc.eq.1) then
                !DL_d = nothing     :: No lower diagonal for first row
                D_d(l_cc)  = 1.0_ip + (Imp_fac)*CenterFac & ! This is part of the normal stencil
                                     -(Imp_fac)*LeftFac     ! This line is the BC that ensures
                                                            ! that q(ncell)=q(ncell+1) at t1
                                                            ! i.e. no outward flux

                DU_d(l_cc) =        - (Imp_fac)*RightFac

                  ! RHS contains left boundary term
                BC_left_t0 = q_cc(1)   ! This sets a Neuman condition at t0
                !BC_left_t1 ~= q_cc(1) ! This condition is baked into the matrix
                B_d(l_cc)  =   (1.0_ip-Imp_fac)*LeftFac    * BC_left_t0 + &
                     (1.0_ip - (1.0_ip-Imp_fac)*CenterFac) * q_cc(l_cc) + &
                               (1.0_ip-Imp_fac)*RightFac   * q_cc(l_cc+1)


              elseif(l_cc.lt.ncells)then
                DL_d(l_cc-1) =      - (Imp_fac)*LeftFac
                D_d(l_cc)  = 1.0_ip + (Imp_fac)*CenterFac
                DU_d(l_cc)   =      - (Imp_fac)*RightFac

                B_d(l_cc)    = (1.0_ip-Imp_fac)*LeftFac    * q_cc(l_cc-1) + &
                     (1.0_ip - (1.0_ip-Imp_fac)*CenterFac) * q_cc(l_cc  ) + &
                               (1.0_ip-Imp_fac)*RightFac   * q_cc(l_cc+1)
              elseif(l_cc.eq.ncells)then
                DL_d(l_cc-1) =      - (Imp_fac)*LeftFac
                D_d(l_cc)  = 1.0_ip + (Imp_fac)*CenterFac  & ! This is part of the normal stencil
                                     -(Imp_fac)*RightFac     ! This line is the BC that ensures
                                                             ! that q(ncell)=q(ncell+1) at t1
                                                             ! i.e. no outward flux
                !DU_d = nothing   :: No upper diagonal for last row

                  ! RHS contains right boundary term
                BC_right_t0 = q_cc(ncells)    ! This sets a Neuman condition at t0
                !BC_right_t1 ~= q_cc(ncells)  ! This condition is baked into the matrix
                B_d(l_cc)    = (1.0_ip-Imp_fac)*LeftFac    * q_cc(l_cc-1) + &
                     (1.0_ip - (1.0_ip-Imp_fac)*CenterFac) * q_cc(l_cc  ) + &
                               (1.0_ip-Imp_fac)*RightFac   * BC_right_t0
              endif
            enddo

#ifdef CRANKNIC
      ! Note: The only reason not to use Crank-Nicolson is if you
      !       don't have blas and lapack installed.  This pre-proc.
      !       directive allows this section to be turned off.
            if(useVarDiffH.or.IsLatLon)then
              ! This is the call for solving single or double
              ! precision general tridiagonal Ax=b
              ! Note: This is the function to call if ky is spatially
              ! variable
              if(ip.eq.4)then
                DL_s = real(DL_d,kind=4)
                D_s  = real(D_d ,kind=4)
                DU_s = real(DU_d,kind=4)
                B_s  = real(B_d ,kind=4)
                call sgtsv(     &
                      nlineq, &  !i The order of the matrix A.  N >= 0.
                      nrhs,   &  !i The number of right hand sides
                      DL_s,   &  !b array, dimension (N-1) sub-diagonal
                      D_s,    &  !b dimension (N) diagonal
                      DU_s,   &  !b array, dimension (N-1) super-diagonal
                      B_s,    &  !b dimension (LDB,NRHS) On entry, the N by NRHS matrix of right hand side matrix B.
                      ldb,    &  !i The leading dimension of the array B. LDB >= max(1,N)
                      info)      !o
              elseif(ip.eq.8)then
                call dgtsv(     &
                      nlineq, &  !i The order of the matrix A.  N >= 0.
                      nrhs,   &  !i The number of right hand sides
                      DL_d,   &  !b array, dimension (N-1) sub-diagonal
                      D_d,    &  !b dimension (N) diagonal
                      DU_d,   &  !b array, dimension (N-1) super-diagonal
                      B_d,    &  !b dimension (LDB,NRHS) On entry, the N by NRHS matrix of right hand side matrix B.
                      ldb,    &  !i The leading dimension of the array B. LDB >= max(1,N)
                      info)      !o
              endif
            else
              ! This is the call for solving single or double
              ! precision symmetric positive definite tridiagonal Ax=b
              ! Note: A will only be symmetric if ky is homogeneous
              !       This is really not much faster than dgtsv
              if(ip.eq.4)then
                D_s  = real(D_d ,kind=4)
                DU_s = real(DU_d,kind=4)
                B_s  = real(B_d ,kind=4)
                call sptsv(     &
                      nlineq, &  !i The order of the matrix A.  N >= 0.
                      nrhs,   &  !i The number of right hand sides
                      D_s,    &  !b dimension (N) diagonal
                      DU_s,   &  !b array, dimension (N-1) sub or super-diagonal
                      B_s,    &  !b dimension (LDB,NRHS) On entry, the N by NRHS matrix of right hand side matrix B.
                      ldb,    &  !i The leading dimension of the array B. LDB >= max(1,N)
                      info)      !o
              elseif(ip.eq.8)then
                call dptsv(     &
                      nlineq, &  !i The order of the matrix A.  N >= 0.
                      nrhs,   &  !i The number of right hand sides
                      D_d,    &  !b dimension (N) diagonal
                      DU_d,   &  !b array, dimension (N-1) sub or super-diagonal
                      B_d,    &  !b dimension (LDB,NRHS) On entry, the N by NRHS matrix of right hand side matrix B.
                      ldb,    &  !i The leading dimension of the array B. LDB >= max(1,N)
                      info)      !o
              endif
            endif
#endif
            concen_pd(i,1:nymax,k,n,ts1) = B_d

          enddo ! loop over i
        enddo ! loop over k
      enddo ! loop over n
      if (ip.eq.4)then
        deallocate(DL_s,D_s,DU_s,B_s)
      endif
      deallocate(DL_d,D_d,DU_d,B_d)

      concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts0) = &
        concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts1)

      endif

      end subroutine diffCN_y

!******************************************************************************

      subroutine diffCN_z

      ! Implicit Crank-Nicolson diffusion routine
      ! Implements Eq 4.13 of LeVeque02

      use mesh,          only : &
         IsLatLon,nxmax,nymax,nzmax,nsmax,ts0,ts1,kappa_pd,sigma_nz_pd

      use solution,      only : &
         concen_pd,IsAloft,imin,imax,jmin,jmax,kmin,kmax

      use time_data,     only : &
         dt

      implicit none

      integer :: i,j,n
      integer :: l_I    ! This is the interface index along the particular diffusion direction
      integer :: l_cc   ! This is the cell-centered index along the particular diffusion direction
      integer :: ncells

      real(kind=ip) :: BC_left_t0,BC_left_t1
      real(kind=ip) :: BC_right_t0,BC_right_t1
      real(kind=sp),allocatable,dimension(:) :: DL_s,D_s,DU_s,B_s
      real(kind=ip),allocatable,dimension(:) :: DL_d,D_d,DU_d,B_d
      integer :: nlineq,nrhs,ldb,info
      real(kind=ip) :: LeftFac,CenterFac,RightFac

       ! arrays that live on cell-centers: Note that we only have 1 ghost cell
      real(kind=ip),dimension(0:nzmax+1)     :: q_cc      ! concen
      real(kind=ip),dimension(0:nzmax+1)     :: vol_cc    ! volume of cell

       ! arrays that live on cell interfaces
       !  Note: interface I for cell i is at (i-1/2); i.e. the left or negative
       !  side of i
       !        We only need the interfaces up to the boundary of the domain
       !        (not the ghost cells)
           ! These should be only from 1 to ncells + 1
      real(kind=ip),dimension( 1:nzmax+1)     :: sig_I        ! area of interface
      real(kind=ip),dimension( 1:nzmax+1)     :: ds_I    ! length measure used at interface
      real(kind=ip),dimension( 1:nzmax+1)     :: vavg_I
      real(kind=ip),dimension( 1:nzmax+1)     :: kavg_I
      real(kind=ip),dimension( 1:nzmax+1)     :: k_ds_I  ! k/ds
      !real(kind=ip),dimension( 1:nzmax+1)     :: ksig2_vol_I  ! k*sig*sig/volavg 

      integer :: rmin, rmax     ! min and max indices of the row

#ifdef CRANKNIC
      ! Note: The only reason not to use Crank-Nicolson is if you
      !       don't have blas and lapack installed.  This pre-proc.
      !       directive allows this section to be turned off.
      INTERFACE
        subroutine sgtsv(N,NRHS,DL,D,DU,B,LDB,INFO)
          integer                         ,intent(in)    :: N
          integer                         ,intent(in)    :: NRHS
          real(kind=4),dimension(N-1)     ,intent(inout) :: DL
          real(kind=4),dimension(N)       ,intent(inout) :: D
          real(kind=4),dimension(N-1)     ,intent(inout) :: DU
          integer                         ,intent(in)    :: LDB
          real(kind=4),dimension(ldb,NRHS),intent(inout) :: B
          integer                         ,intent(out)   :: INFO
        end subroutine
        subroutine dgtsv(N,NRHS,DL,D,DU,B,LDB,INFO)
          integer                         ,intent(in)    :: N
          integer                         ,intent(in)    :: NRHS
          real(kind=8),dimension(N-1)     ,intent(inout) :: DL
          real(kind=8),dimension(N)       ,intent(inout) :: D
          real(kind=8),dimension(N-1)     ,intent(inout) :: DU
          integer                         ,intent(in)    :: LDB
          real(kind=8),dimension(ldb,NRHS),intent(inout) :: B
          integer                         ,intent(out)   :: INFO
        end subroutine
        subroutine sptsv(N,NRHS,D,E,B,LDB,INFO)
          integer                         ,intent(in)    :: N
          integer                         ,intent(in)    :: NRHS
          real(kind=4),dimension(N)       ,intent(inout) :: D
          real(kind=4),dimension(N-1)     ,intent(inout) :: E
          integer                         ,intent(in)    :: LDB
          real(kind=4),dimension(ldb,NRHS),intent(inout) :: B
          integer                         ,intent(out)   :: INFO
        end subroutine
        subroutine dptsv(N,NRHS,D,E,B,LDB,INFO)
          integer                         ,intent(in)    :: N
          integer                         ,intent(in)    :: NRHS
          real(kind=8),dimension(N)       ,intent(inout) :: D
          real(kind=8),dimension(N-1)     ,intent(inout) :: E
          integer                         ,intent(in)    :: LDB
          real(kind=8),dimension(ldb,NRHS),intent(inout) :: B
          integer                         ,intent(out)   :: INFO
        end subroutine
      END INTERFACE
#endif

      ! We are diffusing in z so set the length of the cell list accordingly
      !rmin = kmin
      !rmax = kmax
      ! Since sub-grid for diffusion has not been tested, use the full length
      rmin = 1
      rmax = nzmax
      ncells = rmax - rmin + 1

      concen_pd(:,:,:,:,ts1) = 0.0_ip
      ! Neuman boundary conditions in z
        !***  Bottom/Top (Z)
      concen_pd(:,:,     -1,:,ts0) = concen_pd(:,:,    1,:,ts0)
      concen_pd(:,:,      0,:,ts0) = concen_pd(:,:,    1,:,ts0)
      concen_pd(:,:,nzmax+1,:,ts0) = concen_pd(:,:,nzmax,:,ts0)
      concen_pd(:,:,nzmax+2,:,ts0) = concen_pd(:,:,nzmax,:,ts0)

      if(ncells.gt.1)then

      nlineq = ncells
      ldb = nlineq  ! leading dimension of b is num of equations
      if (ip.eq.4)then
        allocate(DL_s(nlineq-1));
        allocate(D_s(nlineq));
        allocate(DU_s(nlineq-1));
        allocate(B_s(nlineq));
      endif
      allocate(DL_d(nlineq-1));
      allocate(D_d(nlineq));
      allocate(DU_d(nlineq-1));
      allocate(B_d(nlineq));

      do n=1,nsmax
        if(.not.IsAloft(n)) cycle

        do j=jmin,jmax
          do i=imin,imax

            ! solve the problem in z for each y
            vol_cc(rmin-1:rmin-1+ncells+1) = kappa_pd(    i,j,rmin-1:rmin-1+ncells+1)
            q_cc(  rmin-1:rmin-1+ncells+1) = concen_pd(   i,j,rmin-1:rmin-1+ncells+1,n,ts0)
            sig_I( rmin  :rmin-1+ncells+1) = sigma_nz_pd( i,j,rmin  :rmin-1+ncells+1)
            vavg_I(rmin  :rmin-1+ncells+1) = 0.5_ip*(vol_cc(  rmin-1:rmin-1+ncells  ) + &
                                                     vol_cc(  rmin  :rmin-1+ncells+1))
            kavg_I(rmin  :rmin-1+ncells+1) = 0.5_ip*(kz(i,j,  rmin-1:rmin-1+ncells  ) + &
                                                     kz(i,j,  rmin  :rmin-1+ncells+1))
            ds_I(  rmin  :rmin-1+ncells+1) = sig_I(           rmin  :rmin-1+ncells+1) / &
                                            vavg_I(           rmin  :rmin-1+ncells+1)
            k_ds_I(rmin  :rmin-1+ncells+1) = kavg_I(          rmin  :rmin-1+ncells+1)* &
                                                    ds_I(     rmin  :rmin-1+ncells+1)

            ! Note:  nrhs will be 1 in this case, but we can apply this
            ! to all rows at once using nrhs = i*j*n
            nrhs = 1
            ! Loop over all cells in this z-row
            do l_cc=1,ncells
              l_I = l_cc  ! Interface ID refers to the Left side of the cell (k-1/2)
              LeftFac    = dt*k_ds_I(l_I  )*sig_I(l_I  )/vol_cc(l_cc)
              RightFac   = dt*k_ds_I(l_I+1)*sig_I(l_I+1)/vol_cc(l_cc)
              CenterFac = LeftFac + RightFac

              D_d(l_cc)  = 1.0_ip + (Imp_fac)*CenterFac
              if(l_cc.eq.1) then
                !DL_d = nothing     :: No lower diagonal for first row
                D_d(l_cc)  = 1.0_ip + (Imp_fac)*CenterFac & ! This is part of the normal stencil
                                     -(Imp_fac)*LeftFac     ! This line is the BC that ensures
                                                            ! that q(ncell)=q(ncell+1) at t1
                                                            ! i.e. no outward flux

                DU_d(l_cc) =        - (Imp_fac)*RightFac

                  ! RHS contains left boundary term
                BC_left_t0 = q_cc(1)   ! This sets a Neuman condition at t0
                !BC_left_t1 ~= q_cc(1) ! This condition is baked into the matrix

                B_d(l_cc)  =        (1.0_ip-Imp_fac)*LeftFac    * BC_left_t0 + &
                          (1.0_ip - (1.0_ip-Imp_fac)*CenterFac) * q_cc(l_cc) + &
                                    (1.0_ip-Imp_fac)*RightFac   * q_cc(l_cc+1)


              elseif(l_cc.lt.ncells)then
                DL_d(l_cc-1) =        - (Imp_fac)*LeftFac
                D_d(l_cc)    = 1.0_ip + (Imp_fac)*CenterFac
                DU_d(l_cc)   =        - (Imp_fac)*RightFac

                B_d(l_cc)    =           (1.0_ip-Imp_fac)*LeftFac    * q_cc(l_cc-1) + &
                               (1.0_ip - (1.0_ip-Imp_fac)*CenterFac) * q_cc(l_cc  ) + &
                                         (1.0_ip-Imp_fac)*RightFac   * q_cc(l_cc+1)
              elseif(l_cc.eq.ncells)then
                DL_d(l_cc-1) =         - (Imp_fac)*LeftFac
                D_d(l_cc)  = 1.0_ip + (Imp_fac)*CenterFac  & ! This is part of the normal stencil
                                     -(Imp_fac)*RightFac     ! This line is the BC that ensures
                                                             ! that q(ncell)=q(ncell+1) at t1
                                                             ! i.e. no outward flux
                !DU_d = nothing   :: No upper diagonal for last row

                  ! RHS contains right boundary term
                BC_right_t0 = q_cc(ncells)    ! This sets a Neuman condition at t0
                !BC_right_t1 ~= q_cc(ncells)  ! This condition is baked into the matrix
                B_d(l_cc)    =           (1.0_ip-Imp_fac)*LeftFac    * q_cc(l_cc-1) + &
                               (1.0_ip - (1.0_ip-Imp_fac)*CenterFac) * q_cc(l_cc  ) + &
                                         (1.0_ip-Imp_fac)*RightFac   * BC_right_t0
              endif
            enddo
#ifdef CRANKNIC
      ! Note: The only reason not to use Crank-Nicolson is if you
      !       don't have blas and lapack installed.  This pre-proc.
      !       directive allows this section to be turned off.
            if(useVarDiffV.or.IsLatLon)then
              ! This is the call for solving single or double
              ! precision general tridiagonal Ax=b
              ! Note: This is the function to call if kz is spatially
              ! variable
              if(ip.eq.4)then
                DL_s = real(DL_d,kind=4)
                D_s  = real(D_d ,kind=4)
                DU_s = real(DU_d,kind=4)
                B_s  = real(B_d ,kind=4)
                call sgtsv(     &
                      nlineq, &  !i The order of the matrix A.  N >= 0.
                      nrhs,   &  !i The number of right hand sides
                      DL_s,   &  !b array, dimension (N-1) sub-diagonal
                      D_s,    &  !b dimension (N) diagonal
                      DU_s,   &  !b array, dimension (N-1) super-diagonal
                      B_s,    &  !b dimension (LDB,NRHS) On entry, the N by NRHS matrix of right hand side matrix B.
                      ldb,    &  !i The leading dimension of the array B. LDB >= max(1,N)
                      info)      !o
              elseif(ip.eq.8)then
                call dgtsv(     &
                      nlineq, &  !i The order of the matrix A.  N >= 0.
                      nrhs,   &  !i The number of right hand sides
                      DL_d,   &  !b array, dimension (N-1) sub-diagonal
                      D_d,    &  !b dimension (N) diagonal
                      DU_d,   &  !b array, dimension (N-1) super-diagonal
                      B_d,    &  !b dimension (LDB,NRHS) On entry, the N by NRHS matrix of right hand side matrix B.
                      ldb,    &  !i The leading dimension of the array B. LDB >= max(1,N)
                      info)      !o
              endif
            else
              ! This is the call for solving single or double
              ! precision symmetric positive definite tridiagonal Ax=b
              ! Note: A will only be symmetric if kz is homogeneous
              !       This is really not much faster than dgtsv
              if(ip.eq.4)then
                D_s  = real(D_d ,kind=4)
                DU_s = real(DU_d,kind=4)
                B_s  = real(B_d ,kind=4)
                call sptsv(     &
                      nlineq, &  !i The order of the matrix A.  N >= 0.
                      nrhs,   &  !i The number of right hand sides
                      D_s,    &  !b dimension (N) diagonal
                      DU_s,   &  !b array, dimension (N-1) sub or super-diagonal
                      B_s,    &  !b dimension (LDB,NRHS) On entry, the N by NRHS matrix of right hand side matrix B.
                      ldb,    &  !i The leading dimension of the array B. LDB >= max(1,N)
                      info)      !o
              elseif(ip.eq.8)then
                call dptsv(     &
                      nlineq, &  !i The order of the matrix A.  N >= 0.
                      nrhs,   &  !i The number of right hand sides
                      D_d,    &  !b dimension (N) diagonal
                      DU_d,   &  !b array, dimension (N-1) sub or super-diagonal
                      B_d,    &  !b dimension (LDB,NRHS) On entry, the N by NRHS matrix of right hand side matrix B.
                      ldb,    &  !i The leading dimension of the array B. LDB >= max(1,N)
                      info)      !o
              endif
            endif
#endif
            concen_pd(i,j,rmin:rmin-1+ncells,n,ts1) = B_d

!            if(i.eq.179.and.j.eq.102)then
!            if(i.eq.5.and.j.eq.5)then
!              do l_cc=rmin,rmin-1+ncells
!                write(*,*)l_cc,q_cc(l_cc),concen_pd(i,j,l_cc,n,ts1)
!              enddo
!              stop 5
!            endif

          enddo ! loop over j
        enddo ! loop over i
      enddo ! loop over n

      if (ip.eq.4)then
        deallocate(DL_s,D_s,DU_s,B_s)
      endif
      deallocate(DL_d,D_d,DU_d,B_d)

      concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts0) = &
        concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts1)

      endif

      end subroutine diffCN_z

      end module Diffusion

