      module Diffusion

      use precis_param

      use io_units

      use global_param,  only : &
         useVarDiffH,useVarDiffV

      implicit none

      real(kind=ip) :: diffusivity_horz    ! horizontal diffusion coefficient (m2/s)
      real(kind=ip) :: diffusivity_vert    ! vertical diffusion coefficient (m2/s)

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
         IsLatLon,nxmax,nymax,nzmax,nsmax,ts0,ts1,kappa_pd,sigma_nx_pd

      use solution,      only : &
         concen_pd,IsAloft

      use time_data,     only : &
         dt

      implicit none

      integer       :: j,k,n  ! These are the indeces mapping to the global arrays
      integer       :: l        ! This is the index along the particular diffusion direction
      integer       :: i_cc
      integer       :: ncells
      integer       :: idx_dum

       ! arrays that live on cell-centers: Note that we have 2 ghost cells
      real(kind=ip),dimension(-1:nxmax+2)               :: update_cc
      real(kind=ip),dimension(-1:nxmax+2)               :: q_cc      ! concen

      real(kind=ip),dimension( 0:nxmax+2)     :: dq_I
      real(kind=ip),dimension( 0:nxmax+2)     :: k_ds2_I
      real(kind=ip) :: ds

      !integer OMP_GET_MAX_THREADS
      !integer OMP_GET_NUM_THREADS
      !integer OMP_GET_THREAD_NUM
      !integer :: nthreads,thread_num
      !logical :: OMP_get_nested

      ! We are diffusing in x so set the length of the cell list accordingly
      ncells = nxmax
      concen_pd(:,:,:,:,ts1) = 0.0_ip
      do n=1,nsmax
        if(.not.IsAloft(n)) cycle
      !!!$OMP PARALLEL DO &
      !!!$OMP DEFAULT(NONE) &
      !!!$OMP SHARED(n,nymax,nzmax,ncells,nsmax,dt,concen_pd,kappa_pd,&
      !!!$OMP sigma_nx_pd,kx,&
      !!!$OMP IsPeriodic),&
      !!!$OMP PRIVATE(l,j,k,q_cc,update_cc,ds,k_ds2_I,&
      !!!$OMP dq_I,i_cc,k_ds2_I)&
      !!!$OMP collapse(2)
        do k=1,nzmax
          do j=1,nymax
            ! Initialize cell-centered values for this x-row
            ! Note: ghost cells should contain q_cc values at edge (Neumann)
            q_cc(-1:ncells+2) = concen_pd(-1:ncells+2,j,k,n,ts0)
            dq_I( 0:ncells+2) = q_cc(0:ncells+2) - q_cc(-1:ncells+1)
              ! Loop over interfaces and get geometry term
            do l=1,ncells+1
                ! ds is the 1/dx for this cell along x
              ds=sigma_nx_pd(l,j,k)/kappa_pd(l,j,k)
              k_ds2_I(l) = 0.5_ip*(kx(l-1,j,k)+kx(l,j,k)) * ds * ds
            enddo
              ! Loop over cells and update
            do i_cc=1,ncells
                ! Eq 4.11 LeVeque02
              update_cc(i_cc) = dt*(k_ds2_I(i_cc+1)*dq_I(i_cc+1) - &
                                    k_ds2_I(i_cc  )*dq_I(i_cc  ))
            enddo ! loop over l (cell centers)

            concen_pd(1:ncells,j,k,n,ts1) = concen_pd(1:ncells,j,k,n,ts0) + &
               update_cc(1:ncells)

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

      ! Explicit diffusion routine. author RP Denlinger

      use mesh,          only : &
         IsLatLon,nxmax,nymax,nzmax,nsmax,ts0,ts1,kappa_pd,sigma_ny_pd

      use solution,      only : &
         concen_pd,IsAloft

      use time_data,     only : &
         dt

      implicit none

      integer       :: i,k,n  ! These are the indeces mapping to the global arrays
      integer       :: l        ! This is the index along the particular diffusion direction
      integer       :: i_cc
      integer       :: ncells
      integer       :: idx_dum

       ! arrays that live on cell-centers: Note that we have 2 ghost cells
      real(kind=ip),dimension(-1:nymax+2)               :: update_cc
      real(kind=ip),dimension(-1:nymax+2)               :: q_cc      ! concen

      real(kind=ip),dimension( 0:nymax+2)     :: dq_I
      real(kind=ip),dimension( 0:nymax+2)     :: k_ds2_I
      real(kind=ip) :: ds

      !integer OMP_GET_MAX_THREADS
      !integer OMP_GET_NUM_THREADS
      !integer OMP_GET_THREAD_NUM
      !integer :: nthreads,thread_num
      !logical :: OMP_get_nested

      ! We are diffusing in x so set the length of the cell list accordingly
      ncells = nymax
      concen_pd(:,:,:,:,ts1) = 0.0_ip
      do n=1,nsmax
        if(.not.IsAloft(n)) cycle
      !!!$OMP PARALLEL DO &
      !!!$OMP DEFAULT(NONE) &
      !!!$OMP SHARED(n,nxmax,nzmax,ncells,nsmax,dt,concen_pd,kappa_pd,&
      !!!$OMP sigma_ny_pd,ky,&
      !!!$OMP IsPeriodic),&
      !!!$OMP PRIVATE(l,j,k,q_cc,update_cc,ds,k_ds2_I,&
      !!!$OMP dq_I,i_cc,k_ds2_I)&
      !!!$OMP collapse(2)
        do k=1,nzmax
          do i=1,nxmax
            ! Initialize cell-centered values for this y-row
            ! Note: ghost cells should contain q_cc values at edge (Neumann)
            q_cc(-1:ncells+2) = concen_pd(i,-1:ncells+2,k,n,ts0)
            dq_I( 0:ncells+2) = q_cc(0:ncells+2) - q_cc(-1:ncells+1)
              ! Loop over interfaces and get geometry term
            do l=1,ncells+1
                ! ds is the 1/dx for this cell along x
              ds=sigma_ny_pd(i,l,k)/kappa_pd(i,l,k)
              k_ds2_I(l) = 0.5_ip*(ky(i,l-1,k)+ky(i,l,k)) * ds * ds
            enddo
              ! Loop over cells and update
            do i_cc=1,ncells
                ! Eq 4.11 LeVeque02
              update_cc(i_cc) = dt*(k_ds2_I(i_cc+1)*dq_I(i_cc+1) - &
                                    k_ds2_I(i_cc  )*dq_I(i_cc  ))
            enddo ! loop over l (cell centers)

            concen_pd(i,1:ncells,k,n,ts1) = concen_pd(i,1:ncells,k,n,ts0) + &
               update_cc(1:ncells)

          enddo ! 
        enddo ! loop over j=1,nymax
     !!!! !$OMP END PARALLEL do

      enddo ! loop over idx_dum

      concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts0) = &
        concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts1)

      end subroutine diff_y

!******************************************************************************

      subroutine diff_z

      ! Explicit diffusion routine. author RP Denlinger

      use mesh,          only : &
         IsLatLon,nxmax,nymax,nzmax,nsmax,ts0,ts1,kappa_pd,sigma_nz_pd

      use solution,      only : &
         concen_pd,IsAloft

      use time_data,     only : &
         dt

      implicit none

      integer       :: i,j,n  ! These are the indeces mapping to the global arrays
      integer       :: l        ! This is the index along the particular diffusion direction
      integer       :: i_cc
      integer       :: ncells
      integer       :: idx_dum

       ! arrays that live on cell-centers: Note that we have 2 ghost cells
      real(kind=ip),dimension(-1:nzmax+2)               :: update_cc
      real(kind=ip),dimension(-1:nzmax+2)               :: q_cc      ! concen

      real(kind=ip),dimension( 0:nzmax+2)     :: dq_I
      real(kind=ip),dimension( 0:nzmax+2)     :: k_ds2_I
      real(kind=ip) :: ds

      !integer OMP_GET_MAX_THREADS
      !integer OMP_GET_NUM_THREADS
      !integer OMP_GET_THREAD_NUM
      !integer :: nthreads,thread_num
      !logical :: OMP_get_nested

      ! We are diffusing in x so set the length of the cell list accordingly
      ncells = nzmax
      concen_pd(:,:,:,:,ts1) = 0.0_ip
      do n=1,nsmax
        if(.not.IsAloft(n)) cycle
      !!!$OMP PARALLEL DO &
      !!!$OMP DEFAULT(NONE) &
      !!!$OMP SHARED(n,nymax,nzmax,ncells,nsmax,dt,concen_pd,kappa_pd,&
      !!!$OMP sigma_nz_pd,kz,&
      !!!$OMP IsPeriodic),&
      !!!$OMP PRIVATE(l,j,k,q_cc,update_cc,ds,k_ds2_I,&
      !!!$OMP dq_I,i_cc,k_ds2_I)&
      !!!$OMP collapse(2)
        do j=1,nymax
          do i=1,nxmax
            ! Initialize cell-centered values for this z-row
            ! Note: ghost cells should contain q_cc values at edge (Neumann)
            q_cc(-1:ncells+2) = concen_pd(i,j,-1:ncells+2,n,ts0)
            dq_I( 0:ncells+2) = q_cc(0:ncells+2) - q_cc(-1:ncells+1)
              ! Loop over interfaces and get geometry term
            do l=1,ncells+1
                ! ds is the 1/dx for this cell along x
              ds=sigma_nz_pd(i,j,l)/kappa_pd(i,j,l)
              k_ds2_I(l) = 0.5_ip*(kz(i,j,l-1)+kz(i,j,l)) * ds * ds
            enddo
              ! Loop over cells and update
            do i_cc=1,ncells
                ! Eq 4.11 LeVeque02
              update_cc(i_cc) = dt*(k_ds2_I(i_cc+1)*dq_I(i_cc+1) - &
                                    k_ds2_I(i_cc  )*dq_I(i_cc  ))
            enddo ! loop over l (cell centers)

            concen_pd(i,j,1:ncells,n,ts1) = concen_pd(i,j,1:ncells,n,ts0) + &
               update_cc(1:ncells)

          enddo ! 
        enddo ! loop over j=1,nymax
     !!!! !$OMP END PARALLEL do

      enddo ! loop over idx_dum



!      integer :: i,j,k,n
!      real(kind=ip) :: k0,k1,k2,km1,km0,dq0,dq1
!      real(kind=ip) :: kap2,sm1,sp1
!
!      concen_pd(:,:,:,:,ts1) = 0.0_ip
!
!      do n=1,nsmax
!        if(.not.IsAloft(n)) cycle
!
!        do j=1,nymax
!          do i=1,nxmax
!            do k=1,nzmax
!              k2 = kz(i,j,k+1)
!              k1 = kz(i,j,k  )
!              k0 = kz(i,j,k-1)
!              km1 = 0.5_ip*(k1+k2)
!              km0 = 0.5_ip*(k1+k0)
!              if (k.eq.1) then
!                ! Homogeneous Neumann conditions at bottom for zero
!                ! diffusive flux
!                dq0 = 0.0_ip
!              else
!                dq0 = concen_pd(i,j,k  ,n,ts0)-concen_pd(i,j,k-1,n,ts0)
!              endif
!              dq1 = concen_pd(i,j,k+1,n,ts0)-concen_pd(i,j,k  ,n,ts0)
!                ! Eq 4.11 LeVeque02
!              if(IsLatLon) then
!                !concen(i,j,k,n,ts1) = concen(i,j,k,n,ts0) + dt_dz2(k)*(dq1*km1-dq0*km0)
!                kap2 = kappa_pd(i,j,k)*kappa_pd(i,j,k)
!                sm1  = sigma_nz_pd(i,j,k-1)
!                sp1  = sigma_nz_pd(i,j,k  )
!                concen_pd(i,j,k,n,ts1) = concen_pd(i,j,k,n,ts0) + dt/kap2 * &
!                                       (sp1*sp1*dq1*km1 - sm1*sm1*dq0*km0)
!              else
!                concen_pd(i,j,k,n,ts1) = concen_pd(i,j,k,n,ts0) + dtodzdz*(dq1*km1-dq0*km0)
!              endif
!            enddo
!          enddo ! loop over i
!        enddo ! loop over j
!
!      enddo ! loop over n

      concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts0) = &
        concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts1)

      end subroutine diff_z

!******************************************************************************

      subroutine diffCN_x

      ! Implicit Crank-Nicolson diffusion routine
      ! Implements Eq 4.13 of LeVeque02

      use mesh,          only : &
         IsLatLon,nxmax,nymax,nzmax,nsmax,ts0,ts1,kappa_pd,sigma_nx_pd

      use solution,      only : &
         concen_pd,IsAloft

      use time_data,     only : &
         dt,dtodxdx

      implicit none

      integer :: i,j,k,n
      real(kind=ip) :: k0,k1,k2,km1,km0
      real(kind=ip) :: km12,r,BC_left,BC_right
        ! It would probably be better to have these on a stack
      real(kind=sp),allocatable,dimension(:) :: DL_s,D_s,DU_s,B_s
      real(kind=dp),allocatable,dimension(:) :: DL_d,D_d,DU_d,B_d
      integer :: nlineq,nrhs,ldb,info
      real(kind=ip) :: sm1,sp1

#ifdef CRANKNIC
      ! Note: The only reason not to use Crank-Nicolson is if you
      !       don't have blas and lapack installed.  This pre-proc.
      !       directive allows this section to be turnes off.
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
      concen_pd(:,:,:,:,ts1) = 0.0_ip

      if(nxmax.gt.1)then

      nlineq = nxmax
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

      r = 0.5_ip*dtodxdx
      do n=1,nsmax
        if(.not.IsAloft(n)) cycle

        do k=1,nzmax
          do j=1,nymax
            ! solve the problem in x for each y
            ! Note:  nrhs will be 1 in this case, but we can apply this
            ! to all rows at once using nrhs = j*k*n
            nrhs = 1
            do i=1,nxmax
              k2 = kx(i+1,j,k)
              k1 = kx(i  ,j,k)
              k0 = kx(i-1,j,k)
              if(IsLatLon)then
                !r = 0.5_ip*dtode2(i,j,k)
                ! for LatLon, r is time/vol2 with the km's containing a surface
                ! area term * diffusivity
                r = 0.5_ip*dt/(kappa_pd(i,j,k)*kappa_pd(i,j,k))
                sm1  = sigma_nx_pd(i-1,j,k)
                sp1  = sigma_nx_pd(i  ,j,k)
                km1 = 0.5_ip*(k1+k2)*sp1*sp1
                km0 = 0.5_ip*(k1+k0)*sm1*sm1
                km12= km0 + km1
              else
                ! for cartesian, r is time/length2 and the km's just are
                ! diffusivities
                km1 = 0.5_ip*(k1+k2)
                km0 = 0.5_ip*(k1+k0)
                km12= km0 + km1
              endif

              D_d(i)  = (1.0_ip + r * km12)
              If(i.eq.1) then
                  ! No lower diagonal for first row
                D_d(i)  = (1.0_ip + r * km12)
                DU_d(i) =         - r * km1
                  ! RHS contains left boundary term
                BC_left = concen_pd(0,j,k,n,ts0)
                B_d(i)  =           r * km0   * 2.0_ip * BC_left    + &
                          (1.0_ip - r * km12) * concen_pd(i  ,j,k,n,ts0) +  &
                                    r * km1   * concen_pd(i+1,j,k,n,ts0)
              elseif(i.lt.nxmax)then
                DL_d(i-1) =         - r * km0
                D_d(i)    = (1.0_ip + r * km12)
                DU_d(i)   =         - r * km1
                B_d(i)    =           r * km0   * concen_pd(i-1,j,k,n,ts0) + &
                            (1.0_ip - r * km12) * concen_pd(i  ,j,k,n,ts0) + &
                                      r * km1   * concen_pd(i+1,j,k,n,ts0)

              elseif(i.eq.nxmax)then
                DL_d(i-1) = -r * km0
                D_d(i)  = (1.0_ip + r*km12)
                  ! No upper diagonal for last row
                  ! RHS contains right boundary term
                BC_right = concen_pd(nxmax+1,j,k,n,ts0)
                B_d(i)    =           r * km0   * concen_pd(i-1,j,k,n,ts0) + &
                            (1.0_ip - r * km12) * concen_pd(i  ,j,k,n,ts0) + &
                                      r * km1   * 2.0_ip * BC_right

              endif
 
            enddo
#ifdef CRANKNIC
      ! Note: The only reason not to use Crank-Nicolson is if you
      !       don't have blas and lapack installed.  This pre-proc.
      !       directive allows this section to be turned off.
            if(useVarDiffH.or.IsLatLon)then
              ! This is the call for solving single or double
              ! precision general tridiagonal Ax=b
              ! Note: This is the function to call if kx is spatially
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
              ! Note: A will only be symmetric if kx is homogeneous
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
         concen_pd,IsAloft

      use time_data,     only : &
         dt,dtodydy

      implicit none

      integer :: i,j,k,n
      real(kind=ip) :: k0,k1,k2,km1,km0
      real(kind=ip) :: km12,r,BC_left,BC_right
      real(kind=sp),allocatable,dimension(:) :: DL_s,D_s,DU_s,B_s
      real(kind=dp),allocatable,dimension(:) :: DL_d,D_d,DU_d,B_d
      integer :: nlineq,nrhs,ldb,info
      real(kind=ip) :: sm1,sp1

#ifdef CRANKNIC
      ! Note: The only reason not to use Crank-Nicolson is if you
      !       don't have blas and lapack installed.  This pre-proc.
      !       directive allows this section to be turnes off.
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
      concen_pd(:,:,:,:,ts1) = 0.0_ip

      if(nymax.gt.1)then

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

      r = 0.5_ip*dtodydy
      do n=1,nsmax
        if(.not.IsAloft(n)) cycle

        do k=1,nzmax
          do i=1,nxmax
            ! solve the problem in y for each x
            ! Note:  nrhs will be 1 in this case, but we can apply this
            ! to all rows at once using nrhs = i*k*n
            nrhs = 1
            do j=1,nymax
              k2 = ky(i,j+1,k)
              k1 = ky(i,j  ,k)
              k0 = ky(i,j-1,k)
              if(IsLatLon)then
                !r = 0.5_ip*dtodn2(i,j,k)
                ! for LatLon, r is time/vol2 with the km's containing a surface
                ! area term * diffusivity
                r = 0.5_ip*dt/(kappa_pd(i,j,k)*kappa_pd(i,j,k))
                sm1  = sigma_ny_pd(i,j-1,k)
                sp1  = sigma_ny_pd(i,j  ,k)
                km1 = 0.5_ip*(k1+k2)*sp1*sp1
                km0 = 0.5_ip*(k1+k0)*sm1*sm1
                km12= km0 + km1
              else
                ! for cartesian, r is time/length2 and the km's just are
                ! diffusivities
                km1 = 0.5_ip*(k1+k2)
                km0 = 0.5_ip*(k1+k0)
                km12= km0 + km1
              endif

              D_d(j)  = (1.0_ip + r * km12)
              If(j.eq.1) then
                  ! No lower diagonal for first row
                D_d(j)  = (1.0_ip + r * km12)
                DU_d(j) =         - r * km1
                  ! RHS contains left boundary term
                BC_left = concen_pd(i,0,k,n,ts0)
                B_d(j) =            r * km0   * 2.0_ip * BC_left    + &
                          (1.0_ip - r * km12) * concen_pd(i,j  ,k,n,ts0) +  &
                                    r * km1   * concen_pd(i,j+1,k,n,ts0)
              elseif(j.lt.nymax)then
                DL_d(j-1) =         - r * km0
                D_d(j)    = (1.0_ip + r * km12)
                DU_d(j)   =         - r * km1
                B_d(j)    =           r * km0   * concen_pd(i,j-1,k,n,ts0) + &
                            (1.0_ip - r * km12) * concen_pd(i,j  ,k,n,ts0) + &
                                      r * km1   * concen_pd(i,j+1,k,n,ts0)

              elseif(j.eq.nymax)then
                DL_d(j-1) = -r * km0
                D_d(j)  = (1.0_ip + r*km12)
                  ! No upper diagonal for last row
                  ! RHS contains right boundary term
                BC_right = concen_pd(i,nymax+1,k,n,ts0)
                B_d(j)    =           r * km0   * concen_pd(i,j-1,k,n,ts0) + &
                            (1.0_ip - r * km12) * concen_pd(i,j  ,k,n,ts0) + &
                                      r * km1   * 2.0_ip * BC_right

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
         concen_pd,IsAloft

      use time_data,     only : &
         dt,dtodzdz

      implicit none

      integer :: i,j,k,n
      real(kind=ip) :: k0,k1,k2,km1,km0
      real(kind=ip) :: km12,r,BC_left,BC_right
      real(kind=sp),allocatable,dimension(:) :: DL_s,D_s,DU_s,B_s
      real(kind=ip),allocatable,dimension(:) :: DL_d,D_d,DU_d,B_d
      integer :: nlineq,nrhs,ldb,info
      real(kind=ip) :: sm1,sp1

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

      concen_pd(:,:,:,:,ts1) = 0.0_ip

      if(nzmax.gt.1)then

      nlineq = nzmax
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

      r = 0.5_ip*dtodzdz
      do n=1,nsmax
        if(.not.IsAloft(n)) cycle

        do i=1,nxmax
          do j=1,nymax
            ! solve the problem in z for each y
            ! Note:  nrhs will be 1 in this case, but we can apply this
            ! to all rows at once using nrhs = i*j*n
            nrhs = 1
            do k=1,nzmax
              k2 = kz(i,j,k+1)
              k1 = kz(i,j,k  )
              k0 = kz(i,j,k-1)
              if(IsLatLon)then
                !r = 0.5_ip*dt_dz2(k)
                ! for LatLon, r is time/vol2 with the km's containing a surface
                ! area term * diffusivity
                r = 0.5_ip*dt/(kappa_pd(i,j,k)*kappa_pd(i,j,k))
                sm1  = sigma_nz_pd(i,j,k-1)
                sp1  = sigma_nz_pd(i,j,k  )
                km1 = 0.5_ip*(k1+k2)*sp1*sp1
                km0 = 0.5_ip*(k1+k0)*sm1*sm1
                km12= km0 + km1
              else
                ! for cartesian, r is time/length2 and the km's just are
                ! diffusivities
                km1 = 0.5_ip*(k1+k2)
                km0 = 0.5_ip*(k1+k0)
                km12= km0 + km1
              endif

              D_d(k)  = (1.0_ip + r * km12)
              If(k.eq.1) then
                  ! No lower diagonal for first row
                D_d(k)  = (1.0_ip + r * km12)
                DU_d(k) =         - r * km1
                  ! RHS contains left boundary term
                BC_left = concen_pd(i,j,k,n,ts0)
                B_d(k)  =           r * km0   * 2.0_ip * BC_left    + &
                          (1.0_ip - r * km12) * concen_pd(i,j,k  ,n,ts0) +  &
                                    r * km1   * concen_pd(i,j,k+1,n,ts0)
              elseif(k.lt.nzmax)then
                DL_d(k-1) =         - r * km0
                D_d(k)    = (1.0_ip + r * km12)
                DU_d(k)   =         - r * km1
                B_d(k)    =           r * km0   * concen_pd(i,j,k-1,n,ts0) + &
                            (1.0_ip - r * km12) * concen_pd(i,j,k  ,n,ts0) + &
                                      r * km1   * concen_pd(i,j,k+1,n,ts0)

              elseif(k.eq.nzmax)then
                DL_d(k-1) = -r * km0
                D_d(k)  = (1.0_ip + r*km12)
                  ! No upper diagonal for last row
                  ! RHS contains right boundary term
                BC_right = concen_pd(i,j,nzmax+1,n,ts0)
                B_d(k)    =           r * km0   * concen_pd(i,j,k-1,n,ts0) + &
                            (1.0_ip - r * km12) * concen_pd(i,j,k  ,n,ts0) + &
                                      r * km1   * 2.0_ip * BC_right

              endif
 
            enddo
#ifdef CRANKNIC
      ! Note: The only reason not to use Crank-Nicolson is if you
      !       don't have blas and lapack installed.  This pre-proc.
      !       directive allows this section to be turnes off.
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
            concen_pd(i,j,1:nzmax,n,ts1) = B_d

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

