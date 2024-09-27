!##############################################################################
!
!  Tephra module
!
!  This module contains all the variables and subroutines needed for managing
!  the tephra species of the concentration array. Note that there may be
!  additional bins that are not managed by this module (i.e. n_gs_max<nsmax).
!
!      subroutine Allocate_Tephra
!      subroutine Allocate_Tephra_Met
!      subroutine Deallocate_Tephra
!      subroutine Deallocate_Tephra_Met
!      subroutine Set_Vf_Meso
!      subroutine Calculate_Tephra_Shape
!      subroutine Sort_Tephra_Size
!      subroutine partition_gsbins
!      subroutine Prune_GS
!      function vset_WH
!      function vset_WH_slip
!      function vset_WH_PCM
!      function vset_Gans
!      function vset_Gans_slip
!      function vset_Stokes_slip
!
!##############################################################################

      module Tephra

      use precis_param

      use io_units

      use global_param,  only : &
         useCalcFallVel,useLogNormGSbins,PI,GRAV

      implicit none

        ! Set everything to private by default
      private

        ! Publicly available subroutines/functions
      public Allocate_Tephra,           &
             Deallocate_Tephra,         &
             Allocate_Tephra_Met,       &
             Deallocate_Tephra_Met,     &
             Calculate_Tephra_Shape,    &
             Sort_Tephra_Size,          &
             Set_Vf_Meso,               &
             Prune_GS,                  &
             vset_WH                        ! This is only public for test cases

        ! Publicly available variables

      real(kind=ip),public :: MagmaDensity   = 2500.0_ip  !density of magma, in kg/m3
      real(kind=ip),public :: DepositDensity = 1000.0_ip  !deposit density, in kg/m3
      real(kind=ip),public :: LAM_GS_THRESH  =  250.0_ip  ! Invokes Cslip once effect is 1%
                                     !=  125.0_ip  ! Invokes Cslip once effect is 2%
                                     !=   50.0_ip  ! Invokes Cslip once effect is 5%
      real(kind=ip),public :: AIRBORNE_THRESH = 1.0e-3_ip ! Mass threshold for flagging bin as empty (kg)

      integer,public :: n_gs_max                      ! # size classes of particles 
      integer,public :: n_gs_aloft                    ! max gs bin still aloft

              ! Fall velocity model
              !   0 = No fall (just tracer)
              !   1 = Wilson and Huang
              !   2 = Wilson and Huang + Cunningham slip
              !   3 = Wilson and Huang + Mod by Pfeiffer Et al.
              !   4 = Ganser
              !   5 = Ganser + Cunningham slip
              !   6 = Stokes flow for spherical particles + slip
      integer,public                                    :: FV_ID
              ! Shape parameter
              !   1 = Wilson and Huang: F = (b+c)/2a maybe also with G = c/b
              !   2 = Sphericity
      integer,public                                    :: Shape_ID

#ifdef USEPOINTERS
      real(kind=ip), dimension(:)  ,pointer,public  :: Tephra_v_s     =>null()    ! Settling vel (m/s)
      real(kind=ip), dimension(:)  ,pointer,public  :: Tephra_gsdiam  =>null()    ! Grain-size diameter 
                                                                                  !  (read in mm, stored in m)
      real(kind=ip), dimension(:)  ,pointer,public  :: Tephra_bin_mass=>null()    ! mass    (kg)
      real(kind=ip), dimension(:)  ,pointer,public  :: Tephra_rho_m   =>null()    ! density (kg/m3)
      real(kind=ip), dimension(:)  ,pointer,public  :: Tephra_gsF     =>null()    ! Grain-size shape (b+c)/2a
      real(kind=ip), dimension(:)  ,pointer,public  :: Tephra_gsG     =>null()    ! Grain-size shape c/b
      real(kind=ip), dimension(:)  ,pointer,public  :: Tephra_gsPhi   =>null()    ! Grain-size shape sphericity
      real(kind=ip), dimension(:,:),pointer         :: Tephra_gsF_fac =>null()    ! Precalculated shape factors
                                                                                  !  i=1 WH Stokes fac
                                                                                  !  i=2 WH Newton fac
                                                                                  !  i=3 Gans Stokes fac
                                                                                  !  i=4 Gans Newton fac
                                                                                  !  i=5 slip adjustment to diameter

      !real(kind=ip),dimension(:,:,:)  ,pointer :: DepositGranularity => null() ! accumulated ash mass on ground
      real(kind=sp),dimension(:,:,:,:),pointer,public :: vf_meso_last_step_MetP_sp => null()
      real(kind=sp),dimension(:,:,:,:),pointer,public :: vf_meso_next_step_MetP_sp => null()
#else
      real(kind=ip), dimension(:)  ,allocatable,public  :: Tephra_v_s         ! Settling vel (m/s)
      real(kind=ip), dimension(:)  ,allocatable,public  :: Tephra_gsdiam      ! Grain-size diameter 
                                                                              !  (read in mm, stored in m)
      real(kind=ip), dimension(:)  ,allocatable,public  :: Tephra_bin_mass    ! mass    (kg)
      real(kind=ip), dimension(:)  ,allocatable,public  :: Tephra_rho_m       ! density (kg/m3)
      real(kind=ip), dimension(:)  ,allocatable,public  :: Tephra_gsF         ! Grain-size shape (b+c)/2a
      real(kind=ip), dimension(:)  ,allocatable,public  :: Tephra_gsG         ! Grain-size shape c/b
      real(kind=ip), dimension(:)  ,allocatable,public  :: Tephra_gsPhi       ! Grain-size shape sphericity
      real(kind=ip), dimension(:,:),allocatable         :: Tephra_gsF_fac     ! Precalculated shape factors
                                                                              !  i=1 WH Stokes fac
                                                                              !  i=2 WH Newton fac
                                                                              !  i=3 Gans Stokes fac
                                                                              !  i=4 Gans Newton fac
                                                                              !  i=5 slip adjustment to diameter
      !real(kind=ip),dimension(:,:,:)  ,allocatable :: DepositGranularity ! accumulated ash mass on ground 
      real(kind=sp),dimension(:,:,:,:),allocatable,public :: vf_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:,:),allocatable,public :: vf_meso_next_step_MetP_sp
#endif

      real(kind=ip),public :: phi_mean
      real(kind=ip),public :: phi_stddev
      real(kind=ip),parameter :: vset_ConvCrit = 0.001_ip

      contains
      !------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Allocate_Tephra
!
!  Called from: Read_Control_File
!  Arguments:
!    none
!
!  This subroutine allocates several variables that describe the grainsize bins.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Allocate_Tephra

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Allocate_Tephra"
      endif;enddo

      allocate(Tephra_v_s(n_gs_max));       Tephra_v_s      = 0.0_ip
      allocate(Tephra_gsdiam(n_gs_max));    Tephra_gsdiam   = 0.0_ip
      allocate(Tephra_bin_mass(n_gs_max));  Tephra_bin_mass = 0.0_ip
      allocate(Tephra_rho_m(n_gs_max));     Tephra_rho_m    = 0.0_ip
      allocate(Tephra_gsF(n_gs_max));       Tephra_gsF      = 0.0_ip
      allocate(Tephra_gsG(n_gs_max));       Tephra_gsG      = 0.0_ip
      allocate(Tephra_gsPhi(n_gs_max));     Tephra_gsPhi    = 0.0_ip
      allocate(Tephra_gsF_fac(n_gs_max,5)); Tephra_gsF_fac  = 0.0_ip

      end subroutine Allocate_Tephra

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Allocate_Tephra_Met
!
!  Called from: Ash3d.F90
!  Arguments:
!    none
!
!  This subroutine allocates the fall-velocity variables on the NWP grid.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Allocate_Tephra_Met

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Allocate_Tephra_Met"
      endif;enddo

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"--------------------------------------------------"
        write(outlog(io),*)"---------- ALLOCATE_TEPHRA_MET -------------------"
        write(outlog(io),*)"--------------------------------------------------"
      endif;enddo

      allocate(vf_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet,n_gs_max))
      allocate(vf_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet,n_gs_max))

      end subroutine Allocate_Tephra_Met

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Deallocate_Tephra
!
!  Called from: dealloc_arrays
!  Arguments:
!    none
!
!  This subroutine deallocates the variables allocated in Allocate_Tephra.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Deallocate_Tephra

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Deallocate_Tephra"
      endif;enddo

#ifdef USEPOINTERS
      if(associated(Tephra_v_s))        deallocate(Tephra_v_s)
      if(associated(Tephra_gsdiam))     deallocate(Tephra_gsdiam)
      if(associated(Tephra_bin_mass))   deallocate(Tephra_bin_mass)
      if(associated(Tephra_rho_m))      deallocate(Tephra_rho_m)
      if(associated(Tephra_gsF))        deallocate(Tephra_gsF)
      if(associated(Tephra_gsG))        deallocate(Tephra_gsG)
      if(associated(Tephra_gsPhi))      deallocate(Tephra_gsPhi)
      if(associated(Tephra_gsF_fac))    deallocate(Tephra_gsF_fac)
#else
      if(allocated(Tephra_v_s))        deallocate(Tephra_v_s)
      if(allocated(Tephra_gsdiam))     deallocate(Tephra_gsdiam)
      if(allocated(Tephra_bin_mass))   deallocate(Tephra_bin_mass)
      if(allocated(Tephra_rho_m))      deallocate(Tephra_rho_m)
      if(allocated(Tephra_gsF))        deallocate(Tephra_gsF)
      if(allocated(Tephra_gsG))        deallocate(Tephra_gsG)
      if(allocated(Tephra_gsPhi))      deallocate(Tephra_gsPhi)
      if(allocated(Tephra_gsF_fac))    deallocate(Tephra_gsF_fac)
#endif

      end subroutine Deallocate_Tephra

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Deallocate_Tephra_Met
!
!  Called from: dealloc_arrays
!  Arguments:
!    none
!
!  This subroutine deallocates the variable allocated in Deallocate_Tephra_Met
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Deallocate_Tephra_Met

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Deallocate_Tephra_Met"
      endif;enddo

#ifdef USEPOINTERS
      if(associated(vf_meso_last_step_MetP_sp)) deallocate(vf_meso_last_step_MetP_sp)
      if(associated(vf_meso_next_step_MetP_sp)) deallocate(vf_meso_next_step_MetP_sp)
#else
      if(allocated(vf_meso_last_step_MetP_sp)) deallocate(vf_meso_last_step_MetP_sp)
      if(allocated(vf_meso_next_step_MetP_sp)) deallocate(vf_meso_next_step_MetP_sp)
#endif

      end subroutine Deallocate_Tephra_Met

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Set_Vf_Meso(Load_MesoSteps)
!
!  Called from: Read_NextMesoStep
!  Arguments:
!    Load_MesoSteps = logical flag indicating that the next NWP step needs to
!                     be calculated.
!
!  This subroutine populates the fall velocity arrays on the MetP grid.  This
!  subroutine is only called if a variable fall velocity model is used.  For
!  each node, the atmospheric properties are loaded, then a function call to
!  one of the v_s models, such as vset_WH.  Finally, these fall velocites on
!  the metP grid are interpolated onto the compH grid at the NWP time step
!  via MR_Regrid_MetP_to_CompH.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Set_Vf_Meso(Load_MesoSteps)

      use solution,      only : &
         IsAloft

      use wind_grid,     only : &
         vf_meso_next_step_sp

      use atmosphere,    only : &
         AirDens_meso_last_step_MetP_sp,AirDens_meso_next_step_MetP_sp,&
         AirVisc_meso_last_step_MetP_sp,AirVisc_meso_next_step_MetP_sp,&
         AirLamb_meso_last_step_MetP_sp,AirLamb_meso_next_step_MetP_sp

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet,MR_dum3d_metP,MR_dum3d_compH,MR_iMetStep_Now,&
           MR_Regrid_MetP_to_CompH

      logical      ,intent(in) :: Load_MesoSteps

      real(kind=sp) :: vf_sp(n_gs_max)
      integer :: i,j,k,is,isize
      real(kind=ip) :: dens,visc,lambda,Kna
      real(kind=ip) :: v_grav_set
      logical,save  :: first_time = .true.

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Set_Vf_Meso"
      endif;enddo

      if(Load_MesoSteps)then
        vf_meso_next_step_sp(:,:,:,:) = 0.0_sp
        if(first_time)then
          is = 1
          first_time = .false.
        else
          is = 2
        endif
        do i=1,nx_submet
         do j=1,ny_submet
            do k=1,np_fullmet
              ! Note: we solve for the fall velocities on the p-levels since
              ! that is where the atmospheric variables live, then we
              ! interpolate onto the compH grid
               if(is.eq.1)then
                 dens   = AirDens_meso_last_step_MetP_sp(i,j,k) ! kg/m^3
                 visc   = AirVisc_meso_last_step_MetP_sp(i,j,k) ! kg/(m s)
                 lambda = AirLamb_meso_last_step_MetP_sp(i,j,k)
               else
                 dens   = AirDens_meso_next_step_MetP_sp(i,j,k) ! kg/m^3
                 visc   = AirVisc_meso_next_step_MetP_sp(i,j,k) ! kg/(m s)
                 lambda = AirLamb_meso_next_step_MetP_sp(i,j,k)
               endif
                 ! Initialize fall velocities too because we will be
                 ! skipping over some slots
               vf_sp(:) = 0.0_ip
               select case (FV_ID)
                 ! Get settling velocity in m/s (will be converted to
                 ! km/hr and negated later)
               case (0)  ! No fall
                 vf_sp(:) = 0.0_ip
               case (1)  ! Wilson and Huang
                 do isize=1,n_gs_max
                   if(.not.IsAloft(isize)) cycle
                   v_grav_set = vset_WH(dens,Tephra_rho_m(isize),visc, &
                                 Tephra_gsdiam(isize),Tephra_gsF_fac(isize,1),&
                                 Tephra_gsF_fac(isize,2))
                   vf_sp(isize) = real(v_grav_set,kind=sp)
                 enddo
               case (2)  ! Wilson and Huang + Cunningham slip
                 do isize=1,n_gs_max
                   if(.not.IsAloft(isize)) cycle
                   Kna = 2.0_ip*lambda/(Tephra_gsdiam(isize)*Tephra_gsF_fac(isize,5))
                   ! The non-continuum effect are > 1% when
                   ! gs<250*lam_col_windp(k)
                   if(Tephra_gsdiam(isize).gt.LAM_GS_THRESH*lambda)then
                     v_grav_set = vset_WH(dens,Tephra_rho_m(isize),visc, &
                                   Tephra_gsdiam(isize),Tephra_gsF_fac(isize,1),Tephra_gsF_fac(isize,2))
                   else
                     v_grav_set = vset_WH_slip(dens,Tephra_rho_m(isize),visc, &
                                   Tephra_gsdiam(isize),Tephra_gsF_fac(isize,1),Tephra_gsF_fac(isize,2), &
                                   Kna)
                   endif
                   vf_sp(isize) = real(v_grav_set,kind=sp)
                 enddo
               case (3)  ! Wilson and Huang + Mod by Pfeiffer Et al.
                 do isize=1,n_gs_max
                   if(.not.IsAloft(isize)) cycle
                   v_grav_set = vset_WH_PCM(dens,Tephra_rho_m(isize),visc, &
                                 Tephra_gsdiam(isize),Tephra_gsF_fac(isize,1),Tephra_gsF_fac(isize,2))
                   vf_sp(isize) = real(v_grav_set,kind=sp)
                 enddo
               case (4)  ! Ganser
                 do isize=1,n_gs_max
                   if(.not.IsAloft(isize)) cycle
                   v_grav_set = vset_Gans(dens,Tephra_rho_m(isize),&
                                          visc,Tephra_gsdiam(isize), &
                                          Tephra_gsF_fac(isize,3),Tephra_gsF_fac(isize,4))
                   vf_sp(isize) = real(v_grav_set,kind=sp)
                 enddo
               case (5)  ! Ganser + Cunningham slip
                 do isize=1,n_gs_max
                   if(.not.IsAloft(isize)) cycle
                   Kna = 2.0_ip*lambda/(Tephra_gsdiam(isize)*Tephra_gsF_fac(isize,5))
                   ! The non-continuum effect are > 1% when
                   ! gs<250*lam_col_windp(k)
                   if(Tephra_gsdiam(isize).gt.LAM_GS_THRESH*lambda)then
                     v_grav_set = vset_Gans(dens,Tephra_rho_m(isize),&
                                            visc,Tephra_gsdiam(isize), &
                                            Tephra_gsF_fac(isize,3),Tephra_gsF_fac(isize,4))
                   else
                     v_grav_set = vset_Gans_slip(dens,Tephra_rho_m(isize),&
                                            visc,Tephra_gsdiam(isize), &
                                            Tephra_gsF_fac(isize,3),Tephra_gsF_fac(isize,4),&
                                            kna)
                   endif
                   vf_sp(isize) = real(v_grav_set,kind=sp)
                 enddo
               case (6)  ! Stokes flow for spherical particles + slip
                 do isize=1,n_gs_max
                   if(.not.IsAloft(isize)) cycle
                   Kna = 2.0_ip*lambda/(Tephra_gsdiam(isize))
                   v_grav_set = vset_Stokes_slip(Tephra_rho_m(isize),visc,Tephra_gsdiam(isize),Kna)
                   vf_sp(isize) = real(v_grav_set,kind=sp)
                 enddo
               case default ! Wilson and Huang
                 do isize=1,n_gs_max
                   if(.not.IsAloft(isize)) cycle
                   v_grav_set = vset_WH(dens,Tephra_rho_m(isize),visc, &
                                 Tephra_gsdiam(isize),Tephra_gsF_fac(isize,1),Tephra_gsF_fac(isize,2))
                   vf_sp(isize) = real(v_grav_set,kind=sp)
                 enddo
               end select
               vf_meso_next_step_MetP_sp(i,j,k,:) = vf_sp(1:n_gs_max)
            enddo !k
          enddo !j
        enddo !i

        ! Now need to interpolate vs_meso_last_step_MetP_sp onto vf_meso_1_sp
        do isize=1,n_gs_max
          if(.not.IsAloft(isize)) cycle
          MR_dum3d_metP(:,:,:) = vf_meso_next_step_MetP_sp(:,:,:,isize)
          call MR_Regrid_MetP_to_CompH(MR_iMetStep_Now+1)
          vf_meso_next_step_sp(:,:,:,isize) = MR_dum3d_compH(:,:,:)*(-1.0_sp)
        enddo
      endif

      end subroutine Set_Vf_Meso

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Calculate_Tephra_Shape
!
!  Called from: Read_Control_File
!  Arguments:
!    none
!
!  This subroutine calculates various aspect of the tephra shape needed for
!  the Wilson/Huang, Ganser or Cunningham slip models.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Calculate_Tephra_Shape

      integer :: isize,j
      !real(kind=ip) :: ellipse_area,ellipse_vol,equiv_rad
      !real(kind=ip) :: diamN ! diamter of sphere with equivalent projected area
      !real(kind=ip) :: diamS ! diamter of sphere with equivalent surface area
      !real(kind=ip) :: diamV ! diamter of sphere with equivalent volume
      real(kind=ip) :: p_exp
      !real(kind=ip) :: sphere_area
!      real(kind=ip) :: phi_sphere

      real(kind=ip) :: Dahneke_LD(13), Dahneke_RL(13),onF
      real(kind=ip) :: tmp_b, tmp_c
      integer       :: Dahni

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Calculate_Tephra_Shape"
      endif;enddo

!     Calculating terms for Cunningham slip coefficients for non-spherical
!     particles
      Dahneke_LD( 1) =   1.0_ip
      Dahneke_LD( 2) =   2.0_ip
      Dahneke_LD( 3) =   3.0_ip
      Dahneke_LD( 4) =   4.0_ip
      Dahneke_LD( 5) =   6.0_ip
      Dahneke_LD( 6) =   8.0_ip
      Dahneke_LD( 7) =  10.0_ip
      Dahneke_LD( 8) =  15.0_ip
      Dahneke_LD( 9) =  20.0_ip
      Dahneke_LD(10) =  30.0_ip
      Dahneke_LD(11) =  50.0_ip
      Dahneke_LD(12) =  75.0_ip
      Dahneke_LD(13) = 100.0_ip
      Dahneke_RL( 1) = 1.00_ip
      Dahneke_RL( 2) = 1.23_ip
      Dahneke_RL( 3) = 1.39_ip
      Dahneke_RL( 4) = 1.51_ip
      Dahneke_RL( 5) = 1.72_ip
      Dahneke_RL( 6) = 1.87_ip
      Dahneke_RL( 7) = 2.00_ip
      Dahneke_RL( 8) = 2.25_ip
      Dahneke_RL( 9) = 2.43_ip
      Dahneke_RL(10) = 2.69_ip
      Dahneke_RL(11) = 3.02_ip
      Dahneke_RL(12) = 3.28_ip
      Dahneke_RL(13) = 3.47_ip

      do isize=1,n_gs_max
        ! Precalculate the WilsonHuang factors based on particle shapes
        Tephra_gsF_fac(isize,1) = Tephra_gsF(isize)**(-0.828_ip)   ! WH Stokes factor
        Tephra_gsF_fac(isize,2) = sqrt(1.07_ip-Tephra_gsF(isize))  ! WH Newton factor

        ! and precalculate the K1 and K2 if we use the Ganser model
        if(Shape_ID.eq.1)then
          ! First, if we are using F and G for shape, get sphericity.
          ! Note: sphericity is the ratio of the surface area of a sphere with
          !       equivalent volume to the actual surface area of the particle
          ! following http://en.wikipedia.org/wiki/Ellipsoid
          p_exp = 1.6075_ip
          tmp_b = 2.0_ip*Tephra_gsF(isize)/(1+Tephra_gsG(isize))
          tmp_c = 2.0_ip*Tephra_gsF(isize)*Tephra_gsG(isize)/(1+Tephra_gsG(isize))

          Tephra_gsPhi(isize) = (tmp_b*tmp_c)**(2.0_ip/3.0_ip) * &
                        ((tmp_b**p_exp + &
                          tmp_c**p_exp + &
                         (tmp_b*tmp_c)**p_exp)/3.0_ip)**(-1.0_ip/p_exp)
        elseif(Shape_ID.eq.2)then
          ! Shape factor is given as sphericity
          ! Assume prolate ellipsoids with B=C (i.e. G=1.0)
          ! Need to calculate F for logging
          Tephra_gsG(isize) = 1.0_ip
          ! Looks like the solution if 3*phi^(-p) = F^(-p/3)*(2+F)
          ! Probably not worth calculating here.
          Tephra_gsF(isize) = 1.0_ip
        endif

        !ellipse_vol  = 4.0_ip*PI*Tephra_gsF(i)*Tephra_gsF(i)/3.0_ip
        !equiv_rad    = (ellipse_vol*3.0_ip/(4.0_ip*PI))**(1.0_ip/3.0_ip)
        !sphere_area  = 4.0_ip*PI*equiv_rad*equiv_rad
        !ellipse_area = 4.0_ip*PI*(((Tephra_gsF(i)*Tephra_gsF(i))**p_exp + &
        !                    2.0_ip*(Tephra_gsF(i))**p_exp)/3.0_ip)**(1.0_ip/p_exp)
        !phi_sphere = sphere_area/ellipse_area

          ! Calculate K1 for Ganser model (Table 7 Isometric)
        Tephra_gsF_fac(isize,3) = 3.0_ip/(1.0_ip+2.0_ip*Tephra_gsPhi(isize)**(-0.5_ip))
          ! Calculate K2 for Ganser model (Table 7 Isometric)
        Tephra_gsF_fac(isize,4) = 1.84148_ip*(-log10(Tephra_gsPhi(isize)))**0.5743_ip
        Tephra_gsF_fac(isize,4) = 10.0_ip**Tephra_gsF_fac(isize,4);
          ! Now calculate the aspherical correction to the Cunningham
          ! slip correction using Table 1 from Dahneke, Aerosol Sci, v4p163
          ! 1973.
        onF = 1.0_ip/Tephra_gsF(isize)
        Dahni = 0
        do j = 1,12
          if(onF.ge.Dahneke_LD(j).and.onF.lt.Dahneke_LD(j+1)) Dahni = j
        enddo
        if(Dahni.eq.0)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"WARNING: Could not determine slip correction factor."
            write(outlog(io),*)"         Setting factor to 1.0"
          endif;enddo
          Tephra_gsF_fac(isize,5) = 1.0_ip
        else
          Tephra_gsF_fac(isize,5) = Dahneke_RL(Dahni) + (onF-Dahneke_LD(Dahni))* &
                         (Dahneke_RL(Dahni+1)-Dahneke_RL(Dahni))/     &
                          (Dahneke_LD(Dahni+1)-Dahneke_LD(Dahni))
        endif
      enddo

      end subroutine Calculate_Tephra_Shape

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Sort_Tephra_Size
!
!  Called from: Read_Control_File
!  Arguments:
!    none
!
!  This subroutine sorts the list of grainsizes as read by the control file
!  according to size with the smallest in the first slot.  This is a necessary
!  step when calculating a log-normal grainsize distribution.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Sort_Tephra_Size

      integer       :: i
      integer       :: isize
      real(kind=ip) :: tmp2,tmp3,tmp4,tmp5,tmp6,tmp7
      real(kind=ip) :: temp_a(5)

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Sort_Tephra_Size"
      endif;enddo

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*) 'WARNING: Sorting grain-size bins by size'
      endif;enddo

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*) 'GSD before sorting:'
      endif;enddo
      do isize=1,n_gs_max
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)&
             'isize = ',isize,', gsdiam = ',real(Tephra_gsdiam(isize),kind=sp),&
             ', rho_m = ',real(Tephra_rho_m(isize),kind=sp)
        endif;enddo
      enddo

        ! Using insertion sort on p321 of Numerical Recipes
      do isize=2,n_gs_max
        tmp2   = Tephra_gsdiam(isize)
        tmp3   = Tephra_bin_mass(isize)
        tmp4   = Tephra_rho_m(isize)
        tmp5   = Tephra_v_s(isize)
        tmp6   = Tephra_gsF(isize)
        tmp7   = Tephra_gsG(isize)
        temp_a = Tephra_gsF_fac(isize,:)
        do i=isize-1,1,-1
            ! sort on grain-size
          if (Tephra_gsdiam(i).le.tmp2) goto 101
          Tephra_gsdiam(i+1)    = Tephra_gsdiam(i)
          Tephra_bin_mass(i+1)  = Tephra_bin_mass(i)
          Tephra_rho_m(i+1)     = Tephra_rho_m(i)
          Tephra_v_s(i+1)       = Tephra_v_s(i)
          Tephra_gsF(i+1)       = Tephra_gsF(i)
          Tephra_gsG(i+1)       = Tephra_gsG(i)
          Tephra_gsF_fac(i+1,:) = Tephra_gsF_fac(i,:)
        enddo
        i=0
 101    Tephra_gsdiam(i+1)   = tmp2
        Tephra_bin_mass(i+1) = tmp3
        Tephra_rho_m(i+1)    = tmp4
        Tephra_v_s(i+1)      = tmp5
        Tephra_gsF(i+1)      = tmp6
        Tephra_gsG(i+1)      = tmp7
        Tephra_gsF_fac(i+1,:)= temp_a
      enddo

      if (useLogNormGSbins) call partition_gsbins(phi_mean,phi_stddev)

      end subroutine Sort_Tephra_Size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  partition_gsbins(mu,sigma)
!
!  Called from: Sort_Tephra_Size
!  Arguments:
!    mu    = grainsize distribution average phi
!    sigma = grainsize distribution standard deviation in phi
!
!  This subroutine will create a log-normal distribution of grainsizes (normal
!  in the variable phi) according to the input parameters for the mean and
!  standard distribution in phi-space (mu, sigma).  The phi-value of bin
!  boudaries is assumed to be the average of neighboring phi values, so this
!  works best if bins are equally spaced in phi and that the bins comfortably
!  span the distribution.  Values for the upper and lower bins are calculated
!  by integration the tails to the far boundary phi-value.  This log-normal
!  distribution will only be applied to the mass-fraction not already specified.
!  If no mass is specified in bins, then the total distribution will be this
!  log-normal distribution.  If some bin have mass already specified, this
!  log-normal distribution will supplement.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine partition_gsbins(mu,sigma)

      integer :: isize

      real(kind=ip) :: mu,sigma

      real(kind=ip),dimension(n_gs_max)   :: phi
      real(kind=ip),dimension(n_gs_max-1) :: phi_boundaries
      real(kind=ip),dimension(n_gs_max)   :: suppl_frac
      integer       :: mid_bin
      real(kind=ip) :: mid_bin_neg_half
      real(kind=ip) :: mid_bin_pos_half
      real(kind=ip) :: erf_at_a,erf_at_b
      real(kind=ip) :: fac1,frac_to_distrib

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine partition_gsbins"
      endif;enddo

      suppl_frac = 0.0_ip
      mid_bin = 1
      fac1 = 0.5_ip

      !Note: This subroutine assumes that the grainsmax bins have been
      !      sorted from smallest to largest

      ! We have n_gs_max bins with real grainsmax info.  The last bin read
      ! in (init_n_gs_max = n_gs_max+1) is
      ! the remaining fraction of the total volume that will be
      ! redistributed among the first n_gs_max bins with a normal
      ! distribution in phi.

      if(sigma.le.0.0_ip)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: StdDev <= 0.0 for supplimental GS"
          write(errlog(io),*)"       distribution"
        endif;enddo
        stop 1
      endif

      ! Get the phi of the other bins
      do isize = 1,n_gs_max
        phi(isize) = -log(Tephra_gsdiam(isize))/log(2.0_ip)
      enddo

      ! boundary between bins is the average phi
      do isize = 1,n_gs_max-1
        phi_boundaries(isize) = (phi(isize)+phi(isize+1))*0.5_ip
      enddo

      if (n_gs_max.eq.0)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Must have n_gs_max.ge.1"
        endif;enddo
        stop 1
      else if (n_gs_max.eq.1)then
        suppl_frac(1)=1.0_ip
      else
        ! Convert boundaries between grainsmaxs given to locations on a
        ! standard normal Gaussian, given mu and sigma
        do isize = 1,n_gs_max-1
          phi_boundaries(isize) = (phi_boundaries(isize) - mu)/sigma
        enddo
        if(phi_boundaries(1).lt.0.0_ip)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: Mean is not within GS-distribution."
            write(errlog(io),*)&
               "        Bin    :       gsdiam        :",&
               "        phi    :       right bin boundary (phi)"
            do isize = 1,n_gs_max
              write(errlog(io),*)isize,real(Tephra_gsdiam(isize),kind=sp),&
                          real(phi(isize),kind=sp),&
                          real(phi_boundaries(isize),kind=sp)
            enddo
            write(errlog(io),*)"Mean = ",mu
          endif;enddo
          stop 1
        endif
        if(n_gs_max.eq.2)then
          ! We're dividing the surplus into two bins
           erf_at_a = 0.5_ip - fac1*erf(abs(phi_boundaries(1)))
           if(phi_boundaries(1).lt.0.0_ip)then
             suppl_frac(1) = erf_at_a
             suppl_frac(2) = 1.0_ip - erf_at_a
           else
             suppl_frac(2) = erf_at_a
             suppl_frac(1) = 1.0_ip - erf_at_a
           endif
        else if(n_gs_max.ge.3)then
          do isize = 2,n_gs_max-1
            ! Find the bin that contains the mean of the distribution
            if(phi_boundaries(isize).le.0.0_ip.and.phi_boundaries(isize-1).ge.0.0_ip)then
              mid_bin = isize
            endif
          enddo
        endif

        ! Integrating the positive half (lower index) of phi
        !   integrate from zero to the first boundary
        mid_bin_pos_half = fac1*erf(abs(phi_boundaries(mid_bin-1)))
        suppl_frac(mid_bin) = mid_bin_pos_half
        ! Now integrate the remaining bins on the positive side of mu
        do isize = mid_bin-1,2,-1
          erf_at_a = fac1*erf(abs(phi_boundaries(isize)))
          erf_at_b = fac1*erf(abs(phi_boundaries(isize-1)))
          suppl_frac(isize) = erf_at_b - erf_at_a
        enddo
        ! The last positive bin is from phi_boundaries(1) to Inf
        suppl_frac(1) = 0.5_ip - fac1*erf(abs(phi_boundaries(1)))

        ! Now integrate the negative side
        mid_bin_neg_half = fac1*erf(abs(phi_boundaries(mid_bin)))
        suppl_frac(mid_bin) = mid_bin_neg_half
        ! Now integrate the remaining bins on the negative side of mu
        do isize = mid_bin+1,n_gs_max-1
          erf_at_a = fac1*erf(abs(phi_boundaries(isize-1)))
          erf_at_b = fac1*erf(abs(phi_boundaries(isize)))
          suppl_frac(isize) = erf_at_b - erf_at_a
        enddo
          ! The last negative bin is from phi_boundaries(1) to -Inf
        suppl_frac(n_gs_max) = 0.5_ip - fac1*erf(abs(phi_boundaries(n_gs_max-1)))
          ! Reconstruct middle bin from the two one-sided integrations
        suppl_frac(mid_bin) = mid_bin_pos_half + mid_bin_neg_half
      endif

      frac_to_distrib = 1.0_ip - sum(Tephra_bin_mass(1:n_gs_max))
      do isize = 1,n_gs_max
        Tephra_bin_mass(isize) = Tephra_bin_mass(isize) + suppl_frac(isize)*frac_to_distrib
      enddo

      end subroutine partition_gsbins

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Prune_GS
!
!  Called from: Ash3d.F90
!  Arguments:
!    none
!
!  This subroutine the total mass aloft in each grainsize bin and compare it
!  to the value AIRBORNE_THRESH.  If the mass drops below that value, that
!  bin is flagged as fully deposited (or otherwise flushed out) and removed
!  from further calculations.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Prune_GS

      use solution,      only : &
           mass_aloft,IsAloft

      integer :: isize

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Prune_GS"
      endif;enddo

      ! Loop through all the grain sizes and check if any have completely
      ! flushed out.  If so, then 
      ! modify the max index of the grain-size loop (n_gs_aloft)
      ! The threshold for collapsing GS array is if there is less than a gram in
      ! that size bin aloft
      ! Note: mass_aloft is calculated for all species in Output_Vars:Calc_AshVol_Aloft
      n_gs_aloft = 0
      do isize = 1,n_gs_max
        if(IsAloft(isize).and. &                     ! if bin is currently flagged as aloft
           mass_aloft(isize).lt.AIRBORNE_THRESH)then ! but the mass is less than the thresh
          IsAloft(isize) = .false.
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Grainsize bin ",isize," has fully deposited or left the domain."
          endif;enddo
        else
          n_gs_aloft = n_gs_aloft + 1
        endif
      enddo

      end subroutine Prune_GS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  vset_WH(rho_air,rho_m,eta,diam,Ffac1,Ffac2)
!
!  Called from: Set_Vf_Meso
!  Arguments:
!    rho_air= density of air (kg/m^3)
!    rho_m  = density of particle (kg/m^3)
!    eta    = viscosity of air (kg/(m s))
!    diam   = diameter of particle (m)
!    Ffac1  = WH precalculated factor 1 = F**(-0.828)
!    Ffac2  = WH precalculated factor 2 = sqrt(1.07-F)
!
!  This function returns the terminal fall velocity of a prolate ellipsoidal particle
!  given the density and viscosity of air along with the density, diameter and
!  shape of the particle.  Velocity is returned in m/s.  The drag coefficient
!  is that of Wilson and Huang which results in a quadratic expression that
!  allows direct calculation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function vset_WH(rho_air,rho_m,eta,diam,Ffac1,Ffac2)

      real(kind=ip) :: vset_WH   ! Settling velocity in m/s
      real(kind=ip) :: rho_air   ! density of air in kg/m3
      real(kind=ip) :: rho_m     ! density of the particle in km/m3
      real(kind=ip) :: eta       ! dynamic viscosity of air in (kg/(m s))
      real(kind=ip) :: diam      ! diameter of the particle in m
      real(kind=ip) :: Ffac1     ! = F**(-0.828)
      real(kind=ip) :: Ffac2     ! = sqrt(1.07-F)

      real(kind=ip) :: vset2

      ! Wilson and Huang (1979)
      !  EPSL v. 44, pp. 311-324: https://doi.org/10.1016/0012-821X(79)90179-1
      ! This model gives the drag coefficient
      ! of tephra as a function of its shape by defining tephra fragments as 
      ! ellipsoidal with axes a, b, and c. they also define F=(b+c)/2a, and give the
      ! drag coefficient as:
      !   Cd = (24/Re)*F**(-0.828) + 2*sqrt(1.07-F)
      ! They have a table in which a, b, and c were measured for pumice, glass, and
      ! feldspar fragments. The average values of F for these three fragment types
      ! ranges from 0.39 to 0.42.

      ! F = 0.4
      ! Ffac1 = 2.13547419198160     ! = F**(-0.828)  where F=0.4
      ! Ffac1 = 1.34072838100802     ! = F**(-0.32)  where F=0.4; Suzuki's mod
      ! Ffac2 = 0.818535277187245    ! = sqrt(1.07-F) where F=0.4

      ! F = 0.6666
      ! Ffac1 = 1.39896599613964      ! = F**(-0.828)  where F=0.6666
      ! Ffac1 = 1.13854602830969      ! = F**(-0.32)  where F=0.6666; Suzuki's mod
      ! Ffac2 = 0.635137780328017     ! = sqrt(1.07-F) where F=0.6666

      ! The following roots are solved by the maxima script:
      !  load(f90)$
      !  eq_1 : Vs*Vs = (4*diam*rho_m*g/(3*Cd*rho_air))$
      !  eq_2 : Cd = (24/Re)*Ffac1 + 2*Ffac2$
      !  eq_3 : Re = (rho_air*Vs*diam)/eta$
      !  eq_4:subst(eq_3,eq_2)$
      !  eq_5:subst(%,eq_1)$
      !  VsSoln : solve(%,Vs);
      !  f90(VsSoln[2]);    /* We want the second (positive) root */

      !vset1 =-(sqrt(6.)*sqrt(diam**3.*Ffac2*g*rho_air*rho_m + &
      !         54.*eta**2.*Ffac1**2.)+18.*eta*Ffac1)/(diam*Ffac2*rho_air)/3.0
      !vset2 = (sqrt(6.)*sqrt(diam**3.*Ffac2*g*rho_air*rho_m + &
      !         54.*eta**2.*Ffac1**2.)-18.*eta*Ffac1)/(diam*Ffac2*rho_air)/3.0

      do io=1,2;if(VB(io).le.verbosity_debug2)then
        write(outlog(io),*)"     Entered function vset_WH"
      endif;enddo

      ! Here's expanding out exponents and precalculating roots
      vset2 = (2.44948974278318_ip*sqrt(diam*diam*diam*Ffac2*GRAV*rho_air*rho_m + &
               54.0_ip*eta*eta*Ffac1*Ffac1)-18.0_ip*eta*Ffac1)/(diam*Ffac2*rho_air)/3.0_ip
      vset_WH = vset2

      end function vset_WH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  vset_WH_slip(rho_air,rho_m,eta,diam,Ffac1,Ffac2,Kna)
!
!  Called from: Set_Vf_Meso
!  Arguments:
!    rho_air= density of air (kg/m^3)
!    rho_m  = density of particle (kg/m^3)
!    eta    = viscosity of air (kg/(m s))
!    diam   = diameter of particle (m)
!    Ffac1  = WH precalculated factor 1 = F**(-0.828)
!    Ffac2  = WH precalculated factor 2 = sqrt(1.07-F)
!    Kna    = Knudsen number (dimensionless)
!
!  This function returns the terminal fall velocity of a prolate ellipsoidal particle
!  given the density and viscosity of air along with the density, diameter and
!  shape of the particle and Knudsen number.  The Knudsen number is important
!  at high altitudes with low air density and a large mean free path.  In these
!  environments, Cunningham slip corrections should be applied to the fall model.
!  Velocity is returned in m/s.  The drag coefficient is that of Wilson and Huang
!  but with the slip correction.  The result is no longer quadratic and must
!  be solved iteritively with Cd and Vs updated until Vs does not change by more than
!  a small fraction vset_ConvCrit (~0.001).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function vset_WH_slip(rho_air,rho_m,eta,diam,Ffac1,Ffac2,Kna)
      ! Modification to Wilson and Huang's with Cunningham slip

      real(kind=ip) :: vset_WH_slip      ! Settling velocity in m/s
      real(kind=ip) :: rho_air   ! density of air in kg/m3
      real(kind=ip) :: rho_m     ! density of the particle in km/m3
      real(kind=ip) :: eta       ! dynamic viscosity of air in (kg/(m s))
      real(kind=ip) :: diam      ! diameter of the particle in m
      real(kind=ip) :: Ffac1     ! = F**(-0.828)
      real(kind=ip) :: Ffac2     ! = sqrt(1.07-F)
      real(kind=ip) :: Kna       ! adjusted Knudsen number 

      real(kind=ip) :: vnew, vold, Re            ! old and new settling velocity
      real(kind=ip) :: Cd                        ! drag coefficient
      real(kind=ip) :: Cslip                     ! Slip correction

      do io=1,2;if(VB(io).le.verbosity_debug2)then
        write(outlog(io),*)"     Entered function vset_WH_slip"
      endif;enddo

      vold = 1.0_ip                              ! assume an initial settling velocity of 1 m/s

      Re = rho_air*vold*diam/eta                 ! Eq. 4  of Wilson79

      ! Eq. 9.34 of Seinfeld and Pandis
      Cslip = 1.0_ip+ Kna * (1.257_ip + 0.4_ip*exp(-1.1_ip/Kna))

      ! Get initial Cd from Wilson79, Eq 12
      Cd = (24.0_ip/Re)*Ffac1 + Ffac2

      vnew = sqrt((4.0_ip*rho_m*diam*GRAV)/(3.0_ip*rho_air*Cd))   ! Eq. 15 of WoodsBursik91

      do while ((abs(vold-vnew)/vnew).gt.vset_ConvCrit)
        vold = vnew
        Re = rho_air*vold*diam/eta
        Cd = (24.0_ip/Re)*Ffac1 + Ffac2
        Cd = Cd/Cslip
        vnew = sqrt((4.0_ip*rho_m*diam*GRAV)/(3.0_ip*rho_air*Cd))
      enddo

      vset_WH_slip = vnew

      return

      end function vset_WH_slip

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  vset_WH_PCM(rho_air,rho_m,eta,diam,Ffac1,Ffac2)
!
!  Called from: Set_Vf_Meso
!  Arguments:
!    rho_air= density of air (kg/m^3)
!    rho_m  = density of particle (kg/m^3)
!    eta    = viscosity of air (kg/(m s))
!    diam   = diameter of particle (m)
!    Ffac1  = WH precalculated factor 1 = F**(-0.828)
!    Ffac2  = WH precalculated factor 2 = sqrt(1.07-F)
!
!  This function returns the terminal fall velocity of a prolate ellipsoidal particle
!  given the density and viscosity of air along with the density, diameter and
!  shape of the particle.  Velocity is returned in m/s.  The drag coefficient
!  is that of Wilson and Huang with the modification suggested by Pfeiffer, Costa,
!  and Macedonio.  The result is no longer quadratic and must be solved iteritively
!  with Cd and Vs updated until Vs does not change by more than a small fraction
!  vset_ConvCrit (~0.001).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function vset_WH_PCM(rho_air,rho_m,eta,diam,Ffac1,Ffac2)
      ! Modification to Wilson and Huang's model suggested by:
      ! T. Pfeiffer and A. Costa and G. Macedonio, JVGR, v140n4p273 2005
      ! DOI:10.1016/j.jvolgeores.2004.09.001

      real(kind=ip) :: vset_WH_PCM      ! Settling velocity in m/s
      real(kind=ip) :: rho_air   ! density of air in kg/m3
      real(kind=ip) :: rho_m     ! density of the particle in km/m3
      real(kind=ip) :: eta       ! dynamic viscosity of air in (kg/(m s))
      real(kind=ip) :: diam      ! diameter of the particle in m
      real(kind=ip) :: Ffac1     ! = F**(-0.828)
      real(kind=ip) :: Ffac2     ! = sqrt(1.07-F)

      real(kind=ip) :: vnew, vold, Re            ! old and new settling velocity
      real(kind=ip) :: Cd                        ! drag coefficient

      real(kind=ip) :: Cd100

      do io=1,2;if(VB(io).le.verbosity_debug2)then
        write(outlog(io),*)"     Entered function vset_WH_PCM"
      endif;enddo

      Cd100 = (24.0_ip/100.0_ip)*Ffac1 + Ffac2

      vold = 1.0_ip                              ! assume an initial settling velocity of 1 m/s

      Re = rho_air*vold*diam/eta                 ! Eq. 4  of Wilson79

      ! Get initial Cd from Wilson79, Eq 12
      Cd = (24.0_ip/Re)*Ffac1 + Ffac2
      vnew = sqrt((4.0_ip*rho_m*diam*GRAV)/(3.0_ip*rho_air*Cd))   ! Eq. 15 of WoodsBursik91

      do while ((abs(vold-vnew)/vnew).gt.vset_ConvCrit)
        vold = vnew
        Re = rho_air*vold*diam/eta
        if (Re.lt.100.0_ip)then
          Cd = (24.0_ip/Re)*Ffac1 + Ffac2
        elseif (Re.gt.1000.0_ip)then
          Cd = 1.0_ip
        else
          ! interpolate between values at Re=1000 and 100
          Cd = 1.0_ip-(1.0_ip-Cd100)/900.0_ip *(1000.0_ip-Re)
        endif
        vnew = sqrt((4.0_ip*rho_m*diam*GRAV)/(3.0_ip*rho_air*Cd))
      enddo

      vset_WH_PCM = vnew

      return

      end function vset_WH_PCM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  vset_Gans(rho_air,rho_m,eta,diam,K1,K2)
!
!  Called from: Set_Vf_Meso
!  Arguments:
!    rho_air= density of air (kg/m^3)
!    rho_m  = density of particle (kg/m^3)
!    eta    = viscosity of air (kg/(m s))
!    diam   = diameter of particle (m)
!    K1     = precalculated Stokes shape factor
!    K2     = precalculated Newtons shape factor
!
!  This function returns the terminal fall velocity of a general ellipsoidal particle
!  (either prolate or oblate) given the density and viscosity of air along with the
!  density, diameter and shape of the particle.  Velocity is returned in m/s.
!  The drag coefficient is that of Ganser (1993) which uses two shape factors to
!  characterize the type of ellipsoid.  As in the Wilson and Huang models the
!  drag coefficient and fall velocity must be solved iteritively with Cd and Vs
!  updated until Vs does not change by more than a small fraction vset_ConvCrit (~0.001).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function vset_Gans(rho_air,rho_m,eta,diam,K1,K2)
      ! Fall velocity as calculated from
      ! Ganser, Powder Tech., v77,2p143, 1993
      ! DOI:10.1016/0032-5910(93)80051-B

      real(kind=ip) :: vset_Gans ! Settling velocity in m/s
      real(kind=ip) :: rho_air   ! density of air in kg/m3
      real(kind=ip) :: rho_m     ! density of the particle in km/m3
      real(kind=ip) :: eta       ! dynamic viscosity of air in (kg/(m s))
      real(kind=ip) :: diam      ! diameter of the particle in m
      real(kind=ip) :: K1,K2     ! Stokes shape factor and Newtons shape factor

      real(kind=ip) :: vnew, vold, Re            ! old and new settling velocity
      real(kind=ip) :: Cd                        ! drag coefficient

      real(kind=ip) :: Cd1,Cd2

      do io=1,2;if(VB(io).le.verbosity_debug2)then
        write(outlog(io),*)"     Entered function vset_Gans"
      endif;enddo

      vold = 1.0_ip                              ! assume an initial settling velocity of 1 m/s

      Re = rho_air*vold*diam/eta                 ! Eq. 4  of Wilson79

      Cd1 = 1.0_ip + 0.1118_ip * (Re*K1*K2)**0.6567_ip
      Cd1 = Cd1 * 24.0_ip /(K1*Re)
      Cd2 = 0.4305_ip*K2/(1.0_ip+3305.0_ip/(Re*K1*K2))
      Cd = Cd1+Cd2;
      vnew = sqrt((4.0_ip*rho_m*diam*GRAV)/(3.0_ip*rho_air*Cd))   ! Eq. 15 of WoodsBursik91

      do while ((abs(vold-vnew)/vnew).gt.vset_ConvCrit)
        vold = vnew
        Re = rho_air*vold*diam/eta
        Cd1 = 1.0_ip+0.1118_ip*(Re*K1*K2)**0.6567_ip
        Cd1 = Cd1 * 24.0_ip /(K1*Re)
        Cd2 = 0.4305_ip*K2/(1.0_ip+3305.0_ip/(Re*K1*K2))
        Cd = Cd1+Cd2;
        vnew = sqrt((4.0_ip*rho_m*diam*GRAV)/(3.0_ip*rho_air*Cd))
      enddo

      vset_Gans = vnew

      return

      end function vset_Gans

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  vset_Gans_slip(rho_air,rho_m,eta,diam,K1,K2,Kna)
!
!  Called from: Set_Vf_Meso
!  Arguments:
!    rho_air= density of air (kg/m^3)
!    rho_m  = density of particle (kg/m^3)
!    eta    = viscosity of air (kg/(m s))
!    diam   = diameter of particle (m)
!    K1     = precalculated Stokes shape factor
!    K2     = precalculated Newtons shape factor
!    Kna    = Knudsen number (dimensionless)
!
!  This function returns the terminal fall velocity of a general ellipsoidal particle
!  (either prolate or oblate) given the density and viscosity of air along with the
!  density, diameter, shape of the particle, and the Knudsen number.  Velocity is returned
!  in m/s. The drag coefficient is that of Ganser (1993) which uses two shape factors to
!  characterize the type of ellipsoid.  As in the Wilson and Huang models the
!  drag coefficient and fall velocity must be solved iteritively with Cd and Vs
!  updated until Vs does not change by more than a small fraction vset_ConvCrit (~0.001).
!  The Knudsen number is important at high altitudes with low air density and a large
!  mean free path.  In these environments, Cunningham slip corrections should be
!  applied to the fall model. This solution is only valid for very small particles.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function vset_Gans_slip(rho_air,rho_m,eta,diam,K1,K2,Kna)
      ! Fall velocity as calculated from
      ! Ganser, Powder Tech., v77,2p143, 1993
      ! DOI:10.1016/0032-5910(93)80051-B

      real(kind=ip) :: vset_Gans_slip ! Settling velocity in m/s
      real(kind=ip) :: rho_air        ! density of air in kg/m3
      real(kind=ip) :: rho_m          ! density of the particle in km/m3
      real(kind=ip) :: eta            ! dynamic viscosity of air in (kg/(m s))
      real(kind=ip) :: diam           ! diameter of the particle in m
      real(kind=ip) :: K1,K2          ! Stokes shape factor and Newtons shape factor
      real(kind=ip) :: Kna            ! adjusted Knudsen number

      real(kind=ip) :: Cslip                     ! drag coefficient

      real(kind=ip) :: vnew, vold, Re            ! old and new settling velocity
      real(kind=ip) :: Cd                        ! drag coefficient

      real(kind=ip) :: Cd1,Cd2

      do io=1,2;if(VB(io).le.verbosity_debug2)then
        write(outlog(io),*)"     Entered function vset_Gans_slip"
      endif;enddo

      vold = 1.0_ip                              ! assume an initial settling velocity of 1 m/s

      Re = rho_air*vold*diam/eta                 ! Eq. 4  of Wilson79

      ! Eq. 9.34 of Seinfeld and Pandis
      Cslip = 1.0_ip+ Kna * (1.257_ip + 0.4_ip*exp(-1.1_ip/Kna))

      Cd1 = 1.0_ip + 0.1118_ip * (Re*K1*K2)**0.6567_ip
      Cd1 = Cd1 * 24.0_ip /(K1*Re)
      Cd2 = 0.4305_ip*K2/(1.0_ip+3305.0_ip/(Re*K1*K2))
      Cd = Cd1+Cd2;
      vnew = sqrt((4.0_ip*rho_m*diam*GRAV)/(3.0_ip*rho_air*Cd))   ! Eq. 15 of WoodsBursik91

      do while ((abs(vold-vnew)/vnew).gt.vset_ConvCrit)
        vold = vnew
        Re = rho_air*vold*diam/eta
        Cd1 = 1.0_ip+0.1118_ip*(Re*K1*K2)**0.6567_ip
        Cd1 = Cd1 * 24.0_ip /(K1*Re)
        Cd2 = 0.4305_ip*K2/(1.0_ip+3305.0_ip/(Re*K1*K2))
        Cd = Cd1+Cd2;
        Cd = Cd/Cslip
        vnew = sqrt((4.0_ip*rho_m*diam*GRAV)/(3.0_ip*rho_air*Cd))
      enddo

      vset_Gans_slip = vnew

      return

      end function vset_Gans_slip

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  vset_Stokes_slip(rho_m,eta,diam,Kna)
!
!  Called from: Set_Vf_Meso
!  Arguments:
!    rho_m  = density of particle (kg/m^3)
!    eta    = viscosity of air (kg/(m s))
!    diam   = diameter of particle (m)
!    Kna    = Knudsen number (dimensionless)
!
!  This function returns the Stokes fall velocity of a spherical particle
!  given the viscosity of air along with the density and diameter of the
!  particle and Knudsen number.  The Knudsen number is important at high
!  altitudes with low air density and a large mean free path.  In these
!  environments, Cunningham slip corrections should be applied to the fall model.
!  Velocity is returned in m/s. This solution is only valid for very small
!  particles.  Fall velocity is calculated explicitly and no iteration is
!  needed.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function vset_Stokes_slip(rho_m,eta,diam,Kna)
      ! Fall velocity as calculated from Stokes flow plus Cunningham slip

      real(kind=ip) :: vset_Stokes_slip  ! Settling velocity in m/s
      real(kind=ip) :: rho_m             ! density of the particle in km/m3
      real(kind=ip) :: eta               ! dynamic viscosity of air in (kg/(m s))
      real(kind=ip) :: diam              ! diameter of the particle in m
      real(kind=ip) :: Kna               ! adjusted Knudsen number

      real(kind=ip) :: Cslip                     ! drag coefficient

      do io=1,2;if(VB(io).le.verbosity_debug2)then
        write(outlog(io),*)"     Entered function vset_Stokes_slip"
      endif;enddo

      ! Eq. 9.34 of Seinfeld and Pandis
      Cslip = 1.0_ip+ Kna * (1.257_ip + 0.4_ip*exp(-1.1_ip/Kna))

      vset_Stokes_slip = diam*diam*rho_m*GRAV*Cslip/(18.0_ip*eta)

      return

      end function vset_Stokes_slip

      end module Tephra

!##############################################################################
