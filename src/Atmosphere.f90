!##############################################################################
!
!  Atmosphere module
!
!  This module contains the variables storing atmospheric variables such as
!  density, viscosity, mean-free-path, as well as physical constants related
!  to gas equations.
!
!      subroutine Allocate_Atmosphere_Met
!      subroutine Deallocate_Atmosphere_Met
!      subroutine Set_Atmosphere_Meso(Load_MesoSteps,Interval_Frac,first_time)
!      function Dens_IdealGasLaw(pres,temp)
!      function Visc_Sutherland(temp)
!      function lambda_MeanFreePath(visc,pres,temp)
!
!##############################################################################

      module Atmosphere

      use precis_param

      use io_units

      implicit none

        ! Set everything to private by default
      private

        ! Publicly available subroutines/functions
      public Allocate_Atmosphere_Met,   &
             Deallocate_Atmosphere_Met, &
             Set_Atmosphere_Meso

        ! Publicly available variables
        ! Define atmospheric variables only needed on the native Met grid
        !  Physical Properties of air
        !  meso means grid of the met. file
#ifdef USEPOINTERS
        ! Denisity is in kg/m^3
        ! Viscosity is in kg/(m s)
        ! Mean free path (lambda) is in m
      real(kind=sp),dimension(:,:,:),pointer,public :: AirDens_meso_last_step_MetP_sp => null()
      real(kind=sp),dimension(:,:,:),pointer,public :: AirVisc_meso_last_step_MetP_sp => null()
      real(kind=sp),dimension(:,:,:),pointer,public :: AirLamb_meso_last_step_MetP_sp => null()
      real(kind=sp),dimension(:,:,:),pointer,public :: AirDens_meso_next_step_MetP_sp => null()
      real(kind=sp),dimension(:,:,:),pointer,public :: AirVisc_meso_next_step_MetP_sp => null()
      real(kind=sp),dimension(:,:,:),pointer,public :: AirLamb_meso_next_step_MetP_sp => null()
      real(kind=sp),dimension(:,:,:),pointer,public :: AirTemp_meso_last_step_MetP_sp => null()
      real(kind=sp),dimension(:,:,:),pointer,public :: AirRelH_meso_last_step_MetP_sp => null()
      real(kind=sp),dimension(:,:,:),pointer,public :: AirSH_meso_last_step_MetP_sp   => null()
      real(kind=sp),dimension(:,:,:),pointer,public :: AirTemp_meso_next_step_MetP_sp => null()
      real(kind=sp),dimension(:,:,:),pointer,public :: AirRelH_meso_next_step_MetP_sp => null()
      real(kind=sp),dimension(:,:,:),pointer,public :: AirSH_meso_next_step_MetP_sp   => null()
#else
      real(kind=sp),dimension(:,:,:),allocatable,public :: AirDens_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:),allocatable,public :: AirVisc_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:),allocatable,public :: AirLamb_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:),allocatable,public :: AirDens_meso_next_step_MetP_sp
      real(kind=sp),dimension(:,:,:),allocatable,public :: AirVisc_meso_next_step_MetP_sp
      real(kind=sp),dimension(:,:,:),allocatable,public :: AirLamb_meso_next_step_MetP_sp
      real(kind=sp),dimension(:,:,:),allocatable,public :: AirTemp_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:),allocatable,public :: AirRelH_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:),allocatable,public :: AirSH_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:),allocatable,public :: AirTemp_meso_next_step_MetP_sp
      real(kind=sp),dimension(:,:,:),allocatable,public :: AirRelH_meso_next_step_MetP_sp
      real(kind=sp),dimension(:,:,:),allocatable,public :: AirSH_meso_next_step_MetP_sp
#endif

      real(kind=ip), parameter,public :: R_GAS_DRYAIR = 286.98_ip       ! Specific gas constant of R=286.98 J /(kg K)
      real(kind=ip), parameter,public :: R_GAS_IDEAL  = 8.3144621_ip    ! Ideal gas constant (J /(kg K))
      real(kind=ip), parameter,public :: CP_AIR       = 1.004e3_ip      ! Specific heat capacity at p (J /kg K)
      real(kind=ip), parameter,public :: MB_DRY_AIR   = 0.028966_ip     ! Molecular weight of dry air in kg/mol
      real(kind=ip), parameter,public :: BoltzK       = 1.380658e-23_ip ! Boltzmann's constant kg m2 s-2 K-1 molec-1

      ! HFS : Probably delete this
      ! And since we are looking more carefully at atmospheric conditions,
      ! allocate the winds on the MetP grid
      !real(kind=sp),dimension(:,:,:),allocatable :: vx_meso_last_step_MetP_sp
      !real(kind=sp),dimension(:,:,:),allocatable :: vy_meso_last_step_MetP_sp
      !real(kind=sp),dimension(:,:,:),allocatable :: vz_meso_last_step_MetP_sp
      !real(kind=sp),dimension(:,:,:),allocatable :: vx_next_last_step_MetP_sp
      !real(kind=sp),dimension(:,:,:),allocatable :: vy_nxet_last_step_MetP_sp
      !real(kind=sp),dimension(:,:,:),allocatable :: vz_nxet_last_step_MetP_sp

      contains
      !------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Allocate_Atmosphere_Met
!
!  Called from: Main ash3d (Ash3d.F90)
!  Arguments:
!    none
!
!  Allocates all atmospheric variables on the subset of the met grid with
!  p and the vertical axis
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Allocate_Atmosphere_Met

      use global_param,  only : &
         useMoistureVars

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet

      do io=1,2;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)"--------------------------------------------------"
        write(outlog(io),*)"---------- ALLOCATE_ATMOSPHERE_MET ---------------"
        write(outlog(io),*)"--------------------------------------------------"
      endif;enddo

      ! Note: all these arrays are allocated on np_fullmet levels even though
      !       np_fullmet_T or np_fullmet_RH might be less.  These will contain
      !       the full MR_dum3d_MetP array and MetReader will populate the missing
      !       values as best it can (SH,RH = 0 above levels given, T interpolated)

#ifdef USEPOINTERS
      allocate(AirTemp_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(AirDens_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(AirVisc_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(AirLamb_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))

      if(useMoistureVars)THEN
        allocate(AirRelH_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
        allocate(AirSH_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      endif

      allocate(AirTemp_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(AirDens_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(AirVisc_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(AirLamb_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if (useMoistureVars) THEN
        allocate(AirRelH_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
        allocate(AirSH_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      endif
#else
      if(.not.allocated(AirTemp_meso_last_step_MetP_sp))&
        allocate(AirTemp_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if(.not.allocated(AirDens_meso_last_step_MetP_sp))&
        allocate(AirDens_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if(.not.allocated(AirVisc_meso_last_step_MetP_sp))&
        allocate(AirVisc_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if(.not.allocated(AirLamb_meso_last_step_MetP_sp))&
        allocate(AirLamb_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if(useMoistureVars)THEN
        if(.not.allocated(AirRelH_meso_last_step_MetP_sp))&
          allocate(AirRelH_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.allocated(AirSH_meso_last_step_MetP_sp))&
          allocate(AirSH_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      endif
        
      if(.not.allocated(AirTemp_meso_next_step_MetP_sp))&
        allocate(AirTemp_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if(.not.allocated(AirDens_meso_next_step_MetP_sp))&
        allocate(AirDens_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if(.not.allocated(AirVisc_meso_next_step_MetP_sp))&
        allocate(AirVisc_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if(.not.allocated(AirLamb_meso_next_step_MetP_sp))&
        allocate(AirLamb_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if(useMoistureVars)THEN
        if(.not.allocated(AirRelH_meso_next_step_MetP_sp))&
          allocate(AirRelH_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.allocated(AirSH_meso_next_step_MetP_sp))&
          allocate(AirSH_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      endif
#endif

      end subroutine Allocate_Atmosphere_Met

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Deallocate_Atmosphere_Met
!
!  Called from: dealloc_arrays
!  Arguments:
!    none
!
!  Deallocates all atmospheric variables if they were allocated
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Deallocate_Atmosphere_Met

      use global_param,  only : &
         useMoistureVars

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Deallocate_Atmosphere_Met"
      endif;enddo

#ifdef USEPOINTERS
      if(associated(AirTemp_meso_last_step_MetP_sp))deallocate(AirTemp_meso_last_step_MetP_sp)
      if(associated(AirDens_meso_last_step_MetP_sp))deallocate(AirDens_meso_last_step_MetP_sp)
      if(associated(AirVisc_meso_last_step_MetP_sp))deallocate(AirVisc_meso_last_step_MetP_sp)
      if(associated(AirLamb_meso_last_step_MetP_sp))deallocate(AirLamb_meso_last_step_MetP_sp)
      if(useMoistureVars)THEN
        if(associated(AirRelH_meso_last_step_MetP_sp))deallocate(AirRelH_meso_last_step_MetP_sp)
        if(associated(AirSH_meso_last_step_MetP_sp))deallocate(AirSH_meso_last_step_MetP_sp)
      endif

      if(associated(AirTemp_meso_next_step_MetP_sp))deallocate(AirTemp_meso_next_step_MetP_sp)
      if(associated(AirDens_meso_next_step_MetP_sp))deallocate(AirDens_meso_next_step_MetP_sp)
      if(associated(AirVisc_meso_next_step_MetP_sp))deallocate(AirVisc_meso_next_step_MetP_sp)
      if(associated(AirLamb_meso_next_step_MetP_sp))deallocate(AirLamb_meso_next_step_MetP_sp)
      if(useMoistureVars)THEN
        if(associated(AirRelH_meso_next_step_MetP_sp))deallocate(AirRelH_meso_next_step_MetP_sp)
        if(associated(AirSH_meso_next_step_MetP_sp))deallocate(AirSH_meso_next_step_MetP_sp)
      endif
#else
      if(allocated(AirTemp_meso_last_step_MetP_sp))deallocate(AirTemp_meso_last_step_MetP_sp)
      if(allocated(AirDens_meso_last_step_MetP_sp))deallocate(AirDens_meso_last_step_MetP_sp)
      if(allocated(AirVisc_meso_last_step_MetP_sp))deallocate(AirVisc_meso_last_step_MetP_sp)
      if(allocated(AirLamb_meso_last_step_MetP_sp))deallocate(AirLamb_meso_last_step_MetP_sp)
      if(useMoistureVars)THEN
        if(allocated(AirRelH_meso_last_step_MetP_sp))deallocate(AirRelH_meso_last_step_MetP_sp)
        if(allocated(AirSH_meso_last_step_MetP_sp))deallocate(AirSH_meso_last_step_MetP_sp)
      endif

      if(allocated(AirTemp_meso_next_step_MetP_sp))deallocate(AirTemp_meso_next_step_MetP_sp)
      if(allocated(AirDens_meso_next_step_MetP_sp))deallocate(AirDens_meso_next_step_MetP_sp)
      if(allocated(AirVisc_meso_next_step_MetP_sp))deallocate(AirVisc_meso_next_step_MetP_sp)
      if(allocated(AirLamb_meso_next_step_MetP_sp))deallocate(AirLamb_meso_next_step_MetP_sp)
      if(useMoistureVars)THEN
        if(allocated(AirRelH_meso_next_step_MetP_sp))deallocate(AirRelH_meso_next_step_MetP_sp)
        if(allocated(AirSH_meso_next_step_MetP_sp))deallocate(AirSH_meso_next_step_MetP_sp)
      endif
#endif

      end subroutine Deallocate_Atmosphere_Met

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Set_Atmosphere_Meso(Load_MesoSteps,Interval_Frac,first_time)
!
!  Called from: MesoInterpolater
!  Arguments:
!    Load_MesoSteps = logical; triggers loading the next step
!    Interval_Frac  = fraction of the time between last and next met steps
!    first_time     = logical
!
!  Most []_Meso() subroutine interpolates data onto the current time and computational
!  grid.  For the atmospheric data, we only use data on the met steps and on the 
!  subgrid of the NWP grid.  These are used for calculating fall velocities at these
!  met grid nodes, which are then interpolated onto the computational grid. So this
!  subroutine is only called when a new met step needs to be loaded, populating density,
!  viscosity and mean-free path on those met nodes.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Set_Atmosphere_Meso(Load_MesoSteps,Interval_Frac,first_time)
      !Fclaw subroutine Set_Atmosphere_Meso(Load_MesoSteps,Interval_Frac,first_time)

      use global_param,  only : &
         useMoistureVars

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet,p_fullmet_sp,MR_dum3d_MetP,MR_iMetStep_Now,&
         Met_var_IsAvailable, &
           MR_Read_3d_MetP_Variable

      logical      ,intent(in) :: Load_MesoSteps
      real(kind=ip),intent(in) :: Interval_Frac     ! This is a placeholder in cases where
                                                    ! atmospheric variables are interpolated
                                                    ! onto local time steps
      logical      ,intent(in) :: first_time

      integer :: ivar
      integer :: i,j,k
      real(kind=sp) :: temp,pres

      !logical :: first_time
      ! Fclaw
      !logical,save :: first_time = .true.

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Set_Atmosphere_Meso"
      endif;enddo

      if(Load_MesoSteps)THEN
        if(first_time)THEN
          ivar = 5 ! Temperature
          call MR_Read_3d_MetP_Variable(ivar,MR_iMetStep_Now)
            AirTemp_meso_next_step_MetP_sp = MR_dum3d_MetP
          if(useMoistureVars)THEN
            if(Met_var_IsAvailable(31))then
              ! Read Specific humidity if it is available
              ivar = 31
              call MR_Read_3d_MetP_Variable(ivar,MR_iMetStep_Now)
              AirSH_meso_last_step_MetP_sp = MR_dum3d_MetP
            elseif(Met_var_IsAvailable(30))then
              ! Otherwise, read Rel Humidity and convert
              ivar = 30
              call MR_Read_3d_MetP_Variable(ivar,MR_iMetStep_Now)
                ! need to convert RH to SH
              !AirSH_meso_last_step_MetP_sp = MR_dum3d_MetP
              do io=1,2;if(VB(io).le.verbosity_info)then
                write(outlog(io),*)"WARNING: Specific Humidity requested but is unavailable."
                write(outlog(io),*)"         Setting to zero."
              endif;enddo
              AirSH_meso_last_step_MetP_sp = 0.0_sp
            else
              do io=1,2;if(VB(io).le.verbosity_error)then
                write(errlog(io),*)"ERROR: Neither SH nor RH are available"
              endif;enddo
              stop 1
            endif
          endif
          ! Fclaw  : delete this line
          !first_time = .false.
        endif ! first_time
        AirTemp_meso_last_step_MetP_sp = AirTemp_meso_next_step_MetP_sp
        if(useMoistureVars)THEN
          AirRelH_meso_last_step_MetP_sp = AirRelH_meso_next_step_MetP_sp
          AirSH_meso_last_step_MetP_sp   = AirSH_meso_next_step_MetP_sp
        endif

        ivar = 5 ! Temperature
        call MR_Read_3d_MetP_Variable(ivar,MR_iMetStep_Now+1)
        AirTemp_meso_next_step_MetP_sp = MR_dum3d_MetP
        if(useMoistureVars)THEN
          if(Met_var_IsAvailable(31))then
            ! Read Specific humidity if it is available
            ivar = 31
            call MR_Read_3d_MetP_Variable(ivar,MR_iMetStep_Now+1)
            AirSH_meso_next_step_MetP_sp = MR_dum3d_MetP
          elseif(Met_var_IsAvailable(30))then
            ! Otherwise, read Rel Humidity and convert
            ivar = 30
            call MR_Read_3d_MetP_Variable(ivar,MR_iMetStep_Now+1)
                ! need to convert RH to SH
            !AirSH_meso_next_step_MetP_sp = MR_dum3d_MetP
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"WARNING: Specific Humidity requested but is unavailable."
              write(outlog(io),*)"         Setting to zero."
            endif;enddo
            AirSH_meso_last_step_MetP_sp = 0.0_sp
          else
            do io=1,2;if(VB(io).le.verbosity_error)then  
              write(errlog(io),*)"ERROR: Neither SH nor RH are available"
            endif;enddo
            stop 1
          endif
        endif
        do k=1,np_fullmet
          ! Note: this needs to be fixed for WRF data
          pres = p_fullmet_sp(k)
          do i=1,nx_submet
            do j=1,ny_submet
              temp = AirTemp_meso_last_step_MetP_sp(i,j,k)
              AirDens_meso_last_step_MetP_sp(i,j,k) = &
                Dens_IdealGasLaw(pres,temp)
              AirVisc_meso_last_step_MetP_sp(i,j,k) = &
                Visc_Sutherland(temp)
              AirLamb_meso_last_step_MetP_sp(i,j,k) = &
                lambda_MeanFreePath(AirVisc_meso_last_step_MetP_sp(i,j,k),&
                                    pres,temp)

              temp = AirTemp_meso_next_step_MetP_sp(i,j,k)
              AirDens_meso_next_step_MetP_sp(i,j,k) = &
                Dens_IdealGasLaw(pres,temp)
              AirVisc_meso_next_step_MetP_sp(i,j,k) = &
                Visc_Sutherland(temp)
              AirLamb_meso_next_step_MetP_sp(i,j,k) = &
                lambda_MeanFreePath(AirVisc_meso_next_step_MetP_sp(i,j,k),&
                                    pres,temp)
            enddo
          enddo
        enddo
      else
        ! Currently, this subroutine is only called if Load_MesoSteps=.true.
        ! so we shouldn't be here
        do io=1,2;if(VB(io).le.verbosity_error)then              
          write(errlog(io),*)"Calling Set_Atmosphere_Meso outside of a Load_MesoSteps=.true."
          write(errlog(io),*)"case for Interval_Frac = ",Interval_Frac
          write(errlog(io),*)"This is a place-holder for interpolating tempertures to the"
          write(errlog(io),*)"current time.  Not yet implemented."
        endif;enddo
        stop 1
        !Temperature(:,:,:) = real( Temp_meso_last_step_sp(:,:,:),kind=ip) + &
        !                     real((Temp_meso_next_step_sp(:,:,:) - &
        !                           Temp_meso_last_step_sp(:,:,:)),kind=ip) * &
        !                     Interval_Frac
      endif

      end subroutine Set_Atmosphere_Meso

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Dens_IdealGasLaw(pres,temp)
!
!  Called from: Set_Atmosphere_Meso
!  Arguments:
!    pres = pressure (Pa)
!    temp = temperature (K)
!
!  Function that calculates air density give pressure and temperature using
!  the ideal gas law.  Density is returned in kg/m^3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function Dens_IdealGasLaw(pres,temp)

      real(kind=sp) :: Dens_IdealGasLaw
      real(kind=sp) :: pres,temp

      do io=1,2;if(VB(io).le.verbosity_debug2)then
        write(outlog(io),*)"     Entered function Dens_IdealGasLaw"
      endif;enddo

        ! Get the density (kg/m^3) of dry air via the ideal gas law and
        ! Specific gas constant of R=286.98 J /(kg K)
      Dens_IdealGasLaw = pres/(real(temp*R_GAS_DRYAIR,kind=sp))

      return

      end function Dens_IdealGasLaw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Visc_Sutherland(temp)
!
!  Called from: Set_Atmosphere_Meso
!  Arguments:
!    temp = temperature (K)
!
!  Function that calculates air viscosity given temperature in K via
!  Sutherland's equation.  Viscosity is returned in kg/(m s)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function Visc_Sutherland(temp)

      real(kind=sp) :: Visc_Sutherland
      real(kind=sp) :: temp

      do io=1,2;if(VB(io).le.verbosity_debug2)then
        write(outlog(io),*)"     Entered function Visc_Sutherland"
      endif;enddo

            ! Get the dynamic viscosity (kg/(m s)) of air via Sutherland's
            ! equation (Jacobson05 p. 102 Eq 4.54))
      Visc_Sutherland = 1.8325e-5_sp * (4.1616e2_sp/(temp+120.0_sp)) &
                             * (temp/296.16_sp)**1.5_sp
      return

      end function Visc_Sutherland

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  lambda_MeanFreePath(visc,pres,temp)
!
!  Called from: Set_Atmosphere_Meso
!  Arguments:
!    visc = viscosity kg/(m s)
!    pres = pressure (Pa)
!    temp = temperature (K)
!
!  Function that calculates the mean-free path of air given air viscosity,
!  pressure and temperature using Eq. 9.6 of Seinfeld and Pandis. Output
!  units are in m.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function lambda_MeanFreePath(visc,pres,temp)

      use global_param,  only : &
         PI

      real(kind=sp) :: lambda_MeanFreePath
      real(kind=sp) :: visc,pres,temp

      do io=1,2;if(VB(io).le.verbosity_debug2)then
        write(outlog(io),*)"     Entered function lambda_MeanFreePath"
      endif;enddo

        ! Mean-free-path of dry air : Eq. 9.6 of Seinfeld and Pandis
      lambda_MeanFreePath = (2.0_sp*visc)/ &
             (pres*sqrt(8.0_sp*real(MB_DRY_AIR/(PI*R_GAS_IDEAL*temp),kind=sp)))

      return

      end function lambda_MeanFreePath

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module Atmosphere

!##############################################################################
