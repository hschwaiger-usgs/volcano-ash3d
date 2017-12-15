      module Atmosphere

      use precis_param

      use io_units

      real(kind=ip), parameter :: R_GAS_IDEAL  = 8.3144621_ip ! Ideal gas constant (J /(kg K))
      real(kind=ip), parameter :: R_GAS_DRYAIR = 286.98_ip    ! Specific gas constant of R=286.98 J /(kg K)
      real(kind=ip), parameter :: CP_AIR       = 1.004e3_ip   ! Specific heat capacity at p (J /kg K)
      real(kind=ip), parameter :: MB_DRY_AIR   = 0.028966_ip  ! Molecular weight of dry air in kg/mol
      real(kind=ip), parameter :: BoltzK       = 1.380658e-23_ip ! Boltzmann's constant kg m2 s-2 K-1 molec-1

      ! Define atmospheric variables only needed on the native Met grid
      !  Physical Properties of air
      !  meso means grid of the met. file
#ifdef USEPOINTERS
      real(kind=sp),dimension(:,:,:),pointer :: AirTemp_meso_last_step_MetP_sp => null()
      real(kind=sp),dimension(:,:,:),pointer :: AirDens_meso_last_step_MetP_sp => null()
      real(kind=sp),dimension(:,:,:),pointer :: AirVisc_meso_last_step_MetP_sp => null()
      real(kind=sp),dimension(:,:,:),pointer :: AirLamb_meso_last_step_MetP_sp => null()
      real(kind=sp),dimension(:,:,:),pointer :: AirRelH_meso_last_step_MetP_sp => null()
      real(kind=sp),dimension(:,:,:),pointer :: AirSH_meso_last_step_MetP_sp   => null()
      real(kind=sp),dimension(:,:,:),pointer :: AirTemp_meso_next_step_MetP_sp => null()
      real(kind=sp),dimension(:,:,:),pointer :: AirDens_meso_next_step_MetP_sp => null()
      real(kind=sp),dimension(:,:,:),pointer :: AirVisc_meso_next_step_MetP_sp => null()
      real(kind=sp),dimension(:,:,:),pointer :: AirLamb_meso_next_step_MetP_sp => null()
      real(kind=sp),dimension(:,:,:),pointer :: AirRelH_meso_next_step_MetP_sp => null()
      real(kind=sp),dimension(:,:,:),pointer :: AirSH_meso_next_step_MetP_sp   => null()
#else
      real(kind=sp),dimension(:,:,:),allocatable :: AirTemp_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:),allocatable :: AirDens_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:),allocatable :: AirVisc_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:),allocatable :: AirLamb_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:),allocatable :: AirRelH_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:),allocatable :: AirSH_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:),allocatable :: AirTemp_meso_next_step_MetP_sp
      real(kind=sp),dimension(:,:,:),allocatable :: AirDens_meso_next_step_MetP_sp
      real(kind=sp),dimension(:,:,:),allocatable :: AirVisc_meso_next_step_MetP_sp
      real(kind=sp),dimension(:,:,:),allocatable :: AirLamb_meso_next_step_MetP_sp
      real(kind=sp),dimension(:,:,:),allocatable :: AirRelH_meso_next_step_MetP_sp
      real(kind=sp),dimension(:,:,:),allocatable :: AirSH_meso_next_step_MetP_sp
#endif

      ! And since we are looking more carefully at atmospheric conditions,
      ! allocate the winds on the MetP grid
      !real(kind=sp),dimension(:,:,:),allocatable :: vx_meso_last_step_MetP_sp
      !real(kind=sp),dimension(:,:,:),allocatable :: vy_meso_last_step_MetP_sp
      !real(kind=sp),dimension(:,:,:),allocatable :: vz_meso_last_step_MetP_sp
      !real(kind=sp),dimension(:,:,:),allocatable :: vx_next_last_step_MetP_sp
      !real(kind=sp),dimension(:,:,:),allocatable :: vy_nxet_last_step_MetP_sp
      !real(kind=sp),dimension(:,:,:),allocatable :: vz_nxet_last_step_MetP_sp

      contains

!******************************************************************************

      subroutine Allocate_Atmosphere_Met

      use global_param,  only : &
         useMoistureVars

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet,np_fullmet_RH

      implicit none
      
      write(global_production,*)"--------------------------------------------------"
      write(global_production,*)"---------- ALLOCATE_ATMOSPHERE_MET ---------------"
      write(global_production,*)"--------------------------------------------------"

#ifdef USEPOINTERS
      if(.not.associated(AirTemp_meso_last_step_MetP_sp))&
        allocate(AirTemp_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if(.not.associated(AirDens_meso_last_step_MetP_sp))&
        allocate(AirDens_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if(.not.associated(AirVisc_meso_last_step_MetP_sp))&
        allocate(AirVisc_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if(.not.associated(AirLamb_meso_last_step_MetP_sp))&
        allocate(AirLamb_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if(useMoistureVars)THEN
        if(.not.associated(AirRelH_meso_last_step_MetP_sp))&
          allocate(AirRelH_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.associated(AirSH_meso_last_step_MetP_sp))&
          allocate(AirSH_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      endif

      if(.not.associated(AirTemp_meso_next_step_MetP_sp))&
        allocate(AirTemp_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if(.not.associated(AirDens_meso_next_step_MetP_sp))&
        allocate(AirDens_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if(.not.associated(AirVisc_meso_next_step_MetP_sp))&
        allocate(AirVisc_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if(.not.associated(AirLamb_meso_next_step_MetP_sp))&
        allocate(AirLamb_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if(useMoistureVars)THEN
        if(.not.associated(AirRelH_meso_next_step_MetP_sp))&
          allocate(AirRelH_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet_RH))
        if(.not.associated(AirSH_meso_next_step_MetP_sp))&
          allocate(AirSH_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet_RH))
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
          allocate(AirRelH_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet_RH))
        if(.not.allocated(AirSH_meso_next_step_MetP_sp))&
          allocate(AirSH_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet_RH))
      endif
#endif

      end subroutine Allocate_Atmosphere_Met

!******************************************************************************

      subroutine Deallocate_Atmosphere_Met

      use global_param,  only : &
         useMoistureVars

      implicit none

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

!******************************************************************************

      subroutine Set_Atmosphere_Meso(Load_MesoSteps,Interval_Frac,first_time)
      !Fclaw subroutine Set_Atmosphere_Meso(Load_MesoSteps,Interval_Frac,first_time)

      use global_param,  only : &
         useMoistureVars

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet,p_fullmet_sp,MR_dum3d_MetP,MR_iMetStep_Now,&
         Met_var_IsAvailable, &
           MR_Read_3d_MetP_Variable

      implicit none

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
                !HFS : need to convert RH to SH
              !AirSH_meso_last_step_MetP_sp = MR_dum3d_MetP
              AirSH_meso_last_step_MetP_sp = 0.0_sp
            else
              write(global_info,*)"ERROR: Neither SH nor RH are available"
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
                !HFS : need to convert RH to SH
            !AirSH_meso_next_step_MetP_sp = MR_dum3d_MetP
            AirSH_meso_last_step_MetP_sp = 0.0_sp
          else
            write(global_info,*)"ERROR: Neither SH nor RH are available"
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

      endif

      !Temperature(:,:,:) = real( Temp_meso_last_step_sp(:,:,:),kind=ip) + &
      !                     real((Temp_meso_next_step_sp(:,:,:) - &
      !                           Temp_meso_last_step_sp(:,:,:)),kind=ip) * &
      !                     Interval_Frac

      end subroutine Set_Atmosphere_Meso

!******************************************************************************

      function Dens_IdealGasLaw(pres,temp)

      implicit none

      real(kind=sp) :: Dens_IdealGasLaw
      real(kind=sp) :: pres,temp

        ! Get the density (kg/m^3) of dry air via the ideal gas law and
        ! Specific gas constant of R=286.98 J /(kg K)
      Dens_IdealGasLaw = pres/(real(R_GAS_DRYAIR,kind=sp)*temp)

      return

      end function Dens_IdealGasLaw

!******************************************************************************

      function Visc_Sutherland(temp)

      implicit none

      real(kind=sp) :: Visc_Sutherland
      real(kind=sp) :: temp

            ! Get the dynamic viscosity (kg/(m s)) of air via Sutherland's
            ! equation (Jacobson05 p. 102 Eq 4.54))
      Visc_Sutherland = 1.8325e-5_sp * (4.1616e2_sp/(temp+120.0_sp)) &
                             * (temp/296.16_sp)**1.5_sp
      return

      end function Visc_Sutherland

!******************************************************************************

      function lambda_MeanFreePath(visc,pres,temp)

      use global_param,  only : &
         PI

      implicit none

      real(kind=sp) :: lambda_MeanFreePath
      real(kind=sp) :: visc,pres,temp

        ! Mean-free-path of dry air : Eq. 9.6 of Seinfeld and Pandis
      lambda_MeanFreePath = (2.0_sp*visc)/pres * &
             sqrt(8.0_sp*real(MB_DRY_AIR/(PI*R_GAS_IDEAL),kind=sp)/temp)

      return

      end function lambda_MeanFreePath

!******************************************************************************

      end module Atmosphere
