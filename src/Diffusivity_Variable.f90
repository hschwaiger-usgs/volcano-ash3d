!##############################################################################
!
!  Diffusivity_Variable module
!
!  This module calculates the horizontal and vertical diffusivities as a
!  function of the local meteorological conditions.
!
!  Horizontal diffusivity is calculated using spatial derivatives of the
!  horizontal velocities using either the model of Smagorinsky (1963) or
!  from Pielke (1974). This LES approach uses the area of the NWP cells
!  to scale.
!
!  Vertical diffusivity is calculated using vertical derivatives of horizontal
!  velocities, scaled by a function of the stability parameter Ri (itself a
!  function of horizontal velocities and potential temperature). Additionally,
!  a separate calculation of diffusivities in the planetary boundary layer can
!  be invoked. This requires Ri as well, but also the thickness of the planetary
!  boundary layer (provided by the NWP, calculated from Ri, T inversion, or
!  Eckman Layer thickness), and the surface friction velocity. Several
!  options for characterizing the diffusivity in this boundary layer are
!  available through the specification of the similarity function Phi, as
!  well as the functional form of the extension from the surface layer throughout
!  the boundary layer.
!
!      subroutine input_data_VarDiff
!      subroutine Allocate_VarDiff_Met
!      subroutine Prep_output_VarDiff
!      subroutine Deallocate_VarDiff_Met
!      subroutine Eddy_diff
!      subroutine Calc_Vert_Diff
!      subroutine Set_VarDiffH_Meso
!      subroutine Set_VarDiffV_Meso
!      subroutine Calc_Ri
!      subroutine Calc_SurfaceRoughnessLength
!      subroutine Calc_SurfaceFrictionVelocity
!      subroutine Calc_Monin_Length
!      subroutine Calc_PBLH
!      function Fc_Louis
!      function Fc_Jac
!      function Fc_Betts
!      function Fc_Hong
!      function Fc_Collins
!      function MixLen
!      function Phi_WindShear_Similarity
!      function Psi_WindShear_Similarity
!
!
! Line 1: Optional module specifier
! OPTMOD=VARDIFF
! Line 2: indicates whether or not to write VarDiff variables to the output file
!  yes
! Line 3: Horizontal diff
!  yes 1 500.0     # 1=const ; value in m2/s
!  yes 2 0.2       # 2=Smagorinsky ; C
!  yes 3 0.2       # 3=Pielke      ; C
! Line 4: Vertical diff
!  yes
! Line 5: Kv BL parameters
!  0 B, alpha, beta, gamma, pexp
!  1 500.0         # BL model 1=const ; value
!  2               #          2=no BL, only free-air throughout
!  3 [1,2]         #          3=Troen and Mahrt
!  4 [1,2]         #          4=Ulke
!  5               #          5=Shir / Businger,Ayer
! Line 6: Kv Free-air parameters
!  1 500.0         # Free-Air 1=const ; value
!  2               #          2=F(Ri)=Louis 1979
!  3               #          3=F(Ri)=Stull 1988
!  4               #          4=F(Ri)=Betts 1996
!  5               #          5=F(Ri)=Hong 1996
!  6               #          6=F(Ri,z)=Collins 2004
! Line 7: 
!0.4                         # vonKarman
! Line 8: 
!30.0                        # LambdaC
! Line 9:
!0.25                        # RI_CRIT
! 
!##############################################################################

      module Diffusivity_Variable

      use precis_param

      use io_units

      use Diffusion,     only : &
         diffusivity_horz,diffusivity_vert

      integer :: Kh_model_ID     ! [1] = Smagorinsky (1963); 2 = Pielke (1974)
      integer :: Phi_model_ID    ! 
      integer :: KvBL_model_ID   ! 
      integer :: KvBL_MomHeat    !  
      integer :: KvFA_model_ID   ! 

      !  These are the parameters that control the diffusivity calculations
      !    C from Smagorinsky model of horizontal diffusivity
      real(kind=ip) :: KH_SmagC     ! Smagorinsky (1993) constant for LES horizontal diffusivity (0.2 - 0.9)
      real(kind=ip),parameter :: MAX_LES_LengthScale2 = 100.0_ip ! Maximum area that will be used for scaling
      !    These next three are needed for the vertical diffusivity
      real(kind=ip) :: vonKarman    ! von Karman constant (around 0.4)
      real(kind=ip) :: LambdaC      ! Asymptotic length scale (around 30-150 m)
      real(kind=ip) :: RI_CRIT      ! Critical Richardson number (0.25)

      real(kind=ip) :: USTAR_MIN = 0.1_ip   ! Minimum Friction Velocity (m/s)

      !    These are the values controlling the stability function Phi (lots of models out there)
      !   source              alpha   beta   gamma
      ! Businger-Dyer (1971)  -1/4    5.0    -16.0
      ! Carl (1973)           -1/3    5.0    -15.0
      ! Businger-Arya (1974)          4.7
      ! Troen-Mahrt (1986)    -1/3    4.7     -7.0 ** Default
      ! Ulke (2000)           -1/2    9.2    -13.0
      real(kind=ip) :: phi_prefac = 1.0_ip       ! typically 1 for momentum and Pr (<1) for heat
      real(kind=ip) :: phi_alpha = -0.33333_ip   ! Exponent in unstable term
      real(kind=ip) :: phi_beta  =  4.7_ip       ! Coefficient in stable term (pretty much always 4.7->5.2
      real(kind=ip) :: phi_gamma = -7.0_ip       ! Coefficient in unstable term
      integer       :: PBL_exp_int = 1
      real(kind=ip) :: PBL_exp

      logical       :: useBoundaryLayer = .true.
      real(kind=ip) :: diffusivity_BL            ! This is used if Kv is constant in BL
      ! Set the number of output variables for this module
      ! This depends on settings from the input block
      logical :: use_Output_Vars_VarDiff       = .true.
      integer, parameter :: nvar_User2d_static_XY_VarDiff = 0
      integer            :: nvar_User2d_XY_VarDiff        = 0 ! If using Kz, then =2 : Pblh, Ust
      integer, parameter :: nvar_User3d_XYGs_VarDiff      = 0
      integer            :: nvar_User3d_XYZ_VarDiff       = 0 ! If using Kh, then =1 khorz; if also Kz, then =3 kvert, Ri
      integer, parameter :: nvar_User4d_XYZGs_VarDiff     = 0
      integer, parameter,public :: nvar_User_charlines_VarDiff   = 9 ! number of line of the special block of control file

      character(len=30),dimension(:),allocatable :: temp_2d_name_VarDiff
      character(len=30),dimension(:),allocatable :: temp_2d_unit_VarDiff
      character(len=30),dimension(:),allocatable :: temp_2d_lname_VarDiff
      real(kind=op),    dimension(:),allocatable :: temp_2d_MissVal_VarDiff
      real(kind=op),    dimension(:),allocatable :: temp_2d_FillVal_VarDiff

      character(len=30),dimension(:),allocatable :: temp_3d_name_VarDiff
      character(len=30),dimension(:),allocatable :: temp_3d_unit_VarDiff
      character(len=30),dimension(:),allocatable :: temp_3d_lname_VarDiff
      real(kind=op),    dimension(:),allocatable :: temp_3d_MissVal_VarDiff
      real(kind=op),    dimension(:),allocatable :: temp_3d_FillVal_VarDiff
      character(len=80),dimension(nvar_User_charlines_VarDiff),public :: var_User_charlines_VarDiff = ''

      ! These are used to keep track of which index in the global list, this
      ! modules output vars corespond to
      integer :: indx_User2d_static_XY_VarDiff
      integer :: indx_User2d_XY_VarDiff
      integer :: indx_User3d_XYGs_VarDiff
      integer :: indx_User3d_XYZ_VarDiff
      integer :: indx_User4d_XYZGs_VarDiff
      integer :: indx_User_charlines_VarDiff

      real(kind=ip) :: LES_L2ScaleCoeff

      ! 3d Variables needed on MetP grid
!      real(kind=sp),dimension(:,:,:)  ,allocatable :: dVel_dz_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: du_dx_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: du_dy_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: dv_dx_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: dv_dy_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: dV_dz_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Khz_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Khz_meso_next_step_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Kv_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Kv_meso_next_step_MetP_sp

      ! 2d variables needed at meso steps on Met grid
      real(kind=sp),dimension(:,:)    ,allocatable :: SurfRoughLen_Met_sp

      ! Variables needed on Comp grid (kx,y,z are already allocated)
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Khz_meso_last_step_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Khz_meso_next_step_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Kv_meso_last_step_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Kv_meso_next_step_sp

      ! Both Khz and Kv need U and V values on MetP grid so store local copies
      ! Note: The core Ash3d code reads directly into the computational grid, bypassing
      !       the MetP grid
      real(kind=sp),dimension(:,:,:)  ,allocatable :: vx_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: vy_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: vx_meso_next_step_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: vy_meso_next_step_MetP_sp

      contains
      !------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  input_data_VarDiff
!
!  Called from: Ash3d.f90
!  Arguments: None
!    none
!
!  This subroutine reopens the Ash3d control file, searches for the optional
!  module identifier OPTMOD=VARDIFF then parses the block for user-specified
!  Variable Diffisivity options/parameters.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine input_data_VarDiff

      use global_param,  only : &
         nmods,useTemperature,useVarDiffH,useVarDiffV

      use io_data,       only : &
         infile

      use MetReader,     only : &
         MR_Save_Velocities

      implicit none

      character(len=3 )  :: answer
      character(len=80)  :: linebuffer080
      integer            :: ios,ios2,ioerr
      character(len=20)  :: mod_name
      integer            :: substr_pos
      real(kind=ip)      :: tmp

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine input_data_VarDiff"
      endif;enddo

      open(unit=10,file=infile,status='old',err=1900)

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"    Searching for OPTMOD=VARDIFF"
      endif;enddo
      nmods = 0
      read(10,'(a80)',iostat=ios)linebuffer080
      do while(ios.eq.0)
        read(10,'(a80)',iostat=ios)linebuffer080

        substr_pos = index(linebuffer080,'OPTMOD')
        if(substr_pos.eq.1)then
          ! found an optional module
          !  Parse for the keyword
          read(linebuffer080,1104)mod_name
          if(adjustl(trim(mod_name)).eq.'VARDIFF')then
            exit
          endif
        endif
1104    format(7x,a20)
      enddo
      ! Found Line 1:
      var_User_charlines_VarDiff(1) = trim(adjustl(linebuffer080))

      useVarDiffH = .false.
      useVarDiffV = .false.
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"    Continue reading input file for VarDiff block"
      endif;enddo

      ! Check if variables from the VarDiff module will be written to output file
      read(10,'(a80)',iostat=ios,err=2010)linebuffer080
      ! Line 2:
      var_User_charlines_VarDiff(2) = trim(adjustl(linebuffer080))
      read(linebuffer080,'(a3)',err=2011) answer
      if (answer.eq.'yes') then
        use_Output_Vars_VarDiff = .true.
      else
        use_Output_Vars_VarDiff = .false.
      endif

      !Check if we're going to use variable diffusivity
      read(10,'(a80)',iostat=ios,err=2010)linebuffer080
      ! Line 3:
      var_User_charlines_VarDiff(3) = trim(adjustl(linebuffer080))
      read(linebuffer080,'(a3)',err=2011) answer
      if (answer.eq.'yes') then
        useVarDiffH = .true.  ! might be changed back below if we are holding Kh constant

        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"    Horizontal variable diffusivity:  ON"
        endif;enddo
        ! Try to read the horizontal model ID
        read(linebuffer080(4:),*,iostat=ios)Kh_model_ID,tmp
        if(ios.eq.0)then
          if(Kh_model_ID.eq.1)then
            useVarDiffH = .false.
            diffusivity_horz = tmp
            ! Error-checking diffusivity
            if(diffusivity_horz.lt.0.0_ip)then
              do io=1,2;if(VB(io).le.verbosity_error)then
                write(errlog(io),*)"ERROR: diffusivity_horz must be >0"
              endif;enddo
              stop 1
            endif
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"        Kh model ID       = 1: Constant"
              write(outlog(io),*)"        with Kh (in m2/s) = ",real(diffusivity_horz,kind=4)
            endif;enddo
          elseif(Kh_model_ID.eq.2)then
            KH_SmagC = tmp
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"        Kh model ID  = 2: Smagorinsky (1963)"
              write(outlog(io),*)"              with C = ",real(KH_SmagC,kind=4)
            endif;enddo
          elseif(Kh_model_ID.eq.3)then
            KH_SmagC = tmp
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"        Kh model ID  = 3: Pielke (1974)"
              write(outlog(io),*)"              with C = ",real(KH_SmagC,kind=4)
            endif;enddo
          else
            KH_SmagC = 0.2_ip
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"        Horizontal diffusivity model ID not recognized."
              write(outlog(io),*)"        Using Kh model ID = 2: Smagorinsky (1963)"
              write(outlog(io),*)"                   with C = ",real(KH_SmagC,kind=4)
            endif;enddo
          endif
        else
          KH_SmagC = 0.2_ip
          Kh_model_ID = 2
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"        Error reading Kh model ID and/or parameter."
            write(outlog(io),*)"        Offending line:",linebuffer080
            write(outlog(io),*)"        Horizontal diffusivity model ID not recognized."
            write(outlog(io),*)"        Using Kh model ID = 2: Smagorinsky (1963)"
            write(outlog(io),*)"                   with C = ",real(KH_SmagC,kind=4)
          endif;enddo
        endif
      elseif(answer(1:2).eq.'no') then
        useVarDiffH = .false.
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"    Horizontal variable diffusivity:  OFF"
        endif;enddo
      else
        goto 2011
      endif

      ! Vertical variable diffusivity
      read(10,'(a80)',iostat=ios,err=2010)linebuffer080
      ! Line 4:
      var_User_charlines_VarDiff(4) = trim(adjustl(linebuffer080))
      read(linebuffer080,'(a3)',err=2011) answer
      if (answer.eq.'yes') then
        useVarDiffV    = .true.
        useTemperature = .true.
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"    Vertical variable diffusivity:    ON"
        endif;enddo
      elseif(answer(1:2).eq.'no') then
        useVarDiffV = .false.
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"    Vertical variable diffusivity:    OFF"
        endif;enddo
      else
        goto 2011
      endif
      if(useVarDiffV)then
        ! Need to read two more lines defining first the Boundary Layer model, then the Free-Air model
        read(10,'(a80)',iostat=ios,err=2010)linebuffer080
        ! Line 5:
        var_User_charlines_VarDiff(5) = trim(adjustl(linebuffer080))
        read(linebuffer080,*,iostat=ios)KvBL_model_ID
        if(ios.ne.0)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: ",&
                  'Boundary Layer model ID not provided.'
          endif;enddo
          stop 1
        endif
        if(KvBL_model_ID.eq.1)then
          ! Diffusivity is constant in the BL
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using constant vertical diffusivity in the boundary layer."
          endif;enddo
          ! Try to read the diffusivity value
          read(linebuffer080,*,iostat=ios)KvBL_model_ID,diffusivity_BL
          if(ios.ne.0)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: ",&
                    'Constant diffusivity in BL specified, but diffusivity value not provided.'
            endif;enddo
            stop 1
          else
            ! Error-checking diffusivity
            if(diffusivity_BL.lt.0.0_ip)then
              do io=1,2;if(VB(io).le.verbosity_error)then
                write(errlog(io),*)"ERROR: diffusivity_BL must be >0"
              endif;enddo
              stop 1
            endif
          endif
        elseif(KvBL_model_ID.eq.2)then
          ! Boundary Layer is turned off and vertical diffusivity will just be from the free-air equation
          useBoundaryLayer = .false.
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using free-air expression for vertical diffusivity throughout domain."
          endif;enddo
        elseif(KvBL_model_ID.eq.3)then
          ! Model from Troen and Mahrt, 1973
          phi_prefac = 1.0_ip
          phi_alpha = -0.33333_ip   ! Exponent in unstable term
          phi_beta  =  4.7_ip       ! Coefficient in stable term (pretty much always 4.7->5.2
          phi_gamma = -7.0_ip       ! Coefficient in unstable term
          PBL_exp_int = 2
          PBL_exp     = real(PBL_exp_int,kind=ip)
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using boundary layer vertical diffusivity as outlined by Troen and Mahrt (1973)."
          endif;enddo
        elseif(KvBL_model_ID.eq.4)then
          ! Model from Ulke (2000)
          ! Try for KvBL_MomHeat specifying Momentum vs Heat
          read(linebuffer080,*,iostat=ios2)KvBL_model_ID,KvBL_MomHeat
          if(ios2.eq.0)then
            if(KvBL_MomHeat.ne.2)then
              ! 2 specifies heat version; if it is anything but, then default to 1
              KvBL_MomHeat = 1
            endif
          else
            KvBL_MomHeat = 1
          endif
          if(KvBL_MomHeat.eq.1)then
            ! Momentum values
            phi_prefac = 1.0_ip
            phi_alpha = -0.25_ip  ! Exponent in unstable term
            phi_beta  =  6.9_ip   ! Coefficient in stable term (pretty much always 4.7->5.2
            phi_gamma = -22.0_ip  ! Coefficient in unstable term
            PBL_exp_int = 1
            PBL_exp     = real(PBL_exp_int,kind=ip)
          else
            ! Heat values
            phi_prefac = 1.0_ip
            phi_alpha = -0.5_ip   ! Exponent in unstable term
            phi_beta  =  9.2_ip   ! Coefficient in stable term (pretty much always 4.7->5.2
            phi_gamma = -13.0_ip  ! Coefficient in unstable term
            PBL_exp_int = 1
            PBL_exp     = real(PBL_exp_int,kind=ip)
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using boundary layer vertical diffusivity as outlined by Ulke (2000)."
          endif;enddo
        elseif(KvBL_model_ID.eq.5)then
          ! Model from Shir / Businger,Ayer outlined in Seinfeld and Pandis Eq 18.125
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using boundary layer vertical diffusivity as outlined by in Seinfeld and Pandis."
          endif;enddo
          ! Try for KvBL_MomHeat specifying Momentum vs Heat
          read(linebuffer080,*,iostat=ios2)KvBL_model_ID,KvBL_MomHeat
          if(ios2.eq.0)then
            if(KvBL_MomHeat.ne.2)then
              ! 2 specifies heat version; if it is anything but, then default to 1
              KvBL_MomHeat = 1
            endif
          else
            KvBL_MomHeat = 1
          endif
          if(KvBL_MomHeat.eq.1)then
            ! Momentum values
            phi_prefac = 1.0_ip   
            phi_alpha = -0.25_ip  ! Exponent in unstable term
            phi_beta  =  5.0_ip   ! Coefficient in stable term (pretty much always 4.7->5.2
            phi_gamma = -16.0_ip  ! Coefficient in unstable term
            PBL_exp_int = 0
          else
            ! Heat values
            ! Momentum values
            phi_prefac = 0.74_ip
            phi_alpha = -0.5_ip   ! Exponent in unstable term
            phi_beta  =  4.7_ip   ! Coefficient in stable term (pretty much always 4.7->5.2
            phi_gamma = -9.0_ip   ! Coefficient in unstable term
            PBL_exp_int = 0
          endif
        elseif(KvBL_model_ID.eq.0)then
          ! Coefficients of phi model specified in control file
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using boundary layer vertical diffusivity with generalized phi."
          endif;enddo
          ! Try to read all the parameters
          read(linebuffer080,*,iostat=ios)KvBL_model_ID,phi_prefac,phi_alpha,phi_beta,phi_gamma,PBL_exp_int
          if(ios.ne.0)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: ",&
                    'Custom BL model selected (KvBL_model_ID=0), but expected parameters could not be read.'
              write(errlog(io),*)"Expected format = "
              write(errlog(io),*)" KvBL_model_ID phi_prefac phi_alpha phi_beta phi_gamma PBL_exp_int"
            endif;enddo
            stop 1
          endif
          ! Some error-checking on values
          if(phi_prefac.le.0.0_ip.or.phi_prefac.gt.1.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: phi_prefac must be in range (0,1]"
            endif;enddo
            stop 1
          endif
          if(phi_alpha.ge.0.0_ip.or.phi_alpha.lt.-1.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: phi_alpha must be in range [-1,0)"
            endif;enddo
            stop 1
          endif
          if(phi_beta.lt.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: phi_beta must be >0"
            endif;enddo
            stop 1
          endif
          if(phi_gamma.gt.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: phi_gamma must be <0"
            endif;enddo
            stop 1
          endif
          if(PBL_exp_int.lt.0)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: PBL_exp_int must be a positive integer"
              write(errlog(io),*)"         or = 0 to signify useing an exponential taper/"
            endif;enddo
            stop 1
          endif
        else
          KvBL_model_ID = 2
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using boundary layer vertical diffusivity as outlined by Troen and Mahrt (1973)."
          endif;enddo
        endif
        if(KvBL_model_ID.eq.0.or.&
           KvBL_model_ID.eq.3.or.&
           KvBL_model_ID.eq.4.or.&
           KvBL_model_ID.eq.5)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using Phi with the following paramters:"
            write(outlog(io),*)"      phi_prefac  = ",real(phi_prefac,kind=4)
            write(outlog(io),*)"      phi_alpha   = ",real(phi_alpha,kind=4)
            write(outlog(io),*)"      phi_beta    = ",real(phi_beta,kind=4)
            write(outlog(io),*)"      phi_gamma   = ",real(phi_gamma,kind=4)
            write(outlog(io),*)"      PBL_exp_int = ",PBL_exp_int
          endif;enddo
        endif
        ! Free-air model: essentially, which scaling factor F(Ri) to use.
        read(10,'(a80)',iostat=ios,err=2010)linebuffer080
        ! Line 6:
        var_User_charlines_VarDiff(6) = trim(adjustl(linebuffer080))
        read(linebuffer080,*,iostat=ios)KvFA_model_ID
        if(KvFA_model_ID.eq.1)then
          ! Diffusivity is constant in the FA; read the value
          read(linebuffer080,*,iostat=ios)KvFA_model_ID,diffusivity_vert
          if(ios.eq.0)then
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"    Using constant vertical diffusivity above the boundary layer."
            endif;enddo
          else
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: Constant diffusivity specified in free-atmosphere, but no"
              write(errlog(io),*)"       diffusivity given."
            endif;enddo
            stop 1
          endif
          ! Error-checking diffusivity
          if(diffusivity_vert.lt.0.0_ip)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: diffusivity_FA must be >0"
            endif;enddo
            stop 1
          endif
        elseif(KvFA_model_ID.eq.2)then
          ! Mixing length model with Fc from Louis (1979)
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using free-air vertical diffusivity with stability function from Louis (1979)."
          endif;enddo
        elseif(KvFA_model_ID.eq.3)then
          ! Mixing length model with Fc from Stull (1988)
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using free-air vertical diffusivity with stability function from Stull (1988)"
          endif;enddo
        elseif(KvFA_model_ID.eq.4)then
          ! Mixing length model with Fc from Betts (1996)
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using free-air vertical diffusivity with stability function from Betts et al (1996)"
          endif;enddo
        elseif(KvFA_model_ID.eq.5)then
          ! Mixing length model with Fc from Hong (1996)
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using free-air vertical diffusivity with stability function from Hong and Pan (1996)"
          endif;enddo
        elseif(KvFA_model_ID.eq.6)then
          ! Mixing length model with Fc from Collins et al. (2004)
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using free-air vertical diffusivity with stability function from Collins et al (2004)"
          endif;enddo
        else
          KvFA_model_ID = 4
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    WARNING: free-air model not read correctly. Setting to default."
            write(outlog(io),*)"    Using free-air vertical diffusivity with stability function from Betts et al (1996)"
          endif;enddo
        endif

        if(useBoundaryLayer)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"      The planetary boundary layer will be identified, if possible, first by"
            write(outlog(io),*)"      trying to read PBLH from the Met files. If that is unavailable, then"
            write(outlog(io),*)"      Ash3d will search for a low-level temperature inversion, inspect Ri(z)"
            write(outlog(io),*)"      relative to Ri_crit, and calculate the Eckman layer thickness from the"
            write(outlog(io),*)"      latitude and friction velocity. If the atmospheric data is too coarse"
            write(outlog(io),*)"      to determine the boundary layer, the free-air mixing-length will be used"
            write(outlog(io),*)"      for vertical diffusivity calculations throughout the domain."
          endif;enddo
        endif

      endif

      if (useVarDiffH.or.useVarDiffV) then
        ! Check if we're using variable diffusivity, then get the constants
        read(10,'(a80)',iostat=ios,err=2010)linebuffer080
        ! Line 7:
        var_User_charlines_VarDiff(7) = trim(adjustl(linebuffer080))
        read(linebuffer080,*,iostat=ioerr) vonKarman
        read(10,'(a80)',iostat=ios,err=2010)linebuffer080
        ! Line 8:
        var_User_charlines_VarDiff(8) = trim(adjustl(linebuffer080))
        read(linebuffer080,*,iostat=ioerr) LambdaC
        read(10,'(a80)',iostat=ios,err=2010)linebuffer080
        ! Line 9:
        var_User_charlines_VarDiff(9) = trim(adjustl(linebuffer080))
        read(linebuffer080,*,iostat=ioerr) RI_CRIT

        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"von Karman constant = ",real(vonKarman,kind=4)
          write(outlog(io),*)"      Mixing Length = ",real(LambdaC,kind=4)
          write(outlog(io),*)"        Critical Ri = ",real(RI_CRIT,kind=4)
        endif;enddo

        ! We will want to reuse velocities on the MetP grid for this module
        MR_Save_Velocities = .true.
      endif

      ! Now set up output variable options
      if(use_Output_Vars_VarDiff.and.useVarDiffH)then
        nvar_User3d_XYZ_VarDiff = nvar_User3d_XYZ_VarDiff + 1  ! for Kh
      endif

      if(use_Output_Vars_VarDiff.and.useVarDiffV)then
        nvar_User2d_XY_VarDiff  = nvar_User2d_XY_VarDiff  + 1  ! for Pbl
        nvar_User2d_XY_VarDiff  = nvar_User2d_XY_VarDiff  + 1  ! for U*
        nvar_User3d_XYZ_VarDiff = nvar_User3d_XYZ_VarDiff + 1  ! for Kv
        nvar_User3d_XYZ_VarDiff = nvar_User3d_XYZ_VarDiff + 1  ! for Ri
      endif

2010  continue
      close(10)

      return

1900  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error: cannot find input file: ',infile
        write(errlog(io),*)  'Program stopped'
      endif;enddo
      stop 1

2011  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*) 'Error reading whether to use variable diffusivity.'
        write(errlog(io),*) 'Answer must be yes or no.'
        write(errlog(io),*) 'You gave:',linebuffer080
        write(errlog(io),*) 'Program stopped'
      endif;enddo
      stop 1

      end subroutine input_data_VarDiff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Allocate_VarDiff_Met
!
!  Called from: Ash3d.F90
!  Arguments:
!    none
!
!  This subroutine allocates all the variables needed for calculating Kv and Kh
!  on the MetP grid. This includes some variables that are part of the Atmosphere
!  module, and might have been allocated elsewhere, but are needed here.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Allocate_VarDiff_Met

      use global_param,  only : &
         PI,useVarDiffH,useVarDiffV

      use io_data,       only : &
         nvar_User2d_static_XY,nvar_User3d_XYGs,nvar_User2d_XY,&
         nvar_User4d_XYZGs,nvar_User3d_XYZ,nvar_User_charlines

      use mesh,          only : &
         nxmax,nymax,nzmax

      use Atmosphere,    only : &
         Ri_meso_last_step_MetP_sp,     Ri_meso_next_step_MetP_sp,       &
         PBLH_meso_last_step_Met_sp,    PBLH_meso_next_step_Met_sp,      &
         L_MonOb_meso_last_step_Met_sp, L_MonOb_meso_next_step_Met_sp,   &
         FricVel_meso_last_step_Met_sp, FricVel_meso_next_step_Met_sp

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet

      implicit none

      integer :: i

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"--------------------------------------------------"
        write(outlog(io),*)"---------- ALLOCATE_VARDIFF_MET ------------------"
        write(outlog(io),*)"--------------------------------------------------"
      endif;enddo

#ifdef USEPOINTERS
      if(useVarDiffH)then
        if(.not.associated(du_dx_MetP_sp))   &
                  allocate(du_dx_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.associated(du_dy_MetP_sp))   &
                  allocate(du_dy_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.associated(dv_dx_MetP_sp))   &
                  allocate(dv_dx_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.associated(dv_dx_MetP_sp))   &
                  allocate(dv_dy_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.associated(Khz_meso_last_step_sp))   &
                  allocate(Khz_meso_last_step_sp(nxmax,nymax,nzmax))
        if(.not.associated(Khz_meso_next_step_sp))   &
                  allocate(Khz_meso_next_step_sp(nxmax,nymax,nzmax))
      endif
      if(useVarDiffV)then
        if(.not.associated(dV_dz_MetP_sp))   &
                  allocate(dV_dz_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.associated(SurfRoughLen_Met_sp)) &
                  allocate(SurfRoughLen_Met_sp(nx_submet,ny_submet))
        if(.not.associated(Ri_meso_last_step_MetP_sp))  &
                  allocate(Ri_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.associated(Ri_meso_next_step_MetP_sp))  &
                  allocate(Ri_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.associated(Khz_meso_last_step_MetP_sp)) &
                  allocate(Khz_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.associated(Khz_meso_next_step_MetP_sp)) &
                  allocate(Khz_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.associated(Kv_meso_last_step_MetP_sp))  &
                  allocate(Kv_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.associated(Kv_meso_next_step_MetP_sp))  &
                  allocate(Kv_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.associated(PBLH_meso_last_step_Met_sp)) &
                  allocate(PBLH_meso_last_step_Met_sp(nx_submet,ny_submet))
        if(.not.associated(PBLH_meso_next_step_Met_sp)) &
                  allocate(PBLH_meso_next_step_Met_sp(nx_submet,ny_submet))
        if(.not.associated(L_MonOb_meso_last_step_Met_sp)) &
                  allocate(L_MonOb_meso_last_step_Met_sp(nx_submet,ny_submet))
        if(.not.associated(L_MonOb_meso_next_step_Met_sp)) &
                  allocate(L_MonOb_meso_next_step_Met_sp(nx_submet,ny_submet))
        if(.not.associated(FricVel_meso_last_step_Met_sp)) &
                  allocate(FricVel_meso_last_step_Met_sp(nx_submet,ny_submet))
        if(.not.associated(FricVel_meso_next_step_Met_sp)) &
                  allocate(FricVel_meso_next_step_Met_sp(nx_submet,ny_submet))
        if(.not.associated(Kv_meso_last_step_sp))    &
                  allocate(Kv_meso_last_step_sp(nxmax,nymax,nzmax))
        if(.not.associated(Kv_meso_next_step_sp))    &
                  allocate(Kv_meso_next_step_sp(nxmax,nymax,nzmax))
      endif

      if(.not.associated(vx_meso_last_step_MetP_sp))   &
              associated allocate(vx_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if(.not.associated(vy_meso_last_step_MetP_sp))   &
              associated allocate(vy_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if(.not.associated(vx_meso_next_step_MetP_sp))   &
              associated allocate(vx_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if(.not.associated(vy_meso_next_step_MetP_sp))   &
               allocate(vy_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
#else
      if(useVarDiffH)then
        if(.not.allocated(du_dx_MetP_sp))   &
                 allocate(du_dx_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.allocated(du_dy_MetP_sp))   &
                 allocate(du_dy_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.allocated(dv_dx_MetP_sp))   &
                 allocate(dv_dx_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.allocated(dv_dx_MetP_sp))   &
                 allocate(dv_dy_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.allocated(dV_dz_MetP_sp))   &
                 allocate(dV_dz_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.allocated(Khz_meso_last_step_sp))   &
                 allocate(Khz_meso_last_step_sp(nxmax,nymax,nzmax))
        if(.not.allocated(Khz_meso_next_step_sp))   &
                 allocate(Khz_meso_next_step_sp(nxmax,nymax,nzmax))
      endif
      if(useVarDiffV)then
        if(.not.allocated(SurfRoughLen_Met_sp)) &
                 allocate(SurfRoughLen_Met_sp(nx_submet,ny_submet))
        if(.not.allocated(Ri_meso_last_step_MetP_sp))  &
                 allocate(Ri_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.allocated(Ri_meso_next_step_MetP_sp))  &
                 allocate(Ri_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.allocated(Khz_meso_last_step_MetP_sp)) &
                 allocate(Khz_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.allocated(Khz_meso_next_step_MetP_sp)) &
                 allocate(Khz_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.allocated(Kv_meso_last_step_MetP_sp))  &
                 allocate(Kv_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.allocated(Kv_meso_next_step_MetP_sp))  &
                 allocate(Kv_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
        if(.not.allocated(PBLH_meso_last_step_Met_sp)) &
                 allocate(PBLH_meso_last_step_Met_sp(nx_submet,ny_submet))
        if(.not.allocated(PBLH_meso_next_step_Met_sp)) &
                 allocate(PBLH_meso_next_step_Met_sp(nx_submet,ny_submet))
        if(.not.allocated(L_MonOb_meso_last_step_Met_sp)) &
                 allocate(L_MonOb_meso_last_step_Met_sp(nx_submet,ny_submet))
        if(.not.allocated(L_MonOb_meso_next_step_Met_sp)) &
                 allocate(L_MonOb_meso_next_step_Met_sp(nx_submet,ny_submet))
        if(.not.allocated(FricVel_meso_last_step_Met_sp)) &
                 allocate(FricVel_meso_last_step_Met_sp(nx_submet,ny_submet))
        if(.not.allocated(FricVel_meso_next_step_Met_sp)) &
                 allocate(FricVel_meso_next_step_Met_sp(nx_submet,ny_submet))
        if(.not.allocated(Kv_meso_last_step_sp))    &
                 allocate(Kv_meso_last_step_sp(nxmax,nymax,nzmax))
        if(.not.allocated(Kv_meso_next_step_sp))    &
                 allocate(Kv_meso_next_step_sp(nxmax,nymax,nzmax))
      endif

      if(.not.allocated(vx_meso_last_step_MetP_sp))   &
               allocate(vx_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if(.not.allocated(vy_meso_last_step_MetP_sp))   &
               allocate(vy_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if(.not.allocated(vx_meso_next_step_MetP_sp))   &
               allocate(vx_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      if(.not.allocated(vy_meso_next_step_MetP_sp))   &
               allocate(vy_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
#endif

        ! Precalculate the LES term
      LES_L2ScaleCoeff = (KH_SmagC*KH_SmagC/(PI*PI))

      ! Set the start indecies
      indx_User2d_static_XY_VarDiff = nvar_User2d_static_XY
      indx_User2d_XY_VarDiff        = nvar_User2d_XY
      indx_User3d_XYGs_VarDiff      = nvar_User3d_XYGs
      indx_User3d_XYZ_VarDiff       = nvar_User3d_XYZ
      indx_User4d_XYZGs_VarDiff     = nvar_User4d_XYZGs
      indx_User_charlines_VarDiff   = nvar_User_charlines

      ! Allocate output variables if needed
      if(.not.allocated(temp_2d_name_VarDiff))    &
               allocate(temp_2d_name_VarDiff(nvar_User2d_XY_VarDiff))
      if(.not.allocated(temp_2d_unit_VarDiff))    &
               allocate(temp_2d_unit_VarDiff(nvar_User2d_XY_VarDiff))
      if(.not.allocated(temp_2d_lname_VarDiff))   &
               allocate(temp_2d_lname_VarDiff(nvar_User2d_XY_VarDiff))
      if(.not.allocated(temp_2d_MissVal_VarDiff)) &
               allocate(temp_2d_MissVal_VarDiff(nvar_User2d_XY_VarDiff))
      if(.not.allocated(temp_2d_FillVal_VarDiff)) &
               allocate(temp_2d_FillVal_VarDiff(nvar_User2d_XY_VarDiff))

      if(.not.allocated(temp_3d_name_VarDiff))    &
               allocate(temp_3d_name_VarDiff(nvar_User3d_XYZ_VarDiff))
      if(.not.allocated(temp_3d_unit_VarDiff))    &
               allocate(temp_3d_unit_VarDiff(nvar_User3d_XYZ_VarDiff))
      if(.not.allocated(temp_3d_lname_VarDiff))   &
               allocate(temp_3d_lname_VarDiff(nvar_User3d_XYZ_VarDiff))
      if(.not.allocated(temp_3d_MissVal_VarDiff)) &
               allocate(temp_3d_MissVal_VarDiff(nvar_User3d_XYZ_VarDiff))
      if(.not.allocated(temp_3d_FillVal_VarDiff)) &
               allocate(temp_3d_FillVal_VarDiff(nvar_User3d_XYZ_VarDiff))

      i = 0 
      if(use_Output_Vars_VarDiff.and.useVarDiffH)then
        i = 1
        temp_3d_name_VarDiff(i) = "Kh"
        temp_3d_lname_VarDiff(i) = "Diffusivity_Horizontal"
        temp_3d_unit_VarDiff(i) = "m2/s"
        temp_3d_MissVal_VarDiff(i) = -9999.0_op
        temp_3d_FillVal_VarDiff(i) = -9999.0_op
      endif

      if(use_Output_Vars_VarDiff.and.useVarDiffV)then
        temp_2d_name_VarDiff(1) = "PBLH"
        temp_2d_lname_VarDiff(1) = "Planetary Boundary Layer H"
        temp_2d_unit_VarDiff(1) = "km"
        temp_2d_MissVal_VarDiff(1) = -9999.0_op
        temp_2d_FillVal_VarDiff(1) = -9999.0_op

        temp_2d_name_VarDiff(2) = "Ust"
        temp_2d_lname_VarDiff(2) = "Friction Velocity"
        temp_2d_unit_VarDiff(2) = "m/s"
        temp_2d_MissVal_VarDiff(2) = -9999.0_op
        temp_2d_FillVal_VarDiff(2) = -9999.0_op

        i = i + 1
        temp_3d_name_VarDiff(i) = "Kv"
        temp_3d_lname_VarDiff(i) = "Diffusivity_Vertical"
        temp_3d_unit_VarDiff(i) = "m2/s"
        temp_3d_MissVal_VarDiff(i) = -9999.0_op
        temp_3d_FillVal_VarDiff(i) = -9999.0_op

        i = i + 1
        temp_3d_name_VarDiff(i) = "Ri"
        temp_3d_lname_VarDiff(i) = "Gradient_Richardson_Number"
        temp_3d_unit_VarDiff(i) = "none"
        temp_3d_MissVal_VarDiff(i) = -9999.0_op
        temp_3d_FillVal_VarDiff(i) = -9999.0_op
      endif

      nvar_User2d_static_XY = nvar_User2d_static_XY + nvar_User2d_static_XY_VarDiff
      nvar_User2d_XY        = nvar_User2d_XY        + nvar_User2d_XY_VarDiff
      nvar_User3d_XYGs      = nvar_User3d_XYGs      + nvar_User3d_XYGs_VarDiff
      nvar_User3d_XYZ       = nvar_User3d_XYZ       + nvar_User3d_XYZ_VarDiff
      nvar_User4d_XYZGs     = nvar_User4d_XYZGs     + nvar_User4d_XYZGs_VarDiff
      nvar_User_charlines   = nvar_User_charlines   + nvar_User_charlines_VarDiff

      end subroutine Allocate_VarDiff_Met

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Prep_output_VarDiff
!
!  Called from: Ash3d.F90
!  Arguments:
!    none
!
!  This subroutine fills the module output variables. The output variables that
!  are part of the atmospheric stability metrics (Ustar, PBLH, Ri) are only on 
!  the Met steps since they are only needed on those intervals to get Kv, Kh.
!  Kv and Kh are interpolated to each time-step so will vary smoothly.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Prep_output_VarDiff

      use global_param,  only : &
         useVarDiffH,KM_2_M,M2PS_2_KM2PHR

      use mesh,          only : &
         nxmax,nymax,nzmax !,lon_cc_pd,lat_cc_pd

      use Diffusion,     only : &
         kx,kz

      use Output_Vars,   only : &
         var_User2d_XY_name,var_User2d_XY_unit,var_User2d_XY_lname,            &
         var_User2d_XY_MissVal,var_User2d_XY_FillVal,var_User2d_XY,            &
         var_User3d_XYZ_name,var_User3d_XYZ_unit,var_User3d_XYZ_lname,         &
         var_User3d_XYZ_MissVal,var_User3d_XYZ_FillVal,var_User3d_XYZ,         &
         var_User_charlines

      use Atmosphere,    only : &
         FricVel_meso_last_step_Met_sp, &
         PBLH_meso_last_step_Met_sp, &
         Ri_meso_last_step_MetP_sp, &
         SolZen_meso_next_step_Met_sp, &
           Set_SolarZenith

      use MetReader,     only : &
         MR_dum2d_met,MR_dum2d_comp,MR_dum3d_compH,MR_dum3d_metP,MR_iMetStep_Now,&
           MR_Regrid_MetP_to_CompH,&
           MR_Regrid_Met2d_to_Comp2D

      implicit none

      integer :: i,ii,indx

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Prep_output_VarDiff"
      endif;enddo

      ! Might have to build in some logic for Kh vs Kv

      do i=1,nvar_User2d_XY_VarDiff
        ! We could put a conditional block here, but we would not enter this
        ! do loop if we didn't plan to write out Planetary Boundary Layer Height
        indx = indx_User2d_XY_VarDiff+i
        var_User2d_XY_name(indx)   = temp_2d_name_VarDiff(i)
        var_User2d_XY_unit(indx)   = temp_2d_unit_VarDiff(i)
        var_User2d_XY_lname(indx)  = temp_2d_lname_VarDiff(i)
        var_User2d_XY_MissVal(indx)= temp_2d_MissVal_VarDiff(i)
        var_User2d_XY_FillVal(indx)= temp_2d_FillVal_VarDiff(i)
        if(i.eq.1)then
           ! Now resample onto computational grid
          MR_dum2d_met = PBLH_meso_last_step_Met_sp/real(KM_2_M,kind=sp)
          call MR_Regrid_Met2d_to_Comp2D
          var_User2d_XY(1:nxmax,1:nymax,indx) = MR_dum2d_comp(1:nxmax,1:nymax)

          ! Uncomment these lines to export Solar Zenith to the PBLH slot
          !call Set_SolarZenith(MR_MetStep_Hour_since_baseyear(MR_iMetStep_Now))
          !MR_dum2d_met = SolZen_meso_next_step_Met_sp
          !call MR_Regrid_Met2d_to_Comp2D
          !var_User2d_XY(1:nxmax,1:nymax,indx) = MR_dum2d_comp(1:nxmax,1:nymax)

        elseif(i.eq.2)then
           ! Now resample onto computational grid
          MR_dum2d_met = FricVel_meso_last_step_Met_sp
          call MR_Regrid_Met2d_to_Comp2D
          var_User2d_XY(1:nxmax,1:nymax,indx) = MR_dum2d_comp(1:nxmax,1:nymax)
        endif
      enddo

      do i=1,nvar_User3d_XYZ_VarDiff
        indx = indx_User3d_XYZ_VarDiff+i
        var_User3d_XYZ_name(indx)   = temp_3d_name_VarDiff(i)
        var_User3d_XYZ_unit(indx)   = temp_3d_unit_VarDiff(i)
        var_User3d_XYZ_lname(indx)  = temp_3d_lname_VarDiff(i)
        var_User3d_XYZ_MissVal(indx)= temp_3d_MissVal_VarDiff(i)
        var_User3d_XYZ_FillVal(indx)= temp_3d_FillVal_VarDiff(i)
        if(use_Output_Vars_VarDiff.and.useVarDiffH)then
          ii = 1
        else
          ii = 0
        endif
        if(i.eq.ii  )then
          ! Horizontal diffusivity is already on the comp grid
          ! This branch is unused if useVarDiffH = .false. since ii=0
          ! Note that the native units are km2/hr, but we need to convert to m2/s
          var_User3d_XYZ(1:nxmax,1:nymax,1:nzmax,indx) = kx(1:nxmax,1:nymax,1:nzmax)/M2PS_2_KM2PHR
        endif
        if(i.eq.ii+1)then
          ! Vertical diffusivity is already on the comp grid
          !var_User3d_XYZ(1:nxmax,1:nymax,1:nzmax,indx) = kz(1:nxmax,1:nymax,1:nzmax)*KM_2_M*KM_2_M/HR_2_S
          var_User3d_XYZ(1:nxmax,1:nymax,1:nzmax,indx) = kz(1:nxmax,1:nymax,1:nzmax)/M2PS_2_KM2PHR
        endif
        if(i.eq.ii+2)then
           ! Ri just exists on the MetP grid for output
           ! Now resample onto computational grid
          MR_dum3d_metP = Ri_meso_last_step_MetP_sp
          call MR_Regrid_MetP_to_CompH(MR_iMetStep_Now)
          var_User3d_XYZ(1:nxmax,1:nymax,1:nzmax,indx) = MR_dum3d_compH(1:nxmax,1:nymax,1:nzmax)
        endif
      enddo

      do i=1,nvar_User_charlines_VarDiff
        indx = indx_User_charlines_VarDiff+i
        var_User_charlines(indx) = var_User_charlines_VarDiff(i)
      enddo

      end subroutine Prep_output_VarDiff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Deallocate_VarDiff_Met
!
!  Called from: Ash3d.F90
!  Arguments:
!    none
!
!  This subroutine deallocates the variables allocated in Allocate_VarDiff_Met
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Deallocate_VarDiff_Met

      use Atmosphere,    only : &
         Ri_meso_last_step_MetP_sp,Ri_meso_next_step_MetP_sp, &
         PBLH_meso_last_step_Met_sp,PBLH_meso_next_step_Met_sp,&
         L_MonOb_meso_last_step_Met_sp,L_MonOb_meso_next_step_Met_sp,&
         FricVel_meso_last_step_Met_sp,FricVel_meso_next_step_Met_sp

      implicit none

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Deallocate_VarDiff_Met"
      endif;enddo

!      if(allocated(dVel_dz_MetP_sp))               deallocate(dVel_dz_MetP_sp)
      if(allocated(du_dx_MetP_sp))                 deallocate(du_dx_MetP_sp)
      if(allocated(du_dy_MetP_sp))                 deallocate(du_dy_MetP_sp)
      if(allocated(dv_dx_MetP_sp))                 deallocate(dv_dx_MetP_sp)
      if(allocated(dv_dy_MetP_sp))                 deallocate(dv_dy_MetP_sp)
      if(allocated(dV_dz_MetP_sp))                 deallocate(dV_dz_MetP_sp)
      if(allocated(SurfRoughLen_Met_sp))           deallocate(SurfRoughLen_Met_sp)
      if(allocated(Ri_meso_last_step_MetP_sp))     deallocate(Ri_meso_last_step_MetP_sp)
      if(allocated(Ri_meso_next_step_MetP_sp))     deallocate(Ri_meso_next_step_MetP_sp)
      if(allocated(Khz_meso_last_step_MetP_sp))    deallocate(Khz_meso_last_step_MetP_sp)
      if(allocated(Khz_meso_next_step_MetP_sp))    deallocate(Khz_meso_next_step_MetP_sp)
      if(allocated(Kv_meso_last_step_MetP_sp))     deallocate(Kv_meso_last_step_MetP_sp)
      if(allocated(Kv_meso_next_step_MetP_sp))     deallocate(Kv_meso_next_step_MetP_sp)
      if(allocated(PBLH_meso_last_step_Met_sp))    deallocate(PBLH_meso_last_step_Met_sp)
      if(allocated(PBLH_meso_next_step_Met_sp))    deallocate(PBLH_meso_next_step_Met_sp)
      if(allocated(L_MonOb_meso_last_step_Met_sp)) deallocate(L_MonOb_meso_last_step_Met_sp)
      if(allocated(L_MonOb_meso_next_step_Met_sp)) deallocate(L_MonOb_meso_next_step_Met_sp)
      if(allocated(FricVel_meso_last_step_Met_sp)) deallocate(FricVel_meso_last_step_Met_sp)
      if(allocated(FricVel_meso_next_step_Met_sp)) deallocate(FricVel_meso_next_step_Met_sp)

      if(allocated(Khz_meso_last_step_sp))         deallocate(Khz_meso_last_step_sp)
      if(allocated(Khz_meso_next_step_sp))         deallocate(Khz_meso_next_step_sp)
      if(allocated(Kv_meso_last_step_sp))          deallocate(Kv_meso_last_step_sp)
      if(allocated(Kv_meso_next_step_sp))          deallocate(Kv_meso_next_step_sp)

      if(allocated(vx_meso_last_step_MetP_sp))     deallocate(vx_meso_last_step_MetP_sp)
      if(allocated(vy_meso_last_step_MetP_sp))     deallocate(vy_meso_last_step_MetP_sp)
      if(allocated(vx_meso_next_step_MetP_sp))     deallocate(vx_meso_next_step_MetP_sp)
      if(allocated(vy_meso_next_step_MetP_sp))     deallocate(vy_meso_next_step_MetP_sp)

      end subroutine Deallocate_VarDiff_Met

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Eddy_diff
!
!  Called from: Set_VarDiffH_Meso
!  Arguments:
!    none
!
!  This subroutine calculates the horizontal diffusivity based on the specified
!  model at the 'next' Met step.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Eddy_diff

      use global_param,  only : &
         HR_2_S,M2PS_2_KM2PHR

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet,MR_sigma_nz_submet

      implicit none

      integer :: i,j,k

      real(kind=ip) :: E11,E12,E21,E22
      real(kind=ip) :: D2_tension,D2_strain
      real(kind=ip) :: LES_TimeScale
      real(kind=ip) :: LES_LengthScale2

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Eddy_diff"
      endif;enddo

      if(Kh_model_ID.eq.1)then
        Khz_meso_next_step_MetP_sp(:,:,:) = real(diffusivity_vert*M2PS_2_KM2PHR,kind=sp)
      else
        do i=1,nx_submet
          do j=1,ny_submet
            do k=1,np_fullmet
  
        ! Smagorinsky LES horizontal eddy diffusivity is proportional
        ! to sqrt((E12+E21)^2 + (E11-E22)^2) where E is the velocity gradient
        ! tensor (just in x and y)
        ! This is following the description in
        ! Griffies and Hallberg, MWR, 2000 doi:10.1175/1520-0493(2000)128<2935:BFWASL>2.0.CO;2
  
          ! spatial derivatives of velocity (in 1/s)
        E11 = du_dx_MetP_sp(i,j,k)
        E12 = du_dy_MetP_sp(i,j,k)
        E21 = dv_dx_MetP_sp(i,j,k)
        E22 = dv_dy_MetP_sp(i,j,k)
  
        D2_strain  = (E12+E21)**2.0_ip
        if(Kh_model_ID.eq.2)then
          D2_tension = (E11-E22)**2.0_ip          ! Smagorinsky (1963, 1993)
        elseif(Kh_model_ID.eq.3)then
          D2_tension = 0.5_ip*(E11*E11+E22*E22)   ! Pielke (1974)
        else
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)  'error: currently only Smagorinsky and Pielke models allowed.'
          endif;enddo
          stop 1
        endif
          ! in units of 1/s
        LES_TimeScale = sqrt(D2_tension+D2_strain)
          ! in units of 1/hr
        LES_TimeScale = LES_TimeScale * HR_2_S
          ! length scale^2 in km^2
        LES_LengthScale2 = min(MR_sigma_nz_submet(i,j),MAX_LES_LengthScale2)
  
          ! Diffusivity in km2/hr
        Khz_meso_next_step_MetP_sp(i,j,k) = real(KH_SmagC*LES_LengthScale2*LES_TimeScale,kind=sp)
              enddo !k
          enddo !j
        enddo !i
      endif

      return

      end subroutine Eddy_diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Calc_Vert_Diff
!
!  Called from: Set_VarDiffV_Meso
!  Arguments: last_or_next
!    none
!
!  This subroutine loops through all the nodes of the NWP subgrid and calculates
!  the vertical and horizontal diffusivities based on the atmospheric layer and
!  on the model specifications. This subroutine will work on either the next
!  or the last Met step.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Calc_Vert_Diff(last_or_next)

      use global_param,  only : &
         KM_2_M,HR_2_S,KM_2_M

      use Atmosphere,    only : &
         Ri_meso_last_step_MetP_sp,Ri_meso_next_step_MetP_sp, &
         PBLH_meso_last_step_Met_sp,PBLH_meso_next_step_Met_sp,&
         FricVel_meso_last_step_Met_sp,FricVel_meso_next_step_Met_sp

      use MetReader,     only : &
        nx_submet,ny_submet,np_fullmet,MR_geoH_metP_last,MR_geoH_metP_next

      implicit none

      integer, intent(in) :: last_or_next

      integer :: i,j,k
      real(kind=ip) :: z0
      real(kind=ip) :: Ri_col(np_fullmet)
      real(kind=ip) :: z_col(np_fullmet)
      real(kind=ip) :: dv_dz_col(np_fullmet)
      real(kind=ip) :: FricVel
      real(kind=ip) :: Kv_col(np_fullmet)
      real(kind=ip) :: PBLz
      real(kind=ip) :: Phi
      real(kind=ip) :: Fc
      real(kind=ip) :: PBL_profile_fac
      real(kind=ip) :: Kv_FreeAir, Kv_BL
      real(kind=ip) :: Lc
      real(kind=ip) :: zeta

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Calc_Vert_Diff"
      endif;enddo

      ! Even if we are not using a BL, we need Ri
      do i=1,nx_submet
        do j=1,ny_submet
          ! For the calculations, we need:
          !  PBLz, L_MonOb, FricVel, Ri, dv_dz, and z
          ! We need to set a minimum for z0 since some NWP files set this to 0
          ! Table 6.1 of Panofsky and Dutton giv 1.0e-4 for water.
          z0 = max(SurfRoughLen_Met_sp(i,j),1.0e-3_sp) ! We need to set a minimum for z0
          if(last_or_next.eq.0)then
              Ri_col(:) = real(Ri_meso_last_step_MetP_sp(i,j,:),kind=ip)    ! dimensionless
               z_col(:) = real(MR_geoH_metP_last(i,j,:),kind=ip)*KM_2_M     ! m
           dv_dz_col(:) = real(dV_dz_MetP_sp(i,j,:),kind=ip)                ! 1/s
                PBLz    = real(PBLH_meso_last_step_Met_sp(i,j),kind=ip)     ! m
!             L_MonOb    = real(L_MonOb_meso_last_step_Met_sp(i,j),kind=ip)  ! m
             FricVel    = real(FricVel_meso_last_step_Met_sp(i,j),kind=ip)  ! m/s
          else
              Ri_col(:) = real(Ri_meso_next_step_MetP_sp(i,j,:),kind=ip)    ! dimensionless
               z_col(:) = real(MR_geoH_metP_next(i,j,:),kind=ip)*KM_2_M     ! m
           dv_dz_col(:) = real(dV_dz_MetP_sp(i,j,:),kind=ip)                ! 1/s
                PBLz    = real(PBLH_meso_next_step_Met_sp(i,j),kind=ip)     ! m
!             L_MonOb    = real(L_MonOb_meso_next_step_Met_sp(i,j),kind=ip)  ! m
             FricVel    = real(FricVel_meso_next_step_Met_sp(i,j),kind=ip)  ! m/s
          endif

          do k = np_fullmet,1,-1
            ! Determine which form of Kv based on height relative to
            ! atmospheric boundary layer
            Kv_BL      = 0.0_ip
            Kv_FreeAir = 0.0_ip
            if(z_col(k).le.0.0_sp)then
                ! If point is at a negative gpm, then assign the kz from the
                ! node above
              Kv_BL = Kv_col(k+1)
            else
                ! In free atmosphere above the PBL, use Prandtl's mixing
                ! length theory for thermally stratified atmosphere.
                ! We calculate this term for all cases and update the PBL zone only if
                ! we exceed that calculated from mixing-length theory.
                !  First get mixing length scale
                !    There are several ways to parameterize the mixing
                !    length (Randerson, p155, 1984; Monin and Yaglom, v1,
                !    p409. Collins et al, NCAR TN-464, 2004, eq. 4.461)
              Lc = MixLen(real(z_col(k),kind=ip))
    
                ! calculate eq 8
                ! The Ri-term seems to zero out anything above the PBL
                ! since Ri is too high
              if(KvFA_model_ID.eq.1)then
                Kv_FreeAir = diffusivity_vert
              elseif(KvFA_model_ID.eq.2)then
                ! use scaling from Louis 1979
                Fc = Fc_Louis(Ri_col(k),z_col(k)/z0)
              elseif(KvFA_model_ID.eq.3)then
                ! use scaling from Jacobson (given earlier in Stull 1988)
                Fc = Fc_Jac(Ri_col(k))
              elseif(KvFA_model_ID.eq.4)then
                ! use scaling from Betts 1996
                Fc = Fc_Betts(Ri_col(k))
              elseif(KvFA_model_ID.eq.5)then
                ! use scaling from Hong (1996)
                Fc = Fc_Hong(Ri_col(k))
              elseif(KvFA_model_ID.eq.6)then
                ! use scaling from Collins (2004)
                Fc = Fc_Collins(Ri_col(k))
              endif
              Kv_FreeAir = Lc*Lc*abs(dv_dz_col(k))*Fc

              if(useBoundaryLayer)then
                if(KvBL_model_ID.eq.1.and.z_col(k).lt.PBLz)then
                  ! Kv constant and specified in BL
                  Kv_BL = diffusivity_BL
                elseif(KvBL_model_ID.eq.2)then
                  ! No BL model; set to zero and Kv will default to Free-Air
                  Kv_BL = 0.0_ip
                else
                  ! For this case, we need the stability function and the PBL taper
                  PBL_profile_fac = 0.0_ip
                  if(PBL_exp_int.eq.0.and.z_col(k).lt.3.0_ip*PBLz)then
                    ! PBL_exp_int = 0 indicates that we will be using the exponential taper
                    ! either because KvBL_model_ID = 5 (Businger-Arya) or 0 (custom) with exponential
                    ! selected. In either case, extend Kv zone above the PBLz (~3x)
                      PBL_profile_fac = exp(-2.0_ip*z_col(k)/PBLz);
                  elseif(z_col(k).lt.PBLz)then
                    ! All polynomial tapers are only valid below the PBLz

                    ! Within the PBL, use similarity theory
                    if(KvBL_model_ID.eq.3.or.KvBL_model_ID.eq.4.or.&
                      (KvBL_model_ID.eq.0.and.PBL_exp_int.gt.0))then
                      ! These are two polynomial models to taper profile for Kv between 0 and PBL
                      ! 3=> Troen-Mahrt: PBL_exp=2; quadratic taper
                      ! 4=> Ulke:        PBL_exp=1; linear taper
                      PBL_profile_fac = (1.0_sp-z_col(k)/PBLz)**PBL_exp
                    endif
                  endif

                  ! Now calculate stability function; get zeta using formula from Panofsky Eq. 6.7.1-2
                  if(Ri_col(k).lt.0.0_ip)then
                    zeta = Ri_col(k)
                  else
                    zeta = Ri_col(k)/(1.0_ip - 5.0_ip*Ri_col(k))
                  endif
                  Phi = Phi_WindShear_Similarity(zeta)
                  ! Kv from similarity theory (Eq. 8.48 of Jacobson)
                  Kv_BL = z_col(k)*vonKarman*FricVel*PBL_profile_fac/Phi
                endif
              endif ! useBoundaryLayer
            endif ! test on if this is a valid zone for Kv
   
            ! assign to array and convert from m2/s to km2/hr
            Kv_col(k) = max(Kv_BL,Kv_FreeAir) * HR_2_S/KM_2_M/KM_2_M
          enddo

          if(last_or_next.eq.0)then
            Kv_meso_last_step_MetP_sp(i,j,:) = real(Kv_col(:),kind=sp)
          else
            Kv_meso_next_step_MetP_sp(i,j,:) = real(Kv_col(:),kind=sp)
          endif

        enddo
      enddo

      return

      end subroutine Calc_Vert_Diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Set_VarDiffH_Meso
!
!  Called from: Ash3d.F90
!  Arguments: 
!    Load_MesoSteps,Interval_Frac
!
!  This subroutine calls MetReader routines to calculate velocity gradiants, then
!  calls Eddy_diff to get Kh at the Met steps. Finally, it interpolates Kh from the
!  Met steps to kx,ky on the time step in question.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Set_VarDiffH_Meso(Load_MesoSteps,Interval_Frac)

      use global_param,  only : &
         useDiffusion

      use mesh,          only : &
         nxmax,nymax,nzmax

      use Diffusion,     only : &
         kx,ky

      use MetReader,     only : &
         MR_dum3d_compH,MR_vx_metP_last,MR_dum3d_metP,MR_dum3d2_metP,MR_iMetStep_Now,&
         MR_vy_metP_last,MR_vy_metP_next,MR_vx_metP_next,&
           MR_DelMetP_Dx,&
           MR_DelMetP_Dy,&
           MR_Regrid_MetP_to_CompH

      implicit none

      logical      ,intent(in) :: Load_MesoSteps
      real(kind=dp),intent(in) :: Interval_Frac

      logical,save  :: first_time = .true.
      real(kind=sp) :: M_2_KM = 1.0e-3_sp

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Set_VarDiffH_Meso"
      endif;enddo

      if(.not.useDiffusion)useDiffusion = .true.

      if(Load_MesoSteps)then
        if(first_time)then
          ! Need to fill _last_step_sp
          !  First fill next step so that outside this 'first_time' loop, the
          !  'next' can be copied to the 'last'
          ! Load U winds on MetP
          vx_meso_next_step_MetP_sp = MR_vx_metP_last
          MR_dum3d_metP             = MR_vx_metP_last
            ! Now differentiate
            ! Note: velocities are in m/s, but the dx and dy are in km
            !       We want du_dx to be in 1/s
          call MR_DelMetP_Dx
          du_dx_MetP_sp = MR_dum3d2_metP * M_2_KM
          call MR_DelMetP_Dy
          du_dy_MetP_sp = MR_dum3d2_metP * M_2_KM

          ! Load V winds on MetP
          vy_meso_next_step_MetP_sp = MR_vy_metP_last
          MR_dum3d_metP             = MR_vy_metP_last
            ! Now differentiate
            ! Again, velocities are in m/s, but the dx and dy are in km
          call MR_DelMetP_Dx
          dv_dx_MetP_sp = MR_dum3d2_metP * M_2_KM
          call MR_DelMetP_Dy
          dv_dy_MetP_sp = MR_dum3d2_metP * M_2_KM
          call Eddy_diff  ! this sets Khz in km2/hr
           ! Now resample onto computational grid
          MR_dum3d_metP = Khz_meso_next_step_MetP_sp
          call MR_Regrid_MetP_to_CompH(MR_iMetStep_Now)
          Khz_meso_next_step_sp = MR_dum3d_compH
          first_time = .false.
        endif ! first_time
        Khz_meso_last_step_MetP_sp = Khz_meso_next_step_MetP_sp
        Khz_meso_last_step_sp      = Khz_meso_next_step_sp
        vx_meso_last_step_MetP_sp  = vx_meso_next_step_MetP_sp
        vy_meso_last_step_MetP_sp  = vy_meso_next_step_MetP_sp

        ! Need to fill _next_step_sp
        ! Load U winds on MetP
        !ivar = 2 ! U winds
        !call MR_Read_3d_MetP_Variable(ivar,MR_iMetStep_Now+1)
        vx_meso_next_step_MetP_sp = MR_vx_metP_next
        MR_dum3d_metP             = MR_vx_metP_next
          ! Now differentiate
        call MR_DelMetP_Dx
        du_dx_MetP_sp = MR_dum3d2_metP * M_2_KM
        call MR_DelMetP_Dy
        du_dy_MetP_sp = MR_dum3d2_metP * M_2_KM
        ! Load V winds on MetP
        vy_meso_next_step_MetP_sp = MR_vy_metP_next
        MR_dum3d_metP             = MR_vy_metP_next
          ! Now differentiate
        call MR_DelMetP_Dx
        dv_dx_MetP_sp = MR_dum3d2_metP * M_2_KM
        call MR_DelMetP_Dy
        dv_dy_MetP_sp = MR_dum3d2_metP * M_2_KM
        call Eddy_diff  ! this sets Khz in km2/hr
         ! Now resample onto computational grid
        MR_dum3d_metP = Khz_meso_next_step_MetP_sp
        call MR_Regrid_MetP_to_CompH(MR_iMetStep_Now+1)
        Khz_meso_next_step_sp = MR_dum3d_compH
      endif

      kx(1:nxmax,1:nymax,1:nzmax) = real( Khz_meso_last_step_sp(:,:,:),kind=ip) + &
                                    real((Khz_meso_next_step_sp(:,:,:) - &
                                          Khz_meso_last_step_sp(:,:,:)),kind=ip) * &
                                    real(Interval_Frac,kind=ip)
      ky = kx

      ! Set boundary kx and ky
        ! Bottom
      kx(0:nxmax+1,0:nymax+1,0) = kx(0:nxmax+1,0:nymax+1,1)
      ky(0:nxmax+1,0:nymax+1,0) = ky(0:nxmax+1,0:nymax+1,1)
        ! Top
      kx(0:nxmax+1,0:nymax+1,nzmax+1) = kx(0:nxmax+1,0:nymax+1,nzmax)
      ky(0:nxmax+1,0:nymax+1,nzmax+1) = ky(0:nxmax+1,0:nymax+1,nzmax)
        ! Left (West)
      kx(0,0:nymax+1,0:nzmax+1) = kx(1,0:nymax+1,0:nzmax+1)
      ky(0,0:nymax+1,0:nzmax+1) = ky(1,0:nymax+1,0:nzmax+1)
        ! Right (East)
      kx(nxmax+1,0:nymax+1,0:nzmax+1) = kx(nxmax,0:nymax+1,0:nzmax+1)
      ky(nxmax+1,0:nymax+1,0:nzmax+1) = ky(nxmax,0:nymax+1,0:nzmax+1)
        ! -y (South)
      kx(0:nxmax+1,0,0:nzmax+1) = kx(0:nxmax+1,1,0:nzmax+1)
      ky(0:nxmax+1,0,0:nzmax+1) = ky(0:nxmax+1,1,0:nzmax+1)
        ! +y (North)
      kx(0:nxmax+1,nymax+1,0:nzmax+1) = kx(0:nxmax+1,nymax,0:nzmax+1)
      ky(0:nxmax+1,nymax+1,0:nzmax+1) = ky(0:nxmax+1,nymax,0:nzmax+1)

      end subroutine Set_VarDiffH_Meso

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Set_VarDiffV_Meso
!
!  Called from: Ash3d.F90
!  Arguments: 
!    Load_MesoSteps,Interval_Frac
!
!  This subroutine calculates kz for the current time. Kv_meso_ is calculated
!  on the MetP grid from the virtual potential temperature, Ri, Monin-Length,
!  Ustar, and the planetary boundary layer height. Function calls are used to
!  populate those values on the NWP steps so that Kv_meso can be calculated,
!  then interpolated onto the compH grid at the current time.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Set_VarDiffV_Meso(Load_MesoSteps,Interval_Frac)

      use global_param,  only : &
         M2PS_2_KM2PHR,useDiffusion

      use mesh,          only : &
         nxmax,nymax,nzmax

      use Diffusion,     only : &
         kz

      use Atmosphere,    only : &
         Set_VirtPotenTemp, &
         Ri_meso_last_step_MetP_sp,Ri_meso_next_step_MetP_sp, &
         FricVel_meso_last_step_Met_sp,FricVel_meso_next_step_Met_sp,  &
         PBLH_meso_last_step_Met_sp,    PBLH_meso_next_step_Met_sp,      &
         L_MonOb_meso_last_step_Met_sp, L_MonOb_meso_next_step_Met_sp

      use MetReader,     only : &
         MR_iMetStep_Now,MR_dum3d_MetP,MR_dum3d_compH,&
           MR_Regrid_MetP_to_CompH

      implicit none

      logical      ,intent(in) :: Load_MesoSteps
      real(kind=dp),intent(in) :: Interval_Frac

      logical,save :: first_time = .true.

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Set_VarDiffV_Meso"
      endif;enddo

      ! To set the vertical diffusivity, we need to:
      !  1. Calculate the Richardson Number on MetP grid
      !  2. Calculate friction velocity (if not provided)
      !  3. Calculate boundary layer lengths
      !       Atmos. Boundary Layer Height (if not provided by Met file)
      !       surface layer thickness
      !       Monin-Obukhov Length
      !  4. Calculate Kv(Ri,u*,L,PBLz)
      ! We will calculate these values on the MetP grid, then interpolate Kv on
      ! the compH and then onto the current time.  
      ! Note: these are all non-linear functions so the better approach would be
      ! to evaluate everything on the computational grid at each time, but this
      ! is probably overkill.

      if(.not.useDiffusion)useDiffusion = .true.

      if(Load_MesoSteps)then
        if(first_time)then
          call Calc_SurfaceRoughnessLength
          !  Populate values for the 'last' step
          call Set_VirtPotenTemp(.true.)         ! True indicates load prestep
          call Calc_Ri(.true.)
          call Calc_Monin_Length(.true.)
          call Calc_SurfaceFrictionVelocity(.true.)
          call Calc_PBLH(.true.)

          call Calc_Vert_Diff(0)
          MR_dum3d_MetP = Kv_meso_last_step_MetP_sp
          call MR_Regrid_MetP_to_CompH(MR_iMetStep_Now)
          Kv_meso_last_step_sp = MR_dum3d_compH

          first_time = .false.
        else
          ! If we've already filled 'next', copy 'next' to 'last'
          Ri_meso_last_step_MetP_sp     = Ri_meso_next_step_MetP_sp
          L_MonOb_meso_last_step_Met_sp = L_MonOb_meso_next_step_Met_sp
          FricVel_meso_last_step_Met_sp = FricVel_meso_next_step_Met_sp
          PBLH_meso_last_step_Met_sp    = PBLH_meso_next_step_Met_sp
          Kv_meso_last_step_sp          = Kv_meso_next_step_sp
        endif ! first_time

          ! Populate Ri for the 'next' step
        call Set_VirtPotenTemp
        call Calc_Ri                              ! sets Ri_meso_next_step_MetP_sp
        call Calc_Monin_Length                    !  and L_MonOb_meso_next_step_Met_sp
        call Calc_SurfaceFrictionVelocity         ! sets FricVel_meso_next_step_Met_sp
        call Calc_PBLH

        call Calc_Vert_Diff(1)
        MR_dum3d_MetP = Kv_meso_next_step_MetP_sp
        call MR_Regrid_MetP_to_CompH(MR_iMetStep_Now+1)
        Kv_meso_next_step_sp = MR_dum3d_compH

      endif

      kz(1:nxmax,1:nymax,1:nzmax) = real( Kv_meso_last_step_sp(:,:,:),kind=ip) + &
                                    real((Kv_meso_next_step_sp(:,:,:) - &
                                          Kv_meso_last_step_sp(:,:,:)),kind=ip) * &
                                    real(Interval_Frac,kind=ip) * M2PS_2_KM2PHR
      ! Set boundary kz
        ! Bottom
      kz(0:nxmax+1,0:nymax+1,0) = kz(0:nxmax+1,0:nymax+1,1)
        ! Top
      kz(0:nxmax+1,0:nymax+1,nzmax+1) = kz(0:nxmax+1,0:nymax+1,nzmax)
        ! Left (West)
      kz(0,0:nymax+1,0:nzmax+1) = kz(1,0:nymax+1,0:nzmax+1)
        ! Right (East)
      kz(nxmax+1,0:nymax+1,0:nzmax+1) = kz(nxmax,0:nymax+1,0:nzmax+1)
        ! -y (South)
      kz(0:nxmax+1,0,0:nzmax+1) = kz(0:nxmax+1,1,0:nzmax+1)
        ! +y (North)
      kz(0:nxmax+1,nymax+1,0:nzmax+1) = kz(0:nxmax+1,nymax,0:nzmax+1)

      end subroutine Set_VarDiffV_Meso

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Calc_Ri(Load_Prestep)
!
!  Called from: Set_VarDiffV_Meso
!  Arguments:
!    Load_Prestep   = logical, optional; triggers loading 'last' only
!
!  This subroutine calculates the Richardson number, given the virtual potential
!  temperature. The 'next' data array is populated unless 'Load_Prestep' is given
!  and is .true.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Calc_Ri(Load_Prestep)

      use global_param,  only : &
         GRAV,KM_2_M

      use Atmosphere,    only : &
         Ri_meso_last_step_MetP_sp,Ri_meso_next_step_MetP_sp,&
         AirVPTemp_meso_last_step_MetP_sp,AirVPTemp_meso_next_step_MetP_sp

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet,p_fullmet_sp,&
         MR_geoH_metP_last,MR_geoH_MetP_next

      implicit none

      logical, intent(in), optional :: Load_Prestep

      real(kind=ip),parameter :: MIN_DVDZ = 3.0e-2_ip ! The minimum vertical shear assumed
                                                      ! No min leads to singular Ri
                                                      ! This value is assumed based on comparisons to MERRA

      real(kind=ip),dimension(:),allocatable :: z ! in m
      real(kind=ip),dimension(:),allocatable :: u ! in m/s
      real(kind=ip),dimension(:),allocatable :: v ! in m/s
      real(kind=ip),dimension(:),allocatable :: p ! in Pa
      real(kind=ip),dimension(:),allocatable :: Tpoten

      integer       :: i,j,k,k1,k2
      real(kind=ip) :: refP
      real(kind=ip) :: del_z
      real(kind=ip) :: dudz,dvdz,dtdz
      real(kind=ip) :: dveldz2
      real(kind=ip) :: temp_term,mech_term
      real(kind=ip) :: Ri

      logical      :: first_time

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Calc_Ri"
      endif;enddo

      if(present(Load_Prestep))then
        first_time = Load_Prestep
      else
        first_time = .false.
      endif

      allocate(z(np_fullmet))
      allocate(u(np_fullmet))
      allocate(v(np_fullmet))
      allocate(p(np_fullmet))
      allocate(Tpoten(np_fullmet))

      refP = 1.0e5_ip   ! reference pressure for potential temperature

      p(1:np_fullmet) = p_fullmet_sp(1:np_fullmet)
      if(first_time)then
        do i=1,nx_submet
          do j=1,ny_submet
            z(1:np_fullmet) = MR_geoH_metP_last(i,j,1:np_fullmet) * KM_2_M
            u(1:np_fullmet) = vx_meso_last_step_MetP_sp(i,j,1:np_fullmet)
            v(1:np_fullmet) = vy_meso_last_step_MetP_sp(i,j,1:np_fullmet)
            Tpoten(1:np_fullmet) = AirVPTemp_meso_last_step_MetP_sp(i,j,1:np_fullmet)
            do k=1,np_fullmet
              ! We need vertical derivatives of theta_v and u
              !  Use one-sided differences
              if(k.lt.np_fullmet)then
                k1 = k
                k2 = k+1
              elseif(k.eq.np_fullmet)then
                k1 = k-1
                k2 = k
              endif
              del_z  = z(k2)- z(k1)
              dtdz   = (Tpoten(k2)-Tpoten(k1)) / del_z
              dudz   = (u(k2)-u(k1)) / del_z
              dvdz   = (v(k2)-v(k1)) / del_z

              ! Only need magnitudes and nothing too close to 0 (leads to singular Ri)
              dudz   = max(abs(dudz),MIN_DVDZ)
              dvdz   = max(abs(dvdz),MIN_DVDZ)

              temp_term = dtdz/Tpoten(k)
              dveldz2   = dudz*dudz + dvdz*dvdz
              mech_term = dveldz2/GRAV
              Ri = real(temp_term / mech_term,kind=sp)

              ! Log this term since we will need it later when calculating free-air Kz
              dV_dz_MetP_sp(i,j,k) = real(sqrt(dveldz2),kind=sp)

              Ri_meso_last_step_MetP_sp(i,j,k) = real(Ri,kind=sp)
            enddo ! k
          enddo ! j
        enddo ! i
      else
        do i=1,nx_submet
          do j=1,ny_submet
            z(1:np_fullmet) = MR_geoH_MetP_next(i,j,1:np_fullmet) * KM_2_M
            u(1:np_fullmet) = vx_meso_next_step_MetP_sp(i,j,1:np_fullmet)
            v(1:np_fullmet) = vy_meso_next_step_MetP_sp(i,j,1:np_fullmet)
            Tpoten(1:np_fullmet) = AirVPTemp_meso_next_step_MetP_sp(i,j,1:np_fullmet)
            do k=1,np_fullmet
              ! We need vertical derivatives of theta_v and u
              !  Use one-sided differences
              if(k.lt.np_fullmet)then
                k1 = k
                k2 = k+1
              elseif(k.eq.np_fullmet)then
                k1 = k-1
                k2 = k
              endif
              del_z  = z(k2)- z(k1)
              dtdz   = (Tpoten(k2)-Tpoten(k1)) / del_z
              dudz   = (u(k2)-u(k1)) / del_z
              dvdz   = (v(k2)-v(k1)) / del_z

              ! Only need magnitudes and nothing too close to 0 (leads to singular Ri)
              dudz   = max(abs(dudz),MIN_DVDZ)
              dvdz   = max(abs(dvdz),MIN_DVDZ)

              temp_term = dtdz/Tpoten(k)
              dveldz2   = dudz*dudz + dvdz*dvdz
              mech_term = dveldz2/GRAV
              Ri = real(temp_term / mech_term,kind=sp)

              ! Log this term since we will need it later when calculating free-air Kz
              dV_dz_MetP_sp(i,j,k) = real(sqrt(dveldz2),kind=sp)

              Ri_meso_next_step_MetP_sp(i,j,k) = real(Ri,kind=sp)
            enddo ! k
          enddo ! j
        enddo ! i
      endif

      end subroutine Calc_Ri

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Calc_SurfaceRoughnessLength
!
!  Called from: Set_VarDiffV_Meso
!  Arguments: None
!    none
!
!  This subroutine loads the surface roughness variable from the first NWP file,
!  if available, and initializes to 0.1 if the variable is not provided by
!  the NWP file.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Calc_SurfaceRoughnessLength

      use MetReader,     only : &
         Met_var_IsAvailable,MR_iMetStep_Now,MR_dum2d_Met,&
         nx_submet,ny_submet, &
           MR_Read_2d_Met_Variable

      integer :: ivar

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Calc_SurfaceRoughnessLength"
      endif;enddo

      ! Check if the windfile being used provides surface roughness
      ivar = 17 ! Surface_roughness_surface
      if(Met_var_IsAvailable(ivar))then
        ! Surface roughness is provided, read it from the met file
        call MR_Read_2d_Met_Variable(ivar,MR_iMetStep_Now)
        SurfRoughLen_Met_sp(1:nx_submet,1:ny_submet)  = MR_dum2d_Met(1:nx_submet,1:ny_submet)
      else
        ! Set SurfRoughLen_Met_sp by assumption
          SurfRoughLen_Met_sp(1:nx_submet,1:ny_submet)  = 0.1_sp
      endif

      end subroutine Calc_SurfaceRoughnessLength

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Calc_SurfaceFrictionVelocity(Load_Prestep)
!
!  Called from: Set_VarDiffV_Meso
!  Arguments:
!    Load_Prestep   = logical, optional; triggers loading 'last' only
!
!  This subroutine sets the surface friction velocity on the Met grid for the
!  'last' (if Load_Prestep is provided) and 'next' data arrays. If Ustar is
!  provided by the NWP file, it is read and checked for validity. If it is not
!  provided by the NWP file (or if is contains invalid values), then Ustar is
!  calculated using the lowesti-level velocities available (10m or highest p).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Calc_SurfaceFrictionVelocity(Load_Prestep)

      use global_param,  only : &
         EPS_SMALL,KM_2_M

      use Atmosphere,    only : &
         FricVel_meso_last_step_Met_sp,FricVel_meso_next_step_Met_sp, &
         L_MonOb_meso_last_step_Met_sp,L_MonOb_meso_next_step_Met_sp

      use MetReader,     only : &
         Met_var_IsAvailable,MR_iMetStep_Now,np_fullmet,MR_iMetStep_Now,MR_Topo_met,&
         MR_geoH_metP_last,MR_geoH_metP_next,MR_dum2d_Met,nx_submet,ny_submet,&
           MR_Read_2d_Met_Variable

      implicit none

      logical, intent(in), optional :: Load_Prestep

      real(kind=ip) :: U_mag
      real(kind=ip) :: denom1,denom2
      real(kind=ip) :: z0
      real(kind=ip) :: L_MonOb
      real(kind=ip) :: zonL
      real(kind=sp),dimension(:,:),allocatable :: SurfVelx_meso_Met_sp
      real(kind=sp),dimension(:,:),allocatable :: SurfVely_meso_Met_sp
      real(kind=sp),dimension(:,:),allocatable :: SurfVelh_meso_Met_sp
      integer :: i,j,k
      integer :: ivar
      logical :: FV_override = .false.

      logical      :: first_time

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Calc_SurfaceFrictionVelocity"
      endif;enddo

      if(present(Load_Prestep))then
        first_time = Load_Prestep
      else
        first_time = .false.
      endif

      ! Check if the windfile being used provides friction velocity
      ivar = 13 ! Friction_velocity_surface
      if(Met_var_IsAvailable(ivar))then
        if(first_time)then
          ! Friction velocity is provided, read it from the met file
          call MR_Read_2d_Met_Variable(ivar,MR_iMetStep_Now)
            FricVel_meso_last_step_Met_sp = MR_dum2d_Met
          if(maxval(MR_dum2d_Met(:,:)).lt.EPS_SMALL)then
            ! variable was present in file, but filled with nonsense
            FV_override = .true.
          endif
        else
          call MR_Read_2d_Met_Variable(ivar,MR_iMetStep_Now+1)
            FricVel_meso_next_step_Met_sp = MR_dum2d_Met
          if(maxval(MR_dum2d_Met(:,:)).lt.EPS_SMALL)then
            ! variable was present in file, but filled with nonsense
            FV_override = .true.
          endif
        endif
      endif
      if(first_time)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          if(.not.Met_var_IsAvailable(ivar).or.FV_override)then
            write(outlog(io),*)"     Friction Velocity will be calculated."
          else
            write(outlog(io),*)"     Friction Velocity will be loaded from NWP file."
          endif
        endif;enddo
      endif

      if(.not.Met_var_IsAvailable(ivar).or.FV_override)then
        ! friction velocity is not provided by the met file (or is not valid)
        ! Calculate it ourselves
        ! Get surface friction velocity using Panofsky/Dutton p376

        allocate(SurfVelx_meso_Met_sp(nx_submet,ny_submet))
        allocate(SurfVely_meso_Met_sp(nx_submet,ny_submet))
        allocate(SurfVelh_meso_Met_sp(nx_submet,ny_submet))

        ! First check if we can get the 10m wind (both x and y)
        if(Met_var_IsAvailable(11).and.Met_var_IsAvailable(12))then
          if(first_time)then
            call MR_Read_2d_Met_Variable(11,MR_iMetStep_Now)
              SurfVelx_meso_Met_sp = MR_dum2d_Met
            call MR_Read_2d_Met_Variable(12,MR_iMetStep_Now)
              SurfVely_meso_Met_sp = MR_dum2d_Met
          else
            call MR_Read_2d_Met_Variable(11,MR_iMetStep_Now+1)
              SurfVelx_meso_Met_sp = MR_dum2d_Met
            call MR_Read_2d_Met_Variable(12,MR_iMetStep_Now+1)
              SurfVely_meso_Met_sp = MR_dum2d_Met
          endif
          SurfVelh_meso_Met_sp = 10.0_sp
          if(first_time)then
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"       Found surface winds in NWP file. Using these for"
              write(outlog(io),*)"       U* calculation."
            endif;enddo
          endif
        else
          ! If the 10m winds are not available, then use the lower levels of the 3d winds
          if(first_time)then
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"       No surface winds found in NWP file. Using lowest-level"
              write(outlog(io),*)"       winds for U* calculation."
            endif;enddo
          endif
          do i=1,nx_submet
            do j=1,ny_submet
              z0 = SurfRoughLen_Met_sp(i,j)
              if(first_time)then
                do k=1,np_fullmet
                  if(MR_geoH_metP_last(i,j,k)*KM_2_M.gt.1000.0_ip*MR_Topo_met(i,j)+z0)then
                    exit
                  endif
                enddo
                SurfVelx_meso_Met_sp(i,j) = vx_meso_last_step_MetP_sp(i,j,k)
                SurfVely_meso_Met_sp(i,j) = vy_meso_last_step_MetP_sp(i,j,k)
                SurfVelh_meso_Met_sp(i,j) = real((MR_geoH_metP_last(i,j,k)-MR_Topo_met(i,j))*KM_2_M,kind=sp)
              else
                do k=1,np_fullmet
                  if(MR_geoH_metP_next(i,j,k)*KM_2_M.gt.1000.0_ip*MR_Topo_met(i,j)+z0)then
                    exit
                  endif
                enddo
                SurfVelx_meso_Met_sp(i,j) = vx_meso_next_step_MetP_sp(i,j,k)
                SurfVely_meso_Met_sp(i,j) = vy_meso_next_step_MetP_sp(i,j,k)
                SurfVelh_meso_Met_sp(i,j) = real((MR_geoH_metP_next(i,j,k)-MR_Topo_met(i,j))*KM_2_M,kind=sp)
              endif
            enddo
          enddo
        endif

        ! Now we have Vx, Vy, and H needed for calculating Ust, either from 10m data or lower-level winds
        do i=1,nx_submet
          do j=1,ny_submet
            z0 = SurfRoughLen_Met_sp(i,j)
            U_mag = sqrt(SurfVelx_meso_Met_sp(i,j)**2.0_sp + &
                         SurfVely_meso_Met_sp(i,j)**2.0_sp)! / MPS_2_KMPHR
            denom1 = log(SurfVelh_meso_Met_sp(i,j)/z0)
            if(first_time)then
              L_MonOb   = real(L_MonOb_meso_last_step_Met_sp(i,j),kind=ip)
              zonL = SurfVelh_meso_Met_sp(i,j)/L_MonOb
              denom2 = Phi_WindShear_Similarity(zonL)
              FricVel_meso_last_step_Met_sp(i,j) = real(U_mag*vonKarman/(denom1+denom2),kind=sp)
            else
              L_MonOb   = real(L_MonOb_meso_next_step_Met_sp(i,j),kind=ip)
              zonL = SurfVelh_meso_Met_sp(i,j)/L_MonOb
              denom2 = Phi_WindShear_Similarity(zonL)
              FricVel_meso_next_step_Met_sp(i,j) = real(U_mag*vonKarman/(denom1+denom2),kind=sp)
            endif
          enddo
        enddo

        deallocate(SurfVelx_meso_Met_sp)
        deallocate(SurfVely_meso_Met_sp)
        deallocate(SurfVelh_meso_Met_sp)

      endif

      end subroutine Calc_SurfaceFrictionVelocity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Calc_PBLH(Load_Prestep)
!
!  Called from: Set_VarDiffV_Meso
!  Arguments:
!    Load_Prestep   = logical, optional; triggers loading 'last' only
!
!  This subroutine sets the planetary boundary layer height for the Met grid
!  at the 'last' (if Load_Prestep is provided) and 'next' data arrays. If
!  PBLH is not provided by the NWP file, then it is calculated by examining
!  the Eckman thickness (requires lat and Ustar), the temperature profile
!  for the lowest inversion, and the Richardson number.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Calc_PBLH(Load_Prestep)

      use global_param,  only : &
         EPS_SMALL,KM_2_M,DEG2RAD

      use Atmosphere,    only : &
         AirVPTemp_meso_last_step_MetP_sp,AirVPTemp_meso_next_step_MetP_sp, &
         Ri_meso_last_step_MetP_sp,Ri_meso_next_step_MetP_sp, &
         PBLH_meso_last_step_Met_sp,PBLH_meso_next_step_Met_sp ,&
         L_MonOb_meso_last_step_Met_sp,L_MonOb_meso_next_step_Met_sp, &
         FricVel_meso_last_step_Met_sp,FricVel_meso_next_step_Met_sp

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet,Met_var_IsAvailable,&
         MR_geoH_metP_last,MR_geoH_metP_next,MR_dum2d_Met,MR_iMetStep_Now,&
         MR_Have_LL_mapping,MR_xy2ll_ylat,IsLatLon_MetGrid,y_submet_sp,&
           MR_Read_2d_Met_Variable,MR_Set_LL_mapping

      implicit none

      logical, intent(in), optional :: Load_Prestep

      integer :: ivar
      integer :: i,j,k,kk
      real(kind=ip) :: denom,tmp
      real(kind=ip) :: Ri_col(np_fullmet)
      real(kind=ip) :: z_col(np_fullmet)
      real(kind=ip) :: vpt_col(np_fullmet)
      real(kind=ip) :: lapse
      real(kind=ip) :: lat
      real(kind=ip) :: L_MonOb
      real(kind=ip) :: Ust
      real(kind=ip) :: EckF
      real(kind=ip) :: PBLz
      real(kind=ip) :: cn,cs
      real(kind=ip),dimension(:,:,:),allocatable :: PBLtmp
      logical :: PBL_override = .false.

      logical      :: first_time

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Calc_PBLH"
      endif;enddo

      if(present(Load_Prestep))then
        first_time = Load_Prestep
      else
        first_time = .false.
      endif

      ! Check if the windfile being used provides PBLH
      ivar = 10 ! Planetary Boundary Level Height
      if(Met_var_IsAvailable(ivar))then
        ! PBLH is provided, read it from the met file
        if(first_time)then
          call MR_Read_2d_Met_Variable(ivar,MR_iMetStep_Now)
            PBLH_meso_last_step_Met_sp = MR_dum2d_Met
          if(maxval(MR_dum2d_Met(:,:)).lt.0.0_sp)then
            ! variable was present in file, but filled with nonsense
            PBL_override = .true.
          endif
        else
          call MR_Read_2d_Met_Variable(ivar,MR_iMetStep_Now+1)
            PBLH_meso_next_step_Met_sp = MR_dum2d_Met
          if(maxval(MR_dum2d_Met(:,:)).lt.0.0_sp)then
            ! variable was present in file, but filled with nonsense
            PBL_override = .true.
          endif
        endif
      endif
      if(first_time)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          if(.not.Met_var_IsAvailable(ivar).or.PBL_override)then
            write(outlog(io),*)"     Planetary Boundary Layer will be calculated."
          else
            write(outlog(io),*)"     Planetary Boundary Layer will be loaded from NWP file."
          endif
        endif;enddo
      endif

      if(.not.Met_var_IsAvailable(ivar).or.PBL_override)then
        ! If PBLH is not provided by the NWP file or is corrupted, then
        ! we need to calculate PBLH internally.  There are many more involved
        ! methods of determining the PBL, but flagging the first temperature
        ! inversion, using Ustar to determine the Eckman layer and using Ri with
        ! Ri_crit seem to work fairly well.

        ! From Sugiyama and Nasstrom, 1999
        ! https://digital.library.unt.edu/ark:/67531/metadc740014/m2/1/high_res_d/8191.pdf

        ! Calculate three measures of the PBL:
        !   (1) First, look for critical temperature inversion
        !   (2) Second approach uses the Eckman layer
        !   (3) Third approach used Ri with a critical threshold (~0.5)
        allocate(PBLtmp(nx_submet,ny_submet,3))
        PBLtmp(1:nx_submet,1:ny_submet,1:3) = 0.0_ip

        ! Eckman layer approach requires knowing the latitude; get that if neded
        if(.not.IsLatLon_MetGrid.and..not.MR_Have_LL_mapping)then
          call MR_Set_LL_mapping
        endif

        do i=1,nx_submet
          do j=1,ny_submet
            if(first_time)then
              Ri_col(:) = Ri_meso_last_step_MetP_sp(i,j,:)
              z_col(:)  = MR_geoH_metP_last(i,j,:)*KM_2_M
              vpt_col(:)= AirVPTemp_meso_last_step_MetP_sp(i,j,:)
              L_MonOb   = real(L_MonOb_meso_last_step_Met_sp(i,j),kind=ip)
              Ust       = real(FricVel_meso_last_step_Met_sp(i,j),kind=ip)
            else
              Ri_col(:) = Ri_meso_next_step_MetP_sp(i,j,:)
              z_col(:)  = MR_geoH_metP_next(i,j,:)*KM_2_M
              vpt_col(:)= AirVPTemp_meso_next_step_MetP_sp(i,j,:)
              L_MonOb   = real(L_MonOb_meso_next_step_Met_sp(i,j),kind=ip)
              Ust       = real(FricVel_meso_next_step_Met_sp(i,j),kind=ip)
            endif

            ! (1) temperature inversion
            kk = 0
            do k = 2,np_fullmet
              lapse = -(vpt_col(k) - vpt_col(k-1)) / &
                       (z_col(k)   - z_col(k-1))
              if(lapse.le.-0.005_ip)then
                kk = k
                exit
              endif
            enddo
            if(kk.gt.0)then
              PBLz = z_col(kk) - 2.0_ip/lapse
            else
              PBLz = EPS_SMALL
            endif
            PBLtmp(i,j,1) = PBLz

            ! (2) Eckman thickness
            !     This always has a non-zero solution
            if(IsLatLon_MetGrid)then
              lat = real(y_submet_sp(j),kind=ip)
            else
              lat = real(MR_xy2ll_ylat(i,j),kind=ip)
            endif
            lat = max(20.0_ip,abs(lat));
            EckF= 2.0_ip*7.292e-5_ip*sin(lat*DEG2RAD);
            tmp = abs(Ust/EckF/L_MonOb);
            if(tmp.lt.4.0_ip)then
              cn = 0.2;
              PBLtmp(i,j,2) = cn*Ust/abs(EckF);
            else
              cs = 0.4;
              PBLtmp(i,j,2) = cs*sqrt(abs(Ust*L_MonOb/EckF));
            endif

            ! (3) evaluate Ri crit
              ! Initialize boundary layer height to sea level
            PBLz = EPS_SMALL
            do k = 2,np_fullmet-1
              if(Ri_col(k).gt.RI_CRIT.and.Ri_col(k-1).le.RI_CRIT &
                 .and.z_col(k).lt.3000.0_ip)then ! We need this upper limit of 3km to avoid
                                                 ! missing the PBL and flagging the tropopause

                ! This height is above the PBL; interpolate back to
                ! k-1 to get PBLz
                if(abs(Ri_col(k)-Ri_col(k-1)).lt.EPS_SMALL)cycle
                denom = (Ri_col(k)-Ri_col(k-1))
                if(abs(denom).lt.1.0e-3_ip)then
                  tmp = 1.0_ip
                else
                  tmp = (RI_CRIT-Ri_col(k-1)) / denom
                  tmp = min(tmp,1.0_ip)
                endif
                PBLz = z_col(k-1)+tmp*(z_col(k)-z_col(k-1))
                ! if we have a PBLz, then exit the do loop
                exit
              endif
            enddo
    
            ! Make sure that PBLz is not negative
            PBLtmp(i,j,3) = max(PBLz,EPS_SMALL)

            ! Set PBLH to temperature inversion height
            !PBLz = PBLtmp(i,j,1)
            ! Set PBLH to Eckman thickness
            !PBLz = PBLtmp(i,j,2)
            ! Set PBLH to height of Ri_crit
            !PBLz = PBLtmp(i,j,3)

            PBLz = maxval(PBLtmp(i,j,1:3))

            if(first_time)then
              PBLH_meso_last_step_Met_sp(i,j) = real(PBLz,kind=sp)
            else
              PBLH_meso_next_step_Met_sp(i,j) = real(PBLz,kind=sp)
            endif
          enddo
        enddo
        deallocate(PBLtmp)
      endif

      return

      end subroutine Calc_PBLH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Calc_Monin_Length(Load_Prestep)
!
!  Called from: Set_VarDiffV_Meso
!  Arguments:
!    Load_Prestep   = logical, optional; triggers loading 'last' only
!
!  This subroutine Monin-Obukhov length to be used in surface-layer calculations.
!  The equation for L can be used throughout, but we are only interested in the
!  lowest layer (surface) so this subroutine returns a 2d array, either 'last'
!  (if Load_Prestep is provided) or the 'next' data array.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Calc_Monin_Length(Load_Prestep)

      use global_param,  only : &
         EPS_SMALL,KM_2_M

      use Atmosphere,    only : &
         Ri_meso_last_step_MetP_sp,Ri_meso_next_step_MetP_sp,&
         L_MonOb_meso_last_step_Met_sp,L_MonOb_meso_next_step_Met_sp

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet,MR_geoH_metP_last,MR_geoH_metP_next,&
           MR_Read_2d_Met_Variable

      implicit none

      logical, intent(in), optional :: Load_Prestep

      integer :: i,j,k_L
      real(kind=ip) :: Ri
      real(kind=ip) :: Ri_col(np_fullmet)
      real(kind=ip) :: z_col(np_fullmet)
      real(kind=ip) :: L_MonOb

      logical      :: first_time

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Calc_Monin_Length"
      endif;enddo

      if(present(Load_Prestep))then
        first_time = Load_Prestep
      else
        first_time = .false.
      endif

      ! Get Monin-Obukhov length from the
      ! Businger-Dyer-Pandolfo empirical result 
      ! using z and Ri at k=2
        ! Eq 6.7.1 and 6.7 2 of "Atmospheric Turbulence";
        ! Panofsky and Dutton,1984
        ! also Eq 11.24 of "Introduction to Micrometeorology";
        ! Arya, 1988
      if(first_time)then
        do i=1,nx_submet
          do j=1,ny_submet
            Ri_col(:) = Ri_meso_last_step_MetP_sp(i,j,:)
            z_col(:)  = MR_geoH_metP_last(i,j,:)*KM_2_M
            ! Pick the bottom (non-zero) z
            do k_L=1,np_fullmet
              if(z_col(k_L).gt.0.0_ip)then
                exit
              endif
            enddo
            ! For the purpuse of calculating L, don't let Ri get too close to 0
            Ri = sign(max(abs(Ri_col(k_L)),1.0e-2_ip),Ri_col(k_L))
            if(Ri_col(k_L).lt.RI_CRIT)then
              ! Test for special case of neutrally stable case, L->Inf ; set to 1 km
              if (abs(Ri_col(k_L)).lt.EPS_SMALL)then
                L_MonOb = sign(abs(L_MonOb),100.0_ip)
              else
                ! Unstable (negative L)
                L_MonOb = z_col(k_L)/Ri
              endif
            elseif(Ri_col(k_L).gt.RI_CRIT)then
                ! Stable (positive L)
              L_MonOb = z_col(k_L)/Ri * (1.0_ip - 5.0_ip*Ri)
            endif
            L_MonOb = sign(min(abs(L_MonOb),100.0_ip),L_MonOb)

            L_MonOb_meso_last_step_Met_sp(i,j) = real(L_MonOb,kind=sp)
          enddo
        enddo
      else
        do i=1,nx_submet
          do j=1,ny_submet
            Ri_col(:) = Ri_meso_next_step_MetP_sp(i,j,:)
            z_col(:)  = MR_geoH_metP_next(i,j,:)*KM_2_M
            ! Pick the bottom (non-zero) z
            do k_L=1,np_fullmet
              if(z_col(k_L).gt.0.0_ip)then
                exit
              endif
            enddo
            ! For the purpuse of calculating L, don't let Ri get too close to 0
            Ri = sign(max(abs(Ri_col(k_L)),1.0e-2_ip),Ri_col(k_L))
            if(Ri_col(k_L).lt.RI_CRIT)then
              ! Test for special case of neutrally stable case, L->Inf ; set to 1 km
              if (abs(Ri_col(k_L)).lt.EPS_SMALL)then
                L_MonOb = sign(abs(L_MonOb),100.0_ip)
              else
                ! Unstable (negative L)
                L_MonOb = z_col(k_L)/Ri
              endif
            elseif(Ri_col(k_L).gt.RI_CRIT)then
                ! Stable (positive L)
              L_MonOb = z_col(k_L)/Ri * (1.0_ip - 5.0_ip*Ri)
            endif
            L_MonOb = sign(min(abs(L_MonOb),100.0_ip),L_MonOb)

            L_MonOb_meso_next_step_Met_sp(i,j) = real(L_MonOb,kind=sp)
          enddo
        enddo
      endif

      return

      end subroutine Calc_Monin_Length

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Fc_Louis
!
!  Called from: Calc_Vert_Diff
!  Arguments:
!    Ri    : Richardson #
!    zonz0 : z / z0
!
!  This function is the dimensionless stability function for free-air calculations
!  of Kv.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function Fc_Louis(Ri,zonz0)

      ! Stability function for vertical diffusion in the free atmosphere above the PBL
      ! Described in Louis 1979

      implicit none

      real(kind=ip) :: Fc_Louis ! dimensionless
      real(kind=ip) :: Ri       ! dimensionless
      real(kind=ip) :: zonz0    ! dimensionless

      real(kind=ip) :: a,b,c,bprime

      do io=1,2;if(VB(io).le.verbosity_debug2)then
        write(outlog(io),*)"     Entered function Fc_Louis"
      endif;enddo

      a = vonKarman/log(zonz0)    ! Eq. 13
      b = 9.4_ip
      bprime = 0.5_ip*b
      c = 7.4*a*a*b*sqrt(zonz0)   ! Eq. 20

      if(Ri.le.0.0_ip)then
          ! Unstable atmosphere
          ! Eq. 14
        Fc_Louis = 1.0_ip - b*Ri/(1.0_ip+c*sqrt(abs(Ri)))
      else
          ! Stable atmosphere
          ! Eq. 15
        Fc_Louis = (1.0_ip + bprime*Ri)**(-2.0_ip)
      endif

      return

      end function Fc_Louis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Fc_Jac
!
!  Called from: Calc_Vert_Diff
!  Arguments:
!    Ri    : Richardson #
!
!  This function is the dimensionless stability function for free-air calculations
!  of Kv.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function Fc_Jac(Ri)

      ! Stability function for vertical diffusion in the free atmosphere above the PBL
      ! Described in Jacobson (p251 Eq.8.70)
      ! An earlier reference is Stull (1988) Table 6.4 of Chap 6 Turbulence Closure

      implicit none

      real(kind=ip) :: Fc_Jac ! dimensionless
      real(kind=ip) :: Ri ! dimensionless

      do io=1,2;if(VB(io).le.verbosity_debug2)then
        write(outlog(io),*)"     Entered function Fc_Jac"
      endif;enddo

      if(Ri.le.0.0_ip)then
          ! Unstable atmosphere
          ! Tab 6.4 of Stull (unstable, line 2)
        Fc_Jac = sqrt(1.0_ip-18.0_ip*Ri)
      elseif(Ri.ge.0.0_ip.and.Ri.le.RI_CRIT)then
          ! Weakly unstable atmosphere
          ! Jacobson Eq. 8.70 or Table 6.4 of Stull (unstable, line 1)
        Fc_Jac = (RI_CRIT-Ri)/RI_CRIT
      else
          ! Stable atmosphere
        Fc_Jac = 0.0_ip
      endif

      return

      end function Fc_Jac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Fc_Betts
!
!  Called from: Calc_Vert_Diff
!  Arguments:
!    Ri    : Richardson #
!
!  This function is the dimensionless stability function for free-air calculations
!  of Kv.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function Fc_Betts(Ri)

      ! Stability function for vertical diffusion in the free atmosphere above the PBL
      ! Originally from Betts et al, 1996

      implicit none

      real(kind=ip) :: Fc_Betts ! dimensionless
      real(kind=ip) :: Ri ! dimensionless

      do io=1,2;if(VB(io).le.verbosity_debug2)then
        write(outlog(io),*)"     Entered function Fc_Betts"
      endif;enddo

      if(Ri.le.0.0_ip)then
          ! Unstable atmosphere
          ! Eq. A20b
        Fc_Betts = 1.0_ip-8.0_ip*Ri/(1.0_ip+1.746_ip*sqrt(-Ri))
      else
          ! Stable atmosphere
          ! Eq. A18
        Fc_Betts = (1.0_ip + 5.0_ip*Ri)**(-2.0_ip)
      endif

      return

      end function Fc_Betts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Fc_Hong
!
!  Called from: Calc_Vert_Diff
!  Arguments:
!    Ri    : Richardson #
!
!  This function is the dimensionless stability function for free-air calculations
!  of Kv.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function Fc_Hong(Ri)

      ! Stability function for vertical diffusion in the free atmosphere above the PBL
      ! Originally from Hong and Pan, 1996

      implicit none

      real(kind=ip) :: Fc_Hong ! dimensionless
      real(kind=ip) :: Ri ! dimensionless

      do io=1,2;if(VB(io).le.verbosity_debug2)then
        write(outlog(io),*)"     Entered function Fc_Hong"
      endif;enddo

      if(Ri.le.0.0_ip)then
          ! Unstable atmosphere
          ! Eq. 5
        Fc_Hong = (1.0_ip-1.6_ip*Ri)**(-0.25_ip)
      else
          ! Stable atmosphere
          ! Eq. 13
        Fc_Hong = exp(-8.5_ip*Ri) + 0.15_ip/(Ri+3.0_ip)
      endif

      return

      end function Fc_Hong

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Fc_Collins
!
!  Called from: Calc_Vert_Diff
!  Arguments:
!    Ri    : Richardson #
!
!  This function is the dimensionless stability function for free-air calculations
!  of Kv.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function Fc_Collins(Ri)

      ! Stability function for vertical diffusion in the free atmosphere above the PBL
      ! Originally from Collins et al, NCAR TN-464, 2004
      ! http://www.cesm.ucar.edu/models/atm-cam/docs/description/description.pdf

      implicit none

      real(kind=ip) :: Fc_Collins ! dimensionless
      real(kind=ip) :: Ri ! dimensionless

      do io=1,2;if(VB(io).le.verbosity_debug2)then
        write(outlog(io),*)"     Entered function Fc_Collins"
      endif;enddo

      if(Ri.le.0.0_ip)then
          ! Unstable atmosphere
          ! Eq. 4.464  (Also Tab 6.4 of Stull)
        Fc_Collins = sqrt(1.0_ip-18.0_ip*Ri)
      else
          ! Stable atmosphere
          ! Eq. 4.465 
        Fc_Collins = 1.0_ip/(1.0_ip+10.0_ip*Ri*(1.0_ip+8.0_ip*Ri))
      endif

      return

      end function Fc_Collins

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Fc_PMB
!
!  Called from: No called (since it isn't working), but would be called from Calc_Vert_Diff
!  Arguments:
!    Ri    : Richardson #
!    ml    : mixing length
!    z     : altitude
!
!  This function is the dimensionless stability function for free-air calculations
!  of Kv.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function Fc_PMB(Ri,ml,z)

      ! Stability function for vertical diffusion
      ! see Piedelievre, Jean Philippe, Lue Musson-Genon, François Bompay, 1990:
      ! MEDIA—An Eulerian Model of Atmospheric Dispersion: First Validation on
      ! the Chernobyl Release. J. Appl. Meteor., 29, 1205–1220.
      ! doi: http://dx.doi.org/10.1175/1520-0450(1990)029<1205:MEMOAD>2.0.CO;2 
      ! This is the model used by MDLP0 (also used in EMERRAUDE or PERIDOT).
      ! It also actually used a subsequent F form from Louis (1980).

      implicit none

      real(kind=ip) :: Fc_PMB  ! dimensionless
      real(kind=ip) :: Ri      ! dimensionless
      real(kind=ip) :: ml      ! m
      real(kind=ip) :: z       ! m
        ! These are kept written as separate constants for consistancy with the
        ! equation in the paper above
      real(kind=ip),parameter :: b = 5.0_ip
      real(kind=ip),parameter :: c = 5.0_ip
      real(kind=ip),parameter :: d = 5.0_ip
      real(kind=ip),parameter :: f = 5.19615242270663_ip ! = sqrt(27)

      do io=1,2;if(VB(io).le.verbosity_debug2)then
        write(outlog(io),*)"     Entered function Fc_PMB"
      endif;enddo

      if(Ri.le.0.0_ip)then
          ! Unstable atmosphere
        Fc_PMB = 1.0_ip/(1.0_ip+3.0_ip*b*c*(ml*ml*sqrt(Ri)/(z*z*f)))
      else
          ! Stable atmosphere
        Fc_PMB = 1.0_ip/(1.0_ip+3.0_ip*b*Ri*sqrt(1.0_ip+d*Ri))
      endif

      return

      end function Fc_PMB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  MixLen
!
!  Called from: Calc_Vert_Diff
!  Arguments:
!    z   : altitude
!
!  This function calculates the mixing legth for free-air Kv using the asymptotic
!  mixing length LambdaC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function MixLen(z)

      ! Returns the mixing length for Prandtl's turbulent diffusion
      ! See Collins et al, NCAR TN-464, 2004, eq. 4.461
      ! http://www.cesm.ucar.edu/models/atm-cam/docs/description/description.pdf
      ! This is originally from Blackadar (1962)
      ! Note that LambdaC is the free-atmospheric mixing length and is set in
      ! the optional module block.
      !  Jacobson (p251) gives the range of 70-200m
      !  Louis (2000) found LambdaC to be 100m
      !  Collins used 30m in the CAM3 model
      !  Piedelievere et al (1990) used 150 m in the MEDIA model

      implicit none

      real(kind=ip) :: MixLen  ! m
      real(kind=ip) :: z       ! m

      do io=1,2;if(VB(io).le.verbosity_debug2)then
        write(outlog(io),*)"     Entered function MixLen"
      endif;enddo

      !MixLen = 1.0_ip/(1.0_ip/(z*vonKarman) + 1.0_ip/LambdaC)
      MixLen = (z*vonKarman)/(1.0_ip+(z*vonKarman/LambdaC))

      return

      end function MixLen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Phi_WindShear_Similarity
!
!  Called from: Calc_Vert_Diff
!  Arguments:
!    z_on_L
!
!  This function calculates the stability function from Monin similarity theory
!  for surface layers. It assumes that the function is a general form with
!  parameters alpha, beta, gamma, set in the input block.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function Phi_WindShear_Similarity(z_on_L)

      ! Generalized wind shear similarity function of z/L
      ! Coefficients and exponents are module variables

      implicit none

      real(kind=ip) :: Phi_WindShear_Similarity

      real(kind=ip) :: z_on_L

      do io=1,2;if(VB(io).le.verbosity_debug2)then
        write(outlog(io),*)"     Entered function Phi_WindShear_Similarity"
      endif;enddo

      if(z_on_L.le.0.0_ip)then
          ! Unstable
        Phi_WindShear_Similarity = (1.0_ip + phi_gamma*z_on_L)**phi_alpha
      else
          ! Stable
        Phi_WindShear_Similarity = (1.0_ip + phi_beta*z_on_L)
      endif

      return

      end function Phi_WindShear_Similarity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Psi_WindShear_Similarity
!
!  Called from: Currently not called
!  Arguments:
!    z_on_L
!
!  This function calculates the Psi from Monin similarity theory.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function Psi_WindShear_Similarity(z_on_L)

      ! Generalized wind shear similarity function of z/L (integral of phi)

      use global_param,  only : &
         PI

      implicit none

      real(kind=ip) :: Psi_WindShear_Similarity

      real(kind=ip) :: z_on_L

      real(kind=ip) :: x,tmp1,tmp2

      do io=1,2;if(VB(io).le.verbosity_debug2)then
        write(outlog(io),*)"     Entered function Psi_WindShear_Similarity"
      endif;enddo

      if(z_on_L.le.0.0_ip)then
          ! Unstable
        x = (1.0_ip-16.0_ip*z_on_L)**0.25_ip
        tmp1 = 0.5_ip*(1.0_ip+x*x);
        tmp2 = (0.5_ip*(1.0_ip+x))**2.0_ip;
        Psi_WindShear_Similarity = log(tmp1*tmp2) - 2.0_ip*atan(x) + 0.5_ip*PI;
      else
          ! Stable
        Psi_WindShear_Similarity = -5.0_ip*z_on_L
      endif

      return

      end function Psi_WindShear_Similarity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module Diffusivity_Variable

!##############################################################################
