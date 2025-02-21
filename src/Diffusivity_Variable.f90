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
!  to scale 
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
!      function Fc
!      function Fc_PMB
!      function MixLen
!      function Phi_WindShear_Similarity
!      function Psi_WindShear_Similarity
!
!
! OPTMOD=VARDIFF
! Horizontal diff
!  yes 1 500.0     # 1=const ; value in m2/s
!  yes 2 0.2       # 2=Smagorinsky ; C
!  yes 3 0.2       # 3=Pielke      ; C
! Vertical diff
!  yes 
!  1 500.0         # BL model 1=const ; value
!  2               #          2=Troen and Mahrt
!  3               #          3=Ulke
!  4               #          4=Shir / Businger,Ayer
!  1 500.0         # Free-Air 1=const ; value
!  2               #          2=Jac
!  3               #          3=Collin
!  4               #          4=Piedelievre
!0.4                         # vonKarman
!30.0                        # LambdaC
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
      integer :: KvBL_model_ID
      integer :: KvFA_model_ID

      !  These are the parameters that control the diffusivity calculations
      !    C from Smagorinsky model of horizontal diffusivity
      real(kind=ip) :: KH_SmagC     ! Smagorinsky (1993) constant for LES horizontal diffusivity (0.2 - 0.9)
      !    These next three are needed for the vertical diffusivity
      real(kind=ip) :: vonKarman    ! von Karman constant (around 0.4)
      real(kind=ip) :: LambdaC      ! Asymptotic length scale (around 30 m)
      real(kind=ip) :: RI_CRIT      ! Critical Richardson number (0.25)

      !    These are the values controlling the stability function Phi (lots of models out there)
      !   source              alpha   beta   gamma
      ! Businger-Dyer (1971)  -1/4    5.0    -16.0
      ! Carl (1973)           -1/3    5.0    -15.0
      ! Troen-Mahrt (1986)    -1/3    4.7     -7.0
      ! Ulke (2000)           -1/2    9.2    -13.0
      real(kind=ip) :: phi_alpha = -0.33333_ip   ! Exponent in unstable term
      real(kind=ip) :: phi_beta  =  4.7_ip       ! Coefficient in stable term (pretty much always 4.7->5.2
      real(kind=ip) :: phi_gamma = -7.0_ip       ! Coefficient in unstable term

      real(kind=ip) :: PBL_exp = 1.0_ip

      real(kind=ip) :: diffusivity_BL
      ! Set the number of output variables for this module
      ! This depends on settings from the input block
      logical :: use_Output_Vars_VarDiff       = .true.
      integer :: nvar_User2d_static_XY_VarDiff = 0
      integer :: nvar_User2d_XY_VarDiff        = 0 ! If using Kz, then =1 : Pbl
      integer :: nvar_User3d_XYGs_VarDiff      = 0
      integer :: nvar_User3d_XYZ_VarDiff       = 0 ! If using Kh, then =1 khorz; if also Kz, then =3 kvert, Ri
      integer :: nvar_User4d_XYZGs_VarDiff     = 0

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

      ! These are used to keep track of which index in the global list, this
      ! modules output vars corespond to
      integer :: indx_User2d_static_XY_VarDiff
      integer :: indx_User2d_XY_VarDiff
      integer :: indx_User3d_XYGs_VarDiff
      integer :: indx_User3d_XYZ_VarDiff
      integer :: indx_User4d_XYZGs_VarDiff

      real(kind=ip) :: LES_L2ScaleCoeff

      ! Note: RoughLen_z can be related to Land use
      !    From Stohl et al, ACP, v5n9p2461, 2005 Table 3
      !         Grassland       :: 0.10
      !         Arable land     :: 0.15
      !         Permanent crops :: 0.30
      !         Forest          :: 0.60
      !         Inland water    :: Charnock
      !         Urban areas     :: 0.70
      !         Other           :: 0.10
      !         Ocean           :: Charnock
      !   Note: Charnok relation is z_0 = a(u_star^2/g) with a~ 0.018
      !    From Stohl et al, Tech Note FLEXPART 8.2 :: surfdata.t
      !         landuse   comment                               z0        glcf
      !         --------------------------------------------------------
      !          1 Urban land                                   0.7       13
      !          2 Agricultural land                            0.1       11
      !          3 Range land                                   0.1       10
      !          4 Deciduous forest                             1.         3,4
      !          5 Coniferous forest                            1.         1,2
      !          6 Mixed forest including wetland               0.7        5
      !          7 water, both salt and fresh                   0.001      0
      !          8 barren land mostly desert                    0.01      12
      !          9 nonforested wetland                          0.1
      !         10 mixed agricultural and range land            0.1        6,7,8,9
      !         11 rocky open areas with low grow shrubs        0.05 
      !         12 snow and ice                                 0.001
      !         13 rainforest                                   1.
      ! For resuspension cases, friction velocity is needed at every cell and so
      ! will the z0 (RoughLen_z).  Vr (U10 and V10) need to be calculated of
      ! read for the wind grid and regridded to the comp grid.
      
      ! 3d Variables needed on MetP grid
      real(kind=sp),dimension(:,:,:)  ,allocatable :: dVel_dz_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: du_dx_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: du_dy_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: dv_dx_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: dv_dy_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: dV_dz_MetP_sp
      real(kind=sp),dimension(:,:)    ,allocatable :: SurfRoughLen_Met_sp

      ! HFS: Consider moving Ri, PBLH, L_MonOb, FricVel, TropoH, SurfRoughLen, displacement height to Atmosphere
      !      These would still be allocated here if needed here

        ! and at both meso steps (also MetP)
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Ri_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Ri_meso_next_step_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Khz_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Khz_meso_next_step_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Kv_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Kv_meso_next_step_MetP_sp

      ! 2d variables needed at meso steps on Met grid
      real(kind=sp),dimension(:,:)    ,allocatable :: PBLH_meso_last_step_Met_sp
      real(kind=sp),dimension(:,:)    ,allocatable :: PBLH_meso_next_step_Met_sp
      real(kind=sp),dimension(:,:)    ,allocatable :: L_MonOb_meso_last_step_Met_sp
      real(kind=sp),dimension(:,:)    ,allocatable :: L_MonOb_meso_next_step_Met_sp
      real(kind=sp),dimension(:,:)    ,allocatable :: FricVel_meso_last_step_Met_sp
      real(kind=sp),dimension(:,:)    ,allocatable :: FricVel_meso_next_step_Met_sp
      real(kind=sp),dimension(:,:)    ,allocatable :: v10x_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:)    ,allocatable :: v10y_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:)    ,allocatable :: v10x_meso_next_step_MetP_sp
      real(kind=sp),dimension(:,:)    ,allocatable :: v10y_meso_next_step_MetP_sp

      ! Variables needed on Comp grid (kx,y,z are already allocated)
      real(kind=sp),dimension(:,:)    ,allocatable :: FricVel_meso_last_step_sp
      real(kind=sp),dimension(:,:)    ,allocatable :: FricVel_meso_next_step_sp
      real(kind=ip),dimension(:,:)    ,allocatable :: FricVel_ip
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Khz_meso_last_step_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Khz_meso_next_step_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Kv_meso_last_step_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: Kv_meso_next_step_sp

      ! Both Khz adn Kv need U and V values on MetP grid so store local copies
      ! Note: The core Ash3d code reads directly into the computational grid
      real(kind=sp),dimension(:,:,:)  ,allocatable :: vx_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: vy_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: vx_meso_next_step_MetP_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable :: vy_meso_next_step_MetP_sp

      contains

!******************************************************************************

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
      integer            :: ios,ioerr
      character(len=20)  :: mod_name
      integer            :: substr_pos
      real(kind=ip)      :: tmp

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

      useVarDiffH = .false.
      useVarDiffV = .false.
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"    Continue reading input file for VarDiff block"
      endif;enddo

      !Check if we're going to use variable diffusivity
      read(10,'(a80)',iostat=ios,err=2010)linebuffer080
      read(linebuffer080,'(a3)',err=2011) answer
      if (answer.eq.'yes') then
        useVarDiffH = .true.  ! might be changed back below if we are holding Kh constant:w

        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"    Using horizontal variable diffusivity"
        endif;enddo
        ! Try to read the horizontal model ID
        read(linebuffer080(4:),*,iostat=ios)Kh_model_ID,tmp
        if(ios.eq.0)then
          if(Kh_model_ID.eq.1)then
            useVarDiffH = .false.
            diffusivity_horz = tmp
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"    Horizontal diffusivity model ID = 1: Constant"
              write(outlog(io),*)"                            with Kh = ",diffusivity_horz
            endif;enddo
          elseif(Kh_model_ID.eq.2)then
            KH_SmagC = tmp
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"    Horizontal diffusivity model ID = 2: Smagorinsky (1963)"
              write(outlog(io),*)"                             with C = ",KH_SmagC
            endif;enddo
          elseif(Kh_model_ID.eq.3)then
            KH_SmagC = tmp
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"    Horizontal diffusivity model ID = 3: Pielke (1974)"
              write(outlog(io),*)"                             with C = ",KH_SmagC
            endif;enddo
          else
            KH_SmagC = 0.2_ip
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"    Horizontal diffusivity model ID not recognized."
              write(outlog(io),*)"    Using model ID = 2: Smagorinsky (1963)"
              write(outlog(io),*)"                             with C = ",KH_SmagC
            endif;enddo
          endif
        else
          Kh_model_ID = 2
        endif
      elseif(answer(1:2).eq.'no') then
        useVarDiffH = .false.
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"    Not using horizontal variable diffusivity"
        endif;enddo
      else
        goto 2011
      endif
      read(10,'(a80)',iostat=ios,err=2010)linebuffer080
      read(linebuffer080,'(a3)',err=2011) answer
      if (answer.eq.'yes') then
        useVarDiffV = .true.
        useTemperature = .true.
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"    Using vertical variable diffusivity"
        endif;enddo
      elseif(answer(1:2).eq.'no') then
        useVarDiffV = .false.
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"    Not using vertical variable diffusivity"
        endif;enddo
      else
        goto 2011
      endif
      if(useVarDiffV)then
        ! Need to read two more lines defining first the Boundary Layer model, then the Free-Air model
        read(10,'(a80)',iostat=ios,err=2010)linebuffer080
        read(linebuffer080,*,iostat=ios)KvBL_model_ID
        if(KvBL_model_ID.eq.1)then
          ! Diffusivity is constant in the BL; read the value
          read(linebuffer080,*,iostat=ios)KvBL_model_ID,diffusivity_BL
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using constant vertical diffusivity in the boundary layer."
          endif;enddo
        elseif(KvBL_model_ID.eq.2)then
          ! Model from Troen and Mahrt, 1973
          phi_alpha = -0.33333_ip   ! Exponent in unstable term
          phi_beta  =  4.7_ip       ! Coefficient in stable term (pretty much always 4.7->5.2
          phi_gamma = -7.0_ip       ! Coefficient in unstable term
          PBL_exp   = 2.0_ip
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using boundary layer vertical diffusivity as outlined by Troen and Mahrt (1973)."
          endif;enddo
        elseif(KvBL_model_ID.eq.3)then
          ! Model from Ulke (2000)
          phi_alpha = -0.5_ip   ! Exponent in unstable term
          phi_beta  =  9.2_ip       ! Coefficient in stable term (pretty much always 4.7->5.2
          phi_gamma = -13.0_ip       ! Coefficient in unstable term
          PBL_exp   = 1.0_ip
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using boundary layer vertical diffusivity as outlined by Ulke (2000)."
          endif;enddo
        elseif(KvBL_model_ID.eq.4)then
          ! Model from Shir / Businger,Ayer outlined in Seinfeld and Pandis
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using boundary layer vertical diffusivity as outlined by in Seinfeld and Pandis."
          endif;enddo
        else
          KvBL_model_ID = 2
         do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using boundary layer vertical diffusivity as outlined by Troen and Mahrt (1973)."
          endif;enddo
        endif

        read(linebuffer080(4:),*,iostat=ios)KvFA_model_ID
        if(KvFA_model_ID.eq.1)then
          ! Diffusivity is constant in the BL; read the value
          read(linebuffer080,*,iostat=ios)KvFA_model_ID,diffusivity_vert
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using constant vertical diffusivity above the boundary layer."
          endif;enddo
        elseif(KvFA_model_ID.eq.2)then
          ! Mixing length model with Fc from Jacobson
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using free air vertical diffusivity with stability function from Jacobson."
          endif;enddo
        elseif(KvFA_model_ID.eq.3)then
          ! Mixing length model with Fc from Collin et al
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using free air vertical diffusivity with stability function from Collin et al."
          endif;enddo
        elseif(KvFA_model_ID.eq.4)then
          ! Mixing length model with Fc from Piedelievre et al
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using free air vertical diffusivity with stability function from Piedelievre et al."
          endif;enddo
        else
          KvFA_model_ID = 3
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"    Using free air vertical diffusivity with stability function from Collin et al."
          endif;enddo
        endif


      endif

      if (useVarDiffH.or.useVarDiffV) then
        ! Check if we're using variable diffusivity, then get the constants
        read(10,'(a80)',iostat=ios,err=2010)linebuffer080
        read(linebuffer080,*,iostat=ioerr) vonKarman
        read(10,'(a80)',iostat=ios,err=2010)linebuffer080
        read(linebuffer080,*,iostat=ioerr) LambdaC
        read(10,'(a80)',iostat=ios,err=2010)linebuffer080
        read(linebuffer080,*,iostat=ioerr) RI_CRIT

        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)vonKarman
          write(outlog(io),*)LambdaC
          write(outlog(io),*)RI_CRIT
        endif;enddo

        ! We will want to reuse velocities on the metP grid for this module
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

!******************************************************************************

!******************************************************************************

      subroutine Allocate_VarDiff_Met

      use global_param,  only : &
         PI,useVarDiffH,useVarDiffV

      use io_data,       only : &
         nvar_User2d_static_XY,nvar_User3d_XYGs,nvar_User2d_XY,&
         nvar_User4d_XYZGs,nvar_User3d_XYZ

      use mesh,          only : &
         nxmax,nymax,nzmax

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet

      implicit none

      integer :: i

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"--------------------------------------------------"
        write(outlog(io),*)"---------- ALLOCATE_VARDIFF_MET ------------------"
        write(outlog(io),*)"--------------------------------------------------"
      endif;enddo

      ! HFS only allocate the bits needed for Kh vs Kv

      allocate(dVel_dz_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(du_dx_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(du_dy_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(dv_dx_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(dv_dy_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(dV_dz_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(SurfRoughLen_Met_sp(nx_submet,ny_submet))
      allocate(Ri_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(Ri_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(Khz_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(Khz_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(Kv_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(Kv_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(PBLH_meso_last_step_Met_sp(nx_submet,ny_submet))
      allocate(PBLH_meso_next_step_Met_sp(nx_submet,ny_submet))
      allocate(L_MonOb_meso_last_step_Met_sp(nx_submet,ny_submet))
      allocate(L_MonOb_meso_next_step_Met_sp(nx_submet,ny_submet))
      allocate(FricVel_meso_last_step_Met_sp(nx_submet,ny_submet))
      allocate(FricVel_meso_next_step_Met_sp(nx_submet,ny_submet))

      allocate(FricVel_meso_last_step_sp(nxmax,nymax))
      allocate(FricVel_meso_next_step_sp(nxmax,nymax))
      allocate(FricVel_ip(nxmax,nymax))

      allocate(Khz_meso_last_step_sp(nxmax,nymax,nzmax))
      allocate(Khz_meso_next_step_sp(nxmax,nymax,nzmax))
      allocate(Kv_meso_last_step_sp(nxmax,nymax,nzmax))
      allocate(Kv_meso_next_step_sp(nxmax,nymax,nzmax))

      allocate(vx_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(vy_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(vx_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(vy_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))

        ! Precalculate the LES term
      LES_L2ScaleCoeff = (KH_SmagC*KH_SmagC/(PI*PI))

      ! Set the start indecies
      indx_User2d_static_XY_VarDiff = nvar_User2d_static_XY
      indx_User2d_XY_VarDiff        = nvar_User2d_XY
      indx_User3d_XYGs_VarDiff      = nvar_User3d_XYGs
      indx_User3d_XYZ_VarDiff       = nvar_User3d_XYZ
      indx_User4d_XYZGs_VarDiff     = nvar_User4d_XYZGs

      ! Allocate output variables if needed
      if(.not.allocated(temp_2d_name_VarDiff))allocate(temp_2d_name_VarDiff(nvar_User2d_XY_VarDiff))
      if(.not.allocated(temp_2d_unit_VarDiff))allocate(temp_2d_unit_VarDiff(nvar_User2d_XY_VarDiff))
      if(.not.allocated(temp_2d_lname_VarDiff))allocate(temp_2d_lname_VarDiff(nvar_User2d_XY_VarDiff))
      if(.not.allocated(temp_2d_MissVal_VarDiff))allocate(temp_2d_MissVal_VarDiff(nvar_User2d_XY_VarDiff))
      if(.not.allocated(temp_2d_FillVal_VarDiff))allocate(temp_2d_FillVal_VarDiff(nvar_User2d_XY_VarDiff))

      if(.not.allocated(temp_3d_name_VarDiff))allocate(temp_3d_name_VarDiff(nvar_User3d_XYZ_VarDiff))
      if(.not.allocated(temp_3d_unit_VarDiff))allocate(temp_3d_unit_VarDiff(nvar_User3d_XYZ_VarDiff))
      if(.not.allocated(temp_3d_lname_VarDiff))allocate(temp_3d_lname_VarDiff(nvar_User3d_XYZ_VarDiff))
      if(.not.allocated(temp_3d_MissVal_VarDiff))allocate(temp_3d_MissVal_VarDiff(nvar_User3d_XYZ_VarDiff))
      if(.not.allocated(temp_3d_FillVal_VarDiff))allocate(temp_3d_FillVal_VarDiff(nvar_User3d_XYZ_VarDiff))

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
        temp_2d_lname_VarDiff(1) = "Planetary Boundary Layer Height"
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

      end subroutine Allocate_VarDiff_Met

!******************************************************************************

      subroutine Prep_output_VarDiff

      use global_param,  only : &
         useVarDiffH,KM_2_M,HR_2_S,DEG2RAD,M2PS_2_KM2PHR

      use mesh,          only : &
         nxmax,nymax,nzmax,lon_cc_pd,lat_cc_pd

      use Diffusion,     only : &
         kx,kz

      use Output_Vars,   only : &
         var_User2d_XY_name,var_User2d_XY_unit,var_User2d_XY_lname,&
         var_User2d_XY_MissVal,var_User2d_XY_FillVal,var_User2d_XY, &
         var_User3d_XYZ_name,var_User3d_XYZ_unit,var_User3d_XYZ_lname,&
         var_User3d_XYZ_MissVal,var_User3d_XYZ_FillVal,var_User3d_XYZ

      use MetReader,     only : &
         MR_dum2d_met,MR_dum2d_comp,MR_dum3d_compH,MR_dum3d_metP,MR_iMetStep_Now,&
           MR_Regrid_MetP_to_CompH,&
           MR_Regrid_Met2d_to_Comp2D

      use Atmosphere,    only : &
           solar_zenith

      use time_data,     only : &
         time,SimStartHour,Simtime_in_hours,BaseYear,useLeap

      implicit none

      integer :: i,ii,indx

      integer :: iii,jjj,kkk,hh,mm
      real(kind=ip) :: tmp
      real(kind=8) :: hour
      integer :: jday
      integer :: iyear,imonth,iday,idoy

      INTERFACE
        integer function HS_DayOfYear(HoursSince,byear,useLeaps)
          real(kind=8)               ::  HoursSince
          integer                    ::  byear
          logical                    ::  useLeaps
        end function HS_DayOfYear
        real(kind=8) function HS_HourOfDay(HoursSince,byear,useLeaps)
          real(kind=8)               ::  HoursSince
          integer                    ::  byear
          logical                    ::  useLeaps
        end function HS_HourOfDay
      END INTERFACE

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
          MR_dum2d_met = PBLH_meso_next_step_Met_sp/KM_2_M
          call MR_Regrid_Met2d_to_Comp2D
          var_User2d_XY(1:nxmax,1:nymax,indx) = MR_dum2d_comp(1:nxmax,1:nymax)
        elseif(i.eq.2)then
           ! Now resample onto computational grid
          MR_dum2d_met = FricVel_meso_next_step_Met_sp
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
          !var_User3d_XYZ(1:nxmax,1:nymax,1:nzmax,indx) = kx(1:nxmax,1:nymax,1:nzmax)*KM_2_M*KM_2_M/HR_2_S
          var_User3d_XYZ(1:nxmax,1:nymax,1:nzmax,indx) = kx(1:nxmax,1:nymax,1:nzmax)/M2PS_2_KM2PHR

!          jday = HS_DayOfYear(SimStartHour+time,BaseYear,useLeap)
!          hour = HS_HourOfDay(SimStartHour+time,BaseYear,useLeap)
!          hh   = floor(hour)
!          mm   = floor((hour-hh)*60.0_ip)
!          do iii=1,nxmax
!            do jjj=1,nymax
!              tmp = min(solar_zenith(lon_cc_pd(iii),lat_cc_pd(jjj),jday,hh,mm),90.0_ip)
!              var_User3d_XYZ(iii,jjj,1:nzmax,indx) = 361.0*cos(tmp*DEG2RAD)
!              var_User3d_XYZ(iii,jjj,1:nzmax,indx) = L_MonOb_meso_last_step_Met_sp(i)
!            enddo
!          enddo
           ! Now resample L_Mon onto computational grid
!          MR_dum2d_met = L_MonOb_meso_last_step_Met_sp
!          call MR_Regrid_Met2d_to_Comp2D
!          do kkk=1,nzmax
!            var_User3d_XYZ(1:nxmax,1:nymax,kkk,indx) = MR_dum2d_comp(1:nxmax,1:nymax)
!          enddo
        endif
        if(i.eq.ii+1)then
          ! Vertical diffusivity is already on the comp grid
          var_User3d_XYZ(1:nxmax,1:nymax,1:nzmax,indx) = kz(1:nxmax,1:nymax,1:nzmax)*KM_2_M*KM_2_M/HR_2_S
        endif
        if(i.eq.ii+2)then
           ! Ri just exists on the MetP grid for output
           ! Now resample onto computational grid
          MR_dum3d_metP = Ri_meso_next_step_MetP_sp
          call MR_Regrid_MetP_to_CompH(MR_iMetStep_Now)
          var_User3d_XYZ(1:nxmax,1:nymax,1:nzmax,indx) = MR_dum3d_compH(1:nxmax,1:nymax,1:nzmax)
        endif
      enddo

      end subroutine Prep_output_VarDiff

!******************************************************************************

      subroutine Deallocate_VarDiff_Met

      implicit none

      if(allocated(dVel_dz_MetP_sp))               deallocate(dVel_dz_MetP_sp)
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

      if(allocated(FricVel_meso_last_step_sp))     deallocate(FricVel_meso_last_step_sp)
      if(allocated(FricVel_meso_next_step_sp))     deallocate(FricVel_meso_next_step_sp)
      if(allocated(FricVel_ip))                    deallocate(FricVel_ip)

      if(allocated(Khz_meso_last_step_sp))         deallocate(Khz_meso_last_step_sp)
      if(allocated(Khz_meso_next_step_sp))         deallocate(Khz_meso_next_step_sp)
      if(allocated(Kv_meso_last_step_sp))          deallocate(Kv_meso_last_step_sp)
      if(allocated(Kv_meso_next_step_sp))          deallocate(Kv_meso_next_step_sp)

      if(allocated(vx_meso_last_step_MetP_sp))     deallocate(vx_meso_last_step_MetP_sp)
      if(allocated(vy_meso_last_step_MetP_sp))     deallocate(vy_meso_last_step_MetP_sp)
      if(allocated(vx_meso_next_step_MetP_sp))     deallocate(vx_meso_next_step_MetP_sp)
      if(allocated(vy_meso_next_step_MetP_sp))     deallocate(vy_meso_next_step_MetP_sp)

      end subroutine Deallocate_VarDiff_Met

!******************************************************************************
!******************************************************************************

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
      real(kind=ip) :: LES_LengthScale

      if(Kh_model_ID.eq.1)then
        Khz_meso_next_step_MetP_sp(:,:,:) = diffusivity_vert*M2PS_2_KM2PHR
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
        LES_LengthScale = MR_sigma_nz_submet(i,j)
  
          ! Diffusivity in km2/hr
        Khz_meso_next_step_MetP_sp(i,j,k) = real(KH_SmagC*LES_LengthScale*LES_TimeScale,kind=sp)
  
            enddo !k
          enddo !j
        enddo !i
      endif

      return

      end subroutine Eddy_diff

!******************************************************************************

      subroutine Calc_Vert_Diff(last_or_next)

      use global_param,  only : &
         KM_2_M,HR_2_S,KM_2_M,DEG2RAD

      use MetReader,     only : &
        nx_submet,ny_submet,np_fullmet,MR_geoH_metP_last,MR_geoH_metP_next,y_submet_sp,&
        MR_xy2ll_ylat,IsLatLon_MetGrid

      implicit none

      integer, intent(in) :: last_or_next

      integer :: i,j,k
      real(kind=ip) :: Ri_col(np_fullmet)
      real(kind=ip) :: z_col(np_fullmet)
      real(kind=ip) :: dv_dz_col(np_fullmet)
      real(kind=ip) :: FricVel
      real(kind=ip) :: Kv_col(np_fullmet)
      real(kind=ip) :: L_MonOb
      real(kind=ip) :: PBLz
      real(kind=ip) :: Phi
      real(kind=ip) :: PBL_profile_fac
      real(kind=ip) :: Kv_FreeAir, Kv_BL
      real(kind=ip) :: Lc
      real(kind=ip) :: EckF
      real(kind=ip) :: lat

      do i=1,nx_submet
        do j=1,ny_submet
          ! For the calculations, we need:
          !  PBLz, L_MonOb, FricVel, Ri, dv_dz, and z
          if(last_or_next.eq.0)then
              Ri_col(:) = real(Ri_meso_last_step_MetP_sp(i,j,:),kind=ip)    ! dimensionless
               z_col(:) = real(MR_geoH_metP_last(i,j,:),kind=ip)*KM_2_M     ! m
           dv_dz_col(:) = real(dV_dz_MetP_sp(i,j,:),kind=ip)                ! 1/s
                PBLz    = real(PBLH_meso_last_step_Met_sp(i,j),kind=ip)     ! m
             L_MonOb    = real(L_MonOb_meso_last_step_Met_sp(i,j),kind=ip)  ! m
             FricVel    = real(FricVel_meso_last_step_Met_sp(i,j),kind=ip)  ! m/s
          else
              Ri_col(:) = real(Ri_meso_next_step_MetP_sp(i,j,:),kind=ip)    ! dimensionless
               z_col(:) = real(MR_geoH_metP_next(i,j,:),kind=ip)*KM_2_M     ! m
           dv_dz_col(:) = real(dV_dz_MetP_sp(i,j,:),kind=ip)                ! 1/s
                PBLz    = real(PBLH_meso_next_step_Met_sp(i,j),kind=ip)     ! m
             L_MonOb    = real(L_MonOb_meso_next_step_Met_sp(i,j),kind=ip)  ! m
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
              Kv_FreeAir = Lc*Lc*abs(dv_dz_col(k))!*Fc(Ri_col_windp(k))

              if(z_col(k).lt.PBLz)then
                ! Within the PBL, use similarity theory
                  ! if PBL_exp=1; linear taper profile factor for Kv between 0 and PBL
                PBL_profile_fac = (1.0_sp-z_col(k)/PBLz)**PBL_exp

                if(IsLatLon_MetGrid)then
                  lat = y_submet_sp(j)
                else
                  lat = MR_xy2ll_ylat(i,j)
                endif
                lat = max(20.0_ip,abs(lat));
                EckF= 2.0_ip*7.292e-5_ip*sin(lat*DEG2RAD);
!                PBL_profile_fac = exp(-8.0_ip*EckF*z_col(k)/FricVel);

                Phi = Phi_WindShear_Similarity(z_col(k)/L_MonOb)
                ! Kv from similarity theory (Eq. 8.48 of Jacobson)
                Kv_BL = z_col(k)*vonKarman*FricVel*PBL_profile_fac/Phi
              endif

            endif
    
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

!******************************************************************************

      subroutine Set_VarDiffH_Meso(Load_MesoSteps,Interval_Frac)

      use mesh,          only : &
         nxmax,nymax,nzmax,dx,dy,IsLatLon,sigma_nz_pd

      use Diffusion,     only : &
         kx,ky

      use MetReader,     only : &
         MR_dum3d_compH,MR_vx_metP_last,MR_dum3d_metP,MR_dum3d2_metP,MR_iMetStep_Now,&
         MR_vy_metP_last,MR_vy_metP_next,MR_vx_metP_next,&
         nx_submet,ny_submet,&
           MR_DelMetP_Dx,&
           MR_DelMetP_Dy,&
           MR_Regrid_MetP_to_CompH

      implicit none

      logical      ,intent(in) :: Load_MesoSteps
      real(kind=ip),intent(in) :: Interval_Frac

      logical,save  :: first_time = .true.
      real(kind=sp) :: M_2_KM = 1.0e-3_sp

      integer :: i

      if(Load_MesoSteps)then
        if(first_time)then
          ! Need to fill _last_step_sp
          !  First fill next step so that outside this 'first_time' loop, the
          !  'next' can be copied to the 'last'
          ! Load U winds on MetP
          !ivar = 2 ! U winds
          !call MR_Read_3d_MetP_Variable(ivar,MR_iMetStep_Now)
          !vx_meso_next_step_MetP_sp = MR_dum3d_metP
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
          !ivar = 3 ! V winds
          !call MR_Read_3d_MetP_Variable(ivar,MR_iMetStep_Now)
          !vy_meso_next_step_MetP_sp = MR_dum3d_metP
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
        !ivar = 3 ! V winds
        !call MR_Read_3d_MetP_Variable(ivar,MR_iMetStep_Now+1)
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
                                              Interval_Frac
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

!******************************************************************************

      subroutine Set_VarDiffV_Meso(Load_MesoSteps,Interval_Frac)

      use global_param,  only : &
         M2PS_2_KM2PHR

      use mesh,          only : &
         nxmax,nymax,nzmax

      use Diffusion,     only : &
         kz

      use Atmosphere,    only : &
           Set_VirtPotenTemp

      use MetReader,     only : &
         MR_iMetStep_Now,MR_dum3d_MetP,MR_dum3d_compH,&
           MR_Regrid_MetP_to_CompH

      implicit none

      logical      ,intent(in) :: Load_MesoSteps
      real(kind=ip),intent(in) :: Interval_Frac

      logical,save :: first_time = .true.

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

      if(Load_MesoSteps)then
        if(first_time)then
          !  Populate values for the 'last' step
          call Set_VirtPotenTemp(0)
          call Calc_Ri(0)
          call Calc_Monin_Length(0)
          call Calc_SurfaceRoughnessLength
          call Calc_SurfaceFrictionVelocity(0)
          call Calc_PBLH(0)

          call Calc_Vert_Diff(0)
          MR_dum3d_MetP = Kv_meso_last_step_MetP_sp
          call MR_Regrid_MetP_to_CompH(MR_iMetStep_Now)
          Kv_meso_last_step_sp = MR_dum3d_compH

          ! The calls above with parameter 0 sets _last_ directly
          Ri_meso_next_step_MetP_sp     = Ri_meso_last_step_MetP_sp
          L_MonOb_meso_next_step_Met_sp = L_MonOb_meso_last_step_Met_sp
          FricVel_meso_next_step_Met_sp = FricVel_meso_last_step_Met_sp
          PBLH_meso_next_step_Met_sp    = PBLH_meso_last_step_Met_sp
          Kv_meso_next_step_sp          = Kv_meso_last_step_sp

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
        call Set_VirtPotenTemp(1)
        call Calc_Ri(1)                           ! sets Ri_meso_next_step_MetP_sp
        call Calc_Monin_Length(1)
        call Calc_SurfaceFrictionVelocity(1)      ! sets FricVel_meso_next_step_Met_sp
        call Calc_PBLH(1)
                                                  !  and L_MonOb_meso_next_step_Met_sp
        call Calc_Vert_Diff(1)
        MR_dum3d_MetP = Kv_meso_next_step_MetP_sp
        call MR_Regrid_MetP_to_CompH(MR_iMetStep_Now+1)
        Kv_meso_next_step_sp = MR_dum3d_compH

      endif

      kz(1:nxmax,1:nymax,1:nzmax) = real(Kv_meso_last_step_sp(:,:,:),kind=ip) + &
                                    real((Kv_meso_next_step_sp(:,:,:) - &
                                          Kv_meso_last_step_sp(:,:,:)),kind=ip) * &
                                              Interval_Frac * M2PS_2_KM2PHR
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

!******************************************************************************
!******************************************************************************

      subroutine Calc_Ri(last_or_next)

      use global_param,  only : &
         GRAV,KM_2_M,useMoistureVars,EPS_SMALL

      use Atmosphere,    only : &
         AirSH_meso_last_step_MetP_sp,AirSH_meso_next_step_MetP_sp,&
         AirTemp_meso_last_step_MetP_sp,AirTemp_meso_next_step_MetP_sp,&
         AirVPTemp_meso_last_step_MetP_sp,AirVPTemp_meso_next_step_MetP_sp,&
         R_GAS_DRYAIR,CP_AIR,R_GAS_WATVAP

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet,p_fullmet_sp,&
         x_submet_sp,y_submet_sp,&
         MR_geoH_metP_last,MR_geoH_MetP_next

      implicit none

      integer, intent(in) :: last_or_next

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
      real(kind=ip) :: mixrat
      real(kind=ip) :: del_z
      real(kind=ip) :: dudz,dvdz,dtdz
      real(kind=ip) :: dveldz2
      real(kind=ip) :: temp_term,mech_term
      real(kind=ip) :: Ri

      allocate(z(np_fullmet))
      allocate(u(np_fullmet))
      allocate(v(np_fullmet))
      allocate(p(np_fullmet))
      allocate(Tpoten(np_fullmet))

      refP = 1.0e5_ip   ! reference pressure for potential temperature

      p(1:np_fullmet) = p_fullmet_sp(1:np_fullmet)
      do i=1,nx_submet
        do j=1,ny_submet
          
          if(last_or_next.eq.0)then
            z(1:np_fullmet) = MR_geoH_metP_last(i,j,1:np_fullmet) * KM_2_M
            u(1:np_fullmet) = vx_meso_last_step_MetP_sp(i,j,1:np_fullmet)
            v(1:np_fullmet) = vy_meso_last_step_MetP_sp(i,j,1:np_fullmet)
            Tpoten(1:np_fullmet) = AirVPTemp_meso_last_step_MetP_sp(i,j,1:np_fullmet)
          else
            z(1:np_fullmet) = MR_geoH_MetP_next(i,j,1:np_fullmet) * KM_2_M
            u(1:np_fullmet) = vx_meso_next_step_MetP_sp(i,j,1:np_fullmet)
            v(1:np_fullmet) = vy_meso_next_step_MetP_sp(i,j,1:np_fullmet)
            Tpoten(1:np_fullmet) = AirVPTemp_meso_next_step_MetP_sp(i,j,1:np_fullmet)
          endif

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

            if(last_or_next.eq.0)then
              Ri_meso_last_step_MetP_sp(i,j,k) = Ri
            else
              Ri_meso_next_step_MetP_sp(i,j,k) = Ri
            endif

          enddo ! k
        enddo ! j
      enddo ! i

      end subroutine Calc_Ri

!******************************************************************************

      subroutine Calc_SurfaceRoughnessLength

      use MetReader,     only : &
         Met_var_IsAvailable,MR_iMetStep_Now,MR_dum2d_Met,&
         nx_submet,ny_submet, &
           MR_Read_2d_Met_Variable

      integer :: ivar

      ! Check if the windfile being used provides surface roughness
      ivar = 17 ! Surface_roughness_surface
      if(Met_var_IsAvailable(ivar))then
        ! Surface roughness is provided, read it from the met file
        call MR_Read_2d_Met_Variable(ivar,MR_iMetStep_Now)
        SurfRoughLen_Met_sp(1:nx_submet,1:ny_submet)  = MR_dum2d_Met(1:nx_submet,1:ny_submet)

      !elseif(useLandCover)then
      !  ! Set SurfRoughLen_Met_sp from Land use classification

      !  if(LandCover_Format.eq.1)then
      !    ! 1 degree resolution
      !    if(LC_grid)!0       Water
      !    !1       Broadleaf Evergreen Forest
      !    !2       Coniferous Evergreen Forest and Woodland
      !    !3       High Latitude Deciduous Forest and Woodland
      !    !4       Tundra
      !    !5       Mixed Coniferous Forest and Woodland
      !    !6       Wooded Grassland
      !    !7       Grassland
      !    !8       Bare Ground
      !    !9       Shrubs and Bare Ground
      !    !10      Cultivated Crops
      !    !11      Broadleaf Deciduous Forest and Woodland
      !    !12      Data Unavailable
      !  elseif(LandCover_Format.eq.2)then
      !    ! 8 km resolution
      !    !0       Water (and Goode's interrupted space)
      !    !1       Evergreen Needleleaf Forest
      !    !2       Evergreen Broadleaf Foreset
      !    !3       Deciduous Needleleaf Forest
      !    !4       Deciduous Broadleaf Forest
      !    !5       Mixed Forest
      !    !6       Woodland
      !    !7       Wooded Grassland
      !    !8       Closed Shrubland
      !    !9       Open Shrubland
      !    !10      Grassland
      !    !11      Cropland
      !    !12      Bare Ground
      !    !13      Urban and Built-up
      !  elseif(LandCover_Format.eq.3)then
      !    ! 1 km resolution
      !    !0       Water (and Goode's interrupted space)
      !    !1       Evergreen Needleleaf Forest
      !    !2       Evergreen Broadleaf Foreset
      !    !3       Deciduous Needleleaf Forest
      !    !4       Deciduous Broadleaf Forest
      !    !5       Mixed Forest
      !    !6       Woodland
      !    !7       Wooded Grassland
      !    !8       Closed Shrubland
      !    !9       Open Shrubland
      !    !10      Grassland
      !    !11      Cropland
      !    !12      Bare Ground
      !    !13      Urban and Built-up
      !  endif
      else
        ! Set SurfRoughLen_Met_sp by assumption
          SurfRoughLen_Met_sp(1:nx_submet,1:ny_submet)  = 0.1_sp
      endif

      end subroutine Calc_SurfaceRoughnessLength

!******************************************************************************

      subroutine Calc_SurfaceFrictionVelocity(last_or_next)

      use global_param,  only : &
         EPS_SMALL,KM_2_M,MPS_2_KMPHR

      use MetReader,     only : &
         Met_var_IsAvailable,MR_iMetStep_Now,np_fullmet,MR_iMetStep_Now,MR_Topo_met,&
         MR_geoH_metP_last,MR_geoH_metP_next,MR_dum2d_Met,nx_submet,ny_submet,&
           MR_Read_2d_Met_Variable

      implicit none

      integer, intent(in) :: last_or_next

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

      ! Check if the windfile being used provides friction velocity
      ivar = 13 ! Friction_velocity_surface
      if(Met_var_IsAvailable(ivar))then
        ! Friction velocity is provided, read it from the met file
        if(last_or_next.eq.0)then
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

      if(.not.Met_var_IsAvailable(ivar).or.FV_override)then
        ! friction velocity is not provided by the met file (or is not valid)
        ! Calculate it ourselves
        ! Get surface friction velocity using Panofsky/Dutton p376

        allocate(SurfVelx_meso_Met_sp(nx_submet,ny_submet))
        allocate(SurfVely_meso_Met_sp(nx_submet,ny_submet))
        allocate(SurfVelh_meso_Met_sp(nx_submet,ny_submet))

        ! First check if we can get the 10m wind (both x and y)
        if(Met_var_IsAvailable(11).and.Met_var_IsAvailable(12))then
          if(last_or_next.eq.0)then
            call MR_Read_2d_Met_Variable(11,MR_iMetStep_Now)
              SurfVelx_meso_Met_sp = MR_dum2d_Met
!            if(maxval(MR_dum2d_Met(:,:)).lt.EPS_SMALL)then
!              ! variable was present in file, but filled with nonsense
!              FV_override = .true.
!            endif
            call MR_Read_2d_Met_Variable(12,MR_iMetStep_Now)
              SurfVely_meso_Met_sp = MR_dum2d_Met
          else
            call MR_Read_2d_Met_Variable(11,MR_iMetStep_Now+1)
              SurfVelx_meso_Met_sp = MR_dum2d_Met
!            if(maxval(MR_dum2d_Met(:,:)).lt.EPS_SMALL)then
!              ! variable was present in file, but filled with nonsense
!              FV_override = .true.
!            endif
            call MR_Read_2d_Met_Variable(12,MR_iMetStep_Now+1)
              SurfVely_meso_Met_sp = MR_dum2d_Met
          endif
          SurfVelh_meso_Met_sp = 10.0_sp
        else
          ! If the 10m winds are not available, then use the lower levels of the 3d winds
          do i=1,nx_submet
            do j=1,ny_submet
              z0 = SurfRoughLen_Met_sp(i,j)
              if(last_or_next.eq.0)then
                do k=1,np_fullmet
                  if(MR_geoH_metP_last(i,j,k)*KM_2_M.gt.1000.0_ip*MR_Topo_met(i,j)+z0)then
                    exit
                  endif
                enddo
                SurfVelx_meso_Met_sp(i,j) = vx_meso_last_step_MetP_sp(i,j,k)
                SurfVely_meso_Met_sp(i,j) = vy_meso_last_step_MetP_sp(i,j,k)
                SurfVelh_meso_Met_sp(i,j) = (MR_geoH_metP_last(i,j,k)-MR_Topo_met(i,j))*KM_2_M
              else
                do k=1,np_fullmet
                  if(MR_geoH_metP_next(i,j,k)*KM_2_M.gt.1000.0_ip*MR_Topo_met(i,j)+z0)then
                    exit
                  endif
                enddo
                SurfVelx_meso_Met_sp(i,j) = vx_meso_next_step_MetP_sp(i,j,k)
                SurfVely_meso_Met_sp(i,j) = vy_meso_next_step_MetP_sp(i,j,k)
                SurfVelh_meso_Met_sp(i,j) = (MR_geoH_metP_next(i,j,k)-MR_Topo_met(i,j))*KM_2_M
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
            if(last_or_next.eq.0)then
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

!******************************************************************************

      subroutine Calc_PBLH(last_or_next)

      use global_param,  only : &
         EPS_SMALL,KM_2_M,DEG2RAD

      use Atmosphere,    only : &
         AirVPTemp_meso_last_step_MetP_sp,AirVPTemp_meso_next_step_MetP_sp

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet,Met_var_IsAvailable,&
         MR_geoH_metP_last,MR_geoH_metP_next,MR_dum2d_Met,MR_iMetStep_Now,&
         MR_Have_LL_mapping,MR_xy2ll_ylat,IsLatLon_MetGrid,y_submet_sp,&
           MR_Read_2d_Met_Variable,MR_Set_LL_mapping

      implicit none

      integer, intent(in) :: last_or_next

      integer :: ivar
      integer :: i,j,k,k_L,kk
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

      ! Check if the windfile being used provides PBLH
      ivar = 10 ! Planetary Boundary Level Height
      if(Met_var_IsAvailable(ivar))then
        ! PBLH is provided, read it from the met file
        if(last_or_next.eq.0)then
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

      PBL_override = .true.

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
            if(last_or_next.eq.0)then
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
              lat = y_submet_sp(j)
            else
              lat = MR_xy2ll_ylat(i,j)
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

            if(last_or_next.eq.0)then
              PBLH_meso_last_step_Met_sp(i,j) = real(PBLz,kind=sp)
            else
              PBLH_meso_next_step_Met_sp(i,j) = real(PBLz,kind=sp)
            endif
          enddo
        enddo
        deallocate(PBLtmp)
      endif


      end subroutine Calc_PBLH

!!******************************************************************************


!******************************************************************************

      subroutine Calc_Monin_Length(last_or_next)

      use global_param,  only : &
         EPS_SMALL,KM_2_M

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet,Met_var_IsAvailable,&
         MR_geoH_metP_last,MR_geoH_metP_next,MR_dum2d_Met,MR_iMetStep_Now,&
           MR_Read_2d_Met_Variable

      implicit none

      integer, intent(in) :: last_or_next

      integer :: ivar
      integer :: i,j,k,k_L
      real(kind=ip) :: denom,tmp
      real(kind=ip) :: Ri
      real(kind=ip) :: Ri_col(np_fullmet)
      real(kind=ip) :: z_col(np_fullmet)
      real(kind=ip) :: L_MonOb

      ! Get Monin-Obukhov length from the
      ! Businger-Dyer-Pandolfo empirical result 
      ! using z and Ri at k=2
        ! Eq 6.7.1 and 6.7 2 of "Atmospheric Turbulence";
        ! Panofsky and Dutton,1984
        ! also Eq 11.24 of "Introduction to Micrometeorology";
        ! Arya, 1988
      do i=1,nx_submet
        do j=1,ny_submet
          if(last_or_next.eq.0)then
            Ri_col(:) = Ri_meso_last_step_MetP_sp(i,j,:)
            z_col(:)  = MR_geoH_metP_last(i,j,:)*KM_2_M
          else
            Ri_col(:) = Ri_meso_next_step_MetP_sp(i,j,:)
            z_col(:)  = MR_geoH_metP_next(i,j,:)*KM_2_M
          endif

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

          if(last_or_next.eq.0)then
            L_MonOb_meso_last_step_Met_sp(i,j) = real(L_MonOb,kind=sp)
          else
            L_MonOb_meso_next_step_Met_sp(i,j) = real(L_MonOb,kind=sp)
          endif
        enddo
      enddo

      end subroutine Calc_Monin_Length
!
!!******************************************************************************

      function Fc(Ri)

      ! Stability function for vertical diffusion in the free atmosphere above the PBL
      ! Originally from Collins et al, NCAR TN-464, 2004
      ! http://www.cesm.ucar.edu/models/atm-cam/docs/description/description.pdf

      implicit none

      real(kind=ip) :: Fc ! dimensionless
      real(kind=ip) :: Ri ! dimensionless

      if(Ri.ge.0.0_ip)then
          ! Eq. 4.465  : Stable atmosphere
        Fc = 1.0_ip/(1.0_ip+10.0_ip*Ri*(1.0_ip+8.0_ip*Ri))
      else
          ! Eq. 4.464  : unstable
        Fc = sqrt(1.0_ip-18.0_ip*Ri)
      endif

      return

      end function Fc

!!******************************************************************************

      function Fc_PMB(Ri,ml,z)

      ! Stability function for vertical diffusion
      ! see Piedelievre, Jean Philippe, Lue Musson-Genon, François Bompay, 1990:
      ! MEDIA—An Eulerian Model of Atmospheric Dispersion: First Validation on
      ! the Chernobyl Release. J. Appl. Meteor., 29, 1205–1220.
      ! doi: http://dx.doi.org/10.1175/1520-0450(1990)029<1205:MEMOAD>2.0.CO;2 
      ! This is the model used by MDLP0 (also used in EMERRAUDE or PERIDOT)

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

      if(Ri.ge.0.0_ip)then
        Fc_PMB = 1.0_ip/(1.0_ip+3.0_ip*b*Ri*sqrt(1.0_ip+d*Ri))
      else
        Fc_PMB = 1.0_ip/(1.0_ip+3.0_ip*b*c*(ml*ml*sqrt(Ri)/(z*z*f)))
      endif

      return

      end function Fc_PMB

!!******************************************************************************

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

      !MixLen = 1.0_ip/(1.0_ip/(z*vonKarman) + 1.0_ip/LambdaC)
      MixLen = (z*vonKarman)/(1.0_ip+(z*vonKarman/LambdaC))

      return

      end function MixLen

!******************************************************************************

      function Phi_WindShear_Similarity(z_on_L)

      ! Generalized wind shear similarity function of z/L
      ! Coefficients and exponents are module variables

      implicit none

      real(kind=ip) :: Phi_WindShear_Similarity

      real(kind=ip) :: z_on_L

      if(z_on_L.le.0.0_ip)then
          ! Unstable
        Phi_WindShear_Similarity = (1.0_ip + phi_gamma*z_on_L)**phi_alpha
      else
          ! Stable
        Phi_WindShear_Similarity = (1.0_ip + phi_beta*z_on_L)
      endif

      return

      end function Phi_WindShear_Similarity

!******************************************************************************

      function Psi_WindShear_Similarity(z_on_L)

      ! Generalized wind shear similarity function of z/L (integral of phi)

      use global_param,  only : &
         PI

      implicit none

      real(kind=ip) :: Psi_WindShear_Similarity

      real(kind=ip) :: z_on_L

      real(kind=ip) :: x,tmp1,tmp2

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

!******************************************************************************

      end module Diffusivity_Variable
