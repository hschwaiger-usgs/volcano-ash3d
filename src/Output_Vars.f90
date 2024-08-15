!##############################################################################
!
!  Output_Vars module
!
!  This module manages data, settings and the calculation of output variables.
!  
!
!    area(lat, lon)
!    depotime(lat, lon)
!    ash_arrival_time(lat, lon)
!    depothickFin(lat, lon)
!    
!    
!    depothick(t, lat, lon)
!    ashcon_max(t, lat, lon)
!    cloud_height(t, lat, lon)
!    cloud_load(t, lat, lon) 
!    cloud_mask(t, lat, lon)
!    cloud_bottom(t, lat, lon)
!    
!    radar_reflectivity(t, z, lat, lon)
!    depocon(t, bn, lat, lon)
!    ashcon(t, bn, z, lat, lon)

!
!      subroutine Allocate_Output_Vars
!      subroutine Allocate_NTime
!      subroutine Allocate_Profile
!      subroutine Allocate_Output_UserVars
!      subroutine Deallocate_Output_Vars
!      subroutine Deallocate_NTime
!      subroutine Deallocate_Profile
!      subroutine Deallocate_Output_UserVars
!      subroutine Set_OutVar_ContourLevel
!      subroutine AshThicknessCalculator
!      subroutine AshTotalCalculator
!      subroutine dbZCalculator      
!      subroutine ConcentrationCalculator      
!      subroutine CloudAreaCalculator      
!      subroutine Gen_Output_Vars
!      subroutine Calc_AshVol_Aloft
!      subroutine Calc_vprofile
!      subroutine Calc_AshVol_Deposit
!      subroutine Calc_AshVol_Outflow
!      subroutine FirstAsh
!
!##############################################################################

      module Output_Vars

      use precis_param

      use io_units

      use global_param,  only : &
         PI, M_2_MM, KM2_2_M2, KG_2_MG, KM3_2_M3, EPS_SMALL

      use time_data,     only : &
         time_native

      implicit none

        ! Set everything to private by default
      private

        ! Publicly available subroutines/functions
      public Allocate_Output_Vars,      &
             Deallocate_Output_Vars,    &
             Allocate_NTime,            &
             Deallocate_NTime,          &
             Allocate_Profile,          &
             Deallocate_Profile,        &
             Allocate_Output_UserVars,  &
             Deallocate_Output_UserVars,&
             Set_OutVar_ContourLevel,   &
             AshTotalCalculator,        &
             dbZCalculator,             &
             Gen_Output_Vars,           &
             Calc_AshVol_Aloft,         &
             Calc_vprofile,             &
             Calc_AshVol_Deposit,       &
             Calc_AshVol_Outflow,       &
             FirstAsh

        ! PRODUCT ID: KML var ID:: variable                     : units  : dims        : Output kml file name
        !----------------------------------------------------------------------------------------------------------
        ! iprod =  1:ivar =     :: full concentration array     : kg/km3 : (x,y,z,t,g) :
        ! iprod =  2:ivar =     :: deposit granularity          : kg/m2  : (x,y,t,g)   :
        ! iprod =  3:ivar =   7 :: deposit                      : mm     : (x,y,t)     : deposit_thickness_mm.kml
        ! iprod =  4:ivar =   8 :: deposit (NWS)                : inches : (x,y,t)     : deposit_thickness_inches.kml
        ! iprod =  5:ivar =  *7 :: deposit final                : mm     : (x,y)       :
        ! iprod =  6:ivar =  *8 :: deposit final (NWS)          : inches : (x,y)       :
        ! iprod =  7:ivar =   9 :: ashfall arrival time         : hours  : (x,y)       : ashfall_arrivaltimes_hours.kml
        ! iprod =  8:ivar =     :: ashfall at airports          : mm     : (pt,tn)     : ash_arrivaltimes_airports.kml

        ! iprod =  9:ivar =   1 :: ash-cloud concentration      : mg/m3  : (x,y,t)     : CloudConcentration.kml
        ! iprod = 10:ivar =   2 :: ash-cloud height (top)       : km     : (x,y,t)     : CloudHeight.kml
        ! iprod = 11:ivar =   3 :: ash-cloud height (bot)       : km     : (x,y,t)     : CloudBottom.kml
        ! iprod = 12:ivar =   4 :: ash-cloud load               : T/km2  : (x,y,t)     : CloudLoad.kml
        ! iprod = 13:ivar =   6 :: ash-cloud radar reflectivity : dBz    : (x,y,t)     : reflectivity.kml
        ! iprod = 14:ivar =   5 :: ash-cloud arrival time       : hours  : (x,y)       : cloud_arrivaltimes_hours.kml

        ! iprod = 15:ivar =  10 :: topography                   : km     : (x,y)       : Topography.kml
        ! iprod = 16:ivar =     :: profiles concentration?      :

        ! Publicly available variables
        !  Variable threshold values
        !  The deposit and cloudload thresholds are special in that they are used
        !  for delineating the deposit and cloud masks respectively.  These masks
        !  are written to the output netcdf file if post-processing required
        !  thresholding (if flooded contours are to be suppressed below the
        !  threshold, for example).
      real(kind=ip),public    :: DEPO_THRESH           = 1.0e-2_ip  ! threshold deposit thickness (mm)
      real(kind=ip),public    :: CLOUDLOAD_THRESH      = 2.0e-1_ip  ! threshold cloud load (t/km2)
                                       ! 0.2 T/km2 is roughly the detection
                                       ! limit of Pavolonis's SEVIRI satellite retrievals

      real(kind=ip),public    :: DEPRATE_THRESH        = 1.0e-2_ip  ! threshold deposition rate (mm/hr)
      real(kind=ip),public    :: CLOUDCON_THRESH       = 1.0e-3_ip  ! threshold cloud concentration (kg/km3) for output
      real(kind=ip),public    :: CLOUDCON_GRID_THRESH  = 1.0e-7_ip  ! threshold cloud concentration (kg/km3) for subgrid

      real(kind=ip),public    :: THICKNESS_THRESH = 1.0e-2_ip  ! threshold thickness for start of deposition (mm)
      real(kind=ip),public    :: DBZ_THRESH       =-2.0e+1_ip  ! threshold dbZ

      ! These are the initialized values
      real(kind=op),public    :: DepositThickness_FillValue   = -9999.0_op
      real(kind=op),public    :: MaxConcentration_FillValue   = -9999.0_op
      real(kind=op),public    :: DepArrivalTime_FillValue     = -9999.0_op
      real(kind=op),public    :: CloudArrivalTime_FillValue   = -9999.0_op
      real(kind=op),public    :: CloudLoad_FillValue          = -9999.0_op
      real(kind=op),public    :: MaxHeight_FillValue          = -9999.0_op
      real(kind=op),public    :: MinHeight_FillValue          = -9999.0_op
      real(kind=op),public    :: dbZCol_FillValue             = -9999.0_op
      real(kind=op),public    :: pr_ash_FillValue             = -9999.0_op

      logical,public          :: Calculated_Cloud_Load
      logical,public          :: Calculated_AshThickness
        ! Set this parameter if you want to include velocities in the output file
      logical,public          :: useWindVars  = .false.

        ! Set this to true if you want the extra output variables defined in the
        ! optional modules
      logical, parameter,public :: USE_OPTMOD_VARS  = .true.

        ! Set this parameter if you want to export additional variables
        ! to the netcdf file
      !logical, parameter :: USE_ADDITIONAL_VARS  = .false.
      !logical, parameter :: log_2d1_tmp          = .false.
      !logical, parameter :: log_3d1_tmp          = .false.
      !logical, parameter :: log_4d1_tmp          = .false.

        ! This variable is set to true indicating that the output file should include
        ! the standard derived variables in the output file.
      logical,public :: useOutprodVars = .true.

        ! This variable will be set to false if you do not want raw concentration
        ! values exported (only derived products and deposits) if indicated
        ! in the input file on block 5/line 15
        ! yes 2   # Write out 3-D ash concentration at specified times? / [output code: 1=2d+concen,2=2d only]
        ! Can also be reset in RESETPARAMS
      logical,public :: useRestartVars = .false.

      real(kind=ip),public :: CloudArea                ! area of ash cloud at a given time
      real(kind=ip),public :: LoadVal(5)               ! 5 threshold values for area calculations
      real(kind=ip),public :: CloudLoadArea(5)         ! Corresponding areas where ash cloud exceeds LoadVal(i)
      real(kind=ip),public :: DepositAreaCovered       ! area covered by ash deposit

      integer      ,public :: iplotpref  = 0           ! Used in post-processing to specify plotting package

      ! Contour colors and levels
      logical                                     ,public:: ContourFilled = .false. ! T if using filled contours, F if lines
      integer                                     ,public:: nConLev
      integer,parameter                           ,public:: CONTOUR_MAXCURVES  = 40
      integer,parameter                           ,public:: CONTOUR_MAXPOINTS  = 1000
        ! User-specified contour interval and colors
      logical                                     ,public:: Con_Cust   = .false.    ! T if using a custom set of contours
      integer                                     ,public:: Con_Cust_N

#ifdef USEPOINTERS
      real(kind=ip),dimension(:)    ,pointer ,public:: ContourLev         => null()
      real(kind=ip),dimension(:,:,:),pointer ,public:: ContourDataX       => null() ! x curve data with dims: ilev, icurve, ipnt
      real(kind=ip),dimension(:,:,:),pointer ,public:: ContourDataY       => null() ! x curve data with dims: ilev, icurve, ipnt
      integer      ,dimension(:)    ,pointer ,public:: ContourDataNcurves => null() ! num of curves for each level (some = 0)
      integer      ,dimension(:,:)  ,pointer ,public:: ContourDataNpoints => null() ! num of pts for ilev and icurve
        ! User-specified contour interval and colors
      integer      ,dimension(:,:)  ,pointer,public:: Con_Cust_RGB        => null()
      real(kind=ip),dimension(:)    ,pointer,public:: Con_Cust_Lev        => null()
#else
      real(kind=ip),dimension(:)    ,allocatable  ,public:: ContourLev
      real(kind=ip),dimension(:,:,:),allocatable  ,public:: ContourDataX        ! x curve data with dims: ilev, icurve, ipnt
      real(kind=ip),dimension(:,:,:),allocatable  ,public:: ContourDataY        ! x curve data with dims: ilev, icurve, ipnt
      integer      ,dimension(:)    ,allocatable  ,public:: ContourDataNcurves  ! num of curves for each level (some = 0)
      integer      ,dimension(:,:)  ,allocatable  ,public:: ContourDataNpoints  ! num of pts for ilev and icurve
        ! User-specified contour interval and colors
      integer      ,dimension(:,:)  ,allocatable  ,public:: Con_Cust_RGB
      real(kind=ip),dimension(:)    ,allocatable  ,public:: Con_Cust_Lev
#endif



        ! Fixed size arrays for output products
      integer,parameter                           ,public:: Con_DepThick_mm_N   = 10
      integer      ,dimension(Con_DepThick_mm_N,3),public:: Con_DepThick_mm_RGB
      real(kind=ip),dimension(Con_DepThick_mm_N)  ,public:: Con_DepThick_mm_Lev
      integer,parameter                           ,public:: Con_DepThick_in_N   = 5
      integer      ,dimension(Con_DepThick_in_N,3),public:: Con_DepThick_in_RGB
      real(kind=ip),dimension(Con_DepThick_in_N)  ,public:: Con_DepThick_in_Lev
      integer,parameter                           ,public:: Con_DepTime_N       = 9
      integer      ,dimension(Con_DepTime_N,3)    ,public:: Con_DepTime_RGB
      real(kind=ip),dimension(Con_DepTime_N)      ,public:: Con_DepTime_Lev
      integer,parameter                           ,public:: Con_CloudCon_N      = 9
      integer      ,dimension(Con_CloudCon_N,3)   ,public:: Con_CloudCon_RGB
      real(kind=ip),dimension(Con_CloudCon_N)     ,public:: Con_CloudCon_Lev
      integer,parameter                           ,public:: Con_CloudTop_N      = 9
      integer      ,dimension(Con_CloudTop_N,3)   ,public:: Con_CloudTop_RGB
      real(kind=ip),dimension(Con_CloudTop_N)     ,public:: Con_CloudTop_Lev
      integer,parameter                           ,public:: Con_CloudBot_N      = 9
      integer      ,dimension(Con_CloudBot_N,3)   ,public:: Con_CloudBot_RGB
      real(kind=ip),dimension(Con_CloudBot_N)     ,public:: Con_CloudBot_Lev
      integer,parameter                           ,public:: Con_CloudLoad_N     = 9
      integer      ,dimension(Con_CloudLoad_N,3)  ,public:: Con_CloudLoad_RGB
      real(kind=ip),dimension(Con_CloudLoad_N)    ,public:: Con_CloudLoad_Lev
      integer,parameter                           ,public:: Con_CloudRef_N      = 9
      integer      ,dimension(Con_CloudRef_N,3)   ,public:: Con_CloudRef_RGB
      real(kind=ip),dimension(Con_CloudRef_N)     ,public:: Con_CloudRef_Lev
      integer,parameter                           ,public:: Con_CloudTime_N     = 9
      integer      ,dimension(Con_CloudTime_N,3)  ,public:: Con_CloudTime_RGB
      real(kind=ip),dimension(Con_CloudTime_N)    ,public:: Con_CloudTime_Lev

#ifdef USEPOINTERS
        ! 2-D variables (in x,y)
      logical,       dimension(:,:)  ,pointer,public :: Mask_Cloud       => null()
      logical,       dimension(:,:)  ,pointer,public :: Mask_Deposit     => null()
      real(kind=ip), dimension(:,:)  ,pointer,public :: DepositThickness => null() ! accumulated ash thickness on ground (mm)
      real(kind=ip), dimension(:,:)  ,pointer,public :: MaxConcentration => null() ! max concentration in the cloud at any i,j node (mg/m3)
      real(kind=dp), dimension(:,:)  ,pointer,public :: DepArrivalTime   => null() ! (hours)
      real(kind=dp), dimension(:,:)  ,pointer,public :: CloudArrivalTime => null() ! (hours)
      real(kind=ip), dimension(:,:)  ,pointer,public :: CloudLoad        => null() ! Ash load in cloud, (tonnes/km2)
      real(kind=ip), dimension(:,:)  ,pointer,public :: MaxHeight        => null() ! maximum cloud height (km)
      real(kind=ip), dimension(:,:)  ,pointer,public :: MinHeight        => null() ! cloud bottom height (km)
      real(kind=ip), dimension(:,:)  ,pointer,public :: dbZCol           => null() ! max reflectivity in a vertical column (dB)
      real(kind=op), dimension(:,:,:),pointer,public :: ashcon_tot       => null() ! Total ash concentration (3d)
      real(kind=ip), dimension(:,:,:),pointer,public :: pr_ash           => null() ! concentration profile
      real(kind=ip), dimension(:,:)  ,pointer        :: CloudLoadLast    => null() ! Ash load at last time step, (tonnes/km2)

        ! 3-D variables
        !   (in x,y,z)
      real(kind=ip), dimension(:,:,:),pointer,public :: dbZ => null()               ! radar reflectivty at time t (dbZ)

        ! User-defined 2-D static variables (in x,y)
      character(len=30), dimension(:),  pointer,public :: var_User2d_static_XY_name    => null()
      character(len=30), dimension(:),  pointer,public :: var_User2d_static_XY_unit    => null()
      character(len=30), dimension(:),  pointer,public :: var_User2d_static_XY_lname   => null()
      real(kind=op),     dimension(:),  pointer,public :: var_User2d_static_XY_MissVal => null()
      real(kind=op),     dimension(:),  pointer,public :: var_User2d_static_XY_FillVal => null()
      real(kind=op), dimension(:,:,:),  pointer,public :: var_User2d_static_XY         => null()
        ! User-defined 2-D variables (in x,y)
      character(len=30), dimension(:),  pointer,public :: var_User2d_XY_name           => null()
      character(len=30), dimension(:),  pointer,public :: var_User2d_XY_unit           => null()
      character(len=30), dimension(:),  pointer,public :: var_User2d_XY_lname          => null()
      real(kind=op),     dimension(:),  pointer,public :: var_User2d_XY_MissVal        => null()
      real(kind=op),     dimension(:),  pointer,public :: var_User2d_XY_FillVal        => null()
      real(kind=op), dimension(:,:,:),  pointer,public :: var_User2d_XY                => null()
        ! User-defined 3-D variables (in x,y,gs)
      character(len=30), dimension(:),  pointer,public :: var_User3d_XYGs_name         => null()
      character(len=30), dimension(:),  pointer,public :: var_User3d_XYGs_unit         => null()
      character(len=30), dimension(:),  pointer,public :: var_User3d_XYGs_lname        => null()
      real(kind=op),     dimension(:),  pointer,public :: var_User3d_XYGs_MissVal      => null()
      real(kind=op),     dimension(:),  pointer,public :: var_User3d_XYGs_FillVal      => null()
      real(kind=op), dimension(:,:,:,:),pointer,public :: var_User3d_XYGs              => null()
        ! User-defined 3-D variables (in x,y,z)
      character(len=30), dimension(:),  pointer,public :: var_User3d_XYZ_name          => null()
      character(len=30), dimension(:),  pointer,public :: var_User3d_XYZ_unit          => null()
      character(len=30), dimension(:),  pointer,public :: var_User3d_XYZ_lname         => null()
      real(kind=op),     dimension(:),  pointer,public :: var_User3d_XYZ_MissVal       => null()
      real(kind=op),     dimension(:),  pointer,public :: var_User3d_XYZ_FillVal       => null()
      real(kind=ip), dimension(:,:,:,:),pointer,public :: var_User3d_XYZ               => null()
        ! User-defined 4-D variables (in x,y,z,gs)
      character(len=30), dimension(:),    pointer,public :: var_User4d_XYZGs_name      => null()
      character(len=30), dimension(:),    pointer,public :: var_User4d_XYZGs_unit      => null()
      character(len=30), dimension(:),    pointer,public :: var_User4d_XYZGs_lname     => null()
      real(kind=op),     dimension(:),    pointer,public :: var_User4d_XYZGs_MissVal   => null()
      real(kind=op),     dimension(:),    pointer,public :: var_User4d_XYZGs_FillVal   => null()
      real(kind=op), dimension(:,:,:,:,:),pointer,public :: var_User4d_XYZGs           => null()
#else
        ! 2-D variables (in x,y)
      logical,       dimension(:,:)  ,allocatable,public :: Mask_Cloud
      logical,       dimension(:,:)  ,allocatable,public :: Mask_Deposit
      real(kind=ip), dimension(:,:)  ,allocatable,public :: DepositThickness   ! accumulated ash thickness on ground in mm (x,y)
      real(kind=ip), dimension(:,:)  ,allocatable,public :: MaxConcentration   ! max concentration in the cloud at any i,j node (mg/m3)
      real(kind=dp), dimension(:,:)  ,allocatable,public :: DepArrivalTime     ! (hours)
      real(kind=dp), dimension(:,:)  ,allocatable,public :: CloudArrivalTime   ! (hours)
      real(kind=ip), dimension(:,:)  ,allocatable,public :: CloudLoad          ! Ash load in cloud, tonnes/km2
      real(kind=ip), dimension(:,:)  ,allocatable,public :: MaxHeight          ! maximum cloud height
      real(kind=ip), dimension(:,:)  ,allocatable,public :: MinHeight          ! cloud bottom height
      real(kind=ip), dimension(:,:)  ,allocatable,public :: dbZCol             ! max reflectivity in a vertical column
      real(kind=op), dimension(:,:,:),allocatable,public :: ashcon_tot         ! Total ash concentration (3d)
      real(kind=ip), dimension(:,:,:),allocatable,public :: pr_ash             ! concentration profile
      real(kind=ip), dimension(:,:)  ,allocatable        :: CloudLoadLast      ! Ash load at last time step, tonnes/km2

        ! 3-D variables
        !   (in x,y,z)
      real(kind=ip), dimension(:,:,:),allocatable,public :: dbZ                ! radar reflectivty at time t (dbZ)

        ! User-defined 2-D static variables (in x,y)
      character(len=30), dimension(:),    allocatable,public :: var_User2d_static_XY_name
      character(len=30), dimension(:),    allocatable,public :: var_User2d_static_XY_unit
      character(len=30), dimension(:),    allocatable,public :: var_User2d_static_XY_lname
      real(kind=op),     dimension(:),    allocatable,public :: var_User2d_static_XY_MissVal
      real(kind=op),     dimension(:),    allocatable,public :: var_User2d_static_XY_FillVal
      real(kind=op), dimension(:,:,:),    allocatable,public :: var_User2d_static_XY
        ! User-defined 2-D variables (in x,y)
      character(len=30), dimension(:),    allocatable,public :: var_User2d_XY_name
      character(len=30), dimension(:),    allocatable,public :: var_User2d_XY_unit
      character(len=30), dimension(:),    allocatable,public :: var_User2d_XY_lname
      real(kind=op),     dimension(:),    allocatable,public :: var_User2d_XY_MissVal
      real(kind=op),     dimension(:),    allocatable,public :: var_User2d_XY_FillVal
      real(kind=op), dimension(:,:,:),    allocatable,public :: var_User2d_XY
        ! User-defined 3-D variables (in x,y,gs)
      character(len=30), dimension(:),    allocatable,public :: var_User3d_XYGs_name
      character(len=30), dimension(:),    allocatable,public :: var_User3d_XYGs_unit
      character(len=30), dimension(:),    allocatable,public :: var_User3d_XYGs_lname
      real(kind=op),     dimension(:),    allocatable,public :: var_User3d_XYGs_MissVal
      real(kind=op),     dimension(:),    allocatable,public :: var_User3d_XYGs_FillVal
      real(kind=op), dimension(:,:,:,:),  allocatable,public :: var_User3d_XYGs
        ! User-defined 3-D variables (in x,y,z)
      character(len=30), dimension(:),    allocatable,public :: var_User3d_XYZ_name
      character(len=30), dimension(:),    allocatable,public :: var_User3d_XYZ_unit
      character(len=30), dimension(:),    allocatable,public :: var_User3d_XYZ_lname
      real(kind=op),     dimension(:),    allocatable,public :: var_User3d_XYZ_MissVal
      real(kind=op),     dimension(:),    allocatable,public :: var_User3d_XYZ_FillVal
      real(kind=ip), dimension(:,:,:,:),  allocatable,public :: var_User3d_XYZ
        ! User-defined 4-D variables (in x,y,z,gs)
      character(len=30), dimension(:),    allocatable,public :: var_User4d_XYZGs_name
      character(len=30), dimension(:),    allocatable,public :: var_User4d_XYZGs_unit
      character(len=30), dimension(:),    allocatable,public :: var_User4d_XYZGs_lname
      real(kind=op),     dimension(:),    allocatable,public :: var_User4d_XYZGs_MissVal
      real(kind=op),     dimension(:),    allocatable,public :: var_User4d_XYZGs_FillVal
      real(kind=op), dimension(:,:,:,:,:),allocatable,public :: var_User4d_XYZGs
#endif

      character(len=30)                       ,public :: Extra2dVarName
      real(kind=ip),dimension(:,:),allocatable,public :: Extra2dVar

      contains
      !------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Allocate_Output_Vars
!
!  Called from: alloc_arrays
!  Arguments:
!    none
!
!  This subroutine allocated the arrays that store the output products, mostly
!  2-d variables.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Allocate_Output_Vars

      use mesh,          only : &
         nxmax,nymax,nzmax

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Allocate_Output_Vars"
      endif;enddo

      allocate(Mask_Cloud(nxmax,nymax))
      Mask_Cloud = .false.
      allocate(Mask_Deposit(nxmax,nymax))
      Mask_Deposit = .false.
      allocate(DepositThickness(nxmax,nymax))
      DepositThickness = 0.0_ip
      allocate(MaxConcentration(nxmax,nymax))
      MaxConcentration = 0.0_ip
      allocate(DepArrivalTime(nxmax,nymax))
      DepArrivalTime = DepArrivalTime_FillValue
      allocate(CloudArrivalTime(nxmax,nymax))                    ! time of arrival of ash cloud
      CloudArrivalTime = CloudArrivalTime_FillValue
      allocate(CloudLoad(nxmax,nymax))
      CloudLoad    = 0.0_ip
      allocate(CloudLoadLast(nxmax,nymax))
      CloudLoadLast = 0.0_ip
      allocate(MaxHeight(nxmax,nymax))                           ! maximum height (top) of ash cloud
      MaxHeight = 0.0_ip
      allocate(MinHeight(nxmax,nymax))                           ! minimum height (bottom) of ash cloud
      MinHeight = 0.0_ip
      allocate(dbZCol(nxmax,nymax))                              ! reflectivity in a column of nodes
      dbZCol = 0.0_ip
      ! ashcon_tot is allocated/deallocated as needed since it can be a bit big
      allocate(dbZ(nxmax,nymax,nzmax))                              ! radar reflectivity (dbZ)
      dbZ = 0.0_ip

      end subroutine Allocate_Output_Vars

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Allocate_NTime(nt)
!
!  Called from: Ash3d.F90
!  Arguments:
!    nt  = total anticipated time steps
!
!  This subroutine allocated the variable for storing the time steps for the
!  simulation.  These are the full steps, not just the output or log steps.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Allocate_NTime(nt)

      integer,intent(in) :: nt

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Allocate_NTime"
      endif;enddo

      allocate(time_native(nt)); time_native = 0.0_dp

      end subroutine Allocate_NTime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Allocate_Profile(nz,nt,nv)
!
!  Called from: Ash3d.F90
!  Arguments:
!    nz  = z length of vertical profile
!    nt  = total anticipated time steps
!    nv  = number of vertical profiles
!
!  This subroutine allocates the variable for storing profile data, if profile
!  data output is requested.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Allocate_Profile(nz,nt,nv)

      integer,intent(in) :: nz
      integer,intent(in) :: nt
      integer,intent(in) :: nv

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Allocate_Profile"
      endif;enddo

      if(nv.gt.0)then
        allocate(pr_ash(nz,nt,nv))                       ! vertical ash profile
        pr_ash = 0.0_op
      endif

      end subroutine Allocate_Profile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Allocate_Output_UserVars(nx,ny,nz,ns)
!
!  Called from: Ash3d.F90
!  Arguments:
!    nx  = x length of computational grid
!    ny  = y length of computational grid
!    nz  = z length of computational grid
!    ns  = number of bin characterizing grainsize/species
!
!  This subroutine allocates all the variables needed for user-defined output
!  products.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Allocate_Output_UserVars(nx,ny,nz,ns)

      use io_data,       only : &
         nvar_User2d_static_XY,nvar_User2d_XY,&
         nvar_User3d_XYGs,nvar_User3d_XYZ,    &
         nvar_User4d_XYZGs

      integer,intent(in) :: nx
      integer,intent(in) :: ny
      integer,intent(in) :: nz
      integer,intent(in) :: ns

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Allocate_Output_UserVars"
      endif;enddo

        ! User-defined 2-D static variables (in x,y)
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Allocating var_User2d_static_XY: ",nx,ny,nvar_User2d_static_XY
      endif;enddo
      allocate(var_User2d_static_XY_name(nvar_User2d_static_XY));   var_User2d_static_XY_name   = ''
      allocate(var_User2d_static_XY_unit(nvar_User2d_static_XY));   var_User2d_static_XY_unit   = ''
      allocate(var_User2d_static_XY_lname(nvar_User2d_static_XY));  var_User2d_static_XY_lname  = ''
      allocate(var_User2d_static_XY_MissVal(nvar_User2d_static_XY));var_User2d_static_XY_MissVal= 0.0_op
      allocate(var_User2d_static_XY_FillVal(nvar_User2d_static_XY));var_User2d_static_XY_FillVal= 0.0_op
      allocate(var_User2d_static_XY(nx,ny,nvar_User2d_static_XY));  var_User2d_static_XY        = 0.0_op
        ! User-defined 2-D variables (in x,y)
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Allocating var_User2d_XY: ",nx,ny,nvar_User2d_XY
      endif;enddo
      allocate(var_User2d_XY_name(nvar_User2d_XY));    var_User2d_XY_name   = ''
      allocate(var_User2d_XY_unit(nvar_User2d_XY));    var_User2d_XY_unit   = ''
      allocate(var_User2d_XY_lname(nvar_User2d_XY));   var_User2d_XY_lname  = ''
      allocate(var_User2d_XY_MissVal(nvar_User2d_XY)); var_User2d_XY_MissVal= 0.0_op
      allocate(var_User2d_XY_FillVal(nvar_User2d_XY)); var_User2d_XY_FillVal= 0.0_op
      allocate(var_User2d_XY(nx,ny,nvar_User2d_XY));   var_User2d_XY        = 0.0_op
        ! User-defined 3-D variables (in x,y,gs)
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Allocating var_User3d_XYGs: ",nx,ny,ns,nvar_User3d_XYGs
      endif;enddo
      allocate(var_User3d_XYGs_name(nvar_User3d_XYGs));     var_User3d_XYGs_name   = ''
      allocate(var_User3d_XYGs_unit(nvar_User3d_XYGs));     var_User3d_XYGs_unit   = ''
      allocate(var_User3d_XYGs_lname(nvar_User3d_XYGs));    var_User3d_XYGs_lname  = ''
      allocate(var_User3d_XYGs_MissVal(nvar_User3d_XYGs));  var_User3d_XYGs_MissVal= 0.0_op
      allocate(var_User3d_XYGs_FillVal(nvar_User3d_XYGs));  var_User3d_XYGs_FillVal= 0.0_op
      allocate(var_User3d_XYGs(nx,ny,ns,nvar_User3d_XYGs)); var_User3d_XYGs        = 0.0_op
        ! User-defined 3-D variables (in x,y,z)
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Allocating var_User3d_XYZ: ",nx,ny,nz,nvar_User3d_XYZ
      endif;enddo
      allocate(var_User3d_XYZ_name(nvar_User3d_XYZ));     var_User3d_XYZ_name   = ''
      allocate(var_User3d_XYZ_unit(nvar_User3d_XYZ));     var_User3d_XYZ_unit   = ''
      allocate(var_User3d_XYZ_lname(nvar_User3d_XYZ));    var_User3d_XYZ_lname  = ''
      allocate(var_User3d_XYZ_MissVal(nvar_User3d_XYZ));  var_User3d_XYZ_MissVal= 0.0_op
      allocate(var_User3d_XYZ_FillVal(nvar_User3d_XYZ));  var_User3d_XYZ_FillVal= 0.0_op
      allocate(var_User3d_XYZ(nx,ny,nz,nvar_User3d_XYZ)); var_User3d_XYZ        = 0.0_op
        ! User-defined 4-D variables (in x,y,z,gs)
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Allocating var_User4d_XYZGs: ",nx,ny,ns,nvar_User4d_XYZGs
      endif;enddo
      allocate(var_User4d_XYZGs_name(nvar_User4d_XYZGs));        var_User4d_XYZGs_name   = ''
      allocate(var_User4d_XYZGs_unit(nvar_User4d_XYZGs));        var_User4d_XYZGs_unit   = ''
      allocate(var_User4d_XYZGs_lname(nvar_User4d_XYZGs));       var_User4d_XYZGs_lname  = ''
      allocate(var_User4d_XYZGs_MissVal(nvar_User4d_XYZGs));     var_User4d_XYZGs_MissVal= 0.0_op
      allocate(var_User4d_XYZGs_FillVal(nvar_User4d_XYZGs));     var_User4d_XYZGs_FillVal= 0.0_op
      allocate(var_User4d_XYZGs(nx,ny,nz,ns,nvar_User4d_XYZGs)); var_User4d_XYZGs        = 0.0_op

      end subroutine Allocate_Output_UserVars

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Deallocate_Output_Vars()
!
!  Called from: dealloc_arrays
!  Arguments:
!    none
!
!  This subroutine deallocates the variables allocated in Allocate_Output_Vars.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Deallocate_Output_Vars()

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Deallocate_Output_Vars"
      endif;enddo

#ifdef USEPOINTERS
      if(associated(Mask_Cloud))       deallocate(Mask_Cloud)
      if(associated(Mask_Deposit))     deallocate(Mask_Deposit)
      if(associated(DepositThickness)) deallocate(DepositThickness)
      if(associated(MaxConcentration)) deallocate(MaxConcentration)
      if(associated(DepArrivalTime))   deallocate(DepArrivalTime)
      if(associated(CloudArrivalTime)) deallocate(CloudArrivalTime)
      if(associated(CloudLoad))        deallocate(CloudLoad)
      if(associated(CloudLoadLast))    deallocate(CloudLoadLast)
      if(associated(MaxHeight))        deallocate(MaxHeight)
      if(associated(MinHeight))        deallocate(MinHeight)
      if(associated(dbZCol))           deallocate(dbZCol)
      if(associated(dbZ))              deallocate(dbZ)
      if(associated(Con_Cust_RGB))     deallocate(Con_Cust_RGB)
      if(associated(Con_Cust_Lev))     deallocate(Con_Cust_Lev)
#else
      if(allocated(Mask_Cloud))       deallocate(Mask_Cloud)
      if(allocated(Mask_Deposit))     deallocate(Mask_Deposit)
      if(allocated(DepositThickness)) deallocate(DepositThickness)
      if(allocated(MaxConcentration)) deallocate(MaxConcentration)
      if(allocated(DepArrivalTime))   deallocate(DepArrivalTime)
      if(allocated(CloudArrivalTime)) deallocate(CloudArrivalTime)
      if(allocated(CloudLoad))        deallocate(CloudLoad)
      if(allocated(CloudLoadLast))    deallocate(CloudLoadLast)
      if(allocated(MaxHeight))        deallocate(MaxHeight)
      if(allocated(MinHeight))        deallocate(MinHeight)
      if(allocated(dbZCol))           deallocate(dbZCol)
      if(allocated(dbZ))              deallocate(dbZ)
      if(allocated(Con_Cust_RGB))     deallocate(Con_Cust_RGB)
      if(allocated(Con_Cust_Lev))     deallocate(Con_Cust_Lev)
#endif

      end subroutine Deallocate_Output_Vars

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Deallocate_NTime
!
!  Called from: dealloc_arrays
!  Arguments:
!    none
!
!  This subroutine deallocates the variables allocated in Allocate_NTime.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Deallocate_NTime

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Deallocate_NTime"
      endif;enddo

#ifdef USEPOINTERS
      if(associated(time_native))deallocate(time_native)
#else
      if(allocated(time_native))deallocate(time_native)
#endif

      end subroutine Deallocate_NTime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Deallocate_Profile
!
!  Called from: dealloc_arrays
!  Arguments:
!    none
!
!  This subroutine deallocates the variables allocated in Allocate_Profile.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Deallocate_Profile

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Deallocate_Profile"
      endif;enddo

#ifdef USEPOINTERS
      if(associated(pr_ash))deallocate(pr_ash)                       ! vertical ash profile
#else
      if(allocated(pr_ash))deallocate(pr_ash)                       ! vertical ash profile
#endif

      end subroutine Deallocate_Profile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Deallocate_Output_UserVars()
!
!  Called from: dealloc_arrays
!  Arguments:
!    none
!
!  This subroutine deallocates the variables allocated in Allocate_Output_UserVars.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Deallocate_Output_UserVars()

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Deallocate_Output_UserVars"
      endif;enddo

#ifdef USEPOINTERS
      if(associated(var_User2d_static_XY_name))    deallocate(var_User2d_static_XY_name)
      if(associated(var_User2d_static_XY_unit))    deallocate(var_User2d_static_XY_unit)
      if(associated(var_User2d_static_XY_lname))   deallocate(var_User2d_static_XY_lname)
      if(associated(var_User2d_static_XY_MissVal)) deallocate(var_User2d_static_XY_MissVal)
      if(associated(var_User2d_static_XY_FillVal)) deallocate(var_User2d_static_XY_FillVal)
      if(associated(var_User2d_static_XY))         deallocate(var_User2d_static_XY)
      if(associated(var_User2d_XY_name))           deallocate(var_User2d_XY_name)
      if(associated(var_User2d_XY_unit))           deallocate(var_User2d_XY_unit)
      if(associated(var_User2d_XY_lname))          deallocate(var_User2d_XY_lname)
      if(associated(var_User2d_XY_MissVal))        deallocate(var_User2d_XY_MissVal)
      if(associated(var_User2d_XY_FillVal))        deallocate(var_User2d_XY_FillVal)
      if(associated(var_User2d_XY))                deallocate(var_User2d_XY)
      if(associated(var_User3d_XYGs_name))         deallocate(var_User3d_XYGs_name)
      if(associated(var_User3d_XYGs_unit))         deallocate(var_User3d_XYGs_unit)
      if(associated(var_User3d_XYGs_lname))        deallocate(var_User3d_XYGs_lname)
      if(associated(var_User3d_XYGs_MissVal))      deallocate(var_User3d_XYGs_MissVal)
      if(associated(var_User3d_XYGs_FillVal))      deallocate(var_User3d_XYGs_FillVal)
      if(associated(var_User3d_XYGs))              deallocate(var_User3d_XYGs)
      if(associated(var_User3d_XYZ_name))          deallocate(var_User3d_XYZ_name)
      if(associated(var_User3d_XYZ_unit))          deallocate(var_User3d_XYZ_unit)
      if(associated(var_User3d_XYZ_lname))         deallocate(var_User3d_XYZ_lname)
      if(associated(var_User3d_XYZ_MissVal))       deallocate(var_User3d_XYZ_MissVal)
      if(associated(var_User3d_XYZ_FillVal))       deallocate(var_User3d_XYZ_FillVal)
      if(associated(var_User3d_XYZ))               deallocate(var_User3d_XYZ)
      if(associated(var_User4d_XYZGs_name))        deallocate(var_User4d_XYZGs_name)
      if(associated(var_User4d_XYZGs_unit))        deallocate(var_User4d_XYZGs_unit)
      if(associated(var_User4d_XYZGs_lname))       deallocate(var_User4d_XYZGs_lname)
      if(associated(var_User4d_XYZGs_MissVal))     deallocate(var_User4d_XYZGs_MissVal)
      if(associated(var_User4d_XYZGs_FillVal))     deallocate(var_User4d_XYZGs_FillVal)
      if(associated(var_User4d_XYZGs))             deallocate(var_User4d_XYZGs)
#else
      if(allocated(var_User2d_static_XY_name))    deallocate(var_User2d_static_XY_name)
      if(allocated(var_User2d_static_XY_unit))    deallocate(var_User2d_static_XY_unit)
      if(allocated(var_User2d_static_XY_lname))   deallocate(var_User2d_static_XY_lname)
      if(allocated(var_User2d_static_XY_MissVal)) deallocate(var_User2d_static_XY_MissVal)
      if(allocated(var_User2d_static_XY_FillVal)) deallocate(var_User2d_static_XY_FillVal)
      if(allocated(var_User2d_static_XY))         deallocate(var_User2d_static_XY)
      if(allocated(var_User2d_XY_name))           deallocate(var_User2d_XY_name)
      if(allocated(var_User2d_XY_unit))           deallocate(var_User2d_XY_unit)
      if(allocated(var_User2d_XY_lname))          deallocate(var_User2d_XY_lname)
      if(allocated(var_User2d_XY_MissVal))        deallocate(var_User2d_XY_MissVal)
      if(allocated(var_User2d_XY_FillVal))        deallocate(var_User2d_XY_FillVal)
      if(allocated(var_User2d_XY))                deallocate(var_User2d_XY)
      if(allocated(var_User3d_XYGs_name))         deallocate(var_User3d_XYGs_name)
      if(allocated(var_User3d_XYGs_unit))         deallocate(var_User3d_XYGs_unit)
      if(allocated(var_User3d_XYGs_lname))        deallocate(var_User3d_XYGs_lname)
      if(allocated(var_User3d_XYGs_MissVal))      deallocate(var_User3d_XYGs_MissVal)
      if(allocated(var_User3d_XYGs_FillVal))      deallocate(var_User3d_XYGs_FillVal)
      if(allocated(var_User3d_XYGs))              deallocate(var_User3d_XYGs)
      if(allocated(var_User3d_XYZ_name))          deallocate(var_User3d_XYZ_name)
      if(allocated(var_User3d_XYZ_unit))          deallocate(var_User3d_XYZ_unit)
      if(allocated(var_User3d_XYZ_lname))         deallocate(var_User3d_XYZ_lname)
      if(allocated(var_User3d_XYZ_MissVal))       deallocate(var_User3d_XYZ_MissVal)
      if(allocated(var_User3d_XYZ_FillVal))       deallocate(var_User3d_XYZ_FillVal)
      if(allocated(var_User3d_XYZ))               deallocate(var_User3d_XYZ)
      if(allocated(var_User4d_XYZGs_name))        deallocate(var_User4d_XYZGs_name)
      if(allocated(var_User4d_XYZGs_unit))        deallocate(var_User4d_XYZGs_unit)
      if(allocated(var_User4d_XYZGs_lname))       deallocate(var_User4d_XYZGs_lname)
      if(allocated(var_User4d_XYZGs_MissVal))     deallocate(var_User4d_XYZGs_MissVal)
      if(allocated(var_User4d_XYZGs_FillVal))     deallocate(var_User4d_XYZGs_FillVal)
      if(allocated(var_User4d_XYZGs))             deallocate(var_User4d_XYZGs)
#endif

      end subroutine Deallocate_Output_UserVars

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Set_OutVar_ContourLevel
!
!  Called from: NC_Read_Output_Products
!  Arguments:
!    none
!
!  This subroutine fills several variables that are needed when plotting
!  maps in Ash3d_PostProc.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Set_OutVar_ContourLevel

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Set_OutVar_ContourLevel"
      endif;enddo

      ! Set contour levels and colors
      !  Deposit Thickness (mm)
      ! Recall that Con_DepThick_mm_N   = 10
      Con_DepThick_mm_Lev = (/0.01_ip, 0.03_ip, 0.1_ip, 0.3_ip, 1.0_ip, &
                             3.0_ip, 10.0_ip, 30.0_ip, 100.0_ip, 300.0_ip /)
      Con_DepThick_mm_RGB( 1,1:3) = (/ 214,222,105 /)
      Con_DepThick_mm_RGB( 2,1:3) = (/ 249,167,113 /)
      Con_DepThick_mm_RGB( 3,1:3) = (/ 128,  0,128 /)
      Con_DepThick_mm_RGB( 4,1:3) = (/   0,  0,255 /)
      Con_DepThick_mm_RGB( 5,1:3) = (/   0,128,255 /)
      Con_DepThick_mm_RGB( 6,1:3) = (/   0,255,128 /)
      Con_DepThick_mm_RGB( 7,1:3) = (/ 195,195,  0 /)
      Con_DepThick_mm_RGB( 8,1:3) = (/ 255,128,  0 /)
      Con_DepThick_mm_RGB( 9,1:3) = (/ 255,  0,  0 /)
      Con_DepThick_mm_RGB(10,1:3) = (/ 128,  0,  0 /)

      !  Deposit Thickness (in)
      ! Recall that Con_DepThick_in_N   = 5
      Con_DepThick_in_Lev = (/0.004_ip, 0.031_ip, 0.25_ip, 1.0_ip,   4.0_ip /)
      Con_DepThick_in_RGB( 1,1:3) = (/ 255,  0,  0 /)
      Con_DepThick_in_RGB( 2,1:3) = (/   0,  0,255 /)
      Con_DepThick_in_RGB( 3,1:3) = (/   0,183,255 /)
      Con_DepThick_in_RGB( 4,1:3) = (/ 255,  0,255 /)
      Con_DepThick_in_RGB( 5,1:3) = (/   0, 51, 51 /)

      !  Deposit (ashfall) arrival time (hours)
      ! Recall that Con_DepTime_N   = 9
      Con_DepTime_Lev = (/0.0_ip, 3.0_ip, 6.0_ip, 9.0_ip, 12.0_ip, 15.0_ip, &
                         18.0_ip, 24.0_ip, 36.0_ip  /)
      Con_DepTime_RGB( 1,1:3) = (/ 255,  0,  0 /)  ! ff0000
      Con_DepTime_RGB( 2,1:3) = (/ 255,128,255 /)  ! ff8000
      Con_DepTime_RGB( 3,1:3) = (/ 255,255,  0 /)  ! ffff00
      Con_DepTime_RGB( 4,1:3) = (/ 128,255,128 /)  ! 80ff80
      Con_DepTime_RGB( 5,1:3) = (/   0,255,255 /)  ! 00ffff
      Con_DepTime_RGB( 6,1:3) = (/   0,128,255 /)  ! 0080ff
      Con_DepTime_RGB( 7,1:3) = (/   0,  0,255 /)  ! 0000ff
      Con_DepTime_RGB( 8,1:3) = (/   0,  0,128 /)  ! 000080
      Con_DepTime_RGB( 9,1:3) = (/   0,  0,  0 /)  ! 000000


      !  Maximum ash concentration (mg/m3)
      ! Recall that Con_CloudCon_N   = 9
      Con_CloudCon_Lev = (/0.1_ip, 0.3_ip, 1.0_ip, 3.0_ip, 10.0_ip, 30.0_ip,&
                         100.0_ip, 300.0_ip, 1000.0_ip /)
        ! This is what is in the KML file, but it is too light for a white background
      Con_CloudCon_RGB( 1,1:3) = (/ 204,204,255 /)  ! ffe5e5
      Con_CloudCon_RGB( 2,1:3) = (/ 178,178,255 /)  ! ffcccc
      Con_CloudCon_RGB( 3,1:3) = (/ 153,153,255 /)  ! ffb2b2
      Con_CloudCon_RGB( 4,1:3) = (/ 255,153,255 /)  ! ff99ff
      Con_CloudCon_RGB( 5,1:3) = (/ 255,126,255 /)  ! ff7fff
      Con_CloudCon_RGB( 6,1:3) = (/ 255,102,255 /)  ! ff66ff
      Con_CloudCon_RGB( 7,1:3) = (/ 255, 76,255 /)  ! ff4cff
      Con_CloudCon_RGB( 8,1:3) = (/ 255, 51,255 /)  ! ff33ff
      Con_CloudCon_RGB( 9,1:3) = (/ 255, 51,255 /)  ! ff33ff
        ! Using the same colors as depotime
      Con_CloudCon_RGB( 1:9,1:3) = Con_DepTime_RGB( 1:9,1:3)

      ! Cloud Height Top (km)
      ! Recall that Con_CloudTop_N   = 9
      Con_CloudTop_Lev = (/ 0.24_ip, 3.0_ip, 6.0_ip, 10.0_ip, 13.0_ip, 16.0_ip,&
                         20.0_ip, 25.0_ip, 30.0_ip /)
      Con_CloudTop_RGB( 1,1:3) = (/ 128,  0,128 /)  ! 800080
      Con_CloudTop_RGB( 2,1:3) = (/ 255,  0,  0 /)  ! ff0000
      Con_CloudTop_RGB( 3,1:3) = (/ 255,128,255 /)  ! ff8000
      Con_CloudTop_RGB( 4,1:3) = (/ 255,255,  0 /)  ! ffff00
      Con_CloudTop_RGB( 5,1:3) = (/ 128,255,128 /)  ! 80ff80
      Con_CloudTop_RGB( 6,1:3) = (/   0,255,255 /)  ! 00ffff
      Con_CloudTop_RGB( 7,1:3) = (/   0,128,255 /)  ! 0080ff
      Con_CloudTop_RGB( 8,1:3) = (/   0,  0,255 /)  ! 0000ff
      Con_CloudTop_RGB( 9,1:3) = (/   0,  0,128 /)  ! 000080

      ! Cloud Height Bot (km)
      ! Recall that Con_CloudBot_N   = 9
      Con_CloudBot_Lev(1:Con_CloudBot_N) = Con_CloudTop_Lev(1:Con_CloudBot_N)
      Con_CloudBot_RGB( 1:9,1:3) = Con_CloudTop_RGB( 1:9,1:3)

      ! Cloud Load (T/km2)
      ! Recall that Con_CloudLoad_N   = 9
      Con_CloudLoad_Lev = (/ 0.2_ip, 1.0_ip, 2.0_ip, 5.0_ip, 10.0_ip, 30.0_ip,&
                         100.0_ip, 300.0_ip, 1000.0_ip /)
      Con_CloudLoad_RGB( 1:9,1:3) = Con_CloudTop_RGB( 1:9,1:3)

      ! Cloud Reflectivity (dBz)
      ! Recall that Con_CloudRef_N   = 9
      Con_CloudRef_Lev = (/-20.0_ip, -10.0_ip, 0.0_ip, 10.0_ip, 20.0_ip, 30.0_ip,&
                         40.0_ip, 50.0_ip, 60.0_ip /)
      Con_CloudRef_RGB( 1:9,1:3) = Con_CloudTop_RGB( 1:9,1:3)

      ! Cloud Arrival Time (hours)
      ! Recall that Con_CloudTime_N   = 9
      Con_CloudTime_Lev(1:Con_CloudTime_N) = Con_DepTime_Lev(1:Con_CloudTime_N)
      Con_CloudTime_RGB( 1:9,1:3) = Con_CloudTop_RGB( 1:9,1:3)

      end subroutine Set_OutVar_ContourLevel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  AshThicknessCalculator
!
!  Called from: Gen_Output_Vars and FirstAsh
!  Arguments:
!    none
!
!  This subroutine calculates the thickness (mm) of the ash deposit by summing
!  the mass of all the tephra bins and converting to thickness via DepositDensity.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine AshThicknessCalculator

      use mesh,          only : &
         nxmax,nymax,dz_vec_pd,sigma_nz_pd

      use solution,      only : &
         DepositGranularity

      use Tephra,        only : &
         DepositDensity,n_gs_max

      integer :: i,j

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine AshThicknessCalculator"
      endif;enddo

      !calculate deposit thickness in mm, and area covered
      Mask_Deposit(:,:)         = .false.
      DepositThickness(:,:)     = 0.0_ip
      DepositAreaCovered        = 0.0_ip

      if(n_gs_max.gt.0)then
        do i=1,nxmax
          do j=1,nymax
            DepositThickness(i,j) = sum(DepositGranularity(i,j,1:n_gs_max)) * &  ! in kg/km^3
                                    dz_vec_pd(1)                            / &  ! convert to kg/km^2
                                    KM2_2_M2                                / &  ! from kg/km^2 to kg/m^2
                                    DepositDensity                          * &  ! from kg/m^2 to m
                                    M_2_MM                                       ! from m to mm
            if (DepositThickness(i,j).gt.DEPO_THRESH)then
              Mask_Deposit(i,j) = .true.
              DepositAreaCovered = DepositAreaCovered + sigma_nz_pd(i,j,1)
            endif
          enddo
        enddo
      endif

      Calculated_AshThickness = .true.

      end subroutine AshThicknessCalculator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  AshTotalCalculator
!
!  Called from: output_results
!  Arguments:
!    none
!
!  This subroutine calculates the total concentration of ash at every point in
!  the domain by summing the concentrations of each tephra bin.  This is only
!  used in the writing out of 3d ASCII and BINARY values.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine AshTotalCalculator

      use mesh,          only : &
         nxmax,nymax,nzmax,ts1

      use solution,      only : &
         concen_pd

      use Tephra,        only : &
         n_gs_max

      integer :: isize

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine AshTotalCalculator"
      endif;enddo

      if(n_gs_max.gt.0)then
        ashcon_tot = 0.0_op
        do isize=1,n_gs_max
          ashcon_tot(1:nxmax,1:nymax,1:nzmax) =  &
           ashcon_tot(1:nxmax,1:nymax,1:nzmax) + &
           real(concen_pd(1:nxmax,1:nymax,1:nzmax,isize,ts1),kind=op)
        enddo
      endif

      end subroutine AshTotalCalculator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  dbZCalculator
!
!  Called from: NC_append_to_netcdf
!  Arguments:
!    none
!
!  This subroutine calculates the radar reflectivity factor (Z) in dBZ.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine dbZCalculator

      use global_param,  only : &
          EPS_TINY,M_2_MM

      use mesh,          only : &
         nzmax,ts1

      use solution,      only : &
         concen_pd,imin,imax,jmin,jmax,kmin,kmax

      use Tephra,        only : &
         n_gs_max,Tephra_gsdiam,Tephra_rho_m

      integer :: i,j,k,isize
      real(kind=ip) :: NumDens          !number densities (#/m3) of particles
      real(kind=ip) :: zcol             !z value of cell
      real(kind=ip) :: tmp

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine dbZCalculator"
      endif;enddo

      dbZCol(:,:) = dbZCol_FillValue
      dbZ(:,:,:)  = dbZCol_FillValue

      ! Calculate particle collision rate between two particle sizes
      ! Note: this requires that only two particle sizes be used as input

      if(n_gs_max.gt.0)then
        do i=imin,imax
          do j=jmin,jmax
            if (CloudLoad(i,j).lt.CLOUDLOAD_THRESH) cycle
            do k=kmin,kmax
              zcol = 0.0_ip
              do isize=1,n_gs_max
                ! convert concentration (kg/km3) to number density (#/m3)
                NumDens = concen_pd(i,j,k,isize,ts1) / &
                            (Tephra_rho_m(isize)*PI*Tephra_gsdiam(isize)**3.0_ip/6.0_ip) / &
                            KM3_2_M3                                  ! particles/m3
                zcol    = zcol + NumDens*(Tephra_gsdiam(isize)*M_2_MM)**6.0_ip
              enddo
              if(zcol.lt.EPS_TINY)then
                dbZ(i,j,k) = dbZCol_FillValue
              else
                tmp = 10.0_ip*log10(zcol)
                if(tmp.gt.DBZ_THRESH)then
                  dbZ(i,j,k) = tmp
                else
                  dbZ(i,j,k) = dbZCol_FillValue
                endif
              endif
            enddo
            dbZCol(i,j) = maxval(dbZ(i,j,1:nzmax))
          enddo
        enddo
      endif

      end subroutine dbZCalculator      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  ConcentrationCalculator
!
!  Called from: Gen_Output_Vars and FirstAsh
!  Arguments:
!    none
!
!  This subroutine several output variables derived from the airborne ash
!  concentration: CloudLoad (vertically integrated value), 
!  MaxConcentration (peak value of all z), Min and Max cloud height, the
!  cloud mask (logical variable) and the CloudArea of the cloud mask.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ConcentrationCalculator

      use mesh,          only : &
         nxmax,nymax,nzmax,dz_vec_pd,z_cc_pd,ts1,sigma_nz_pd

      use solution,      only : &
         concen_pd,imin,imax,jmin,jmax

      use Tephra,        only : &
         n_gs_max

      integer :: i,j,k
      real(kind=ip) :: CellArea

      ! Both these concentration variables are in the 'natural' units of kg/km3
      real(kind=ip),dimension(nzmax) :: TotalConcentration ! concentration from all grain sizes as a vertical column
      real(kind=ip)                  :: MaxTotalConcentration
 
      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine ConcentrationCalculator"
      endif;enddo

      ! calculate cloud concentration,
      ! cloud height, and cloud bottom

      CloudArea                         = 0.0_ip
      CloudLoad(1:nxmax,1:nymax)        = 0.0_ip
      MaxConcentration(1:nxmax,1:nymax) = 0.0_ip
      MaxHeight(1:nxmax,1:nymax)        = 0.0_ip
      MinHeight(1:nxmax,1:nymax)        = 100.0_ip

      Mask_Cloud(:,:) = .false.

      if(n_gs_max.gt.0)then
        do i=imin,imax
          do j=jmin,jmax
            CellArea = sigma_nz_pd(i,j,1)
            TotalConcentration = 0.0_ip
            do k=1,nzmax
               ! Increment the cloud load for this column
              CloudLoad(i,j) = CloudLoad(i,j) + &
                                sum(concen_pd(i,j,k,1:n_gs_max,ts1)) * & ! in kg/km^3
                                   dz_vec_pd(k)                      / & ! convert to kg/km^2
                                1.0e3_ip                                 ! tonnes/km2
               ! Pull out column total concentration
              TotalConcentration(k) = sum(concen_pd(i,j,k,1:n_gs_max,ts1))

               ! Note max height
              if(TotalConcentration(k).ge.CLOUDCON_THRESH)then
                ! this is overwritten if the k+1 is also identified as an ash cloud
                MaxHeight(i,j)=z_cc_pd(k)+0.5_ip*dz_vec_pd(k)
              endif
            enddo

             ! Now go from the top down and find the bottom
            do k=nzmax,1,-1
             if(TotalConcentration(k).le.CLOUDCON_THRESH)then
                MinHeight(i,j)=z_cc_pd(k)-0.5_ip*dz_vec_pd(k)
              endif
            enddo

            MaxTotalConcentration = maxval(TotalConcentration(:))

            ! Generate the cloud mask base on integrated cloud load
            if(CloudLoad(i,j).ge.CLOUDLOAD_THRESH)then
              Mask_Cloud(i,j) = .true.
            else
              Mask_Cloud(i,j) = .false.
            endif

            ! Finally, modify the output variables based on the mask
            if (Mask_Cloud(i,j))then
              ! This is a cloud, accumulate area and error-check heights
              if(MinHeight(i,j).lt.0.0_ip.and. &
                 MaxHeight(i,j).gt.0.0_ip) &
                   MinHeight(i,j) = 0.0_ip
              ! Double-check that min doesn't exceed max
              MinHeight(i,j)=min(MaxHeight(i,j),MinHeight(i,j))

              CloudArea = CloudArea + sigma_nz_pd(i,j,1)
                ! Set the output variable for max concentration in mg/m3
              MaxConcentration(i,j) = MaxTotalConcentration * KG_2_MG/KM3_2_M3
            else
              ! Not an ash cloud column so write 0.0
              ! This will be set to FillValue in netcdf when writing out
              CloudLoad(i,j)        = 0.0_ip
              MaxConcentration(i,j) = 0.0_ip
              MaxHeight(i,j)        = 0.0_ip
              MinHeight(i,j)        = 0.0_ip
            endif

          enddo
        enddo
      endif

      Calculated_Cloud_Load = .true.

      end subroutine ConcentrationCalculator      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  CloudAreaCalculator
!
!  Called from: Gen_Output_Vars
!  Arguments:
!    none
!
!  This subroutine calculates the area of the cloud with a cloud load that
!  exceeds various threshholds.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CloudAreaCalculator

      use mesh,          only : &
         sigma_nz_pd

      use solution,      only : &
         imin,imax,jmin,jmax

      use Tephra,        only : &
         n_gs_max

      integer :: i,j,k

      real(kind=ip) :: CellArea

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine CloudAreaCalculator"
      endif;enddo

      ! calculate cloud concentration
      LoadVal(1)         = 0.24_ip
      LoadVal(2)         = 1.0_ip
      LoadVal(3)         = 2.0_ip
      LoadVal(4)         = 4.0_ip
      LoadVal(5)         = 6.0_ip
      CloudLoadArea(1:5) = 0.0_ip

      if(n_gs_max.gt.0)then
        do i=imin,imax
          do j=jmin,jmax
            CellArea = sigma_nz_pd(i,j,1)
            do k=1,5
              if (CloudLoad(i,j).gt.LoadVal(k)) &
                CloudLoadArea(k) = CloudLoadArea(k) + CellArea
            enddo
          enddo
        enddo
      endif

      end subroutine CloudAreaCalculator      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Gen_Output_Vars
!
!  Called from: Ash3d.F90
!  Arguments:
!    none
!
!  This subroutine calls routines to calculate deposit thickness, ash cloud
!  variables, and deposit/cloud masks.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Gen_Output_Vars

      use io_data,       only : &
         Called_Gen_Output_Vars

      use mesh,          only : &
         nxmax, nymax

      use time_data,     only : &
         time

      integer :: i,j

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Gen_Output_Vars"
      endif;enddo

      Mask_Cloud   = .false.
      Mask_Deposit = .false.

      if(.not.Calculated_AshThickness)then
        ! Need to calculate deposit thickness since this is a proxy for
        ! where the deposit is located and where FillValues should used
        call AshThicknessCalculator
      endif

      if(.not.Calculated_Cloud_Load)then
        ! Need to calculate cloud load since this is a proxy for
        ! where the cloud is located and where FillValues should used
        call ConcentrationCalculator
      endif

      call CloudAreaCalculator

      ! Mark the arrival time of any new deposit
      do i=1,nxmax
        do j=1,nymax
          if(Mask_Deposit(i,j).and.DepArrivalTime(i,j).lt.0.0_ip)then
            DepArrivalTime(i,j)=time
          endif
          if((Mask_Cloud(i,j)).and.(CloudArrivalTime(i,j).lt.0.0_ip))then
            CloudArrivalTime(i,j)=time
          endif
        enddo
      enddo

      Called_Gen_Output_Vars = .true.

      end subroutine Gen_Output_Vars

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Calc_AshVol_Aloft(vol)
!
!  Called from: TimeStepTotals
!  Arguments:
!    vol  = output variable for the total volume of ash aloft
!
!  This subroutine calculates the total mass aloft and converts it to vol in
!  km^3 DRE for logging.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Calc_AshVol_Aloft(vol)

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,kappa_pd,ts0

      use solution,      only : &
         concen_pd,mass_aloft

      use Tephra,        only : &
         MagmaDensity,n_gs_max

      real(kind=ip),intent(out) :: vol ! Total volume of ash still airborne
      integer :: isize

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Calc_AshVol_Aloft"
      endif;enddo

      vol = 0.0_ip

      ! First calculate mass of all species
      do isize=1,nsmax
        mass_aloft(isize) = sum(concen_pd(1:nxmax,1:nymax,1:nzmax,isize,ts0) *   & ! in kg/km^3
                                 kappa_pd(1:nxmax,1:nymax,1:nzmax))                ! convert to kg
      enddo

      ! Now loop over just the tephra bins (first n_gs_max bins) and
      ! calculate tephra volume
      if(n_gs_max.gt.0)then
        do isize=1,n_gs_max
            ! Increment total ash in air
          vol = vol + mass_aloft(isize)               /   & ! in kg
                      MagmaDensity                    /   & ! convert to m3
                      KM3_2_M3                              ! convert to km3

        enddo
      endif

      end subroutine Calc_AshVol_Aloft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Calc_vprofile(itime)
!
!  Called from: Ash3d.F90
!  Arguments:
!    itime  =  index of time_native (full time array)
!
!  This subroutine sums the concentrations over each tephra bin above the 
!  profile site and saves is to pr_ash in mg/m3.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Calc_vprofile(itime)

      use global_param,  only : &
         KG_2_MG,KM3_2_M3

      use io_data,       only : &
         nvprofiles,i_vprofile,j_vprofile

      use mesh,          only : &
         nzmax,ts1

      use solution,      only : &
         concen_pd

      use time_data,     only : &
         time,ntmax

      use Tephra,         only : &
         n_gs_max

      integer, intent(in) :: itime

      integer :: i,k
      real(kind=ip) :: totalash

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Calc_vprofile"
      endif;enddo

      if(itime.gt.ntmax)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: itime is greater than ntmax"
          write(errlog(io),*)"       cannot write to profile"
        endif;enddo
        return
      else
        time_native(itime) = time
      endif

      do i=1,nvprofiles
        ! Get the total ash aloft in the coloumn at this point in kg/km3
        totalash = sum(concen_pd(i_vprofile(i),j_vprofile(i),1:nzmax,1:n_gs_max,ts1))
        ! don't write if there's no ash
        if(totalash.lt.CLOUDCON_THRESH) cycle
        do k=1,nzmax
          pr_ash(k,itime,i) = sum(concen_pd(i_vprofile(i),j_vprofile(i),k,1:n_gs_max,ts1)) &
                              * KG_2_MG / KM3_2_M3 !convert from kg/km3 to mg/m3
        enddo
      enddo

      return

      end subroutine Calc_vprofile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Calc_AshVol_Deposit(vol)
!
!  Called from: TimeStepTotals
!  Arguments:
!    vol  = output variable for the total volume of ash in the deposit
!
!  This subroutine calculates the total mass of the deposit and converts it to
!  vol in km^3 DRE for logging.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Calc_AshVol_Deposit(vol)

      use mesh,          only : &
         nxmax,nymax,kappa_pd

      use solution,      only : &
         DepositGranularity

      use Tephra,        only : &
         MagmaDensity,n_gs_max

      real(kind=ip),intent(out) :: vol ! Total volume of ash in deposit
      integer :: isize

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Calc_AshVol_Deposit"
      endif;enddo

      vol = 0.0_ip

      if(n_gs_max.gt.0)then
        do isize=1,n_gs_max
          vol = vol + sum(DepositGranularity(1:nxmax,1:nymax,isize) *   & ! in kg/km^3
                           kappa_pd(1:nxmax,1:nymax,0))          /   & ! convert to kg
                      MagmaDensity                      /   & ! convert to m3
                      KM3_2_M3                                ! convert to km3
        enddo
      endif

      end subroutine Calc_AshVol_Deposit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Calc_AshVol_Outflow(vol)
!
!  Called from: TimeStepTotals
!  Arguments:
!    vol  = output variable for the total volume of ash that exited the domain
!
!  This subroutine calculates the total mass that left the domain and converts
!  it to vol in km^3 DRE for logging.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Calc_AshVol_Outflow(vol)

      use mesh,          only : &
         nxmax,nymax,nzmax,kappa_pd

      use solution,      only : &
         outflow_xz1_pd,outflow_xz2_pd,outflow_yz1_pd,outflow_yz2_pd,outflow_xy2_pd

      use Tephra,        only : &
         MagmaDensity,n_gs_max

      real(kind=ip),intent(out) :: vol ! Total volume of ash flowing out of gird

      integer :: isize

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Calc_AshVol_Outflow"
      endif;enddo

      vol = 0.0_ip

      if(n_gs_max.gt.0)then
        do isize=1,n_gs_max
          vol = vol + (                                &
                        sum(outflow_yz1_pd(        1:nymax,1:nzmax,isize)*   &
                                  kappa_pd(      0,1:nymax,1:nzmax)) +       &
                        sum(outflow_yz2_pd(        1:nymax,1:nzmax,isize)*   &
                                  kappa_pd(nxmax+1,1:nymax,1:nzmax)) +       &
                        sum(outflow_xz1_pd(1:nxmax,        1:nzmax,isize)   *   &
                                  kappa_pd(1:nxmax,        0,1:nzmax)) +          &
                        sum(outflow_xz2_pd(1:nxmax,        1:nzmax,isize)  *   &
                                  kappa_pd(1:nxmax,nymax+1,1:nzmax)) +         &
                        sum(outflow_xy2_pd(1:nxmax,1:nymax        ,isize)   *   &
                                  kappa_pd(1:nxmax,1:nymax,  nzmax+1)) )     /   & ! convert to kg
                        MagmaDensity                            /   & ! convert to m3
                        KM3_2_M3                                      ! convert to km3

        enddo
      endif

      end subroutine Calc_AshVol_Outflow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  FirstAsh
!
!  Called from: Ash3d.F90
!  Arguments:
!    none
!
!  This subroutine determines whether the ash has yet hit any airports.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine FirstAsh

      use time_data,     only : &
         time, dt

      use Airports,      only : &
         Airport_CloudHere,Airport_thickness,Airport_depRate,Airport_AshDuration,&
         Airport_thicknessLast,Airport_depRateLast,Airport_CloudHereLast,&
         Airport_CloudDuration,Airport_CloudArrived,Airport_CloudArrivalTime,&
         Airport_AshArrived,Airport_AshArrivalTime,nairports,Airport_i,Airport_j,&
           bilinear_thickness

      integer :: i

      ! Requires that AshThicknessCalculator and ConcentrationCalculator were called
      if(.not.Calculated_Cloud_Load)then
        ! Need to calculate cloud load since this is a proxy for
        ! where the cloud is located and where FillValues should used
        call ConcentrationCalculator
      endif
      if(.not.Calculated_AshThickness)then
        ! Need to calculate deposit thickness since this is a proxy for
        ! where the deposit is located and where FillValues should used
        call AshThicknessCalculator
      endif

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine FirstAsh"
      endif;enddo

        !record thickness & cloud concentration in last time step      
      Airport_CloudHereLast  = Airport_CloudHere
      Airport_ThicknessLast  = Airport_thickness
      Airport_depRateLast    = Airport_depRate
      do i=1,nairports
        Airport_CloudHere(i)  = CloudLoad(Airport_i(i),Airport_j(i))
        Airport_thickness(i)  = bilinear_thickness(i,DepositThickness)
        Airport_depRate(i)    = (Airport_thickness(i)-Airport_thicknessLast(i)) / &
                                real(dt,kind=ip)       !dep. rate, mm/hr

        !For airports where ash has arrived . . .
        !mark fall duration if ash has stopped falling
        if (Airport_AshArrived(i).eqv..true.) then
          if ((Airport_depRate(i).lt.DEPRATE_THRESH).and.&
              (Airport_depRateLast(i).ge.DEPRATE_THRESH)) then
            Airport_AshDuration(i) = time-Airport_AshArrivalTime(i)
          endif
        endif

        !mark cloud duration if cloud has passed
        if (Airport_CloudArrived(i).eqv..true.) then
          if ((Airport_CloudHere(i).le.CLOUDLOAD_THRESH).and.&
              (Airport_CloudHereLast(i).gt.CLOUDLOAD_THRESH)) then
            Airport_CloudDuration(i) = time-Airport_CloudArrivalTime(i)
          endif
        endif

        !For airports where ash has not yet arrived . . .
        !if ash load>CLOUDLOAD_THRESH T/km2 , call it "arrived"
        if ((Airport_CloudArrived(i).eqv..false.).and.&
            (Airport_CloudHere(i).gt.CLOUDLOAD_THRESH)) then
          Airport_CloudArrived(i) = .true.
          Airport_CloudArrivalTime(i) = time
        endif
        !if ash thickness>THICKNESS_THRESH mm, call it "arrived"
        if ((Airport_AshArrived(i).eqv..false.).and.&
            (Airport_thickness(i).gt.THICKNESS_THRESH)) then
          Airport_AshArrived(i) = .true.
          Airport_AshArrivalTime(i) = time
          ! Some cases with high eruptive volume might have a deposit that
          ! arrived before the 'cloud' is triggered.  Force the cloud to be
          ! marked
          if (Airport_CloudArrived(i).eqv..false.) then
            Airport_CloudArrived(i) = .true.
            Airport_CloudArrivalTime(i) = time
          endif
        endif
      enddo

      return

      end subroutine FirstAsh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module Output_Vars
!##############################################################################
