! All the modules in this file just define the essential variables grouped by
! function
!   precis_param
!   io_units
!   global_param
!   io_data
!   mesh
!   solution
!   time_data
!   time_data
!   wind_grid
! Because these are just data modules, all variables are set to public by default.
! A few of these modules have allocate/deallocate subroutines.

!##############################################################################
!
!  precis_param module
!
!  This module just defines the precision parameters used throughout the program.
!
!##############################################################################

      module precis_param

      implicit none

        ! Set everything to public by default
      public

        ! The first two precision parameters can be changed to meet your needs
      integer, parameter,public :: op         = 4 ! Output precision
      integer, parameter,public :: ip         = 8 ! Internal precision

        ! These single and double precision parameters must not be changed from
        ! 4 and 8 or the reading of the windfiles will fail.
      !integer, parameter,public :: sp         = 4 ! single precision
      !integer, parameter,public :: dp         = 8 ! double precision
      integer, parameter,public :: sp = selected_real_kind( 6,   37) ! single precision
      integer, parameter,public :: dp = selected_real_kind(15,  307) ! double precision
      integer, parameter,public :: qp = selected_real_kind(33, 4931) ! quad precision

      end module precis_param

!##############################################################################
!
!  io_units module
!
!  This module defines the input/output variables as well as verbosity
!  considerations.
!
!##############################################################################

      module io_units

      ! This module requires Fortran 2003 or later
      use iso_fortran_env, only : &
         input_unit,output_unit,error_unit

      implicit none

        ! Set everything to public by default
      public

      ! Initilize these with the defaults, but we will reset these in the subroutine
      integer,parameter :: global_log = 9    ! Log file will be unit 9
      integer :: nio  = 1                    ! number of output streams (stdout and logfile)
      integer :: io                          ! index over out-streams
      integer,dimension(2) :: outlog = (/output_unit,global_log/)
      integer,dimension(2) :: errlog = (/ error_unit,global_log/)
      integer :: errcode = 1 

      character(9) :: logfile = 'Ash3d.lst'  ! This is the default Ash3d logfile

      ! default verbosity is 3
      ! This will write everything from verbosity_log to verbosity_error to both stdout and stdlog
      ! Note that this may be reset in Set_OS_Env is ASH3DVERB is set.
      ! You could set the stdout and log to different velocity levels here is desired.
      integer,dimension(2) :: VB = (/3,3/) ! Verbosity level for stdout and logfile, respectively
      character(10)        :: vlevel       ! Text description of verbosity level

      ! These verbosity levels are for harmonizing with forestclaw
      ! Only write statements verbosity level >= VB(1) will be written
      integer,parameter :: verbosity_debug2       = 1  ! Additional debugging information only written to stdout
      integer,parameter :: verbosity_debug1       = 2  ! Debugging information only written to stdout
      integer,parameter :: verbosity_log          = 3  ! Time step information (this is the limit for writing to logfile)
      integer,parameter :: verbosity_info         = 4  ! Additional information on run set up and shutdown
      integer,parameter :: verbosity_statistics   = 5  ! Details on health of run (timing, mass conservation)
      integer,parameter :: verbosity_production   = 6  ! Major program flow info only
      integer,parameter :: verbosity_essential    = 7  ! Only start up and shutdown messages
      integer,parameter :: verbosity_error        = 8  ! No logging to stdout, only stderr (and logfile)
      integer,parameter :: verbosity_silent       = 9  ! No logging to stdout,stderr. Logfile written as normal
      integer,parameter :: verbosity_dark         = 10 ! No logging to stdout,stderr or logfile

      end module io_units

!##############################################################################
!
!  global_param module
!
!  This module sets many fixed parameters used in the run such as unit conversions
!  and some physical constants.  Also included are many variables that control the
!  activation or deactivation of routines such as horizontal/vertical advection,
!  diffusion, or how vertical velocities are calculated.
!
!##############################################################################

      module global_param

      use precis_param

      implicit none

        ! Set everything to public by default
      public

#include "Ash3d_version.h"  ! contain the git commit id

      character(len=8)  :: version           =  ' 1.0  '  ! The Ash3d version number

      real(kind=ip), parameter :: EPS_SMALL  = 1.0e-7_ip       ! Small number
      real(kind=ip), parameter :: EPS_TINY   = 1.0e-12_ip      ! Very small number
      real(kind=ip), parameter :: PI         = 3.141592653589793_ip
      real(kind=ip), parameter :: DEG2RAD    = 1.7453292519943295e-2_ip
      real(kind=ip), parameter :: DEG2KMLAT  = 111.0_ip        ! km/degree latitude 
      real(kind=ip), parameter :: DEG2KMLON  = 111.321_ip      ! km/degree longitude at equator
      real(kind=ip), parameter :: EPS_THRESH = 1.0e-10_ip      ! Threshold for Riemann solver
                                                               ! 1 kg/km3=0.001 mg/m3
      ! Unit conversions
      real(kind=ip), parameter :: KM_2_M     = 1.0e3_ip        ! km to m
      real(kind=ip), parameter :: M_2_MM     = 1.0e3_ip        ! m to mm
      real(kind=ip), parameter :: MM_2_IN    = 3.937e-2_ip     ! mm to inch
      real(kind=ip), parameter :: KM2_2_M2   = 1.0e6_ip        ! km^2 to m^2
      real(kind=ip), parameter :: KM3_2_M3   = 1.0e9_ip        ! km^3 to m^3
      real(kind=ip), parameter :: KG_2_MG    = 1.0e6_ip        ! kg to mg
      real(kind=ip), parameter :: MPS_2_KMPHR= 3.6_ip          ! m/s to km/hr
      real(kind=ip), parameter :: HR_2_S     = 3600.0_ip       ! hour to seconds


      real(kind=ip) :: GRAV       = 9.81_ip     ! Gravitational acceleration m/s^2
      real(kind=ip) :: RAD_EARTH  = 6371.229_ip ! Radius of Earth in km (used for cell
                                                ! geometry calculations: area, volume)
                                                !  Note: a particular projection might
                                                !        use a different radius

      integer,       parameter :: MAXNUM_OPTMODs   = 10   ! used just to preallocate block array
      character(len=20),dimension(MAXNUM_OPTMODs) :: OPTMOD_names
      integer                  :: nmods

      ! Some variables determined by preprocessor flags at compilation time
#ifdef LIM_NONE
      character(len=10)        :: limiter = 'No'
#endif
#ifdef LIM_LAXWEN
      character(len=10)        :: limiter = 'LaxWendrof'
#endif
#ifdef LIM_BW
      character(len=10)        :: limiter = 'BeamWarm'
#endif
#ifdef LIM_FROMM
      character(len=10)        :: limiter = 'Fromm'
#endif
#ifdef LIM_MINMOD
      character(len=10)        :: limiter = 'Minmod'
#endif
#ifdef LIM_SUPERBEE
      character(len=10)        :: limiter = 'Superbee'
#endif
#ifdef LIM_MC
      character(len=10)        :: limiter = 'MC'
#endif

#ifdef CRANKNIC
      logical, parameter       :: useCN           = .true.
#endif
#ifdef EXPLDIFF
      logical, parameter       :: useCN           = .false.
#endif

        ! These should not be changed unless needed for testing
      logical, parameter       :: useDS            = .true.  ! Dimension splitting v.s. something else
      logical, parameter       :: useVertAdvect    = .true.  ! Turns on/off vert. advection
      logical, parameter       :: useHorzAdvect    = .true.  ! Turns on/off horz. advection

        ! This will be reset based on d_coeff (.true. if d_coeff<0.0)
      logical                  :: useDiffusion     = .false. ! Reset in Read_Control_File by d_coeff
      logical                  :: useVarDiffH      = .false. ! Turned on in variable diffusion optmod
      logical                  :: useVarDiffV      = .false. ! Turned on in variable diffusion optmod

        ! These are determined when reading Reading Block 7: Grain Size Groups
      logical                  :: useCalcFallVel   = .false. ! Turned on in Read_Control_File if needed
      logical                  :: useTemperature   = .false. 
      logical                  :: useLogNormGSbins = .false.


      ! The variables below can be reset via OPTMOD=RESETPARAMS
        ! Vertical velocities come from Vertical_Velocity_Pressure (in Pa s)
        ! This can be converted to m/s by dividing by dp/dz.  We have two
        ! ways we can calculate dp/dz: using -rho g, or calculating a finite
        ! difference approximation using p and GPH variable
      logical                  :: useVz_rhoG      = .true.   ! using  -rho g
      !logical                  :: useVz_rhoG      = .false. ! using  finite-differences

        ! Only load temperature and water content data if needed
        ! This must be turned on in optional modules if needed
      logical                  :: useMoistureVars    = .false.

      real(kind=ip)            :: CFL = 0.80_ip       ! courant number
                                                      ! Note: CFL can be reset via environment
                                                      !       variables or via the input file
      real(kind=dp)            :: DT_MIN = 1.0e-5_dp  ! Minimum DT in hours
      real(kind=dp)            :: DT_MAX = 1.0e0_dp   ! Maximum DT in hours

      ! Stop conditions
      !  1 = check if amount aloft is too little
      !        aloft_vol/tot_vol.lt.(1.0_ip-StopValue)
      !  2 = check if time is past sim end
      !        time.ge.Simtime_in_hours
      !  3 = check if there is ash aloft (this might be turned off for certain sources)
      !        n_gs_aloft.eq.0
      !  4 = check on mass balance
      !        MassConsErr.gt.1.0e-3_ip
      !  5 = check on negative volumes
      !        (dep_vol.lt.-1.0_ip*EPS_SMALL).or.&
      !        (aloft_vol.lt.-1.0_ip*EPS_SMALL).or.&
      !        (outflow_vol.lt.-1.0_ip*EPS_SMALL).or.&
      !        (SourceCumulativeVol.lt.-1.0_ip*EPS_SMALL)

      logical, dimension(5) :: StopConditions  = .false.  ! Various conditions that force the run to stop
      logical, dimension(5) :: CheckConditions = .true.   ! Which conditions to check

      !  These are reset in Set_OS_Env
      integer   :: OS_TYPE                        ! 1=linux, 2=apple, 3=windows
      logical   :: IsLitEnd                       ! little-endian-ness; set in Set_OS_Env
      logical   :: IsLinux    = .true.
      logical   :: IsWindows  = .false.
      logical   :: IsMacOS    = .false.
      character (len=7)    :: OS_Flavor
      character (len=2)    :: DirPrefix  = 'c:'
      character (len=1)    :: DirDelim   = '/'
      character (len=255)  :: os_full_command_line
      character (len=32)   :: os_user
      character (len=50)   :: os_host
      character (len=255)  :: os_cwd

      end module global_param

!##############################################################################
!
!  io_data module
!
!  This module stores variables associated with input and out.  Much of the
!  information from the input control file is stored in these variables to be
!  logged in the output netcdf file (if netcdf is requested).  Also, many
!  logical variables are set in Read_Control_File that specify whether point,
!  profile, 3d output are requested, whether KML or ASCII files are requested,
!  etc.  
!
!##############################################################################

      module io_data
      
      use precis_param

      implicit none

        ! Set everything to public by default
      public

      integer            :: log_step = 1

      integer            :: iout3d          ! index for output timestep of 3d/2d data
      integer            :: ioutputFormat   ! determines the format of the output
                                            ! (1=ASCII, 2=raw binary, 3=NetCDF)
      character (len=130):: Ash3dHome       ! path to Ash3d installation
      character (len=130):: infile          ! input file name
      character (len=50) :: outfile
      logical            :: LoadConcen
      character (len=80) :: concenfile
      integer            :: init_tstep
      character (len=130):: cdf_title
      character (len=80) :: cdf_institution
      character (len=80) :: cdf_source
      character (len=80) :: cdf_history
      character (len=80) :: cdf_references
      character (len=80) :: cdf_run_class  ! Forecast, Hypothetical, Analysis
      character (len=80) :: cdf_url
      character (len=80) :: cdf_comment
      character (len=80) :: cdf_conventions
      character (len=80) :: cdf_b1l1 !character strings containing parameters for netcdf file
      character (len=80) :: cdf_b1l2
      character (len=80) :: cdf_b1l3
      character (len=80) :: cdf_b1l4
      character (len=80) :: cdf_b1l5
      character (len=80) :: cdf_b1l6
      character (len=80) :: cdf_b1l7
      character (len=80) :: cdf_b1l8
      character (len=80) :: cdf_b1l9
      character (len=80) :: cdf_b3l1
      character (len=80) :: cdf_b3l2
      character (len=80) :: cdf_b3l3
      character (len=80) :: cdf_b3l4
      character (len=80) :: cdf_b3l5
      character (len=80) :: cdf_b4l1
      character (len=80) :: cdf_b4l2
      character (len=80) :: cdf_b4l3
      character (len=80) :: cdf_b4l4
      character (len=80) :: cdf_b4l5
      character (len=80) :: cdf_b4l6
      character (len=80) :: cdf_b4l7
      character (len=80) :: cdf_b4l8
      character (len=80) :: cdf_b4l9
      character (len=80) :: cdf_b4l10      
      character (len=80) :: cdf_b4l11
      character (len=80) :: cdf_b4l12
      character (len=80) :: cdf_b4l13
      character (len=80) :: cdf_b4l14
      character (len=80) :: cdf_b4l15
      character (len=80) :: cdf_b4l16
      character (len=80) :: cdf_b4l17
      character (len=80) :: cdf_b4l18
      character (len=80) :: cdf_b6l1
      character (len=80) :: cdf_b6l2
      character (len=80) :: cdf_b6l3
      character (len=80) :: cdf_b6l4
      character (len=80) :: cdf_b6l5
!
      logical            :: WriteCloudConcentration_ASCII ! .true. if cloud top files are to be written out
      logical            :: WriteCloudConcentration_KML
      logical            :: WriteCloudHeight_ASCII        ! .true. if kml file of cloud height (km) is to be written out
      logical            :: WriteCloudHeight_KML
      logical            :: WriteCloudLoad_ASCII          ! .true. if kml file of cloud load (T/km2) is to be written out
      logical            :: WriteCloudLoad_KML
      logical            :: WriteCloudTime_ASCII          ! .true. if time of cloud arrival is to be written out
      logical            :: WriteCloudTime_KML
      logical            :: WriteReflectivity_ASCII       ! .true. if radar dbZ is to be calculated
      logical            :: WriteReflectivity_KML
      logical            :: WriteDepositFinal_ASCII       ! .true. if final deposit file is to be written out
      logical            :: WriteDepositFinal_KML
      logical            :: WriteDepositTS_ASCII          ! .true. if time series of deposit files is to be written out
      logical            :: WriteDepositTS_KML
      logical            :: WriteDepositTime_ASCII        ! .true. if time of first ash is to be written out
      logical            :: WriteDepositTime_KML
      logical            :: WriteAirportFile_ASCII        ! .true. if ash arrival times at airports is to be written out
      logical            :: WriteAirportFile_KML
      logical            :: Write_PT_Data                 ! .true. if either of the above is true (writes to netcdf)
      logical            :: Write_PR_Data                 ! .true. if writing profile data
      logical            :: ReadExtAirportFile            ! .true. if external airport file is to be read 
      logical            :: AppendExtAirportFile          ! .true. if external airports in external file are appended

      logical            :: Write3dFiles            ! .true. if 3d files are to be written
      logical            :: WriteGSD                ! .true. if grain-size distribution is to be written to airport file
      logical            :: isFinal_TS              ! .true. if we're writing out the final deposit file
      
      logical            :: Output_every_TS
      logical            :: Output_at_WriteTimes
      logical            :: Output_at_logsteps
      logical            :: Called_Gen_Output_Vars
!
      character (len=30):: VolcanoName       !name of the volcano, from the ESP input file
      integer           :: iTimeNext         !index value of next time step to write
      real(kind=ip)     :: WriteInterval     !time between file writing, used only if nWriteTimes=-1
      real(kind=ip)     :: NextWriteTime     !time to write the next file
      character (len=1) :: OutputStep_Marker !=* if data were written out since the last log_step
      integer           :: nWriteTimes       !number of deposit files to write
      real(kind=ip), dimension(:), allocatable :: WriteTimes      !times (hrs after first eruption start) to write out files

      integer,parameter :: MAXPROFILES = 100
      integer           :: nvprofiles                                        ! number of vertical profiles to write out
      integer,          dimension(:), allocatable :: i_vprofile, j_vprofile  ! i and j values of vertical profiles
      real(kind=ip),    dimension(:), allocatable :: x_vprofile, y_vprofile  ! x and y of vertical profiles
      character(len=50),dimension(:), allocatable :: Site_vprofile           ! name of profile location

      ! These specify the number of user-defined variables.  This is mostly for optional modules.
      integer            :: nvar_User2d_static_XY  = 0
      integer            :: nvar_User2d_XY         = 0
      integer            :: nvar_User3d_XYGs       = 0
      integer            :: nvar_User3d_XYZ        = 0
      integer            :: nvar_User4d_XYZGs      = 0

      contains

      !------------------------------------------------------------------------
      !subroutine Allocate_io_data
      !  ! All allocatable arrays from this module are allocated in
      !  !   Input_Data.f90:Read_Control_File
      !end subroutine Allocate_io_data
      !------------------------------------------------------------------------

      subroutine Deallocate_io_data

      if(allocated(WriteTimes))   deallocate(WriteTimes)
      if(allocated(i_vprofile))   deallocate(i_vprofile)
      if(allocated(j_vprofile))   deallocate(j_vprofile)
      if(allocated(x_vprofile))   deallocate(x_vprofile)
      if(allocated(y_vprofile))   deallocate(y_vprofile)

      end subroutine Deallocate_io_data
      !------------------------------------------------------------------------

      end module io_data

!##############################################################################
!
!  mesh module
!
!  This module contains variables that define all aspect of the computational
!  mesh.
!
!##############################################################################

      module mesh

      use precis_param
  
      implicit none

        ! Set everything to public by default
      public

      integer, parameter :: ts0 = 0
      integer, parameter :: ts1 = 1

      integer            :: ivent,jvent                    ! ij coordinates of volcano
      logical            :: IsLatLon

      ! projection parameters of the computational (Ash3d) mesh.  This might be
      ! different from the projection of the NWP files, or it might be ignored
      ! if the Ash3d mesh is lon/lat
      integer            :: A3d_iprojflag
      real(kind=dp)      :: A3d_k0_scale
      real(kind=dp)      :: A3d_phi0
      real(kind=dp)      :: A3d_lam0
      real(kind=dp)      :: A3d_lam1,A3d_phi1
      real(kind=dp)      :: A3d_lam2,A3d_phi2
      real(kind=dp)      :: A3d_Re

      logical            :: IsPeriodic   = .false.
      real(kind=ip)      :: ZPADDING     = 1.3_ip
      character(len=7)   :: VarDzType
      real(kind=ip)      :: dz_const                        ! z nodal spacing (always km)
      real(kind=ip),dimension(:)    ,allocatable :: z_vec_init

      ! Dimensional parameters in km, used if IsLatLon=.False.        
      real(kind=ip)      :: gridwidth_x, gridwidth_y  ! Dimensions (in km) of the grid
      real(kind=ip)      :: xLL,xUR,yLL,yUR           ! lower-left,upper-right points of grid
      real(kind=ip)      :: dx, dy                    ! horizontal cell sizees (km)

      integer :: nxmax      ! number of nodes in x
      integer :: nymax      ! number of nodes in y
      integer :: nzmax      ! number of nodes in z
      integer :: nsmax      ! total number of species tracked in concen
                            !  i.e. all ash bins + anything else (aggs, water, chem)

      ! Variables that are allocated based on the computational grid
      ! (nxmax,nymax)
        ! The (projected) computational grid back-projected onto lat/lon
      real(kind=ip),dimension(:,:)  ,allocatable :: xy2ll_xlon
      real(kind=ip),dimension(:,:)  ,allocatable :: xy2ll_ylat

! *****************************************************************************
!     Dimensional parameters in degrees, used if IsLatLon=.True.
      real(kind=ip)      :: latLL,lonLL,latUR,lonUR   ! lat/lon of LL corner
      real(kind=ip)      :: gridwidth_e, gridwidth_n            ! Dimensions (in km) of the grid
      real(kind=ip)      :: de, dn                    !nodal spacing east & north, degrees
      real(kind=ip)      :: de_km, dn_km              !nodal spacing, km, at volcano
#ifdef USEPOINTERS
      real(kind=ip),dimension(:)    ,pointer :: dz_vec_pd ! used for variable dz cases
      real(kind=ip),dimension(:)    ,pointer :: x_cc_pd ! x_component of cell centers
      real(kind=ip),dimension(:)    ,pointer :: y_cc_pd ! y_component of cell centers
      real(kind=ip),dimension(:)    ,pointer :: z_cc_pd ! z_component of cell centers
      real(kind=ip),dimension(:)    ,pointer :: z_lb_pd ! z_component of cell lower-boundary
      real(kind=ip),dimension(:,:,:),pointer :: kappa_pd !volume of each node in km3

      real(kind=ip),dimension(:)    ,pointer :: lat_cc_pd   => null() ! lat of i,j cell centers
      real(kind=ip),dimension(:)    ,pointer :: lon_cc_pd   => null() ! lon of i,j cell centers
      real(kind=ip),dimension(:,:,:),pointer :: sigma_nx_pd => null() ! area of x face at i-1/2,j,k
      real(kind=ip),dimension(:,:,:),pointer :: sigma_ny_pd => null() ! area of y face at i,j-1/2,k
      real(kind=ip),dimension(:,:,:),pointer :: sigma_nz_pd => null() ! area of z face at i,j,k-1/2
#else
      real(kind=ip),dimension(:)    ,allocatable :: dz_vec_pd ! used for variable dz cases
      real(kind=ip),dimension(:)    ,allocatable :: x_cc_pd ! x_component of cell centers
      real(kind=ip),dimension(:)    ,allocatable :: y_cc_pd ! y_component of cell centers
      real(kind=ip),dimension(:)    ,allocatable :: z_cc_pd ! z_component of cell centers
      real(kind=ip),dimension(:)    ,allocatable :: z_lb_pd ! z_component of cell lower-boundary
      real(kind=ip),dimension(:,:,:),allocatable :: kappa_pd !volume of each node in km3
      real(kind=ip),dimension(:)    ,allocatable :: lat_cc_pd ! lat of i,j cell centers
      real(kind=ip),dimension(:)    ,allocatable :: lon_cc_pd ! lon of i,j cell centers
      real(kind=ip),dimension(:,:,:),allocatable :: sigma_nx_pd ! area of x face at i-1/2,j,k
      real(kind=ip),dimension(:,:,:),allocatable :: sigma_ny_pd ! area of y face at i,j-1/2,k
      real(kind=ip),dimension(:,:,:),allocatable :: sigma_nz_pd ! area of z face at i,j,k-1/2
#endif

! ********************************************************************************* 

      contains

      !------------------------------------------------------------------------
      subroutine Allocate_mesh
        ! Some allocatable arrays from this module are allocated in
        !   Input_Data.f90:Read_Control_File
        !    z_vec_init
      if (IsLatLon) then
        allocate(lon_cc_pd(-1:nxmax+2));                   lon_cc_pd   = 0.0_ip
        allocate(lat_cc_pd(-1:nymax+2));                   lat_cc_pd   = 0.0_ip
      else
        allocate(x_cc_pd(-1:nxmax+2));                         x_cc_pd = 0.0_ip
        allocate(y_cc_pd(-1:nymax+2));                         y_cc_pd = 0.0_ip
      endif

      allocate( kappa_pd(-1:nxmax+2,-1:nymax+2,-1:nzmax+2));  kappa_pd = 1.0_ip
      allocate(dz_vec_pd(-1:nzmax+2));                       dz_vec_pd = 0.0_ip
      allocate(  z_cc_pd(-1:nzmax+2));                         z_cc_pd = 0.0_ip
      allocate(  z_lb_pd(-1:nzmax+2));                         z_lb_pd = 0.0_ip

      allocate(sigma_nx_pd(-1:nxmax+2,-1:nymax+2,-1:nzmax+2));  sigma_nx_pd = 1.0_ip
      allocate(sigma_ny_pd(-1:nxmax+2,-1:nymax+2,-1:nzmax+2));  sigma_ny_pd = 1.0_ip
      allocate(sigma_nz_pd(-1:nxmax+2,-1:nymax+2,-1:nzmax+2));  sigma_ny_pd = 1.0_ip

      end subroutine Allocate_mesh
      !------------------------------------------------------------------------

      subroutine Deallocate_mesh

#ifdef USEPOINTERS
      if(associated(dz_vec_pd))     deallocate(dz_vec_pd)
      if(associated(x_cc_pd))       deallocate(x_cc_pd)
      if(associated(y_cc_pd))       deallocate(y_cc_pd)
      if(associated(z_cc_pd))       deallocate(z_cc_pd)
      if(associated(z_lb_pd))       deallocate(z_lb_pd)
      if(associated(kappa_pd))      deallocate(kappa_pd)
      if(associated(lon_cc_pd))     deallocate(lon_cc_pd)
      if(associated(lat_cc_pd))     deallocate(lat_cc_pd)
      if(associated(sigma_nx_pd))   deallocate(sigma_nx_pd)
      if(associated(sigma_ny_pd))   deallocate(sigma_ny_pd)
      if(associated(sigma_nz_pd))   deallocate(sigma_nz_pd)
#else
      if(allocated(dz_vec_pd))     deallocate(dz_vec_pd)
      if(allocated(x_cc_pd))       deallocate(x_cc_pd)
      if(allocated(y_cc_pd))       deallocate(y_cc_pd)
      if(allocated(z_cc_pd))       deallocate(z_cc_pd)
      if(allocated(z_lb_pd))       deallocate(z_lb_pd)
      if(allocated(kappa_pd))      deallocate(kappa_pd)
      if(allocated(sigma_nx_pd))   deallocate(sigma_nx_pd)
      if(allocated(sigma_ny_pd))   deallocate(sigma_ny_pd)
      if(allocated(sigma_nz_pd))   deallocate(sigma_nz_pd)
      if(allocated(lon_cc_pd))     deallocate(lon_cc_pd)
      if(allocated(lat_cc_pd))     deallocate(lat_cc_pd)
#endif
      if(allocated(z_vec_init))    deallocate(z_vec_init)

      end subroutine Deallocate_mesh
      !------------------------------------------------------------------------

      end module mesh

!##############################################################################
!
!  solution module
!
!  This module stores all the variables associated with the PDE such as concentration,
!  velocities (on Ash3d grid), outflow, as well as the aspect of the solution needed for 
!  evaluating stop conditions.
!
!##############################################################################

      module solution

      use precis_param

      implicit none

        ! Set everything to public by default
      public

#ifdef USEPOINTERS
      real(kind=ip),dimension(:,:,:)    ,pointer :: vx_pd => null() ! u (E) component of wind
      real(kind=ip),dimension(:,:,:)    ,pointer :: vy_pd => null() ! v (N) component of wind
      real(kind=ip),dimension(:,:,:)    ,pointer :: vz_pd => null() ! w (up) component of wind
      real(kind=ip),dimension(:,:,:,:)  ,pointer :: vf_pd => null() ! fall velocity (x,y,z,gs) (positive upward)

      real(kind=ip),dimension(:,:,:,:,:),pointer :: concen_pd      => null() !ash concentration in x,y,z,gs_bin,time
      real(kind=ip),dimension(:,:,:)    ,pointer :: outflow_xz1_pd => null()
      real(kind=ip),dimension(:,:,:)    ,pointer :: outflow_xz2_pd => null()
      real(kind=ip),dimension(:,:,:)    ,pointer :: outflow_yz1_pd => null()
      real(kind=ip),dimension(:,:,:)    ,pointer :: outflow_yz2_pd => null()
      real(kind=ip),dimension(:,:,:)    ,pointer :: outflow_xy1_pd => null()
      real(kind=ip),dimension(:,:,:)    ,pointer :: outflow_xy2_pd => null()
      real(kind=ip),dimension(:,:,:)    ,pointer :: DepositGranularity => null() ! accumulated ash mass on ground
#else
      real(kind=ip),dimension(:,:,:)    ,allocatable :: vx_pd ! u (E) component of wind
      real(kind=ip),dimension(:,:,:)    ,allocatable :: vy_pd ! v (N) component of wind
      real(kind=ip),dimension(:,:,:)    ,allocatable :: vz_pd ! w (up) component of wind
      real(kind=ip),dimension(:,:,:,:)  ,allocatable :: vf_pd ! fall velocity (x,y,z,gs) (positive upward)

      real(kind=ip),dimension(:,:,:,:,:),allocatable :: concen_pd       ! ash concentration in x,y,z,gs_bin,time
      real(kind=ip),dimension(:,:,:)    ,allocatable :: outflow_xz1_pd  ! outflow concentration in x,z,gs
      real(kind=ip),dimension(:,:,:)    ,allocatable :: outflow_xz2_pd  ! outflow concentration in x,z,gs
      real(kind=ip),dimension(:,:,:)    ,allocatable :: outflow_yz1_pd  ! outflow concentration in y,z,gs
      real(kind=ip),dimension(:,:,:)    ,allocatable :: outflow_yz2_pd  ! outflow concentration in y,z,gs
      real(kind=ip),dimension(:,:,:)    ,allocatable :: outflow_xy1_pd  ! outflow concentration in x,y,gs
      real(kind=ip),dimension(:,:,:)    ,allocatable :: outflow_xy2_pd  ! outflow concentration in x,y,gs
      real(kind=ip),dimension(:,:,:)    ,allocatable :: DepositGranularity ! accumulated ash mass on ground 
#endif
      real(kind=ip)      :: dep_percent_accumulated
      real(kind=ip)      :: aloft_percent_remaining
      real(kind=ip)      :: StopValue    !program stops when percent_accumulated>StopValue
      real(kind=ip)      :: dep_vol,aloft_vol,outflow_vol,tot_vol
      real(kind=ip)      :: SourceCumulativeVol

      real(kind=ip), dimension(:), allocatable :: mass_aloft

      integer,       dimension(:), allocatable :: SpeciesID  ! 1 = ash gs bin
                                                             ! 2 = aggregate
                                                             ! 3 = chem species
      integer,       dimension(:), allocatable :: SpeciesSubID ! categorization within the class

      real(kind=ip), dimension(:), allocatable :: v_s         ! Settling vel 
      real(kind=ip), dimension(:), allocatable :: gsdiam      ! diameter (m)
      real(kind=ip), dimension(:), allocatable :: bin_mass    ! mass
      real(kind=ip), dimension(:), allocatable :: rho_m       ! density (kg/m3)

      logical,       dimension(:), allocatable :: IsAloft    ! T/F indicator remaining airborn concentration

        ! These are the max/min indices of the ash cloud used if FAST_SUBGRID is used
      integer :: imin,imax
      integer :: jmin,jmax
      integer :: kmin,kmax

      contains

      !------------------------------------------------------------------------
      subroutine Allocate_solution

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,ts0,ts1

      allocate(vx_pd(-1:nxmax+2,-1:nymax+2,-1:nzmax+2));         vx_pd = 0.0_ip 
      allocate(vy_pd(-1:nxmax+2,-1:nymax+2,-1:nzmax+2));         vy_pd = 0.0_ip
      allocate(vz_pd(-1:nxmax+2,-1:nymax+2,-1:nzmax+2));         vz_pd = 0.0_ip
      allocate(vf_pd(-1:nxmax+2,-1:nymax+2,-1:nzmax+2,1:nsmax)); vf_pd = 0.0_ip

      allocate(concen_pd(-1:nxmax+2,-1:nymax+2,-1:nzmax+2,1:nsmax,ts0:ts1)); concen_pd = 0.0_ip
      allocate(outflow_xz1_pd(-1:nxmax+2,-1:nzmax+2,1:nsmax)); outflow_xz1_pd = 0.0_ip
      allocate(outflow_xz2_pd(-1:nxmax+2,-1:nzmax+2,1:nsmax)); outflow_xz2_pd = 0.0_ip
      allocate(outflow_yz1_pd(-1:nymax+2,-1:nzmax+2,1:nsmax)); outflow_yz1_pd = 0.0_ip
      allocate(outflow_yz2_pd(-1:nymax+2,-1:nzmax+2,1:nsmax)); outflow_yz2_pd = 0.0_ip
      allocate(outflow_xy1_pd(-1:nxmax+2,-1:nymax+2,1:nsmax)); outflow_xy1_pd = 0.0_ip
      allocate(outflow_xy2_pd(-1:nxmax+2,-1:nymax+2,1:nsmax)); outflow_xy2_pd = 0.0_ip

      ! DepositGranularity should probably be a part of the Tephra
      ! module with trailing dimension of n_gs_max
      allocate(DepositGranularity(nxmax,nymax,nsmax)); DepositGranularity = 0.0_ip

      if (.not. allocated(mass_aloft)) then
        allocate(mass_aloft(1:nsmax)); 
        mass_aloft = 0.0_ip
        allocate(SpeciesID(1:nsmax));  SpeciesID = 1  ! Initialize everything to ash
                                                      ! If nsmax>n_gs_max, then the
                                                      ! extra bins will need to be
                                                      ! flagged in the custom source modules

        allocate(SpeciesSubID(1:nsmax));  SpeciesSubID = 0
        allocate(     v_s(1:nsmax)); v_s      = 0.0_ip
        allocate(  gsdiam(1:nsmax)); gsdiam   = 0.0_ip
        allocate(bin_mass(1:nsmax)); bin_mass = 0.0_ip
        allocate(   rho_m(1:nsmax)); rho_m    = 0.0_ip

        allocate(IsAloft(1:nsmax));   IsAloft = .true.
      endif

      end subroutine Allocate_solution
      !------------------------------------------------------------------------

      subroutine Deallocate_solution

#ifdef USEPOINTERS
      if(associated(vx_pd))              deallocate(vx_pd)
      if(associated(vy_pd))              deallocate(vy_pd)
      if(associated(vz_pd))              deallocate(vz_pd)
      if(associated(vf_pd))              deallocate(vf_pd)
      if(associated(concen_pd))          deallocate(concen_pd)
      if(associated(outflow_xz1_pd))     deallocate(outflow_xz1_pd)
      if(associated(outflow_xz2_pd))     deallocate(outflow_xz2_pd)
      if(associated(outflow_yz1_pd))     deallocate(outflow_yz1_pd)
      if(associated(outflow_yz2_pd))     deallocate(outflow_yz2_pd)
      if(associated(outflow_xy1_pd))     deallocate(outflow_xy1_pd)
      if(associated(outflow_xy2_pd))     deallocate(outflow_xy2_pd)
      if(associated(DepositGranularity)) deallocate(DepositGranularity)
#else
      if(allocated(vx_pd))              deallocate(vx_pd)
      if(allocated(vy_pd))              deallocate(vy_pd)
      if(allocated(vz_pd))              deallocate(vz_pd)
      if(allocated(vf_pd))              deallocate(vf_pd)
      if(allocated(concen_pd))          deallocate(concen_pd)
      if(allocated(outflow_xz1_pd))     deallocate(outflow_xz1_pd)
      if(allocated(outflow_xz2_pd))     deallocate(outflow_xz2_pd)
      if(allocated(outflow_yz1_pd))     deallocate(outflow_yz1_pd)
      if(allocated(outflow_yz2_pd))     deallocate(outflow_yz2_pd)
      if(allocated(outflow_xy1_pd))     deallocate(outflow_xy1_pd)
      if(allocated(outflow_xy2_pd))     deallocate(outflow_xy2_pd)
      if(allocated(DepositGranularity)) deallocate(DepositGranularity)
#endif
      if(allocated(mass_aloft))   deallocate(mass_aloft)
      if(allocated(SpeciesID))    deallocate(SpeciesID)
      if(allocated(SpeciesSubID)) deallocate(SpeciesSubID)
      if(allocated(v_s))          deallocate(v_s)
      if(allocated(gsdiam))    deallocate(gsdiam)
      if(allocated(bin_mass))  deallocate(bin_mass)
      if(allocated(rho_m))     deallocate(rho_m)
      if(allocated(IsAloft))   deallocate(IsAloft)

      end subroutine Deallocate_solution
      !------------------------------------------------------------------------

      end module solution

!##############################################################################
!
!  time_data module
!
!  This module contains all the variables associated with time
!
!##############################################################################

      module time_data

      use precis_param

      implicit none

        ! Set everything to public by default
      public

        ! If Ash3d needs to be run with wind data from a different time than the
        ! output (e.g. using a 1992 Spurr case for a present-day excercise) set
        ! this parameter for the output products
      real(kind=dp),parameter :: OutputOffset    = 0.0_dp

      integer            :: BaseYear = 1900
      logical            :: useLeap  = .true.

      ! Variables for the wall time this simulation is launched
      !  relative to BaseYear in module time_data
      character(len=13)  :: RunStartHour_ch  ! start hour of model run
      integer            :: RunStartYear
      integer            :: RunStartMonth
      integer            :: RunStartDay
      integer            :: RunStartHr
      integer            :: RunStartMinute

      ! Note: Hours should must be stored at kind=8 since forecast runs to work with HoursSince
      !       If kind=4 were used, the current year with a BaseYear=0 will cause
      !       overflows and tricky failures
      real(kind=dp)      :: Simtime_in_hours ! simulated time for ash cloud transport      
      real(kind=dp)      :: SimStartHour     ! Simulation start time, in hours since 1900
      real(kind=dp)      :: time             ! physical time simulated by this model

      real(kind=ip)      :: t0,t1,t2        ! CPU time indicators
      real(kind=ip)      :: dt_ip
      real(kind=dp)      :: dt              ! dt used for actual integration
      real(kind=dp)      :: dt_meso_last    ! dt as calculated from meso_last
      real(kind=dp)      :: dt_meso_next    ! dt as calculated from meso_next

      ! Some stings that hold time data used in output files
      !character(len=17)  :: os_time_log
      character(len=20)  :: os_time_log
      character(len=20)  :: xmlSimStartTime                     !start time of simulation in xml format
      character(len=20)  :: xmlTimeSpanStart
      character(len=20)  :: xmlTimeSpanEnd    !time periods written to kml files

      integer :: ntmax ! The maximum anticipated number of steps
      real(kind=dp),dimension(:),allocatable :: time_native
      ! No allocatable arrays to allocate or deallocate

      end module time_data

!##############################################################################
!
!  wind_grid module
!
!  This module contains the Ash3d copies of the wind data on the grid of the
!  NWP files and the computational grid copies at the NWP time intervals.
!
!##############################################################################

      module wind_grid

      use precis_param

      implicit none

        ! Set everything to public by default
      public

        ! toggle used for moving pointers between last and next
      integer            :: Meso_toggle
          ! Each *_meso_[1,2]_sp holds a regridded variable at the last and next
          ! wind file time steps.  The *_meso_[last,next]_sp are pointers that
          ! point to the correct memory locations.
          ! These exist on the computational (Ash3d) grid.
          ! For the fall velocity, we use named arrays (not pointers)
#ifdef USEPOINTERS
      real(kind=sp),dimension(:,:,:)  ,pointer          :: vx_meso_last_step_sp => null()
      real(kind=sp),dimension(:,:,:)  ,pointer          :: vx_meso_next_step_sp => null()
      real(kind=sp),dimension(:,:,:)  ,pointer          :: vy_meso_last_step_sp => null()
      real(kind=sp),dimension(:,:,:)  ,pointer          :: vy_meso_next_step_sp => null()
      real(kind=sp),dimension(:,:,:)  ,pointer          :: vz_meso_last_step_sp => null()
      real(kind=sp),dimension(:,:,:)  ,pointer          :: vz_meso_next_step_sp => null()
      real(kind=sp),dimension(:,:,:)  ,pointer          :: vx_meso_1_sp => null()
      real(kind=sp),dimension(:,:,:)  ,pointer          :: vx_meso_2_sp => null()
      real(kind=sp),dimension(:,:,:)  ,pointer          :: vy_meso_1_sp => null()
      real(kind=sp),dimension(:,:,:)  ,pointer          :: vy_meso_2_sp => null()
      real(kind=sp),dimension(:,:,:)  ,pointer          :: vz_meso_1_sp => null()
      real(kind=sp),dimension(:,:,:)  ,pointer          :: vz_meso_2_sp => null()
      real(kind=sp),dimension(:,:,:,:),pointer          :: vf_meso_last_step_sp => null()
      real(kind=sp),dimension(:,:,:,:),pointer          :: vf_meso_next_step_sp => null()
#else
      real(kind=sp),dimension(:,:,:)  ,allocatable      :: vx_meso_last_step_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable      :: vx_meso_next_step_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable      :: vy_meso_last_step_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable      :: vy_meso_next_step_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable      :: vz_meso_last_step_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable      :: vz_meso_next_step_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable      :: vx_meso_1_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable      :: vx_meso_2_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable      :: vy_meso_1_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable      :: vy_meso_2_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable      :: vz_meso_1_sp
      real(kind=sp),dimension(:,:,:)  ,allocatable      :: vz_meso_2_sp
      real(kind=sp),dimension(:,:,:,:),allocatable      :: vf_meso_last_step_sp 
      real(kind=sp),dimension(:,:,:,:),allocatable      :: vf_meso_next_step_sp
#endif

      contains

      !------------------------------------------------------------------------
      subroutine Allocate_wind_grid

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax

      allocate(vx_meso_last_step_sp(nxmax,nymax,nzmax))
      allocate(vx_meso_next_step_sp(nxmax,nymax,nzmax))
      allocate(vy_meso_last_step_sp(nxmax,nymax,nzmax))
      allocate(vy_meso_next_step_sp(nxmax,nymax,nzmax))
      allocate(vz_meso_last_step_sp(nxmax,nymax,nzmax))
      allocate(vz_meso_next_step_sp(nxmax,nymax,nzmax))

      allocate(vx_meso_1_sp(nxmax,nymax,nzmax))
      allocate(vx_meso_2_sp(nxmax,nymax,nzmax))
      allocate(vy_meso_1_sp(nxmax,nymax,nzmax))
      allocate(vy_meso_2_sp(nxmax,nymax,nzmax))
      allocate(vz_meso_1_sp(nxmax,nymax,nzmax))
      allocate(vz_meso_2_sp(nxmax,nymax,nzmax))

      ! Space for the fall velocity of each species at met steps
      allocate(vf_meso_last_step_sp(nxmax,nymax,nzmax,nsmax)); vf_meso_last_step_sp = 0.0_sp
      allocate(vf_meso_next_step_sp(nxmax,nymax,nzmax,nsmax)); vf_meso_next_step_sp = 0.0_sp

      end subroutine Allocate_wind_grid
      !------------------------------------------------------------------------

      subroutine Deallocate_wind_grid


#ifdef USEPOINTERS
      if(associated(vx_meso_last_step_sp)) deallocate(vx_meso_last_step_sp)
      if(associated(vx_meso_next_step_sp)) deallocate(vx_meso_next_step_sp)
      if(associated(vy_meso_last_step_sp)) deallocate(vy_meso_last_step_sp)
      if(associated(vy_meso_next_step_sp)) deallocate(vy_meso_next_step_sp)
      if(associated(vz_meso_last_step_sp)) deallocate(vz_meso_last_step_sp)
      if(associated(vz_meso_next_step_sp)) deallocate(vz_meso_next_step_sp)

      if(associated(vx_meso_1_sp))   deallocate(vx_meso_1_sp)
      if(associated(vx_meso_2_sp))   deallocate(vx_meso_2_sp)
      if(associated(vy_meso_1_sp))   deallocate(vy_meso_1_sp)
      if(associated(vy_meso_2_sp))   deallocate(vy_meso_2_sp)
      if(associated(vz_meso_1_sp))   deallocate(vz_meso_1_sp)
      if(associated(vz_meso_2_sp))   deallocate(vz_meso_2_sp)

      if(associated(vf_meso_last_step_sp)) deallocate(vf_meso_last_step_sp)
      if(associated(vf_meso_next_step_sp)) deallocate(vf_meso_next_step_sp)
#else
      if(allocated(vx_meso_last_step_sp)) deallocate(vx_meso_last_step_sp)
      if(allocated(vx_meso_next_step_sp)) deallocate(vx_meso_next_step_sp)
      if(allocated(vy_meso_last_step_sp)) deallocate(vy_meso_last_step_sp)
      if(allocated(vy_meso_next_step_sp)) deallocate(vy_meso_next_step_sp)
      if(allocated(vz_meso_last_step_sp)) deallocate(vz_meso_last_step_sp)
      if(allocated(vz_meso_next_step_sp)) deallocate(vz_meso_next_step_sp)

      if(allocated(vx_meso_1_sp))   deallocate(vx_meso_1_sp)
      if(allocated(vx_meso_2_sp))   deallocate(vx_meso_2_sp)
      if(allocated(vy_meso_1_sp))   deallocate(vy_meso_1_sp)
      if(allocated(vy_meso_2_sp))   deallocate(vy_meso_2_sp)
      if(allocated(vz_meso_1_sp))   deallocate(vz_meso_1_sp)
      if(allocated(vz_meso_2_sp))   deallocate(vz_meso_2_sp)

      if(allocated(vf_meso_last_step_sp)) deallocate(vf_meso_last_step_sp)
      if(allocated(vf_meso_next_step_sp)) deallocate(vf_meso_next_step_sp)
#endif

      end subroutine Deallocate_wind_grid
      !------------------------------------------------------------------------

      end module wind_grid

!##############################################################################
