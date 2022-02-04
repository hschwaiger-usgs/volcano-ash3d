      module Ash3d_Netcdf

      use precis_param

      use io_units

      use netcdf

      implicit none

      integer :: ncid
      integer :: t_dim_id     = 0 ! Time
      integer :: x_dim_id     = 0 ! X
      integer :: y_dim_id     = 0 ! Y
      integer :: z_dim_id     = 0 ! Z
      integer :: bn_dim_id    = 0 ! Full generalized species class ID 1:nsmax
      integer :: er_dim_id    = 0 ! eruption number
      integer :: wf_dim_id    = 0 ! windfile number
      integer :: sl_dim_id    = 0 ! string length
      integer :: pt_dim_id    = 0 ! point output number (airport or POI)
      integer :: pr_dim_id    = 0 ! profile output number
      integer :: tn_dim_id    = 0 ! time on native grid

      integer :: t_len
      integer :: x_len
      integer :: y_len
      integer :: z_len
      integer :: bn_len
      integer :: pt_len
      integer :: pr_len
      integer :: tn_len

      !integer :: pj_dim_id    = 0 ! projection parameter dimension
      !  Coordinate variables
      integer :: t_var_id              = 0 ! Time
      integer :: x_var_id              = 0 ! X-distance
      integer :: y_var_id              = 0 ! Y-distance
      integer :: z_var_id              = 0 ! Z-distance
      integer :: bn_var_id             = 0 ! Species class ID
      integer :: er_var_id             = 0 ! eruption index
      integer :: wf_var_id             = 0 ! wind file index
      integer :: pt_var_id             = 0 ! point (airport/POI) index
      integer :: pr_var_id             = 0 ! profile output index
      integer :: tn_var_id             = 0 ! time (native)

      integer :: proj_var_id           = 0 ! Projection
      integer :: spec_var_id           = 0 ! Species class ID
      integer :: subspec_var_id        = 0 ! Species sub-class ID

      integer :: vx_var_id             = 0 ! Vx
      integer :: vy_var_id             = 0 ! Vy
      integer :: vz_var_id             = 0 ! Vz
      integer :: vf_var_id             = 0 ! Vf
      integer :: ashcon_var_id         = 0 ! Ash concentration
      integer :: depocon_var_id        = 0 ! Deposit mass/area
      !integer :: gencon_var_id         = 0 ! General concentration for any
      !slices of concen abov
      integer :: ashconMax_var_id      = 0 ! Max Ash concentration in column
      integer :: ashheight_var_id      = 0 ! Height of top of ash cloud
      integer :: ashload_var_id        = 0 ! Vert. integrated load of ash cloud
      integer :: radrefl_var_id        = 0 ! Radar reflectivity in dbZ
      integer :: depothick_var_id      = 0 ! Total deposit thickness
      integer :: depotime_var_id       = 0 ! Deposit arrival time
      integer :: ashcloudtime_var_id   = 0 ! Cloud arrival time
      integer :: ashcloudBot_var_id    = 0 ! Height of bottom of ash cloud

      integer :: area_var_id           = 0 ! area of cell (km^2)
      integer :: gssd_var_id           = 0 ! Grain diameter
      integer :: gsmf_var_id           = 0 ! Grain mass fraction
      integer :: gsdens_var_id         = 0 ! Grain density
      integer :: er_stime_var_id       = 0 ! eruption start time
      integer :: er_duration_var_id    = 0 ! eruption duration
      integer :: er_plumeheight_var_id = 0 ! eruption plume height
      integer :: er_volume_var_id      = 0 ! eruption volume
      integer :: wf_name_var_id        = 0 ! wind file name

      ! Airport / POI variables
      integer :: pt_x_var_id             = 0 ! x coordinate of point
      integer :: pt_y_var_id             = 0 ! y coordinate of point
      integer :: pt_code_var_id          = 0 ! Airport code (or other 3-char label)
      integer :: pt_name_var_id          = 0 ! Airport name
      integer :: pt_asharrival_var_id    = 0 ! Arrival time of ashfall
      integer :: pt_ashduration_var_id   = 0 ! Duration of ashfall
      integer :: pt_cloudarrival_var_id  = 0 ! Arrival time of ash cloud
      integer :: pt_cloudduration_var_id = 0 ! Duration of ash cloud
      integer :: pt_ashthickness_var_id  = 0 ! TS of ash accumulation

      ! Vertical profile variable
      integer :: pr_ash_var_id           = 0 ! concentration profile

      integer :: temp1_2d_var_id !,temp2_2d_var_id,temp3_2d_var_id,temp4_2d_var_id
      integer :: temp1_3d_var_id !,temp2_3d_var_id,temp3_3d_var_id,temp4_3d_var_id
      integer :: temp1_4d_var_id !,temp2_4d_var_id!,temp3_4d_var_id!,temp4_4d_var_id
      integer :: ns_extra

      ! work spaces
      real(kind=op)                                 :: dumscal_out
      integer,       dimension(:)      ,allocatable :: dum1dint_out
      real(kind=op), dimension(:)      ,allocatable :: dum1d_out
      real(kind=op), dimension(:,:)    ,allocatable :: dum2d_out
      real(kind=op), dimension(:,:,:)  ,allocatable :: dum3d_out
      character                                     :: dumchar

      real(kind=op), dimension(:,:,:,:),allocatable :: ashcon
      real(kind=op), dimension(:,:,:)  ,allocatable :: depocon

      contains

!##############################################################################
!
!    create_netcdf_file
!
!##############################################################################

      subroutine create_netcdf_file

      use global_param,  only : &
         EPS_SMALL,EPS_TINY,KM2_2_M2,useCalcFallVel,VERB

      use io_data,       only : &
         nvprofiles, &
         cdf_b1l1,cdf_b1l2,cdf_b1l3,cdf_b1l4,cdf_b1l5,cdf_b1l6,cdf_b1l7,cdf_b1l8,cdf_b1l9,&
         cdf_b3l1,cdf_b3l2,cdf_b3l3,cdf_b3l4,cdf_b3l5,cdf_b4l1,cdf_b4l2,cdf_b4l3,cdf_b4l4,&
         cdf_b4l5,cdf_b4l6,cdf_b4l7,cdf_b4l8,cdf_b4l9,cdf_b4l10,cdf_b4l11,cdf_b6l1,cdf_b6l2,&
         cdf_b6l3,cdf_b6l4,cdf_conventions,&
         cdf_comment,cdf_title,cdf_institution,cdf_source,cdf_history,cdf_references,outfile,&
         nvar_User2d_static_XY,nvar_User2d_XY,nvar_User3d_XYGs,nvar_User3d_XYZ,&
         nvar_User4d_XYZGs

      use Output_Vars,   only : &
         DepositThickness_FillValue,MaxConcentration_FillValue,DepArrivalTime_FillValue,&
         CloudArrivalTime_FillValue,CloudLoad_FillValue,MaxHeight_FillValue,&
         MinHeight_FillValue,dbZCol_FillValue,DepositThickness_FillValue,&
         CloudArrivalTime_FillValue, &
         var_User2d_static_XY_name,var_User2d_static_XY_unit,var_User2d_static_XY_lname,&
         var_User2d_static_XY_MissVal,var_User2d_static_XY_FillVal,var_User2d_static_XY,&
         var_User2d_XY_name,var_User2d_XY_unit,var_User2d_XY_lname,&
         var_User2d_XY_MissVal,var_User2d_XY_FillVal,var_User2d_XY,&
         var_User3d_XYGs_name,var_User3d_XYGs_unit,var_User3d_XYGs_lname,&
         var_User3d_XYGs_MissVal,var_User3d_XYGs_FillVal,var_User3d_XYGs,&
         var_User3d_XYZ_name,var_User3d_XYZ_unit,var_User3d_XYZ_lname,&
         var_User3d_XYZ_MissVal,var_User3d_XYZ_FillVal,var_User3d_XYZ,&
         var_User4d_XYZGs_name,var_User4d_XYZGs_unit,var_User4d_XYZGs_lname,&
         var_User4d_XYZGs_MissVal,var_User4d_XYZGs_FillVal,var_User4d_XYZGs,&
         DBZ_THRESH,USE_OPTMOD_VARS,USE_RESTART_VARS,&
         USE_OUTPROD_VARS,USE_WIND_VARS,DepositThickness,DepArrivalTime,CloudArrivalTime,&
         MaxConcentration,MaxHeight,CloudLoad,dbZ,MinHeight,Mask_Cloud,Mask_Deposit,&
           dbZCalculator

      use Airports,      only : &
         nairports,Airport_Code,Airport_Name,Airport_x,Airport_y,&
         Airport_Latitude,Airport_Longitude,Airport_Thickness_TS

      use Tephra,        only : &
         n_gs_max,Tephra_gsdiam,Tephra_bin_mass,Tephra_rho_m

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,x_cc_pd,y_cc_pd,z_cc_pd,lon_cc_pd,lat_cc_pd,&
         sigma_nz_pd,dx,dy,dz_vec_pd,IsLatLon,ts1,&
         A3d_iprojflag,A3d_k0_scale,A3d_phi0,A3d_lam0,A3d_lam1,A3d_phi1,A3d_lam2,&
         A3d_phi2,A3d_radius_earth

      use solution,      only : &
          vx_pd,vy_pd,vz_pd,vf_pd,concen_pd,DepositGranularity,SpeciesID,SpeciesSubID

      use time_data,     only : &
          BaseYear,useLeap,cdf_time_log,time,SimStartHour,xmlSimStartTime,OutputOffset

      use Source,        only : &
         neruptions,e_Volume,e_Duration,e_StartTime,e_PlumeHeight

      use MetReader,     only : &
         MR_iwindfiles,MR_windfiles,MR_MetStep_Hour_since_baseyear

      implicit none

      integer :: nSTAT
!      integer :: ncid
!      integer :: t_dim_id     = 0 ! Time
!      integer :: x_dim_id     = 0 ! X
!      integer :: y_dim_id     = 0 ! Y
!      integer :: z_dim_id     = 0 ! Z
!      integer :: bn_dim_id    = 0 ! Full generalized species class ID       1:nsmax
!      integer :: er_dim_id    = 0 ! eruption number
!      integer :: wf_dim_id    = 0 ! windfile number
!      integer :: sl_dim_id    = 0 ! string length
!      integer :: pt_dim_id    = 0 ! point output number (airport or POI)
!      integer :: pr_dim_id    = 0 ! profile output number
!      integer :: tn_dim_id    = 0 ! time on native grid
!
!      !integer :: pj_dim_id    = 0 ! projection parameter dimension
!      !  Coordinate variables
!      integer :: t_var_id              = 0 ! Time
!      integer :: x_var_id              = 0 ! X-distance
!      integer :: y_var_id              = 0 ! Y-distance
!      integer :: z_var_id              = 0 ! Z-distance
!      integer :: bn_var_id             = 0 ! Species class ID
!      integer :: er_var_id             = 0 ! eruption index
!      integer :: wf_var_id             = 0 ! wind file index
!      integer :: pt_var_id             = 0 ! point (airport/POI) index
!      integer :: pr_var_id             = 0 ! profile output index
!      integer :: tn_var_id             = 0 ! time (native)
!
!      integer :: proj_var_id           = 0 ! Projection
!      integer :: spec_var_id           = 0 ! Species class ID
!      integer :: subspec_var_id        = 0 ! Species sub-class ID
!
!      integer :: vx_var_id             = 0 ! Vx
!      integer :: vy_var_id             = 0 ! Vy
!      integer :: vz_var_id             = 0 ! Vz
!      integer :: vf_var_id             = 0 ! Vf
!      integer :: ashcon_var_id         = 0 ! Ash concentration
!      integer :: depocon_var_id        = 0 ! Deposit mass/area
!      !integer :: gencon_var_id         = 0 ! General concentration for any slices of concen above the # of GS
!      integer :: ashconMax_var_id      = 0 ! Max Ash concentration in column
!      integer :: ashheight_var_id      = 0 ! Height of top of ash cloud
!      integer :: ashload_var_id        = 0 ! Vert. integrated load of ash cloud
!      integer :: radrefl_var_id        = 0 ! Radar reflectivity in dbZ
!      integer :: depothick_var_id      = 0 ! Total deposit thickness
!      integer :: depotime_var_id       = 0 ! Deposit arrival time
!      integer :: ashcloudtime_var_id   = 0 ! Cloud arrival time
!      integer :: ashcloudBot_var_id    = 0 ! Height of bottom of ash cloud
!
!      integer :: area_var_id           = 0 ! area of cell (km^2)
!      integer :: gssd_var_id           = 0 ! Grain diameter
!      integer :: gsmf_var_id           = 0 ! Grain mass fraction
!      integer :: gsdens_var_id         = 0 ! Grain density
!      integer :: er_stime_var_id       = 0 ! eruption start time
!      integer :: er_duration_var_id    = 0 ! eruption duration
!      integer :: er_plumeheight_var_id = 0 ! eruption plume height
!      integer :: er_volume_var_id      = 0 ! eruption volume
!      integer :: wf_name_var_id        = 0 ! wind file name
!
!      ! Airport / POI variables
!      integer :: pt_x_var_id             = 0 ! x coordinate of point
!      integer :: pt_y_var_id             = 0 ! y coordinate of point
!      integer :: pt_code_var_id          = 0 ! Airport code (or other 3-char label)
!      integer :: pt_name_var_id          = 0 ! Airport name
!      integer :: pt_asharrival_var_id    = 0 ! Arrival time of ashfall
!      integer :: pt_ashduration_var_id   = 0 ! Duration of ashfall
!      integer :: pt_cloudarrival_var_id  = 0 ! Arrival time of ash cloud
!      integer :: pt_cloudduration_var_id = 0 ! Duration of ash cloud
!      integer :: pt_ashthickness_var_id  = 0 ! TS of ash accumulation
!
!      ! Vertical profile variable
!      integer :: pr_ash_var_id           = 0 ! concentration profile
!
!      integer :: temp1_2d_var_id !,temp2_2d_var_id,temp3_2d_var_id,temp4_2d_var_id
!      integer :: temp1_3d_var_id !,temp2_3d_var_id,temp3_3d_var_id,temp4_3d_var_id
!      integer :: temp1_4d_var_id !,temp2_4d_var_id!,temp3_4d_var_id!,temp4_4d_var_id
!      integer :: ns_extra

!      character(len=30),dimension(5) :: temp_2d_name
!      character(len=30),dimension(5) :: temp_2d_unit
!      character(len=30),dimension(5) :: temp_3d_name
!      character(len=30),dimension(5) :: temp_3d_unit
!      character(len=30),dimension(5) :: temp_4d_name
!      character(len=30),dimension(5) :: temp_4d_unit
      
      character(len=32)              :: time_units

      character (len=20)  :: cdf_user
      character (len=50)  :: cdf_host
      character (len=255) :: cdf_cwd
      character (len=20)  :: cdf_WindStartTime

      ! Since output precision might be different from input precision,
      ! we need to allocate to correct memory space
!      real(kind=op)                                 :: dumscal_out
!      integer,       dimension(:)      ,allocatable :: dum1dint_out
!      real(kind=op), dimension(:)      ,allocatable :: dum1d_out
!      real(kind=op), dimension(:,:)    ,allocatable :: dum2d_out
!      real(kind=op), dimension(:,:,:)  ,allocatable :: dum3d_out
!      character                                     :: dumchar

!      real(kind=op), dimension(:,:,:,:),allocatable :: ashcon
!      real(kind=op), dimension(:,:,:)  ,allocatable :: depocon
      !real(kind=op), dimension(:,:,:,:),allocatable :: gencon

      character(len=3) ,dimension(11) :: dim_names
      character(len=30),dimension(11) :: dim_lnames
      character(len=30),dimension(40) :: var_lnames
      character (len=13)         :: reftimestr
      character(len=130):: lllinebuffer

      integer :: NCversion
      integer :: NCsubversion
      integer :: i,j,k,n
      integer :: ivar
      integer,dimension(5) :: chunksizes5

      INTERFACE
        character (len=13) function HS_yyyymmddhhmm_since(HoursSince,byear,useLeaps)
          real(kind=8)               ::  HoursSince
          integer                    ::  byear
          logical                    ::  useLeaps
        end function HS_yyyymmddhhmm_since
        character (len=20) function HS_xmltime(HoursSince,byear,useLeaps)
          real(kind=8)              :: HoursSince
          integer                   :: byear
          logical                   :: useLeaps
        end function HS_xmltime
      END INTERFACE

      if(VERB.gt.1)write(global_info,*)"Inside create_netcdf_file"

      allocate(dum2d_out(nxmax,nymax))
      allocate(dum3d_out(nxmax,nymax,nzmax))

      ! This space is also used for other 4-d variables
      allocate(ashcon(nxmax,nymax,nzmax,nsmax))
      allocate(depocon(nxmax,nymax,nsmax))
      if(nsmax.gt.n_gs_max)then
        ! These are some non-ash species
        ns_extra = nsmax-n_gs_max
      !  allocate(gencon(nxmax,nymax,nzmax,ns_extra))
      else
        ns_extra = 0
      endif

      ! COARDS dimension definition standard
      !  see http://ferret.wrc.noaa.gov/noaa_coop/coop_cdf_profile.html
      ! Order of dimesions is important, but the names are not
      ! standardized; i.e. "t" could be "time" or "record" or "date",
      ! "z" could be "height" or "elev" or "level", etc.
      ! Note: These dimesion names will also be used to define the
      ! coordinate variables.
      dim_names(1) = "t"
      dim_lnames(1) = "Time"
      dim_names(2) = "z"
      dim_lnames(2) = "Elevation"
      if (IsLatLon) then
        dim_names(3) = "lat"
        dim_lnames(3) = "Latitude"
        dim_names(4) = "lon"
        dim_lnames(4) = "Longitude"
      else
        dim_names(3) = "y"
        dim_lnames(3) = "Y (North)"
        dim_names(4) = "x"
        dim_lnames(4) = "X (East)"
      endif
      dim_names(5)  = "bn" ! grainsize bin for 1:n_gs_max, then generalized bin up to nsmax
      dim_lnames(5) = "Bin index; Grainsize,spec. ID"
      dim_names(6)  = "er" ! eruption index
      dim_lnames(6) = "Eruption number"
      dim_names(7)  = "wf" ! windfile index
      dim_lnames(7) = "Wind file number"
      dim_names(8)  = "sl" ! string length
      dim_lnames(8) = "string length"
      dim_names(9)  = "pt" ! point index for airport or POI points
      dim_lnames(9) = "Airport/POI point number"
      dim_names(10) = "pr" ! profile index
      dim_lnames(10)= "Profile number"
      dim_names(11) = "tn" ! time (native)
      dim_lnames(11)= "Time native"

      
      var_lnames(8) = "Wind velocity (x)"
      var_lnames(9) = "Wind velocity (y)"
      var_lnames(10) = "Wind velocity (z)"
      var_lnames(11) = "Ash/species Concentration"
      var_lnames(12) = "Mass/area of ash deposit"
      var_lnames(13) = "grain diameter"
      var_lnames(14) = "Mass fraction"
      var_lnames(15) = "Eruption start time"
      var_lnames(16) = "Eruption duration"
      var_lnames(17) = "Eruption plume height"
      var_lnames(18) = "Eruption volume"
      var_lnames(19) = "Wind file name"

      var_lnames(20) = "Grain density"
      var_lnames(21) = "Deposit arrival time"
      var_lnames(22) = "Elevation"

      var_lnames(30) = "Deposit thickness"
      var_lnames(31) = "Airborne ash arrival time"
      var_lnames(32) = "Ash Concentration (Max)"
      var_lnames(33) = "Cloud Top Height"
      var_lnames(34) = "Cloud Load"
      var_lnames(35) = "Radar Reflectivity"
      var_lnames(36) = "Cloud Bottom Height"

      var_lnames(38) = "Species class ID"
      var_lnames(39) = "Species sub-class ID"

      ! Declare header info
      call GetLog(cdf_user)
      call Hostnm(cdf_host)
      call getcwd(cdf_cwd)
      cdf_cwd = trim(cdf_cwd)

      !get date and time of start of first wind file
      cdf_WindStartTime = HS_xmltime(MR_MetStep_Hour_since_baseyear(1),BaseYear,useLeap)
      ! Create and open netcdf file
      if(VERB.gt.1)write(global_info,*)"Creating netcdf file"
      lllinebuffer = trim(nf90_inq_libvers())
      read(lllinebuffer,'(i1,a1,i1)')NCversion,dumchar,NCsubversion
      write(global_info,*)"Netcdf library version = ",NCversion
      !NCversion = 3
      if(NCversion.eq.4)then
        nSTAT = nf90_create(outfile,nf90_netcdf4,ncid,           &
                            cache_nelems = 1000, &
                            cache_size = 32000000)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"create outfile v4:")
      else
        nSTAT = nf90_create(outfile,nf90_clobber, ncid)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"create outfile classic:")
      endif

      ! Fill in header info
      if(VERB.gt.1)write(global_info,*)"Filling in header info"
      nSTAT = nf90_put_att(ncid,nf90_global,"title",cdf_title)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att title:")

      nSTAT = nf90_put_att(ncid,nf90_global,"institution",cdf_institution)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"tt institution:")

      nSTAT = nf90_put_att(ncid,nf90_global,"source",cdf_source)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att source:")

      cdf_history=cdf_time_log
      nSTAT = nf90_put_att(ncid,nf90_global,"history",cdf_history)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att history:")

      nSTAT = nf90_put_att(ncid,nf90_global,"references",cdf_references)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att references:")

      nSTAT = nf90_put_att(ncid,nf90_global,"Conventions",cdf_conventions)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Conventions:")

      nSTAT = nf90_put_att(ncid,nf90_global,"user",cdf_user)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att user:")
      nSTAT = nf90_put_att(ncid,nf90_global,"date",cdf_time_log)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att date:")
      nSTAT = nf90_put_att(ncid,nf90_global,"NWPStartTime",cdf_WindStartTime)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att NWPStartTime:")
      nSTAT = nf90_put_att(ncid,nf90_global,"host",cdf_host)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att host:")
      nSTAT = nf90_put_att(ncid,nf90_global,"CWD",cdf_cwd)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att CWD:")
      nSTAT = nf90_put_att(ncid,nf90_global,"comment",cdf_comment)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att comment:")
        ! Add lines copied from the input file
      nSTAT = nf90_put_att(ncid,nf90_global,"b1l1",cdf_b1l1)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b1l1:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b1l2",cdf_b1l2)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b1l2:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b1l3",cdf_b1l3)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b1l3:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b1l4",cdf_b1l4)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b1l4:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b1l5",cdf_b1l5)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b1l5:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b1l6",cdf_b1l6)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b1l6:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b1l7",cdf_b1l7)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b1l7:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b1l8",cdf_b1l8)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b1l8:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b1l9",cdf_b1l9)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b1l9:")

      nSTAT = nf90_put_att(ncid,nf90_global,"b3l1",cdf_b3l1)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b3l1:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b3l2",cdf_b3l2)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b3l2:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b3l3",cdf_b3l3)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b3l3:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b3l4",cdf_b3l4)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b3l4:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b3l5",cdf_b3l5)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b3l5:")

      nSTAT = nf90_put_att(ncid,nf90_global,"b4l1",cdf_b4l1)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b4l1:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b4l2",cdf_b4l2)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b4l2:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b4l3",cdf_b4l3)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b4l3:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b4l4",cdf_b4l4)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b4l4:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b4l5",cdf_b4l5)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b4l5:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b4l6",cdf_b4l6)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b4l6:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b4l7",cdf_b4l7)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b4l7:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b4l8",cdf_b4l8)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b4l8:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b4l9",cdf_b4l9)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b4l9:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b4l10",cdf_b4l10)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b4l10:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b4l11",cdf_b4l11)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b4l11:")

      nSTAT = nf90_put_att(ncid,nf90_global,"b6l1",cdf_b6l1)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b6l1:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b6l2",cdf_b6l2)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b6l2:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b6l3",cdf_b6l3)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b6l3:")
      nSTAT = nf90_put_att(ncid,nf90_global,"b6l4",cdf_b6l4)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b6l4:")

      ! Define dimensions
        ! t,z,y,x
        !  or time, elev, lon, lat
        !  or record, level, y, x
        ! and bn (particle bin)
        ! er (eruption index)
        ! wf (wind file index)
        ! sl (string length for storing windfile names)
        ! pt (point output index: e.g. airport/POI)
        ! pr (profile output index)
        ! tn (time on the native time grid)
      if(VERB.gt.1)write(global_info,*)"Defining dimensions"
      ! t
      nSTAT = nf90_def_dim(ncid,dim_names(1),nf90_unlimited,t_dim_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_dim t:")
      ! z
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_dim z:")
      nSTAT = nf90_def_dim(ncid,dim_names(2),nzmax,z_dim_id)
      ! y
      nSTAT = nf90_def_dim(ncid,dim_names(3),nymax,y_dim_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_dim y:")
      ! x
      nSTAT = nf90_def_dim(ncid,dim_names(4),nxmax,x_dim_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_dim x:")
      ! bn
      ! Note we specify the full nsmax, not just n_gs_max
      nSTAT = nf90_def_dim(ncid,dim_names(5),nsmax,bn_dim_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_dim bn:")
      ! er
      nSTAT = nf90_def_dim(ncid,dim_names(6),neruptions,er_dim_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_dim er:")
      ! wf
      nSTAT = nf90_def_dim(ncid,dim_names(7),MR_iwindfiles,wf_dim_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_dim wf:")
      ! sl
      nSTAT = nf90_def_dim(ncid,dim_names(8),130,sl_dim_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_dim sl:")
      ! pt
      if (nairports.gt.0)then
        nSTAT = nf90_def_dim(ncid,dim_names(9),nairports,pt_dim_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_dim pt:")
      endif

      if(NCversion.eq.4)then
        ! pr
        if (nvprofiles.gt.0)then
          nSTAT = nf90_def_dim(ncid,dim_names(10),nvprofiles,pr_dim_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_dim pr:")
        endif
        ! tn
        nSTAT = nf90_def_dim(ncid,dim_names(11),nf90_unlimited,tn_dim_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_dim tn:")
      endif

      ! Define coordinate variables
        ! X,Y,Z,time,bn,er,wf,sl,pt

      ! Define time-dependent variables
        ! Vx,Vy,Vz
        ! ashconc,depdepth
      ! Define other variables
        ! gs_setvel, gs_massfrac, gs_dens
         ! Time
      if(VERB.gt.1)write(global_info,*)"Defining coordinate variables"
      if(VERB.gt.1)write(global_info,*)"     Time",dim_names(1)
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,dim_names(1),&
                             nf90_double,& 
                             (/t_dim_id/),&
                             t_var_id)
      else
        nSTAT = nf90_def_var(ncid,dim_names(1),&
                             nf90_float,&
                             (/t_dim_id/), &
                             t_var_id)
      endif
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var t:")
      nSTAT = nf90_put_att(ncid,t_var_id,"long_name",dim_lnames(1))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att t long_name:")
      write(time_units,4313) xmlSimStartTime
4313  format('hours since ',a20)
      nSTAT = nf90_put_att(ncid,t_var_id,"units",time_units)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att units t:")
      nSTAT = nf90_put_att(ncid,t_var_id,"standard_name","time")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att t standard_name:")
      reftimestr = HS_yyyymmddhhmm_since(SimStartHour+OutputOffset,&
                                         BaseYear,useLeap)

      nSTAT = nf90_put_att(ncid,t_var_id,"ReferenceTime",reftimestr)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att t ReferenceTime:")

         ! Z
      if(VERB.gt.1)write(global_info,*)"     Z",dim_names(2)
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,dim_names(2),&
                             nf90_double,&
                             (/z_dim_id/),&
                             z_var_id)
      else
        nSTAT = nf90_def_var(ncid,dim_names(2),&
                             nf90_float,&
                             (/z_dim_id/), &
                             z_var_id)
      endif
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var z")
      nSTAT = nf90_put_att(ncid,z_var_id,"long_name",dim_lnames(2))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att z long_name")
      nSTAT = nf90_put_att(ncid,z_var_id,"units","km")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att z units")
      nSTAT = nf90_put_att(ncid,z_var_id,"standard_name","altitude")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att z standard_name")
      nSTAT = nf90_put_att(ncid,z_var_id,"positive","up")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att z positive")

         ! Y
      if(VERB.gt.1)write(global_info,*)"     Y",dim_names(3)
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,dim_names(3),&
                             nf90_double,&
                             (/y_dim_id/),&
                             y_var_id)
      else
        nSTAT = nf90_def_var(ncid,dim_names(3),&
                             nf90_float,&
                             (/y_dim_id/), &
                             y_var_id)
      endif
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var y")
      nSTAT = nf90_put_att(ncid,y_var_id,"long_name",dim_lnames(3))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att y long_name")
      if (IsLatLon) then
        nSTAT = nf90_put_att(ncid,y_var_id,"units","degrees_north")
        nSTAT = nf90_put_att(ncid,y_var_id,"standard_name","latitude")
      else
        nSTAT = nf90_put_att(ncid,y_var_id,"units","km")
        nSTAT = nf90_put_att(ncid,y_var_id,"standard_name","projection_y_coordinate")
      endif
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att y units standard_name")

         ! X
      if(VERB.gt.1)write(global_info,*)"     X",dim_names(4)
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,dim_names(4),&
                             nf90_double,&
                             (/x_dim_id/),&
                             x_var_id)
      else
        nSTAT = nf90_def_var(ncid,dim_names(4),&
                             nf90_float,&
                             (/x_dim_id/), &
                             x_var_id)
      endif
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var x")
      nSTAT = nf90_put_att(ncid,x_var_id,"long_name",dim_lnames(4))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att x long_name")
      if (IsLatLon) then
        nSTAT = nf90_put_att(ncid,x_var_id,"units","degrees_east")
        nSTAT = nf90_put_att(ncid,x_var_id,"standard_name","longitude")
      else
        nSTAT = nf90_put_att(ncid,x_var_id,"units","km")
        nSTAT = nf90_put_att(ncid,x_var_id,"standard_name","projection_x_coordinate")
      endif
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att x units standard_name")

         ! BN (Grain size bin ID)
      if(VERB.gt.1)write(global_info,*)"    Bin",dim_names(5)
      nSTAT = nf90_def_var(ncid,dim_names(5),&
                           nf90_int,&
                           (/bn_dim_id/), &
                           bn_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var bn")
      nSTAT = nf90_put_att(ncid,bn_var_id,"long_name",dim_lnames(5))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att bn long_name")
      nSTAT = nf90_put_att(ncid,bn_var_id,"units","index")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att bn units")
      nSTAT = nf90_put_att(ncid,bn_var_id,"Comment",&
                                          "This is just an index")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att bn Comment")

         ! ER (Eruption index)
      if(VERB.gt.1)write(*,*)"     ER",dim_names(6)
      nSTAT = nf90_def_var(ncid,dim_names(6),&
                           nf90_int,&
                           (/er_dim_id/),&
                           er_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var er")
      nSTAT = nf90_put_att(ncid,er_var_id,"long_name",dim_lnames(6))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att er long_name")
      nSTAT = nf90_put_att(ncid,er_var_id,"units","index")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att er units")

         ! WF (Wind file index)
      if(VERB.gt.1)write(*,*)"     WF",dim_names(7)
      nSTAT = nf90_def_var(ncid,dim_names(7),&
                           nf90_int,&
                           (/wf_dim_id/),&
                           wf_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var wf")
      nSTAT = nf90_put_att(ncid,wf_var_id,"long_name",dim_lnames(7))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att wf long_name")
      nSTAT = nf90_put_att(ncid,wf_var_id,"units","index")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att wf units")

        ! sl (string length for storing windfile names)
      ! We don't really need to explicitly have a variable for this one

        ! pt (point airport/POI index)
      if (nairports.gt.0)then
        if(VERB.gt.1)write(*,*)"     PT",dim_names(9)
        nSTAT = nf90_def_var(ncid,dim_names(9),&
                             nf90_int,&
                             (/pt_dim_id/),&
                             pt_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt")
        nSTAT = nf90_put_att(ncid,pt_var_id,"long_name",dim_lnames(9))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pt long_name")
        nSTAT = nf90_put_att(ncid,pt_var_id,"units","index")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pt units")
      endif

        ! pr (profile output index)
      if (nvprofiles.gt.0)then
        if(VERB.gt.1)write(*,*)"     PR",dim_names(10)
        nSTAT = nf90_def_var(ncid,dim_names(10),&
                             nf90_int,&
                             (/pr_dim_id/),&
                             pr_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pr")
        nSTAT = nf90_put_att(ncid,pr_var_id,"long_name",dim_lnames(10))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pr long_name")
        nSTAT = nf90_put_att(ncid,pr_var_id,"units","index")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pr units")

        ! tn (time native)
        ! we can only have a second unlimited dimension with NC version 4
        if(VERB.gt.1)write(*,*)"     TN",dim_names(11)
        nSTAT = nf90_def_var(ncid,dim_names(11),&
                             nf90_int,&
                             (/tn_dim_id/),&
                             tn_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var tn")
        nSTAT = nf90_put_att(ncid,tn_var_id,"long_name",dim_lnames(11))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att tn long_name")
        nSTAT = nf90_put_att(ncid,tn_var_id,"units","index")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att tn units")
      endif
      !   End of dimension variables
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! write projection as a variable
      if (IsLatLon) then
        if(VERB.gt.1)write(*,*)"     LatLon_Projection"
        nSTAT = nf90_def_var(ncid,"LatLon_Projection",&
                             nf90_int,&
                             proj_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var LatLon_Projection")
        nSTAT = nf90_put_att(ncid,proj_var_id,"grid_mapping_name","latitude_longitude")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att LatLon_Projection grid_mapping_name")
        nSTAT = nf90_put_att(ncid,proj_var_id,"semi_major_axis",A3d_radius_earth*1000.0_ip)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att LatLon_Projection semi_major_axis")
        nSTAT = nf90_put_att(ncid,proj_var_id,"inverse_flattening","0")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att LatLon_Projection inverse_flattening")
      else
        select case (A3d_iprojflag)
        case(0)
          ! Non-geographic projection, (x,y) only
          if(VERB.gt.1)write(*,*)"     Non-geographic"
          nSTAT = nf90_def_var(ncid,"Non-geographic",&
                               nf90_int,&
                               proj_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var Non-geographic")
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                               "grid_mapping_name", &
                               "Non-geographic")
        case(1)
          ! Polar stereographic
          if(VERB.gt.1)write(*,*)"     Polar_Stereographic"
          nSTAT = nf90_def_var(ncid,"Polar_Stereographic",&
                               nf90_int,&
                               proj_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var Polar_Stereographic")
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                               "grid_mapping_name", &
                               "polar_stereographic")
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,&
                              "put_att Polar_Stereographic grid_mapping_name")
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                              "longitude_of_projection_origin",A3d_lam0)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1, &
                              "put_att Polar_Stereographic longitude_of_projection_origin")
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                              "straight_vertical_longitude_from_pole",A3d_lam1)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1, &
                              "put_att Polar_Stereographic straight_vertical_longitude_from_pole")
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                              "scale_factor_at_projection_origin",A3d_k0_scale)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1, &
                              "put_att Polar_Stereographic scale_factor_at_projection_origin")
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                              "latitude_of_projection_origin",A3d_phi0)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1, &
                              "put_att Polar_Stereographic latitude_of_projection_origin")
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                              "earth_radius",A3d_radius_earth*1000.0_ip)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1, &
                              "put_att Polar_Stereographic earth_radius")
        case(2)
          ! Albers Equal Area
          if(VERB.gt.1)write(*,*)"     Albers Equal Area"
          nSTAT = nf90_def_var(ncid,"Albers_Equal_Area",&
                               nf90_int,&
                               proj_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var Albers_Equal_Area")
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                               "grid_mapping_name", &
                               "albers_conical_equal_area")
          ! standard_parallel - There may be 1 or 2 values.
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                               "standard_parallel",A3d_phi0)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,&
                             "put_att Albers_Equal_Area standard_parallel")
          ! longitude_of_central_meridian
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                               "longitude_of_central_meridian",A3d_lam0)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1, &
                             "put_att Albers_Equal_Area longitude_of_central_meridian")
          ! false_easting
          ! false_northing 

          nSTAT = nf90_put_att(ncid,proj_var_id,&
                              "earth_radius",A3d_radius_earth*1000.0_ip)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Albers_Equal_Area earth_radius")
        case(3)
          ! UTM

        case(4)
          ! Lambert conformal conic 
          if(VERB.gt.1)write(*,*)"     Lambert_Conformal"
          nSTAT = nf90_def_var(ncid,"Lambert_Conformal",&
                               nf90_int,&
                               proj_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var Lambert_Conformal")
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                               "grid_mapping_name", &
                               "lambert_conformal_conic")
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1, &
                             "put_att Lambert_Conformal grid_mapping_name")
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                               "standard_parallel",A3d_phi0)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1, &
                             "put_att Lambert_Conformal standard_parallel")
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                               "longitude_of_central_meridian",A3d_lam0)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1, &
                             "put_att Lambert_Conformal longitude_of_central_meridian")
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                              "latitude_of_projection_origin",A3d_phi1)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1, &
                             "put_att Lambert_Conformal latitude_of_projection_origin")
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                              "earth_radius",A3d_radius_earth*1000.0_ip)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1, &
                             "put_att Lambert_Conformal earth_radius")
        case(5)
          ! Mercator
        !        Mercator:grid_mapping_name = "mercator" ;
        !        Mercator:standard_parallel = 20. ;
        !        Mercator:longitude_of_projection_origin = 198.475006103516 ;
          if(VERB.gt.1)write(*,*)"     Mercator"
          nSTAT = nf90_def_var(ncid,"Mercator",&
                               nf90_int,&
                               proj_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var Mercator")
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                               "grid_mapping_name", &
                               "mercator")
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1, &
                             "put_att Mercator grid_mapping_name")
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                              "standard_parallel",A3d_phi0)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1, &
                             "put_att Mercator standard_parallel")
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                              "longitude_of_projection_origin",A3d_lam0)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1, &
                             "put_att Mercator longitude_of_projection_origin")

        end select
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !   Now a few other variables that are a function of BN
         ! Species class ID
      if(VERB.gt.1)write(global_info,*)"     SC : Species Class"
      nSTAT = nf90_def_var(ncid,"spec_class",&
                           nf90_int,&
                           (/bn_dim_id/), &
                           spec_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var spec_class")
      nSTAT = nf90_put_att(ncid,spec_var_id,"long_name",var_lnames(38))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att spec_class long_name")
      nSTAT = nf90_put_att(ncid,spec_var_id,"units","index")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att spec_class units")
      nSTAT = nf90_put_att(ncid,spec_var_id,"Comment","1=ash")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att spec_class Comment")

         ! Species sub-class ID
      if(VERB.gt.1)write(global_info,*)"     SSC : Species Sub-class"
      nSTAT = nf90_def_var(ncid,"spec_subclass",&
                           nf90_int,&
                           (/bn_dim_id/), &
                           subspec_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var spec_subclass")
      nSTAT = nf90_put_att(ncid,subspec_var_id,"long_name",var_lnames(39))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att spec_subclass long_name")
      nSTAT = nf90_put_att(ncid,subspec_var_id,"units","ID")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att spec_subclass units")
      nSTAT = nf90_put_att(ncid,subspec_var_id,"Comment",&
       "Non-ash species code")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att spec_subclass Comment")

         ! Grain-size diameter
      if(VERB.gt.1)write(global_info,*)"     GS, grain-diameter"
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,"gs_diameter",&
                             nf90_double,&
                             (/bn_dim_id/),&
                             gssd_var_id)
      else
        nSTAT = nf90_def_var(ncid,"gs_diameter",&
                             nf90_float,&
                             (/bn_dim_id/), &
                             gssd_var_id)
      endif
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var gs_diameter")
      nSTAT = nf90_put_att(ncid,gssd_var_id,"long_name",var_lnames(13))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att gs_diameter long_name")
      nSTAT = nf90_put_att(ncid,gssd_var_id,"units","mm")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att gs_diameter units")

         ! gs_massfrac (Mass fraction of grain size)
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,"gs_massfrac",&
                             nf90_double,&
                             (/bn_dim_id/),&
                             gsmf_var_id)
      else
        nSTAT = nf90_def_var(ncid,"gs_massfrac",&
                             nf90_float,&
                             (/bn_dim_id/), &
                             gsmf_var_id)
      endif
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var gs_massfrac")
      nSTAT = nf90_put_att(ncid,gsmf_var_id,"long_name",var_lnames(14))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att gs_massfrac long_name")
      nSTAT = nf90_put_att(ncid,gsmf_var_id,"units","fraction")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att gs_massfrac units")

         ! gs_dens (Density of grain)
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,"gs_dens",&
                             nf90_double,&
                             (/bn_dim_id/),&
                             gsdens_var_id)
      else
        nSTAT = nf90_def_var(ncid,"gs_dens",&
                             nf90_float,&
                             (/bn_dim_id/), &
                             gsdens_var_id)
      endif
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var gs_dens")
      nSTAT = nf90_put_att(ncid,gsdens_var_id,"long_name",var_lnames(20))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att gs_dens long_name")
      nSTAT = nf90_put_att(ncid,gsdens_var_id,"units","kg/m3")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att gs_dens units")

      !   Now a few other variables that are a function of ER
         ! er_stime (Start time of eruption)
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,"er_stime",&
                             nf90_double,&
                             (/er_dim_id/),&
                             er_stime_var_id)
      else
        nSTAT = nf90_def_var(ncid,"er_stime",&
                             nf90_float,&
                             (/er_dim_id/), &
                             er_stime_var_id)
      endif
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var er_stime")
      nSTAT = nf90_put_att(ncid,er_stime_var_id,"long_name",var_lnames(15))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att er_stime long_name")
      nSTAT = nf90_put_att(ncid,er_stime_var_id,"units", &
                           "hours since 1900")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att er_stime units")

         ! er_duration (Duration of eruption)
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,"er_duration",&
                             nf90_double,&
                             (/er_dim_id/),&
                             er_duration_var_id)
      else
        nSTAT = nf90_def_var(ncid,"er_duration",&
                             nf90_float,&
                             (/er_dim_id/), &
                             er_duration_var_id)
      endif
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var er_duration")
      nSTAT = nf90_put_att(ncid,er_duration_var_id,"long_name",var_lnames(16))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att er_duration long_name")
      nSTAT = nf90_put_att(ncid,er_duration_var_id,"units", "hours")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att er_duration units")

         ! er_plumeheight (Plume height of eruption)
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,"er_plumeheight",&
                             nf90_double,&
                             (/er_dim_id/),&
                             er_plumeheight_var_id)
      else
        nSTAT = nf90_def_var(ncid,"er_plumeheight",&
                             nf90_float,&
                             (/er_dim_id/), &
                             er_plumeheight_var_id)
      endif
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var er_plumeheight")
      nSTAT = nf90_put_att(ncid,er_plumeheight_var_id,"long_name",var_lnames(17))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att er_plumeheight long_name")
      nSTAT = nf90_put_att(ncid,er_plumeheight_var_id,"units", "km")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att er_plumeheight units")

         ! er_volume (Volume of eruption)
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,"er_volume",&
                             nf90_double,&
                             (/er_dim_id/),&
                             er_volume_var_id)
      else
        nSTAT = nf90_def_var(ncid,"er_volume",&
                             nf90_float,&
                             (/er_dim_id/), &
                             er_volume_var_id)
      endif
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var er_volume")
      nSTAT = nf90_put_att(ncid,er_volume_var_id,"long_name",var_lnames(18))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att er_volume long_name")
      nSTAT = nf90_put_att(ncid,er_volume_var_id,"units", &
                           "km3")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att er_volume units")

         ! Now define the other (non-time-dependent) variables
         ! wf_name (Name of windfile)
      nSTAT = nf90_def_var(ncid,"wf_name",&
                           nf90_char, &
                           (/sl_dim_id,wf_dim_id/),&
                           wf_name_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var wf_name")
      nSTAT = nf90_put_att(ncid,wf_name_var_id,"long_name",var_lnames(19))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att wf_name long_name")
      nSTAT = nf90_put_att(ncid,wf_name_var_id,"units", "string")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att wf_name units")

         ! Cell area
      if(VERB.gt.1)write(global_info,*)"     area"
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,"area",&
                             nf90_double, &
                             (/x_dim_id,y_dim_id/),                &
                             area_var_id)
      else
        nSTAT = nf90_def_var(ncid,"area",&
                             nf90_float,  &
                             (/x_dim_id,y_dim_id/),                &
                             area_var_id)
      endif
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var area")
      nSTAT = nf90_put_att(ncid,area_var_id,"long_name","Cell area")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att area long_name")
      nSTAT = nf90_put_att(ncid,area_var_id,"units","km2")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att area units")

      if(VERB.gt.1)write(global_info,*)"Defining variables"

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(USE_WIND_VARS)then
         ! Now define the time-dependent variables
         ! Vx
        if(VERB.gt.1)write(global_info,*)"     vx",var_lnames(8)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"vx",&
                               nf90_double,  &
                               (/x_dim_id,y_dim_id,z_dim_id,t_dim_id/), &
                               vx_var_id)
        else
          nSTAT = nf90_def_var(ncid,"vx",&
                               nf90_float,   &
                               (/x_dim_id,y_dim_id,z_dim_id,t_dim_id/), &
                               vx_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var vx")
        nSTAT = nf90_put_att(ncid,vx_var_id,"long_name",var_lnames(8))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att vx long_name")
        nSTAT = nf90_put_att(ncid,vx_var_id,"units","km/hr")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att vx units")
           ! Vy
        if(VERB.gt.1)write(global_info,*)"     vy",var_lnames(9)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"vy",&
                               nf90_double,  &
                               (/x_dim_id,y_dim_id,z_dim_id,t_dim_id/), &
                               vy_var_id)
        else
          nSTAT = nf90_def_var(ncid,"vy",&
                               nf90_float,   &
                               (/x_dim_id,y_dim_id,z_dim_id,t_dim_id/), &
                               vy_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var vy")
        nSTAT = nf90_put_att(ncid,vy_var_id,"long_name",var_lnames(9))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att vy long_name")
        nSTAT = nf90_put_att(ncid,vy_var_id,"units","km/hr")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att vy units")
           ! Vz
        if(VERB.gt.1)write(global_info,*)"     vz",var_lnames(10)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"vz",&
                               nf90_double,  &
                               (/x_dim_id,y_dim_id,z_dim_id,t_dim_id/), &
                               vz_var_id)
        else
          nSTAT = nf90_def_var(ncid,"vz",&
                               nf90_float,   &
                               (/x_dim_id,y_dim_id,z_dim_id,t_dim_id/), &
                               vz_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var vz")
        nSTAT = nf90_put_att(ncid,vz_var_id,"long_name","Wind velocity (z)")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att vz long_name")
        nSTAT = nf90_put_att(ncid,vz_var_id,"units","km/hr")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att vz units")

           ! Vf
        if(VERB.gt.1)write(global_info,*)"     vf","Fall Velocity"
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"vf",&
                               nf90_double,  &
                               (/x_dim_id,y_dim_id,z_dim_id,bn_dim_id,t_dim_id/), &
                               vf_var_id)
        else
          nSTAT = nf90_def_var(ncid,"vf",&
                               nf90_float,   &
                               (/x_dim_id,y_dim_id,z_dim_id,bn_dim_id,t_dim_id/), &
                               vf_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var vf")
        nSTAT = nf90_put_att(ncid,vf_var_id,"long_name","Fall Velocity")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att vf long_name")
        nSTAT = nf90_put_att(ncid,vf_var_id,"units","km/hr")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att vf units")
      endif  ! USE_WIND_VARS
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(USE_RESTART_VARS)then
         ! Full Ash Concentration
        if(VERB.gt.1)write(global_info,*)"     ashcon",var_lnames(11)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"ashcon",&
                               nf90_double, &
                               (/x_dim_id,y_dim_id,z_dim_id,bn_dim_id,t_dim_id/),    &
                               ashcon_var_id)
        else
          nSTAT = nf90_def_var(ncid,"ashcon",&
                               nf90_float,  &
                               (/x_dim_id,y_dim_id,z_dim_id,bn_dim_id,t_dim_id/),    &
                                ashcon_var_id)
        endif
        if(NCversion.eq.4)then
          chunksizes5 = (/nxmax, nymax, nzmax, nsmax, 1/)
          nSTAT = nf90_def_var_chunking(ncid, ashcon_var_id, &
                                        0, &
                                        chunksizes5)
          nSTAT = nf90_def_var_deflate(ncid, ashcon_var_id,          &
                                  shuffle = 1,                &
                                  deflate = 5,                &
                                  deflate_level = 5  )
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var ashcon")
        nSTAT = nf90_put_att(ncid,ashcon_var_id,"long_name",var_lnames(11))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att ashcon long_name")
        nSTAT = nf90_put_att(ncid,ashcon_var_id,"units","kg/km^3")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att ashcon units")
        nSTAT = nf90_put_att(ncid,ashcon_var_id,&
                 "missing_value", MaxConcentration_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att ashcon missing_value")
        nSTAT = nf90_put_att(ncid,ashcon_var_id,"_FillValue",MaxConcentration_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att ashcon _FillValue")
      endif ! USE_RESTART_VARS
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Derived variables
      if(USE_OUTPROD_VARS)then
        ! DEPOSIT (as function of grainsmax)
           ! Concentration of deposit by grain size(z=0 plane of tot_ashcon)
        if(VERB.gt.1)write(global_info,*)"     depocon",var_lnames(12)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"depocon",&
                               nf90_double, &
                               (/x_dim_id,y_dim_id,bn_dim_id,t_dim_id/),                &
                               depocon_var_id)
        else
          nSTAT = nf90_def_var(ncid,"depocon",&
                               nf90_float,  &
                               (/x_dim_id,y_dim_id,bn_dim_id,t_dim_id/),                &
                               depocon_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var depocon")
        nSTAT = nf90_put_att(ncid,depocon_var_id,"long_name",var_lnames(12))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att depocon long_name")
        nSTAT = nf90_put_att(ncid,depocon_var_id,"units","kg/m2")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att depocon units")
        nSTAT = nf90_put_att(ncid,depocon_var_id,&
                 "missing_value", -9999.0)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att depocon missing_value")
        nSTAT = nf90_put_att(ncid,depocon_var_id,"_FillValue",-9999.0_op)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att depocon _FillValue")
  
        ! DEPOSIT thickness
        if(VERB.gt.1)write(global_info,*)"     depothick",var_lnames(30)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"depothick",&
                               nf90_double, &
                               (/x_dim_id,y_dim_id,t_dim_id/),                &
                               depothick_var_id)
        else
          nSTAT = nf90_def_var(ncid,"depothick",&
                               nf90_float,  &
                               (/x_dim_id,y_dim_id,t_dim_id/),                &
                               depothick_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var depothick")
        nSTAT = nf90_put_att(ncid,depothick_var_id,"long_name",var_lnames(30))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att depothick long_name")
        nSTAT = nf90_put_att(ncid,depothick_var_id,"units","mm")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att depothick units")
        nSTAT = nf90_put_att(ncid,depothick_var_id,&
                 "missing_value", DepositThickness_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att depothick missing_value")
        nSTAT = nf90_put_att(ncid,depothick_var_id,"_FillValue",DepositThickness_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att depothick _FillValue")
 
           ! Arrival time of deposit
        if(VERB.gt.1)write(global_info,*)"     depotime",var_lnames(21)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"depotime",&
                               nf90_double, &
                               (/x_dim_id,y_dim_id/),                &
                               depotime_var_id)
        else
          nSTAT = nf90_def_var(ncid,"depotime",&
                               nf90_float,  &
                               (/x_dim_id,y_dim_id/),                &
                               depotime_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var depotime")
        nSTAT = nf90_put_att(ncid,depotime_var_id,"long_name",var_lnames(21))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var depotime long_name")
        nSTAT = nf90_put_att(ncid,depotime_var_id,"units","hr")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var depotime units")
        nSTAT = nf90_put_att(ncid,depotime_var_id,&
                   "missing_value", DepArrivalTime_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var depotime missing_value")
        nSTAT = nf90_put_att(ncid,depotime_var_id,"_FillValue",DepArrivalTime_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var depotime _FillValue")
 
           ! Arrival time of airborne ash cloud
        if(VERB.gt.1)write(global_info,*)"     ash_arrival_time",var_lnames(31)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"ash_arrival_time",&
                               nf90_double, &
                               (/x_dim_id,y_dim_id/),                &
                               ashcloudtime_var_id)
        else
          nSTAT = nf90_def_var(ncid,"ash_arrival_time",&
                               nf90_float,  &
                               (/x_dim_id,y_dim_id/),                &
                               ashcloudtime_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var ash_arrival_time")
        nSTAT = nf90_put_att(ncid,ashcloudtime_var_id,"long_name",var_lnames(31))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att ash_arrival_time long_name")
        nSTAT = nf90_put_att(ncid,ashcloudtime_var_id,"units","hr")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att ash_arrival_time units")
        nSTAT = nf90_put_att(ncid,ashcloudtime_var_id,&
                 "missing_value", CloudArrivalTime_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att ash_arrival_time missing_value")
        nSTAT = nf90_put_att(ncid,ashcloudtime_var_id,"_FillValue",CloudArrivalTime_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att ash_arrival_time _FillValue")
  
           ! 2d ash concentration (MaxConcentration)
        if(VERB.gt.1)write(global_info,*)"     ashcon_max",var_lnames(32)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"ashcon_max",&
                               nf90_double, &
                               (/x_dim_id,y_dim_id,t_dim_id/),                &
                               ashconMax_var_id)
        else
          nSTAT = nf90_def_var(ncid,"ashcon_max",&
                               nf90_float,  &
                               (/x_dim_id,y_dim_id,t_dim_id/),                &
                               ashconMax_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var ashcon_max")
        nSTAT = nf90_put_att(ncid,ashconMax_var_id,"long_name",var_lnames(32))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att ashcon_max long_name")
        nSTAT = nf90_put_att(ncid,ashconMax_var_id,"units","mg/m3")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att ashcon_max units")
        nSTAT = nf90_put_att(ncid,ashconMax_var_id,&
                 "missing_value", MaxConcentration_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att ashcon_max missing_value")
        nSTAT = nf90_put_att(ncid,ashconMax_var_id,"_FillValue",MaxConcentration_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att ashcon_max _FillValue")
  
           ! 2d ash cloud height (MaxHeight)
        if(VERB.gt.1)write(global_info,*)"     cloud_height",var_lnames(33)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"cloud_height",&
                               nf90_double, &
                               (/x_dim_id,y_dim_id,t_dim_id/),                &
                               ashheight_var_id)
        else
          nSTAT = nf90_def_var(ncid,"cloud_height",&
                               nf90_float,  &
                               (/x_dim_id,y_dim_id,t_dim_id/),                &
                               ashheight_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var cloud_height")
        nSTAT = nf90_put_att(ncid,ashheight_var_id,"long_name",var_lnames(33))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att cloud_height long_name")
        nSTAT = nf90_put_att(ncid,ashheight_var_id,"units","km") ! Note the canonical_units are m
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att cloud_height units")
        nSTAT = nf90_put_att(ncid,ashheight_var_id,"standard_name",&
                             "geopotential_height_at_volcanic_ash_cloud_top")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att cloud_height standard_name")
        nSTAT = nf90_put_att(ncid,ashheight_var_id,&
                 "missing_value", MaxHeight_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att cloud_height missing_value")
        nSTAT = nf90_put_att(ncid,ashheight_var_id,"_FillValue",MaxHeight_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att cloud_height _FillValue")
  
           ! 2d ash cloud load (CloudLoad)
        if(VERB.gt.1)write(global_info,*)"     cloud_load",var_lnames(34)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"cloud_load",&
                               nf90_double, &
                               (/x_dim_id,y_dim_id,t_dim_id/),                &
                               ashload_var_id)
        else
          nSTAT = nf90_def_var(ncid,"cloud_load",&
                               nf90_float,  &
                               (/x_dim_id,y_dim_id,t_dim_id/),                &
                               ashload_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var cloud_load")
        nSTAT = nf90_put_att(ncid,ashload_var_id,"long_name",var_lnames(34))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att cloud_load long_name")
        nSTAT = nf90_put_att(ncid,ashload_var_id,"units","T/km2") ! Note the canonical_units are kg m-2
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att cloud_load units")
        nSTAT = nf90_put_att(ncid,ashload_var_id,"standard_name",&
                             "atmosphere_mass_content_of_volcanic_ash")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att cloud_load standard_name")
        nSTAT = nf90_put_att(ncid,ashload_var_id,&
                 "missing_value", CloudLoad_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att cloud_load missing_value")
        nSTAT = nf90_put_att(ncid,ashload_var_id,"_FillValue",CloudLoad_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att cloud_load _FillValue")
  
           ! 3d radar reflectivity
        if(VERB.gt.1)write(global_info,*)"     radar_reflectivity",var_lnames(35)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"radar_reflectivity",&
                               nf90_double, &
                               (/x_dim_id,y_dim_id,z_dim_id,t_dim_id/),                &
                               radrefl_var_id)
        else
          nSTAT = nf90_def_var(ncid,"radar_reflectivity",&
                               nf90_float,  &
                               (/x_dim_id,y_dim_id,z_dim_id,t_dim_id/),                &
                               radrefl_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var radar_reflectivity")
        nSTAT = nf90_put_att(ncid,radrefl_var_id,"long_name",var_lnames(35))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att radar_reflectivity long_name")
        nSTAT = nf90_put_att(ncid,radrefl_var_id,"units","dbZ")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att radar_reflectivity units")
        nSTAT = nf90_put_att(ncid,radrefl_var_id,&
                 "missing_value", dbZCol_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att radar_reflectivity missing_value")
        nSTAT = nf90_put_att(ncid,radrefl_var_id,"_FillValue",dbZCol_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att radar_reflectivity _FillValue")
  
           ! 2d ash cloud bottom 
        if(VERB.gt.1)write(global_info,*)"     cloud_bottom",var_lnames(36)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"cloud_bottom",&
                               nf90_double, &
                               (/x_dim_id,y_dim_id,t_dim_id/),                &
                               ashcloudBot_var_id)
        else
          nSTAT = nf90_def_var(ncid,"cloud_bottom",&
                               nf90_float,  &
                               (/x_dim_id,y_dim_id,t_dim_id/),                &
                               ashcloudBot_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var cloud_bottom")
        nSTAT = nf90_put_att(ncid,ashcloudBot_var_id,"long_name",var_lnames(36))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att cloud_bottom long_name")
        nSTAT = nf90_put_att(ncid,ashcloudBot_var_id,"units","km")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att cloud_bottom units")
        nSTAT = nf90_put_att(ncid,ashcloudBot_var_id,&
                 "missing_value", MinHeight_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att cloud_bottom missing_value")
        nSTAT = nf90_put_att(ncid,ashcloudBot_var_id,"_FillValue",MinHeight_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att cloud_bottom _FillValue")

      endif ! USE_OUTPROD_VARS
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Point output (Airport/POI)
      if (nairports.gt.0)then
        ! x coordinate of point
        if(VERB.gt.1)write(global_info,*)"     pt_x"
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"pt_x",&
                               nf90_double,   &
                               (/pt_dim_id/), &
                               pt_x_var_id)
        else
          nSTAT = nf90_def_var(ncid,"pt_x",&
                               nf90_float,  &
                               (/pt_dim_id/),                &
                               pt_x_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_x")
        nSTAT = nf90_put_att(ncid,pt_x_var_id,"long_name","x coordinate of point")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pt_x long_name")
        if(IsLatLon)then
           nSTAT = nf90_put_att(ncid,pt_x_var_id,"units","degrees_east")
        else
          nSTAT = nf90_put_att(ncid,pt_x_var_id,"units","km")
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pt_x units")

        ! y coordinate of point
        if(VERB.gt.1)write(global_info,*)"     pt_y"
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"pt_y",&
                               nf90_double,   &
                               (/pt_dim_id/), &
                               pt_y_var_id)
        else
          nSTAT = nf90_def_var(ncid,"pt_y",&
                               nf90_float,  &
                               (/pt_dim_id/),                &
                               pt_y_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_y")
        nSTAT = nf90_put_att(ncid,pt_y_var_id,"long_name","y coordinate of point")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pt_y long_name")
        if(IsLatLon)then
           nSTAT = nf90_put_att(ncid,pt_y_var_id,"units","degrees_north")
        else
          nSTAT = nf90_put_att(ncid,pt_y_var_id,"units","km")
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pt_y units")

        ! Point/Airport code (3-char label)
        if(VERB.gt.1)write(global_info,*)"     pt_code"
        nSTAT = nf90_def_var(ncid,"pt_code",&
                             nf90_char,   &
                             (/sl_dim_id,pt_dim_id/),&
                             pt_code_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_code")
        nSTAT = nf90_put_att(ncid,pt_code_var_id,"long_name","3 character point label")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pt_code long_name")
        nSTAT = nf90_put_att(ncid,pt_code_var_id,"units","text")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pt_code units")

        ! Point/Airport name
        if(VERB.gt.1)write(global_info,*)"     pt_name"
        nSTAT = nf90_def_var(ncid,"pt_name",&
                             nf90_char,   &
                             (/sl_dim_id,pt_dim_id/),&
                             pt_name_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_name")
        nSTAT = nf90_put_att(ncid,pt_name_var_id,"long_name","Point name")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pt_name long_name")
        nSTAT = nf90_put_att(ncid,pt_name_var_id,"units","text")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pt_name units")

        ! Point/Airport ashfall arrival time
        if(VERB.gt.1)write(global_info,*)"     pt_depotime"
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"pt_depotime",&
                               nf90_double, &
                               (/pt_dim_id/),                &
                               pt_asharrival_var_id)
        else
          nSTAT = nf90_def_var(ncid,"pt_depotime",&
                               nf90_float,  &
                               (/pt_dim_id/),                &
                               pt_asharrival_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_depotime")
        nSTAT = nf90_put_att(ncid,pt_asharrival_var_id,"long_name","Ashfall arrival time")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_asharrival_var_id long_name")
        nSTAT = nf90_put_att(ncid,pt_asharrival_var_id,"units","hr")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_asharrival_var_id units")
        nSTAT = nf90_put_att(ncid,pt_asharrival_var_id,&
                   "missing_value", DepArrivalTime_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_asharrival_var_id missing_value")
        nSTAT = nf90_put_att(ncid,pt_asharrival_var_id,"_FillValue",DepArrivalTime_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_asharrival_var_id _FillValue")

        ! Point/Airport ashfall duration
        if(VERB.gt.1)write(global_info,*)"     pt_depodur"
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"pt_depodur",&
                               nf90_double, &
                               (/pt_dim_id/),                &
                               pt_ashduration_var_id)
        else
          nSTAT = nf90_def_var(ncid,"pt_depodur",&
                               nf90_float,  &
                               (/pt_dim_id/),                &
                               pt_ashduration_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_depotime")
        nSTAT = nf90_put_att(ncid,pt_ashduration_var_id,"long_name","Ashfall duration")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_ashduration_var_id long_name")
        nSTAT = nf90_put_att(ncid,pt_ashduration_var_id,"units","hr")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_ashduration_var_id units")
        nSTAT = nf90_put_att(ncid,pt_ashduration_var_id,&
                   "missing_value", DepArrivalTime_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_ashduration_var_id missing_value")
        nSTAT = nf90_put_att(ncid,pt_ashduration_var_id,"_FillValue",DepArrivalTime_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_ashduration_var_id _FillValue")

        ! Point/Airport ash cloud arrival time
        if(VERB.gt.1)write(global_info,*)"     pt_cloud_arrival"
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"pt_cloud_arrival",&
                               nf90_double, &
                               (/pt_dim_id/),                &
                               pt_cloudarrival_var_id)
        else
          nSTAT = nf90_def_var(ncid,"pt_cloud_arrival",&
                               nf90_float,  &
                               (/pt_dim_id/),                &
                               pt_cloudarrival_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_cloud_arrival")
        nSTAT = nf90_put_att(ncid,pt_cloudarrival_var_id,"long_name","Ash cloud arrival time")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_cloudarrival_var_id long_name")
        nSTAT = nf90_put_att(ncid,pt_cloudarrival_var_id,"units","hr")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_cloudarrival_var_id units")
        nSTAT = nf90_put_att(ncid,pt_cloudarrival_var_id,&
                   "missing_value", DepArrivalTime_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_cloudarrival_var_id missing_value")
        nSTAT = nf90_put_att(ncid,pt_cloudarrival_var_id,"_FillValue",DepArrivalTime_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_cloudarrival_var_id _FillValue")

        ! Point/Airport ash cloud duration
        if(VERB.gt.1)write(global_info,*)"     pt_cloud_dur"
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"pt_cloud_dur",&
                               nf90_double, &
                               (/pt_dim_id/),                &
                               pt_cloudduration_var_id)
        else
          nSTAT = nf90_def_var(ncid,"pt_cloud_dur",&
                               nf90_float,  &
                               (/pt_dim_id/),                &
                               pt_cloudduration_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_cloud_dur")
        nSTAT = nf90_put_att(ncid,pt_cloudduration_var_id,"long_name","Ash cloud duration")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_cloudduration_var_id long_name")
        nSTAT = nf90_put_att(ncid,pt_cloudduration_var_id,"units","hr")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_cloudduration_var_id units")
        nSTAT = nf90_put_att(ncid,pt_cloudduration_var_id,&
                   "missing_value", DepArrivalTime_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_cloudduration_var_id missing_value")
        nSTAT = nf90_put_att(ncid,pt_cloudduration_var_id,"_FillValue",DepArrivalTime_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_cloudduration_var_id _FillValue")

        ! Point/Airport ashfall accumulation
        if(VERB.gt.1)write(global_info,*)"     pt_depothick"
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"pt_depothick",&
                               nf90_double, &
                               (/pt_dim_id,t_dim_id/),                &
                               pt_ashthickness_var_id)
        else
          nSTAT = nf90_def_var(ncid,"pt_depothick",&
                               nf90_float,  &
                               (/pt_dim_id,t_dim_id/),                &
                               pt_ashthickness_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_depothick")
        nSTAT = nf90_put_att(ncid,pt_ashthickness_var_id,"long_name","Ash fall thickness")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_ashthickness_var_id long_name")
        nSTAT = nf90_put_att(ncid,pt_ashthickness_var_id,"units","mm")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_ashthickness_var_id units")
        nSTAT = nf90_put_att(ncid,pt_ashthickness_var_id,&
                   "missing_value", DepArrivalTime_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_ashthickness_var_id missing_value")
        nSTAT = nf90_put_att(ncid,pt_ashthickness_var_id,"_FillValue",DepArrivalTime_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_ashthickness_var_id _FillValue")

      endif

      if(NCversion.eq.4.and.nvprofiles.gt.0)then
        ! PR
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"pr_ash",&
                               nf90_double, &
                               (/z_dim_id,tn_dim_id,pr_dim_id/),                &
                               pr_ash_var_id)
        else
          nSTAT = nf90_def_var(ncid,"pr_ash",&
                               nf90_float,  &
                               (/z_dim_id,tn_dim_id,pr_dim_id/),                &
                               pr_ash_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pr_ash")
        nSTAT = nf90_put_att(ncid,pr_ash_var_id,"long_name",var_lnames(32))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pr_ash long_name")
        nSTAT = nf90_put_att(ncid,pr_ash_var_id,"units","mg/m3")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pr_ash units")
        nSTAT = nf90_put_att(ncid,pr_ash_var_id,&
                 "missing_value", MaxConcentration_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pr_ash missing_value")
        nSTAT = nf90_put_att(ncid,pr_ash_var_id,"_FillValue",MaxConcentration_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pr_ash _FillValue")
      endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(USE_OPTMOD_VARS)then
        if(VERB.gt.1)write(global_info,*)"     USE_OPTMOD_VARS"
        ! Define User-specified 2-d static variables
        if(nvar_User2d_static_XY.gt.0)then
          do ivar=1,nvar_User2d_static_XY
            if(op.eq.8)then
              nSTAT = nf90_def_var(ncid,var_User2d_static_XY_name(ivar),&
                                   nf90_double,  &
                                   (/x_dim_id,y_dim_id/), &
                                   temp1_2d_var_id)
            else
              nSTAT = nf90_def_var(ncid,var_User2d_static_XY_name(ivar),&
                                   nf90_float,  &
                                   (/x_dim_id,y_dim_id/), &
                                   temp1_2d_var_id)
            endif
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var XYs")
            nSTAT = nf90_put_att(ncid,temp1_2d_var_id,"long_name",&
                                 var_User2d_static_XY_lname(ivar))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att XYs long_name")
            nSTAT = nf90_put_att(ncid,temp1_2d_var_id,"units",&
                                 var_User2d_static_XY_unit(ivar))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att XYs units")
            nSTAT = nf90_put_att(ncid,temp1_2d_var_id,"missing_value",&
                                 var_User2d_static_XY_MissVal(ivar))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att XYs missing_value")
            nSTAT = nf90_put_att(ncid,temp1_2d_var_id,"_FillValue",&
                                 var_User2d_static_XY_FillVal(ivar))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att XYs _FillValue")
          enddo
        endif

        ! Define User-specified 2-d transient variables
        if(nvar_User2d_XY.gt.0)then
          do ivar=1,nvar_User2d_XY
            if(op.eq.8)then
              nSTAT = nf90_def_var(ncid,var_User2d_XY_name(ivar),&
                                   nf90_double,  &
                                   (/x_dim_id,y_dim_id,t_dim_id/), &
                                   temp1_2d_var_id)
            else
              nSTAT = nf90_def_var(ncid,var_User2d_XY_name(ivar),&
                                   nf90_float,  &
                                   (/x_dim_id,y_dim_id,t_dim_id/), &
                                   temp1_2d_var_id)
            endif
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var XY")
            nSTAT = nf90_put_att(ncid,temp1_2d_var_id,"long_name",&
                                 var_User2d_XY_lname(ivar))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att XY long_name")
            nSTAT = nf90_put_att(ncid,temp1_2d_var_id,'units',&
                                 var_User2d_XY_unit(ivar))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att XY units")
            nSTAT = nf90_put_att(ncid,temp1_2d_var_id,"missing_value",&
                                 var_User2d_XY_MissVal(ivar))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att XY missing_value")
            nSTAT = nf90_put_att(ncid,temp1_2d_var_id,"_FillValue",&
                                 var_User2d_XY_FillVal(ivar))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att XY _FillValue")
          enddo
        endif

        !! Define User-specified 3-d transient variables in x,y,gs
        if(nvar_User3d_XYGs.gt.0)then
          do ivar=1,nvar_User3d_XYGs
            if(op.eq.8)then
              nSTAT = nf90_def_var(ncid,var_User3d_XYGs_name(ivar),&
                                   nf90_double,  &
                                   (/x_dim_id,y_dim_id,bn_dim_id,t_dim_id/), &
                                   temp1_3d_var_id)
            else
              nSTAT = nf90_def_var(ncid,var_User3d_XYGs_name(ivar),&
                                   nf90_float,  &
                                   (/x_dim_id,y_dim_id,bn_dim_id,t_dim_id/), &
                                   temp1_3d_var_id)
            endif
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var XYGs")
            nSTAT = nf90_put_att(ncid,temp1_3d_var_id,"long_name",&
                                 var_User3d_XYGs_lname(ivar))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att XYGs long_name")
            nSTAT = nf90_put_att(ncid,temp1_3d_var_id,'units',&
                                 var_User3d_XYGs_unit(ivar))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att XYGs units")
            nSTAT = nf90_put_att(ncid,temp1_3d_var_id,"missing_value",&
                                 var_User3d_XYGs_MissVal(ivar))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att XYGs missing_value")
            nSTAT = nf90_put_att(ncid,temp1_3d_var_id,"_FillValue",&
                                 var_User3d_XYGs_FillVal(ivar))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att XYGs _FillValue")
          enddo
        endif

        !! Define User-specified 3-d transient variables in x,y,z
        if(nvar_User3d_XYZ.gt.0)then
          do ivar=1,nvar_User3d_XYZ
            if(op.eq.8)then
              nSTAT = nf90_def_var(ncid,var_User3d_XYZ_name(ivar),&
                                   nf90_double,  &
                                   (/x_dim_id,y_dim_id,z_dim_id,t_dim_id/), &
                                   temp1_3d_var_id)
            else
              nSTAT = nf90_def_var(ncid,var_User3d_XYZ_name(ivar),&
                                   nf90_float,  &
                                   (/x_dim_id,y_dim_id,z_dim_id,t_dim_id/), &
                                   temp1_3d_var_id)
            endif
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var XYZ")
            nSTAT = nf90_put_att(ncid,temp1_3d_var_id,"long_name",&
                                 var_User3d_XYZ_lname(ivar))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att XYZ long_name")
            nSTAT = nf90_put_att(ncid,temp1_3d_var_id,'units',&
                                 var_User3d_XYZ_unit(ivar))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att XYZ units")
            nSTAT = nf90_put_att(ncid,temp1_3d_var_id,"missing_value",&
                                 var_User3d_XYZ_MissVal(ivar))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att XYZ missing_value")
            nSTAT = nf90_put_att(ncid,temp1_3d_var_id,"_FillValue",&
                                 var_User3d_XYZ_FillVal(ivar))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att XYZ _FillValue")
          enddo
        endif

        ! Define User-specified 4-d transient variables in x,y,z,gs
        if(nvar_User4d_XYZGs.gt.0)then
          do ivar=1,nvar_User4d_XYZGs
            if(op.eq.8)then
              nSTAT = nf90_def_var(ncid,var_User4d_XYZGs_name(ivar),&
                                   nf90_double,  &
                                   (/x_dim_id,y_dim_id,z_dim_id,bn_dim_id,t_dim_id/), &
                                   temp1_4d_var_id)
            else
              nSTAT = nf90_def_var(ncid,var_User4d_XYZGs_name(ivar),&
                                   nf90_float, &
                                   (/x_dim_id,y_dim_id,z_dim_id,bn_dim_id,t_dim_id/), &
                                   temp1_4d_var_id)
            endif
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var XYZGs")
            nSTAT = nf90_put_att(ncid,temp1_4d_var_id,"long_name",&
                                 var_User4d_XYZGs_lname(ivar))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att XYZGs long_name")
            nSTAT = nf90_put_att(ncid,temp1_4d_var_id,'units',&
                                 var_User4d_XYZGs_unit(ivar))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att XYZGs units")
            nSTAT = nf90_put_att(ncid,temp1_4d_var_id,"missing_value",&
                                 var_User4d_XYZGs_MissVal(ivar))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att XYZGs missing_value")
            nSTAT = nf90_put_att(ncid,temp1_4d_var_id,"_FillValue",&
                                 var_User4d_XYZGs_FillVal(ivar))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att XYZGs _FillValue")
          enddo
        endif

      endif  ! USE_OPTMOD_VARS`
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                    Leaving define mode.                               !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      nSTAT = nf90_enddef(ncid)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"enddef")
      if(VERB.gt.1)write(*,*)"Leaving define mode"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Fill variables with initial values
        ! Time
      if(VERB.gt.1)write(global_info,*)"     Fill time"
      dumscal_out = real(time,kind=op)
      nSTAT=nf90_put_var(ncid,t_var_id,dumscal_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var time")
        ! Z
      if(VERB.gt.1)write(global_info,*)"     Fill Z"
      allocate(dum1d_out(nzmax))
      dum1d_out(:) = real(z_cc_pd(1:nzmax),kind=op)
      nSTAT=nf90_put_var(ncid,z_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var z")
      deallocate(dum1d_out)
        ! Y
      if(VERB.gt.1)write(global_info,*)"     Fill Y"
      allocate(dum1d_out(nymax))
      if (IsLatLon) then
        dum1d_out(1:nymax) = real(lat_cc_pd(1:nymax),kind=op)
      else
        dum1d_out(1:nymax) = real(y_cc_pd(1:nymax),kind=op)
      endif
      nSTAT=nf90_put_var(ncid,y_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var y")
      deallocate(dum1d_out)
        ! X
      if(VERB.gt.1)write(global_info,*)"     Fill X"
      allocate(dum1d_out(nxmax))
      if (IsLatLon) then
        if(lon_cc_pd(nxmax).gt.360.0_op)then
          dum1d_out(1:nxmax) = real(lon_cc_pd(1:nxmax),kind=op)-360.0_op
        else
          dum1d_out(1:nxmax) = real(lon_cc_pd(1:nxmax),kind=op)
        endif
      else
        dum1d_out(1:nxmax) = real(x_cc_pd(1:nxmax),kind=op)
      endif
      nSTAT=nf90_put_var(ncid,x_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var x")
      deallocate(dum1d_out)
        ! BN (Grain size bin ID)
      if(VERB.gt.1)write(global_info,*)"     Fill BN"
      allocate(dum1dint_out(nsmax))
      do i=1,nsmax  ! This is just an index
        dum1dint_out(i) = i
      enddo
      nSTAT=nf90_put_var(ncid,bn_var_id,dum1dint_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var bn")
      deallocate(dum1dint_out)
        ! ER
      if(VERB.gt.1)write(global_info,*)"     Fill ER"
      allocate(dum1dint_out(neruptions))
      ! This is variable associated with the dimension for eruptions ID
      ! This only contains the index starting with 1
      do i=1,neruptions
        dum1dint_out(i) = i
      enddo
      nSTAT=nf90_put_var(ncid,er_var_id,dum1dint_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var er")
      deallocate(dum1dint_out)
        ! WS
      if(VERB.gt.1)write(global_info,*)"     Fill WF"
      allocate(dum1dint_out(MR_iwindfiles))
      ! This is variable associated with the dimension for windfile ID
      ! This only contains the index starting with 1
      do i=1,MR_iwindfiles
        dum1dint_out(i) = i
      enddo
      nSTAT=nf90_put_var(ncid,wf_var_id,dum1dint_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var wf")
      deallocate(dum1dint_out)
        ! PT
      if (nairports.gt.0)then
        if(VERB.gt.1)write(global_info,*)"     Fill PT"
        allocate(dum1dint_out(nairports))
        ! This is variable associated with the dimension for point output (airport/POI)
        ! This only contains the index starting with 1
        do i=1,nairports
          dum1dint_out(i) = i
        enddo
        nSTAT=nf90_put_var(ncid,pt_var_id,dum1dint_out,(/1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pt")
        deallocate(dum1dint_out)
      endif

      if(NCversion.eq.4.and.nvprofiles.gt.0)then
          ! PR
        if(VERB.gt.1)write(global_info,*)"     Fill PR"
        allocate(dum1dint_out(nvprofiles))
        ! This is variable associated with the dimension for profile output
        ! This only contains the index starting with 1
        do i=1,nvprofiles
          dum1dint_out(i) = i
        enddo
        nSTAT=nf90_put_var(ncid,pr_var_id,dum1dint_out,(/1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pr")
        deallocate(dum1dint_out)
          ! TN
        !if(VERB.gt.1)write(global_info,*)"     Fill TN"
        !dumscal_out = real(time,kind=op)
        !nSTAT=nf90_put_var(ncid,tn_var_id,dumscal_out,(/1/))
        !if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var tn")
      endif

      !   End of filling dimension variables
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !   Now fill a few other variables that are a function of BN
         ! Species class ID
      if(VERB.gt.1)write(global_info,*)"     Fill Species class ID"
      allocate(dum1dint_out(nsmax))
      dum1dint_out(1:nsmax) = SpeciesID(1:nsmax)
      nSTAT=nf90_put_var(ncid,spec_var_id,dum1dint_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var spec_class")
      deallocate(dum1dint_out)
         ! Species sub-class ID
      if(VERB.gt.1)write(global_info,*)"     Fill Species sub-class ID"
      allocate(dum1dint_out(nsmax))
      dum1dint_out(1:nsmax) = SpeciesSubID(1:nsmax)
      nSTAT=nf90_put_var(ncid,subspec_var_id,dum1dint_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var spec_subclass")
      deallocate(dum1dint_out)
         ! Grain-size diameter
      if(VERB.gt.1)write(global_info,*)"     Fill GS Diameter"
      allocate(dum1d_out(nsmax))
      dum1d_out = 0.0_op
      if(useCalcFallVel)then
        dum1d_out(1:n_gs_max) = real(Tephra_gsdiam(1:n_gs_max),kind=op)
      else
        do i=1,n_gs_max
          dum1d_out(i) = real(i,kind=op)
        enddo
      endif
      nSTAT=nf90_put_var(ncid,gssd_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var gs_diameter")
         ! gs_massfrac (Mass fraction of grain size)
      if(VERB.gt.1)write(global_info,*)"     Fill GS MassFrac"
      dum1d_out = 0.0_op
      dum1d_out(1:n_gs_max) = real(Tephra_bin_mass(1:n_gs_max),kind=op)
      nSTAT=nf90_put_var(ncid,gsmf_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var gs_massfrac")
       ! gs_dens (Density of grain)
      if(useCalcFallVel)then
        dum1d_out(1:n_gs_max) = real(Tephra_rho_m(1:n_gs_max),kind=op)
      else
        dum1d_out = 0.0_op
      endif
      nSTAT=nf90_put_var(ncid,gsdens_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var gs_dens")
      deallocate(dum1d_out)

      !   Now fill a few other variables that are a function of ER
        ! er_stime (Start time of eruption)
      allocate(dum1d_out(neruptions))
      dum1d_out = real(e_StartTime + SimStartHour,kind=op)
      nSTAT=nf90_put_var(ncid,er_stime_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var er_stime")
        ! er_duration (Duration of eruption)
      dum1d_out = real(e_Duration,kind=op)
      nSTAT=nf90_put_var(ncid,er_duration_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var er_duration")
        ! er_plumeheight (Plume height of eruption)
      dum1d_out = real(e_PlumeHeight,kind=op)
      nSTAT=nf90_put_var(ncid,er_plumeheight_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var er_plumeheight")
        ! er_volume (Volume of eruption)
      dum1d_out = real(e_Volume,kind=op)
      nSTAT=nf90_put_var(ncid,er_volume_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var er_volume")
      deallocate(dum1d_out)

         ! Now fill the other (non-time-dependent) variables
         ! wf_name (Name of windfile)
      do i=1,MR_iwindfiles
        write(lllinebuffer,'(130a)')MR_windfiles(i)
        do j=1,130
          nSTAT=nf90_put_var(ncid,wf_name_var_id,lllinebuffer(j:j),(/j,i/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var wf_name")
        enddo
      enddo

         ! Cell area
      if(VERB.gt.1)write(global_info,*)"     Fill area"
      if (IsLatLon) then
        dum2d_out(1:nxmax,1:nymax) = real(sigma_nz_pd(1:nxmax,1:nymax,0),kind=op)
      else
        dum2d_out(1:nxmax,1:nymax) = real(dx*dy,kind=op)
      endif
      nSTAT=nf90_put_var(ncid,area_var_id,dum2d_out)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var area")

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(USE_WIND_VARS)then
          ! Vz
        if(VERB.gt.1)write(global_info,*)"     Fill Vz"
        dum3d_out(1:nxmax,1:nymax,1:nzmax) = real(vz_pd(1:nxmax,1:nymax,1:nzmax),kind=op)
        nSTAT=nf90_put_var(ncid,vz_var_id,dum3d_out,(/1,1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var vz")
          ! Vy
        if(VERB.gt.1)write(global_info,*)"     Fill Vy"
        dum3d_out(1:nxmax,1:nymax,1:nzmax) = real(vy_pd(1:nxmax,1:nymax,1:nzmax),kind=op)
        nSTAT=nf90_put_var(ncid,vy_var_id,dum3d_out,(/1,1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var vy")
          ! Vx
        if(VERB.gt.1)write(global_info,*)"     Fill Vx"
        dum3d_out(1:nxmax,1:nymax,1:nzmax) = real(vx_pd(1:nxmax,1:nymax,1:nzmax),kind=op)
        nSTAT=nf90_put_var(ncid,vx_var_id,dum3d_out,(/1,1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var vx")
          ! Vf (Note, using ashcon as a dummy variable for this 4-d output)
        ashcon(1:nxmax,1:nymax,1:nzmax,1:nsmax) = &
             real(vf_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax),kind=op)
        nSTAT=nf90_put_var(ncid,vf_var_id,ashcon,(/1,1,1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var vf")
      endif  ! USE_WIND_VARS
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(USE_RESTART_VARS)then
        if(VERB.gt.1)write(global_info,*)"     Fill ashcon"
          ! ashcon
          ! netCDF standard requires that the unlimited dimension (time)
          ! be the most slowly varying, so dum4d_out unfortunately has a
          ! different shape than concen
        do i = 1,nsmax
          ashcon(1:nxmax,1:nymax,1:nzmax,i) = &
             real(concen_pd(1:nxmax,1:nymax,1:nzmax,i,ts1),kind=op)
        enddo
        nSTAT=nf90_put_var(ncid,ashcon_var_id,ashcon,(/1,1,1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var ashcon")
      endif  ! USE_RESTART_VARS
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(USE_OUTPROD_VARS)then
          ! depocon
        if(VERB.gt.1)write(global_info,*)"     Fill depocon"
        depocon = 0.0_op
          ! netCDF standard requires that the unlimited dimension (time)
          ! be the most slowly varying, so dum3d_out unfortunately has a
          ! different shape than depocon
        do i = 1,n_gs_max
            ! Here's the conversion to kg/m^2 from kg/km^2
          depocon(1:nxmax,1:nymax,i) = real(DepositGranularity(1:nxmax,1:nymax,i) * &
                                         dz_vec_pd(0)/KM2_2_M2,kind=op)
        enddo
        do n=1,n_gs_max
          do i=1,nxmax
            do j=1,nymax
              if(depocon(i,j,n).le.EPS_SMALL)depocon(i,j,n)=0.0_op
            enddo
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,depocon_var_id,depocon,(/1,1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var depocon")
  
        if(VERB.gt.1)write(global_info,*)"     Calling dbZCalculator"
        call dbZCalculator            ! get radar reflectivity
  
        ! depothick
        if(VERB.gt.1)write(global_info,*)"     Fill depothick"
        dum2d_out(:,:) = DepositThickness_FillValue
        do i=1,nxmax
          do j=1,nymax
            if(Mask_Deposit(i,j))&
                dum2d_out(i,j) = real(DepositThickness(i,j),kind=op)
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,depothick_var_id,dum2d_out,(/1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var depothick")
 
          ! depotime
        if(VERB.gt.1)write(global_info,*)"     Fill depotime"
        dum2d_out(:,:) = DepArrivalTime_FillValue
        do i=1,nxmax
          do j=1,nymax
            if(DepArrivalTime(i,j).ge.0.0_ip)&
                dum2d_out(i,j)=real(DepArrivalTime(i,j),kind=op)
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,depotime_var_id,dum2d_out,(/1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var depotime")

          ! ashtime
        if(VERB.gt.1)write(global_info,*)"     Fill ashtime"
        dum2d_out(:,:) = CloudArrivalTime_FillValue
        do i=1,nxmax
          do j=1,nymax
            if(CloudArrivalTime(i,j).ge.0.0)&
                dum2d_out(i,j)=real(CloudArrivalTime(i,j),kind=op)
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,ashcloudtime_var_id,dum2d_out,(/1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var ashcloudtime_")
  
        ! ashconMax
        if(VERB.gt.1)write(global_info,*)"     Fill ashconMax"
        dum2d_out(:,:) = real(MaxConcentration_FillValue,kind=op)
        do i=1,nxmax
          do j=1,nymax
            if(Mask_Cloud(i,j))then
              dumscal_out=real(MaxConcentration(i,j)/1.0e3_ip,kind=op)
              dum2d_out(i,j)=dumscal_out
            endif
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,ashconMax_var_id,dum2d_out,(/1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var ashconMax")
  
        ! ash cloud_height (top)
        if(VERB.gt.1)write(global_info,*)"     Fill ash_height"
        dum2d_out(:,:) = real(MaxHeight_FillValue,kind=op)
        do i=1,nxmax
          do j=1,nymax
            !if(MaxHeight(i,j).gt.0.0_ip)then
            if(Mask_Cloud(i,j))then
              dumscal_out=real(MaxHeight(i,j),kind=op)
              dum2d_out(i,j) = dumscal_out
            endif
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,ashheight_var_id,dum2d_out,(/1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var ashheight")

        ! ash-load
        if(VERB.gt.1)write(global_info,*)"     Fill ashload"
        dum2d_out(:,:) = real(CloudLoad_FillValue,kind=op)
        do i=1,nxmax
          do j=1,nymax
            if(Mask_Cloud(i,j))then
              dumscal_out=real(CloudLoad(i,j),kind=op)
              dum2d_out(i,j)=dumscal_out
            endif
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,ashload_var_id,dum2d_out,(/1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var ashload")

        ! radar-reflectivity
        if(VERB.gt.1)write(global_info,*)"     Fill dbZ"
        dum3d_out(:,:,:) = real(DBZ_THRESH,kind=op)
        do i=1,nxmax
          do j=1,nymax
            do k=1,nzmax
              if(dbZ(i,j,k).ge.DBZ_THRESH)then
                dumscal_out=real(dbZ(i,j,k),kind=op)
                dum3d_out(i,j,k)=dumscal_out
              endif
            enddo
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,radrefl_var_id,dum3d_out,(/1,1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var reflectivity")

        ! ash cloud_bottom
        if(VERB.gt.1)write(global_info,*)"     Fill ash height min"
        dum2d_out(:,:) = real(MinHeight_FillValue,kind=op)
        do i=1,nxmax
          do j=1,nymax
            if(Mask_Cloud(i,j))then
              dumscal_out=real(MinHeight(i,j),kind=op)
              dum2d_out(i,j) = dumscal_out
            endif
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,ashcloudBot_var_id,dum2d_out,(/1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var cloud-bottom")

      endif ! USE_OUTPROD_VARS

      if (nairports.gt.0)then
        ! These are variable associated with the dimension for point output
        ! (airport/POI)
        if(VERB.gt.1)write(global_info,*)"     Fill PT variables"

        ! x coordinate of point
        allocate(dum1d_out(nairports))
        dum1d_out(1:nairports) = real(Airport_Longitude(1:nairports),kind=op)
        nSTAT=nf90_put_var(ncid,pt_x_var_id,dum1d_out,(/1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pt_x")
        deallocate(dum1d_out)

        ! y coordinate of point
        allocate(dum1d_out(nairports))
        dum1d_out(1:nairports) = real(Airport_Latitude(1:nairports),kind=op)
        nSTAT=nf90_put_var(ncid,pt_y_var_id,dum1d_out,(/1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pt_y")
        deallocate(dum1d_out)

        ! Point/Airport code (3-char label)
        do i=1,nairports
          write(lllinebuffer,'(130a)')Airport_Code(i)
          do j=1,130
            nSTAT=nf90_put_var(ncid,pt_code_var_id,lllinebuffer(j:j),(/j,i/))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pt_code")
          enddo
        enddo

        ! Point/Airport name
        do i=1,nairports
          write(lllinebuffer,'(130a)')Airport_Name(i)
          do j=1,130
            nSTAT=nf90_put_var(ncid,pt_name_var_id,lllinebuffer(j:j),(/j,i/))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pt_name")
          enddo
        enddo

        ! Time-series of ash accumulation at points
        allocate(dum1d_out(nairports))
        dum1d_out(1:nairports) = real(Airport_Thickness_TS(1:nairports,1),kind=op)
        nSTAT=nf90_put_var(ncid,pt_ashthickness_var_id,dum1d_out,(/1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pt_depothick")
        deallocate(dum1d_out)

      endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(USE_OPTMOD_VARS)then
        if(VERB.gt.1)write(global_info,*)"     Fill OPTMOD_VARS"

        ! Fill User-specified 2-d static variables
        if(nvar_User2d_static_XY.gt.0)then
          do ivar=1,nvar_User2d_static_XY
            dum2d_out(1:nxmax,1:nymax) = real(var_User2d_static_XY(1:nxmax,1:nymax,ivar),kind=op)
            nSTAT = nf90_inq_varid(ncid,var_User2d_static_XY_name(ivar),temp1_2d_var_id)
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid static_XY")
            nSTAT = nf90_put_var(ncid,temp1_2d_var_id,dum2d_out,(/1,1/))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var static_XY")
          enddo
        endif

        ! Fill User-specified 2-d transient variables
        if(nvar_User2d_XY.gt.0)then
          do ivar=1,nvar_User2d_XY
            dum2d_out(:,:) = real(var_User2d_XY(:,:,ivar),kind=op)
            nSTAT = nf90_inq_varid(ncid,var_User2d_XY_name(ivar),temp1_2d_var_id)
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_var XY")
            nSTAT = nf90_put_var(ncid,temp1_2d_var_id,dum2d_out,(/1,1,1/))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var XY")
          enddo
        endif

        ! Fill User-specified 3-d transient variables in x,y,gs
        if(nvar_User3d_XYGs.gt.0)then
          do ivar=1,nvar_User3d_XYGs
            depocon(:,:,:) = real(var_User3d_XYGs(:,:,:,ivar),kind=op)
            nSTAT = nf90_inq_varid(ncid,var_User3d_XYGs_name(ivar),temp1_3d_var_id)
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_var XYGs")
            nSTAT = nf90_put_var(ncid,temp1_3d_var_id,depocon,(/1,1,1,1/))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var XYGs")
          enddo
        endif

        ! Fill User-specified 3-d transient variables in x,y,z
        if(nvar_User3d_XYZ.gt.0)then
          do ivar=1,nvar_User3d_XYZ
            dum3d_out(:,:,:) = real(var_User3d_XYZ(:,:,:,ivar),kind=op)
            nSTAT = nf90_inq_varid(ncid,var_User3d_XYZ_name(ivar),temp1_3d_var_id)
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_var XYZ")
            nSTAT = nf90_put_var(ncid,temp1_3d_var_id,dum3d_out,(/1,1,1,1/))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var XYZ")
          enddo
        endif

        ! Fill User-specified 4-d transient variables in x,y,z,gs
        if(nvar_User4d_XYZGs.gt.0)then
          do ivar=1,nvar_User4d_XYZGs
            ashcon(:,:,:,:) = real(var_User4d_XYZGs(1:nxmax,1:nymax,1:nzmax,1:nsmax,ivar),kind=op)
            nSTAT = nf90_inq_varid(ncid,var_User4d_XYZGs_name(ivar),temp1_4d_var_id)
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_var XYZGs")
            nSTAT = nf90_put_var(ncid,temp1_4d_var_id,ashcon,(/1,1,1,1,1/))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var XYZGs")
          enddo
        endif

      endif  ! USE_OPTMOD_VARS
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Close file
      nSTAT = nf90_close(ncid)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"nf90_close")

      deallocate(dum2d_out,dum3d_out)
      deallocate(ashcon)
      deallocate(depocon)
      if(VERB.gt.1)write(global_info,*)"Created netcdf file"

      end subroutine create_netcdf_file
!##############################################################################


!##############################################################################
!
!    append_to_netcdf
!
!##############################################################################
      subroutine append_to_netcdf

      use global_param,  only : &
         EPS_SMALL,EPS_TINY,VERB,KM2_2_M2

      use io_data,       only : &
         iout3d,nvar_User2d_XY,nvar_User3d_XYGs,nvar_User3d_XYZ,nvar_User4d_XYZGs,&
         outfile,isFinal_TS,nvprofiles

      use Output_Vars,   only : &
         CloudLoad_FillValue,dbZCol_FillValue,MaxConcentration_FillValue,&
         MaxHeight_FillValue,MinHeight_FillValue,Mask_Cloud,Mask_Deposit,&
         DepositThickness_FillValue,CloudArrivalTime_FillValue,DepArrivalTime_FillValue,&
         var_User2d_XY_name,var_User2d_XY,var_User3d_XYGs_name,var_User3d_XYGs,&
         var_User3d_XYZ_name,var_User3d_XYZ,var_User4d_XYZGs_name,var_User4d_XYZGs,&
         DBZ_THRESH,USE_OPTMOD_VARS,USE_RESTART_VARS,&
         USE_OUTPROD_VARS,USE_WIND_VARS,DepositThickness,DepArrivalTime,CloudArrivalTime,&
         MaxConcentration,MaxHeight,CloudLoad,dbZ,MinHeight,pr_ash,&
           dbZCalculator

      use Tephra,        only : &
         n_gs_max

      use Airports,      only : &
         nairports,Airport_Thickness_TS,Airport_AshDuration,Airport_CloudArrivalTime,&
         Airport_CloudDuration,Airport_AshArrivalTime

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,ts1,dz_vec_pd

      use solution,      only : &
         vx_pd,vy_pd,vz_pd,vf_pd,concen_pd,DepositGranularity

      use time_data,     only : &
         time,ntmax,time_native

      implicit none

      integer :: i,j,k,n

      integer :: nSTAT
      integer :: ivar
!      integer :: ncid
!      integer :: t_var_id          = 0
!      integer :: tn_var_id         = 0
!      integer :: vx_var_id         = 0
!      integer :: vy_var_id         = 0
!      integer :: vz_var_id         = 0
!      integer :: vf_var_id         = 0
!      integer :: ashcon_var_id     = 0
!      integer :: depocon_var_id    = 0
!
!      integer :: ashconMax_var_id      = 0 ! Max Ash concentration in column
!      integer :: ashheight_var_id      = 0 ! Height of top of ash cloud
!      integer :: ashload_var_id        = 0 ! Vert. integrated load of ash cloud
!      integer :: radrefl_var_id        = 0 ! Radar reflectivity in dbZ
!      integer depothick_var_id     = 0 ! Total deposit (summed over gs)
!      integer :: depotime_var_id       = 0 ! Deposit arrival time
!      integer :: ashcloudtime_var_id   = 0 ! Cloud arrival time
!      integer :: ashcloudBot_var_id    = 0 ! Height of bottom of ash cloud
!      integer :: pt_ashthickness_var_id = 0 ! TS of ashfall
!      integer :: pt_cloudarrival_var_id = 0
!      integer :: pt_cloudduration_var_id= 0
!      integer :: pt_asharrival_var_id   = 0
!      integer :: pt_ashduration_var_id  = 0
!      integer :: pr_ash_var_id          = 0
!
!      integer :: temp1_2d_var_id !,temp2_2d_var_id,temp3_2d_var_id,temp4_2d_var_id
!      integer :: temp1_3d_var_id !,temp2_3d_var_id,temp3_3d_var_id,temp4_3d_var_id
!      integer :: temp1_4d_var_id !,temp2_4d_var_id!,temp3_4d_var_id!,temp4_4d_var_id
!      integer :: ns_extra
!
!      real(kind=op)                                 :: dumscal_out
!      real(kind=op), dimension(:)      ,allocatable :: dum1d_out
!      real(kind=op), dimension(:,:)    ,allocatable :: dum2d_out
!      real(kind=op), dimension(:,:,:)  ,allocatable :: dum3d_out
!
!      real(kind=op), dimension(:,:,:,:),allocatable :: ashcon
!      real(kind=op), dimension(:,:,:)  ,allocatable :: depocon

      if(VERB.gt.2)write(global_info,*)"Allocating output vars"
      !allocate(dum2d_out(nxmax,nymax))
      !allocate(dum3d_out(nxmax,nymax,nzmax))

      allocate(ashcon(nxmax,nymax,nzmax,nsmax))
      allocate(depocon(nxmax,nymax,nsmax))
      if(nsmax.gt.n_gs_max)then
        ! These are some non-ash species
        ns_extra = nsmax-n_gs_max
      else
        ns_extra = 0
      endif

      ! Open netcdf file for writing
      if(VERB.gt.2)write(global_info,*)"opening netcdf file"
      nSTAT=nf90_open(outfile,nf90_write, ncid)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"nf90_open")

      if(isFinal_TS)then
          ! depotime
        allocate(dum2d_out(nxmax,nymax))
        nSTAT = nf90_inq_varid(ncid,"depotime",depotime_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid depotime")
        if(VERB.gt.2)write(global_info,*)"  Writing depotime"
        dum2d_out(:,:) = DepArrivalTime_FillValue
        do i=1,nxmax
          do j=1,nymax
            if(DepArrivalTime(i,j).ge.0.0_ip)&
                dum2d_out(i,j)=real(DepArrivalTime(i,j),kind=op)
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,depotime_var_id,dum2d_out,(/1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var depotime")
        deallocate(dum2d_out)

          ! ashtime
        allocate(dum2d_out(nxmax,nymax))
        nSTAT = nf90_inq_varid(ncid,"ash_arrival_time",ashcloudtime_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid ash_arrival_time")
        if(VERB.gt.2)write(global_info,*)"  Writing ashtime"
        dum2d_out(:,:) = CloudArrivalTime_FillValue
        do i=1,nxmax
          do j=1,nymax
            if(CloudArrivalTime(i,j).ge.0.0_ip)&
                dum2d_out(i,j)=real(CloudArrivalTime(i,j),kind=op)
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,ashcloudtime_var_id,dum2d_out,(/1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var ash_arrival_time")
        deallocate(dum2d_out)

        if (nairports.gt.0)then

          ! Arrival time of ashfall
          nSTAT = nf90_inq_varid(ncid,"pt_depotime",pt_asharrival_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid pt_depotime")
          allocate(dum1d_out(nairports))
          dum1d_out(1:nairports) = real(Airport_AshArrivalTime(1:nairports),kind=op)
          nSTAT=nf90_put_var(ncid,pt_asharrival_var_id,dum1d_out,(/1/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pt_depotime")
          deallocate(dum1d_out)

          ! Duration of ashfall
          nSTAT = nf90_inq_varid(ncid,"pt_depodur",pt_ashduration_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid pt_depodur")
          allocate(dum1d_out(nairports))
          dum1d_out(1:nairports) = real(Airport_AshDuration(1:nairports),kind=op)
          nSTAT=nf90_put_var(ncid,pt_ashduration_var_id,dum1d_out,(/1/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pt_depodur")
          deallocate(dum1d_out)

          ! Arrival time of ash cloud
          nSTAT = nf90_inq_varid(ncid,"pt_cloud_arrival",pt_cloudarrival_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid pt_cloud_arrival")
          allocate(dum1d_out(nairports))
          dum1d_out(1:nairports) = real(Airport_CloudArrivalTime(1:nairports),kind=op)
          nSTAT=nf90_put_var(ncid,pt_cloudarrival_var_id,dum1d_out,(/1/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pt_cloud_arrival")
          deallocate(dum1d_out)

          ! Duration of ash cloud
          nSTAT = nf90_inq_varid(ncid,"pt_cloud_dur",pt_cloudduration_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid pt_cloud_dur")
          allocate(dum1d_out(nairports))
          dum1d_out(1:nairports) = real(Airport_CloudDuration(1:nairports),kind=op)
          nSTAT=nf90_put_var(ncid,pt_cloudduration_var_id,dum1d_out,(/1/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pt_cloud_dur")
          deallocate(dum1d_out)

        endif

        if(nvprofiles.gt.0)then
          ! PR
          ! write out native time
          nSTAT = nf90_inq_varid(ncid,"tn",tn_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid tn")
          allocate(dum1d_out(ntmax))
          dum1d_out(1:ntmax) = real(time_native(1:ntmax),kind=op)
          nSTAT=nf90_put_var(ncid,tn_var_id,dum1d_out,(/1/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var tn")
          deallocate(dum1d_out)

          ! write out profiles
          nSTAT = nf90_inq_varid(ncid,"pr_ash",pr_ash_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid pr_ash")
          allocate(dum3d_out(nzmax,ntmax,nvprofiles))
          dum3d_out(1:nzmax,1:ntmax,1:nvprofiles) = &
                    real(pr_ash(1:nzmax,1:ntmax,1:nvprofiles),kind=op)
          nSTAT=nf90_put_var(ncid,pr_ash_var_id,dum3d_out,(/1,1,1/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pr_ash")
          deallocate(dum3d_out)
          write(*,*)"WROTE :",nvprofiles,nzmax,ntmax
        endif

      endif

      if(.not.isFinal_TS)then
        ! Get variable ids
        nSTAT = nf90_inq_varid(ncid,"t",t_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid t")
  
        if(VERB.gt.2)write(global_info,*)"Got var IDs, now writing data"
        ! Write data
        ! Time
        if(VERB.gt.2)write(global_info,*)"  Writing Time"
        dumscal_out = real(time,kind=op)
        nSTAT=nf90_put_var(ncid,t_var_id,dumscal_out,(/iout3d/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var t")
  
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(USE_WIND_VARS)then
          allocate(dum3d_out(nxmax,nymax,nzmax))
            ! Vz
          nSTAT = nf90_inq_varid(ncid,"vz",vz_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid vz")
          if(VERB.gt.2)write(global_info,*)"  Writing Vz"
          dum3d_out(1:nxmax,1:nymax,1:nzmax) = real(vz_pd(1:nxmax,1:nymax,1:nzmax),kind=op)
          nSTAT=nf90_put_var(ncid,vz_var_id,dum3d_out,(/1,1,1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var vz")
            ! Vy
          nSTAT = nf90_inq_varid(ncid,"vy",vy_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid vy")
          if(VERB.gt.2)write(global_info,*)"  Writing Vy"
          dum3d_out(1:nxmax,1:nymax,1:nzmax) = real(vy_pd(1:nxmax,1:nymax,1:nzmax),kind=op)
          nSTAT=nf90_put_var(ncid,vy_var_id,dum3d_out,(/1,1,1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var vy")
            ! Vx
          nSTAT = nf90_inq_varid(ncid,"vx",vx_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid vx")
          if(VERB.gt.2)write(global_info,*)"  Writing Vx"
          dum3d_out(1:nxmax,1:nymax,1:nzmax) = real(vx_pd(1:nxmax,1:nymax,1:nzmax),kind=op)
          nSTAT=nf90_put_var(ncid,vx_var_id,dum3d_out,(/1,1,1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var vx")
          deallocate(dum3d_out)
            ! Vf
          nSTAT = nf90_inq_varid(ncid,"vf",vf_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid vf")
          ashcon(1:nxmax,1:nymax,1:nzmax,1:nsmax) = real(vf_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax),kind=op)
          nSTAT=nf90_put_var(ncid,vf_var_id,ashcon,(/1,1,1,1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var vf")
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(USE_RESTART_VARS)then
            ! ashcon
            ! netCDF standard requires that the unlimited dimension (time)
            ! be the most slowly varying, so ashcon unfortunately has a
            ! different shape than concen
          nSTAT = nf90_inq_varid(ncid,"ashcon",ashcon_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid ashcon")
          if(VERB.gt.2)write(global_info,*)"  Writing ashcon"
          ashcon = 0.0_op
          ashcon(1:nxmax,1:nymax,1:nzmax,1:nsmax) = real(concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts1),kind=op)
          nSTAT=nf90_put_var(ncid,ashcon_var_id,ashcon,(/1,1,1,1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var ashcon")
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(USE_OUTPROD_VARS)then
  
            ! depocon
          if(VERB.gt.1)write(global_info,*)"     Fill depocon"
          nSTAT = nf90_inq_varid(ncid,"depocon",depocon_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid depocon")
          depocon = 0.0_op
            ! netCDF standard requires that the unlimited dimension (time)
            ! be the most slowly varying, so dum3d_out unfortunately has a
            ! different shape than depocon
          do i = 1,n_gs_max
              ! Here's the conversion to kg/m^2 from kg/km^2
            depocon(1:nxmax,1:nymax,i) = real(DepositGranularity(1:nxmax,1:nymax,i) * &
                                           dz_vec_pd(0)/KM2_2_M2,kind=op)
          enddo
          do n=1,n_gs_max
            do i=1,nxmax
              do j=1,nymax
                if(depocon(i,j,n).le.EPS_SMALL)depocon(i,j,n)=0.0_op
              enddo
            enddo
          enddo
          nSTAT=nf90_put_var(ncid,depocon_var_id,depocon,(/1,1,1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var depocon")
    
          if(VERB.gt.2)write(global_info,*)"Calling dbZCalculator"
          call dbZCalculator            ! get radar reflectivity
    
          ! depothick
          allocate(dum2d_out(nxmax,nymax))
          nSTAT = nf90_inq_varid(ncid,"depothick",depothick_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid depothick")
          if(VERB.gt.2)write(global_info,*)"  Writing depothick"
          dum2d_out(:,:) = DepositThickness_FillValue
          do i=1,nxmax
            do j=1,nymax
              if(Mask_Deposit(i,j))&
                  dum2d_out(i,j)=real(DepositThickness(i,j),kind=op)
            enddo
          enddo
          nSTAT=nf90_put_var(ncid,depothick_var_id,dum2d_out,(/1,1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var depothick")
          deallocate(dum2d_out)

          ! ashconMax
          allocate(dum2d_out(nxmax,nymax))
          dum2d_out(:,:) = MaxConcentration_FillValue
          nSTAT = nf90_inq_varid(ncid,"ashcon_max",ashconMax_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid ashcon_max")
          if(VERB.gt.2)write(global_info,*)"  Writing ashconMax"
          dum2d_out(:,:) = real(MaxConcentration,kind=op)
          do i=1,nxmax
            do j=1,nymax
              if(Mask_Cloud(i,j))then
                dumscal_out=real(MaxConcentration(i,j)/1.0e3_ip,kind=op)
                dum2d_out(i,j)=dumscal_out
              endif
            enddo
          enddo
          nSTAT=nf90_put_var(ncid,ashconMax_var_id,dum2d_out,(/1,1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var ashcon_max")
          deallocate(dum2d_out)
   
          ! ash cloud_height
          allocate(dum2d_out(nxmax,nymax))
          nSTAT = nf90_inq_varid(ncid,"cloud_height",ashheight_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid cloud_height")
          if(VERB.gt.2)write(global_info,*)"  Writing ash height"
          dum2d_out(:,:) = real(MaxHeight_FillValue,kind=op)
          do i=1,nxmax
            do j=1,nymax
              if(Mask_Cloud(i,j))then
                dumscal_out=real(MaxHeight(i,j),kind=op)
                dum2d_out(i,j) = dumscal_out
              endif
            enddo
          enddo
          nSTAT=nf90_put_var(ncid,ashheight_var_id,dum2d_out,(/1,1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var cloud_height")
          deallocate(dum2d_out)

          ! ash-load
          allocate(dum2d_out(nxmax,nymax))
          nSTAT = nf90_inq_varid(ncid,"cloud_load",ashload_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid cloud_load")
          if(VERB.gt.2)write(global_info,*)"  Writing ash-load"
          dum2d_out(:,:) = real(CloudLoad_FillValue,kind=op)
          do i=1,nxmax
            do j=1,nymax
              if(Mask_Cloud(i,j))then
                dumscal_out=real(CloudLoad(i,j),kind=op)
                dum2d_out(i,j) = dumscal_out
              endif
            enddo
          enddo
          nSTAT=nf90_put_var(ncid,ashload_var_id,dum2d_out,(/1,1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var cloud_load")
          deallocate(dum2d_out)

          ! radar reflectivity
          allocate(dum3d_out(nxmax,nymax,nzmax))
          nSTAT = nf90_inq_varid(ncid,"radar_reflectivity",radrefl_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid radar_reflectivity")
          if(VERB.gt.2)write(global_info,*)"  Writing radar reflectivity"
          dum3d_out(:,:,:) = real(dbZCol_FillValue,kind=op)
          do i=1,nxmax
            do j=1,nymax
              do k=1,nzmax
                if(dbZ(i,j,k).ge.DBZ_THRESH)then
                  dumscal_out=real(dbZ(i,j,k),kind=op)
                  dum3d_out(i,j,k)=dumscal_out
                endif
              enddo
            enddo
          enddo
          nSTAT=nf90_put_var(ncid,radrefl_var_id,dum3d_out,(/1,1,1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var radar_reflectivity")
          deallocate(dum3d_out)

          ! ash cloud_bottom
          allocate(dum2d_out(nxmax,nymax))
          nSTAT = nf90_inq_varid(ncid,"cloud_bottom",ashcloudBot_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid cloud_bottom")
          if(VERB.gt.2)write(global_info,*)"  Writing ash-height (bottom)"
          dum2d_out(:,:) = real(MinHeight_FillValue,kind=op)
          do i=1,nxmax
            do j=1,nymax
              if(Mask_Cloud(i,j))then
                dumscal_out=real(MinHeight(i,j),kind=op)
                dum2d_out(i,j) = dumscal_out
              endif
            enddo
          enddo
          nSTAT=nf90_put_var(ncid,ashcloudBot_var_id,dum2d_out,(/1,1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var cloud_bottom")
          deallocate(dum2d_out)

        endif ! USE_OUTPROD_VARS

        if (nairports.gt.0)then
          ! Time-series of ash accumulation at points
          nSTAT = nf90_inq_varid(ncid,"pt_depothick",pt_ashthickness_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid pt_depothick")
          if(VERB.gt.2)write(global_info,*)"  Writing time-series of ashfall at points"
          allocate(dum1d_out(nairports))
          dum1d_out(1:nairports) = real(Airport_Thickness_TS(1:nairports,iout3d),kind=op)
          nSTAT=nf90_put_var(ncid,pt_ashthickness_var_id,dum1d_out,(/1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pt_depothick")
          deallocate(dum1d_out)
        endif
      endif ! .not.isFinal_TS
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(USE_OPTMOD_VARS)then
        if(isFinal_TS)then
          ! User-specified 2-d static variables were filled in the
          ! create_netcdf subroutine
          write(*,*)"No final variables for optional moduals"
        endif

        if(.not.isFinal_TS)then
          ! Fill User-specified 2-d transient variables
          if(nvar_User2d_XY.gt.0)then
            allocate(dum2d_out(nxmax,nymax))
            do ivar=1,nvar_User2d_XY
              dum2d_out(:,:) = real(var_User2d_XY(:,:,ivar),kind=op)
              nSTAT = nf90_inq_varid(ncid,var_User2d_XY_name(ivar),temp1_2d_var_id)
              if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid XY")
              nSTAT = nf90_put_var(ncid,temp1_2d_var_id,dum2d_out,(/1,1,iout3d/))
              if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var XY")
            enddo
            deallocate(dum2d_out)
          endif

          ! Fill User-specified 3-d transient variables in x,y,gs
          if(nvar_User3d_XYGs.gt.0)then
            do ivar=1,nvar_User3d_XYGs
              depocon(:,:,:) = real(var_User3d_XYGs(:,:,:,ivar),kind=op)
              nSTAT = nf90_inq_varid(ncid,var_User3d_XYGs_name(ivar),temp1_3d_var_id)
              if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid XYGs")
              nSTAT = nf90_put_var(ncid,temp1_3d_var_id,depocon,(/1,1,1,iout3d/))
              if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var XYGs")
            enddo
          endif

          ! Fill User-specified 3-d transient variables in x,y,z
          if(nvar_User3d_XYZ.gt.0)then
            allocate(dum3d_out(nxmax,nymax,nzmax))
            do ivar=1,nvar_User3d_XYZ
              dum3d_out(:,:,:) = real(var_User3d_XYZ(:,:,:,ivar),kind=op)
              nSTAT = nf90_inq_varid(ncid,var_User3d_XYZ_name(ivar),temp1_3d_var_id)
              if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid XYZ")
              nSTAT = nf90_put_var(ncid,temp1_3d_var_id,dum3d_out,(/1,1,1,iout3d/))
              if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var XYZ")
            enddo
            deallocate(dum3d_out)
          endif

          ! Fill User-specified 4-d transient variables in x,y,z,gs
          if(nvar_User4d_XYZGs.gt.0)then
            do ivar=1,nvar_User4d_XYZGs
              ashcon(:,:,:,:) = real(var_User4d_XYZGs(1:nxmax,1:nymax,1:nzmax,1:nsmax,ivar),kind=op)
              nSTAT = nf90_inq_varid(ncid,var_User4d_XYZGs_name(ivar),temp1_4d_var_id)
              if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid XYZGs")
              nSTAT = nf90_put_var(ncid,temp1_4d_var_id,ashcon,(/1,1,1,1,iout3d/))
              if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var XYZGs")
            enddo
          endif
        endif
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Close file
      if(VERB.gt.2)write(global_info,*)"closing file"
      nSTAT = nf90_close(ncid)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"nf90_close")

      if(VERB.gt.2)write(global_info,*)"Deallocating"
      deallocate(ashcon)
      deallocate(depocon)

      end subroutine append_to_netcdf

!##############################################################################
!
!    NC_RestartFile_ReadTimes
!
!##############################################################################

      subroutine NC_RestartFile_ReadTimes

      use io_data,           only : &
         concenfile,init_tstep

      use time_data,         only : &
         time

      implicit none

      integer :: nSTAT
      integer :: ncid

      integer :: t_dim_id
      integer :: t_var_id
      integer :: t_len
      integer :: it
      !real(kind=op)                                 :: dumscal_out
      real(kind=op), allocatable, dimension(:) :: t_list

      ! Open netcdf file for writing
      nSTAT=nf90_open(concenfile,nf90_nowrite,ncid)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"nSTAT=nf90_open")

      ! Identify dimension for time (and note size)
      nSTAT = nf90_inq_dimid(ncid,"t",t_dim_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_dimid time")
      nSTAT = nf90_Inquire_Dimension(ncid,t_dim_id,len=t_len)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"Inquire_Dimension time")

      ! Get variable ids
      nSTAT = nf90_inq_varid(ncid,"t",t_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid t")

      allocate(t_list(t_len))
      nSTAT = nf90_get_var(ncid,t_var_id,t_list)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var t")

      write(*,*)"  Step :  time"
      do it = 1,t_len
        write(*,*)it,t_list(it)
      enddo

      write(global_info,*)'Enter timestep for initialization'
      read(5,*) init_tstep

      nSTAT=nf90_get_var(ncid,t_var_id,dumscal_out,(/init_tstep/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var t")
      time = real(dumscal_out,kind=ip)

      deallocate(t_list)

      end subroutine NC_RestartFile_ReadTimes

!##############################################################################
!
!    NC_RestartFile_LoadConcen
!
!##############################################################################

      subroutine NC_RestartFile_LoadConcen

      use global_param,      only : &
         KM2_2_M2

      use io_data,           only : &
         concenfile,init_tstep

      !use Tephra,        only : &
      !   n_gs_max

      use mesh,              only : &
         nxmax,nymax,nzmax,nsmax,dz_vec_pd,ts0,ts1

      use solution,          only : &
         concen_pd

      implicit none

      integer :: i
      integer :: nSTAT
!      integer :: ncid

!      integer :: t_var_id          = 0
!      integer :: vx_var_id         = 0
!      integer :: vy_var_id         = 0
!      integer :: vz_var_id         = 0
!      integer :: ashcon_var_id     = 0
!      integer :: depocon_var_id  = 0

!      real(kind=op)                                 :: dumscal_out
!      real(kind=op), dimension(:,:)    ,allocatable :: dum2d_out
!      real(kind=op), dimension(:,:,:)  ,allocatable :: dum3d_out

!      real(kind=op), dimension(:,:,:,:),allocatable :: ashcon
!      real(kind=op), dimension(:,:,:)  ,allocatable :: depocon

      allocate(dum2d_out(nxmax,nymax))
      allocate(dum3d_out(nxmax,nymax,nzmax))
      allocate(ashcon(nxmax,nymax,nzmax,nsmax))
      allocate(depocon(nxmax,nymax,nsmax))

      write(global_info,*)"WARNING "
      write(global_info,*)"Input file is not currently verified "
      write(global_info,*)" with previous run."

      ! Open netcdf file for writing
      nSTAT=nf90_open(concenfile,nf90_nowrite,ncid)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"nf90_open")

      ! Get variable ids
      nSTAT = nf90_inq_varid(ncid,"t",t_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid t")
      nSTAT = nf90_inq_varid(ncid,"vx",vx_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid vx")
      nSTAT = nf90_inq_varid(ncid,"vy",vy_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid vy")
      nSTAT = nf90_inq_varid(ncid,"vz",vz_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid vz")
      nSTAT = nf90_inq_varid(ncid,"ashcon",ashcon_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid ashcon")
      nSTAT = nf90_inq_varid(ncid,"depocon",depocon_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid depocon")

      nSTAT=nf90_get_var(ncid,ashcon_var_id,ashcon,  &
               start = (/1,1,1,1,init_tstep/),       &
               count = (/nxmax,nymax,nzmax,nsmax,1/))

      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var ashcon")

      do i = 1,nsmax
        concen_pd(1:nxmax,1:nymax,1:nzmax,i,ts1) = real(ashcon(:,:,:,i),kind=ip)
      enddo

      nSTAT=nf90_get_var(ncid,depocon_var_id,depocon,&
               start = (/1,1,1,init_tstep/),       &
               count = (/nxmax,nymax,nsmax,1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var depocon")

      do i = 1,nsmax
          ! Here's the conversion from kg/m^2
        concen_pd(1:nxmax,1:nymax,0,i,ts1) = real(depocon(:,:,i),kind=ip)/(dz_vec_pd(0)/KM2_2_M2)
      enddo
      concen_pd(:,:,:,:,ts0) = concen_pd(:,:,:,:,ts1)

      ! Close file
      nSTAT = nf90_close(ncid)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"nf90_close")

      deallocate(dum2d_out,dum3d_out)
      deallocate(ashcon)
      deallocate(depocon)

      write(global_info,*)"Read concentrations from time ",dumscal_out

      end subroutine NC_RestartFile_LoadConcen


!##############################################################################
!
!    NC_check_status
!
!    nSTAT   = error code returned from netcdf call
!    errcode = user-supplied return value on stopping of code
!    operation = string descriptor of function call causing error
!
!    Error-checking routine for NetCDF function calls.
!    Modeled after a subroutine posted at:
!    https://climate-cms.org/2018/10/12/create-netcdf.html
!
!##############################################################################

      subroutine NC_check_status(nSTAT, errcode, operation)

      use io_units

      implicit none

      integer, intent(in) :: nSTAT
      integer, intent(in) :: errcode
      character(len=*), intent(in) :: operation

      character(len=9) :: severity

      if (errcode.eq.0)then 
        severity = "WARNING: "
       else
        severity = "ERROR:   "
      endif

      if (nSTAT == nf90_noerr) return
      !write(global_essential ,*)severity,errcode,operation,nf90_strerror(nSTAT)
      write(global_log ,*)severity,errcode,operation,nf90_strerror(nSTAT)
      write(global_error ,*)severity,errcode,operation,nf90_strerror(nSTAT)

      ! If user-supplied error code is 0, then consider this a warning,
      ! otherwise do a hard stop
      if (errcode.ne.0) stop

      end subroutine NC_check_status

!##############################################################################
!
!    NC_Read_Output_Products
!
!    if timestep = -1, then use the last step in file
!##############################################################################

      subroutine NC_Read_Output_Products(timestep)

      use io_data,           only : &
         concenfile,init_tstep

      use time_data,         only : &
         time

      use Output_Vars,   only : &
         DepositThickness,DepArrivalTime,CloudArrivalTime,&
         MaxConcentration,MaxHeight,CloudLoad,dbZ,MinHeight,Mask_Cloud,Mask_Deposit

      use Airports,      only : &
         nairports,Airport_Code,Airport_Name,Airport_x,Airport_y,&
         Airport_Latitude,Airport_Longitude,Airport_Thickness_TS

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,x_cc_pd,y_cc_pd,z_cc_pd,lon_cc_pd,lat_cc_pd,&
         dx,dy,dz_vec_pd,IsLatLon

      implicit none

      integer, intent(in), optional :: timestep

      integer :: nSTAT

      integer :: it
      real(kind=op) :: dumscal_out
      real(kind=op), allocatable, dimension(:) :: t_list
      real(kind=op), dimension(:,:)    ,allocatable :: dum2d_out

      write(global_info,*)"Reading NetCDF file ",concenfile
      write(global_info,*)"Found the following dimensions and sizes"

      ! Open netcdf file for writing
      nSTAT=nf90_open(concenfile,nf90_nowrite,ncid)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"nSTAT=nf90_open")

      !!!!!!  TIME  !!!!!!!!!!!
      ! Identify dimension for time (and note size)
      nSTAT = nf90_inq_dimid(ncid,"t",t_dim_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_dimid time")
      nSTAT = nf90_Inquire_Dimension(ncid,t_dim_id,len=t_len)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"Inquire_Dimension time")
      ! Get variable id for this dimension
      nSTAT = nf90_inq_varid(ncid,"t",t_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid t")
      write(*,*)"t",t_len

      !!!!!!  X  !!!!!!!!!!!
      ! Identify dimension for x (and note size)
      nSTAT = nf90_inq_dimid(ncid,"x",x_dim_id)
      !call NC_check_status(nSTAT,0,"inq_dimid x")
      if(nSTAT.ne.0)then
        ! dimension x not found, trying for lon
        !write(*,*)"Could not find dimension x"
        !write(*,*)"Trying for lon"
        nSTAT = nf90_inq_dimid(ncid,"lon",x_dim_id)
        call NC_check_status(nSTAT,1,"inq_dimid x")
        if (nSTAT.eq.0)then
          IsLatLon = .true.
        endif
       else
        IsLatLon = .false.
      endif
      ! continue reading x or lon
      nSTAT = nf90_Inquire_Dimension(ncid,x_dim_id,len=x_len)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"Inquire_Dimension x")
      if(IsLatLon)then
        ! Get variable id for this dimension
        nSTAT = nf90_inq_varid(ncid,"lon",x_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid x")
        if(.not.allocated(lon_cc_pd))then
          nxmax = x_len
          allocate(lon_cc_pd(-1:nxmax+2))
          allocate(dum1d_out(1:nxmax))
          nSTAT=nf90_get_var(ncid,x_var_id,dum1d_out,  &
                   start = (/1/),       &
                   count = (/x_len/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var x")
          lon_cc_pd(1:nxmax) = real(dum1d_out(1:nxmax),kind=ip)
          deallocate(dum1d_out)
        else
          write(*,*)"ERROR: lon_cc_pd already allocated"
        endif
        write(*,*)"lon",x_len
      else
        ! Get variable id for this dimension
        nSTAT = nf90_inq_varid(ncid,"x",x_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid x")
        if(.not.allocated(x_cc_pd))then
          nxmax = x_len
          allocate(x_cc_pd(-1:nxmax+2))
          allocate(dum1d_out(1:nxmax))
          nSTAT=nf90_get_var(ncid,x_var_id,dum1d_out,  &
                   start = (/1/),       &
                   count = (/x_len/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var x")
          x_cc_pd(1:nxmax) = real(dum1d_out(1:nxmax),kind=ip)
          deallocate(dum1d_out)
        else
          write(*,*)"ERROR: x_cc_pd already allocated"
        endif
        write(*,*)"x",x_len
      endif

      !!!!!!  Y  !!!!!!!!!!!
      ! Identify dimension for y (and note size)
      if(IsLatLon)then
        nSTAT = nf90_inq_dimid(ncid,"lat",y_dim_id)
        call NC_check_status(nSTAT,1,"inq_dimid x")
        nSTAT = nf90_Inquire_Dimension(ncid,y_dim_id,len=y_len)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"Inquire_Dimension y")
        ! Get variable id for this dimension
        nSTAT = nf90_inq_varid(ncid,"lat",y_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid y")
        if(.not.allocated(lat_cc_pd))then
          nymax = y_len
          allocate(lat_cc_pd(-1:nymax+2))
          allocate(dum1d_out(1:nymax))
          nSTAT=nf90_get_var(ncid,y_var_id,dum1d_out,  &
                   start = (/1/),       &
                   count = (/y_len/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var y")
          lat_cc_pd(1:nymax) = real(dum1d_out(1:nymax),kind=ip)
          deallocate(dum1d_out)
        else
          write(*,*)"ERROR: lat_cc_pd already allocated"
        endif
        write(*,*)"lat",y_len
      else
        nSTAT = nf90_inq_dimid(ncid,"y",y_dim_id)
        call NC_check_status(nSTAT,1,"inq_dimid x")
        nSTAT = nf90_Inquire_Dimension(ncid,y_dim_id,len=y_len)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"Inquire_Dimension y")
        ! Get variable id for this dimension
        nSTAT = nf90_inq_varid(ncid,"y",y_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid y")
        if(.not.allocated(y_cc_pd))then
          nymax = y_len
          allocate(y_cc_pd(-1:nymax+2))
          allocate(dum1d_out(1:nymax))
          nSTAT=nf90_get_var(ncid,y_var_id,dum1d_out,  &
                   start = (/1/),       &
                   count = (/y_len/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var y")
          y_cc_pd(1:nymax) = real(dum1d_out(1:nymax),kind=ip)
          deallocate(dum1d_out)
        else
          write(*,*)"ERROR: lat_cc_pd already allocated"
        endif
        write(*,*)"y",y_len
      endif

      !!!!!!  Z  !!!!!!!!!!!
      ! Identify dimension for z (and note size)
      nSTAT = nf90_inq_dimid(ncid,"z",z_dim_id)
      call NC_check_status(nSTAT,1,"inq_dimid z")
      nSTAT = nf90_Inquire_Dimension(ncid,z_dim_id,len=z_len)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"Inquire_Dimension z")
      ! Get variable id for this dimension
      nSTAT = nf90_inq_varid(ncid,"z",z_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid z")
      write(*,*)"z",z_len

      !!!!!!  BN  !!!!!!!!!!!
      ! Identify dimension for bn (and note size)
      nSTAT = nf90_inq_dimid(ncid,"bn",bn_dim_id)
      call NC_check_status(nSTAT,1,"inq_dimid bn")
      nSTAT = nf90_Inquire_Dimension(ncid,bn_dim_id,len=bn_len)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"Inquire_Dimension bn")
      ! Get variable id for this dimension
      nSTAT = nf90_inq_varid(ncid,"bn",bn_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid bn")
      write(*,*)"bn",bn_len

      !!!!!!  PT  !!!!!!!!!!!
      ! Identify dimension for pt (and note size)
      nSTAT = nf90_inq_dimid(ncid,"pt",pt_dim_id)
      call NC_check_status(nSTAT,1,"inq_dimid bn")
      nSTAT = nf90_Inquire_Dimension(ncid,pt_dim_id,len=pt_len)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"Inquire_Dimension pt")
      ! Get variable id for this dimension
      nSTAT = nf90_inq_varid(ncid,"pt",pt_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid pt")
      write(*,*)"pt",pt_len

      !!!!!!  PR  !!!!!!!!!!!
      ! Identify dimension for pr (and note size)
      nSTAT = nf90_inq_dimid(ncid,"pr",pr_dim_id)
      call NC_check_status(nSTAT,1,"inq_dimid pr")
      nSTAT = nf90_Inquire_Dimension(ncid,pr_dim_id,len=pr_len)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"Inquire_Dimension pr")
      ! Get variable id for this dimension
      nSTAT = nf90_inq_varid(ncid,"pr",pr_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid pr")
      write(*,*)"pr",pr_len

      !!!!!!  TN !!!!!!!!!!!
      ! Identify dimension for tn(and note size)
      nSTAT = nf90_inq_dimid(ncid,"tn",tn_dim_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_dimid tn")
      nSTAT = nf90_Inquire_Dimension(ncid,tn_dim_id,len=tn_len)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"Inquire_Dimension tn")
      ! Get variable id for this dimension
      nSTAT = nf90_inq_varid(ncid,"tn",tn_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid tn")
      write(*,*)"tn",tn_len

      ! Now get all the other variable info:
      ! Vx
      nSTAT = nf90_inq_varid(ncid,"vx",vx_var_id)
      !if(nSTAT.ne.0)call NC_check_status(nSTAT,0,"inq_varid vx")
      ! Vy
      nSTAT = nf90_inq_varid(ncid,"vy",vy_var_id)
      !if(nSTAT.ne.0)call NC_check_status(nSTAT,0,"inq_varid vy")
      ! Vz
      nSTAT = nf90_inq_varid(ncid,"vz",vz_var_id)
      !if(nSTAT.ne.0)call NC_check_status(nSTAT,0,"inq_varid vz")
      ! ashcon
      nSTAT = nf90_inq_varid(ncid,"ashcon",ashcon_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,0,"inq_varid ashcon")
      ! depocon
      nSTAT = nf90_inq_varid(ncid,"depocon",depocon_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,0,"inq_varid depocon")
      ! ashcon_max
      nSTAT = nf90_inq_varid(ncid,"ashcon_max",ashconMax_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,0,"inq_varid ashconMax")
      ! cloud_height
      nSTAT = nf90_inq_varid(ncid,"cloud_height",ashheight_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,0,"inq_varid cloud_height")
      ! cloud_load
      nSTAT = nf90_inq_varid(ncid,"cloud_load",ashload_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,0,"inq_varid cloud_load")
      ! cloud_bottom
      nSTAT = nf90_inq_varid(ncid,"cloud_bottom",ashcloudBot_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,0,"inq_varid cloud_bottom")
      ! depothick
      nSTAT = nf90_inq_varid(ncid,"depothick",depothick_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,0,"inq_varid depothick")

      ! pt_x
      nSTAT = nf90_inq_varid(ncid,"pt_x",pt_x_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,0,"inq_varid pt_x")
      ! pt_y
      nSTAT = nf90_inq_varid(ncid,"pt_y",pt_y_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,0,"inq_varid pt_y")
      ! pt_code
      nSTAT = nf90_inq_varid(ncid,"pt_code",pt_code_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,0,"inq_varid pt_code")
      ! pt_name
      nSTAT = nf90_inq_varid(ncid,"pt_name",pt_name_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,0,"inq_varid pt_name")
      ! pt_depotime
      nSTAT = nf90_inq_varid(ncid,"pt_depotime",pt_asharrival_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,0,"inq_varid pt_depotime")
      ! pt_depodur
      nSTAT = nf90_inq_varid(ncid,"pt_depodur",pt_ashduration_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,0,"inq_varid pt_depod")
      ! pt_cloud_arrival
      nSTAT = nf90_inq_varid(ncid,"pt_cloud_arrival",pt_cloudarrival_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,0,"inq_varid pt_cloud_arrival")
      ! pt_cloud_dur
      nSTAT = nf90_inq_varid(ncid,"pt_cloud_dur",pt_cloudduration_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,0,"inq_varid pt_cloud_dur")
      ! pt_depothick
      nSTAT = nf90_inq_varid(ncid,"pt_depothick",pt_ashthickness_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,0,"inq_varid pt_depothick")

      allocate(t_list(t_len))
      nSTAT = nf90_get_var(ncid,t_var_id,t_list)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var t")

      if (.not.present(timestep))then
        write(*,*)"Found the following time steps in file:"
        write(*,*)"  Step :  time"
        do it = 1,t_len
          write(*,*)it,t_list(it)
        enddo
        write(global_info,*)'Enter timestep for initialization'
        read(5,*) init_tstep
      else
        if (timestep.eq.-1)then
          init_tstep = t_len
        else
          init_tstep = timestep
        endif
      endif
      if(init_tstep.gt.t_len)then
        write(*,*)"ERROR: Requested time step index is greater than available."
        stop 1
      elseif(init_tstep.lt.0)then
        write(*,*)"ERROR: Requested time step index is invalid."
        stop 1
      endif
      write(*,*)"Requested timestep = ",init_tstep,t_list(init_tstep)


      nSTAT=nf90_get_var(ncid,t_var_id,dumscal_out,(/init_tstep/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var t")
      time = real(dumscal_out,kind=ip)

      allocate(dum2d_out(x_len,y_len))
      nSTAT=nf90_get_var(ncid,depothick_var_id,dum2d_out,  &
               start = (/1,1,init_tstep/),       &
               count = (/x_len,y_len,1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var depothick")

      allocate(DepositThickness(x_len,y_len))
      DepositThickness = real(dum2d_out,kind=ip)
      deallocate(t_list,dum2d_out)

      end subroutine NC_Read_Output_Products

      end module Ash3d_Netcdf

