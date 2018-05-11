!##############################################################################
!
!    create_netcdf_file
!
!##############################################################################

      subroutine create_netcdf_file

      use precis_param

      use io_units

      use global_param,  only : &
         EPS_SMALL,EPS_TINY,KM2_2_M2,useCalcFallVel,VERB

      use io_data,       only : &
         cdf_b1l1,cdf_b1l2,cdf_b1l3,cdf_b1l4,cdf_b1l5,cdf_b1l6,cdf_b1l7,cdf_b1l8,cdf_b1l9,&
         cdf_b3l1,cdf_b3l2,cdf_b3l3,cdf_b3l4,cdf_b3l5,cdf_b4l1,cdf_b4l2,cdf_b4l3,cdf_b4l4,&
         cdf_b4l5,cdf_b4l6,cdf_b4l7,cdf_b4l8,cdf_b4l9,cdf_b4l10,cdf_b4l11,cdf_b6l1,cdf_b6l2,&
         cdf_b6l3,cdf_b6l4,cdf_comment,cdf_title,outfile,&
         nvar_User2d_static_XY,nvar_User2d_XY,nvar_User3d_XYGs,nvar_User3d_XYZ,&
         nvar_User4d_XYZGs

      use Output_Vars,   only : &
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
         CLOUDLOAD_THRESH,DBZ_THRESH,USE_OPTMOD_VARS,USE_RESTART_VARS,&
         USE_OUTPROD_VARS,USE_WIND_VARS,DepositThickness,DepArrivalTime,CloudArrivalTime,&
         MaxConcentration,MaxHeight,CloudLoad,dbZ,MinHeight,&
           dbZCalculator

      use Tephra,        only : &
         n_gs_max,Tephra_gsdiam,Tephra_v_s,Tephra_bin_mass,Tephra_rho_m

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,x_cc_pd,y_cc_pd,z_cc_pd,lon_cc_pd,lat_cc_pd,&
         sigma_nz_pd,dx,dy,dz_vec_pd,IsLatLon,ts1

      use solution,      only : &
          vx_pd,vy_pd,vz_pd,vf_pd,concen_pd,DepositGranularity,SpeciesID,SpeciesSubID

      use time_data,     only : &
          BaseYear,useLeap,cdf_time_log,time,SimStartHour,xmlSimStartTime

      use Source,        only : &
         neruptions,e_Volume,e_Duration,e_StartTime,PlumeHeight

      use MetReader,     only : &
         MR_iwindfiles,MR_windfiles,MR_MetStep_Hour_since_baseyear

      use netcdf

      implicit none

      integer :: nSTAT
      integer :: ncid
      integer :: t_dim_id     = 0 ! Time
      integer :: x_dim_id     = 0 ! X
      integer :: y_dim_id     = 0 ! Y
      integer :: z_dim_id     = 0 ! Z
      integer :: bin_dim_id   = 0 ! Full generalized species class ID       1:nsmax    (replaces gs_dim_id)
      !integer :: gs_dim_id    = 0 ! index just for grain size bins          1:n_gs_max
      !integer :: Xspec_dim_id = 0 ! Index for extra (non-ash) species bins  1:ns_extra

      integer :: er_dim_id    = 0 ! eruption number
      integer :: wf_dim_id    = 0 ! windfile number
      integer :: sl_dim_id    = 0 ! string length
      !integer :: pj_dim_id    = 0 ! projection parameter dimension

      integer :: t_var_id              = 0 ! Time
      integer :: x_var_id              = 0 ! X-distance
      integer :: y_var_id              = 0 ! Y-distance
      integer :: z_var_id              = 0 ! Z-distance
      integer :: bin_var_id            = 0 ! Species class ID
      !integer :: gs_var_id             = 0 ! grain size
      integer :: spec_var_id           = 0 ! Species class ID
      integer :: subspec_var_id        = 0 ! Species sub-class ID
      integer :: er_var_id             = 0 ! eruption index
      integer :: wf_var_id             = 0 ! wind file index
      integer :: vx_var_id             = 0 ! Vx
      integer :: vy_var_id             = 0 ! Vy
      integer :: vz_var_id             = 0 ! Vz
      integer :: vf_var_id             = 0 ! Vf
      integer :: ashcon_var_id         = 0 ! Ash concentration
      integer :: depocon_var_id        = 0 ! Deposit mass/area
      !integer :: gencon_var_id         = 0 ! General concentration for any slices of concen above the # of GS
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

      integer :: temp1_2d_var_id !,temp2_2d_var_id,temp3_2d_var_id,temp4_2d_var_id
      integer :: temp1_3d_var_id !,temp2_3d_var_id,temp3_3d_var_id,temp4_3d_var_id
      integer :: temp1_4d_var_id !,temp2_4d_var_id!,temp3_4d_var_id!,temp4_4d_var_id
      integer :: ns_extra

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
      character (len=20)  :: HS_xmltime


      ! Since output precision might be different from input precision,
      ! we need to allocate to correct memory space
      real(kind=op)                                 :: dumscal_out
      integer,       dimension(:)      ,allocatable :: dum1dint_out
      real(kind=op), dimension(:)      ,allocatable :: dum1d_out
      real(kind=op), dimension(:,:)    ,allocatable :: dum2d_out
      real(kind=op), dimension(:,:,:)  ,allocatable :: dum3d_out

      real(kind=op), dimension(:,:,:,:),allocatable :: ashcon
      real(kind=op), dimension(:,:,:)  ,allocatable :: depocon
      !real(kind=op), dimension(:,:,:,:),allocatable :: gencon

      character(len=3) ,dimension(10) :: dim_names
      character(len=30),dimension(40) :: var_lnames
      character (len=13)         :: reftimestr
      character (len=13)         ::  HS_yyyymmddhhmm_since
      character(len=130):: lllinebuffer

      integer :: i,j,k,n
      integer :: ivar

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
      dim_names(2) = "z"
      if (IsLatLon) then
        dim_names(3) = "lat"
        dim_names(4) = "lon"
      else
        dim_names(3) = "y"
        dim_names(4) = "x"
      endif
      dim_names(5) = "bn" ! grainsize bin for 1:n_gs_max, then generalized bin up to nsmax
      dim_names(6) = "er" ! eruption index
      dim_names(7) = "wf" ! windfile index
      dim_names(8) = "sl" ! string length
      !dim_names(9) = "sc" ! full species class ID
      !dim_names(10)= "xs" ! index for extra species bins

      var_lnames(1) = "Time"
      var_lnames(2) = "Elevation"
      if (IsLatLon) then
        var_lnames(3) = "Latitude"
        var_lnames(4) = "Longitude"
      else
        var_lnames(3) = "Y (North)"
        var_lnames(4) = "X (East)"
      endif
      var_lnames(5) = "Bin index; Grainsize,spec. ID"
      var_lnames(6) = "Eruption number"
      var_lnames(7) = "Wind file number"

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
      nSTAT = nf90_create(outfile,nf90_clobber, ncid)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: create outfile: ', &
                              nf90_strerror(nSTAT)

      ! Fill in header info
      if(VERB.gt.1)write(global_info,*)"Filling in header info"
      nSTAT = nf90_put_att(ncid,nf90_global,"Title",cdf_title)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att title: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"User",cdf_user)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att User: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"Date",cdf_time_log)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Date: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"NWPStartTime",cdf_WindStartTime)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att NWPStartTime: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"Host",cdf_host)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Host: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"CWD",cdf_cwd)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att CWD: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"Comment",cdf_comment)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
        ! Add lines copied from the input file
      nSTAT = nf90_put_att(ncid,nf90_global,"b1l1",cdf_b1l1)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b1l2",cdf_b1l2)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b1l3",cdf_b1l3)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b1l4",cdf_b1l4)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b1l5",cdf_b1l5)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b1l6",cdf_b1l6)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b1l7",cdf_b1l7)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b1l8",cdf_b1l8)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b1l9",cdf_b1l9)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)

      nSTAT = nf90_put_att(ncid,nf90_global,"b3l1",cdf_b3l1)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b3l2",cdf_b3l2)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b3l3",cdf_b3l3)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b3l4",cdf_b3l4)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b3l5",cdf_b3l5)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)

      nSTAT = nf90_put_att(ncid,nf90_global,"b4l1",cdf_b4l1)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b4l2",cdf_b4l2)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b4l3",cdf_b4l3)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b4l4",cdf_b4l4)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b4l5",cdf_b4l5)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b4l6",cdf_b4l6)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b4l7",cdf_b4l7)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b4l8",cdf_b4l8)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b4l9",cdf_b4l9)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b4l10",cdf_b4l10)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b4l11",cdf_b4l11)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)

      nSTAT = nf90_put_att(ncid,nf90_global,"b6l1",cdf_b6l1)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b6l2",cdf_b6l2)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b6l3",cdf_b6l3)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,nf90_global,"b6l4",cdf_b6l4)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att Comment: ',nf90_strerror(nSTAT)

      ! Define dimensions
        ! t,z,y,x
        !  or time, elev, lon, lat
        !  or record, level, y, x
        ! and gs (grain size)
        ! sc (full species class indes)
        ! xs (extra species index)
        ! er (eruption index)
        ! wf (wind file index)
        ! sl (string length for storing windfile names)
      if(VERB.gt.1)write(global_info,*)"Defining dimensions"
      nSTAT = nf90_def_dim(ncid,dim_names(1),nf90_unlimited,t_dim_id)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: ',nf90_strerror(nSTAT)
      nSTAT = nf90_def_dim(ncid,dim_names(2),nzmax,z_dim_id)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: ',nf90_strerror(nSTAT)
      nSTAT = nf90_def_dim(ncid,dim_names(3),nymax,y_dim_id)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: ',nf90_strerror(nSTAT)
      nSTAT = nf90_def_dim(ncid,dim_names(4),nxmax,x_dim_id)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: ',nf90_strerror(nSTAT)
      ! Note we specify the full nxmax, not just n_gs_max
      nSTAT = nf90_def_dim(ncid,dim_names(5),nsmax,bin_dim_id)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: ',nf90_strerror(nSTAT)
      nSTAT = nf90_def_dim(ncid,dim_names(6),neruptions,er_dim_id)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: ',nf90_strerror(nSTAT)
      nSTAT = nf90_def_dim(ncid,dim_names(7),MR_iwindfiles,wf_dim_id)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: ',nf90_strerror(nSTAT)
      nSTAT = nf90_def_dim(ncid,dim_names(8),130,sl_dim_id)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: ',nf90_strerror(nSTAT)

      ! Define coordinate variables
        ! X,Y,Z,time,gs
      ! Define time-dependent variables
        ! Vx,Vy,Vz
        ! ashconc,depdepth
      ! Define other variables
        ! gs_setvel, gs_massfrac, gs_dens
         ! Time
      if(VERB.gt.1)write(global_info,*)"Defining coordinate variables"
      if(VERB.gt.1)write(global_info,*)"     Time",dim_names(1)
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,dim_names(1),nf90_double,(/t_dim_id/),&
                             t_var_id)
      else
        nSTAT = nf90_def_var(ncid,dim_names(1),nf90_float,(/t_dim_id/), &
                             t_var_id)
      endif
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: def_var: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,t_var_id,"long_name",var_lnames(1))
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      write(time_units,4313) xmlSimStartTime
4313  format('hours since ',a20)
      nSTAT = nf90_put_att(ncid,t_var_id,"units",time_units)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      !if(runAsForecast)then
      !  call get_start_time_nc(1,itstart_year,itstart_month, &
      !                           itstart_day,filestart_hour)
      !  offset = hours_since_1900(itstart_year,itstart_month, &
      !                            itstart_day,filestart_hour)
      !  reftimestr = yyyymmddhhmm_since_1900(SimStartHour)
      !else
      !  reftimestr = yyyymmddhhmm_since_1900(e_StartTime(1))
      !endif
      reftimestr = HS_yyyymmddhhmm_since(SimStartHour,BaseYear,useLeap)

      nSTAT = nf90_put_att(ncid,t_var_id,"ReferenceTime",reftimestr)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: ReferenceTime: ',nf90_strerror(nSTAT)

         ! Z
      if(VERB.gt.1)write(global_info,*)"     Z",dim_names(2)
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,dim_names(2),nf90_double,(/z_dim_id/),&
                             z_var_id)
      else
        nSTAT = nf90_def_var(ncid,dim_names(2),nf90_float,(/z_dim_id/), &
                             z_var_id)
      endif
      !write(global_info,*)"z_var_id : ",z_var_id
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: def_var: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,z_var_id,"long_name",var_lnames(2))
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,z_var_id,"units","km")
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)

         ! Y
      if(VERB.gt.1)write(global_info,*)"     Y",dim_names(3)
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,dim_names(3),nf90_double,(/y_dim_id/),&
                             y_var_id)
      else
        nSTAT = nf90_def_var(ncid,dim_names(3),nf90_float,(/y_dim_id/), &
                             y_var_id)
      endif
      !write(global_info,*)"y_var_id : ",y_var_id
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: def_var: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,y_var_id,"long_name",var_lnames(3))
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      if (IsLatLon) then
        nSTAT = nf90_put_att(ncid,y_var_id,"units","degrees")
      else
        nSTAT = nf90_put_att(ncid,y_var_id,"units","km")
      endif
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)

         ! X
      if(VERB.gt.1)write(global_info,*)"     X",dim_names(4)
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,dim_names(4),nf90_double,(/x_dim_id/),&
                             x_var_id)
      else
        nSTAT = nf90_def_var(ncid,dim_names(4),nf90_float,(/x_dim_id/), &
                             x_var_id)
      endif
      !write(global_info,*)"x_var_id : ",x_var_id
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: def_var: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,x_var_id,"long_name",var_lnames(4))
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      if (IsLatLon) then
        nSTAT = nf90_put_att(ncid,x_var_id,"units","degrees")
      else
        nSTAT = nf90_put_att(ncid,x_var_id,"units","km")
      endif
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)

         ! BN (Grain size bin ID)
      if(VERB.gt.1)write(global_info,*)"    Bin",dim_names(5)
      nSTAT = nf90_def_var(ncid,dim_names(5),nf90_int,(/bin_dim_id/), &
                           bin_var_id)
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: def_var: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,bin_var_id,"long_name",var_lnames(5))
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,bin_var_id,"units","mm")
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,bin_var_id,"Comment",&
       "This is just an index")
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)

         ! ER (Eruption index)
      if(VERB.gt.1)write(*,*)"     ER",dim_names(6)
      nSTAT = nf90_def_var(ncid,dim_names(6),nf90_int,(/er_dim_id/),&
                               er_var_id)
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: def_var: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,er_var_id,"long_name",var_lnames(6))
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,er_var_id,"units","index")
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)

         ! WF (Wind file index)
      if(VERB.gt.1)write(*,*)"     WF",dim_names(7)
      nSTAT = nf90_def_var(ncid,dim_names(7),nf90_int,(/wf_dim_id/),&
                               wf_var_id)
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: def_var: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,wf_var_id,"long_name",var_lnames(7))
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,wf_var_id,"units","index")
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)

      !   End of dimension variables
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !   Now a few other variables that are a function of BN
         ! Species class ID
      if(VERB.gt.1)write(global_info,*)"     SC : Species Class"
      nSTAT = nf90_def_var(ncid,"spec_class",nf90_int,(/bin_dim_id/), &
                           spec_var_id)
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: def_var: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,spec_var_id,"long_name",var_lnames(38))
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,spec_var_id,"units","index")
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,spec_var_id,"Comment",&
       "1=ash")
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)

         ! Species sub-class ID
      if(VERB.gt.1)write(global_info,*)"     SSC : Species Sub-class"
      nSTAT = nf90_def_var(ncid,"spec_subclass",nf90_int,(/bin_dim_id/), &
                           subspec_var_id)
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: def_var: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,subspec_var_id,"long_name",var_lnames(39))
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,subspec_var_id,"units","ID")
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,subspec_var_id,"Comment",&
       "Non-ash species code")
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)

         ! Grain-size diameter
      if(VERB.gt.1)write(global_info,*)"     GS, grain-diameter"
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,"gs_diameter",nf90_double,(/bin_dim_id/),&
                             gssd_var_id)
      else
        nSTAT = nf90_def_var(ncid,"gs_diameter",nf90_float,(/bin_dim_id/), &
                             gssd_var_id)
      endif
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: def_var: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,gssd_var_id,"long_name",var_lnames(13))
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,gssd_var_id,"units","mm")
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)

         ! gs_massfrac (Mass fraction of grain size)
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,"gs_massfrac",nf90_double,(/bin_dim_id/),&
                             gsmf_var_id)
      else
        nSTAT = nf90_def_var(ncid,"gs_massfrac",nf90_float,(/bin_dim_id/), &
                             gsmf_var_id)
      endif
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: def_var: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,gsmf_var_id,"long_name",var_lnames(14))
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,gsmf_var_id,"units","fraction")
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)

         ! gs_dens (Density of grain)
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,"gs_dens",nf90_double,(/bin_dim_id/),&
                             gsdens_var_id)
      else
        nSTAT = nf90_def_var(ncid,"gs_dens",nf90_float,(/bin_dim_id/), &
                             gsdens_var_id)
      endif
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: def_var: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,gsdens_var_id,"long_name",var_lnames(20))
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,gsdens_var_id,"units","kg/m3")
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)

      !   Now a few other variables that are a function of ER
         ! er_stime (Start time of eruption)
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,"er_stime",nf90_double,(/er_dim_id/),&
                             er_stime_var_id)
      else
        nSTAT = nf90_def_var(ncid,"er_stime",nf90_float,(/er_dim_id/), &
                             er_stime_var_id)
      endif
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: def_var: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,er_stime_var_id,"long_name",var_lnames(15))
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,er_stime_var_id,"units", &
                           "hours since 1900")
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)

         ! er_duration (Duration of eruption)
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,"er_duration",nf90_double,(/er_dim_id/),&
                             er_duration_var_id)
      else
        nSTAT = nf90_def_var(ncid,"er_duration",nf90_float,(/er_dim_id/), &
                             er_duration_var_id)
      endif
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: def_var: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,er_duration_var_id,"long_name",var_lnames(16))
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,er_duration_var_id,"units", "hours")
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)

         ! er_plumeheight (Plume height of eruption)
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,"er_plumeheight",nf90_double,(/er_dim_id/),&
                             er_plumeheight_var_id)
      else
        nSTAT = nf90_def_var(ncid,"er_plumeheight",nf90_float,(/er_dim_id/), &
                             er_plumeheight_var_id)
      endif
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: def_var: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,er_plumeheight_var_id,"long_name",var_lnames(17))
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,er_plumeheight_var_id,"units", &
                           "km")
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)

         ! er_volume (Volume of eruption)
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,"er_volume",nf90_double,(/er_dim_id/),&
                             er_volume_var_id)
      else
        nSTAT = nf90_def_var(ncid,"er_volume",nf90_float,(/er_dim_id/), &
                             er_volume_var_id)
      endif
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: def_var: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,er_volume_var_id,"long_name",var_lnames(18))
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,er_volume_var_id,"units", &
                           "km3")
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)

         ! Now define the other (non-time-dependent) variables
         ! wf_name (Name of windfile)
      nSTAT = nf90_def_var(ncid,"wf_name",nf90_char, &
                           (/sl_dim_id,wf_dim_id/),wf_name_var_id)
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: def_var: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,wf_name_var_id,"long_name",var_lnames(19))
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,wf_name_var_id,"units", &
                           "string")
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)

         ! Cell area
      if(VERB.gt.1)write(global_info,*)"     area"
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,"area",nf90_double, &
            (/x_dim_id,y_dim_id/),                &
              area_var_id)
      else
        nSTAT = nf90_def_var(ncid,"area",nf90_float,  &
            (/x_dim_id,y_dim_id/),                &
              area_var_id)
      endif
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: def_var cell area: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,area_var_id,"long_name","Cell area")
      if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
      nSTAT = nf90_put_att(ncid,area_var_id,"units","km2")
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: put_att cell area: ',nf90_strerror(nSTAT)

      if(VERB.gt.1)write(global_info,*)"Defining variables"

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(USE_WIND_VARS)then
         ! Now define the time-dependent variables
         ! Vx
        if(VERB.gt.1)write(global_info,*)"     vx",var_lnames(8)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"vx",nf90_double,  &
              (/x_dim_id,y_dim_id,z_dim_id,t_dim_id/), &
                vx_var_id)
        else
          nSTAT = nf90_def_var(ncid,"vx",nf90_float,   &
              (/x_dim_id,y_dim_id,z_dim_id,t_dim_id/), &
                vx_var_id)
        endif
        if(nSTAT.ne.0)write(global_log ,*)'ERROR: def_var vx: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,vx_var_id,"long_name",var_lnames(8))
        if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,vx_var_id,"units","km/hr")
        if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att vx: ',nf90_strerror(nSTAT)
           ! Vy
        if(VERB.gt.1)write(global_info,*)"     vy",var_lnames(9)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"vy",nf90_double,  &
              (/x_dim_id,y_dim_id,z_dim_id,t_dim_id/), &
                vy_var_id)
        else
          nSTAT = nf90_def_var(ncid,"vy",nf90_float,   &
              (/x_dim_id,y_dim_id,z_dim_id,t_dim_id/), &
                vy_var_id)
        endif
        if(nSTAT.ne.0)write(global_log ,*)'ERROR: def_var vy: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,vy_var_id,"long_name",var_lnames(9))
        if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,vy_var_id,"units","km/hr")
        if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att vy: ',nf90_strerror(nSTAT)
           ! Vz
        if(VERB.gt.1)write(global_info,*)"     vz",var_lnames(10)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"vz",nf90_double,  &
              (/x_dim_id,y_dim_id,z_dim_id,t_dim_id/), &
                vz_var_id)
        else
          nSTAT = nf90_def_var(ncid,"vz",nf90_float,   &
              (/x_dim_id,y_dim_id,z_dim_id,t_dim_id/), &
                vz_var_id)
        endif
        if(nSTAT.ne.0)write(global_log ,*)'ERROR: def_var vz: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,vz_var_id,"long_name","Wind velocity (z)")
        if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,vz_var_id,"units","km/hr")
        if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att vz: ',nf90_strerror(nSTAT)

           ! Vf
        if(VERB.gt.1)write(global_info,*)"     vf","Fall Velocity"
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"vf",nf90_double,  &
              (/x_dim_id,y_dim_id,z_dim_id,bin_dim_id,t_dim_id/), &
                vf_var_id)
        else
          nSTAT = nf90_def_var(ncid,"vf",nf90_float,   &
              (/x_dim_id,y_dim_id,z_dim_id,bin_dim_id,t_dim_id/), &
                vf_var_id)
        endif
        if(nSTAT.ne.0)write(global_log ,*)'ERROR: def_var vf: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,vf_var_id,"long_name","Fall Velocity")
        if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,vf_var_id,"units","km/hr")
        if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att vf: ',nf90_strerror(nSTAT)
      endif  ! USE_WIND_VARS
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(USE_RESTART_VARS)then
         ! Full Ash Concentration
        if(VERB.gt.1)write(global_info,*)"     ashcon",var_lnames(11)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"ashcon",nf90_double, &
              (/x_dim_id,y_dim_id,z_dim_id,bin_dim_id,t_dim_id/),    &
                ashcon_var_id)
        else
          nSTAT = nf90_def_var(ncid,"ashcon",nf90_float,  &
              (/x_dim_id,y_dim_id,z_dim_id,bin_dim_id,t_dim_id/),    &
                ashcon_var_id)
        endif
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: def_var ashcon: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashcon_var_id,"long_name",var_lnames(11))
        if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashcon_var_id,"units","kg/km^3")
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att ashcon: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashcon_var_id,&
                 "missing_value", -9999.0)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att ashcon: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashcon_var_id,"_FillValue",-9999.0_op)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att ashcon: ',nf90_strerror(nSTAT)
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
          nSTAT = nf90_def_var(ncid,"depocon",nf90_double, &
              (/x_dim_id,y_dim_id,bin_dim_id,t_dim_id/),                &
                depocon_var_id)
        else
          nSTAT = nf90_def_var(ncid,"depocon",nf90_float,  &
              (/x_dim_id,y_dim_id,bin_dim_id,t_dim_id/),                &
                depocon_var_id)
        endif
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: def_var depocon: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,depocon_var_id,"long_name",var_lnames(12))
        if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,depocon_var_id,"units","kg/m2")
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att depocon: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,depocon_var_id,&
                 "missing_value", -9999.0)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att depocon: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,depocon_var_id,"_FillValue",-9999.0_op)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att depocon: ',nf90_strerror(nSTAT)  
  
        ! DEPOSIT thickness
        if(VERB.gt.1)write(global_info,*)"     depothick",var_lnames(30)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"depothick",nf90_double, &
              (/x_dim_id,y_dim_id,t_dim_id/),                &
                depothick_var_id)
        else
          nSTAT = nf90_def_var(ncid,"depothick",nf90_float,  &
              (/x_dim_id,y_dim_id,t_dim_id/),                &
                depothick_var_id)
        endif
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: def_var depothick: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,depothick_var_id,"long_name",var_lnames(30))
        if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,depothick_var_id,"units","mm")
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att depothick: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,depothick_var_id,&
                 "missing_value", -9999.0)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att depothick: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,depothick_var_id,"_FillValue", -9999.0_op)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att depothick: ',nf90_strerror(nSTAT)
 
           ! Arrival time of deposit
        if(VERB.gt.1)write(global_info,*)"     depotime",var_lnames(21)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"depotime",nf90_double, &
              (/x_dim_id,y_dim_id/),                &
                depotime_var_id)
        else
          nSTAT = nf90_def_var(ncid,"depotime",nf90_float,  &
              (/x_dim_id,y_dim_id/),                &
                depotime_var_id)
        endif
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: def_var depotime: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,depotime_var_id,"long_name",var_lnames(21))
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att: depotime',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,depotime_var_id,"units","hr")
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att depotime: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,depotime_var_id,&
                   "missing_value", -9999.0)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att depotime: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,depotime_var_id,"_FillValue",-9999.0_op)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att depotime: ',nf90_strerror(nSTAT)
 
           ! Arrival time of airborne ash cloud
        if(VERB.gt.1)write(global_info,*)"     ash_arrival_time",var_lnames(31)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"ash_arrival_time",nf90_double, &
              (/x_dim_id,y_dim_id/),                &
                ashcloudtime_var_id)
        else
          nSTAT = nf90_def_var(ncid,"ash_arrival_time",nf90_float,  &
              (/x_dim_id,y_dim_id/),                &
                ashcloudtime_var_id)
        endif
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: def_var ash_arrival_time: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashcloudtime_var_id,"long_name",var_lnames(31))
        if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashcloudtime_var_id,"units","hr")
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att ash_arrival_time: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashcloudtime_var_id,&
                 "missing_value", -9999.0)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att ash_arrival_time: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashcloudtime_var_id,"_FillValue",-9999.0_op)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att ash_arrival_time: ',nf90_strerror(nSTAT)
  
           ! 2d ash concentration (MaxConcentration)
        if(VERB.gt.1)write(global_info,*)"     ashcon_max",var_lnames(32)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"ashcon_max",nf90_double, &
              (/x_dim_id,y_dim_id,t_dim_id/),                &
                ashconMax_var_id)
        else
          nSTAT = nf90_def_var(ncid,"ashcon_max",nf90_float,  &
              (/x_dim_id,y_dim_id,t_dim_id/),                &
                ashconMax_var_id)
        endif
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: def_var ashcon_max: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashconMax_var_id,"long_name",var_lnames(32))
        if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashconMax_var_id,"units","mg/m3")
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att ashcon_max: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashconMax_var_id,&
                 "missing_value", -9999.0)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att ashcon_max: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashconMax_var_id,"_FillValue",-9999.0_op)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att ashcon_max: ',nf90_strerror(nSTAT)
  
           ! 2d ash cloud height (MaxHeight)
        if(VERB.gt.1)write(global_info,*)"     cloud_height",var_lnames(33)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"cloud_height",nf90_double, &
              (/x_dim_id,y_dim_id,t_dim_id/),                &
                ashheight_var_id)
        else
          nSTAT = nf90_def_var(ncid,"cloud_height",nf90_float,  &
              (/x_dim_id,y_dim_id,t_dim_id/),                &
                ashheight_var_id)
        endif
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: def_var cloud_height: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashheight_var_id,"long_name",var_lnames(33))
        if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashheight_var_id,"units","km")
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att cloud_height: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashheight_var_id,&
                 "missing_value", -9999.0)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att cloud_height: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashheight_var_id,"_FillValue",-9999.0_op)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att cloud_height: ',nf90_strerror(nSTAT)
  
           ! 2d ash cloud load (CloudLoad)
        if(VERB.gt.1)write(global_info,*)"     cloud_load",var_lnames(34)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"cloud_load",nf90_double, &
              (/x_dim_id,y_dim_id,t_dim_id/),                &
                ashload_var_id)
        else
          nSTAT = nf90_def_var(ncid,"cloud_load",nf90_float,  &
              (/x_dim_id,y_dim_id,t_dim_id/),                &
                ashload_var_id)
        endif
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: def_var cloud_load: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashload_var_id,"long_name",var_lnames(34))
        if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashload_var_id,"units","T/km2")
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att cloud_load: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashload_var_id,&
                 "missing_value", -9999.0)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att cloud_load: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashload_var_id,"_FillValue",-9999.0_op)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att cloud_load: ',nf90_strerror(nSTAT)
  
           ! 3d radar reflectivity
        if(VERB.gt.1)write(global_info,*)"     radar_reflectivity",var_lnames(35)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"radar_reflectivity",nf90_double, &
              (/x_dim_id,y_dim_id,z_dim_id,t_dim_id/),                &
                radrefl_var_id)
        else
          nSTAT = nf90_def_var(ncid,"radar_reflectivity",nf90_float,  &
              (/x_dim_id,y_dim_id,z_dim_id,t_dim_id/),                &
                radrefl_var_id)
        endif
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: def_var radar_reflectivity: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,radrefl_var_id,"long_name",var_lnames(35))
        if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,radrefl_var_id,"units","dbZ")
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att radar_reflectivity: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,radrefl_var_id,&
                 "missing_value", -9999.0)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att radar_reflectivity: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,radrefl_var_id,"_FillValue",-9999.0_op)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att radar_reflectivity: ',nf90_strerror(nSTAT)
  
           ! 2d ash cloud bottom 
        if(VERB.gt.1)write(global_info,*)"     cloud_bottom",var_lnames(36)
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"cloud_bottom",nf90_double, &
              (/x_dim_id,y_dim_id,t_dim_id/),                &
                ashcloudBot_var_id)
        else
          nSTAT = nf90_def_var(ncid,"cloud_bottom",nf90_float,  &
              (/x_dim_id,y_dim_id,t_dim_id/),                &
                ashcloudBot_var_id)
        endif
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: def_var cloud_bottom: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashcloudBot_var_id,"long_name",var_lnames(36))
        if(nSTAT.ne.0)write(global_log ,*)'ERROR: put_att: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashcloudBot_var_id,"units","km")
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att cloud_bottom: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashcloudBot_var_id,&
                 "missing_value", -9999.0)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att cloud_bottom: ',nf90_strerror(nSTAT)
        nSTAT = nf90_put_att(ncid,ashcloudBot_var_id,"_FillValue",-9999.0_op)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_att cloud_bottom: ',nf90_strerror(nSTAT)

      endif ! USE_OUTPROD_VARS
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(USE_OPTMOD_VARS)then
        if(VERB.gt.1)write(global_info,*)"     USE_OPTMOD_VARS"
        ! Define User-specified 2-d static variables
        if(nvar_User2d_static_XY.gt.0)then
          do ivar=1,nvar_User2d_static_XY
            if(op.eq.8)then
              nSTAT = nf90_def_var(ncid,var_User2d_static_XY_name(ivar),nf90_double,  &
                  (/x_dim_id,y_dim_id/), temp1_2d_var_id)
            else
              nSTAT = nf90_def_var(ncid,var_User2d_static_XY_name(ivar),nf90_float,  &
                  (/x_dim_id,y_dim_id/), temp1_2d_var_id)
            endif
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: def_var XYs: ',nf90_strerror(nSTAT)
            nSTAT = nf90_put_att(ncid,temp1_2d_var_id,"long_name",&
                                 var_User2d_static_XY_lname(ivar))
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: put_att XYs: ',nf90_strerror(nSTAT)
            nSTAT = nf90_put_att(ncid,temp1_2d_var_id,"units",&
                                 var_User2d_static_XY_unit(ivar))
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: put_att XYs: ',nf90_strerror(nSTAT)
            nSTAT = nf90_put_att(ncid,temp1_2d_var_id,"missing_value",&
                                 var_User2d_static_XY_MissVal(ivar))
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: put_att XYs: ',nf90_strerror(nSTAT)
            nSTAT = nf90_put_att(ncid,temp1_2d_var_id,"_FillValue",&
                                 var_User2d_static_XY_FillVal(ivar))
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: put_att XYs: ',nf90_strerror(nSTAT)
          enddo
        endif

        ! Define User-specified 2-d transient variables
        if(nvar_User2d_XY.gt.0)then
          do ivar=1,nvar_User2d_XY
            if(op.eq.8)then
              nSTAT = nf90_def_var(ncid,var_User2d_XY_name(ivar),nf90_double,  &
                  (/x_dim_id,y_dim_id,t_dim_id/), temp1_2d_var_id)
            else
              nSTAT = nf90_def_var(ncid,var_User2d_XY_name(ivar),nf90_float,  &
                  (/x_dim_id,y_dim_id,t_dim_id/), temp1_2d_var_id)
            endif
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: def_var XY: ',nf90_strerror(nSTAT)
            nSTAT = nf90_put_att(ncid,temp1_2d_var_id,"long_name",&
                                 var_User2d_XY_lname(ivar))
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: put_att XY: ',nf90_strerror(nSTAT)
            nSTAT = nf90_put_att(ncid,temp1_2d_var_id,'units',&
                                 var_User2d_XY_unit(ivar))
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: put_att XY: ',nf90_strerror(nSTAT)
            nSTAT = nf90_put_att(ncid,temp1_2d_var_id,"missing_value",&
                                 var_User2d_XY_MissVal(ivar))
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: put_att XY: ',nf90_strerror(nSTAT)
            nSTAT = nf90_put_att(ncid,temp1_2d_var_id,"_FillValue",&
                                 var_User2d_XY_FillVal(ivar))
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: put_att XY: ',nf90_strerror(nSTAT)
          enddo
        endif

        !! Define User-specified 3-d transient variables in x,y,gs
        if(nvar_User3d_XYGs.gt.0)then
          do ivar=1,nvar_User3d_XYGs
            if(op.eq.8)then
              nSTAT = nf90_def_var(ncid,var_User3d_XYGs_name(ivar),nf90_double,  &
                  (/x_dim_id,y_dim_id,bin_dim_id,t_dim_id/), temp1_3d_var_id)
            else
              nSTAT = nf90_def_var(ncid,var_User3d_XYGs_name(ivar),nf90_float,  &
                  (/x_dim_id,y_dim_id,bin_dim_id,t_dim_id/), temp1_3d_var_id)
            endif
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: def_var XYGs: ',nf90_strerror(nSTAT)
            nSTAT = nf90_put_att(ncid,temp1_3d_var_id,"long_name",&
                                 var_User3d_XYGs_lname(ivar))
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: put_att XYGs: ',nf90_strerror(nSTAT)
            nSTAT = nf90_put_att(ncid,temp1_3d_var_id,'units',&
                                 var_User3d_XYGs_unit(ivar))
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: put_att XYGs: ',nf90_strerror(nSTAT)
            nSTAT = nf90_put_att(ncid,temp1_3d_var_id,"missing_value",&
                                 var_User3d_XYGs_MissVal(ivar))
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: put_att XYGs: ',nf90_strerror(nSTAT)
            nSTAT = nf90_put_att(ncid,temp1_3d_var_id,"_FillValue",&
                                 var_User3d_XYGs_FillVal(ivar))
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: put_att XYGs: ',nf90_strerror(nSTAT)
          enddo
        endif

        !! Define User-specified 3-d transient variables in x,y,z
        if(nvar_User3d_XYZ.gt.0)then
          do ivar=1,nvar_User3d_XYZ
            if(op.eq.8)then
              nSTAT = nf90_def_var(ncid,var_User3d_XYZ_name(ivar),nf90_double,  &
                  (/x_dim_id,y_dim_id,z_dim_id,t_dim_id/), temp1_3d_var_id)
            else
              nSTAT = nf90_def_var(ncid,var_User3d_XYZ_name(ivar),nf90_float,  &
                  (/x_dim_id,y_dim_id,z_dim_id,t_dim_id/), temp1_3d_var_id)
            endif
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: def_var XYZ: ',nf90_strerror(nSTAT)
            nSTAT = nf90_put_att(ncid,temp1_3d_var_id,"long_name",&
                                 var_User3d_XYZ_lname(ivar))
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: put_att XYZ: ',nf90_strerror(nSTAT)
            nSTAT = nf90_put_att(ncid,temp1_3d_var_id,'units',&
                                 var_User3d_XYZ_unit(ivar))
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: put_att XYZ: ',nf90_strerror(nSTAT)
            nSTAT = nf90_put_att(ncid,temp1_3d_var_id,"missing_value",&
                                 var_User3d_XYZ_MissVal(ivar))
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: put_att XYZ: ',nf90_strerror(nSTAT)
            nSTAT = nf90_put_att(ncid,temp1_3d_var_id,"_FillValue",&
                                 var_User3d_XYZ_FillVal(ivar))
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: put_att XYZ: ',nf90_strerror(nSTAT)
          enddo
        endif

        ! Define User-specified 4-d transient variables in x,y,z,gs
        if(nvar_User4d_XYZGs.gt.0)then
          do ivar=1,nvar_User4d_XYZGs
            if(op.eq.8)then
              nSTAT = nf90_def_var(ncid,var_User4d_XYZGs_name(ivar),nf90_double,  &
                  (/x_dim_id,y_dim_id,z_dim_id,bin_dim_id,t_dim_id/), temp1_4d_var_id)
            else
              nSTAT = nf90_def_var(ncid,var_User4d_XYZGs_name(ivar),nf90_float, &
                  (/x_dim_id,y_dim_id,z_dim_id,bin_dim_id,t_dim_id/), temp1_4d_var_id)
            endif
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: def_var XYZGs: ',nf90_strerror(nSTAT)
            nSTAT = nf90_put_att(ncid,temp1_4d_var_id,"long_name",&
                                 var_User4d_XYZGs_lname(ivar))
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: put_att XYZGs: ',nf90_strerror(nSTAT)
            nSTAT = nf90_put_att(ncid,temp1_4d_var_id,'units',&
                                 var_User4d_XYZGs_unit(ivar))
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: put_att XYZGs: ',nf90_strerror(nSTAT)
            nSTAT = nf90_put_att(ncid,temp1_4d_var_id,"missing_value",&
                                 var_User4d_XYZGs_MissVal(ivar))
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: put_att XYZGs: ',nf90_strerror(nSTAT)
            nSTAT = nf90_put_att(ncid,temp1_4d_var_id,"_FillValue",&
                                 var_User4d_XYZGs_FillVal(ivar))
            if(nSTAT.ne.0) &
                write(global_log ,*)'ERROR: put_att XYZGs: ',nf90_strerror(nSTAT)
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
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: nf90_enddef: ',nf90_strerror(nSTAT)
      if(VERB.gt.1)write(*,*)"Leaving define mode"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Fill variables with initial values
        ! Time
      if(VERB.gt.1)write(global_info,*)"     Fill time"
      dumscal_out = real(time,kind=op)
      nSTAT=nf90_put_var(ncid,t_var_id,dumscal_out,(/1/))
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: put_var time: ',nf90_strerror(nSTAT)
        ! Z
      if(VERB.gt.1)write(global_info,*)"     Fill Z"
      !z_out(:) = z(1:nzmax)
      allocate(dum1d_out(nzmax))
      dum1d_out(:) = real(z_cc_pd(1:nzmax),kind=op)
      nSTAT=nf90_put_var(ncid,z_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: put_var z: ',nf90_strerror(nSTAT)
      deallocate(dum1d_out)

        ! Y
      if(VERB.gt.1)write(global_info,*)"     Fill Y"
      !y_out(:) = y(1:nymax)
      allocate(dum1d_out(nymax))
      if (IsLatLon) then
        dum1d_out(1:nymax) = real(lat_cc_pd(1:nymax),kind=op)
      else
        dum1d_out(1:nymax) = real(y_cc_pd(1:nymax),kind=op)
      endif
      nSTAT=nf90_put_var(ncid,y_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: put_var y: ',nf90_strerror(nSTAT)
      deallocate(dum1d_out)
        ! X
      if(VERB.gt.1)write(global_info,*)"     Fill X"
      !x_out(:) = x(1:nxmax)
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
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: put_var x: ',nf90_strerror(nSTAT)
      deallocate(dum1d_out)
        ! BN (Grain size bin ID)
      if(VERB.gt.1)write(global_info,*)"     Fill GS"
      allocate(dum1dint_out(nsmax))
      do i=1,nsmax
        dum1dint_out(i) = i
      enddo
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: put_var gs: ',nf90_strerror(nSTAT)
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
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: put_var er: ',nf90_strerror(nSTAT)
      deallocate(dum1dint_out)

        ! WS
      if(VERB.gt.1)write(global_info,*)"     Fill WS"
      allocate(dum1dint_out(MR_iwindfiles))
      ! This is variable associated with the dimension for windfile ID
      ! This only contains the index starting with 1
      do i=1,MR_iwindfiles
        dum1dint_out(i) = i
      enddo
      nSTAT=nf90_put_var(ncid,wf_var_id,dum1dint_out,(/1/))
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: put_var wd: ',nf90_strerror(nSTAT)
      deallocate(dum1dint_out)
      !   End of filling dimension variables
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !   Now fill a few other variables that are a function of BN
         ! Species class ID
      if(VERB.gt.1)write(global_info,*)"     Fill Species class ID"
      allocate(dum1dint_out(ns_extra))
      dum1dint_out(1:nsmax) = SpeciesID(1:nsmax)
      nSTAT=nf90_put_var(ncid,spec_var_id,dum1dint_out,(/1/))
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: put_var sp: ',nf90_strerror(nSTAT)
      deallocate(dum1dint_out)
         ! Species sub-class ID
      if(VERB.gt.1)write(global_info,*)"     Fill Species sub-class ID"
      allocate(dum1dint_out(ns_extra))
      dum1dint_out(1:nsmax) = SpeciesSubID(1:nsmax)
      nSTAT=nf90_put_var(ncid,subspec_var_id,dum1dint_out,(/1/))
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: put_var sp: ',nf90_strerror(nSTAT)
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
         ! gs_massfrac (Mass fraction of grain size)
      if(VERB.gt.1)write(global_info,*)"     Fill GS MassFrac"
      dum1d_out = 0.0_op
      dum1d_out(1:n_gs_max) = real(Tephra_bin_mass(1:n_gs_max),kind=op)
      nSTAT=nf90_put_var(ncid,gsmf_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: put_var gs_massfrac: ',nf90_strerror(nSTAT)
       ! gs_dens (Density of grain)
      if(useCalcFallVel)then
        dum1d_out(1:n_gs_max) = real(Tephra_rho_m(1:n_gs_max),kind=op)
      else
        dum1d_out = 0.0_op
      endif
      nSTAT=nf90_put_var(ncid,gsdens_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: put_var gs_dens: ',nf90_strerror(nSTAT)
      deallocate(dum1d_out)

      !   Now fill a few other variables that are a function of ER
        ! er_stime (Start time of eruption)
      allocate(dum1d_out(neruptions))
      dum1d_out = real(e_StartTime + SimStartHour,kind=op)
      nSTAT=nf90_put_var(ncid,er_stime_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: put_var er_stime: ',nf90_strerror(nSTAT)
        ! er_duration (Duration of eruption)
      dum1d_out = real(e_Duration,kind=op)
      nSTAT=nf90_put_var(ncid,er_duration_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: put_var er_stime: ',nf90_strerror(nSTAT)
        ! er_plumeheight (Plume height of eruption)
      dum1d_out = real(PlumeHeight,kind=op)
      nSTAT=nf90_put_var(ncid,er_plumeheight_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: put_var er_stime: ',nf90_strerror(nSTAT)
        ! er_volume (Volume of eruption)
      dum1d_out = real(e_Volume,kind=op)
      nSTAT=nf90_put_var(ncid,er_volume_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: put_var er_stime: ',nf90_strerror(nSTAT)
      deallocate(dum1d_out)

         ! Now fill the other (non-time-dependent) variables
         ! wf_name (Name of windfile)
      do i=1,MR_iwindfiles
        write(lllinebuffer,'(130a)')MR_windfiles(i)
        do j=1,130
          nSTAT=nf90_put_var(ncid,wf_name_var_id,lllinebuffer(j:j),(/j,i/))
          if(nSTAT.ne.0) &
            write(global_log ,*)'ERROR: put_var wf_name: ',nf90_strerror(nSTAT)
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
      if(nSTAT.ne.0) &
        write(global_info,*)'ERROR: put_var cell area: ',nf90_strerror(nSTAT)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(USE_WIND_VARS)then
          ! Vz
        if(VERB.gt.1)write(global_info,*)"     Fill Vz"
        dum3d_out(1:nxmax,1:nymax,1:nzmax) = real(vz_pd(1:nxmax,1:nymax,1:nzmax),kind=op)
        !dum3d_out(:,:,:) = vf(1:nxmax,1:nymax,1:nzmax,1)
        nSTAT=nf90_put_var(ncid,vz_var_id,dum3d_out,(/1,1,1,1/))
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_var Vz: ',nf90_strerror(nSTAT)
          ! Vy
        if(VERB.gt.1)write(global_info,*)"     Fill Vy"
        dum3d_out(1:nxmax,1:nymax,1:nzmax) = real(vy_pd(1:nxmax,1:nymax,1:nzmax),kind=op)
        nSTAT=nf90_put_var(ncid,vy_var_id,dum3d_out,(/1,1,1,1/))
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_var Vy: ',nf90_strerror(nSTAT)
          ! Vx
        if(VERB.gt.1)write(global_info,*)"     Fill Vx"
        dum3d_out(1:nxmax,1:nymax,1:nzmax) = real(vx_pd(1:nxmax,1:nymax,1:nzmax),kind=op)
        nSTAT=nf90_put_var(ncid,vx_var_id,dum3d_out,(/1,1,1,1/))
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_var Vx: ',nf90_strerror(nSTAT)
          ! Vf
        ashcon(1:nxmax,1:nymax,1:nzmax,1:nsmax) = &
             real(vf_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax),kind=op)
        nSTAT=nf90_put_var(ncid,vf_var_id,ashcon,(/1,1,1,1,1/))
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_var Vf: ',nf90_strerror(nSTAT)
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
          ashcon(1:nxmax,1:nymax,1:nzmax,i) = real(concen_pd(1:nxmax,1:nymax,1:nzmax,i,ts1),kind=op)
        enddo
        nSTAT=nf90_put_var(ncid,ashcon_var_id,ashcon,(/1,1,1,1,1/))
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_var ashcon: ',nf90_strerror(nSTAT)

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
            ! Here's the previous calculation for depth in meters
          !depocon(:,:,i)=concen_pd(1:nxmax,1:nymax,0,i,ts1)*dz*1.0e-6_op/DepositDensity
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
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_var depocon: ',nf90_strerror(nSTAT)
  
        if(VERB.gt.1)write(global_info,*)"     Calling dbZCalculator"
        call dbZCalculator            ! get radar reflectivity
  
        ! depothick
        if(VERB.gt.1)write(global_info,*)"     Fill depothick"
        dum2d_out(:,:) = 0.0_op
        do i=1,nxmax
          do j=1,nymax
            if(DepositThickness(i,j).ge.0.0_ip)&
                dum2d_out(i,j) = real(DepositThickness(i,j),kind=op)
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,depothick_var_id,dum2d_out,(/1,1,1/))
        if(nSTAT.ne.0) &
           write(global_info,*)'ERROR: put_var depothick: ',nf90_strerror(nSTAT)
 
          ! depotime
        if(VERB.gt.1)write(global_info,*)"     Fill depotime"
        dum2d_out(:,:) = -9999.0_op
        do i=1,nxmax
          do j=1,nymax
            if(DepArrivalTime(i,j).ge.0.0_ip)&
                dum2d_out(i,j)=real(DepArrivalTime(i,j),kind=op)
          enddo
        enddo

        nSTAT=nf90_put_var(ncid,depotime_var_id,dum2d_out,(/1,1/))
        if(nSTAT.ne.0) &
          write(global_info,*)'ERROR: put_var depotime: ',nf90_strerror(nSTAT)

          ! ashtime
        if(VERB.gt.1)write(global_info,*)"     Fill ashtime"
        dum2d_out(:,:) = -9999.0_op
        do i=1,nxmax
          do j=1,nymax
            if(CloudArrivalTime(i,j).ge.0.0)&
                dum2d_out(i,j)=real(CloudArrivalTime(i,j),kind=op)
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,ashcloudtime_var_id,dum2d_out,(/1,1/))
        if(nSTAT.ne.0) &
          write(global_info,*)'ERROR: put_var ashcloudtime_: ',nf90_strerror(nSTAT)
  
        ! ashconMax
        if(VERB.gt.1)write(global_info,*)"     Fill ashconMax"
        dum2d_out(:,:) = 0.0_op
        do i=1,nxmax
          do j=1,nymax
            !if(MaxConcentration(i,j).gt.0.0_ip)write(*,*)MaxConcentration(i,j)
            !if(MaxConcentration(i,j).ge.100.0_ip)then
            if(MaxConcentration(i,j).ge.EPS_TINY)then
              dumscal_out=real(MaxConcentration(i,j)/1.0e3_ip,kind=op)
              dum2d_out(i,j)=dumscal_out
            endif
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,ashconMax_var_id,dum2d_out,(/1,1,1/))
        if(nSTAT.ne.0) &
          write(global_info,*)'ERROR: put_var ashconMax: ',nf90_strerror(nSTAT)
  
        ! ash cloud_height (top)
        if(VERB.gt.1)write(global_info,*)"     Fill ash_height"
        dum2d_out(:,:) = -9999.0_op
        do i=1,nxmax
          do j=1,nymax
            if(MaxHeight(i,j).gt.0.0_ip)then
              dumscal_out=real(MaxHeight(i,j),kind=op)
              dum2d_out(i,j) = dumscal_out
            endif
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,ashheight_var_id,dum2d_out,(/1,1,1/))
        if(nSTAT.ne.0) &
          write(global_info,*)'ERROR: put_var ashheight: ',nf90_strerror(nSTAT)

        ! ash-load
        if(VERB.gt.1)write(global_info,*)"     Fill ash-load"
        dum2d_out(:,:) = 0.0_op
        do i=1,nxmax
          do j=1,nymax
            if(CloudLoad(i,j).ge.CLOUDLOAD_THRESH)then
            !if(CloudLoad(i,j).ge.0.0_ip)then
              dumscal_out=real(CloudLoad(i,j),kind=op)
              dum2d_out(i,j)=dumscal_out
            endif
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,ashload_var_id,dum2d_out,(/1,1,1/))
        if(nSTAT.ne.0) &
          write(global_info,*)'ERROR: put_var ashload: ',nf90_strerror(nSTAT)

        ! radar-reflectivity
        if(VERB.gt.1)write(global_info,*)"     Fill dbZ"
        dum3d_out(:,:,:) = 0.0_op
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
        if(nSTAT.ne.0) &
          write(global_info,*)'ERROR: put_var reflectivity: ',nf90_strerror(nSTAT)

        ! ash cloud_bottom
        if(VERB.gt.1)write(global_info,*)"     Fill ash height min"
        dum2d_out(:,:) = -9999.0_op
        do i=1,nxmax
          do j=1,nymax
            if(MinHeight(i,j).ge.0.0_ip)then
              dumscal_out=real(MinHeight(i,j),kind=op)
              dum2d_out(i,j) = dumscal_out
            endif
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,ashcloudBot_var_id,dum2d_out,(/1,1,1/))
        if(nSTAT.ne.0) &
          write(global_info,*)'ERROR: put_var cloud bottom: ',nf90_strerror(nSTAT)

      endif ! USE_OUTPROD_VARS
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
            if(nSTAT.ne.0) &
               write(global_info,*)'ERROR: put_inq_varid static_XY: ',ivar,nf90_strerror(nSTAT)
            nSTAT = nf90_put_var(ncid,temp1_2d_var_id,dum2d_out,(/1,1/))
            if(nSTAT.ne.0) &
               write(global_info,*)'ERROR: put_var static_XY: ',ivar,nf90_strerror(nSTAT)
          enddo
        endif

        ! Fill User-specified 2-d transient variables
        if(nvar_User2d_XY.gt.0)then
          do ivar=1,nvar_User2d_XY
            dum2d_out(:,:) = real(var_User2d_XY(:,:,ivar),kind=op)
            nSTAT = nf90_inq_varid(ncid,var_User2d_XY_name(ivar),temp1_2d_var_id)
            if(nSTAT.ne.0) &
               write(global_info,*)'ERROR: inq_var XY: ',ivar,nf90_strerror(nSTAT)
            nSTAT = nf90_put_var(ncid,temp1_2d_var_id,dum2d_out,(/1,1,1/))
            if(nSTAT.ne.0) &
               write(global_info,*)'ERROR: put_var XY: ',ivar,nf90_strerror(nSTAT)
          enddo
        endif

        ! Fill User-specified 3-d transient variables in x,y,gs
        if(nvar_User3d_XYGs.gt.0)then
          do ivar=1,nvar_User3d_XYGs
            depocon(:,:,:) = real(var_User3d_XYGs(:,:,:,ivar),kind=op)
            nSTAT = nf90_inq_varid(ncid,var_User3d_XYGs_name(ivar),temp1_3d_var_id)
            if(nSTAT.ne.0) &
               write(global_info,*)'ERROR: inq_var XYGs: ',ivar,nf90_strerror(nSTAT)
            nSTAT = nf90_put_var(ncid,temp1_3d_var_id,depocon,(/1,1,1,1/))
            if(nSTAT.ne.0) &
               write(global_info,*)'ERROR: put_var XYGs: ',ivar,nf90_strerror(nSTAT)
          enddo
        endif

        ! Fill User-specified 3-d transient variables in x,y,z
        if(nvar_User3d_XYZ.gt.0)then
          do ivar=1,nvar_User3d_XYZ
            dum3d_out(:,:,:) = real(var_User3d_XYZ(:,:,:,ivar),kind=op)
            nSTAT = nf90_inq_varid(ncid,var_User3d_XYZ_name(ivar),temp1_3d_var_id)
            if(nSTAT.ne.0) &
               write(global_info,*)'ERROR: inq_var XYZ: ',ivar,nf90_strerror(nSTAT)
            nSTAT = nf90_put_var(ncid,temp1_3d_var_id,dum3d_out,(/1,1,1,1/))
            if(nSTAT.ne.0) &
               write(global_info,*)'ERROR: put_var XYZ: ',ivar,nf90_strerror(nSTAT)
          enddo
        endif

        ! Fill User-specified 4-d transient variables in x,y,z,gs
        if(nvar_User4d_XYZGs.gt.0)then
          do ivar=1,nvar_User4d_XYZGs
            ashcon(:,:,:,:) = real(var_User4d_XYZGs(1:nxmax,1:nymax,1:nzmax,1:nsmax,ivar),kind=op)
            nSTAT = nf90_inq_varid(ncid,var_User4d_XYZGs_name(ivar),temp1_4d_var_id)
            if(nSTAT.ne.0) &
               write(global_info,*)'ERROR: inq_var XYZGs: ',ivar,nf90_strerror(nSTAT)
            nSTAT = nf90_put_var(ncid,temp1_4d_var_id,ashcon,(/1,1,1,1,1/))
            if(nSTAT.ne.0) &
               write(global_info,*)'ERROR: put_var XYZGs: ',ivar,nf90_strerror(nSTAT)
          enddo
        endif

      endif  ! USE_OPTMOD_VARS
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Close file
      nSTAT = nf90_close(ncid)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: close outfile: ',                &
                              nf90_strerror(nSTAT)

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

      use precis_param

      use io_units

      use global_param,  only : &
         EPS_SMALL,EPS_TINY,VERB,KM2_2_M2

      use io_data,       only : &
         io,nvar_User2d_XY,nvar_User3d_XYGs,nvar_User3d_XYZ,nvar_User4d_XYZGs,&
         outfile

      use Output_Vars,   only : &
         var_User2d_XY_name,var_User2d_XY,var_User3d_XYGs_name,var_User3d_XYGs,&
         var_User3d_XYZ_name,var_User3d_XYZ,var_User4d_XYZGs_name,var_User4d_XYZGs,&
         CLOUDLOAD_THRESH,DBZ_THRESH,USE_OPTMOD_VARS,USE_RESTART_VARS,&
         USE_OUTPROD_VARS,USE_WIND_VARS,DepositThickness,DepArrivalTime,CloudArrivalTime,&
         MaxConcentration,MaxHeight,CloudLoad,dbZ,MinHeight,&
           dbZCalculator

      use Tephra,        only : &
         n_gs_max

      use mesh,          only : &
         nxmax,nymax,nzmax,ts1,dz_vec_pd,nsmax

      use solution,      only : &
         vx_pd,vy_pd,vz_pd,vf_pd,concen_pd,DepositGranularity

      use time_data,     only : &
         time

      use netcdf

      implicit none

      integer :: i,j,k,n

      integer :: nSTAT
      integer :: ncid
      integer :: ivar
      integer :: t_var_id          = 0
      integer :: vx_var_id         = 0
      integer :: vy_var_id         = 0
      integer :: vz_var_id         = 0
      integer :: vf_var_id         = 0
      integer :: ashcon_var_id     = 0
      integer :: depocon_var_id    = 0

      integer :: ashconMax_var_id      = 0 ! Max Ash concentration in column
      integer :: ashheight_var_id      = 0 ! Height of top of ash cloud
      integer :: ashload_var_id        = 0 ! Vert. integrated load of ash cloud
      integer :: radrefl_var_id        = 0 ! Radar reflectivity in dbZ
      integer :: depothick_var_id     = 0 ! Total deposit (summed over gs)
      integer :: depotime_var_id       = 0 ! Deposit arrival time
      integer :: ashcloudtime_var_id   = 0 ! Cloud arrival time
      integer :: ashcloudBot_var_id    = 0 ! Height of bottom of ash cloud

      integer :: temp1_2d_var_id !,temp2_2d_var_id,temp3_2d_var_id,temp4_2d_var_id
      integer :: temp1_3d_var_id !,temp2_3d_var_id,temp3_3d_var_id,temp4_3d_var_id
      integer :: temp1_4d_var_id !,temp2_4d_var_id!,temp3_4d_var_id!,temp4_4d_var_id
      integer :: ns_extra

      real(kind=op)                                 :: dumscal_out
      real(kind=op), dimension(:,:)    ,allocatable :: dum2d_out
      real(kind=op), dimension(:,:,:)  ,allocatable :: dum3d_out

      real(kind=op), dimension(:,:,:,:),allocatable :: ashcon
      real(kind=op), dimension(:,:,:)  ,allocatable :: depocon

      if(VERB.gt.2)write(global_info,*)"Allocating output vars"
      allocate(dum2d_out(nxmax,nymax))
      allocate(dum3d_out(nxmax,nymax,nzmax))

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
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: open outfile: ',nf90_strerror(nSTAT)

      ! Get variable ids
      nSTAT = nf90_inq_varid(ncid,"t",t_var_id)
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: inq_varid t: ',nf90_strerror(nSTAT)

      if(VERB.gt.2)write(global_info,*)"Got var IDs, now writing data"
      ! Write data
      ! Time
      if(VERB.gt.2)write(global_info,*)"  Writing Time"
      dumscal_out = real(time,kind=op)
      nSTAT=nf90_put_var(ncid,t_var_id,dumscal_out,(/io/))
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: put_var t: ',nf90_strerror(nSTAT)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(USE_WIND_VARS)then
          ! Vz
        nSTAT = nf90_inq_varid(ncid,"vz",vz_var_id)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: inq_varid vz: ',nf90_strerror(nSTAT)
        if(VERB.gt.2)write(global_info,*)"  Writing Vz"
        dum3d_out(1:nxmax,1:nymax,1:nzmax) = real(vz_pd(1:nxmax,1:nymax,1:nzmax),kind=op)
        nSTAT=nf90_put_var(ncid,vz_var_id,dum3d_out,(/1,1,1,io/))
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_var Vz: ',nf90_strerror(nSTAT)
          ! Vy
        nSTAT = nf90_inq_varid(ncid,"vy",vy_var_id)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: inq_varid vy: ',nf90_strerror(nSTAT)
        if(VERB.gt.2)write(global_info,*)"  Writing Vy"
        dum3d_out(1:nxmax,1:nymax,1:nzmax) = real(vy_pd(1:nxmax,1:nymax,1:nzmax),kind=op)
        nSTAT=nf90_put_var(ncid,vy_var_id,dum3d_out,(/1,1,1,io/))
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_var Vy: ',nf90_strerror(nSTAT)
          ! Vx
        nSTAT = nf90_inq_varid(ncid,"vx",vx_var_id)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: inq_varid vx: ',nf90_strerror(nSTAT)
        if(VERB.gt.2)write(global_info,*)"  Writing Vx"
        dum3d_out(1:nxmax,1:nymax,1:nzmax) = real(vx_pd(1:nxmax,1:nymax,1:nzmax),kind=op)
        nSTAT=nf90_put_var(ncid,vx_var_id,dum3d_out,(/1,1,1,io/))
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_var Vz: ',nf90_strerror(nSTAT)
          ! Vf
        nSTAT = nf90_inq_varid(ncid,"vf",vf_var_id)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: inq_varid vf: ',nf90_strerror(nSTAT)
        ashcon(1:nxmax,1:nymax,1:nzmax,1:nsmax) = real(vf_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax),kind=op)
        nSTAT=nf90_put_var(ncid,vf_var_id,ashcon,(/1,1,1,1,io/))
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
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: inq_varid ashcon: ',nf90_strerror(nSTAT)
        if(VERB.gt.2)write(global_info,*)"  Writing ashcon"
        ashcon = 0.0_op
        ashcon(1:nxmax,1:nymax,1:nzmax,1:nsmax) = real(concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts1),kind=op)
        nSTAT=nf90_put_var(ncid,ashcon_var_id,ashcon,(/1,1,1,1,io/))
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_var ashcon: ',nf90_strerror(nSTAT)
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(USE_OUTPROD_VARS)then

          ! depocon
        if(VERB.gt.1)write(global_info,*)"     Fill depocon"
        nSTAT = nf90_inq_varid(ncid,"depocon",depocon_var_id)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: inq_varid depocon: ',nf90_strerror(nSTAT)
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
        nSTAT=nf90_put_var(ncid,depocon_var_id,depocon,(/1,1,1,io/))
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: put_var depocon: ',nf90_strerror(nSTAT)
  
        if(VERB.gt.2)write(global_info,*)"Calling dbZCalculator"
        call dbZCalculator            ! get radar reflectivity
  
        ! depothick
        nSTAT = nf90_inq_varid(ncid,"depothick",depothick_var_id)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: inq_varid depothick: ',nf90_strerror(nSTAT)
        if(VERB.gt.2)write(global_info,*)"  Writing depothick"
        dum2d_out(:,:) = 0.0_op
        do i=1,nxmax
          do j=1,nymax
            if(DepositThickness(i,j).ge.0.0_ip)&
                dum2d_out(i,j)=real(DepositThickness(i,j),kind=op)
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,depothick_var_id,dum2d_out,(/1,1,io/))
  
          ! depotime
        nSTAT = nf90_inq_varid(ncid,"depotime",depotime_var_id)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: inq_varid depotime: ',nf90_strerror(nSTAT)
        if(VERB.gt.2)write(global_info,*)"  Writing depotime"
        dum2d_out(:,:) = -9999.0_op
        do i=1,nxmax
          do j=1,nymax
            if(DepArrivalTime(i,j).ge.0.0_ip)&
                dum2d_out(i,j)=real(DepArrivalTime(i,j),kind=op)
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,depotime_var_id,dum2d_out,(/1,1/))
        if(nSTAT.ne.0) &
          write(global_info,*)'ERROR: put_var depotime: ',nf90_strerror(nSTAT)
  
          ! ashtime
        nSTAT = nf90_inq_varid(ncid,"ash_arrival_time",ashcloudtime_var_id)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: inq_varid ash_arrival_time : ',nf90_strerror(nSTAT)
        if(VERB.gt.2)write(global_info,*)"  Writing ashtime"
        dum2d_out(:,:) = -9999.0_op
        do i=1,nxmax
          do j=1,nymax
            if(CloudArrivalTime(i,j).ge.0.0_ip)&
                dum2d_out(i,j)=real(CloudArrivalTime(i,j),kind=op)
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,ashcloudtime_var_id,dum2d_out,(/1,1/))
        if(nSTAT.ne.0) &
          write(global_info,*)'ERROR: put_var ashcloudtime_: ',nf90_strerror(nSTAT)
  
        ! ashconMax
        nSTAT = nf90_inq_varid(ncid,"ashcon_max",ashconMax_var_id)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: inq_varid ashcon_max : ',nf90_strerror(nSTAT)
        if(VERB.gt.2)write(global_info,*)"  Writing ashconMax"
        dum2d_out(:,:) = 0.0_op
        do i=1,nxmax
          do j=1,nymax
            if(MaxConcentration(i,j).ge.EPS_TINY)then
              dumscal_out=real(MaxConcentration(i,j)/1.0e03_ip,kind=op)
              dum2d_out(i,j)=dumscal_out
            endif
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,ashconMax_var_id,dum2d_out,(/1,1,io/))
  
        ! ash cloud_height
        nSTAT = nf90_inq_varid(ncid,"cloud_height",ashheight_var_id)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: inq_varid cloud_height : ',nf90_strerror(nSTAT)
        if(VERB.gt.2)write(global_info,*)"  Writing ash height"
        dum2d_out(:,:) = -9999.0_op
        do i=1,nxmax
          do j=1,nymax
            if(MaxHeight(i,j).gt.0.0_ip)then
              dumscal_out=real(MaxHeight(i,j),kind=op)
              dum2d_out(i,j) = dumscal_out
            endif
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,ashheight_var_id,dum2d_out,(/1,1,io/))
  
        ! ash-load
        nSTAT = nf90_inq_varid(ncid,"cloud_load",ashload_var_id)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: inq_varid cloud_load : ',nf90_strerror(nSTAT)
        if(VERB.gt.2)write(global_info,*)"  Writing ash-load"
        dum2d_out(:,:) = -9999.0_op
        do i=1,nxmax
          do j=1,nymax
            if(CloudLoad(i,j).ge.CLOUDLOAD_THRESH)then
              dumscal_out=real(CloudLoad(i,j),kind=op)
              dum2d_out(i,j) = dumscal_out
            endif
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,ashload_var_id,dum2d_out,(/1,1,io/))
  
        ! radar reflectivity
        nSTAT = nf90_inq_varid(ncid,"radar_reflectivity",radrefl_var_id)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: inq_varid radar_reflectivity : ',nf90_strerror(nSTAT)
        if(VERB.gt.2)write(global_info,*)"  Writing radar reflectivity"
        dum3d_out(:,:,:) = 0.0_op
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
        nSTAT=nf90_put_var(ncid,radrefl_var_id,dum3d_out,(/1,1,1,io/))
  
        ! ash cloud_bottom
        nSTAT = nf90_inq_varid(ncid,"cloud_bottom",ashcloudBot_var_id)
        if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: inq_varid cloud_bottom : ',nf90_strerror(nSTAT)
        if(VERB.gt.2)write(global_info,*)"  Writing ash-height (bottom)"
        dum2d_out(:,:) = -9999.0_op
        do i=1,nxmax
          do j=1,nymax
            if(MinHeight(i,j).ge.0.0_ip)then
              dumscal_out=real(MinHeight(i,j),kind=op)
              dum2d_out(i,j) = dumscal_out
            endif
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,ashcloudBot_var_id,dum2d_out,(/1,1,io/))

      endif ! USE_OUTPROD_VARS
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(USE_OPTMOD_VARS)then
        ! User-specified 2-d static variables were filled in the
        ! create_netcdf subroutine

        ! Fill User-specified 2-d transient variables
        if(nvar_User2d_XY.gt.0)then
          do ivar=1,nvar_User2d_XY
            dum2d_out(:,:) = real(var_User2d_XY(:,:,ivar),kind=op)
            nSTAT = nf90_inq_varid(ncid,var_User2d_XY_name(ivar),temp1_2d_var_id)
            nSTAT = nf90_put_var(ncid,temp1_2d_var_id,dum2d_out,(/1,1,io/))
          enddo
        endif

        ! Fill User-specified 3-d transient variables in x,y,gs
        if(nvar_User3d_XYGs.gt.0)then
          do ivar=1,nvar_User3d_XYGs
            depocon(:,:,:) = real(var_User3d_XYGs(:,:,:,ivar),kind=op)
            nSTAT = nf90_inq_varid(ncid,var_User3d_XYGs_name(ivar),temp1_3d_var_id)
            nSTAT = nf90_put_var(ncid,temp1_3d_var_id,depocon,(/1,1,1,io/))
          enddo
        endif

        ! Fill User-specified 3-d transient variables in x,y,z
        if(nvar_User3d_XYZ.gt.0)then
          do ivar=1,nvar_User3d_XYZ
            dum3d_out(:,:,:) = real(var_User3d_XYZ(:,:,:,ivar),kind=op)
            nSTAT = nf90_inq_varid(ncid,var_User3d_XYZ_name(ivar),temp1_3d_var_id)
            nSTAT = nf90_put_var(ncid,temp1_3d_var_id,dum3d_out,(/1,1,1,io/))
          enddo
        endif

        ! Fill User-specified 4-d transient variables in x,y,z,gs
        if(nvar_User4d_XYZGs.gt.0)then
          do ivar=1,nvar_User4d_XYZGs
            ashcon(:,:,:,:) = real(var_User4d_XYZGs(1:nxmax,1:nymax,1:nzmax,1:nsmax,ivar),kind=op)
            nSTAT = nf90_inq_varid(ncid,var_User4d_XYZGs_name(ivar),temp1_4d_var_id)
            nSTAT = nf90_put_var(ncid,temp1_4d_var_id,ashcon,(/1,1,1,1,io/))
          enddo
        endif

      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Close file
      if(VERB.gt.2)write(global_info,*)"closing file"
      nSTAT = nf90_close(ncid)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: close outfile: ',                &
                              nf90_strerror(nSTAT)

      if(VERB.gt.2)write(global_info,*)"Deallocating"
      deallocate(dum2d_out,dum3d_out)
      deallocate(ashcon)
      deallocate(depocon)

      end subroutine append_to_netcdf

!##############################################################################
!
!    NC_RestartFile_ReadTimes
!
!##############################################################################

      subroutine NC_RestartFile_ReadTimes

      use precis_param

      use io_units

      use io_data,           only : &
         concenfile,init_tstep

      use time_data,         only : &
         time

      use netcdf

      implicit none

      integer :: nSTAT
      integer :: ncid

      integer :: t_dim_id
      integer :: t_var_id
      integer :: t_len
      integer :: it
      real(kind=op)                                 :: dumscal_out
      real(kind=op), allocatable, dimension(:) :: t_list

      ! Open netcdf file for writing
      nSTAT=nf90_open(concenfile,NF90_NOWRITE, ncid)
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: open concenfile: ',nf90_strerror(nSTAT)

      ! Identify dimension for time (and note size)
      nSTAT = nf90_inq_dimid(ncid,"t",t_dim_id)
      if(nSTAT.ne.0) &
        write(global_info,*)'ERROR: inq_dimid time: ', &
                            nf90_strerror(nSTAT)
      nSTAT = nf90_Inquire_Dimension(ncid,t_dim_id,len=t_len)
      if(nSTAT.ne.0) &
        write(global_info,*)'ERROR: Inquire_Dimension time: ', &
                             nf90_strerror(nSTAT)

      ! Get variable ids
      nSTAT = nf90_inq_varid(ncid,"t",t_var_id)
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: inq_varid t: ',nf90_strerror(nSTAT)

      allocate(t_list(t_len))
      nSTAT = nf90_get_var(ncid,t_var_id,t_list)

      write(*,*)"  Step :  time"
      do it = 1,t_len
        write(*,*)it,t_list(it)
      enddo

      write(global_info,*)'Enter timestep for initialization'
      read(5,*) init_tstep

      nSTAT=nf90_get_var(ncid,t_var_id,dumscal_out,(/init_tstep/))
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: put_var t: ',nf90_strerror(nSTAT)
      !write(global_info,*)dumscal_out
      time = real(dumscal_out,kind=ip)

      deallocate(t_list)

      end subroutine NC_RestartFile_ReadTimes

!##############################################################################
!
!    NC_RestartFile_LoadConcen
!
!##############################################################################

      subroutine NC_RestartFile_LoadConcen

      use precis_param

      use io_units

      use global_param,      only : &
         KM2_2_M2

      use io_data,           only : &
         concenfile,init_tstep

      use Tephra,        only : &
         n_gs_max

      use mesh,              only : &
         nxmax,nymax,nzmax,nsmax,dz_vec_pd,ts0,ts1

      use solution,          only : &
         concen_pd

      use netcdf

      implicit none

      integer :: i
      integer :: nSTAT
      integer :: ncid

      integer :: t_var_id          = 0
      integer :: vx_var_id         = 0
      integer :: vy_var_id         = 0
      integer :: vz_var_id         = 0
      integer :: ashcon_var_id     = 0
      integer :: depocon_var_id  = 0

      real(kind=op)                                 :: dumscal_out
      real(kind=op), dimension(:,:)    ,allocatable :: dum2d_out
      real(kind=op), dimension(:,:,:)  ,allocatable :: dum3d_out

      real(kind=op), dimension(:,:,:,:),allocatable :: ashcon
      real(kind=op), dimension(:,:,:)  ,allocatable :: depocon

      allocate(dum2d_out(nxmax,nymax))
      allocate(dum3d_out(nxmax,nymax,nzmax))
      allocate(ashcon(nxmax,nymax,nzmax,n_gs_max))
      allocate(depocon(nxmax,nymax,n_gs_max))

      write(global_info,*)"WARNING "
      write(global_info,*)"Input file is not currently verified "
      write(global_info,*)" with previous run."

      ! Open netcdf file for writing
      nSTAT=nf90_open(concenfile,NF90_NOWRITE, ncid)
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: open concenfile: ',nf90_strerror(nSTAT)

      ! Get variable ids
      nSTAT = nf90_inq_varid(ncid,"t",t_var_id)
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: inq_varid t: ',nf90_strerror(nSTAT)
      nSTAT = nf90_inq_varid(ncid,"vx",vx_var_id)
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: inq_varid vx: ',nf90_strerror(nSTAT)
      nSTAT = nf90_inq_varid(ncid,"vy",vy_var_id)
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: inq_varid vy: ',nf90_strerror(nSTAT)
      nSTAT = nf90_inq_varid(ncid,"vz",vz_var_id)
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: inq_varid vz: ',nf90_strerror(nSTAT)
      nSTAT = nf90_inq_varid(ncid,"ashcon",ashcon_var_id)
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: inq_varid ashcon: ',nf90_strerror(nSTAT)
      nSTAT = nf90_inq_varid(ncid,"depocon",depocon_var_id)
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: inq_varid depocon: ',nf90_strerror(nSTAT)

      ! Write data
      ! Time
      !dumscal_out = time
!      nSTAT=nf90_get_var(ncid,t_var_id,dumscal_out,(/init_tstep/))
!      nSTAT=nf90_get_var(ncid,t_var_id,dumscal_out,(/init_tstep/))
!      if(nSTAT.ne.0) &
!        write(global_log ,*)'ERROR: put_var t: ',nf90_strerror(nSTAT)
!      write(global_info,*)dumscal_out
!      time = dumscal_out
        ! ashcon
        ! netCDF standard requires that the unlimited dimension (time)
        ! be the most slowly varying, so ashcon unfortunately has a
        ! different shape than concen
      !ashcon = 0.0_op
      !do i = 1,nsmax
      !  ashcon(:,:,:,i) = concen_pd(1:nxmax,1:nymax,1:nzmax,i,ts1)
      !enddo

      nSTAT=nf90_get_var(ncid,ashcon_var_id,ashcon,  &
               start = (/1,1,1,1,init_tstep/),       &
               count = (/nxmax,nymax,nzmax,n_gs_max,1/))

      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: get_var ashcon: ',nf90_strerror(nSTAT)

      do i = 1,n_gs_max
        concen_pd(1:nxmax,1:nymax,1:nzmax,i,ts1) = real(ashcon(:,:,:,i),kind=ip)
      enddo

      nSTAT=nf90_get_var(ncid,depocon_var_id,depocon,&
               start = (/1,1,1,init_tstep/),       &
               count = (/nxmax,nymax,n_gs_max,1/))
      if(nSTAT.ne.0) &
        write(global_log ,*)'ERROR: put_var depocon: ',nf90_strerror(nSTAT)

      do i = 1,n_gs_max
          ! Here's the conversion from kg/m^2
        concen_pd(1:nxmax,1:nymax,0,i,ts1) = real(depocon(:,:,i),kind=ip)/(dz_vec_pd(0)/KM2_2_M2)
      enddo
      concen_pd(:,:,:,:,ts0) = concen_pd(:,:,:,:,ts1)

      ! Close file
      nSTAT = nf90_close(ncid)
      if(nSTAT.ne.0) &
          write(global_log ,*)'ERROR: close outfile: ',                &
                              nf90_strerror(nSTAT)

      deallocate(dum2d_out,dum3d_out)
      deallocate(ashcon)
      deallocate(depocon)

      write(global_info,*)"Read concentrations from time ",dumscal_out

      end subroutine NC_RestartFile_LoadConcen


