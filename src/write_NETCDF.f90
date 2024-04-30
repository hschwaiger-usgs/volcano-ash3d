!##############################################################################
! 
! Ash3d_Netcdf_IO module
!
! This module manages all input and output of Ash3d data in netcdf format.
!
!      subroutine NC_create_netcdf_file
!      subroutine NC_append_to_netcdf
!      subroutine NC_RestartFile_ReadTimes
!      subroutine NC_RestartFile_LoadConcen
!      subroutine NC_check_status
!      subroutine NC_Read_Output_Products
!
!##############################################################################

      module Ash3d_Netcdf_IO

      use precis_param

      use io_units

      use netcdf

      implicit none

        ! Set everything to private by default
      private

        ! Publicly available subroutines/functions
      public NC_create_netcdf_file,     &
             NC_append_to_netcdf,       &
             NC_RestartFile_ReadTimes,  &
             NC_RestartFile_LoadConcen, &
             NC_Read_Output_Products

        ! Publicly available variables
      integer,public :: tn_len

      integer :: NCversion
      integer :: NCsubversion
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

      !integer :: pj_dim_id    = 0 ! projection parameter dimension
      !  Coordinate variables
      integer :: t_var_id              = 0 ! Time
      integer :: x_var_id              = 0 ! X-distance
      integer :: y_var_id              = 0 ! Y-distance
      integer :: z_var_id              = 0 ! Z-distance
      integer :: s_var_id              = 0 ! shifted (or sigma) altitude
      integer :: bn_var_id             = 0 ! index for species (grain-size bin, gas, water, etc.)
      integer :: er_var_id             = 0 ! eruption index
      integer :: wf_var_id             = 0 ! wind file index
      integer :: pt_var_id             = 0 ! point (airport/POI) index
      integer :: pr_var_id             = 0 ! profile output index
      integer :: tn_var_id             = 0 ! time (native)

      integer :: proj_var_id           = 0 ! Projection
      integer :: FV_var_id             = 0 ! Fall model
      integer :: spec_var_id           = 0 ! Species class ID
      integer :: subspec_var_id        = 0 ! Species sub-class ID

      integer :: vx_var_id             = 0 ! Vx
      integer :: vy_var_id             = 0 ! Vy
      integer :: vz_var_id             = 0 ! Vz
      integer :: vf_var_id             = 0 ! Vf
      integer :: ashcon_var_id         = 0 ! Ash concentration
      integer :: depocon_var_id        = 0 ! Deposit mass/area
      !integer :: gencon_var_id         = 0 ! General concentration for any slices of concen abov
      integer :: cloudmask_var_id      = 0
      integer :: ashconMax_var_id      = 0 ! Max Ash concentration in column
      integer :: ashheight_var_id      = 0 ! Height of top of ash cloud
      integer :: ashload_var_id        = 0 ! Vert. integrated load of ash cloud
      integer :: radrefl_var_id        = 0 ! Radar reflectivity in dbZ
      integer :: depothick_var_id      = 0 ! Total deposit thickness
      integer :: depothickFin_var_id   = 0 ! Final total deposit thickness at simulation end
      integer :: depotime_var_id       = 0 ! Deposit arrival time
      integer :: ashcloudtime_var_id   = 0 ! Cloud arrival time
      integer :: ashcloudBot_var_id    = 0 ! Height of bottom of ash cloud

      integer :: area_var_id           = 0 ! area of cell (km^2)
      integer :: gssd_var_id           = 0 ! Grain diameter
      integer :: gsmf_var_id           = 0 ! Grain mass fraction
      integer :: gsdens_var_id         = 0 ! Grain density
      integer :: gsF_var_id            = 0 ! Grain shape fac F
      integer :: gsG_var_id            = 0 ! Grain shape fac G
      integer :: gsP_var_id            = 0 ! Grain shape fac Phi

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
      integer :: pt_ashthicknessFin_var_id  = 0 ! Final ash thickness at point

      ! Vertical profile variable
      integer :: pr_x_var_id             = 0 ! x coordinate of point
      integer :: pr_y_var_id             = 0 ! y coordinate of point
      integer :: pr_name_var_id          = 0 ! site name of point
      integer :: pr_ash_var_id           = 0 ! concentration profile

      integer :: temp1_2d_var_id !,temp2_2d_var_id,temp3_2d_var_id,temp4_2d_var_id
      integer :: temp1_3d_var_id !,temp2_3d_var_id,temp3_3d_var_id,temp4_3d_var_id
      integer :: temp1_4d_var_id !,temp2_4d_var_id!,temp3_4d_var_id!,temp4_4d_var_id
      integer :: ns_extra

      ! work spaces
      real(kind=op)                                 :: dumscal_out
      integer,       dimension(:)      ,allocatable :: dum1dint_out
      integer,       dimension(:,:)    ,allocatable :: dum2dint_out
      real(kind=op), dimension(:)      ,allocatable :: dum1d_out
      real(kind=op), dimension(:,:)    ,allocatable :: dum2d_out
      real(kind=op), dimension(:,:,:)  ,allocatable :: dum3d_out
      real(kind=dp), dimension(:)      ,allocatable :: dum1d_out_dp
      character                                     :: dumchar
      real(kind=op), dimension(:,:)    ,allocatable :: dum2dPOI_out

      real(kind=op), dimension(:,:,:,:),allocatable :: ashcon
      real(kind=op), dimension(:,:,:)  ,allocatable :: depocon

      contains
      !------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    NC_create_netcdf_file
!
!  Called from: output_results
!  Arguments:
!    none
!
!  This subroutine is called once before time integration starts, creating the output
!  netcdf file, defining all the dimensions and variables.  Variables are filled with
!  initial values.  All details of the control file, and resetable parameters are
!  included as global attributes.  The output netcdf file is closed at the end of this
!  subroutine and reopened as needed in NC_append_to_netcdf.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine NC_create_netcdf_file

      use global_param,  only : &
         EPS_SMALL,KM2_2_M2,M_2_MM,useCalcFallVel,&
         GRAV,CFL,DT_MIN,DT_MAX,RAD_EARTH,Ash3d_GitComID,os_cwd,os_host,os_user,&
         useVz_rhoG

      use io_data,       only : &
         nvprofiles,Site_vprofile,x_vprofile,y_vprofile, &
         cdf_b1l1,cdf_b1l2,cdf_b1l3,cdf_b1l4,cdf_b1l5,cdf_b1l6,cdf_b1l7,cdf_b1l8,cdf_b1l9,&
         cdf_b3l1,cdf_b3l2,cdf_b3l3,cdf_b3l4,cdf_b3l5,cdf_b4l1,cdf_b4l2,cdf_b4l3,cdf_b4l4,&
         cdf_b4l5,cdf_b4l6,cdf_b4l7,cdf_b4l8,cdf_b4l9,cdf_b4l10,cdf_b4l11,cdf_b6l1,cdf_b6l2,&
         cdf_b6l3,cdf_b6l4,cdf_b6l5,cdf_conventions,&
         cdf_comment,cdf_title,cdf_institution,cdf_source,cdf_history,cdf_references,&
         cdf_run_class,cdf_url,infile,concenfile,&
         nvar_User2d_static_XY,nvar_User2d_XY,nvar_User3d_XYGs,nvar_User3d_XYZ,&
         nvar_User4d_XYZGs,Write_PT_Data,Write_PR_Data

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,x_cc_pd,y_cc_pd,z_cc_pd,lon_cc_pd,lat_cc_pd,s_cc_pd,&
         sigma_nz_pd,dx,dy,dz_vec_pd,IsLatLon,ts1,&
         A3d_iprojflag,A3d_k0_scale,A3d_phi0,A3d_lam0,A3d_lam1,A3d_phi1,A3d_lam2,&
         A3d_phi2,A3d_Re,ZPADDING,ZScaling_ID,Ztop

      use solution,      only : &
          vx_pd,vy_pd,vz_pd,vf_pd,concen_pd,DepositGranularity,SpeciesID,SpeciesSubID

      use time_data,     only : &
          BaseYear,useLeap,os_time_log,time,SimStartHour,xmlSimStartTime,OutputOffset

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
         DEPO_THRESH,DEPRATE_THRESH,CLOUDCON_THRESH,CLOUDLOAD_THRESH,&
         THICKNESS_THRESH,DBZ_THRESH,CLOUDCON_GRID_THRESH,&
         DBZ_THRESH,USE_OPTMOD_VARS,useRestartVars,&
         useOutprodVars,useWindVars,DepositThickness,DepArrivalTime,CloudArrivalTime,&
         MaxConcentration,MaxHeight,CloudLoad,dbZ,MinHeight,Mask_Cloud,Mask_Deposit,&
           dbZCalculator

      use Airports,      only : &
         nairports,Airport_Code,Airport_Name,                   &
         Airport_Latitude,Airport_Longitude,Airport_Thickness_TS

      use Tephra,        only : &
         n_gs_max,Tephra_gsdiam,Tephra_bin_mass,Tephra_rho_m,FV_ID,&
         Tephra_gsF,Tephra_gsG,Tephra_gsPhi,Shape_ID,&
         MagmaDensity,DepositDensity,LAM_GS_THRESH,AIRBORNE_THRESH

      use Source,        only : &
         neruptions,e_Volume,e_Duration,e_StartTime,e_PlumeHeight

      use Source_Umbrella, only : &
         VelMod_umb,k_entrainment_umb,lambda_umb,N_BV_umb,SuzK_umb

      use MetReader,     only : &
         MR_iwindfiles,MR_windfiles,MR_GitComID,&
         MR_MetStep_findex,MR_windfile_starthour,Met_proj4,Comp_proj4

      use projection,      only : &
         PJ_GitComID

      integer :: nSTAT

      character(len=32)              :: time_units

      character (len=20)  :: cdf_WindStartTime

      character(len=3 ),dimension(11) :: dim_names
      character(len=30),dimension(11) :: dim_lnames
      character(len=30),dimension(40) :: var_lnames
      character(len=13)  :: reftimestr
      character(len=16)  :: outstring
      character(len= 50) :: linebuffer050
      character(len=130) :: linebuffer130
      integer :: strlen
      integer :: i,j,k,isize
      integer :: ivar
      integer,dimension(5) :: chunksizes5
      logical :: IsThere
      !character(len=3)  :: answer
      integer           :: iostatus
      character(len=120):: iomessage

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

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"     Entered Subroutine NC_create_netcdf_file"
      endif;enddo

      allocate(dum2dint_out(nxmax,nymax))
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
      dim_lnames(5) = "Bin index"
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
      var_lnames(23) = "Grain shape factor F"
      var_lnames(24) = "Grain shape factor G"
      var_lnames(25) = "Grain sphericity"

      var_lnames(30) = "Deposit thickness"
      var_lnames(31) = "Airborne ash arrival time"
      var_lnames(32) = "Ash Concentration (Max)"
      var_lnames(33) = "Cloud Top Height"
      var_lnames(34) = "Cloud Load"
      var_lnames(35) = "Radar Reflectivity"
      var_lnames(36) = "Cloud Bottom Height"
      var_lnames(37) = "Deposit thickness (Final)"

      var_lnames(38) = "Species class ID"
      var_lnames(39) = "Species sub-class ID"
      var_lnames(40) = "Cloud presence flag"

      ! Declare header info

      ! get date and time of start of first wind file
      !cdf_WindStartTime = HS_xmltime(MR_MetStep_Hour_since_baseyear(1),BaseYear,useLeap)
      ! Here we assign the start time of the forecast package
      cdf_WindStartTime = HS_xmltime(MR_windfile_starthour(MR_MetStep_findex(1)),BaseYear,useLeap)

      ! Create and open netcdf file
      inquire(file=concenfile,exist=IsThere)
      if(IsThere)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"Specified output netcdf file already exists."
          write(errlog(io),*)"Output filename requested = ",concenfile
!          write(errlog(io),*)"Would you like to over-write this file? (yes or no)"
        endif;enddo
!        read(input_unit,'(a3)',iostat=iostatus,iomsg=iomessage) answer
!        if (adjustl(trim(answer)).eq.'y'.or.adjustl(trim(answer)).eq.'yes') then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"Over-writing file: ",concenfile
          endif;enddo
!        elseif (adjustl(trim(answer)).eq.'n'.or.adjustl(trim(answer)).eq.'no') then
!          do io=1,2;if(VB(io).le.verbosity_error)then
!            write(errlog(io),*)"Exiting."
!          endif;enddo
!          stop 1
!        else
!          io = 1
!          write(errlog(io),*) 'Sorry, I cannot understand your answer.'
!          write(errlog(io),*) "Expected either 'yes' or 'no', but you provided:",answer
!          stop 1
!        endif

      endif

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Creating netcdf file"
      endif;enddo
      linebuffer130 = trim(nf90_inq_libvers())
      read(linebuffer130,'(i1,a1,i1)',iostat=iostatus,iomsg=iomessage)NCversion,dumchar,NCsubversion
      if(iostatus.ne.0)then
        ! If we couldn't read and verify the library version, assume v3
        NCversion = 3
      endif
#ifdef NC3
      ! library version might be forced to v3 by preprocessor flags
      NCversion = 3
#endif
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Netcdf library version = ",NCversion
      endif;enddo
      if(NCversion.eq.4)then
#ifndef NC3
        nSTAT = nf90_create(concenfile,nf90_netcdf4,ncid,           &
                            cache_nelems = 1000, &
                            cache_size = 32000000)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"create outfile v4:")
#endif
      else
        nSTAT = nf90_create(concenfile,nf90_clobber, ncid)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"create outfile classic:")
      endif

      ! Fill in header info with global attributes
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Filling in header info"
      endif;enddo

      nSTAT = nf90_put_att(ncid,nf90_global,"control_file",infile)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att control_file:")
      nSTAT = nf90_put_att(ncid,nf90_global,"title",cdf_title)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att title:")
      nSTAT = nf90_put_att(ncid,nf90_global,"comment",cdf_comment)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att comment:")
      nSTAT = nf90_put_att(ncid,nf90_global,"run_class",cdf_run_class)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att run_class:")
      nSTAT = nf90_put_att(ncid,nf90_global,"institution",cdf_institution)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"tt institution:")

      nSTAT = nf90_put_att(ncid,nf90_global,"source",cdf_source)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att source:")

      cdf_history=os_time_log
      nSTAT = nf90_put_att(ncid,nf90_global,"history",cdf_history)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att history:")

      nSTAT = nf90_put_att(ncid,nf90_global,"references",cdf_references)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att references:")
      nSTAT = nf90_put_att(ncid,nf90_global,"url",cdf_url)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att url:")

      nSTAT = nf90_put_att(ncid,nf90_global,"Conventions",cdf_conventions)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Conventions:")

      nSTAT = nf90_put_att(ncid,nf90_global,"user",os_user)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att user:")
      nSTAT = nf90_put_att(ncid,nf90_global,"date",os_time_log)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att date:")
      nSTAT = nf90_put_att(ncid,nf90_global,"NWPStartTime",cdf_WindStartTime)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att NWPStartTime:")
      nSTAT = nf90_put_att(ncid,nf90_global,"MepProj4",trim(adjustl(Met_proj4)))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att MepProj4:")
      nSTAT = nf90_put_att(ncid,nf90_global,"CompProj4",trim(adjustl(Comp_proj4)))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att CompProj4:")
      nSTAT = nf90_put_att(ncid,nf90_global,"host",os_host)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att host:")
      nSTAT = nf90_put_att(ncid,nf90_global,"CWD",os_cwd)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att CWD:")
      nSTAT = nf90_put_att(ncid,nf90_global,"BaseYear",BaseYear)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att BaseYear:")
      if(useLeap)then
        nSTAT = nf90_put_att(ncid,nf90_global,"useLeap",1)
      else
        nSTAT = nf90_put_att(ncid,nf90_global,"useLeap",0)
      endif
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att useLeap:")

      nSTAT = nf90_put_att(ncid,nf90_global,"Projection_Git_Commit_ID",PJ_GitComID)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att PJ_GitComID:")
      nSTAT = nf90_put_att(ncid,nf90_global,"Projection_Git_Commit_ID",MR_GitComID)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att MR_GitComID:")
      nSTAT = nf90_put_att(ncid,nf90_global,"Ash3d_Git_Commit_ID",Ash3d_GitComID)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Ash3d_GitComID:")

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
      nSTAT = nf90_put_att(ncid,nf90_global,"b6l5",cdf_b6l5)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment b6l5:")

      ! Now writing out all the resettable parameters
      nSTAT = nf90_put_att(ncid,nf90_global,"MagmaDensity",MagmaDensity)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment MagmaDensity:")
      nSTAT = nf90_put_att(ncid,nf90_global,"DepositDensity",DepositDensity)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment DepositDensity:")
      nSTAT = nf90_put_att(ncid,nf90_global,"LAM_GS_THRESH",LAM_GS_THRESH)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment LAM_GS_THRESH:")
      nSTAT = nf90_put_att(ncid,nf90_global,"AIRBORNE_THRESH",AIRBORNE_THRESH)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment AIRBORNE_THRESH:")
      nSTAT = nf90_put_att(ncid,nf90_global,"GRAV",GRAV)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment GRAV:")
      nSTAT = nf90_put_att(ncid,nf90_global,"RAD_EARTH",RAD_EARTH)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment RAD_EARTH:")
      nSTAT = nf90_put_att(ncid,nf90_global,"CFL",CFL)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment CFL:")
      nSTAT = nf90_put_att(ncid,nf90_global,"DT_MIN",DT_MIN)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment DT_MIN:")
      nSTAT = nf90_put_att(ncid,nf90_global,"DT_MAX",DT_MAX)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment DT_MAX:")
      nSTAT = nf90_put_att(ncid,nf90_global,"ZPADDING",ZPADDING)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment ZPADDING:")
      if(useVz_rhoG)then
        nSTAT = nf90_put_att(ncid,nf90_global,"useVz_rhoG","true")
      else
        nSTAT = nf90_put_att(ncid,nf90_global,"useVz_rhoG","false")
      endif
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment useVz_rhoG:")
      nSTAT = nf90_put_att(ncid,nf90_global,"DEPO_THRESH",DEPO_THRESH)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment DEPO_THRESH:")
      nSTAT = nf90_put_att(ncid,nf90_global,"DEPRATE_THRESH",DEPRATE_THRESH)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment DEPRATE_THRESH:")
      nSTAT = nf90_put_att(ncid,nf90_global,"CLOUDCON_THRESH",CLOUDCON_THRESH)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment CLOUDCON_THRESH:")
      nSTAT = nf90_put_att(ncid,nf90_global,"CLOUDCON_GRID_THRESH",CLOUDCON_GRID_THRESH)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment CLOUDCON_GRID_THRESH:")
      nSTAT = nf90_put_att(ncid,nf90_global,"CLOUDLOAD_THRESH",CLOUDLOAD_THRESH)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment CLOUDLOAD_THRESH:")
      nSTAT = nf90_put_att(ncid,nf90_global,"THICKNESS_THRESH",THICKNESS_THRESH)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment THICKNESS_THRESH:")
      nSTAT = nf90_put_att(ncid,nf90_global,"DBZ_THRESH",DBZ_THRESH)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment DBZ_THRESH:")
      ! Parameters for umbrella clouds
      nSTAT = nf90_put_att(ncid,nf90_global,"VelMod_umb",VelMod_umb)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment VelMod_umb:")
      nSTAT = nf90_put_att(ncid,nf90_global,"lambda_umb",lambda_umb)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment lambda_umb:")
      nSTAT = nf90_put_att(ncid,nf90_global,"N_BV_umb",N_BV_umb)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment N_BV_umb:")
      nSTAT = nf90_put_att(ncid,nf90_global,"k_entrainment_umb",k_entrainment_umb)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment k_entrainment_umb:")
      nSTAT = nf90_put_att(ncid,nf90_global,"SuzK_umb",SuzK_umb)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Comment SuzK_umb:")

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
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Defining dimensions"
      endif;enddo
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
      !if (nairports.gt.0)then
      if (Write_PT_Data)then
        nSTAT = nf90_def_dim(ncid,dim_names(9),nairports,pt_dim_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_dim pt:")
      endif

       ! pr
      if (Write_PR_Data)then
        nSTAT = nf90_def_dim(ncid,dim_names(10),nvprofiles,pr_dim_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_dim pr:")
      endif

      ! Define coordinate variables
        ! X,Y,Z,time,bn,er,wf,sl,pt

      ! Define time-dependent variables
        ! Vx,Vy,Vz
        ! ashconc,depdepth
      ! Define other variables
        ! gs_setvel, gs_massfrac, gs_dens
         ! Time
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Defining coordinate variables"
        write(outlog(io),*)"     Time: ",dim_names(1)
      endif;enddo
      ! Time variables should always be doubles to match with libhourssince
      nSTAT = nf90_def_var(ncid,dim_names(1),&
                           nf90_double,& 
                           (/t_dim_id/),&
                           t_var_id)
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
      do io=1,2;if(VB(io).le.verbosity_info)then
        if (ZScaling_ID.eq.0) then
          write(outlog(io),*)"     Z: ",dim_names(2)
        else
          write(outlog(io),*)"     S: ","s"
        endif
      endif;enddo
      ! This is the normal z=altitude coordinate (assumes topography = 0)
      ! Always written
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
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"     Y: ",dim_names(3)
      endif;enddo
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
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"     X: ",dim_names(4)
      endif;enddo
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

      if (ZScaling_ID.gt.0) then
        ! This branch is for z-shifting/scaling
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"s",&
                               nf90_double,&
                               (/z_dim_id/),&
                               s_var_id)
        else
          nSTAT = nf90_def_var(ncid,"s",&
                               nf90_float,&
                               (/z_dim_id/), &
                               s_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var s")
        if(ZScaling_ID.eq.1)then
          nSTAT = nf90_put_att(ncid,s_var_id,"long_name","shifted-altitude")
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att s long_name")
          nSTAT = nf90_put_att(ncid,s_var_id,"units","km")
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att s units")
          nSTAT = nf90_put_att(ncid,s_var_id,"note","s=z-zsurf")
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att s note")
        elseif(ZScaling_ID.eq.2)then
          nSTAT = nf90_put_att(ncid,s_var_id,"long_name","sigma-altitude")
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att s long_name")
          nSTAT = nf90_put_att(ncid,s_var_id,"units","none")
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att s units")
          nSTAT = nf90_put_att(ncid,s_var_id,"note","s=(z-zsurf)/(top-surf)")
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att s note")
        endif
        nSTAT = nf90_put_att(ncid,s_var_id,"positive","up")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att s positive")
        nSTAT = nf90_put_att(ncid,s_var_id,"ztop",Ztop)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att s ztop")
      endif

         ! BN (Grain size bin ID)
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"     Bin: ",dim_names(5)
      endif;enddo
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
                                          "index for grainsizes, gas, water, etc.")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att bn Comment")

         ! ER (Eruption index)
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"     ER: ",dim_names(6)
      endif;enddo
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
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"     WF: ",dim_names(7)
      endif;enddo
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
      !if (nairports.gt.0)then
      if (Write_PT_Data)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     PT: ",dim_names(9)
       endif;enddo
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

!        ! pr (profile output index)
      if (Write_PR_Data)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     PR: ",dim_names(10)
        endif;enddo
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
        ! We can only have a second unlimited dimension with NC version 4.
        ! With NC version 3, we need to know the dimension length first, so
        ! we define this dimension in NC_append_to_netcdf, but only for the
        ! final step.
      endif
      !   End of dimension variables
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Defining variables"
      endif;enddo

      ! write projection as a variable
      if (IsLatLon) then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     LatLon_Projection"
        endif;enddo
        nSTAT = nf90_def_var(ncid,"LatLon_Projection",&
                             nf90_int,&
                             proj_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var LatLon_Projection")
        nSTAT = nf90_put_att(ncid,proj_var_id,"grid_mapping_name","latitude_longitude")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att LatLon_Projection grid_mapping_name")
        nSTAT = nf90_put_att(ncid,proj_var_id,"semi_major_axis",A3d_Re*1000.0_ip)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att LatLon_Projection semi_major_axis")
        nSTAT = nf90_put_att(ncid,proj_var_id,"inverse_flattening","0")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att LatLon_Projection inverse_flattening")
      else
        select case (A3d_iprojflag)
        case(0)
          ! Non-geographic projection, (x,y) only
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"     Projection : Non-geographic"
          endif;enddo
          nSTAT = nf90_def_var(ncid,"Non-geographic",&
                               nf90_int,&
                               proj_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var Non-geographic")
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                               "grid_mapping_name", &
                               "Non-geographic")
        case(1)
          ! Polar stereographic
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"     Projection : Polar_Stereographic"
          endif;enddo
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
                              "earth_radius",A3d_Re*1000.0_ip)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1, &
                              "put_att Polar_Stereographic earth_radius")
        case(2)
          ! Albers Equal Area
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"     Projection : Albers Equal Area"
          endif;enddo
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
                              "earth_radius",A3d_Re*1000.0_ip)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att Albers_Equal_Area earth_radius")
        case(3)
          ! UTM

        case(4)
          ! Lambert conformal conic 
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"     Projection : Lambert_Conformal"
          endif;enddo
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
                              "earth_radius",A3d_Re*1000.0_ip)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1, &
                             "put_att Lambert_Conformal earth_radius")
        case(5)
          ! Mercator
        !        Mercator:grid_mapping_name = "mercator" ;
        !        Mercator:standard_parallel = 20. ;
        !        Mercator:longitude_of_projection_origin = 198.475006103516 ;
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"     Projection : Mercator"
          endif;enddo
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
        case default ! Just write all projection parameters to file
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"     Projection : Not specified"
          endif;enddo
          nSTAT = nf90_def_var(ncid,"Projection",&
                               nf90_int,&
                               proj_var_id)
          !nSTAT = nf90_put_att(ncid,proj_var_id,&
          !                     "grid_mapping_name", &
          !                     "lambert_conformal_conic")
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                               "parallel0",A3d_phi0)
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                               "parallel1",A3d_phi1)
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                               "parallel2",A3d_phi2)
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                               "meridian0",A3d_lam0)
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                               "meridian1",A3d_lam1)
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                               "meridian2",A3d_lam2)
          nSTAT = nf90_put_att(ncid,proj_var_id,&
                              "earth_radius",A3d_Re*1000.0_ip)
        end select
      endif

      ! Create variable defining the fall velocity model used
      nSTAT = nf90_def_var(ncid,"Fall_Model",&
                           nf90_int,&
                           FV_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var Fall_Model")
      select case (FV_ID)
      case(0)
        ! No fall (just tracer)
        nSTAT = nf90_put_att(ncid,FV_var_id,&
                             "model", &
                             "No fall (just tracer)")
      case(1)
        ! Wilson and Huang
        nSTAT = nf90_put_att(ncid,FV_var_id,&
                             "model", &
                             "Wilson and Huang")
        nSTAT = nf90_put_att(ncid,FV_var_id,&
                             "doi","10.1016/0012-821X(79)90179-1")
      case(2)
        ! Wilson and Huang + Cunningham slip
        nSTAT = nf90_put_att(ncid,FV_var_id,&
                             "model", &
                             "Wilson and Huang + Cunningham slip")
        nSTAT = nf90_put_att(ncid,FV_var_id,&
                             "doi","10.1016/0012-821X(79)90179-1")
      case(3)
        ! Wilson and Huang + Mod by Pfeiffer Et al.
        nSTAT = nf90_put_att(ncid,FV_var_id,&
                             "model", &
                             "Wilson and Huang + Mod by Pfeiffer Et al.")
        nSTAT = nf90_put_att(ncid,FV_var_id,&
                             "doi","10.1016/j.jvolgeores.2004.09.001")
      case(4)
        ! Ganser
        nSTAT = nf90_put_att(ncid,FV_var_id,&
                             "model", &
                             "Ganser")
        nSTAT = nf90_put_att(ncid,FV_var_id,&
                             "doi","10.1016/0032-5910(93)80051-B");
      case(5)
        ! Ganser + Cunningham slip
        nSTAT = nf90_put_att(ncid,FV_var_id,&
                             "model", &
                             "Ganser + Cunningham slip")
        nSTAT = nf90_put_att(ncid,FV_var_id,&
                             "doi","10.1016/0032-5910(93)80051-B");
      case(6)
        ! Stokes flow for spherical particles + slip
        nSTAT = nf90_put_att(ncid,FV_var_id,&
                             "model", &
                             "Stokes,flow for spherical particles + slip")
      case default
        ! Wilson and Huang
        nSTAT = nf90_put_att(ncid,FV_var_id,&
                             "model", &
                             "Wilson and Huang")
        nSTAT = nf90_put_att(ncid,FV_var_id,&
                             "doi","10.1016/0012-821X(79)90179-1")
      end select

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !   Now a few other variables that are a function of BN
         ! Species class ID
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"     SC: Species Class"
      endif;enddo
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
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"     SSC: Species Sub-class"
      endif;enddo
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
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"     gs_diameter: grain diameter"
      endif;enddo
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
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"     gs_massfrac: grain mass fraction"
      endif;enddo
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
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"     gs_dens: grain density"
      endif;enddo
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

         ! gs_F (Shape factor of grain: F = (B+C)/2A)
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"     gs_F: grain shape factor (B+C)/2A"
      endif;enddo
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,"gs_F",&
                             nf90_double,&
                             (/bn_dim_id/),&
                             gsF_var_id)
      else
        nSTAT = nf90_def_var(ncid,"gs_F",&
                             nf90_float,&
                             (/bn_dim_id/), &
                             gsF_var_id)
      endif
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var gs_F")
      nSTAT = nf90_put_att(ncid,gsF_var_id,"long_name",var_lnames(23))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att gs_F long_name")
      nSTAT = nf90_put_att(ncid,gsF_var_id,"units","none")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att gs_F units")
      nSTAT = nf90_put_att(ncid,gsF_var_id,"note","F=(C+B)/2A")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att gs_F note")
      if (Shape_ID.eq.1) then
        nSTAT = nf90_put_att(ncid,gsF_var_id,"source","provided")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att gs_F source")
      else
        nSTAT = nf90_put_att(ncid,gsF_var_id,"source","calculated")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att gs_F source")
      endif

         ! gs_G (Shape factor of grain: G = B/C)
      do io=1,2;if(VB(io).le.verbosity_info)then      
        write(outlog(io),*)"     gs_G: grain shape factor B/C"
      endif;enddo
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,"gs_G",&
                             nf90_double,&
                             (/bn_dim_id/),&
                             gsG_var_id)
      else
        nSTAT = nf90_def_var(ncid,"gs_G",&
                             nf90_float,&
                             (/bn_dim_id/), &
                             gsG_var_id)
      endif
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var gs_G")
      nSTAT = nf90_put_att(ncid,gsG_var_id,"long_name",var_lnames(24))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att gs_G long_name")
      nSTAT = nf90_put_att(ncid,gsG_var_id,"units","none")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att gs_G units")
      nSTAT = nf90_put_att(ncid,gsG_var_id,"note","G=C/B")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att gs_G note")
      if (Shape_ID.eq.1) then
        nSTAT = nf90_put_att(ncid,gsG_var_id,"source","provided")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att gs_F source")
      else
        nSTAT = nf90_put_att(ncid,gsG_var_id,"source","calculated")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att gs_F source")
      endif

         ! gs_Phi (Sphericity of grain: Psi = Area-of-vol.eq.sphere/Area-of-particle)
      do io=1,2;if(VB(io).le.verbosity_info)then      
        write(outlog(io),*)"     gs_Phi: grain sphericity Psi = Area-of-vol.eq.sphere/Area-of-particle"
      endif;enddo
      if(op.eq.8)then
        nSTAT = nf90_def_var(ncid,"gs_Phi",&
                             nf90_double,&
                             (/bn_dim_id/),&
                             gsP_var_id)
      else
        nSTAT = nf90_def_var(ncid,"gs_Phi",&
                             nf90_float,&
                             (/bn_dim_id/), &
                             gsP_var_id)
      endif
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var gs_P")
      nSTAT = nf90_put_att(ncid,gsP_var_id,"long_name",var_lnames(25))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att gs_P long_name")
      nSTAT = nf90_put_att(ncid,gsP_var_id,"units","none")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att gs_P units")
      nSTAT = nf90_put_att(ncid,gsP_var_id,"note","Phi=(Area of vol.equ.sphere)/(Area of particle)")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att gs_P note")
      if (Shape_ID.eq.1) then
        nSTAT = nf90_put_att(ncid,gsP_var_id,"source","calculated")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att gs_P source")
      else
        nSTAT = nf90_put_att(ncid,gsP_var_id,"source","provided")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att gs_P source")
      endif

      !   Now a few other variables that are a function of ER
         ! er_stime (Start time of eruption)
      ! Time variables should always be doubles to match with libhourssince
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"     er_stime: Eruption pulse start time"
      endif;enddo
      nSTAT = nf90_def_var(ncid,"er_stime",&
                           nf90_double,&
                           (/er_dim_id/),&
                           er_stime_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var er_stime")
      nSTAT = nf90_put_att(ncid,er_stime_var_id,"long_name",var_lnames(15))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att er_stime long_name")
      write(outstring,4502)BaseYear
 4502 format('hours since ',i4)
      nSTAT = nf90_put_att(ncid,er_stime_var_id,"units",outstring)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att er_stime units")

         ! er_duration (Duration of eruption)
      ! Time variables should always be doubles to match with libhourssince
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"     er_duration: Eruption pulse duration"
      endif;enddo
      nSTAT = nf90_def_var(ncid,"er_duration",&
                           nf90_double,&
                           (/er_dim_id/),&
                           er_duration_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var er_duration")
      nSTAT = nf90_put_att(ncid,er_duration_var_id,"long_name",var_lnames(16))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att er_duration long_name")
      nSTAT = nf90_put_att(ncid,er_duration_var_id,"units", "hours")
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att er_duration units")

         ! er_plumeheight (Plume height of eruption)
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"     er_plumeheight: Eruption pulse Plume Height"
      endif;enddo
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
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"     er_volume: Eruption pulse volume (DRE)"
      endif;enddo
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
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"     area: cell area"
      endif;enddo
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


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(useWindVars)then
         ! Now define the time-dependent variables
         ! Vx
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     vx: ",var_lnames(8)
        endif;enddo
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
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     vy: ",var_lnames(9)
        endif;enddo
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
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     vz: ",var_lnames(10)
        endif;enddo
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
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     vf: ","Fall Velocity"
        endif;enddo
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
      endif  ! useWindVars
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(useRestartVars)then
         ! Full Ash Concentration
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     ashcon :",var_lnames(11)
        endif;enddo
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
#ifndef NC3
          chunksizes5 = (/nxmax, nymax, nzmax, nsmax, 1/)
          nSTAT = nf90_def_var_chunking(ncid, ashcon_var_id, &
                                        0, &
                                        chunksizes5)
          nSTAT = nf90_def_var_deflate(ncid, ashcon_var_id,          &
                                  shuffle = 1,                &
                                  deflate = 5,                &
                                  deflate_level = 5  )
#endif
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
      endif ! useRestartVars
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Derived variables
      if(useOutprodVars)then
        ! Deposit(as function of grainsmax)
           ! Concentration of deposit by grain size(z=0 plane of tot_ashcon)
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     depocon: ",var_lnames(12)
        endif;enddo
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

        ! Deposit thickness
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     depothick: ",var_lnames(30)
        endif;enddo
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

        ! Depositthickness (Final)
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     depothickFin: ",var_lnames(37)
        endif;enddo
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"depothickFin",&
                               nf90_double, &
                               (/x_dim_id,y_dim_id/),                &
                               depothickFin_var_id)
        else
          nSTAT = nf90_def_var(ncid,"depothickFin",&
                               nf90_float,  &
                               (/x_dim_id,y_dim_id/),                &
                               depothickFin_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var depothickFin")
        nSTAT = nf90_put_att(ncid,depothickFin_var_id,"long_name",var_lnames(30))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att depothickFin long_name")
        nSTAT = nf90_put_att(ncid,depothickFin_var_id,"units","mm")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att depothickFin units")
        nSTAT = nf90_put_att(ncid,depothickFin_var_id,&
                 "missing_value", DepositThickness_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att depothickFin missing_value")
        nSTAT = nf90_put_att(ncid,depothickFin_var_id,"_FillValue",DepositThickness_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att depothickFin _FillValue")

           ! Arrival time of deposit
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     depotime: ",var_lnames(21)
        endif;enddo
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
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     ash_arrival_time: ",var_lnames(31)
        endif;enddo
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
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     ashcon_max: ",var_lnames(32)
        endif;enddo
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
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     cloud_height: ",var_lnames(33)
        endif;enddo
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
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     cloud_load: ",var_lnames(34)
        endif;enddo
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

           ! Cloud Mask (this is an integer array flagging where the CloudLoad is non-zero
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     cloud_mask: 1=volcanic cloud, 0=clear air"
        endif;enddo
        nSTAT = nf90_def_var(ncid,"cloud_mask",&
                             nf90_int,  &
                             (/x_dim_id,y_dim_id,t_dim_id/),                &
                             cloudmask_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var cloud_mask")
        nSTAT = nf90_put_att(ncid,cloudmask_var_id,"long_name",var_lnames(40))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att cloud_mask long_name")
        nSTAT = nf90_put_att(ncid,cloudmask_var_id,"units","0 or 1") !
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att cloud_mask units")

           ! 3d radar reflectivity
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     radar_reflectivity: ",var_lnames(35)
        endif;enddo
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
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     cloud_bottom: ",var_lnames(36)
        endif;enddo
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

      endif ! useOutprodVars
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Point output (Airport/POI)
      !if (nairports.gt.0)then
      if (Write_PT_Data)then
        ! x coordinate of point
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     pt_x: x (or lon) of Airport/POI data"
        endif;enddo
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
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     pt_y: y (or lat) of Airport/POI data"
        endif;enddo
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
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     pt_code: Airport/POI data 3-character code"
        endif;enddo
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
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     pt_name: Airport/POI data location name"
        endif;enddo
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
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     pt_depotime: Airport/POI time of deposit arrival"
        endif;enddo
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
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     pt_depodur: Airport/POI duration of deposition"
        endif;enddo
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
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     pt_cloud_arrival: Airport/POI time of airborne cloud arrival"
        endif;enddo
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
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     pt_cloud_dur: Airport/POI duration of airborne cloud transit"
        endif;enddo
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
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     pt_depothick: Airport/POI time-series of deposit accumulation"
        endif;enddo
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

        ! Point/Airport ashfall thickness
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     pt_depothickFin: Airport/POI final deposit thickness"
        endif;enddo
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"pt_depothickFin",&
                               nf90_double, &
                               (/pt_dim_id/),                &
                               pt_ashthicknessFin_var_id)
        else
          nSTAT = nf90_def_var(ncid,"pt_depothickFin",&
                               nf90_float,  &
                               (/pt_dim_id/),                &
                               pt_ashthicknessFin_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_depothickFin")
        nSTAT = nf90_put_att(ncid,pt_ashthicknessFin_var_id,"long_name","Final ash fall thickness")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_ashthicknessFin_var_id long_name")
        nSTAT = nf90_put_att(ncid,pt_ashthicknessFin_var_id,"units","mm")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_ashthicknessFin_var_id units")
        nSTAT = nf90_put_att(ncid,pt_ashthicknessFin_var_id,&
                   "missing_value", DepArrivalTime_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_ashthicknessFin_var_id missing_value")
        nSTAT = nf90_put_att(ncid,pt_ashthicknessFin_var_id,"_FillValue",DepArrivalTime_FillValue)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pt_ashthicknessFin_var_id _FillValue")

      endif

      if(Write_PR_Data)then
        ! x coordinate of point
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     pr_x: x (or lon) coordinate of vertical profile"
        endif;enddo
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"pr_x",&
                               nf90_double,   &
                               (/pr_dim_id/), &
                               pr_x_var_id)
        else
          nSTAT = nf90_def_var(ncid,"pr_x",&
                               nf90_float,  &
                               (/pr_dim_id/),                &
                               pr_x_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pr_x")
        nSTAT = nf90_put_att(ncid,pr_x_var_id,"long_name","x coordinate of prof point")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pr_x long_name")
        if(IsLatLon)then
           nSTAT = nf90_put_att(ncid,pr_x_var_id,"units","degrees_east")
        else
          nSTAT = nf90_put_att(ncid,pr_x_var_id,"units","km")
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pr_x units")

        ! y coordinate of point
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     pr_y: y (or lat) coordinate of vertical profile"
        endif;enddo
        if(op.eq.8)then
          nSTAT = nf90_def_var(ncid,"pr_y",&
                               nf90_double,   &
                               (/pr_dim_id/), &
                               pr_y_var_id)
        else
          nSTAT = nf90_def_var(ncid,"pr_y",&
                               nf90_float,  &
                               (/pr_dim_id/),                &
                               pr_y_var_id)
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pr_y")
        nSTAT = nf90_put_att(ncid,pr_y_var_id,"long_name","y coordinate of prof point")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pr_y long_name")
        if(IsLatLon)then
           nSTAT = nf90_put_att(ncid,pr_y_var_id,"units","degrees_north")
        else
          nSTAT = nf90_put_att(ncid,pr_y_var_id,"units","km")
        endif
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pr_y units")

        ! Name of point
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     pr_name: Vertical profile name"
        endif;enddo
        nSTAT = nf90_def_var(ncid,"pr_name",&
                             nf90_char,   &
                             (/sl_dim_id,pr_dim_id/), &
                             pr_name_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var pr_name")
        nSTAT = nf90_put_att(ncid,pr_y_var_id,"long_name","Name of prof point")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pr_name long_name")
        nSTAT = nf90_put_att(ncid,pr_name_var_id,"units","text")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pr_name units")

        ! Profile data
        ! Variable for profile data (pr_ash) will be defined in NC_append_to_netcdf
        ! in the finalization step.
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(USE_OPTMOD_VARS)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     USE_OPTMOD_VARS"
        endif;enddo
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
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Leaving define mode"
      endif;enddo
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Fill variables with initial values
        ! Time
      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Fill time"
      endif;enddo
      nSTAT=nf90_put_var(ncid,t_var_id,real(time,kind=dp),(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var time")
        ! Z
      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Fill Z"
      endif;enddo
      allocate(dum1d_out(nzmax))
      dum1d_out(:) = real(z_cc_pd(1:nzmax),kind=op)
      nSTAT=nf90_put_var(ncid,z_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var z")
      deallocate(dum1d_out)
        ! Y
      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Fill Y"
      endif;enddo
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
      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Fill X"
      endif;enddo
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

        ! S
      if (ZScaling_ID.gt.0) then
        do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"     Fill S"
        endif;enddo
        allocate(dum1d_out(nzmax))
        dum1d_out(:) = real(s_cc_pd(1:nzmax),kind=op)
        nSTAT=nf90_put_var(ncid,s_var_id,dum1d_out,(/1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var s")
        deallocate(dum1d_out)
      endif

        ! BN (Grain size bin ID)
      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Fill BN"
      endif;enddo
      allocate(dum1dint_out(nsmax))
      do i=1,nsmax  ! This is just an index
        dum1dint_out(i) = i
      enddo
      nSTAT=nf90_put_var(ncid,bn_var_id,dum1dint_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var bn")
      deallocate(dum1dint_out)
        ! ER
      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Fill ER"
      endif;enddo
      allocate(dum1dint_out(neruptions))
      ! This is variable associated with the dimension for eruptions ID
      ! This only contains the index starting with 1
      do i=1,neruptions
        dum1dint_out(i) = i
      enddo
      nSTAT=nf90_put_var(ncid,er_var_id,dum1dint_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var er")
      deallocate(dum1dint_out)
        ! WF
      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Fill WF"
      endif;enddo
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
      if (Write_PT_Data)then
        do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"     Fill PT"
        endif;enddo
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

      if(Write_PR_Data)then
          ! PR
        do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"     Fill PR"
        endif;enddo
        allocate(dum1dint_out(nvprofiles))
        ! This is variable associated with the dimension for profile output
        ! This only contains the index starting with 1
        do i=1,nvprofiles
          dum1dint_out(i) = i
        enddo
        nSTAT=nf90_put_var(ncid,pr_var_id,dum1dint_out,(/1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pr")
        deallocate(dum1dint_out)
      endif

      !   End of filling dimension variables
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !   Now fill a few other variables that are a function of BN
         ! Species class ID
      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Fill Species class ID"
      endif;enddo
      allocate(dum1dint_out(nsmax))
      dum1dint_out(1:nsmax) = SpeciesID(1:nsmax)
      nSTAT=nf90_put_var(ncid,spec_var_id,dum1dint_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var spec_class")
      deallocate(dum1dint_out)
         ! Species sub-class ID
      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Fill Species sub-class ID"
      endif;enddo
      allocate(dum1dint_out(nsmax))
      dum1dint_out(1:nsmax) = SpeciesSubID(1:nsmax)
      nSTAT=nf90_put_var(ncid,subspec_var_id,dum1dint_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var spec_subclass")
      deallocate(dum1dint_out)
         ! Grain-size diameter
      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Fill GS Diameter"
      endif;enddo
      allocate(dum1d_out(nsmax))
      dum1d_out = 0.0_op
      if(useCalcFallVel)then
        dum1d_out(1:n_gs_max) = real(Tephra_gsdiam(1:n_gs_max)*M_2_MM,kind=op)
      else
        do isize=1,n_gs_max
          dum1d_out(isize) = real(isize,kind=op)
        enddo
      endif
      nSTAT=nf90_put_var(ncid,gssd_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var gs_diameter")
         ! gs_massfrac (Mass fraction of grain size)
      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Fill GS MassFrac"
      endif;enddo
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
       ! gs_F (Shape factor of grain F)
      if(useCalcFallVel)then
        dum1d_out(1:n_gs_max) = real(Tephra_gsF(1:n_gs_max),kind=op)
      else
        dum1d_out = 0.0_op
      endif
      nSTAT=nf90_put_var(ncid,gsF_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var gs_F")

       ! gs_G (Shape factor of grain G
      if(useCalcFallVel)then
        dum1d_out(1:n_gs_max) = real(Tephra_gsG(1:n_gs_max),kind=op)
      else
        dum1d_out = 0.0_op
      endif
      nSTAT=nf90_put_var(ncid,gsG_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var gs_G")
       ! gs_Phi (Shape factor of grain Phi)
      if(useCalcFallVel)then
        dum1d_out(1:n_gs_max) = real(Tephra_gsPhi(1:n_gs_max),kind=op)
      else
        dum1d_out = 0.0_op
      endif
      nSTAT=nf90_put_var(ncid,gsP_var_id,dum1d_out,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var gs_Phi")
      deallocate(dum1d_out)

      !   Now fill a few other variables that are a function of ER
        ! er_stime (Start time of eruption)
      allocate(dum1d_out(neruptions))
      !dum1d_out = real(e_StartTime + SimStartHour,kind=op)
      nSTAT=nf90_put_var(ncid,er_stime_var_id,&
                         real(e_StartTime + SimStartHour,kind=dp) ,(/1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var er_stime")
        ! er_duration (Duration of eruption)
      !dum1d_out = real(e_Duration,kind=op)
      nSTAT=nf90_put_var(ncid,er_duration_var_id,real(e_Duration,kind=dp),(/1/))
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

      ! And a few variables that are not functions of any dimensions
        ! projection flag (name and attributes already defind this)
      nSTAT=nf90_put_var(ncid,proj_var_id,A3d_iprojflag)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var proj_var_id")
        ! Fall Model ID
      nSTAT=nf90_put_var(ncid,FV_var_id,FV_ID)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var FV_var_id")

         ! Now fill the other (non-time-dependent) variables
         ! wf_name (Name of windfile)
      do i=1,MR_iwindfiles
        write(linebuffer130,'(130a)')trim(adjustl(MR_windfiles(i)))
        strlen = len(trim(adjustl(MR_windfiles(i))))
        do j=1,strlen
          nSTAT=nf90_put_var(ncid,wf_name_var_id,linebuffer130(j:j),(/j,i/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var wf_name")
        enddo
      enddo

         ! Cell area
      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Fill area"
      endif;enddo
      if (IsLatLon) then
        dum2d_out(1:nxmax,1:nymax) = real(sigma_nz_pd(1:nxmax,1:nymax,0),kind=op)
      else
        dum2d_out(1:nxmax,1:nymax) = real(dx*dy,kind=op)
      endif
      nSTAT=nf90_put_var(ncid,area_var_id,dum2d_out)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var area")

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(useWindVars)then
          ! Vz
        do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"     Fill Vz"
        endif;enddo
        dum3d_out(1:nxmax,1:nymax,1:nzmax) = real(vz_pd(1:nxmax,1:nymax,1:nzmax),kind=op)
        nSTAT=nf90_put_var(ncid,vz_var_id,dum3d_out,(/1,1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var vz")
          ! Vy
        do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"     Fill Vy"
        endif;enddo
        dum3d_out(1:nxmax,1:nymax,1:nzmax) = real(vy_pd(1:nxmax,1:nymax,1:nzmax),kind=op)
        nSTAT=nf90_put_var(ncid,vy_var_id,dum3d_out,(/1,1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var vy")
          ! Vx
        do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"     Fill Vx"
        endif;enddo
        dum3d_out(1:nxmax,1:nymax,1:nzmax) = real(vx_pd(1:nxmax,1:nymax,1:nzmax),kind=op)
        nSTAT=nf90_put_var(ncid,vx_var_id,dum3d_out,(/1,1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var vx")
          ! Vf (Note, using ashcon as a dummy variable for this 4-d output)
        ashcon(1:nxmax,1:nymax,1:nzmax,1:nsmax) = &
             real(vf_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax),kind=op)
        nSTAT=nf90_put_var(ncid,vf_var_id,ashcon,(/1,1,1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var vf")
      endif  ! useWindVars
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(useRestartVars)then
        do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"     Fill ashcon"
        endif;enddo
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
      endif  ! useRestartVars
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(useOutprodVars)then
          ! depocon
        do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"     Fill depocon"
        endif;enddo
        depocon = 0.0_op
          ! netCDF standard requires that the unlimited dimension (time)
          ! be the most slowly varying, so dum3d_out unfortunately has a
          ! different shape than depocon
        do isize=1,n_gs_max
            ! Here's the conversion to kg/m^2 from kg/km^2
          depocon(1:nxmax,1:nymax,isize) = real(DepositGranularity(1:nxmax,1:nymax,isize) * &
                                         dz_vec_pd(0)/KM2_2_M2,kind=op)
        enddo
        do isize=1,n_gs_max
          do i=1,nxmax
            do j=1,nymax
              if(depocon(i,j,isize).le.EPS_SMALL)depocon(i,j,isize)=0.0_op
            enddo
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,depocon_var_id,depocon,(/1,1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var depocon")

        !do io=1,2;if(VB(io).le.verbosity_info)then
        !  write(outlog(io),*)"     Fill "
        !endif;enddo
        !call dbZCalculator            ! get radar reflectivity

        ! depothick
        do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"     Fill depothick"
        endif;enddo
        dum2d_out(:,:) = real(DepositThickness,kind=op)
        nSTAT=nf90_put_var(ncid,depothick_var_id,dum2d_out,(/1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var depothick")

        ! ashconMax
        do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"     Fill ashconMax"
        endif;enddo
        dum2d_out(:,:) = real(MaxConcentration,kind=op)
        nSTAT=nf90_put_var(ncid,ashconMax_var_id,dum2d_out,(/1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var ashconMax")

        ! ash cloud_height (top)
        do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"     Fill ash_height"
        endif;enddo
        dum2d_out(:,:) = real(MaxHeight_FillValue,kind=op)
        dum2d_out = merge(real(MaxHeight,kind=op),dum2d_out,Mask_Cloud)
        nSTAT=nf90_put_var(ncid,ashheight_var_id,dum2d_out,(/1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var ashheight")

        ! ash-load
        do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"     Fill ashload"
        endif;enddo
        dum2d_out(:,:) = real(CloudLoad,kind=op)
        nSTAT=nf90_put_var(ncid,ashload_var_id,dum2d_out,(/1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var ashload")

        ! cloud-mask
        do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"     Fill cloud mask"
        endif;enddo
        dum2dint_out = merge(1,0,Mask_Cloud)
        nSTAT=nf90_put_var(ncid,cloudmask_var_id,dum2dint_out,(/1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var cloud_mask")

        ! radar-reflectivity
        do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"     Fill dbZ"
        endif;enddo
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
        do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"     Fill ash height min"
        endif;enddo
        dum2d_out(:,:) = real(MinHeight_FillValue,kind=op)
        dum2d_out = merge(real(MinHeight,kind=op),dum2d_out,Mask_Cloud)
        nSTAT=nf90_put_var(ncid,ashcloudBot_var_id,dum2d_out,(/1,1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var cloud-bottom")

        ! Final time outputs 
        !  Note: this may not align with the last output time step since it is
        !  the values at the end of the simulation, not on the last output step
        ! depothickFin
        do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"     Fill depothickFin"
        endif;enddo
        dum2d_out(:,:) = real(DepositThickness,kind=op)
        nSTAT=nf90_put_var(ncid,depothickFin_var_id,dum2d_out,(/1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var depothickFin")

          ! depotime
        do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"     Fill depotime"
        endif;enddo
        dum2d_out(:,:) = DepArrivalTime_FillValue
        dum2d_out = merge(real(DepArrivalTime,kind=op),dum2d_out,Mask_Deposit)
        nSTAT=nf90_put_var(ncid,depotime_var_id,dum2d_out,(/1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var depotime")

          ! ashtime
       do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"     Fill ashtime"
        endif;enddo
        dum2d_out(:,:) = CloudArrivalTime_FillValue
        do i=1,nxmax
          do j=1,nymax
            if(CloudArrivalTime(i,j).ge.0.0)&
                dum2d_out(i,j)=real(CloudArrivalTime(i,j),kind=op)
          enddo
        enddo
        nSTAT=nf90_put_var(ncid,ashcloudtime_var_id,dum2d_out,(/1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var ashcloudtime_")

      endif ! useOutprodVars

      if (Write_PT_Data)then
        ! These are variable associated with the dimension for point output
        ! (airport/POI)
        do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"     Fill PT variables"
        endif;enddo

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
          write(linebuffer130,'(130a)')Airport_Code(i)
          strlen = len(trim(adjustl(Airport_Code(i))))
          do j=1,strlen
            nSTAT=nf90_put_var(ncid,pt_code_var_id,linebuffer130(j:j),(/j,i/))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pt_code")
          enddo
        enddo

        ! Point/Airport name
        do i=1,nairports
          write(linebuffer130,'(130a)')Airport_Name(i)
          strlen = len(adjustl(Airport_Name(i))) ! This should be 130
          do j=1,strlen
            nSTAT=nf90_put_var(ncid,pt_name_var_id,linebuffer130(j:j),(/j,i/))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pt_name")
          enddo
        enddo

        ! Time-series of ash accumulation at points
        ! Since we are creating a new file, we will use step 1 of Airport_Thickness_TS
        allocate(dum2dPOI_out(nairports,1))
        dum2dPOI_out(1:nairports,1) = real(Airport_Thickness_TS(1:nairports,1),kind=op)
        nSTAT=nf90_put_var(ncid,pt_ashthickness_var_id,dum2dPOI_out,(/1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pt_depothick")
        ! And copy this step (the first) to the variable for the final deposit
        nSTAT=nf90_put_var(ncid,pt_ashthicknessFin_var_id,dum2dPOI_out(:,1),(/1/))
        deallocate(dum2dPOI_out)

      endif ! Write_PT_Data

      if(Write_PR_Data)then
        ! Profile data
        ! x coordinate of point
        allocate(dum1d_out(nvprofiles))
        dum1d_out(1:nvprofiles) = real(x_vprofile(1:nvprofiles),kind=op)
        nSTAT=nf90_put_var(ncid,pr_x_var_id,dum1d_out,(/1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pr_x")
        deallocate(dum1d_out)

        ! y coordinate of point
        allocate(dum1d_out(nvprofiles))
        dum1d_out(1:nvprofiles) = real(y_vprofile(1:nvprofiles),kind=op)
        nSTAT=nf90_put_var(ncid,pr_y_var_id,dum1d_out,(/1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pr_y")
        deallocate(dum1d_out)

        ! profile name
        do i=1,nvprofiles
          write(linebuffer050,'(50a)')Site_vprofile(i)
          strlen = len(adjustl(Site_vprofile(i))) ! This should be 50
          do j=1,strlen
            nSTAT=nf90_put_var(ncid,pr_name_var_id,linebuffer050(j:j),(/j,i/))
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pr_name")
          enddo
        enddo
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(USE_OPTMOD_VARS)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     Fill OPTMOD_VARS"
        endif;enddo

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

      deallocate(dum2dint_out,dum2d_out,dum3d_out)
      deallocate(ashcon)
      deallocate(depocon)

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Created netcdf file"
      endif;enddo

      end subroutine NC_create_netcdf_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    NC_append_to_netcdf
!
!  Called from: output_results
!  Arguments:
!    none
!
!  This subroutine is called each time output_results is called when a new output
!  step is reached, writing the next step in the transient variables.  Additionally,
!  this subroutine is called at the end of the simulation after the time integration.
!  Variables for the Airport/POI data are also written each time this is called.
!  Several additional 'Final' variables are written at this point such as the
!  final deposit thickness (depothickFin) and the deposit/ash-cloud arrival
!  times (depotime,ash_arrival_time). Airport/POI variables written in this final
!  step are the arrival and duration for both deposit and ash cloud.  If profile
!  data is tracked, the whole transient variable (time_native) is written in this
!  final step.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine NC_append_to_netcdf

      use global_param,  only : &
         EPS_SMALL,KM2_2_M2

      use io_data,       only : &
         iout3d,nvar_User2d_XY,nvar_User3d_XYGs,nvar_User3d_XYZ,nvar_User4d_XYZGs,&
         concenfile,isFinal_TS,nvprofiles,Write_PT_Data,Write_PR_Data

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,ts1,dz_vec_pd

      use solution,      only : &
         vx_pd,vy_pd,vz_pd,vf_pd,concen_pd,DepositGranularity

      use time_data,     only : &
         time,ntmax,time_native,BaseYear

      use Output_Vars,   only : &
         dbZCol_FillValue,MaxConcentration_FillValue,&
         MaxHeight_FillValue,MinHeight_FillValue,Mask_Cloud,Mask_Deposit,&
         CloudArrivalTime_FillValue,DepArrivalTime_FillValue,&
         var_User2d_XY_name,var_User2d_XY,var_User3d_XYGs_name,var_User3d_XYGs,&
         var_User3d_XYZ_name,var_User3d_XYZ,var_User4d_XYZGs_name,var_User4d_XYZGs,&
         DBZ_THRESH,USE_OPTMOD_VARS,useRestartVars,&
         useOutprodVars,useWindVars,DepositThickness,DepArrivalTime,CloudArrivalTime,&
         MaxConcentration,MaxHeight,CloudLoad,dbZ,MinHeight,pr_ash,&
           dbZCalculator

      use Tephra,        only : &
         n_gs_max

      use Airports,      only : &
         nairports,Airport_Thickness_TS,Airport_AshDuration,Airport_CloudArrivalTime,&
         Airport_CloudDuration,Airport_AshArrivalTime,Airport_thickness

      integer :: i,j,k,isize

      integer :: nSTAT
      integer :: ivar
      character (len=16)         :: outstring

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine NC_append_to_netcdf"
      endif;enddo

      allocate(ashcon(nxmax,nymax,nzmax,nsmax))
      allocate(depocon(nxmax,nymax,nsmax))
      if(nsmax.gt.n_gs_max)then
        ! These are some non-ash species
        ns_extra = nsmax-n_gs_max
      else
        ns_extra = 0
      endif

      ! Open netcdf file for writing
      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"opening netcdf file"
      endif;enddo
      nSTAT=nf90_open(concenfile,nf90_write, ncid)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"nf90_open")

      if(isFinal_TS)then

        ! If this is the final time-step, then enter define mode and add the
        ! native time dimension and associated variables
        nSTAT = nf90_redef(ncid)
        ! tn (time native)
        if(NCversion.eq.4)then
#ifndef NC3
          nSTAT = nf90_def_dim(ncid,"tn",nf90_unlimited,tn_dim_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_dim tn:")
#endif
        else
          nSTAT = nf90_def_dim(ncid,"tn",ntmax,tn_dim_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_dim tn:")
        endif
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"     TN: ","tn"
        endif;enddo
        nSTAT = nf90_def_var(ncid,"tn",&
                             nf90_double,&
                             (/tn_dim_id/),&
                             tn_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"def_var tn")
        nSTAT = nf90_put_att(ncid,tn_var_id,"long_name","Time native")
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att tn long_name")
        write(outstring,4501)BaseYear
 4501   format('hours since ',i4)
        nSTAT = nf90_put_att(ncid,tn_var_id,"units",outstring)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att tn units")

        if (Write_PR_Data)then
          ! Profile data
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
          nSTAT = nf90_put_att(ncid,pr_ash_var_id,"long_name","Ash Concentration (Max)")
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pr_ash long_name")
          nSTAT = nf90_put_att(ncid,pr_ash_var_id,"units","mg/m3")
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pr_ash units")
          nSTAT = nf90_put_att(ncid,pr_ash_var_id,&
                   "missing_value", MaxConcentration_FillValue)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pr_ash missing_value")
          nSTAT = nf90_put_att(ncid,pr_ash_var_id,"_FillValue",MaxConcentration_FillValue)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_att pr_ash _FillValue")

        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                    Leaving define mode.                               !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        nSTAT = nf90_enddef(ncid)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"enddef")
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Leaving define mode"
        endif;enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! depothickFin
        allocate(dum2d_out(nxmax,nymax))
        nSTAT = nf90_inq_varid(ncid,"depothickFin",depothickFin_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid depothickFin")
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"  Writing depothickFin"
        endif;enddo
        dum2d_out(:,:) = real(DepositThickness,kind=op)
        nSTAT=nf90_put_var(ncid,depothickFin_var_id,dum2d_out,(/1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var depothickFin")
        deallocate(dum2d_out)

          ! depotime
        allocate(dum2d_out(nxmax,nymax))
        nSTAT = nf90_inq_varid(ncid,"depotime",depotime_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid depotime")
        do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"  Writing depotime"
        endif;enddo
        dum2d_out(:,:) = DepArrivalTime_FillValue
        dum2d_out = merge(real(DepArrivalTime,kind=op),dum2d_out,Mask_Deposit)
        nSTAT=nf90_put_var(ncid,depotime_var_id,dum2d_out,(/1,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var depotime")
        deallocate(dum2d_out)

          ! ashtime
        allocate(dum2d_out(nxmax,nymax))
        nSTAT = nf90_inq_varid(ncid,"ash_arrival_time",ashcloudtime_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid ash_arrival_time")
        do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"  Writing ashtime"
        endif;enddo
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

        if (Write_PT_Data)then
          allocate(dum1d_out(nairports))

          ! Arrival time of ashfall
          nSTAT = nf90_inq_varid(ncid,"pt_depotime",pt_asharrival_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid pt_depotime")
          dum1d_out(1:nairports) = real(Airport_AshArrivalTime(1:nairports),kind=op)
          nSTAT=nf90_put_var(ncid,pt_asharrival_var_id,dum1d_out,(/1/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pt_depotime")

          ! Duration of ashfall
          nSTAT = nf90_inq_varid(ncid,"pt_depodur",pt_ashduration_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid pt_depodur")
          dum1d_out(1:nairports) = real(Airport_AshDuration(1:nairports),kind=op)
          nSTAT=nf90_put_var(ncid,pt_ashduration_var_id,dum1d_out,(/1/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pt_depodur")

          ! Arrival time of ash cloud
          nSTAT = nf90_inq_varid(ncid,"pt_cloud_arrival",pt_cloudarrival_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid pt_cloud_arrival")
          dum1d_out(1:nairports) = real(Airport_CloudArrivalTime(1:nairports),kind=op)
          nSTAT=nf90_put_var(ncid,pt_cloudarrival_var_id,dum1d_out,(/1/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pt_cloud_arrival")

          ! Duration of ash cloud
          nSTAT = nf90_inq_varid(ncid,"pt_cloud_dur",pt_cloudduration_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid pt_cloud_dur")
          dum1d_out(1:nairports) = real(Airport_CloudDuration(1:nairports),kind=op)
          nSTAT=nf90_put_var(ncid,pt_cloudduration_var_id,dum1d_out,(/1/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pt_cloud_dur")

          ! Final deposit thickness
          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"Trying to find pt_depothickFin"
          endif;enddo
          nSTAT = nf90_inq_varid(ncid,"pt_depothickFin",pt_ashthicknessFin_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid pt_depothickFin")
          dum1d_out(1:nairports) = real(Airport_thickness(1:nairports),kind=op)
          nSTAT=nf90_put_var(ncid,pt_ashthicknessFin_var_id,dum1d_out,(/1/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pt_depothickFin")

          deallocate(dum1d_out)
        endif

        ! write out native time
        nSTAT = nf90_inq_varid(ncid,"tn",tn_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid tn")

        allocate(dum1d_out_dp(ntmax))
        dum1d_out_dp(1:ntmax) = real(time_native(1:ntmax),kind=dp)
        nSTAT=nf90_put_var(ncid,tn_var_id,dum1d_out_dp,(/1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var tn")
        deallocate(dum1d_out_dp)

        if(Write_PR_Data)then
          ! Profile data
          ! write out profiles
          nSTAT = nf90_inq_varid(ncid,"pr_ash",pr_ash_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid pr_ash")
          allocate(dum3d_out(nzmax,ntmax,nvprofiles))
          dum3d_out(1:nzmax,1:ntmax,1:nvprofiles) = &
                    real(pr_ash(1:nzmax,1:ntmax,1:nvprofiles),kind=op)
          nSTAT=nf90_put_var(ncid,pr_ash_var_id,dum3d_out,(/1,1,1/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pr_ash")
          deallocate(dum3d_out)
          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"WROTE :",nvprofiles,nzmax,ntmax
          endif;enddo
        endif

      endif

      if(.not.isFinal_TS)then
        ! Get variable ids
        nSTAT = nf90_inq_varid(ncid,"t",t_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid t")
  
        do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"Got var IDs, now writing data"
        endif;enddo
        ! Write data
        ! Time
        do io=1,2;if(VB(io).le.verbosity_debug1)then
          write(outlog(io),*)"  Writing Time"
        endif;enddo
        dumscal_out = real(time,kind=op)
        nSTAT=nf90_put_var(ncid,t_var_id,dumscal_out,(/iout3d/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var t")

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(useWindVars)then
          allocate(dum3d_out(nxmax,nymax,nzmax))
            ! Vz
          nSTAT = nf90_inq_varid(ncid,"vz",vz_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid vz")
          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"  Writing Vz"
          endif;enddo
          dum3d_out(1:nxmax,1:nymax,1:nzmax) = real(vz_pd(1:nxmax,1:nymax,1:nzmax),kind=op)
          nSTAT=nf90_put_var(ncid,vz_var_id,dum3d_out,(/1,1,1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var vz")
            ! Vy
          nSTAT = nf90_inq_varid(ncid,"vy",vy_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid vy")
          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"  Writing Vy"
          endif;enddo
          dum3d_out(1:nxmax,1:nymax,1:nzmax) = real(vy_pd(1:nxmax,1:nymax,1:nzmax),kind=op)
          nSTAT=nf90_put_var(ncid,vy_var_id,dum3d_out,(/1,1,1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var vy")
            ! Vx
          nSTAT = nf90_inq_varid(ncid,"vx",vx_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid vx")
          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"  Writing Vx"
          endif;enddo
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
        if(useRestartVars)then
            ! ashcon
            ! netCDF standard requires that the unlimited dimension (time)
            ! be the most slowly varying, so ashcon unfortunately has a
            ! different shape than concen
          nSTAT = nf90_inq_varid(ncid,"ashcon",ashcon_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid ashcon")
          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"  Writing ashcon"
          endif;enddo
          ashcon = 0.0_op
          ashcon(1:nxmax,1:nymax,1:nzmax,1:nsmax) = real(concen_pd(1:nxmax,1:nymax,1:nzmax,1:nsmax,ts1),kind=op)
          nSTAT=nf90_put_var(ncid,ashcon_var_id,ashcon,(/1,1,1,1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var ashcon")
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(useOutprodVars)then

            ! depocon
          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"     Writing depocon"
          endif;enddo
          nSTAT = nf90_inq_varid(ncid,"depocon",depocon_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid depocon")
          depocon = 0.0_op
            ! netCDF standard requires that the unlimited dimension (time)
            ! be the most slowly varying, so dum3d_out unfortunately has a
            ! different shape than depocon
          do isize=1,n_gs_max
              ! Here's the conversion to kg/m^2 from kg/km^2
            depocon(1:nxmax,1:nymax,isize) = real(DepositGranularity(1:nxmax,1:nymax,isize) * &
                                           dz_vec_pd(0)/KM2_2_M2,kind=op)
          enddo
          do isize=1,n_gs_max
            do i=1,nxmax
              do j=1,nymax
                if(depocon(i,j,isize).le.EPS_SMALL)depocon(i,j,isize)=0.0_op
              enddo
            enddo
          enddo
          nSTAT=nf90_put_var(ncid,depocon_var_id,depocon,(/1,1,1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var depocon")

          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"Calling dbZCalculator"
          endif;enddo
          call dbZCalculator            ! get radar reflectivity

          ! depothick
          allocate(dum2d_out(nxmax,nymax))
          nSTAT = nf90_inq_varid(ncid,"depothick",depothick_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid depothick")
          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"  Writing depothick"
          endif;enddo
          dum2d_out(:,:) = real(DepositThickness,kind=op)
          nSTAT=nf90_put_var(ncid,depothick_var_id,dum2d_out,(/1,1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var depothick")
          deallocate(dum2d_out)

          ! ashconMax
          allocate(dum2d_out(nxmax,nymax))
          nSTAT = nf90_inq_varid(ncid,"ashcon_max",ashconMax_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid ashcon_max")
          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"  Writing ashconMax"
          endif;enddo
          dum2d_out(:,:) = real(MaxConcentration,kind=op)
          nSTAT=nf90_put_var(ncid,ashconMax_var_id,dum2d_out,(/1,1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var ashcon_max")
          deallocate(dum2d_out)

          ! ash cloud_height
          allocate(dum2d_out(nxmax,nymax))
          nSTAT = nf90_inq_varid(ncid,"cloud_height",ashheight_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid cloud_height")
          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"  Writing ash height"
          endif;enddo
          dum2d_out(:,:) = real(MaxHeight_FillValue,kind=op)
          dum2d_out = merge(real(MaxHeight,kind=op),dum2d_out,Mask_Cloud)
          nSTAT=nf90_put_var(ncid,ashheight_var_id,dum2d_out,(/1,1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var cloud_height")
          deallocate(dum2d_out)

          ! ash-load
          allocate(dum2d_out(nxmax,nymax))
          nSTAT = nf90_inq_varid(ncid,"cloud_load",ashload_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid cloud_load")
          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"  Writing ash-load"
          endif;enddo
          dum2d_out(:,:) = real(CloudLoad,kind=op)
          nSTAT=nf90_put_var(ncid,ashload_var_id,dum2d_out,(/1,1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var cloud_load")
          deallocate(dum2d_out)

          ! cloud-mask
          allocate(dum2dint_out(nxmax,nymax))
          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"     Fill cloud mask"
          endif;enddo
          dum2dint_out = merge(1,0,Mask_Cloud)
          nSTAT=nf90_put_var(ncid,cloudmask_var_id,dum2dint_out,(/1,1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var cloud_mask")
          deallocate(dum2dint_out)

          ! radar reflectivity
          allocate(dum3d_out(nxmax,nymax,nzmax))
          nSTAT = nf90_inq_varid(ncid,"radar_reflectivity",radrefl_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid radar_reflectivity")
          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"  Writing radar reflectivity"
          endif;enddo
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
          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"  Writing ash-height (bottom)"
          endif;enddo
          dum2d_out(:,:) = real(MinHeight_FillValue,kind=op)
          dum2d_out = merge(real(MinHeight,kind=op),dum2d_out,Mask_Cloud)
          nSTAT=nf90_put_var(ncid,ashcloudBot_var_id,dum2d_out,(/1,1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var cloud_bottom")
          deallocate(dum2d_out)

        endif ! useOutprodVars

        if (Write_PT_Data)then
          ! Time-series of ash accumulation at points
          nSTAT = nf90_inq_varid(ncid,"pt_depothick",pt_ashthickness_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid pt_depothick")
          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"  Writing time-series of ashfall at points"
          endif;enddo
          allocate(dum2d_out(nairports,1))
          dum2d_out(1:nairports,1) = real(Airport_Thickness_TS(1:nairports,iout3d),kind=op)
          nSTAT=nf90_put_var(ncid,pt_ashthickness_var_id,dum2d_out,(/1,iout3d/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pt_depothick")
          deallocate(dum2d_out)

          ! Final deposit thickness
          nSTAT = nf90_inq_varid(ncid,"pt_depothickFin",pt_ashthickness_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid pt_depothickFin")
          allocate(dum1d_out(nairports))
          dum1d_out(1:nairports) = real(Airport_thickness(1:nairports),kind=op)
          nSTAT=nf90_put_var(ncid,pt_ashthickness_var_id,dum1d_out,(/1/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var pt_depothickFin")
          deallocate(dum1d_out)

        endif
      endif ! Write_PT_Data
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(USE_OPTMOD_VARS)then
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
      endif ! .not.isFinal_TS
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Close file
      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"closing file"
      endif;enddo
      nSTAT = nf90_close(ncid)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"nf90_close")

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"Deallocating"
      endif;enddo
      deallocate(ashcon)
      deallocate(depocon)

      end subroutine NC_append_to_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    NC_RestartFile_ReadTimes
!
!  Called from: Read_Control_File
!  Arguments:
!    none
!
!  This subroutine is called if a restart run is requested.  The concentration
!  file is opened and the time variable read.  The user if prompted for the
!  start step (init_tstep).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine NC_RestartFile_ReadTimes

      use io_data,           only : &
         concenfile,init_tstep

      use time_data,         only : &
         time

      integer :: nSTAT
      integer :: ncid

      integer :: t_dim_id
      integer :: t_var_id
      integer :: t_len
      integer :: it
      real(kind=op), allocatable, dimension(:) :: t_list
      logical           :: IsThere
      integer           :: iostatus
      character(len=120):: iomessage
      character(len= 50):: linebuffer050 
      character(len= 80):: linebuffer080

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine NC_RestartFile_ReadTimes"
      endif;enddo

      ! Since we haven't opened a logfile yet, only write out to stdout if not a
      ! control file case.
      io = 1

      inquire(file=trim(adjustl(concenfile)),exist=IsThere)
      if(.not.IsThere)then
        write(outlog(io),*)"Concentraion file not found."
        stop 1
      endif

      ! Open netcdf file for reading
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

      write(outlog(io),*)"  Step :  time"
      do it = 1,t_len
        write(outlog(io),*)it,t_list(it)
      enddo

      write(outlog(io),*)'Enter timestep for initialization'
      read(input_unit,*,iostat=iostatus,iomsg=iomessage) linebuffer080
      linebuffer050 = "Reading init_tstep from stdin"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) init_tstep
      linebuffer050 = "Reading init_tstep from linebuffer"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)

      nSTAT=nf90_get_var(ncid,t_var_id,dumscal_out,(/init_tstep/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var t")
      time = real(dumscal_out,kind=ip)

      deallocate(t_list)

      end subroutine NC_RestartFile_ReadTimes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    NC_RestartFile_LoadConcen
!
!  Called from: Ash3d.F90
!  Arguments:
!    none
!
!  This subroutine reads the concentration and deposit arrays from the restart
!  netcdf file for the time step requested in NC_RestartFile_ReadTimes and
!  copies these data to concen_pd and outflow_xy1_pd (mirrored to DepositGranularity)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine NC_RestartFile_LoadConcen

      use io_data,           only : &
         concenfile,init_tstep

      use mesh,              only : &
         nxmax,nymax,nzmax,nsmax,ts0,ts1

      use solution,          only : &
         concen_pd,outflow_xy1_pd,DepositGranularity

      integer :: i
      integer :: nSTAT

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine NC)_RestartFile_LoadConcen"
      endif;enddo

      allocate(dum2d_out(nxmax,nymax))
      allocate(dum3d_out(nxmax,nymax,nzmax))
      allocate(ashcon(nxmax,nymax,nzmax,nsmax))
      allocate(depocon(nxmax,nymax,nsmax))

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"WARNING "
        write(outlog(io),*)"Input file is not currently verified "
        write(outlog(io),*)" with previous run."
      endif;enddo

      ! Open netcdf file for reading
      nSTAT=nf90_open(concenfile,nf90_nowrite,ncid)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"nf90_open")

      ! Get variable ids
      nSTAT = nf90_inq_varid(ncid,"t",t_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid t")
      nSTAT = nf90_inq_varid(ncid,"ashcon",ashcon_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid ashcon")
      nSTAT = nf90_inq_varid(ncid,"depocon",depocon_var_id)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid depocon")

      nSTAT=nf90_get_var(ncid,ashcon_var_id,ashcon,  &
               start = (/1,1,1,1,init_tstep/),       &
               count = (/nxmax,nymax,nzmax,nsmax,1/))

      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var ashcon")

      nSTAT=nf90_get_var(ncid,depocon_var_id,depocon,&
               start = (/1,1,1,init_tstep/),       &
               count = (/nxmax,nymax,nsmax,1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"put_var depocon")

      do i = 1,nsmax
        concen_pd(1:nxmax,1:nymax,1:nzmax,i,ts1) = real(ashcon(:,:,:,i),kind=ip)
        outflow_xy1_pd(1:nxmax,1:nymax,i) = real(depocon(:,:,i),kind=ip)
        DepositGranularity(1:nxmax,1:nymax,i) = outflow_xy1_pd(1:nxmax,1:nymax,i) 
      enddo

      concen_pd(:,:,:,:,ts0) = concen_pd(:,:,:,:,ts1)

      ! Close file
      nSTAT = nf90_close(ncid)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"nf90_close")

      deallocate(dum2d_out,dum3d_out)
      deallocate(ashcon)
      deallocate(depocon)

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Read concentrations from time ",dumscal_out
      endif;enddo

      end subroutine NC_RestartFile_LoadConcen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    NC_check_status
!
!  Called from: within the Ash3d_Netcdf module after each call to the ncf90 interface
!  Arguments:
!    nSTAT     = error code returned from netcdf call
!    icode     = user-supplied return value on stopping of code (0 just issues warning)
!    operation = string descriptor of function call causing error
!
!  Error-checking routine for NetCDF function calls. Modeled after a subroutine
!  posted at:
!    https://climate-cms.org/2018/10/12/create-netcdf.html
!  This subroutine provides a gracefull handling of errors returned from calls
!  to the nf90 interface.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine NC_check_status(nSTAT, icode, operation)

      use io_units

      integer, intent(in) :: nSTAT
      integer, intent(in) :: icode
      character(len=*), intent(in) :: operation

      character(len=9) :: severity

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine NC_check_status"
      endif;enddo

      if (icode.eq.0)then 
        severity = "WARNING: "
       else
        severity = "ERROR:   "
      endif

      if (nSTAT.eq.nf90_noerr) return
      do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io) ,*)severity,icode,operation,":",nf90_strerror(nSTAT)
      endif;enddo

      ! If user-supplied error code is 0, then consider this a warning,
      ! otherwise do a hard stop
      !if (icode.ne.0) stop icode
      if (icode.ne.0) stop 1

      end subroutine NC_check_status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    NC_Read_Output_Products
!
!  Called from: Ash3d_PostProc.f90
!  Arguments:
!    timestep
!
!  This subroutine is used in the post-processing tool to populate the variables
!  for a requested time step.  All the standard output products are read, but
!  only for timestep. If timestep = -1, then the last step in file is used
!  unless there is a 'final' variable as in depothickFin.  All other variables
!  needed for post-processing are also read, such as the grid coordinates,
!  grainsize and eruption variables as well as many global attributes.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine NC_Read_Output_Products(timestep)

      use global_param,  only : &
         GRAV,CFL,DT_MIN,DT_MAX,RAD_EARTH

      use io_data,           only : &
         concenfile,init_tstep,nWriteTimes,WriteTimes,cdf_b1l1,cdf_b1l5,cdf_b3l1, &
         cdf_b1l2,cdf_b3l3,VolcanoName,Write_PT_Data,isFinal_TS,&
         cdf_run_class,cdf_url,cdf_institution,&
         Write_PR_Data,nvprofiles,x_vprofile,y_vprofile,Site_vprofile

      use mesh,          only : &
         nxmax,nymax,nsmax,nzmax,x_cc_pd,y_cc_pd,lon_cc_pd,lat_cc_pd,z_cc_pd, &
         dx,dy,de,dn,dz_const,IsLatLon,latLL,lonLL,latUR,lonUR,xLL,yLL,xUR,yUR,&
         A3d_iprojflag,A3d_k0_scale,A3d_phi0,A3d_lam0,A3d_lam1,A3d_phi1,A3d_lam2,&
         A3d_phi2,A3d_Re,ZPADDING

      use solution,      only : &
         SpeciesID

      use time_data,     only : &
          BaseYear,useLeap,os_time_log,time,time_native,SimStartHour,xmlSimStartTime, &
          Simtime_in_hours,ntmax

      use Source,        only : &
         neruptions,e_Volume,e_Duration,e_StartTime,e_PlumeHeight, &
         lat_volcano,lon_volcano,x_volcano,y_volcano

      use Output_Vars,   only : &
         DepositThickness,DepArrivalTime,CloudArrivalTime,pr_ash,&
         MaxConcentration,MaxHeight,CloudLoad,dbZCol,MinHeight,Mask_Cloud,&
         CLOUDCON_GRID_THRESH,CLOUDCON_THRESH,THICKNESS_THRESH, &
         CLOUDLOAD_THRESH,DBZ_THRESH,DEPO_THRESH,DEPRATE_THRESH,ashcon_tot, &
         useRestartVars, &
           dbZCalculator, &
           Allocate_NTime, &
           Allocate_Profile, &
           Set_OutVar_ContourLevel

      use Airports,      only : &
         nairports,Airport_Code,Airport_Name,&
         Airport_Latitude,Airport_Longitude,Airport_Thickness_TS,Airport_thickness,&
         Airport_AshArrived,Airport_CloudArrivalTime,&
         Airport_CloudDuration,Airport_AshArrivalTime,Airport_AshDuration,&
         Airport_Thickness,Airport_CloudArrived,Airport_i,Airport_j,&
           Allocate_Airports

      use Tephra,        only : &
         n_gs_max,MagmaDensity,DepositDensity,LAM_GS_THRESH,AIRBORNE_THRESH

      use projection,    only : &
         PJ_iprojflag,PJ_k0,PJ_lam0,PJ_lam1,PJ_lam2,PJ_phi0,PJ_phi1,PJ_phi2,PJ_Re,&
           PJ_Set_Proj_Params,PJ_proj_for,PJ_proj_inv

      integer, intent(in), optional :: timestep

      logical,save :: first_time = .true.
      integer :: nSTAT
      integer :: it,i,j,isize
      integer :: var_xtype
      integer :: fop  ! output precision used for Ash3d output file being read in (4 or 8)
      character(len=NF90_MAX_NAME)  :: invar
      real(kind=op) :: dumscal_out
      real(kind=sp), dimension(:),allocatable :: dum1d_sp
      real(kind=dp), dimension(:),allocatable :: dum1d_dp
      character(len=32) :: time_units
      integer           :: iendstr
      integer           :: iostatus
      character(len=120):: iomessage
      character(len= 50):: linebuffer050 
      character(len= 80):: linebuffer080

      integer :: itstart_year,itstart_month,itstart_day
      integer :: itstart_hour,itstart_min,itstart_sec
      real(kind=ip) :: filestart_hour
      integer :: tmp_int
      real(kind=ip) :: lat_in,lon_in
      real(kind=ip) :: xnow,ynow,xout,yout

      INTERFACE
        real(kind=8) function HS_hours_since_baseyear(iyear,imonth,iday,hours,byear,useLeaps)
          integer            :: iyear
          integer            :: imonth
          integer            :: iday
          real(kind=8)       :: hours
          integer            :: byear
          logical            :: useLeaps
        end function HS_hours_since_baseyear
      END INTERFACE

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine NC_Read_Output_Products"
      endif;enddo

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Reading NetCDF file ",concenfile
      endif;enddo

      ! Open netcdf file for reading
      nSTAT=nf90_open(concenfile,nf90_nowrite,ncid)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"nSTAT=nf90_open")

      if(first_time)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Found the following dimensions and sizes"
        endif;enddo

        !!!!!!  TIME  !!!!!!!!!!!
        ! Identify dimension for time (and note size)
        nSTAT = nf90_inq_dimid(ncid,"t",t_dim_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_dimid time")
        nSTAT = nf90_Inquire_Dimension(ncid,t_dim_id,len=t_len)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"Inquire_Dimension time")
        ! Get variable id for this dimension
        nSTAT = nf90_inq_varid(ncid,"t",t_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid t")
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),2501)"  t",t_len
        endif;enddo
 2501 format(6x,a3,5x,i5) 

        !!!!!!  X  !!!!!!!!!!!
        ! Identify dimension for x (and note size)
        nSTAT = nf90_inq_dimid(ncid,"x",x_dim_id)
        !call NC_check_status(nSTAT,0,"inq_dimid x")
        if(nSTAT.ne.0)then
          ! dimension x not found, trying for lon
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
          nSTAT = nf90_inquire_variable(ncid, x_var_id, invar, xtype = var_xtype)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inquire_variable x")
          if(var_xtype.eq.NF90_FLOAT)then
            fop = 4
          elseif(var_xtype.eq.NF90_DOUBLE)then
            fop = 8
          else
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"var type for x is: ",var_xtype
              write(errlog(io),*)"Expecting either NF90_FLOAT or NF90_DOUBLE"
            endif;enddo
            stop 1
          endif
#ifdef USEPOINTERS
          if(.not.associated(lon_cc_pd))then
#else
          if(.not.allocated(lon_cc_pd))then
#endif
            nxmax = x_len
            allocate(lon_cc_pd(-1:nxmax+2))
            if(fop.eq.4)then
              allocate(dum1d_sp(1:nxmax))
              nSTAT=nf90_get_var(ncid,x_var_id,dum1d_sp,  &
                       start = (/1/),       &
                       count = (/x_len/))
              if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var x")
              lon_cc_pd(1:nxmax) = real(dum1d_sp(1:nxmax),kind=ip)
              deallocate(dum1d_sp)
            else
              allocate(dum1d_dp(1:nxmax))
              nSTAT=nf90_get_var(ncid,x_var_id,dum1d_dp,  &
                       start = (/1/),       &
                       count = (/x_len/))
              if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var x")
              lon_cc_pd(1:nxmax) = real(dum1d_dp(1:nxmax),kind=ip)
              deallocate(dum1d_dp)
            endif
          else
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)" lon_cc_pd already allocated"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),2501)"lon",x_len
          endif;enddo
        else
          ! Get projection information
          nSTAT = nf90_get_att(ncid,nf90_global,"b1l2",cdf_b1l2)
          if(nSTAT.ne.0)then
            call NC_check_status(nSTAT,0,"get_att b1l2:")
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"Did not find att b1l2: Projection parameters"
            endif;enddo
          endif
          call PJ_Set_Proj_Params(cdf_b1l2)
          A3d_iprojflag  = PJ_iprojflag
          A3d_k0_scale   = PJ_k0
          A3d_Re         = PJ_Re
          A3d_lam0       = PJ_lam0
          A3d_lam1       = PJ_lam1
          A3d_lam2       = PJ_lam2
          A3d_phi0       = PJ_phi0
          A3d_phi1       = PJ_phi1
          A3d_phi2       = PJ_phi2

          ! Get variable id for this dimension
          nSTAT = nf90_inq_varid(ncid,"x",x_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid x")
          nSTAT = nf90_inquire_variable(ncid, x_var_id, invar, xtype = var_xtype)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inquire_variable x")
          if(var_xtype.eq.NF90_FLOAT)then
            fop = 4
          elseif(var_xtype.eq.NF90_DOUBLE)then
            fop = 8
          else
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"var type for x is: ",var_xtype
              write(errlog(io),*)"Expecting either NF90_FLOAT or NF90_DOUBLE"
            endif;enddo
            stop 1
          endif
          if(fop.ne.op)then
            ! Until we check all input variables, force a hard stop if we have
            ! a mis-match in precision
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: The default real for this file is NF90_DOUBLE"
              write(errlog(io),*)"       but this post-processing program expects op=NF90_FLOAT"
              write(errlog(io),*)"       Please recompile with op=8"
            endif;enddo
            stop 1
          endif
#ifdef USEPOINTERS
          if(.not.associated(x_cc_pd))then
#else
          if(.not.allocated(x_cc_pd))then
#endif
            nxmax = x_len
            allocate(x_cc_pd(-1:nxmax+2))
            if(fop.eq.4)then
              allocate(dum1d_sp(1:nxmax))
              nSTAT=nf90_get_var(ncid,x_var_id,dum1d_sp,  &
                       start = (/1/),       &
                       count = (/x_len/))
              if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var x")
              x_cc_pd(1:nxmax) = real(dum1d_sp(1:nxmax),kind=ip)
              deallocate(dum1d_sp)
            else
              allocate(dum1d_dp(1:nxmax))
              nSTAT=nf90_get_var(ncid,x_var_id,dum1d_dp,  &
                       start = (/1/),       &
                       count = (/x_len/))
              if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var x")
              x_cc_pd(1:nxmax) = real(dum1d_dp(1:nxmax),kind=ip)
              deallocate(dum1d_dp)
            endif
          else
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"x_cc_pd already allocated"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),2501)"  x",x_len
          endif;enddo
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
#ifdef USEPOINTERS
          if(.not.associated(lat_cc_pd))then
#else
          if(.not.allocated(lat_cc_pd))then
#endif
            nymax = y_len
            allocate(lat_cc_pd(-1:nymax+2))
            ! Assume we know fop already from x above
            if(fop.eq.4)then
              allocate(dum1d_sp(1:nymax))
              nSTAT=nf90_get_var(ncid,y_var_id,dum1d_sp,  &
                       start = (/1/),       &
                       count = (/y_len/))
              if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var y")
              lat_cc_pd(1:nymax) = real(dum1d_sp(1:nymax),kind=ip)
              deallocate(dum1d_sp)
            else
              allocate(dum1d_dp(1:nymax))
              nSTAT=nf90_get_var(ncid,y_var_id,dum1d_dp,  &
                       start = (/1/),       &
                       count = (/y_len/))
              if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var y")
              lat_cc_pd(1:nymax) = real(dum1d_dp(1:nymax),kind=ip)
              deallocate(dum1d_dp)
            endif
          else
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"lat_cc_pd already allocated"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),2501)"lat",y_len
          endif;enddo
        else
          nSTAT = nf90_inq_dimid(ncid,"y",y_dim_id)
          call NC_check_status(nSTAT,1,"inq_dimid x")
          nSTAT = nf90_Inquire_Dimension(ncid,y_dim_id,len=y_len)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"Inquire_Dimension y")
          ! Get variable id for this dimension
          nSTAT = nf90_inq_varid(ncid,"y",y_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid y")
#ifdef USEPOINTERS
          if(.not.associated(y_cc_pd))then
#else
          if(.not.allocated(y_cc_pd))then
#endif
            nymax = y_len
            allocate(y_cc_pd(-1:nymax+2))
            ! Assume we know fop already from x above
            if(fop.eq.4)then
              allocate(dum1d_sp(1:nymax))
              nSTAT=nf90_get_var(ncid,y_var_id,dum1d_sp,  &
                       start = (/1/),       &
                       count = (/y_len/))
              if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var y")
              y_cc_pd(1:nymax) = real(dum1d_sp(1:nymax),kind=ip)
              deallocate(dum1d_sp)
            else
              allocate(dum1d_dp(1:nymax))
              nSTAT=nf90_get_var(ncid,y_var_id,dum1d_dp,  &
                       start = (/1/),       &
                       count = (/y_len/))
              if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var y")
              y_cc_pd(1:nymax) = real(dum1d_dp(1:nymax),kind=ip)
              deallocate(dum1d_dp)
            endif
          else
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"lat_cc_pd already allocated"
            endif;enddo
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),2501)"  y",y_len
          endif;enddo
        endif

        if(IsLatLon)then
          de = lon_cc_pd(2)-lon_cc_pd(1)
          dn = lat_cc_pd(2)-lat_cc_pd(1)
          lon_cc_pd(0)  = lon_cc_pd(1) - de
          lon_cc_pd(-1) = lon_cc_pd(0) - de
          lon_cc_pd(nxmax+1) = lon_cc_pd(nxmax  ) + de
          lon_cc_pd(nxmax+2) = lon_cc_pd(nxmax+1) + de
          lat_cc_pd(0)  = lat_cc_pd(1) - dn
          lat_cc_pd(-1) = lat_cc_pd(0) - dn
          lat_cc_pd(nymax+1) = lat_cc_pd(nymax  ) + dn
          lat_cc_pd(nymax+2) = lat_cc_pd(nymax+1) + dn
          lonLL = lon_cc_pd(1)-0.5_ip*de
          latLL = lat_cc_pd(1)-0.5_ip*dn
          lonUR = lon_cc_pd(nxmax)+0.5_ip*de
          latUR = lat_cc_pd(nymax)+0.5_ip*dn
        else
          dx  = x_cc_pd(2)-x_cc_pd(1)
          dy  = y_cc_pd(2)-y_cc_pd(1)
          x_cc_pd(0)  = x_cc_pd(1) - dx
          x_cc_pd(-1) = x_cc_pd(0) - dx
          x_cc_pd(nxmax+1) = x_cc_pd(nxmax  ) + dx
          x_cc_pd(nxmax+2) = x_cc_pd(nxmax+1) + dx
          y_cc_pd(0)  = y_cc_pd(1) - dy
          y_cc_pd(-1) = y_cc_pd(0) - dy
          y_cc_pd(nymax+1) = y_cc_pd(nymax  ) + dy
          y_cc_pd(nymax+2) = y_cc_pd(nymax+1) + dy
          xLL = x_cc_pd(1)-0.5_ip*dx
          yLL = y_cc_pd(1)-0.5_ip*dy
          xUR = x_cc_pd(nxmax)+0.5_ip*dx
          yUR = y_cc_pd(nymax)+0.5_ip*dy
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
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),2501)"  z",z_len
        endif;enddo
        nzmax = z_len 
        allocate(z_cc_pd(-1:nzmax+2))
        ! Assume we know fop already from x above
        if(fop.eq.4)then
          allocate(dum1d_sp(1:nzmax))
          nSTAT=nf90_get_var(ncid,z_var_id,dum1d_sp,  &
                   start = (/1/),       &
                   count = (/z_len/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var z")
          z_cc_pd(1:nzmax) = real(dum1d_sp(1:nzmax),kind=ip)
          deallocate(dum1d_sp)
        else
          allocate(dum1d_dp(1:nzmax))
          nSTAT=nf90_get_var(ncid,z_var_id,dum1d_dp,  &
                   start = (/1/),       &
                   count = (/z_len/))
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var z")
          z_cc_pd(1:nzmax) = real(dum1d_dp(1:nzmax),kind=ip)
          deallocate(dum1d_dp)
        endif
        dz_const  = z_cc_pd(2)-z_cc_pd(1)
        z_cc_pd(0)  = z_cc_pd(1) - dz_const
        z_cc_pd(-1) = z_cc_pd(0) - dz_const
        z_cc_pd(nzmax+1) = z_cc_pd(nzmax  ) + dz_const
        z_cc_pd(nzmax+2) = z_cc_pd(nzmax+1) + dz_const

        !!!!!!  BN  !!!!!!!!!!!
        ! Identify dimension for bn (and note size)
        nSTAT = nf90_inq_dimid(ncid,"bn",bn_dim_id)
        call NC_check_status(nSTAT,1,"inq_dimid bn")
        nSTAT = nf90_Inquire_Dimension(ncid,bn_dim_id,len=bn_len)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"Inquire_Dimension bn")
        ! Get variable id for this dimension
        nSTAT = nf90_inq_varid(ncid,"bn",bn_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid bn")
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),2501) "bn",bn_len
        endif;enddo
        nsmax = bn_len

        !!!!!!  PT  !!!!!!!!!!!
        ! Identify dimension for pt (and note size)
        nSTAT = nf90_inq_dimid(ncid,"pt",pt_dim_id)
        if(nSTAT.ne.0)then
          Write_PT_Data = .false.
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find dim pt: Output file has no point data"
          endif;enddo
          nairports = 0
        else
          Write_PT_Data = .true.
          nSTAT = nf90_Inquire_Dimension(ncid,pt_dim_id,len=pt_len)
          call NC_check_status(nSTAT,1,"Inquire_Dimension pt")
          ! Get variable id for this dimension
          nSTAT = nf90_inq_varid(ncid,"pt",pt_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid pt")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),2501)" pt",pt_len
          endif;enddo
          nairports = pt_len
        endif

        !!!!!!  PR  !!!!!!!!!!!
        ! Identify dimension for pr (and note size)
        nSTAT = nf90_inq_dimid(ncid,"pr",pr_dim_id)
        if(nSTAT.ne.0)then
          ! Since we are only issuing a warning and not a hard stop, reset id
          pr_dim_id = 0
          Write_PR_Data = .false.
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find dim pr: Output file has no profile data"
          endif;enddo
        else
          Write_PR_Data = .true.
          nSTAT = nf90_Inquire_Dimension(ncid,pr_dim_id,len=pr_len)
          call NC_check_status(nSTAT,0,"Inquire_Dimension pr")
          nvprofiles = pr_len
          ! Get variable id for this dimension
          nSTAT = nf90_inq_varid(ncid,"pr",pr_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid pr")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),2501)" pr",nvprofiles
          endif;enddo
        endif

        !!!!!!  TN !!!!!!!!!!!
        ! Identify dimension for tn (and note size)
        nSTAT = nf90_inq_dimid(ncid,"tn",tn_dim_id)
        if(nSTAT.ne.0)then
          ! Since we are only issuing a warning and not a hard stop, reset id
          tn_dim_id = 0
          call NC_check_status(nSTAT,0,"inq_dimid tn")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find dim tn: Output file has no native time data (for profiles)"
          endif;enddo
          tn_len = 0
          ntmax  = 0
        else
          nSTAT = nf90_Inquire_Dimension(ncid,tn_dim_id,len=tn_len)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,0,"Inquire_Dimension tn")
          ! Get variable id for this dimension
          nSTAT = nf90_inq_varid(ncid,"tn",tn_var_id)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,0,"inq_varid tn")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),2501)" tn",tn_len
          endif;enddo
          ntmax = tn_len
          call Allocate_NTime(ntmax)
        endif

        !!!!!!  ER !!!!!!!!!!!
        ! Identify dimension for er (and note size)
        nSTAT = nf90_inq_dimid(ncid,"er",er_dim_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_dimid er")
        nSTAT = nf90_Inquire_Dimension(ncid,er_dim_id,len=neruptions)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"Inquire_Dimension er")
        ! Get variable id for this dimension
        nSTAT = nf90_inq_varid(ncid,"er",er_var_id)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inq_varid er")
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),2501)" er",neruptions
        endif;enddo

        ! Now get the expected global attributes
        nSTAT = nf90_get_att(ncid,nf90_global,"BaseYear",BaseYear)
        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_att BaseYear:")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find att BaseYear: Assuming BaseYear=1900"
          endif;enddo
          BaseYear = 1900
        endif

        nSTAT = nf90_get_att(ncid,nf90_global,"useLeap",tmp_int)
        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_att useLeap:")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find att useLeap: Assuming useLeap=T"
          endif;enddo
          tmp_int = 1
        endif
        if(tmp_int.eq.1)then
          useLeap = .true.
        else
          useLeap = .false.
        endif

        nSTAT = nf90_get_att(ncid,nf90_global,"MagmaDensity",MagmaDensity)
        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_att Comment MagmaDensity:")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find att MagmaDensity: Assuming MagmaDensity=2500.0"
          endif;enddo
          MagmaDensity = 2500.0_ip
        endif

        nSTAT = nf90_get_att(ncid,nf90_global,"DepositDensity",DepositDensity)
        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_att Comment DepositDensity:")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find att DepositDensity: Assuming DepositDensity=1000.0"
          endif;enddo
          DepositDensity = 1000.0_ip
        endif

        nSTAT = nf90_get_att(ncid,nf90_global,"LAM_GS_THRESH",LAM_GS_THRESH)
        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_att Comment LAM_GS_THRESH:")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find att LAM_GS_THRESH: Assuming LAM_GS_THRESH=250.0"
          endif;enddo
          LAM_GS_THRESH = 250.0_ip
        endif

        nSTAT = nf90_get_att(ncid,nf90_global,"AIRBORNE_THRESH",AIRBORNE_THRESH)
        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_att Comment AIRBORNE_THRESH:")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find att AIRBORNE_THRESH: Assuming AIRBORNE_THRESH=1.0e-3"
          endif;enddo
          AIRBORNE_THRESH = 1.0e-3_ip
        endif

        nSTAT = nf90_get_att(ncid,nf90_global,"GRAV",GRAV)
        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_att Comment GRAV:")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find att GRAV: Assuming GRAV=9.81"
          endif;enddo
          GRAV = 9.81_ip
        endif

        nSTAT = nf90_get_att(ncid,nf90_global,"RAD_EARTH",RAD_EARTH)
        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_att Comment RAD_EARTH:")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find att RAD_EARTH: Assuming RAD_EARTH=6371.229"
          endif;enddo
          RAD_EARTH = 6371.229_ip
        endif

        nSTAT = nf90_get_att(ncid,nf90_global,"CFL",CFL)
        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_att Comment CFL:")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find att CFL: Assuming CFL=0.8"
          endif;enddo
          CFL = 0.80_ip
        endif

        nSTAT = nf90_get_att(ncid,nf90_global,"DT_MIN",DT_MIN)
        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_att Comment DT_MIN:")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find att DT_MIN: Assuming DT_MIN=1.0e-5"
          endif;enddo
          DT_MIN = 1.0e-5_ip
        endif

        nSTAT = nf90_get_att(ncid,nf90_global,"DT_MAX",DT_MAX)
        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_att Comment DT_MAX:")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find att DT_MAX: Assuming DT_MAX=1.0"
          endif;enddo
          DT_MAX = 1.0e0_ip
        endif

        nSTAT = nf90_get_att(ncid,nf90_global,"ZPADDING",ZPADDING)
        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_att Comment ZPADDING:")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find att ZPADDING Assuming ZPADDING=1.3"
          endif;enddo
          ZPADDING = 1.3_ip
        endif

        nSTAT = nf90_get_att(ncid,nf90_global,"DEPO_THRESH",DEPO_THRESH)
        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_att Comment DEPO_THRESH:")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find att DEPO_THRESH: Assuming DEPO_THRESH=1.0e-2"
          endif;enddo
          DEPO_THRESH = 1.0e-2_ip
        endif

        nSTAT = nf90_get_att(ncid,nf90_global,"DEPRATE_THRESH",DEPRATE_THRESH)
        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_att Comment DEPRATE_THRESH:")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find att DEPRATE_THRESH: Assuming DEPRATE_THRESH=1.0e-2"
          endif;enddo
          DEPRATE_THRESH = 1.0e-2_ip
        endif

        nSTAT = nf90_get_att(ncid,nf90_global,"CLOUDCON_THRESH",CLOUDCON_THRESH)
        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_att Comment CLOUDCON_THRESH:")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find att CLOUDCON_THRESH: Assuming CLOUDCON_THRESH=1.0e-3"
          endif;enddo
          CLOUDCON_THRESH = 1.0e-3_ip
        endif

        nSTAT = nf90_get_att(ncid,nf90_global,"CLOUDCON_GRID_THRESH",CLOUDCON_GRID_THRESH)
        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_att Comment CLOUDCON_GRID_THRESH:")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find att CLOUDCON_GRID_THRESH: Assuming CLOUDCON_GRID_THRESH=1.0e-7"
          endif;enddo
          CLOUDCON_GRID_THRESH = 1.0e-7_ip
        endif

        nSTAT = nf90_get_att(ncid,nf90_global,"CLOUDLOAD_THRESH",CLOUDLOAD_THRESH)
        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_att Comment CLOUDLOAD_THRESH:")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find att CLOUDLOAD_THRESH: Assuming CLOUDLOAD_THRESH=2.0e-1"
          endif;enddo
          CLOUDLOAD_THRESH = 2.0e-1_ip
        endif

        nSTAT = nf90_get_att(ncid,nf90_global,"THICKNESS_THRESH",THICKNESS_THRESH)
        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_att Comment THICKNESS_THRESH:")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find att THICKNESS_THRESH: Assuming THICKNESS_THRESH=1.0e-2"
          endif;enddo
          THICKNESS_THRESH = 1.0e-2_ip
        endif

        nSTAT = nf90_get_att(ncid,nf90_global,"DBZ_THRESH",DBZ_THRESH)
        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_att Comment DBZ_THRESH:")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find att DBZ_THRESH: Assuming DBZ_THRESH=-2.0e+1"
          endif;enddo
          DBZ_THRESH = -2.0e+1_ip
        endif


        ! Now get all the other variable info:
        ! Species class variable is needed to verify that nsmax = n_gs_max
        nSTAT = nf90_inq_varid(ncid,"spec_class",spec_var_id)
        if(nSTAT.ne.0)then
          ! Since we are only issuing a warning and not a hard stop, reset id
          spec_var_id = 0
          call NC_check_status(nSTAT,0,"inq_varid spec_class")
        endif
        !!!! Time-series vars
        ! Vx
        nSTAT = nf90_inq_varid(ncid,"vx",vx_var_id)
        if(nSTAT.ne.0)then
          vx_var_id = 0
          !call NC_check_status(nSTAT,0,"inq_varid vx")
        endif
        ! Vy
        nSTAT = nf90_inq_varid(ncid,"vy",vy_var_id)
        if(nSTAT.ne.0)then
          vy_var_id = 0
          !call NC_check_status(nSTAT,0,"inq_varid vy")
        endif
        ! Vz
        nSTAT = nf90_inq_varid(ncid,"vz",vz_var_id)
        if(nSTAT.ne.0)then
          vz_var_id = 0
          !call NC_check_status(nSTAT,0,"inq_varid vz")
        endif
        ! ashcon
        nSTAT = nf90_inq_varid(ncid,"ashcon",ashcon_var_id)
        if(nSTAT.ne.0)then
          ashcon_var_id = 0
          call NC_check_status(nSTAT,0,"inq_varid ashcon")
          useRestartVars = .false.
        else
          useRestartVars = .true.
        endif
        ! depocon
        nSTAT = nf90_inq_varid(ncid,"depocon",depocon_var_id)
        if(nSTAT.ne.0)then
          depocon_var_id = 0
          call NC_check_status(nSTAT,0,"inq_varid depocon")
        endif
        ! ashcon_max
        nSTAT = nf90_inq_varid(ncid,"ashcon_max",ashconMax_var_id)
        if(nSTAT.ne.0)then
          ashconMax_var_id = 0
          call NC_check_status(nSTAT,0,"inq_varid ashconMax")
        endif
        ! cloud_height
        nSTAT = nf90_inq_varid(ncid,"cloud_height",ashheight_var_id)
        if(nSTAT.ne.0)then
          ashheight_var_id = 0
          call NC_check_status(nSTAT,0,"inq_varid cloud_height")
        endif
        ! cloud_load
        nSTAT = nf90_inq_varid(ncid,"cloud_load",ashload_var_id)
        if(nSTAT.ne.0)then
          ashload_var_id = 0
          call NC_check_status(nSTAT,0,"inq_varid cloud_load")
        endif
        ! cloud_mask
        nSTAT = nf90_inq_varid(ncid,"cloud_mask",cloudmask_var_id)
        if(nSTAT.ne.0)then
          cloudmask_var_id = 0
          call NC_check_status(nSTAT,0,"inq_varid cloud_mask")
        endif
        ! cloud_bottom
        nSTAT = nf90_inq_varid(ncid,"cloud_bottom",ashcloudBot_var_id)
        if(nSTAT.ne.0)then
          ashcloudBot_var_id = 0
          call NC_check_status(nSTAT,0,"inq_varid cloud_bottom")
        endif
        ! depothick
        nSTAT = nf90_inq_varid(ncid,"depothick",depothick_var_id)
        if(nSTAT.ne.0)then
          depothick_var_id = 0
          call NC_check_status(nSTAT,0,"inq_varid depothick")
        endif
        ! radar_reflectivity
        nSTAT = nf90_inq_varid(ncid,"radar_reflectivity",radrefl_var_id)
        if(nSTAT.ne.0)then
          radrefl_var_id = 0
          call NC_check_status(nSTAT,0,"inq_varid radar_reflectivity")
        endif

        !!!! Final vars (at simulation end)
        ! depothickFin
        nSTAT = nf90_inq_varid(ncid,"depothickFin",depothickFin_var_id)
        if(nSTAT.ne.0)then
          depothickFin_var_id = 0
          call NC_check_status(nSTAT,0,"inq_varid depothickFin") 
        endif
        ! depotime
        nSTAT = nf90_inq_varid(ncid,"depotime",depotime_var_id)
        if(nSTAT.ne.0)then
          depotime_var_id = 0
          call NC_check_status(nSTAT,0,"inq_varid depotime")
        endif
        ! ash_arrival_time
        nSTAT = nf90_inq_varid(ncid,"ash_arrival_time",ashcloudtime_var_id)
        if(nSTAT.ne.0)then
          ashcloudtime_var_id = 0
          call NC_check_status(nSTAT,0,"inq_varid ash_arrival_time")
        endif

        if(Write_PT_Data)then
          ! pt_x
          nSTAT = nf90_inq_varid(ncid,"pt_x",pt_x_var_id)
          if(nSTAT.ne.0)then
            pt_x_var_id = 0
            call NC_check_status(nSTAT,0,"inq_varid pt_x")
          endif
          ! pt_y
          nSTAT = nf90_inq_varid(ncid,"pt_y",pt_y_var_id)
          if(nSTAT.ne.0)then
            pt_y_var_id = 0
            call NC_check_status(nSTAT,0,"inq_varid pt_y")
          endif
          ! pt_code
          nSTAT = nf90_inq_varid(ncid,"pt_code",pt_code_var_id)
          if(nSTAT.ne.0)then
            pt_code_var_id = 0
            call NC_check_status(nSTAT,0,"inq_varid pt_code")
          endif
          ! pt_name
          nSTAT = nf90_inq_varid(ncid,"pt_name",pt_name_var_id)
          if(nSTAT.ne.0)then
            pt_name_var_id = 0
            call NC_check_status(nSTAT,0,"inq_varid pt_name")
          endif
          ! pt_depotime
          nSTAT = nf90_inq_varid(ncid,"pt_depotime",pt_asharrival_var_id)
          if(nSTAT.ne.0)then
            pt_asharrival_var_id = 0
            call NC_check_status(nSTAT,0,"inq_varid pt_depotime")
          endif
          ! pt_depodur
          nSTAT = nf90_inq_varid(ncid,"pt_depodur",pt_ashduration_var_id)
          if(nSTAT.ne.0)then
            pt_ashduration_var_id = 0
            call NC_check_status(nSTAT,0,"inq_varid pt_depod")
          endif
          ! pt_cloud_arrival
          nSTAT = nf90_inq_varid(ncid,"pt_cloud_arrival",pt_cloudarrival_var_id)
          if(nSTAT.ne.0)then
            pt_cloudarrival_var_id = 0
            call NC_check_status(nSTAT,0,"inq_varid pt_cloud_arrival")
          endif
          ! pt_cloud_dur
          nSTAT = nf90_inq_varid(ncid,"pt_cloud_dur",pt_cloudduration_var_id)
          if(nSTAT.ne.0)then
            pt_cloudduration_var_id = 0
            call NC_check_status(nSTAT,0,"inq_varid pt_cloud_dur")
          endif
          ! pt_depothick
          nSTAT = nf90_inq_varid(ncid,"pt_depothick",pt_ashthickness_var_id)
          if(nSTAT.ne.0)then
            pt_ashthickness_var_id = 0
            call NC_check_status(nSTAT,0,"inq_varid pt_depothick")
          endif
          ! pt_depothickFin
          nSTAT = nf90_inq_varid(ncid,"pt_depothickFin",pt_ashthicknessFin_var_id)
          if(nSTAT.ne.0)then
            pt_ashthicknessFin_var_id = 0
            call NC_check_status(nSTAT,0,"inq_varid pt_depothickFin")
          endif
        endif

        ! Load and check species class so we can set n_gs_max
#ifdef USEPOINTERS
        if(.not.associated(SpeciesID))then
#else
        if(.not.allocated(SpeciesID))then
#endif
          allocate(SpeciesID(nsmax))
        endif
        nSTAT = nf90_get_var(ncid,spec_var_id,SpeciesID)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var SpeciesID")
        n_gs_max = 0
        do isize=1,nsmax
          if(SpeciesID(isize).eq.1) n_gs_max = n_gs_max + 1
        enddo

        nWriteTimes = t_len
#ifdef USEPOINTERS
        if(.not.associated(WriteTimes))then
#else
        if(.not.allocated(WriteTimes))then
#endif
          allocate(WriteTimes(nWriteTimes))
        endif
        ! Time should always be written with dp, however some older files used float
        ! Double-check type
        nSTAT = nf90_inquire_variable(ncid, t_var_id, invar, xtype = var_xtype)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inquire_variable t")
        if(var_xtype.eq.NF90_FLOAT)then
          allocate(dum1d_sp(1:nWriteTimes))
          nSTAT = nf90_get_var(ncid,t_var_id,dum1d_sp)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var t")
          WriteTimes(1:nWriteTimes) = real(dum1d_sp(1:nWriteTimes),kind=ip)
          deallocate(dum1d_sp)
        elseif(var_xtype.eq.NF90_DOUBLE)then
          allocate(dum1d_dp(1:nWriteTimes))
          nSTAT = nf90_get_var(ncid,t_var_id,dum1d_dp)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var t")
          WriteTimes(1:nWriteTimes) = real(dum1d_dp(1:nWriteTimes),kind=ip)
          deallocate(dum1d_dp)
        endif
        nSTAT = nf90_get_att(ncid,t_var_id,"units",time_units)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_att units t:")
        read(time_units,4313,iostat=iostatus,iomsg=iomessage) xmlSimStartTime
        linebuffer080 = xmlSimStartTime
        linebuffer050 = "Reading xmlSimStartTime from time var attribute"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
  4313  format(12x,a20)
        read(xmlSimStartTime,4314,iostat=iostatus,iomsg=iomessage)&
                             itstart_year,itstart_month,itstart_day, &
                             itstart_hour,itstart_min,itstart_sec
        linebuffer050 = "Reading xmlSimStartTime from time var attribute"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        filestart_hour = real(itstart_hour,kind=sp) + &
                             real(itstart_min,kind=sp)/60.0_sp      + &
                             real(itstart_sec,kind=sp)/3600.0_sp
  4314  format(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2,1x)
        SimStartHour = real(HS_hours_since_baseyear(itstart_year,itstart_month, &
                               itstart_day,real(filestart_hour,kind=8),&
                               BaseYear,useLeap),kind=4)
        if(Write_PT_Data.and.&
           nairports.gt.0)then
          ! Allocate and fill Airport variables
          call Allocate_Airports(nairports,nWriteTimes)

          ! Read TS of deposit thickness at point locations
          allocate(dum2d_out(1:nairports,1:nWriteTimes))
          nSTAT = nf90_get_var(ncid,pt_ashthickness_var_id,dum2d_out)
          Airport_Thickness_TS(:,:) = real(dum2d_out(:,:),kind=ip)
          deallocate(dum2d_out)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var pt_ashthickness:")

          ! Read final deposit thickness at point locations
          allocate(dum1d_out(1:nairports))
          nSTAT = nf90_get_var(ncid,pt_ashthicknessFin_var_id,dum1d_out)
          Airport_Thickness(:) = real(dum1d_out(:),kind=ip)
          deallocate(dum1d_out)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var pt_depothickFin:")

          ! Read Airport code
          nSTAT = nf90_get_var(ncid,pt_code_var_id,Airport_Code)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var pt_code:")

          ! Read Airport names
          nSTAT = nf90_get_var(ncid,pt_name_var_id,Airport_Name)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var pt_name:")
          if(fop.eq.4)then
            allocate(dum1d_sp(1:nairports))
            nSTAT = nf90_get_var(ncid,pt_x_var_id,dum1d_sp)
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var pt_x:")
            Airport_Longitude(:) = real(dum1d_sp(1:nairports),kind=ip)
            nSTAT = nf90_get_var(ncid,pt_y_var_id,dum1d_sp)
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var pt_y:")
            Airport_Latitude(:) = real(dum1d_sp(1:nairports),kind=ip)
            nSTAT = nf90_get_var(ncid,pt_asharrival_var_id,dum1d_sp)
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var pt_asharrival:")
            Airport_AshArrivalTime(:) = real(dum1d_sp(1:nairports),kind=ip)
            nSTAT = nf90_get_var(ncid,pt_ashduration_var_id,dum1d_sp)
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var pt_depodur:")
            Airport_AshDuration(:) = real(dum1d_sp(1:nairports),kind=ip)
            nSTAT = nf90_get_var(ncid,pt_cloudarrival_var_id,dum1d_sp)
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var pt_cloudarrival:")
            Airport_CloudArrivalTime(:) = real(dum1d_sp(1:nairports),kind=ip)
            nSTAT = nf90_get_var(ncid,pt_cloudduration_var_id,dum1d_sp)
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var pt_cloud_dur:")
            Airport_CloudDuration(:) = real(dum1d_sp(1:nairports),kind=ip)
            deallocate(dum1d_sp)
          else
            allocate(dum1d_dp(1:nairports))
            nSTAT = nf90_get_var(ncid,pt_x_var_id,dum1d_dp)
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var pt_x:")
            Airport_Longitude(:) = real(dum1d_dp(1:nairports),kind=ip)
            nSTAT = nf90_get_var(ncid,pt_y_var_id,dum1d_dp)
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var pt_y:")
            Airport_Latitude(:) = real(dum1d_dp(1:nairports),kind=ip)
            nSTAT = nf90_get_var(ncid,pt_asharrival_var_id,dum1d_dp)
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var pt_asharrival:")
            Airport_AshArrivalTime(:) = real(dum1d_dp(1:nairports),kind=ip)
            nSTAT = nf90_get_var(ncid,pt_ashduration_var_id,dum1d_dp)
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var pt_depodur:")
            Airport_AshDuration(:) = real(dum1d_dp(1:nairports),kind=ip)
            nSTAT = nf90_get_var(ncid,pt_cloudarrival_var_id,dum1d_dp)
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var pt_cloudarrival:")
            Airport_CloudArrivalTime(:) = real(dum1d_dp(1:nairports),kind=ip)
            nSTAT = nf90_get_var(ncid,pt_cloudduration_var_id,dum1d_dp)
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var pt_cloud_dur:")
            Airport_CloudDuration(:) = real(dum1d_dp(1:nairports),kind=ip)
            deallocate(dum1d_dp)
          endif

          if (IsLatLon) then
            do i=1,nairports
              Airport_i(i) = int((Airport_Longitude(i)-lonLL)/de) +1
              Airport_j(i) = int((Airport_Latitude(i)-latLL)/dn) +1
              ! Ash has arrived if the arrival time is positive
              if(Airport_AshArrivalTime(i).gt.0.0_ip)then
                Airport_AshArrived(i) = .true.
              endif
              if(Airport_CloudArrivalTime(i).gt.0.0_ip)then
                Airport_CloudArrived(i) = .true.
              endif
            enddo
          else
            do i=1,nairports
              lon_in = Airport_Longitude(i)
              lat_in = Airport_Latitude(i)
              call PJ_proj_for(lon_in,lat_in, A3d_iprojflag, &
                         A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2,A3d_k0_scale,A3d_Re, &
                         xout,yout)
              xnow = real(xout,kind=ip)
              ynow = real(yout,kind=ip)

              Airport_i(i) = int((xnow-xLL)/dx) +1
              Airport_j(i) = int((ynow-yLL)/dy) +1
              ! Ash has arrived if the arrival time is positive
              if(Airport_AshArrivalTime(i).gt.0.0_ip)then
                Airport_AshArrived(i) = .true.
              endif
              if(Airport_CloudArrivalTime(i).gt.0.0_ip)then
                Airport_CloudArrived(i) = .true.
              endif
            enddo
          endif

        endif

        ! Time is always dp
        ! Note: if the netcdf file does not have time_native as a variable, then
        !       tn_len will be 0
        allocate(dum1d_dp(1:tn_len))
        nSTAT = nf90_get_var(ncid,tn_var_id,dum1d_dp)
!              nSTAT=nf90_get_var(ncid,x_var_id,dum1d_sp,  &
!                       start = (/1/),       &
!                       count = (/x_len/))
!              if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var x")
!              lon_cc_pd(1:nxmax) = real(dum1d_sp(1:nxmax),kind=ip)

        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_var tn")
          dum1d_dp = 0.0_dp
        endif
        time_native(1:tn_len) = real(dum1d_dp(1:tn_len),kind=dp)
        deallocate(dum1d_dp)

        if(Write_PR_Data)then
          ! pr_x
          nSTAT = nf90_inq_varid(ncid,"pr_x",pr_x_var_id)
          if(nSTAT.ne.0)then
            pr_x_var_id = 0
            call NC_check_status(nSTAT,0,"inq_varid pr_x")
          endif
          ! pr_y
          nSTAT = nf90_inq_varid(ncid,"pr_y",pr_y_var_id)
          if(nSTAT.ne.0)then
            pr_y_var_id = 0
            call NC_check_status(nSTAT,0,"inq_varid pr_y")
          endif          
          ! pr_name
          nSTAT = nf90_inq_varid(ncid,"pr_name",pr_name_var_id)
          if(nSTAT.ne.0)then
            pr_name_var_id = 0
            call NC_check_status(nSTAT,0,"inq_varid pr_name")
          endif          
          ! pr_ash
          nSTAT = nf90_inq_varid(ncid,"pr_ash",pr_ash_var_id)
          if(nSTAT.ne.0)then
            pr_ash_var_id = 0
            call NC_check_status(nSTAT,0,"inq_varid pr_ash")
          endif          

          ! Allocate profile variables
          call Allocate_Profile(nzmax,tn_len,nvprofiles)
          allocate(x_vprofile(nvprofiles))
          allocate(y_vprofile(nvprofiles))
          allocate(Site_vprofile(nvprofiles))

          if(fop.eq.4)then
            ! Read x,y at profile point locations
            allocate(dum1d_sp(1:nvprofiles))
            nSTAT = nf90_get_var(ncid,pr_x_var_id,dum1d_sp)
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var pr_x:")
            x_vprofile(:) = real(dum1d_sp(1:nvprofiles),kind=ip)
            nSTAT = nf90_get_var(ncid,pr_y_var_id,dum1d_sp)
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var pr_y:")
            y_vprofile(:) = real(dum1d_sp(1:nvprofiles),kind=ip)
            deallocate(dum1d_sp)
          else
            ! Read x,y at profile point locations
            allocate(dum1d_dp(1:nvprofiles))
            nSTAT = nf90_get_var(ncid,pr_x_var_id,dum1d_dp)
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var pr_x:")
            x_vprofile(:) = real(dum1d_dp(1:nvprofiles),kind=ip)
            nSTAT = nf90_get_var(ncid,pr_y_var_id,dum1d_dp)
            if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var pr_y:")
            y_vprofile(:) = real(dum1d_dp(1:nvprofiles),kind=ip)
            deallocate(dum1d_dp)
          endif

          ! Read vprofile names
          nSTAT = nf90_get_var(ncid,pr_name_var_id,Site_vprofile)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var pr_name:")

          ! Read TS of ash concentration at profile point locations
          allocate(dum3d_out(1:nzmax,1:tn_len,1:nvprofiles))
          nSTAT = nf90_get_var(ncid,pr_ash_var_id,dum3d_out)
          pr_ash(:,:,:) = real(dum3d_out(:,:,:),kind=ip)
          deallocate(dum3d_out)
          if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var pr_ash:")
        endif

        ! Now populate a few of the header values needed for
        ! post-processing/annotation
        ! Get volcano name
        nSTAT = nf90_get_att(ncid,nf90_global,"b1l1",cdf_b1l1)
        iendstr = SCAN(cdf_b1l1, "#")
        if (iendstr.eq.1)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: ","Volcano name cannot start with #"
          endif;enddo
          stop 1
        endif
        VolcanoName = trim(adjustl(cdf_b1l1(1:iendstr-1)))

        ! Get date of run
        nSTAT = nf90_get_att(ncid,nf90_global,"date",os_time_log)
        ! Get windfile info
        nSTAT = nf90_get_att(ncid,nf90_global,"b3l1",cdf_b3l1)
        ! Get Simulation time info
        nSTAT = nf90_get_att(ncid,nf90_global,"b3l3",cdf_b3l3)
        read(cdf_b3l3,*,iostat=iostatus,iomsg=iomessage) Simtime_in_hours        ! simulated transport time
        linebuffer050 = "Reading simtime from cdf_b3l3"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,cdf_b3l3,iomessage)
        ! Get vent location
        nSTAT = nf90_get_att(ncid,nf90_global,"b1l5",cdf_b1l5)
        if (IsLatLon) then                        !get lon_volcano and lat_volcano
          read(cdf_b1l5,*,iostat=iostatus,iomsg=iomessage)lon_volcano, lat_volcano
          linebuffer050 = "Reading vlon,vlat from cdf_b3l5"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,cdf_b1l5,iomessage)
        else
          read(cdf_b1l5,*,iostat=iostatus,iomsg=iomessage)x_volcano, y_volcano
          linebuffer050 = "Reading vx,vy from cdf_b3l5"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,cdf_b1l5,iomessage)
          call PJ_proj_inv(real(x_volcano,kind=dp),real(y_volcano,kind=dp), &
                     A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
                     A3d_k0_scale,A3d_Re, &
                     lon_in,lat_in)
          lon_volcano = real(lon_in,kind=ip)
          lat_volcano = real(lat_in,kind=ip)
        endif

        nSTAT = nf90_get_att(ncid,nf90_global,"institution",cdf_institution)
        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_att institution:")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find att institution: Assuming institution=N/A"
          endif;enddo
          cdf_institution="N/A"
        endif

        nSTAT = nf90_get_att(ncid,nf90_global,"run_class",cdf_run_class)
        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_att run_class:")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find att institution: Assuming run class=Analysis"
          endif;enddo
          cdf_run_class="Analysis"
        endif

        nSTAT = nf90_get_att(ncid,nf90_global,"url",cdf_url)
        if(nSTAT.ne.0)then
          call NC_check_status(nSTAT,0,"get_att url:")
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Did not find att url: Assuming url=https://vsc-ash.wr.usgs.gov/ash3d-gui"
          endif;enddo
          cdf_url="https://vsc-ash.wr.usgs.gov/ash3d-gui"
        endif

        ! Get eruptions ESPs
#ifdef USEPOINTERS
        if(.not.associated(e_StartTime))then
#else
        if(.not.allocated(e_StartTime))then
#endif
          allocate(e_StartTime(neruptions))
          e_StartTime(:) = 0.0_ip
        endif
#ifdef USEPOINTERS
        if(.not.associated(e_Duration))then
#else
        if(.not.allocated(e_Duration))then
#endif
          allocate(e_Duration(neruptions))
          e_Duration(:) = 0.0_ip
        endif
#ifdef USEPOINTERS
        if(.not.associated(e_PlumeHeight))then
#else
        if(.not.allocated(e_PlumeHeight))then
#endif
          allocate(e_PlumeHeight(neruptions))
          e_PlumeHeight(:) = 0.0_ip
        endif
#ifdef USEPOINTERS
        if(.not.associated(e_Volume))then
#else
        if(.not.allocated(e_Volume))then
#endif
          allocate(e_Volume(neruptions))
          e_Volume(:) = 0.0_ip
        endif

        ! Double-checking float vs double for er_stime
        nSTAT = nf90_inq_varid(ncid,"er_stime",er_stime_var_id)
        if(nSTAT.ne.0)then
          er_stime_var_id = 0
          call NC_check_status(nSTAT,1,"inq_varid er_stime")
        endif
        nSTAT = nf90_inquire_variable(ncid, er_stime_var_id, invar, xtype = var_xtype)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inquire_variable er_stime")
        if(var_xtype.eq.NF90_FLOAT)then
          allocate(dum1d_sp(1:neruptions))
          nSTAT=nf90_get_var(ncid,er_stime_var_id,dum1d_sp,(/1/))
          e_StartTime = real(dum1d_sp,kind=ip)
          deallocate(dum1d_sp)
        elseif(var_xtype.eq.NF90_DOUBLE)then
          allocate(dum1d_dp(1:neruptions))
          nSTAT=nf90_get_var(ncid,er_stime_var_id,dum1d_dp,(/1/))
          e_StartTime = real(dum1d_dp,kind=ip)
          deallocate(dum1d_dp)
        endif

        ! Double-checking float vs double for er_duration
        nSTAT = nf90_inq_varid(ncid,"er_duration",er_duration_var_id)
        if(nSTAT.ne.0)then
          er_duration_var_id = 0
          call NC_check_status(nSTAT,1,"inq_varid er_duration")
        endif
        nSTAT = nf90_inquire_variable(ncid, er_duration_var_id, invar, xtype = var_xtype)
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"inquire_variable er_duration")
        if(var_xtype.eq.NF90_FLOAT)then
          allocate(dum1d_sp(1:neruptions))
          nSTAT=nf90_get_var(ncid,er_duration_var_id,dum1d_sp,(/1/))
          e_Duration = real(dum1d_sp,kind=ip)
          deallocate(dum1d_sp)
        elseif(var_xtype.eq.NF90_DOUBLE)then
          allocate(dum1d_dp(1:neruptions))
          nSTAT=nf90_get_var(ncid,er_duration_var_id,dum1d_dp,(/1/))
          e_Duration = real(dum1d_dp,kind=ip)
          deallocate(dum1d_dp)
        endif

        ! Eruption start time from file is in hours since BaseYear
        ! We need to reset e_StartTime to be the offset from the SimStartHour
        SimStartHour = e_StartTime(1)
        e_StartTime(:) = e_StartTime(:) - SimStartHour
        allocate(dum1d_out(1:neruptions))
        nSTAT = nf90_inq_varid(ncid,"er_plumeheight",er_plumeheight_var_id)
        if(nSTAT.ne.0)then
          er_plumeheight_var_id = 0
          call NC_check_status(nSTAT,1,"inq_varid er_plumeheight")
        endif
        nSTAT=nf90_get_var(ncid,er_plumeheight_var_id,dum1d_out,(/1/))
        e_PlumeHeight = real(dum1d_out,kind=ip)

        nSTAT = nf90_inq_varid(ncid,"er_volume",er_volume_var_id)
        if(nSTAT.ne.0)then
          er_volume_var_id = 0
          call NC_check_status(nSTAT,1,"inq_varid er_volume")
        endif
        nSTAT=nf90_get_var(ncid,er_volume_var_id,dum1d_out,(/1/))
        e_Volume = real(dum1d_out,kind=ip)
        deallocate(dum1d_out)

        first_time = .false.

      endif ! first_time

      if (.not.present(timestep))then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Found the following time steps in file:"
          write(outlog(io),*)"  Step :  time"
          do it = 1,nWriteTimes
            write(outlog(io),*)it,WriteTimes(it)
          enddo
          write(outlog(io),*)'Enter timestep for initialization'
        endif;enddo

        read(input_unit,*,iostat=iostatus,iomsg=iomessage) linebuffer080
        linebuffer050 = "Reading init_tstep from stdin"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) init_tstep
        linebuffer050 = "Reading init_tstep from stdin linebuffer"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      else
        if (timestep.eq.-1)then
          isFinal_TS = .true.
          init_tstep = t_len
        else
          init_tstep = timestep
        endif
      endif
      if(init_tstep.gt.t_len)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Requested time step index is greater than available."
        endif;enddo
        stop 1
      elseif(init_tstep.lt.0)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Requested time step index is invalid."
        endif;enddo
        stop 1
      endif
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Requested timestep = ",init_tstep,real(WriteTimes(init_tstep),kind=sp)
      endif;enddo

      nSTAT=nf90_get_var(ncid,t_var_id,dumscal_out,(/init_tstep/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var t")
      time = real(dumscal_out,kind=ip)

      ! Now that we have all the header/dim/var info, we can read the data
      ! for the time step in question

      ! Load all 2-d variables for this time step
      allocate(ashcon(x_len,y_len,z_len,bn_len))
      allocate(dum2d_out(x_len,y_len))
      allocate(dum2dint_out(x_len,y_len))

      ! Full concentration array
      if(useRestartVars)then
        nSTAT=nf90_get_var(ncid,ashcon_var_id,ashcon,  &
                 start = (/1,1,1,1/),       &
                 count = (/x_len,y_len,z_len,bn_len/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var ashcon")
#ifdef USEPOINTERS
        if(.not.associated(ashcon_tot) )then
#else
        if(.not.allocated(ashcon_tot) )then
#endif
          allocate(ashcon_tot(x_len,y_len,z_len))
        endif
        ashcon_tot = 0.0_op
        if(n_gs_max.gt.0)then
          do isize=1,n_gs_max
            ashcon_tot(1:x_len,1:y_len,1:z_len) =  &
             ashcon_tot(1:x_len,1:y_len,1:z_len) + &
             real(ashcon(1:x_len,1:y_len,1:z_len,isize),kind=op)
          enddo
        endif
      endif

      ! Deposit Thickness
      if(isFinal_TS)then
        ! Deposit Thickness (Final)
        if(depothickFin_var_id.eq.0)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: Trying to read depothickFin, but the variable is not available"
          endif;enddo
          stop 1
        endif
        nSTAT=nf90_get_var(ncid,depothickFin_var_id,dum2d_out,  &
                 start = (/1,1/),       &
                 count = (/x_len,y_len/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var depothickFin")
#ifdef USEPOINTERS
        if(.not.associated(DepositThickness))then
#else
        if(.not.allocated(DepositThickness))then
#endif
          allocate(DepositThickness(x_len,y_len))
          DepositThickness(:,:) = 0.0_ip
        endif
        DepositThickness = real(dum2d_out,kind=ip)
      else
        ! Deposit Thickness (Time-series) 
        nSTAT=nf90_get_var(ncid,depothick_var_id,dum2d_out,  &
                 start = (/1,1,init_tstep/),       &
                 count = (/x_len,y_len,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var depothick")
#ifdef USEPOINTERS
        if(.not.associated(DepositThickness))then
#else
        if(.not.allocated(DepositThickness))then
#endif
          allocate(DepositThickness(x_len,y_len))
          DepositThickness(:,:) = 0.0_ip
        endif
        DepositThickness = real(dum2d_out,kind=ip)
      endif

      ! Deposit Arrival Time
      if(depotime_var_id.eq.0)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Trying to read depotime, but the variable is not available"
        endif;enddo
        stop 1
      endif
      nSTAT=nf90_get_var(ncid,depotime_var_id,dum2d_out,  &
               start = (/1,1/),       &
               count = (/x_len,y_len/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var depotime")
#ifdef USEPOINTERS
      if(.not.associated(DepArrivalTime))then
#else
      if(.not.allocated(DepArrivalTime))then
#endif
        allocate(DepArrivalTime(x_len,y_len))
        DepArrivalTime(:,:) = 0.0_ip
      endif
      DepArrivalTime = real(dum2d_out,kind=ip)

      ! Ash Cloud Arrival Time
      if(ashcloudtime_var_id.eq.0)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Trying to read ashcloudtime, but the variable is not available"
        endif;enddo
        stop 1
      endif
      nSTAT=nf90_get_var(ncid,ashcloudtime_var_id,dum2d_out,  &
               start = (/1,1/),       &
               count = (/x_len,y_len/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var ashcloudtime")
#ifdef USEPOINTERS
      if(.not.associated(CloudArrivalTime))then
#else
      if(.not.allocated(CloudArrivalTime))then
#endif
        allocate(CloudArrivalTime(x_len,y_len))
        CloudArrivalTime(:,:) = 0.0_ip
      endif
      CloudArrivalTime = real(dum2d_out,kind=ip)

      ! Ash-cloud concentration
      if(ashconMax_var_id.eq.0)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Trying to read ashconMax, but the variable is not available"
        endif;enddo
        stop 1
      endif
      nSTAT=nf90_get_var(ncid,ashconMax_var_id,dum2d_out,  &
               start = (/1,1,init_tstep/),       &
               count = (/x_len,y_len,1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var ashconMax_var")
#ifdef USEPOINTERS
      if(.not.associated(MaxConcentration))then
#else
      if(.not.allocated(MaxConcentration))then
#endif
        allocate(MaxConcentration(x_len,y_len))
        MaxConcentration(:,:) = 0.0_ip
      endif
      MaxConcentration = real(dum2d_out,kind=ip)

      ! Ash-cloud max height
      if(ashheight_var_id.eq.0)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Trying to read cloud_height, but the variable is not available"
        endif;enddo
        stop 1
      endif
      nSTAT=nf90_get_var(ncid,ashheight_var_id,dum2d_out,  &
               start = (/1,1,init_tstep/),       &
               count = (/x_len,y_len,1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var cloud_height")
#ifdef USEPOINTERS
      if(.not.associated(MaxHeight))then
#else
      if(.not.allocated(MaxHeight))then
#endif
        allocate(MaxHeight(x_len,y_len))
        MaxHeight(:,:) = 0.0_ip
      endif
      MaxHeight = real(dum2d_out,kind=ip)
      ! Ash-cloud bottom
      nSTAT=nf90_get_var(ncid,ashcloudBot_var_id,dum2d_out,  &
               start = (/1,1,init_tstep/),       &
               count = (/x_len,y_len,1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var cloud_bottom")
#ifdef USEPOINTERS
      if(.not.associated(MinHeight))then
#else
      if(.not.allocated(MinHeight))then
#endif
        allocate(MinHeight(x_len,y_len))
        MinHeight(:,:) = 0.0_ip
      endif
      MinHeight = real(dum2d_out,kind=ip)

      ! Cloud-load
      if(ashload_var_id.eq.0)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Trying to read cloud_load, but the variable is not available"
        endif;enddo
        stop 1
      endif
      nSTAT=nf90_get_var(ncid,ashload_var_id,dum2d_out,  &
               start = (/1,1,init_tstep/),       &
               count = (/x_len,y_len,1/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var cloud_load")
#ifdef USEPOINTERS
      if(.not.associated(CloudLoad))then
#else
      if(.not.allocated(CloudLoad))then
#endif
        allocate(CloudLoad(x_len,y_len))
        CloudLoad(:,:) = 0.0_ip
      endif
      CloudLoad = real(dum2d_out,kind=ip)

      ! Cloud-mask
      if(cloudmask_var_id.eq.0)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"WARNING: Trying to read cloud_mask, but the variable is not available"
          write(errlog(io),*)"         Recreating from CloudLoad"
        endif;enddo
        dum2dint_out = 0
        do i=1,x_len
          do j=1,y_len
            if(CloudLoad(i,j).ge.CLOUDLOAD_THRESH)then
              dum2dint_out(i,j) = 1
            endif
          enddo
        enddo
      else
        nSTAT=nf90_get_var(ncid,cloudmask_var_id,dum2dint_out,  &
                 start = (/1,1,init_tstep/),       &
                 count = (/x_len,y_len,1/))
        if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"get_var cloud_mask")
      endif
#ifdef USEPOINTERS
      if(.not.associated(Mask_Cloud))then
#else
      if(.not.allocated(Mask_Cloud))then
#endif
        allocate(Mask_Cloud(x_len,y_len))
        Mask_Cloud(:,:) = .false.
      endif
      do i=1,x_len
        do j=1,y_len
          if(dum2dint_out(i,j).eq.1)then
            Mask_Cloud(i,j) = .true.
          else
            Mask_Cloud(i,j) = .false.
          endif
        enddo
      enddo

      ! Radar_Reflec
      allocate(dum3d_out(x_len,y_len,z_len))
      if(radrefl_var_id.eq.0)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Trying to read radar_reflectivity, but the variable is not available"
        endif;enddo
        stop 1
      endif
      nSTAT=nf90_get_var(ncid,radrefl_var_id,dum3d_out,  &
               start = (/1,1,1,init_tstep/),       &
               count = (/x_len,y_len,z_len/))
      if(nSTAT.ne.0)call NC_check_status(nSTAT,0,"get_var radar_reflectivity")
#ifdef USEPOINTERS
      if(.not.associated(dbZCol))then
#else
      if(.not.allocated(dbZCol))then
#endif
        allocate(dbZCol(x_len,y_len))
        dbZCol(:,:) = 0.0_ip
      endif
      do i=1,x_len
        do j=1,y_len
          dbZCol(i,j) = real(maxval(dum3d_out(i,j,:)),kind=ip)
        enddo
      enddo
      deallocate(dum3d_out)

      ! Load contour levels in case we are post-processing and plotting
      call Set_OutVar_ContourLevel

      ! Cleaning up
#ifdef USEPOINTERS
      if(associated(ashcon_tot))   deallocate(ashcon_tot)
#else
      if(allocated(ashcon_tot))   deallocate(ashcon_tot)
#endif
      if(allocated(ashcon))       deallocate(ashcon)
      if(allocated(dum2d_out))    deallocate(dum2d_out)
      if(allocated(dum2dint_out)) deallocate(dum2dint_out)

      ! Close file
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Closing netCDF file."
      endif;enddo      
      nSTAT = nf90_close(ncid)
      if(nSTAT.ne.0)call NC_check_status(nSTAT,1,"nf90_close")

      end subroutine NC_Read_Output_Products
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module Ash3d_Netcdf_IO

!##############################################################################
