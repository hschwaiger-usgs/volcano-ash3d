
      !subroutine Read_Control_File(fc_inputfile)
      subroutine Read_Control_File

      ! Subroutine that reads ASCII input file and contains error traps for input
      use iso_c_binding
      use precis_param

      use io_units

      use global_param,  only : &
         EPS_SMALL,EPS_TINY,nmods,OPTMOD_names,VERB,limiter,&
         useDS,useTemperature,useCalcFallVel,useVariableGSbins,&
         useDiffusion,useCN,KM3_2_M3

      use io_data,       only : &
         cdf_b1l1,cdf_b1l2,cdf_b1l3,cdf_b1l4,cdf_b1l5,cdf_b1l6,cdf_b1l7,cdf_b1l8,cdf_b1l9,&
         cdf_b3l1,cdf_b3l2,cdf_b3l3,cdf_b3l4,cdf_b3l5,cdf_b4l1,cdf_b4l2,cdf_b4l3,cdf_b4l4,&
         cdf_b4l5,cdf_b4l6,cdf_b4l7,cdf_b4l8,cdf_b4l9,cdf_b4l10,cdf_b4l11,cdf_b4l12,cdf_b4l13,&
         cdf_b4l14,cdf_b4l15,cdf_b4l16,cdf_b4l17,cdf_b4l18,cdf_b6l1,cdf_b6l2,cdf_b6l3,cdf_b6l4,&
         cdf_b6l5,cdf_comment,cdf_title,cdf_institution,cdf_source,cdf_history,cdf_references,&
         outfile,VolcanoName,WriteTimes,nWriteTimes,cdf_conventions,cdf_run_class,cdf_url,&
         x_vprofile,y_vprofile,i_vprofile,j_vprofile,Site_vprofile,&
         concenfile,infile,ioutputFormat,LoadConcen,log_step,NextWriteTime,&
         AppendExtAirportFile,WriteInterval,WriteGSD,WriteDepositTS_KML,WriteDepositTS_ASCII,&
         WriteDepositTime_KML,WriteDepositTime_ASCII,WriteDepositFinal_KML,&
         WriteDepositFinal_ASCII,WriteCloudTime_KML,WriteCloudTime_ASCII,&
         WriteCloudLoad_KML,WriteReflectivity_KML,WriteCloudLoad_ASCII,WriteCloudHeight_KML,&
         WriteCloudHeight_ASCII,WriteCloudConcentration_KML,WriteCloudConcentration_ASCII,&
         WriteAirportFile_KML,WriteAirportFile_ASCII,Write3dFiles,ReadExtAirportFile,&
         Output_every_TS,Output_at_WriteTimes,Output_at_logsteps,nvprofiles,iTimeNext,&
         Write_PT_Data,Write_PR_Data

      use Source,        only : &
         neruptions,e_Duration,e_Volume,e_PlumeHeight,e_prof_Volume,e_prof_dz,&
         MassFlux,e_EndTime,e_prof_MassFlux,e_prof_zpoints,e_StartTime,&
         ESP_duration,ESP_height,ESP_Vol,e_EndTime_final,ibase,itop,&
         lat_volcano,lon_volcano,x_volcano,y_volcano,z_volcano,Suzuki_A,&
         IsCustom_SourceType,SourceType,&
           Allocate_Source_eruption

      use Tephra,        only : &
         DepositDensity,MagmaDensity,Tephra_v_s,Tephra_gsdiam,Tephra_bin_mass,Tephra_rho_m,&
         Tephra_gsF,Tephra_gsG,FV_ID,phi_mean,phi_stddev,n_gs_max,n_gs_aloft,&
           Calculate_Tephra_Shape,&
           Allocate_Tephra, &
           Sort_Tephra_Size

      use mesh,          only : &
         de,dn,dx,dy,z_vec_init,dz_const,nxmax,nymax,nzmax,nsmax,VarDzType,ivent,jvent,&
         gridwidth_e,gridwidth_n,gridwidth_x,gridwidth_y,&
         lonLL,latLL,lonUR,latUR,xLL,yLL,xUR,yUR,&
         A3d_iprojflag,A3d_k0_scale,A3d_phi0,A3d_lam0,A3d_lam1,A3d_phi1,A3d_lam2,&
         A3d_phi2,A3d_radius_earth,IsLatLon,IsPeriodic,ZPADDING

      use solution,      only : &
         StopValue,imin,imax,jmin,jmax,kmin,kmax

      use Output_Vars,   only : &
         USE_OUTPROD_VARS, USE_RESTART_VARS

      use time_data,     only : &
         BaseYear,useLeap,time,SimStartHour,Simtime_in_hours,xmlSimStartTime

      use Airports,      only : &
         AirportInFile,&
           ProjectAirportLocations

      use VotW_ESP,      only : &
           get_ESP

      use Diffusion,     only : &
         diffusivity_horz,diffusivity_vert,&
           Allocate_Diff

      use projection,    only : &
         PJ_iprojflag,PJ_k0,PJ_lam0,PJ_lam1,PJ_lam2,PJ_phi0,PJ_phi1,PJ_phi2,&
         PJ_radius_earth,&
           PJ_Set_Proj_Params

      use MetReader,     only : &
         MR_iwindfiles,MR_windfiles,MR_BaseYear,MR_useLeap,MR_Comp_StartHour,&
         MR_windfiles_GRIB_index,MR_windfiles_Have_GRIB_index,MR_Comp_Time_in_hours,&
         MR_windfile_starthour,MR_windfile_stephour,MR_iHeightHandler,&
         MR_iwf_template,MR_iwind,&
         MR_global_essential,MR_global_production,MR_global_debug,&
         MR_global_info,MR_global_log,MR_global_error, &
           MR_Allocate_FullMetFileList, &
           MR_Read_Met_DimVars

#ifdef USENETCDF
      use Ash3d_Netcdf
#endif

      implicit none

      !character(kind=c_char), dimension(1:130) :: fc_inputfile
      character, dimension(1:130) :: fc_inputfile

      integer           :: i,k,ii
      !integer           :: iargc
      integer           :: nargs          ! number of command-line arguments

      integer, allocatable, dimension(:)       :: iyear  ! time data read from files
      integer, allocatable, dimension(:)       :: imonth
      integer, allocatable, dimension(:)       :: iday
      real(kind=dp), allocatable, dimension(:) :: hour   ! Start time of eruption in
                                                         !  hour (UT)
      character(len=80) :: linebuffer080
      character(len=130):: linebuffer130
      character(len=3)  :: answer
      character(len=6)  :: formatanswer
      character(len=20) :: mod_name
      character(len=20) :: dumstr20

      integer           :: iw,iwf,igrid,idf,iwfiles
      integer           :: ivalue1, ivalue2, ivalue3, ivalue4
      integer           :: status
      integer           :: loc
      integer           :: iform

      character(len=80) :: Comp_projection_line
      integer           :: ilatlonflag
      character         :: testkey,testkey2
      integer           :: iendstr,ios,ioerr,init_n_gs_max
      real(kind=ip)     :: value1, value2, value3, value4, value5
      real(kind=dp)     :: tmp_dp
      !real(kind=dp)     :: StartHour
      !real(kind=dp)     :: RunStartHour    ! Start time of model run, in hours since BaseYear
      real(kind=ip)     :: sum_bins
      character(len=8)  :: volc_code
      real(kind=ip),allocatable,dimension(:) :: dum_prof

      real(kind=ip),allocatable,dimension(:) :: temp_v_s,temp_gsdiam
      real(kind=ip),allocatable,dimension(:) :: temp_bin_mass,temp_rho_m
      real(kind=ip),allocatable,dimension(:) :: temp_gsF,temp_gsG,temp_phi
      real(kind=ip)     :: fracfine = 0.0_ip
      real(kind=ip)     :: CompGrid_height
      real(kind=ip)     :: last_z
      integer           :: nz_init,nsegments
      integer      ,allocatable,dimension(:) :: nz_plin_segments
      real(kind=ip),allocatable,dimension(:) :: dz_plin_segments
      integer           :: substr_pos1
      integer           :: substr_pos2
      logical           :: IsThere
      !character(len=8)  :: version             =  ' 1.0  '
      logical           :: StopWhenDeposited                       ! If true, StopValue=0.99, else StopValue=1e5.
      logical           :: runAsForecast       = .false.           ! This will be changed if year=0
      real(kind=dp)     :: FC_Offset = 0.0_dp

      !! Size matches length of infile (specified in module io_data)
      integer fc_len

      INTERFACE
        subroutine help_input(blockID)
          integer,intent(in) :: blockID
        end subroutine help_input
        subroutine LatLonChecker(latLL,lonLL,lat_volcano,lon_volcano,gridwidth_e,gridwidth_n)
          integer,parameter  :: ip         = 8 ! Internal precision
          real(kind=ip), intent(in)    :: latLL
          real(kind=ip), intent(in)    :: lonLL
          real(kind=ip), intent(in)    :: lat_volcano
          real(kind=ip), intent(inout) :: lon_volcano ! we might need to remap this value
          real(kind=ip), intent(in)    :: gridwidth_e
          real(kind=ip), intent(in)    :: gridwidth_n
        end subroutine LatLonChecker
        subroutine xyChecker(xLL,yLL,dx,dy,x_volcano,y_volcano,gridwidth_x,gridwidth_y)
          integer,parameter  :: ip         = 8 ! Internal precision
          real(kind=ip), intent(in) :: xLL,yLL
          real(kind=ip), intent(in) :: dx,dy
          real(kind=ip), intent(in) :: x_volcano,y_volcano
          real(kind=ip), intent(in) :: gridwidth_x,gridwidth_y
        end subroutine xyChecker
        subroutine vprofchecker(iprof)
          integer, intent(in) :: iprof
        end subroutine vprofchecker
        real(kind=8) function HS_hours_since_baseyear(iyear,imonth,iday,hours,byear,useLeaps)
          integer            :: iyear
          integer            :: imonth
          integer            :: iday
          real(kind=8)       :: hours
          integer            :: byear
          logical            :: useLeaps
        end function HS_hours_since_baseyear
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
        subroutine MR_Set_Gen_Index_GRIB(grib_file)
          character(len=130),intent(in)  :: grib_file
        end subroutine MR_Set_Gen_Index_GRIB
      END INTERFACE


      write(global_production,*)&
       "--------------------------------------------------"
      write(global_production,*)&
       "---------- READ_CONTROL_FILE ---------------------"
      write(global_production,*)&
       "--------------------------------------------------"

      !initialize output
      formatanswer = 'null'
      nWriteTimes  = 0                  !number of output files to write (default=0)
      NextWriteTime = 1.0_ip/EPS_TINY   !Time to write the next file (default = never)

      ! TEST READ COMMAND LINE ARGUMENTS
      nargs = command_argument_count()
      if (nargs.eq.0) then
          ! If no command-line arguments are given, then prompt user
          ! interactively for the command file name and possible a 
          ! restart file
        write(global_info,*)'Enter name of ESP input file:$'
        read(5,*) infile
        write(global_info,*)'Load concentration file?'
        read(5,'(a3)') answer
        if (answer.eq.'yes') then
          LoadConcen = .true.
        elseif (answer.eq.'no') then
          LoadConcen = .false.
        else
          write(global_info,*) 'Sorry, I cannot understand your answer.'
          stop 1
        endif
        if(LoadConcen)then
          write(global_info,*)'Enter name of concentration file'
          read(5,*) concenfile
#ifdef USENETCDF
          call NC_RestartFile_ReadTimes
#else
          write(global_error,*)"ERROR: ",&
           "Loading concentration files requires previous netcdf"
          write(global_error,*)&
           "       output.  This Ash3d executable was not compiled with"
          write(global_error,*)&
           "       netcdf support.  Please recompile Ash3d with"
          write(global_error,*)&
           "       USENETCDF=T, or select another source."
          stop 1
#endif
        endif
      elseif (nargs.ge.1) then
        if (nargs.gt.1) then
          write(global_info,*)&
           "Only one command-line argument is expected."
          write(global_info,*)&
           "Reading first arguement as the control files and"
          write(global_info,*)&
           "disregarding all other arguements."
        endif
          ! If only one argument is given, first test for the '-h' indicating a help
          ! request.
        !call getarg(1,linebuffer130)
        call get_command_argument(1, linebuffer130, status)
        testkey  = linebuffer130(1:1)
        testkey2 = linebuffer130(2:2)
        if(testkey.eq.'-')then
          if(testkey2.eq.'h')then
            ! There should be some parsing of the help parameters
            ! such as
            ! Ash3d -h     or   Ash3d -h run
            !            : gives information on how to execute
            ! call help_run()
            !
            ! Ash3d -h input
            !            : gives information on input file format
            !call help_input(0)
            !
            ! Ash3d -h make
            !            : gives information on the different make options
            ! call help_make()
            !
            ! Ash3d -h wind
            !            : gives information on the different wind
            !            files, where to get them, etc.
            ! call help_wind()
            !

            ! For now, just assume info is wanted on input file format

            ! Call help routine for block=0 (i.e. all blocks of input
            ! file)
            call help_input(0) 
          !else
          !  This branch is reserved for when we can't figure out
          !  what help the user wants
          endif
          stop 1
        else
          ! If the first argument does not begin with '-', then
          ! assume it is the input file name
          read(linebuffer130,*)infile
        endif
      elseif (nargs < 0) then
        !! When code called from ForestClaw, nargs is -1
        fc_len = 0
        do
          if (fc_inputfile(fc_len+1) == C_NULL_CHAR) exit
          fc_len = fc_len + 1
          infile(fc_len:fc_len) = fc_inputfile(fc_len)
        end do
        write(global_info,*) 'Reading input file ''',&
                             infile,''' from ForestClaw'
      endif


      !OPEN AND READ ESP FILE
      write(global_info,3) infile
      write(global_log ,3) infile

      inquire( file=infile, exist=IsThere )
      if(.not.IsThere)then
        write(global_error,*)"ERROR: Cannot file input file"
        stop 1
      endif
      open(unit=10,file=infile,status='old',err=1900)

      !************************************************************************
      ! BLOCK 1: GRID INFO
      ! Start reading the input file assuming there is a variable length
      ! header with each header line flagged by a '#' or '*'
      write(global_info,*)' *******************************************'
      write(global_info,*)' Reading Block 1: Volcano/grid specification'
      write(global_info,*)' *******************************************'
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      do while(testkey.eq.'#'.or.testkey.eq.'*')
         ! Line is a comment, read next line
        read(10,'(a80)')linebuffer080
        read(linebuffer080,*)testkey
      enddo

      ! Block 1 Line 1
      ! Read volcano name
      cdf_b1l1 = linebuffer080
      iendstr = SCAN(linebuffer080, "#")
      if (iendstr.eq.1)then
        write(global_error,*)"ERROR: ",&
         "Volcano name cannot start with #"
        stop 1
      endif
      VolcanoName = trim(adjustl(linebuffer080(1:iendstr-1)))
      ! Check if the volcano name is a text name or a Smithsonian
      ! database ID
      read(VolcanoName,*)testkey
      if(testkey.eq.'0'.or.testkey.eq.'1')then
        ! the 'name' is the CAVW Smithsonian ID
        ! get the source parameters for this volcano
        if(VERB.gt.1)write(global_info,*)"Calling get_ESP"
        volc_code = VolcanoName(1:8)
        call get_ESP(volc_code)
      endif

      write(global_info,*)
      write(global_log ,*)
      write(global_info,37) VolcanoName
      write(global_log ,37) VolcanoName

      ! Block 1 Line 2
      ! Read projection parameters
      read(10,'(a80)') linebuffer080
      cdf_b1l2 = linebuffer080
      Comp_projection_line = linebuffer080
      read(Comp_projection_line,*)ilatlonflag
      if (ilatlonflag.eq.0) then
        ! expecting input variables to be in the same projection as
        ! specified by iprojflag and parameters
        IsLatLon          = .false.
      else
        ! expecting input variables to be in lat/lon
       IsLatLon          = .true.
      endif

      ! Set Projection Parameters
      if (IsLatLon.eqv..false.)then
        call PJ_Set_Proj_Params(Comp_projection_line)
        A3d_iprojflag    = PJ_iprojflag
        A3d_k0_scale     = PJ_k0
        A3d_radius_earth = PJ_radius_earth
        A3d_lam0         = PJ_lam0
        A3d_lam1         = PJ_lam1
        A3d_lam2         = PJ_lam2
        A3d_phi0         = PJ_phi0
        A3d_phi1         = PJ_phi1
        A3d_phi2         = PJ_phi2
      endif

      ! READ BOUNDARIES OF MODEL DOMAIN
      if(IsLatLon)then
        ! If input coordinates are in lat/lon, interpret lines as follows
        ! Block 1 Line 3
        read(10,'(a80)')cdf_b1l3
        read(cdf_b1l3,*,err=1901) lonLL, latLL            ! lat/lon of LL corner
        ! Block 1 Line 4
        read(10,'(a80)')cdf_b1l4
        read(cdf_b1l4,*,err=1902) gridwidth_e, gridwidth_n          ! Dimensions (in degrees) of the grid
        ! Block 1 Line 5
        read(10,'(a80)') cdf_b1l5
        read(cdf_b1l5,*,iostat=ioerr) value1, value2, value3
        if (ioerr.eq.0)then
          lon_volcano = value1
          lat_volcano = value2
          z_volcano   = value3
        else
          lon_volcano = value1
          lat_volcano = value2
          z_volcano   = 0.0_ip
        endif
        ! Block 1 Line 6
        read(10,'(a80)')cdf_b1l6
        read(cdf_b1l6,*,err=1904) de, dn                 ! cell size in degrees 
        
        !Make sure longitudes are between 0 and 360 degrees
        if (lonLL.lt.-360.0_ip) then
          write(global_info,*)&
           "Please give longitude values between -360 and 360."
          write(global_log ,*)&
           "Please give longitude values between -360 and 360."
          stop 1
        endif
        if (lonLL.lt.  0.0_ip) lonLL=lonLL+360.0_ip
        if (lonLL.ge.360.0_ip) lonLL=mod(lonLL,360.0_ip)
        if (lon_volcano.lt.-360.0) then
          write(global_info,*)&
           "Please give longitude values between -360 and 360."
          write(global_log ,*)&
           "Please give longitude values between -360 and 360."
          stop 1
        endif
        if (lon_volcano.lt.  0.0_ip) lon_volcano = lon_volcano+360.0_ip
        if (lon_volcano.ge.360.0_ip) lon_volcano = mod(lon_volcano,360.0_ip)

        If(IsLatLon.and.&
           (gridwidth_e.ge.360.0_ip.or.&
            abs(gridwidth_e-360.0_ip).lt.EPS_TINY))then
          IsPeriodic  = .true.
          lonLL       = 0.0_ip
          gridwidth_e = 360.0_ip
        endif
        lonUR = lonLL + gridwidth_e
        latUR = latLL + gridwidth_n

        write(global_info,*)'lonLL=',real(lonLL,kind=sp)
        write(global_info,*)'lonUR=',real(lonUR,kind=sp)
        write(global_info,*)'latLL=',real(latLL,kind=sp)
        write(global_info,*)'latUR=',real(latUR,kind=sp)
        write(global_info,*)'lon_volcano=',real(lon_volcano,kind=sp)
        write(global_info,*)'lat_volcano=',real(lat_volcano,kind=sp)

        write(global_info,4) lonLL, latLL, gridwidth_e, gridwidth_n, &
                   lon_volcano, lat_volcano
        write(global_info,*) "z_volcano = ",real(z_volcano,kind=sp)," km"
        write(global_log ,4) lonLL, latLL, gridwidth_e, gridwidth_n, &
                   lon_volcano, lat_volcano
        write(global_log ,*) "z_volcano = ",real(z_volcano,kind=sp)," km"
        write(global_info,5) de, dn
        write(global_log ,5) de, dn
        
       !check for errors in input
        call LatLonChecker(latLL,lonLL,lat_volcano,lon_volcano,gridwidth_e,gridwidth_n)
      else  ! IsLatLon
        ! Block 1 Line 3
        read(10,'(a80)')cdf_b1l3
        read(cdf_b1l3,*,err=1901) xLL, yLL                ! LL corner in km 
        ! Block 1 Line 4
        read(10,'(a80)')cdf_b1l4
        read(cdf_b1l4,*,err=1902) gridwidth_x, gridwidth_y         ! width and height of simulation area in km
        xUR = xLL + gridwidth_x
        yUR = yLL + gridwidth_y
        ! Block 1 Line 5
        read(10,'(a80)')cdf_b1l5
        read(cdf_b1l5,*,iostat=ioerr) value1, value2, value3
        if (ioerr.eq.0)then
          x_volcano = value1
          y_volcano = value2
          z_volcano = value3
        else
          x_volcano = value1
          y_volcano = value2
          z_volcano = 0.0_ip
        endif
        ! Block 1 Line 6
        read(10,'(a80)')cdf_b1l6
        read(cdf_b1l6,*,err=1904) dx, dy                 ! cell size in horizontal, vertical, in km
        write(global_info,4) xLL, yLL, gridwidth_x, gridwidth_y, x_volcano,y_volcano  !write out input data
        write(global_log ,4) xLL, yLL, gridwidth_x, gridwidth_y, x_volcano,y_volcano
        write(global_info,*) "z_volcano = ",z_volcano," km"
        write(global_log ,*) "z_volcano = ",z_volcano," km"
        write(global_info,5) dx, dy
        write(global_log ,5) dx, dy
        call xyChecker(xLL,yLL,dx,dy,x_volcano,y_volcano,gridwidth_x,gridwidth_y)
      endif

      ! Block 1 Line 7
      read(10,'(a80)')cdf_b1l7
      read(cdf_b1l7,*,err=5215) dz_const    ! nodal spacing in z (always km)
      VarDzType = "dz_cons"
      write(global_info,43) dz_const
      write(global_log ,43) dz_const
      ! Set up initial z_vector up to 50km or so.  This is to match the variable
      ! dz cases in which the max height is specified.  The computational grid
      ! height will be truncated below to just that needed to cover the plume
      nz_init  = ceiling(100.0_ip/dz_const)+1
      allocate(z_vec_init(0:nz_init))
      z_vec_init = 0.0_ip
      do i=1,nz_init
        z_vec_init(i)=dz_const*(i) ! This the top of cell-boundaries (lower bound at 0)
      enddo
      goto 5220

5215  write(global_info,*)&
       "Could not read dz. Trying to reinterpret as alternate z-spacing"
      write(global_info,*)cdf_b1l7
      read(cdf_b1l7,*,err=1905) VarDzType
      if (VarDzType.eq.'dz_plin')then
        ! Piece-wise linear
        !  Read another line with n-segments, nz1, dz1, nz2, dz2, ...
        write(global_info,*)&
         "z is piecewise linear:  Now reading the segments."
        read(10,'(a80)')cdf_b1l7
        read(cdf_b1l7,*,err=1905) nsegments
        if(nsegments.lt.1)then
          write(global_info,*)"ERROR: ",&
           "nsegments must be positive integer"
          write(global_info,*)&
           "       nsegments = ",nsegments
          stop 1
        endif
        allocate(nz_plin_segments(nsegments))
        allocate(dz_plin_segments(nsegments))
        do i=1,nsegments
          read(10,'(a80)')cdf_b1l7
          read(cdf_b1l7,*,err=1905) nz_plin_segments(i), dz_plin_segments(i)
        enddo
        nz_init = sum(nz_plin_segments(:))
        allocate(z_vec_init(0:nz_init))
        z_vec_init = 0.0_ip
        k = 0
        do i=1,nsegments
          do ii = 1,nz_plin_segments(i)
            if(k.eq.0)then
              last_z = 0.0_ip
            else
              last_z = z_vec_init(k)
            endif
            k=k+1
            z_vec_init(k)=last_z+dz_plin_segments(i)
          enddo
        enddo
        !dz_const = 0.25
      elseif (VarDzType.eq.'dz_clog')then
        write(global_info,*)&
         "Logrithmic dz not yet implemented.  Setting to constant dz=0.25"
        dz_const = 0.25_ip
      elseif (VarDzType.eq.'dz_cust')then
        write(global_info,*)&
         "Custom dz not yet implemented.  Setting to constant dz=0.25"
        dz_const = 0.25_ip
      else
        write(global_info,*)&
         "dz type must be either a number (in km) for constant dz, or"
        write(global_info,*)&
         "dz_plin, dz_clog, or dz_cust for variable dz"
        write(global_info,*)&
         "You entered: ",cdf_b1l7
        write(global_info,*)&
         "Interpreted as: ",VarDzType
        stop 1
      endif

      ! Block 1 Line 8
      !Read this line looking for diffusion coefficient and either a Suzuki constant, 
      !or a plume type ('line', 'point', 'profile', 'umbrella', or 'umbrella_air')
5220  read(10,'(a80)')cdf_b1l8
      read(cdf_b1l8,*,err=5225) diffusivity_horz, Suzuki_A       ! First, try Suzuki coefficient
      SourceType='suzuki'
      !write(global_info,*)"eruption type = suzuki with coefficient of ",Suzuki_A
      goto 5230
      !if the second item is not a number, read SourceType
5225  write(global_info,*)&
       "Source type is not suzuki. Trying to read another standard type"
      read(cdf_b1l8,*) diffusivity_horz, SourceType
      if ((SourceType.eq.'point').or. &
          (SourceType.eq.'Point').or. &
          (SourceType.eq.'POINT')) then
          SourceType='point'
      elseif ((SourceType.eq.'line').or. &
                 (SourceType.eq.'Line').or. &
                 (SourceType.eq.'LINE')) then
          SourceType='line'
      elseif ((SourceType.eq.'profile').or. &
                 (SourceType.eq.'Profile').or. &
                 (SourceType.eq.'PROFILE')) then
          SourceType='profile'
      elseif ((SourceType.eq.'umbrella').or. &
                 (SourceType.eq.'Umbrella').or. &
                 (SourceType.eq.'UMBRELLA')) then
          SourceType='umbrella'
          Suzuki_A = 12.0_ip
      elseif ((SourceType.eq.'umbrella_air').or. &
                 (SourceType.eq.'Umbrella_air').or. &
                 (SourceType.eq.'UMBRELLA_AIR')) then
          !umbrella_air is the same as 'umbrella'
          !but it is assumed to be an airborne run.
          !Thus if gsbins=1, the MER is multiplied by 20
          !to obtain the right rate of umbrella growth.
          SourceType='umbrella_air'
          Suzuki_A = 12.0_ip
      else
        write(global_info,*)&
         "SourceType is not point, line, profile, umbrella or umbrella_air."
        write(global_info,*)&
         "Assuming this is a custom source type."
        write(global_info,*)&
        "For now, just read eruptions start time, duration, and height."
         IsCustom_SourceType = .true.
      endif
      write(global_info,*)"  SourceType = ",SourceType

5230  diffusivity_horz = diffusivity_horz*3.6e-3_ip  !convert diffusion coefficient from m2/s to km2/hr
      diffusivity_vert = diffusivity_horz

      if(abs(diffusivity_horz).lt.EPS_SMALL)then
        useDiffusion = .false.
        write(global_info,*)"Not using turbulent diffusivity."
        write(global_log ,*)"Not using turbulent diffusivity."
      elseif(diffusivity_horz.lt.0.0)then
        write(global_error,*)"ERROR: ",&
         "Diffusivity must be non-negative."
        write(global_log,*)"ERROR: ",&
         "Diffusivity must be non-negative."
        stop 1
      else
        write(global_info,*)"Using constant turbulent diffusivity:  ",&
                  diffusivity_horz/3.6e-3_ip," m2/s"
        write(global_log ,*)"Using constant turbulent diffusivity:  ",&
                  diffusivity_horz/3.6e-3_ip," m2/s"
        useDiffusion = .true.
      endif

      read(10,'(a80)')cdf_b1l9
      read(cdf_b1l9,*,err=1907) neruptions              ! read in number of eruptions or pulses
      write(global_info,*) 'Expecting to read ',neruptions,&
                           ' eruptions lines in Block 2.'
      if (((SourceType.eq.'umbrella').or.(SourceType.eq.'umbrella_air')) &
           .and.(neruptions.gt.1)) then
        write(global_error,*)"ERROR: ",&
         'when SourceType=umbrella, neruptions must equal 1'
        write(global_error,*)&
         'You gave neruptions=',neruptions
        write(global_error,*)&
         'Program stopped'
        write(global_log,*)"ERROR: ",&
         'when SourceType=umbrella, neruptions must equal 1'
        write(global_log ,*)&
         'You gave neruptions=',neruptions
        write(global_log ,*)&
         'Program stopped'
        stop 1
      endif
      !if(dz_const.le.0.0)then
      !  write(global_info,*)"ERROR: dz_const must be positive, not ",dz_const
      !  stop 1
      !endif
      if(SourceType.eq.'suzuki'.and.(Suzuki_A.le.0.0_ip))then
        write(global_error,*)"ERROR: ",&
         "Suzuki_A must be positive, not ",Suzuki_A
        stop 1
      endif
      if(neruptions.le.0)then
        write(global_error,*)"ERROR: ",&
         "neruptions must be positive, not ",neruptions
        stop 1
      endif
      ! END OF BLOCK 1
      !************************************************************************

      ! ALLOCATE ARRAYS OF ERUPTIVE PROPERTIES      
      call Allocate_Source_eruption

      allocate (iyear(neruptions))
      allocate (imonth(neruptions))
      allocate (iday(neruptions))
      allocate (hour(neruptions))
      !************************************************************************
      ! BLOCK 2: ERUPTION PARAMETERS
      ! Again, assuming there is a variable length
      ! header with each header line flagged by a '#' or '*'
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      if (testkey.ne.'#'.and.testkey.ne.'*')then
        write(global_error,*)"ERROR: ",&
         'Expecting a comment line separating blocks.'
        write(global_error,*)&
         '       Check that Block 1 is correct.'
        stop 1
      endif
      do while(testkey.eq.'#'.or.testkey.eq.'*')
         ! Line is a comment, read next line
        read(10,'(a130)')linebuffer130
        read(linebuffer130,*)testkey
      enddo
      write(global_info,*)' *******************************************'
      write(global_info,*)' Reading Block 2: Eruption parameters'
      write(global_info,*)' *******************************************'
      ! BEGIN READING TIMES OF ERUPTIVE PULSES      
      do i=1,neruptions  
        ! Always check if we have overshot the block
        read(linebuffer130,*)testkey
        if (testkey.eq.'#'.or.testkey.eq.'*') then
          write(global_log,*)"ERROR: ",&
           'Trying to read Block 2 and detecting comment line'
          write(global_log,*)&
           '  Eruption ',i,'of',neruptions
          write(global_log,*)&
           '  Offending line: ',linebuffer130
          write(global_error,*)"ERROR: ",&
           'Trying to read Block 2 and detecting comment line'
          write(global_error,*)&
           '  Eruption ',i,'of',neruptions
          write(global_error,*)&
           '  Offending line: ',linebuffer130
          stop 1
        endif
        if(i.eq.1)then
          read(linebuffer130,*,err=1910) iyear(i)
          if(iyear(i).ne.0.and.iyear(i).lt.BaseYear.or.iyear(i)-BaseYear.gt.200)then
            ! Reset BaseYear to the start of the century containing the eruption year
            BaseYear = iyear(i) - mod(iyear(i),100)
            !BaseYear = 1
            write(global_info,*)"WARNING: Resetting BaseYear to ",BaseYear
          endif
          if(iyear(i).eq.0)then  !HFS: KLUDGE-- This should be changed to test for FC or something
                                 !              Start time should also be calculated by Ash3d, not MetReader
            runAsForecast = .true.
            write(global_info,*)"Running as forecast."
          endif
        endif
        if(SourceType.eq.'suzuki'      .or. &
           SourceType.eq.'point'       .or. &
           SourceType.eq.'line'        .or. &
           SourceType.eq.'umbrella'    .or. &
           SourceType.eq.'umbrella_air')then
         !read start time, duration, plume height, volume of each pulse
          read(linebuffer130,*,err=1910) iyear(i),imonth(i),iday(i),hour(i), &
                                e_Duration(i), e_PlumeHeight(i), e_Volume(i)
        elseif(SourceType.eq.'profile')then
          write(global_info,*)"Start reading eruption profiles."
          !read start time, duration, plume height, volume of each pulse
          read(linebuffer130,*,err=1910) iyear(i),imonth(i),iday(i),hour(i), &
                                e_Duration(i), e_PlumeHeight(i), e_prof_dz(i),e_prof_zpoints(i)
          allocate(dum_prof(e_prof_zpoints(i)))
          read(10,*)dum_prof(1:e_prof_zpoints(i))
          e_prof_Volume(i,1:e_prof_zpoints(i))=dum_prof(1:e_prof_zpoints(i))
          deallocate(dum_prof)
          e_Volume(i) = sum(e_prof_Volume(i,:))
        else
          ! This is the custom source.  A special call to a source reader
          ! will need to made from Ash3d_??.F90.  For now, just read the
          ! start time, duration, and plume height
          read(linebuffer130,*,err=1910) iyear(i),imonth(i),iday(i),hour(i),&
                                       e_Duration(i), e_PlumeHeight(i)
          e_Volume(i)    = 0.0_ip
          if(neruptions.gt.1)then
            ! For more than one custom source, the next iteration might cause problems
            ! since a custom source might require multiple input lines per source.
            ! For now, copy slot 1 to all the others and break out of the do loop.
            ! The full source list must be populated by the user-provided custom source
            ! readers.
            iyear(2:neruptions)         = iyear(1)
            imonth(2:neruptions)        = imonth(1)
            iday(2:neruptions)          = iday(1)
            hour(2:neruptions)          = hour(1)
            e_Duration(2:neruptions)    = e_Duration(1)
            e_PlumeHeight(2:neruptions) = e_PlumeHeight(1)
            e_Volume(2:neruptions)      = e_Volume(1)
            exit
          endif

        endif

        ! if we're doing a forecast run, we'll need to add the windfile
        ! reference time to SimStartHour and e_StartTime(i) so
        ! temporarily set year and month to Jan, 1900.  If iday is
        ! considered 'days after start of wind file', then we need to
        ! add 1 so that the hours are calculated properly.
        if(runAsForecast)then
          iyear(i) = BaseYear
          imonth(i) = 1
          iday(i) = iday(i) + 1
          FC_Offset = real(hour(1),sp)
        endif
        if(e_Duration(i).lt.0.0_ip)    e_Duration(i)    = ESP_duration
        if(e_PlumeHeight(i).lt.0.0_ip) e_PlumeHeight(i) = ESP_height
        if(e_Volume(i).lt.0.0_ip)      e_Volume(i)      = ESP_Vol

        ! read next line of input file
        read(10,'(a130)')linebuffer130
      enddo

      !Error trap if more pulses are entered than are specified
      if(IsCustom_SourceType)then
        ! If we are using custom sources, suppress the error-checking
        ! since we don't know the number of lines of input we need and
        ! we will need to loop through here again anyway.
        write(global_info,*)&
          'Skipping error-checking of eruption source lines since '
        write(global_info,*)&
          'the source is a custom type'
         ! For the custom source, we will need to read to the end of block 2
         do while(testkey.ne.'#'.and.testkey.ne.'*')
            ! Line is a comment, read next line
           read(10,'(a80)')linebuffer080
           read(linebuffer080,*)testkey
         enddo

      else
        if (linebuffer130(1:5).ne.'*****') then
          write(global_info,*)&
           'The beginning of the line following the list of', &
           ' eruptive pulses did not'
          write(global_info,*)&
           'start with ''*****''.  Did you enter the correct', &
           '  number of eruptive pulses?'
          write(global_info,*) 'Program stopped.'
          stop 1
        endif
      endif
      ! END OF BLOCK 2
      !************************************************************************

      !************************************************************************
      ! BLOCK 3: WIND PARAMETERS
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      if (.not.IsCustom_SourceType.and.&  ! only perform this check for standard src
          testkey.ne.'#'.and.testkey.ne.'*')then
        write(global_error,*)"ERROR: ",&
         'Expecting a comment line separating blocks.'
        write(global_error,*)&
         '       Check that Block 2 is correct.'
        stop 1
      endif
      do while(testkey.eq.'#'.or.testkey.eq.'*')
         ! Line is a comment, read next line
        read(10,'(a80)')linebuffer080
        read(linebuffer080,*)testkey
      enddo
      write(global_info,*)' *******************************************'
      write(global_info,*)' Reading Block 3: Windfile parameters'
      write(global_info,*)' *******************************************'
      ! Block 3 Line 1
      cdf_b3l1 = linebuffer080
      read(linebuffer080,*,iostat=ioerr) iw,iwf
      idf = 0
      if (ioerr.eq.0)then
        ! Succeeded in reading the two required values, try for three
        read(linebuffer080,*,iostat=ioerr) iw, iwf, ivalue3
        if (ioerr.eq.0)then
          ! Success reading three values, try for four
          igrid = ivalue3
          read(linebuffer080,*,iostat=ioerr) iw, iwf, ivalue3, ivalue4
          if (ioerr.eq.0)then
            ! Success!, set data format (ascii, netcdf, grib)
            idf = ivalue4
          endif
        else
          igrid = 0
        endif
        if(idf.lt.1)then
          ! Data format is not given, assume netcdf unless ascii specified
          if(iw.eq.1.or.iw.eq.2)then
            idf = 1 ! ASCII
          else
            idf = 2 ! Netcdf
          endif
        endif

      else
        ! We need at least two values
        write(global_error,*)"ERROR: ",&
         "could not read iwind, iwindformat"
        stop 1
      endif

      if(iwf.eq.0)then
        ! If iwindformat = 0, then the input file is a not a known format
        ! Read an extra line given the name of a template file.
        if(idf.ne.2)then
          write(global_info,*)&
           "Currently only netcdf reader implemented for templates.",&
           "  Resetting idf to 2"
          idf = 2
        endif
        read(10,'(a80)')linebuffer080
        read(linebuffer080,*)testkey
        if (testkey.eq.'#'.or.testkey.eq.'*') then
          write(global_error,*)"ERROR: ",&
           "Trying to template name and detecting comment line"
          stop 1
        endif
        read(linebuffer080,'(a80)',err=1970) MR_iwf_template
      else
        ! iwf is a known format.
        if(iwf.eq.33)then
            ! The only known format that does not use leap years is the
            ! paleoclimate CAM files.
          useLeap=.false.
        else
          useLeap=.true.
        endif
      endif

      read(10,'(a80)')linebuffer080
      if (testkey.eq.'#'.or.testkey.eq.'*') then
        write(global_error,*)"ERROR: ",&
         "Trying to read MR_iHeightHandler and detecting comment line"
        stop 1
      endif
      ! Block 3 Line 2
      cdf_b3l2 = linebuffer080
      read(linebuffer080,*,err=1932) MR_iHeightHandler ! parameter that determines what to do if the
                                        ! plume height exceeds the wind sounding max. height

      read(10,'(a80)')linebuffer080
      if (testkey.eq.'#'.or.testkey.eq.'*') then
        write(global_error,*)"ERROR: ",&
         "Trying to read Simtime_in_hours and detecting comment line"
        stop 1
      endif
      ! Block 3 Line 3
      cdf_b3l3 = linebuffer080
      read(linebuffer080,*,err=1921) Simtime_in_hours        ! simulated transport time
                                                          ! for ash cloud, in hours

!     Read whether to stop calculation when percent_accumulated>0.99
      read(10,'(a80)') linebuffer080
      if (testkey.eq.'#'.or.testkey.eq.'*') then
        write(global_error,*)"ERROR: ",&
         "Trying to read StopWhenDeposited and detecting comment line"
        stop 1
      endif
      ! Block 3 Line 4
      cdf_b3l4 = linebuffer080
      read(linebuffer080,'(a3)',err=1922) answer
      if (answer.eq.'yes') then
        StopWhenDeposited = .true.
        StopValue = 0.99_ip
       else if (answer(1:2).eq.'no') then
        StopWhenDeposited = .false.
        StopValue = 1.0e2_ip
       else
        goto 1922
      endif

      read(10,'(a80)')linebuffer080
      if (testkey.eq.'#'.or.testkey.eq.'*') then
        write(global_error,*)"ERROR: ",&
         "Trying to read nWindFiles and detecting comment line"
        stop 1
      endif
      ! Block 3 Line 5
      cdf_b3l5 = linebuffer080
      read(linebuffer080,*,err=1923) iwfiles              ! number of wind files to read

      ! Check if there are more lines than expected
      !read(linebuffer080,*)testkey
      !if (testkey.ne.'#'.and.testkey.ne.'*') then
      !  write(global_error,*)'ERROR: Trying to read Block 3 and not detecting the ending comment line'
      !  write(global_error,*)'       Block 3 should look something like this...'
      !  write(global_error,*)'******************* BLOCK 3 ***************************************************'
      !  write(global_error,*)'4 20 4 1            #iwind, iwindFormat, [NCEP grid ID], [data format (netcdf, grib)]'
      !  write(global_error,*)'2                   #iHeightHandler'
      !  write(global_error,*)'48.0                #Simulation time in hours'
      !  write(global_error,*)'yes                 #stop computation when 99% of erupted mass has deposited?'
      !  write(global_error,*)'34                  #nWindFiles, number of gridded wind files (used if iwind>1)'
      !  write(global_error,*)'*******************************************************************************'
      !  write(global_error,*)'  Offending line: ',linebuffer080
      !  write(global_log,*)  'ERROR: Trying to read Block 3 and not detecting the ending comment line'
      !  write(global_log,*)  '  Offending line: ',linebuffer080
      !  stop 1
      !endif
      ! END OF BLOCK 3
      !************************************************************************

      ! Now that we know which calendar we are using (BaseYear, useLeap), now we
      ! can set the HoursSince time for the source terms
      do i=1,neruptions
        !set start time of simulation
        ! Note: This will be updated for forecast runs once we know the start
        !       time of the windfiles
        if (i.eq.1) then
          SimStartHour = HS_hours_since_baseyear(iyear(i),imonth(i),  &
                        iday(i),hour(i),BaseYear,useLeap)
          xmlSimStartTime = HS_xmltime(SimStartHour,BaseYear,useLeap)
          MR_Comp_StartHour     = SimStartHour
          MR_Comp_Time_in_hours = Simtime_in_hours
        endif
        e_StartTime(i) = HS_hours_since_baseyear(iyear(i),imonth(i),  &
                                iday(i),hour(i),BaseYear,useLeap) - SimStartHour
        !error trap if eruptions are not in chronological order
        if(.not.IsCustom_SourceType)then
          ! relax the chronological requirement for custom sources
          if (i.ge.2)then
            !(add 0.001 hours to make sure that rounding error does not cause
            !the program to stop)
            if((e_StartTime(i)+0.001_ip).lt.(e_StartTime(i-1)+e_Duration(i-1))) goto 1912
          endif
        endif
      enddo

      ! write out eruption information
      write(global_info,7) neruptions
      write(global_log ,7) neruptions
      do i=1,neruptions
        write(global_info,8) i, e_PlumeHeight(i), iyear(i), imonth(i), &
                   iday(i), hour(i), e_Duration(i), e_Volume(i)
        write(global_log ,8) i, e_PlumeHeight(i), iyear(i), imonth(i), &
                   iday(i), hour(i), e_Duration(i), e_Volume(i)
        write(global_info,*)"  e_StartTime = ",1,&
                            BaseYear,useLeap,          &
                            real(SimStartHour,kind=4), &
                            real(e_StartTime(i),kind=4)
      if(SourceType.eq.'profile')then
        do ii=1,e_prof_zpoints(i)
          write(global_info,*)"         ",ii,real((ii-1)*e_prof_dz(i),kind=4),&
                                             real(e_prof_Volume(i,ii),kind=4)
        enddo
        write(global_info,*)"         Total Volume for this pulse = ",&
                            real(sum(e_prof_Volume(i,:)),kind=4),"km3 DRE"
      endif
      enddo
      write(global_info,*)"Total volume of all eruptions = ",&
                          real(sum(e_volume),kind=sp),"km3 DRE"
      ! Now that we know the requested dz profile and the plume heights, we can
      ! set up the z-grid for computation
      CompGrid_height = ZPADDING*maxval(e_PlumeHeight(1:neruptions))
!      if(CompGrid_height.lt.10.0_ip*dz_const)then
!        write(global_info,*)&
!         "  Highest source = ", maxval(e_PlumeHeight(1:neruptions))
!        write(global_info,*)&
!         "        ZPADDING = ",ZPADDING
!        write(global_info,*)&
!         " CompGrid_height = ",CompGrid_height
!        write(global_info,*)&
!         "This height for the grid is less than 10 * dz"
!        write(global_info,*)&
!         "Resetting grid height to ",10.0_ip*dz_const
!        CompGrid_height = max(10.0_ip*dz_const,CompGrid_height)
!      endif
      nzmax = 0
      do k = 1,nz_init-1
        if(z_vec_init(k+1).gt.CompGrid_height.and. &
           z_vec_init(k).le.CompGrid_height)then
          nzmax = k
        endif
      enddo
      if(nzmax.eq.0)then
        write(global_error,*)"ERROR: ",&
         "Specified z-grid does not extend high enough"
        write(global_info,*)"        for given plume heights."
        write(global_info,*)"    e_PlumeHeight = ",e_PlumeHeight(1:neruptions)
        write(global_info,*)"         ZPADDING = ",ZPADDING
        write(global_info,*)"  CompGrid_height = ",CompGrid_height
        write(global_info,*)"       z_vec_init = "
        do k = 1,nz_init-1
          write(global_info,*)"                ",k,real(z_vec_init(k),kind=4)
        enddo
        stop 1
      endif

      !for umbrella clouds, 
      itop  = 0
      ibase = 0
      if ((SourceType.eq.'umbrella')    .or. &
          (SourceType.eq.'umbrella_air')) then
        !set nodes at base and top  of cloud
        !itop  = int(e_PlumeHeight(1)/dz)+1           !node at top of plume
        !ibase = int(0.75*e_PlumeHeight(1)/dz)        !node at base
        do k = 1,nzmax
          if(z_vec_init(k)  .ge.e_PlumeHeight(1).and.&
             z_vec_init(k-1).lt.e_PlumeHeight(1))then
            itop = k
          endif
          ! The Suzuki distribution for umbrella clouds is hardwired to k=12
          ! The 0.75 for the plume bottom just ensures that the significant part
          ! of the Suzuki mass loading profile will be within the ibase to itop nodes
          if(z_vec_init(k)  .ge.0.75_ip*e_PlumeHeight(1).and. &
             z_vec_init(k-1).lt.0.75_ip*e_PlumeHeight(1))then
            ibase = k
          endif
        enddo
      endif

      write(global_info,*)
      write(global_log ,*)

      ! Error-check initial wind info and allocate windfile arrays
      MR_BaseYear = BaseYear
      MR_useLeap  = useLeap
      ! Over-ride MetReader output units
      MR_global_essential  = global_essential
      MR_global_production = global_production
      MR_global_debug      = global_debug
      MR_global_info       = global_info
      MR_global_log        = global_log
      MR_global_error      = global_error
      call MR_Allocate_FullMetFileList(iw,iwf,igrid,idf,iwfiles)
      if(MR_useLeap.neqv.useLeap)then
        useLeap  = MR_useLeap
        BaseYear = MR_BaseYear
        write(global_info,*)"Change in calandar; resetting e_StartTime"

          tmp_dp = HS_hours_since_baseyear(iyear(1),imonth(1),  &
                           iday(1),hour(1),BaseYear,useLeap)
          tmp_dp = tmp_dp - SimStartHour   ! Recast tmp_dp as the difference in calandars
          SimStartHour = SimStartHour + tmp_dp
          xmlSimStartTime = HS_xmltime(SimStartHour,BaseYear,useLeap)
      endif
      if (SourceType.eq.'suzuki') then
         write(global_info,6) diffusivity_horz, Suzuki_A, &
                              StopWhenDeposited, Simtime_in_hours
         write(global_log ,6) diffusivity_horz, Suzuki_A, &
                              StopWhenDeposited, Simtime_in_hours
       else
         write(global_info,1438) diffusivity_horz, &
                                 StopWhenDeposited, Simtime_in_hours
         write(global_log ,1438) diffusivity_horz, &
                                 StopWhenDeposited, Simtime_in_hours
      endif
      write(global_info,1439) SourceType
      write(global_log ,1439) SourceType

!     CALCULATE MASS FLUX AND END TIMES OF EACH ERUPTIVE PULSE      
      do i=1,neruptions                            
             !mass flux in kg/hr
        if(SourceType.eq.'suzuki'      .or. &
           SourceType.eq.'point'       .or. &
           SourceType.eq.'line'        .or. &
           SourceType.eq.'umbrella'    .or. &
           SourceType.eq.'umbrella_air')then
          MassFlux(i)  = MagmaDensity * & ! kg/m3
                         e_Volume(i)  * & ! km3
                         KM3_2_M3     / & ! m3/km3
                         e_Duration(i)    ! hours  => kg/hr
          e_EndTime(i) = e_StartTime(i) + e_Duration(i)

          !for debugging >>>>>
          write(6,1023) MagmaDensity, e_Duration(i), MassFlux(i), e_Volume(i)
1023      format('   Magma density (kg/m3) = ',f6.1,', e_Duration(1) (hrs) = ',f6.3,/, &
                 '   Mass flux (kg/hr) = ',e12.4,', Total volume (km3 DRE)=',f8.4,/,'Continue?')
          !read(5,'(a1)') answer
          !if (answer.eq.'n') stop
          !for debugging <<<<
        elseif(SourceType.eq.'profile')then
          e_prof_MassFlux(i,1:e_prof_zpoints(i)) = &
                         MagmaDensity  * &                        ! kg/m3
                         e_prof_Volume(i,1:e_prof_zpoints(i)) * & ! km3
                         KM3_2_M3 / &                             ! m3/km3
                         e_Duration(i)                            ! hours = kg/hr
          MassFlux(i) = sum(e_prof_MassFlux(i,1:e_prof_zpoints(i)))
        else
          ! Custom source, initializing MassFlux and end time
          MassFlux(i)  = 0.0_ip
          e_EndTime(i) = 0.0_ip
        endif

      enddo
      e_EndTime_final = maxval(e_EndTime)  ! this marks the end of all eruptions

!     FIND THE I AND J VALUES OF THE NODE WHERE THE VOLCANO LIES
      if (IsLatLon) then
        ivent = int((lon_volcano-lonLL)/de) + 1
        jvent = int((lat_volcano-latLL)/dn) + 1
      else
        ivent = int((x_volcano-xLL)/dx) + 1
        jvent = int((y_volcano-yLL)/dy) + 1
      endif

      write(global_info,104) ivent,jvent
      write(global_log ,104) ivent,jvent
104   format(4x,'i and j coordinates of volcano:',/, &
             4x,'i=',i4,/, &
             4x,'j=',i4,/)

      ! Now back to reading the input file
      !************************************************************************
      ! BLOCK 4: OUTPUT OPTIONS
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      if (testkey.ne.'#'.and.testkey.ne.'*')then
        write(global_error,*)"ERROR: ",&
         "Expecting a comment line separating blocks."
        write(global_error,*)'       Check that Block 3 is correct.'
        stop 1
      endif
      do while(testkey.eq.'#'.or.testkey.eq.'*')
         ! Line is a comment, read next line
        read(10,'(a80)')linebuffer080
        read(linebuffer080,*)testkey
      enddo
      ! Block 4 Line 1
      write(global_info,*)' *******************************************'
      write(global_info,*)' Reading Block 4: Output options '
      write(global_info,*)' *******************************************'
      WriteDepositFinal_ASCII       = .false.
      WriteDepositFinal_KML         = .false.
      WriteDepositTS_ASCII          = .false.
      WriteDepositTS_KML            = .false.
      WriteCloudConcentration_ASCII = .false.
      WriteCloudConcentration_KML   = .false.
      WriteCloudHeight_ASCII        = .false.
      WriteCloudHeight_KML          = .false.
      WriteCloudLoad_ASCII          = .false.
      WriteCloudLoad_KML            = .false.
      WriteDepositTime_ASCII        = .false.
      WriteDepositTime_KML          = .false.
      WriteCloudTime_ASCII          = .false.
      WriteCloudTime_KML            = .false.
      WriteAirportFile_ASCII        = .false.
      WriteAirportFile_KML          = .false.
      Write_PT_Data                 = .false.
      Write_PR_Data                 = .false.
      cdf_b4l1 = linebuffer080
      ! Read whether to write out final ESRI ASCII deposit file
      read(linebuffer080,*)testkey
      if (testkey.eq.'#'.or.testkey.eq.'*') then
        write(global_error,*)"ERROR: ",&
         "Trying to read WriteDepositFinal_ASCII and detecting",& 
         " comment line"
        stop 1
      endif
      read(linebuffer080,'(a3)',err=1953) answer
      if (answer.eq.'yes') then
        WriteDepositFinal_ASCII = .true.
       else if (answer(1:2).eq.'no') then
        WriteDepositFinal_ASCII = .false.
       else
        goto 1953
      endif

      ! Block 4 Line 2
!     Read whether to write out final KML deposit file
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      if (testkey.eq.'#'.or.testkey.eq.'*') then
        write(global_error,*)"ERROR: ",&
         "Trying to read WriteDepositFinal_KML and detecting", &
         " comment line"
        stop 1
      endif
      read(linebuffer080,'(a3)',err=1954) answer
      cdf_b4l2 = linebuffer080
      if (answer.eq.'yes') then
        WriteDepositFinal_KML = .true.
       else if (answer(1:2).eq.'no') then
        WriteDepositFinal_KML = .false.
       else
        goto 1954
      endif

      ! Block 4 Line 3
!     Read whether to write out ESRI ASCII deposit files at specified times
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      if (testkey.eq.'#'.or.testkey.eq.'*') then
        write(global_error,*)"ERROR: ",&
         "Trying to read WriteDepositTS_ASCII and detecting",&
         " comment line"
        stop 1
      endif
      cdf_b4l3 = linebuffer080
      read(linebuffer080,'(a3)',err=1955) answer
      if (answer.eq.'yes') then
        WriteDepositTS_ASCII = .true.
       else if (answer(1:2).eq.'no') then
        WriteDepositTS_ASCII = .false.
       else
        goto 1955
      endif

      ! Block 4 Line 4
!     Read whether to write out KML deposit files at specified times
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      if (testkey.eq.'#'.or.testkey.eq.'*') then
        write(global_error,*)"ERROR: ",&
         "Trying to read WriteDepositTS_KML and detecting comment line"
        stop 1
      endif
      cdf_b4l4 = linebuffer080
      read(linebuffer080,'(a3)',err=1956) answer
      if (answer.eq.'yes') then
        WriteDepositTS_KML = .true.
       else if (answer(1:2).eq.'no') then
        WriteDepositTS_KML = .false.
       else
        goto 1956
      endif

      ! Block 4 Line 5
!     Read whether to write out ESRI ASCII files of cloud concentration
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      if (testkey.eq.'#'.or.testkey.eq.'*') then
        write(global_error,*)"ERROR: ",&
         "Trying to read WriteCloudConcentration_ASCII and detecting",&
         " comment line"
        stop 1
      endif
      cdf_b4l5 = linebuffer080
      read(linebuffer080,'(a3)',err=1957) answer
      if (answer.eq.'yes') then
        WriteCloudConcentration_ASCII = .true.
       else if (answer(1:2).eq.'no') then
        WriteCloudConcentration_ASCII = .false.
       else
        goto 1957
      endif

      ! Block 4 Line 6
!     Read whether to write out KML files of cloud concentration
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      if (testkey.eq.'#'.or.testkey.eq.'*') then
        write(global_error,*)"ERROR: ",&
         "Trying to read WriteCloudConcentration_KML and detecting",&
         " comment line"
        stop 1
      endif
      cdf_b4l6 = linebuffer080
      read(linebuffer080,'(a3)',err=1958) answer
      if (answer.eq.'yes') then
        WriteCloudConcentration_KML = .true.
       else if (answer(1:2).eq.'no') then
        WriteCloudConcentration_KML = .false.
       else
        goto 1958
      endif

      ! Block 4 Line 7
!     Read whether to write out ASCII files of cloud height
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      if (testkey.eq.'#'.or.testkey.eq.'*') then
        write(global_error,*)"ERROR: ",&
         "Trying to read WriteCloudHeight_ASCII and detecting",&
         " comment line"
        stop 1
      endif
      cdf_b4l7 = linebuffer080
      read(linebuffer080,'(a3)',err=1959) answer
      if (answer.eq.'yes') then
        WriteCloudHeight_ASCII = .true.
       else if (answer(1:2).eq.'no') then
        WriteCloudHeight_ASCII = .false.
       else
        goto 1959
      endif

      ! Block 4 Line 8
!     Read whether to write out KML files of cloud height
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      if (testkey.eq.'#'.or.testkey.eq.'*') then
        write(global_error,*)"ERROR: ",&
         "Trying to read WriteCloudHeight_KML and detecting",&
         " comment line"
        stop 1
      endif
      cdf_b4l8 = linebuffer080
      read(linebuffer080,'(a3)',err=1960) answer
      if (answer.eq.'yes') then
        WriteCloudHeight_KML = .true.
       else if (answer(1:2).eq.'no') then
        WriteCloudHeight_KML = .false.
       else
        goto 1960
      endif

      ! Block 4 Line 9
!     Read whether to write out ASCII files of ashcloud load
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      if (testkey.eq.'#'.or.testkey.eq.'*') then
        write(global_error,*)"ERROR: ",&
         "Trying to read WriteCloudLoad_ASCII and detecting",&
         " comment line"
        stop 1
      endif
      cdf_b4l9 = linebuffer080
      read(linebuffer080,'(a3)',err=19601) answer
      if (answer.eq.'yes') then
        WriteCloudLoad_ASCII = .true.
       else if (answer(1:2).eq.'no') then
        WriteCloudLoad_ASCII = .false.
       else
        goto 19601
      endif

      ! Block 4 Line 10
!     Read whether to write out KML files of ashcloud load
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      if (testkey.eq.'#'.or.testkey.eq.'*') then
        write(global_error,*)"ERROR: ",&
         "Trying to read WriteCloudLoad_KML and detecting comment line"
        stop 1
      endif
      cdf_b4l10 = linebuffer080
      read(linebuffer080,'(a3)',err=1961) answer
      if (answer.eq.'yes') then
        WriteCloudLoad_KML = .true.
       else if (answer(1:2).eq.'no') then
        WriteCloudLoad_KML = .false.
       else
        goto 1961
      endif

      ! Block 4 Line 11
!     Read whether to write out ASCII file of deposit arrival time
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      if (testkey.eq.'#'.or.testkey.eq.'*') then
        write(global_error,*)"ERROR: ",&
         "Trying to read WriteDepositTime_ASCII and detecting",&
         " comment line"
        stop 1
      endif
      cdf_b4l11 = linebuffer080
      read(linebuffer080,'(a3)',err=1962) answer
      if (answer.eq.'yes') then
        WriteDepositTime_ASCII = .true.
       else if (answer(1:2).eq.'no') then
        WriteDepositTime_ASCII = .false.
       else
        goto 1962
      endif

      ! Block 4 Line 12
!     Read whether to write out KML files of deposit arrival time
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      if (testkey.eq.'#'.or.testkey.eq.'*') then
        write(global_error,*)"ERROR: ",&
         "Trying to read WriteDepositTime_KML and detecting",&
         " comment line"
        stop 1
      endif
      cdf_b4l12 = linebuffer080
      read(linebuffer080,'(a3)',err=1963) answer
      if (answer.eq.'yes') then
        WriteDepositTime_KML = .true.
       else if (answer(1:2).eq.'no') then
        WriteDepositTime_KML = .false.
       else
        goto 1963
      endif

      ! Block 4 Line 13
!     Read whether to write out ASCII file of cloud arrival time
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      if (testkey.eq.'#'.or.testkey.eq.'*') then
        write(global_error,*)"ERROR: ",&
         "Trying to read WriteCloudTime_ASCII and detecting",&
         " comment line"
        stop 1
      endif
      cdf_b4l13 = linebuffer080
      read(linebuffer080,'(a3)',err=19621) answer
      if (answer.eq.'yes') then
        WriteCloudTime_ASCII = .true.
       else if (answer(1:2).eq.'no') then
        WriteCloudTime_ASCII = .false.
       else
        goto 19621
      endif

      ! Block 4 Line 14
!     Read whether to write out KML files of cloud arrival time
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      if (testkey.eq.'#'.or.testkey.eq.'*') then
        write(global_error,*)"ERROR: ",&
         "Trying to read WriteCloudTime_KML and detecting comment line"
        stop 1
      endif
      cdf_b4l14 = linebuffer080
      read(linebuffer080,'(a3)',err=19631) answer
      if (answer.eq.'yes') then
        WriteCloudTime_KML = .true.
       else if (answer(1:2).eq.'no') then
        WriteCloudTime_KML = .false.
       else
        goto 19631
      endif

      !If IsLatLon=false, KML files can't be written out until they're re-projected.
      !if ((WriteDepositFinal_KML.or.WriteDepositTS_KML.or.WriteCloudConcentration_KML).and. &
      !    (.not.IsLatLon)) then
      !  write(global_info,38)
      !  write(global_log ,38)
      !  write(global_info,39)
      !  read (5,'(a1)') answer
      !  if (answer.ne.'y') stop 1
      !  WriteCloudConcentration_KML = .false.
      !  WriteDepositFinal_KML    = .false.
      !  WriteDepositTS_KML   = .false.
      !endif
      
      ! Block 4 Line 15
!     Read whether to write out 3D files of ash concentration
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      if (testkey.eq.'#'.or.testkey.eq.'*') then
        write(global_error,*)"ERROR: ",&
         "Trying to read Write3dFiles and detecting comment line"
        stop 1
      endif
      cdf_b4l15 = linebuffer080
      read(linebuffer080,'(a3)',err=1964) answer
      if (answer.eq.'yes') then
        Write3dFiles = .true.
        ! if a consolodate output file will be written, assume both standard
        ! variables and the 3d ash concentrations will be written
        USE_OUTPROD_VARS = .true.
        USE_RESTART_VARS = .true.

        ! Try to read an output code
        loc = index(linebuffer080,'yes')
        dumstr20 = linebuffer080(loc+4:loc+24)
        read(dumstr20,*,iostat=ioerr) iform
        if (ioerr.eq.0)then
          ! Succeeded in reading the format code
          if(iform.eq.1)then
            write(global_info,*)" Successfully read format code=1"
            write(global_info,*)"  Both output products and ash concentrations will be written"
            USE_OUTPROD_VARS = .true.
            USE_RESTART_VARS = .true.
          elseif(iform.eq.2)then
            write(global_info,*)" Successfully read format code=2"
            write(global_info,*)"  Only output products will be written"
            USE_OUTPROD_VARS = .true.
            USE_RESTART_VARS = .false.
          else
            write(global_info,*)" Could not read format code"
            write(global_info,*)"  Assuming both output products and ash concentration will be written"
            USE_OUTPROD_VARS = .true.
            USE_RESTART_VARS = .true.
          endif
        else

        endif
       else if (answer(1:2).eq.'no') then
        Write3dFiles = .false.
       else
        goto 1964
      endif

      ! Block 4 Line 16
!     Read output file format
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      if (testkey.eq.'#'.or.testkey.eq.'*') then
        write(global_error,*)"ERROR: ",&
         "Trying to read output format and detecting comment line"
        stop 1
      endif
      cdf_b4l16 = linebuffer080
      if (Write3dFiles) then
         read(linebuffer080,'(a6)',err=1965) formatanswer
         if (formatanswer(1:5).eq.'ascii') then
            ioutputFormat = 1
          else if (formatanswer(1:6).eq.'binary') then
            ioutputFormat = 2
          else if (formatanswer(1:6).eq.'netcdf') then
            ioutputFormat = 3
          else
           goto 1965
         endif
      endif
      
      ! Block 4 Line 17
!     Read number of files to write out
      read(10,'(a80)') linebuffer080
      read(linebuffer080,*)testkey
      if (testkey.eq.'#'.or.testkey.eq.'*') then
        write(global_error,*)"ERROR: ",&
         "Trying to read nWriteTimes and detecting comment line"
        stop 1
      endif
      cdf_b4l17 = linebuffer080

      ! Block 4 Line 18
      read(10,'(a130)')linebuffer130
      read(linebuffer130,*)testkey
      if (testkey.eq.'#'.or.testkey.eq.'*') then
        write(global_error,*)"ERROR: ",&
         "Trying to read WriteTimes or WriteInterval and detecting",&
         " comment line"
        stop 1
      endif
      cdf_b4l18 = linebuffer130(1:80)
      if (WriteDepositFinal_ASCII      .or. &
          WriteDepositFinal_KML        .or. &
          WriteDepositTS_ASCII         .or. &
          WriteDepositTS_KML           .or. &
          WriteCloudConcentration_ASCII.or. &
          WriteCloudConcentration_KML  .or. &
          Write3dFiles                 .or. &
          WriteCloudHeight_ASCII       .or. &
          WriteCloudHeight_KML         .or. &
          WriteCloudLoad_KML) then
          read(cdf_b4l17,*,err=1966) nWriteTimes
            ! Check how to interpret nWriteTimes
          if (nWriteTimes.gt.0) then
            ! If a positive number, then we're reading an array of times
            allocate(WriteTimes(nWriteTimes))
            read(linebuffer130,*,err=1965) WriteTimes(1:nWriteTimes)
          elseif (nWriteTimes.eq.0) then
            write(global_info,*)&
             "nWriteTimes = 0: Running without output"
          elseif (nWriteTimes.ne.-1) then
            ! If not a positive number, then it should be -1
            ! Report error otherwise
            goto 1966
          else
            !If -1, then read a single WriteTimes and interpret it as a time interval
            read(cdf_b4l18,*,err=1967) WriteInterval
              ! Redefine nWriteTimes since it was read in as -1
            nWriteTimes = int(Simtime_in_hours/WriteInterval)+1
            write(global_info,*)"  WriteInterval    = ",real(WriteInterval,kind=sp)
            write(global_info,*)"  Simtime_in_hours = ",real(Simtime_in_hours,kind=sp)
            write(global_info,*)"  nWriteTimes      = ",nWriteTimes
            allocate(WriteTimes(nWriteTimes))

            forall (i=1:nWriteTimes)  WriteTimes(i) = (i-1)*WriteInterval
            do i=1,nWriteTimes
              WriteTimes(i) = (i-1)*WriteInterval
            enddo
            do i=1,nWriteTimes     !check writetimes for errors
              if (WriteTimes(i).lt.0.0_ip) then  !if the time <0
                  ! Abort the program
                goto 1968
              elseif(i.gt.1)then
                if(WriteTimes(i).lt.Writetimes(i-1))then !if times are not in chronological order
                    ! Abort the program
                  goto 1969
                endif
              elseif (WriteTimes(i).gt.Simtime_in_hours) then   !if some times exceed the simulation time
                write(global_info,32)
                write(global_log ,32)
                nWriteTimes = i-1
                exit
              endif
            enddo
          endif

          if(LoadConcen)then
            ! Find the output time index that is next
            if (time.le.WriteTimes(1))then
              iTimeNext = 1
            else
              do i = 2, nWriteTimes
                if(time.ge.WriteTimes(i-1).and. &
                   time.lt.WriteTimes(i))then
                  iTimeNext = i
                endif
              enddo
            endif
          else
            iTimeNext = 1
          endif
          NextWriteTime = WriteTimes(iTimeNext)
       endif

       !WRITE OUT THE TYPES OF OUTPUT TO BE WRITTEN
        !output options
       write(global_info,33) WriteDepositFinal_ASCII,       &
                   WriteDepositFinal_KML,         &
                   WriteDepositTS_ASCII,          &
                   WriteDepositTS_KML,            &
                   WriteCloudConcentration_ASCII, &
                   WriteCloudConcentration_KML,   &
                   WriteCloudHeight_ASCII,        &
                   WriteCloudHeight_KML,          &
                   WriteCloudLoad_ASCII,          &
                   WriteCloudLoad_KML,            &
                   WriteDepositTime_ASCII,        &
                   WriteDepositTime_KML,          &
                   WriteCloudTime_ASCII,          &
                   WriteCloudTime_KML,            &
                   Write3dFiles,                  &
                   formatanswer,                  &
                   nWriteTimes
       write(global_log ,33) WriteDepositFinal_ASCII,       &
                   WriteDepositFinal_KML,         &
                   WriteDepositTS_ASCII,          &
                   WriteDepositTS_KML,            &
                   WriteCloudConcentration_ASCII, &
                   WriteCloudConcentration_KML,   &
                   WriteCloudHeight_ASCII,        &
                   WriteCloudHeight_KML,          &
                   WriteCloudLoad_ASCII,          &
                   WriteCloudLoad_KML,            &
                   WriteDepositTime_ASCII,        &
                   WriteDepositTime_KML,          &
                   WriteCloudTime_ASCII,          &
                   WriteCloudTime_KML,            &
                   Write3dFiles,                  &
                   formatanswer,                  &
                   nWriteTimes

      if (WriteDepositTS_ASCII          .or. &
          WriteDepositTS_KML            .or. &
          WriteCloudConcentration_ASCII .or. &
          WriteCloudConcentration_KML   .or. &
          WriteCloudHeight_ASCII        .or. &
          WriteCloudHeight_KML          .or. &
          WriteCloudLoad_ASCII          .or. &
          WriteCloudLoad_KML            .or. &
          WriteCloudTime_ASCII          .or. &
          WriteCloudTime_KML            .or. &
          Write3dFiles) then
          write(global_info,34)                                !write out the types of output specified
          write(global_log ,34)
          do i=1,nWriteTimes
            write(global_info,35) WriteTimes(i)                !write out the times when  files will be written
            write(global_log ,35) Writetimes(i)
          enddo
       endif
      ! END OF BLOCK 4
      !************************************************************************

      !************************************************************************
      ! BLOCK 5: INPUT WIND FILES
      if(MR_iwindfiles.gt.0)then
        read(10,'(a130)')linebuffer130
        read(linebuffer130,*)testkey
        if (testkey.ne.'#'.and.testkey.ne.'*')then
          write(global_error,*)"ERROR: ",&
           "Expecting a comment line separating blocks."
          write(global_error,*)'       Check that Block 4 is correct.'
          stop 1
        endif      
        do while(testkey.eq.'#'.or.testkey.eq.'*')
           ! Line is a comment, read next line
          read(10,'(a130)')linebuffer130
          ! Normally, we would read the first character of the string linebuffer130, but
          ! this seems to fail if the first character is '/'
          !read(linebuffer130,*)testkey
          testkey=linebuffer130(1:1)
        enddo
        write(global_info,*)' *****************************************'
        write(global_info,*)' Reading Block 5: Windfile names'
        write(global_info,*)' *****************************************'
        write(global_info,13)
        write(global_log ,13)
          ! Read list of windfiles.
        if(MR_iwind.eq.5)then
          ! For NCEP 2.5 degree (25), NOAA product (27), ERA5 (29), or ERA-20C (30)
          ! just read the path to the files
          read(linebuffer130,'(a130)',err=1970) MR_windfiles(1)
          write(global_info,1034) 1,trim(adjustl(MR_windfiles(1)))
          read(10,'(a130)')linebuffer130
        else
          ! For all other iwf (MR_iwindformats), read the full list
          do i=1,iwfiles
            ! Always check if we have overshot the block
            testkey=linebuffer130(1:1)
            if (testkey.eq.'#'.or.testkey.eq.'*') then
              write(global_log,*)"ERROR: ",&
               "Trying to read Block 5 and detecting comment line"
              write(global_log,*)  '  Windfile ',i,'of',iwfiles
              write(global_log,*)  '  Offending line:',linebuffer130
              write(global_error,*)"ERROR: ",&
               "Trying to read Block 5 and detecting comment line"
              write(global_error,*)'  Windfile ',i,'of',iwfiles
              write(global_error,*)'  Offending line:',linebuffer130
              stop 1
            endif
            read(linebuffer130,'(a130)',err=1970) MR_windfiles(i)
            write(global_info,1034) i,trim(adjustl(MR_windfiles(i)))
            if(idf.eq.3)then
              ! If we are reading grib files, check that the index file has been
              ! generated
#ifndef USEGRIB
              write(global_error,*)"ERROR: ",&
               "This input file specifies that the Met files are"
              write(global_error,*)&
               "       in grib format, but Ash3d has not been compiled"
              write(global_error,*)&
               "       with grib support.  Please recompile with grib"
              write(global_error,*)&
               "       enabled.  MetReader must also support grib."
              stop 1
#else
              MR_windfiles_GRIB_index(i) = trim(adjustl(MR_windfiles(i))) // ".index"
              inquire( file=MR_windfiles_GRIB_index(i), exist=IsThere )
              MR_windfiles_Have_GRIB_index(i) = IsThere
              if(.not.IsThere)then
                ! Grib index file is not there, Try to generate it.
                ! Note, we might have some permission problems here and we should
                ! set up a fail-safe to the cwd or something.
                write(global_info,*)&
                 " Grib index file not found; attempting to create it."
                call MR_Set_Gen_Index_GRIB(MR_windfiles(i))
              endif
#endif
            endif
            read(10,'(a130)')linebuffer130
          enddo
        endif
1034    format(' i=',i3,'  MR_windfiles(i) = ',a)
      else
        write(global_error,*)"ERROR: ",&
         "MR_iwindfiles = 0"
        write(global_error,*)&
         "       Either the number of windfiles specified = 0, or"
        write(global_error,*)&
         "       MR_Allocate_FullMetFileList has not been called."
        stop 1
      endif

      !Error trap if more windfiles are entered than are specified
      if (linebuffer130(1:5).ne.'*****') then
        write(global_info,*)&
         'The beginning of the line following the list of', &
         ' windfiles did not'
        write(global_info,*)&
         'start with ''*****''.  Did you enter the correct', &
         '  number of windfiles?'
        write(global_info,*) 'Program stopped.'
        stop 1
      endif
      ! END OF BLOCK 5
      !************************************************************************

        ! Check for existance and compatibility with simulation time requirements
      call MR_Read_Met_DimVars(iyear(1))
        ! Now that we have the actual times available from the Met files, we can reset
        ! the Simulation Start times for forecast runs
      if(runAsForecast)then
        MR_Comp_StartHour = MR_windfile_starthour(1) + MR_windfile_stephour(1,1) + FC_Offset
        SimStartHour      = MR_Comp_StartHour
        xmlSimStartTime   = HS_xmltime(SimStartHour,BaseYear,useLeap)
      endif

      !************************************************************************
      ! BLOCK 6: AIRPORT FILE
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      if (testkey.ne.'#'.and.testkey.ne.'*')then
        write(global_error,*)"ERROR: ",&
         "Expecting a comment line separating blocks."
        write(global_error,*)'       Check that Block 5 is correct.'
        stop 1
      endif
      do while(testkey.eq.'#'.or.testkey.eq.'*')
         ! Line is a comment, read next line
        read(10,'(a80)')linebuffer080
        read(linebuffer080,*)testkey
      enddo
      write(global_info,*)' *******************************************'
      write(global_info,*)' Reading Block 6: Airport location output'
      write(global_info,*)' *******************************************'
      !Read whether to write out ASCII airport file
      ! Block 6 Line 1
      cdf_b6l1 = linebuffer080
      read(linebuffer080,'(a3)',err=1980) answer
      if (answer.eq.'yes') then
        WriteAirportFile_ASCII = .true.
        Write_PT_Data          = .true.
       else if (answer(1:2).eq.'no') then
        WriteAirportFile_ASCII = .false.
       else
        goto 1980
      endif
      
      !Read whether to write out grain-size distribution to airport file
      read(10,'(a80)') linebuffer080
      ! Block 6 Line 2
      cdf_b6l2 = linebuffer080
      read(linebuffer080,'(a3)',err=1981) answer
      if (answer.eq.'yes') then
        WriteGSD = .true.
       else if (answer(1:2).eq.'no') then
        WriteGSD = .false.
       else
        goto 1981
      endif

      !Read whether to write out kml airport file
      read(10,'(a80)') linebuffer080
      ! Block 6 Line 3
      cdf_b6l3 = linebuffer080
      read(linebuffer080,'(a3)',err=1982) answer
      if (answer.eq.'yes') then
        WriteAirportFile_KML = .true.
        Write_PT_Data        = .true.
       else if (answer(1:2).eq.'no') then
        WriteAirportFile_KML = .false.
       else
        goto 1982
      endif
            
      !Read name of input file containing airport locations
      ! Block 6 Line 4
      read(10,'(a80)') cdf_b6l4
      AirportInFile = cdf_b6l4(1:scan(cdf_b6l4,' ')-1)     !Read to the first blank space      

      !See if we need to read an external airport file
      if ((AirportInFile.ne.'internal').and. &
          (AirportInFile.ne.'')) then
        ReadExtAirportFile=.true.              !read external data, do not append
        if (AirportInFile(1:1).eq.'+') then
          AppendExtAirportFile=.true.       !read and append external data
          AirportInFile = AirportInFile(2:)      !strip off the "plus" at the beginning
        else
          AppendExtAirportFile=.false.       !read and append external data
        endif
        !Make sure the external file exists and can be opened.
        write(global_info,*)&
        'Making sure the external airport file exists and can be opened'
        open(unit=17,file=AirportInFile,status='old',err=19825) !try opening the external file
        close(17)                              !if it opens, close it back up.
      else
         ReadExtAirportFile = .false.              !read external data, do not append
      endif

      !Read whether to project airport coordinates
19828 continue
      ! Block 6 Line 5
      read(10,'(a80)') cdf_b6l5
      read(cdf_b6l5,'(a3)',err=1983) answer
      if (answer.eq.'yes') then
        ProjectAirportLocations = .true.
       else if (answer(1:2).eq.'no') then
        ProjectAirportLocations = .false.
       else
        goto 1983
      endif

      !Write out parameters
      write(global_info,44) ReadExtAirportFile, AppendExtAirportFile, &
                  WriteAirportFile_ASCII, WriteGSD, WriteAirportFile_KML, &
                  ProjectAirportLocations, AirportInFile
      write(global_log ,44) ReadExtAirportFile, AppendExtAirportFile, &
                  WriteAirportFile_ASCII, WriteGSD, WriteAirportFile_KML, &
                  ProjectAirportLocations, AirportInFile
      ! END OF BLOCK 6
      !************************************************************************

      !************************************************************************
      ! BLOCK 7: GRAIN-SIZE BINS, SETTLING VELOCITY
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      if (testkey.ne.'#'.and.testkey.ne.'*')then
        write(global_error,*)"ERROR: ",&
         "Expecting a comment line separating blocks."
        write(global_error,*)'       Check that Block 6 is correct.'
        stop 1
      endif
      do while(testkey.eq.'#'.or.testkey.eq.'*')
         ! Line is a comment, read next line
        read(10,'(a80)')linebuffer080
        read(linebuffer080,*)testkey
      enddo
      write(global_info,*)' *******************************************'
      write(global_info,*)' Reading Block 7: Grain Size Groups'
      write(global_info,*)' *******************************************'
      ! READ GRAIN-SIZE BINS
      ! First, get the number of tephra bins to read
      ! Note: This might be one greater than what is calculated if the last bin
      !       has a negative diameter.  In this case, the remaining mass fraction
      !       neglecting the last bin is distributed over the previous bins with
      !       a gaussian distribution given by a phi_mean and phi_stdev
      !    e.g.  -1 4 2
      !       Also note that the number of tephra bins can be zero if the species
      !       will be defined in optional modules such as gas, aggregates, etc.
      ! Block 6 Line 1
      read(linebuffer080,*,iostat=ioerr) ivalue1
      init_n_gs_max = ivalue1
      read(linebuffer080,*,iostat=ioerr) ivalue1, ivalue2
      ! Assume we can read at least read one value, try for two with the second being
      ! the fall model:
      !  FV_ID = 0 -> No fall (just tracer)
      !          1 -> Wilson and Huang
      !          2 -> Wilson and Huang + Cunningham slip
      !          3 -> Wilson and Huang + Mod by Pfeiffer Et al.
      !          4 -> Ganser
      !          5 -> Stokes flow for spherical particles + slip
      if (ioerr.eq.0)then
        FV_ID = ivalue2
      else
        FV_ID = 1
      endif
      allocate(temp_v_s(init_n_gs_max))
      allocate(temp_gsdiam(init_n_gs_max))
      allocate(temp_bin_mass(init_n_gs_max))
      allocate(temp_rho_m(init_n_gs_max))
      allocate(temp_gsF(init_n_gs_max))
      allocate(temp_gsG(init_n_gs_max))

      if(init_n_gs_max.gt.0)then
        do i=1,init_n_gs_max
          value1 = -1.99_ip
          value2 = -1.99_ip
          value3 = -1.99_ip
          read(10,'(a80)')linebuffer080
          ! Always check if we have overshot the block
          testkey = linebuffer080(1:1)
          if ((testkey.eq.'*').or.(testkey.eq.'#')) then
            write(global_error,*)"ERROR: ",&
             "Error in specifying grain sizes.  You specified ",&
             init_n_gs_max,', sizes,'
            write(global_error,*) 'but only ',i-1,&
             ' size classes were listed in the input file.'
            write(global_error,*)' Offending line: ',linebuffer080
            write(global_error,*) 'Program stopped'
            write(global_error,*)"ERROR: ",&
             "Error in specifying grain sizes.  You specified ",&
             init_n_gs_max,', sizes,'
            write(global_log ,*) 'but only ',i-1,&
             ' size classes were listed in the input file.'
            write(global_log ,*) 'Program stopped'
            stop 1
          endif
          read(linebuffer080,*,iostat=ioerr) value1, value2
          ! Assume we can read at least read two values, try for three
          if (ioerr.eq.0)then
            read(linebuffer080,*,iostat=ioerr) value1, value2, value3
            if (ioerr.eq.0)then
              ! Three values were successfully read, interpret as:
              ! grain-size, mass fraction, density
              ! W&H suggest 800 kg/m3 for d>300um and 2000 for d<88um for pumice
              ! fragments
              useCalcFallVel = .true. 
              useTemperature = .true. ! When calculating Fall Vel. we need T
              temp_gsdiam(i) = value1
              temp_bin_mass(i) = value2
              temp_rho_m(i) = value3
              ! Try for a forth value for shape
              read(linebuffer080,*,iostat=ioerr) value1, value2, value3, value4
              if (ioerr.eq.0)then
                ! Fourth value was successfully read, interpret as W/H shape
                ! parameter
                temp_gsF(i) = value4
                ! Try for a fifth value for second shape parameter for Ganser
                ! model
                read(linebuffer080,*,iostat=ioerr) value1, value2, value3, value4, value5
                if (ioerr.eq.0)then
                  ! Fourth value was successfully read, interpret as Ganser 2nd
                  ! shape parameter
                  temp_gsG(i) = value5
                else
                  temp_gsG(i) = 1.0_ip
                endif
              else
                temp_gsF(i) = 0.44_ip
              endif
                ! Initialize this to zero
              temp_v_s(i) = 0.0_ip
              if(temp_gsdiam(i).lt.0.0_ip)then
                if(i.lt.init_n_gs_max)then
                  write(global_error,*)"ERROR: ",&
                   "diameter must be positive"
                  stop 1
                else
                  phi_mean   = value2
                  phi_stddev = value3
                  write(global_info,*) &
          "Last grain-size bin will be partitioned across all previous."
                  write(global_info,*)"Volume fraction partitioned = ",&
          1.0_ip-sum(temp_bin_mass(1:init_n_gs_max-1))
                  write(global_info,*)&
          "  Assuming remainder is Gaussian in phi"
                  write(global_info,*)&
          "    phi_mean   = ", phi_mean
                  write(global_info,*)&
          "    phi_stddev = ", phi_stddev
                  useVariableGSbins = .true.
                endif
              endif
            else
              ! Only two values were successfully read, interpret with
              ! old format as:
              ! FallVel, mass fraction
              useCalcFallVel = .false.
              temp_v_s(i)     = value1
              if (temp_v_s(i).lt.0.0_ip)then
                write(global_info,*)&
                 "WARNING: fall velocity is negative.  Grains will 'fall' upward"
              endif
              temp_bin_mass(i)= value2
                ! Initialize these
              temp_gsdiam(i)  = 0.1_ip
              temp_rho_m(i)   = 2000.0_ip
              temp_gsF(i)     = 0.44_ip
            endif
          endif
        enddo
        ! Set the number of grain-size bins
        if (useVariableGSbins)then
          ! In this case, the last bin is the remainder to be distributed over the
          ! previous bins.  So we need to decrement the tephra bins by 1
          n_gs_max = init_n_gs_max-1
        else
          n_gs_max = init_n_gs_max
        endif
      endif ! if init_n_gs_max > 0
      ! END OF BLOCK 7
      !************************************************************************

      ! Since this subroutine is called before any optional modules, we can
      ! initialize nsmax to the number of tephra bins
      nsmax      = n_gs_max  ! Total tracked bins
      n_gs_aloft = n_gs_max  ! Number of tephra species aloft

      call Allocate_Tephra
      allocate(temp_phi(n_gs_max))

      Tephra_v_s(1:n_gs_max)      = -1.0_ip * temp_v_s(1:n_gs_max) ! make sure 'fall velocity'
                                                                   ! is in the -z direction
      Tephra_gsdiam(1:n_gs_max)   = temp_gsdiam(1:n_gs_max)
      Tephra_bin_mass(1:n_gs_max) = temp_bin_mass(1:n_gs_max)
      Tephra_rho_m(1:n_gs_max)    = temp_rho_m(1:n_gs_max)
      Tephra_gsF(1:n_gs_max)      = temp_gsF(1:n_gs_max)
      Tephra_gsG(1:n_gs_max)      = temp_gsG(1:n_gs_max)

      deallocate(temp_v_s,temp_gsdiam,temp_bin_mass,temp_rho_m,temp_gsF)

      if(n_gs_max.gt.0)then
        call Calculate_Tephra_Shape

        ! If a log-normal distribution is to be added, make sure the grainsize
        ! bins are sorted by fall velocity (slowest first)
        if (useVariableGSbins) call Sort_Tephra_Size

        temp_phi = -log(Tephra_gsdiam)/log(2.0)

        if (useCalcFallVel) Tephra_gsdiam = Tephra_gsdiam/1000.0_ip   ! convert diameter from mm to m

        ! Find the fraction of fine (<= phi4, 63um)
        fracfine = 0.0_ip
        do i=1,n_gs_max
          if(Tephra_gsdiam(i).lt.6.4e-5_ip)then
            fracfine = fracfine + Tephra_bin_mass(i)
          endif
        enddo

!       Make sure that Tephra_bin_mass>0 for  each size class
        write(global_info,2531)
        write(global_log,2531)
2531    format('Checking to make sure Tephra_bin_mass>0 for all size classes')
        do i=1,n_gs_max
           if (Tephra_bin_mass(i).lt.0.0_ip) then
             write(global_info,25311) i
             write(global_log,25311) i
25311        format('Error: mass of bin ',i2,' is less than zero.  Program stopped')
             stop 1
           end if
        end do
!       Send error message if sum(bin_mass(1:n_gs_max)) does not equal 1
        sum_bins=sum(Tephra_bin_mass(1:n_gs_max))
        if(abs(sum_bins-1.0_ip).gt.0.02_ip) then
          write(global_info,2532)
          write(global_log ,2532)
          do i=1,n_gs_max
            write(global_info,2533) i, Tephra_bin_mass(i)
            write(global_log ,2533) i, Tephra_bin_mass(i)
          enddo
          write(global_info,2534) sum_bins
          write(global_log ,2534) sum_bins
2532      format('Error.  Sum of mass fractions of grain-size bins',/, &
                 'does not equal 1.',/, &
                 'bin   mass')
2533      format(i3,f7.4)
2534      format(3x,f7.4,'  total',/,'Program stopped')
          stop 1
            !If it differs just slightly from 1, adjust automatically
        else if (abs(sum_bins-1.0_ip).gt.0.001_ip) then
          write(global_info,2535) sum_bins
          write(global_log ,2535) sum_bins
2535      format('Warning: sum(bin_mass(1:n_gs_max))=',f10.5,/, &
                 'This differs slightly from 1.0',/, &
                 'adjusting bin masses automatically.')
          do i=1,n_gs_max
            Tephra_bin_mass(i) = Tephra_bin_mass(i)*1.0_ip/sum_bins
          enddo
        endif

!       WRITE OUT GRAIN-SIZE BINS:
        write(global_info,9) n_gs_max, FV_ID                               ! write out grain size bins
        write(global_log ,9) n_gs_max, FV_ID
        if(FV_ID.eq.0)then
          write(global_info,*)"Fall Model = None (tracer)"
          write(global_log ,*)"Fall Model = None (tracer)"
        elseif(FV_ID.eq.1)then
          write(global_info,*)"Fall Model = Wilson and Huang"
          write(global_log ,*)"Fall Model = Wilson and Huang"
        elseif(FV_ID.eq.2)then
          write(global_info,*)"Fall Model = Wilson and Huang + Cunningham slip"
          write(global_log ,*)"Fall Model = Wilson and Huang + Cunningham slip"
        elseif(FV_ID.eq.3)then
          write(global_info,*)"Fall Model = Wilson and Huang + Mod by PCM"
          write(global_log ,*)"Fall Model = Wilson and Huang + Mod by PCM"
        elseif(FV_ID.eq.4)then
          write(global_info,*)"Fall Model = Ganser"
          write(global_log ,*)"Fall Model = Ganser"
        elseif(FV_ID.eq.5)then
          write(global_info,*)"Fall Model = Stokes flow + slip"
          write(global_log ,*)"Fall Model = Stokes flow + slip"
        else
          write(global_info,*)"Default Fall Model = Wilson and Huang"
          write(global_log ,*)"Default Fall Model = Wilson and Huang"
        endif
        if(useCalcFallVel)then
          write(global_info,10)
          write(global_log ,10)
          do i=1,n_gs_max
              !write out diameter in mm, not m
            write(global_info,11) Tephra_bin_mass(i), Tephra_gsdiam(i)*1000.0_ip, Tephra_rho_m(i),&
                          Tephra_gsF(i), Tephra_gsG(i), temp_phi(i)
            write(global_log ,11) Tephra_bin_mass(i), Tephra_gsdiam(i)*1000.0_ip, Tephra_rho_m(i),&
                          Tephra_gsF(i), Tephra_gsG(i), temp_phi(i)
          enddo
        else
          write(global_info,2110)
          write(global_log ,2110)
          do i=1,n_gs_max
            write(global_info,2111) Tephra_bin_mass(i), Tephra_v_s(i)
            write(global_log ,2111) Tephra_bin_mass(i), Tephra_v_s(i)
          enddo
        endif
        write(global_info,*)
        write(global_log ,*)
        
        write(global_info,*)"Using deposit density (kg/m3) of: ",DepositDensity
        write(global_log,*)"Using deposit density (kg/m3) of: ",DepositDensity

        write(global_info,*)
        write(global_log ,*)        

        write(global_info,*)"Using a mass-fraction of fines (< 63um) of : ",real(fracfine,kind=sp)
        write(global_log ,*)"Using a mass-fraction of fines (< 63um) of : ",real(fracfine,kind=sp)
        !if bin masses sum up close to 1, adjust automatically      
        if ((abs(sum(Tephra_bin_mass)-1.0_ip).gt.1.0e-5_ip).and. &  
            (abs(sum(Tephra_bin_mass)-1.0_ip).le.1.0e-02_ip)) then
          write(global_info,1001) sum(Tephra_bin_mass)
          write(global_log ,1001) sum(Tephra_bin_mass)
1001      format(4x,'The sum of the bin masses=',f5.3,&
                 '  Adjusting to equal 1',/, &
                 4x,'New bin masses:')
          Tephra_bin_mass(1:n_gs_max) = Tephra_bin_mass(1:n_gs_max)*1.0_ip/sum(Tephra_bin_mass(1:n_gs_max))

          if(useCalcFallVel)then
            do i=1,n_gs_max
              write(global_info,11) Tephra_bin_mass(i), Tephra_gsdiam(i)*1000.0_ip, Tephra_rho_m(i)       !write out diameter in mm, not m
              write(global_log ,11) Tephra_bin_mass(i), Tephra_gsdiam(i)*1000.0_ip, Tephra_rho_m(i)
            enddo
          else
            do i=1,n_gs_max
              write(global_info,2111) Tephra_bin_mass(i), Tephra_v_s(i)
              write(global_log ,2111) Tephra_bin_mass(i), Tephra_v_s(i)
            enddo
          endif
        else if (abs(sum(Tephra_bin_mass)-1.0_ip).gt.1.0e-2_ip) then   !if bin masses are farther off, stop
          write(global_info,1000)
          write(global_log ,1000)
1000      format(4x,'Error: Sum of the mass fraction of the grain',&
                 ' size bins differs from 1 by more than 0.01',/, &
                 4x,'Program stopped')
        endif
      endif

      !************************************************************************
      ! BLOCK 8: VERTICAL PROFILES
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*)testkey
      if (testkey.ne.'#'.and.testkey.ne.'*')then
        write(global_error,*)"ERROR: ",&
         "Expecting a comment line separating blocks."
        write(global_error,*)'       Check that Block 7 is correct.'
        stop 1
      endif
      do while(testkey.eq.'#'.or.testkey.eq.'*')
         ! Line is a comment, read next line
        read(10,'(a80)')linebuffer080
        read(linebuffer080,*)testkey
      enddo
      write(global_info,*)' *******************************************'
      write(global_info,*)' Reading Block 8: Vertical profile output'
      write(global_info,*)' *******************************************'

      !READ NUMBER OF VERTICAL PROFILES
      write(global_info,*) 'Reading vertical profile information'
      read(linebuffer080,*,err=2000) nvprofiles
      write(global_info,*) 'number of vertical profiles=',nvprofiles
      write(global_log ,*) 'number of vertical profiles=',nvprofiles

      if (nvprofiles.gt.0) then
        Write_PR_Data = .true.
        write(global_info,*)"Allocating profile arrays:",nvprofiles
        allocate(x_vprofile(nvprofiles))
        allocate(y_vprofile(nvprofiles))
        allocate(i_vprofile(nvprofiles))
        allocate(j_vprofile(nvprofiles))
        allocate(Site_vprofile(nvprofiles))
        write(global_info,46)
        write(global_log ,46)
46      format(/,'     vertical profile locations',/, &
                  '          #         x         y     i     j')
        do i=1,nvprofiles
          read(10,'(a80)') linebuffer080
          ! Always check if we have overshot the block
          testkey=linebuffer080(1:1)
          if (testkey.eq.'#'.or.testkey.eq.'*') then
            write(global_log,*)"ERROR: ",&
             "Trying to read Block 8 and detecting comment line"
            write(global_log,*)  '  Vert Prof',i,'of',nvprofiles
            write(global_log,*)  '  Offending line:',linebuffer080
            write(global_error,*)"ERROR: ",&
             "Trying to read Block 8 and detecting comment line"
            write(global_error,*)'  Vert Prof ',i,'of',nvprofiles
            write(global_error,*)'  Offending line:',linebuffer080
            stop 1
          endif

          read(linebuffer080,*,iostat=ioerr) value1, value2
          x_vprofile(i) = value1
          y_vprofile(i) = value2
          write(Site_vprofile(i),'(a14,1x,i3)')"Vertical Prof ",i
          ! Assume we can read at least read two values, try for three
          if (ioerr.eq.0)then
            read(linebuffer080,*,iostat=ioerr) value1, value2, Site_vprofile(i)
            if (ioerr.eq.0)then
              ! HFS do some logic here to check if this is just a comment
              substr_pos1 = index(linebuffer080,trim(adjustl(Site_vprofile(i))))
              substr_pos2 = index(linebuffer080,'#')
              if(substr_pos2.eq.0)then
                ! comment indicator '#' not found, set end of string to length
                substr_pos2 = min(len(linebuffer080),substr_pos1+50)
              endif
              Site_vprofile(i) = trim(adjustl(linebuffer080(substr_pos1:substr_pos2)))
              !write(*,*)"Trying to interpret label:",Site_vprofile(i)
              !Site_vprofile(i) = trim(adjustl(Site_vprofile(i)))
            endif
          else
            write(global_info,*) 'Error in x or y location of a vertical profile.'
            write(global_info,*) 'Answer should be two real numbers.'
            write(global_info,*) 'You gave: ',linebuffer080
            write(global_info,*) 'Program stopped.'
            write(global_log ,*) 'Error in x or y location of a vertical profile.'
            write(global_log ,*) 'Answer should be two real numbers.'
            write(global_log ,*) 'You gave: ',linebuffer080
            write(global_log ,*) 'Program stopped.'
            stop 1
          endif
          call vprofchecker(i)
        enddo
      else
        write(global_info,*) 'No vertical profile specified'
      endif
      ! END OF BLOCK 8
      !************************************************************************

      !************************************************************************
      ! BLOCK 9 (Optional): NETCDF ANNOTATIONS
      read(10,'(a80)',iostat=ios)linebuffer080
      read(linebuffer080,*)testkey
      if (testkey.ne.'#'.and.testkey.ne.'*')then
        write(global_error,*)"ERROR: ",&
         "Expecting a comment line separating blocks."
        write(global_error,*)'       Check that Block 8 is correct.'
        stop 1
      endif      
      do while(ios.eq.0.and.(testkey.eq.'#'.or.testkey.eq.'*'))
         ! Line is a comment, read next line
        read(10,'(a80)',iostat=ios)linebuffer080
        read(linebuffer080,*)testkey
      enddo

      ! Here are the default output file name and comments if Block 9 is not given
      outfile = "3d_tephra_fall.nc"
      cdf_title = infile
      cdf_comment = "None"
      cdf_institution="USGS"
      cdf_source="ash3d v1.0b"
      cdf_run_class="Analysis"
      cdf_url="https://vsc-ash.wr.usgs.gov/ash3d-gui"
      cdf_history=""
      cdf_references="https://vsc-ash.wr.usgs.gov/ash3d-gui"
      cdf_conventions='CF-1.5'
      if(ios.ne.0)then
        write(global_info,*)'  Setting outfile to: 3d_tephra_fall.nc'
        write(global_info,*)'  Setting Title to: ',infile
        write(global_info,*)'  Setting comment to: None'
        write(global_info,*)'  Setting run class to: Analysis'
        write(global_info,*)'  Setting institution to: USGS'
      else
        write(global_info,*)' *****************************************'
        write(global_info,*)' Reading Block 9: Output file / comments'
        write(global_info,*)' *****************************************'
        !Start reading annotation info

        ! First line is the output file name
        read(linebuffer080,*) outfile
        outfile = trim(adjustl(outfile))

        ! Next line is the title of the job
          ! Read title line up until the first '#', then truncate
        read(10,'(a80)',iostat=ios,err=2010)linebuffer080
        iendstr = SCAN(linebuffer080, "#")
        if (iendstr.eq.0)then
             ! '#' not found, just copy linebuffer080 to title
          cdf_title = trim(adjustl(linebuffer080))
        else
            ! clip title at key
          cdf_title = trim(linebuffer080(1:iendstr-1))
        endif

          ! Read comment line up until the first '#', then truncate
        read(10,'(a80)',iostat=ios,err=2010)linebuffer080
        if(ios.ne.0)goto 2010
        iendstr = SCAN(linebuffer080, "#")
        if (iendstr.eq.0)then
             ! '#' not found, just copy linebuffer080 to comment
          cdf_comment = trim(adjustl(linebuffer080))
        else
            ! clip comment at key
          cdf_comment = trim(linebuffer080(1:iendstr-1))
        endif

          ! Read run_class line up until the first '#', then truncate
        read(10,'(a80)',iostat=ios,err=2010)linebuffer080
        if(ios.ne.0)goto 2010
        iendstr = SCAN(linebuffer080, "#")
        if (iendstr.eq.0)then
             ! '#' not found, just copy linebuffer080 to run_class
          cdf_run_class = trim(adjustl(linebuffer080))
        else
            ! clip run class at key
          cdf_run_class = trim(linebuffer080(1:iendstr-1))
        endif

          ! Read institution line up until the first '#', then truncate
        read(10,'(a80)',iostat=ios,err=2010)linebuffer080
        if(ios.ne.0)goto 2010
        iendstr = SCAN(linebuffer080, "#")
        if (iendstr.eq.0)then
             ! '#' not found, just copy linebuffer080 to institution
          cdf_institution = trim(adjustl(linebuffer080))
        else
            ! clip institution at key
          cdf_institution = trim(linebuffer080(1:iendstr-1))
        endif


        ! This is the end of the standard blocks
      endif
      ! END OF BLOCK 9
      !************************************************************************

      !************************************************************************
      ! Searching for optional blocks labled by OPTMOD
      write(global_info,*)' *****************************************'
      write(global_info,*)' Reading Post-Block 9: optional modules   '
      write(global_info,*)' *****************************************'

      write(global_info,*)"Searching for blocks with OPTMOD"
      nmods = 0
      read(10,'(a80)',iostat=ios)linebuffer080
      do while(ios.eq.0)
        substr_pos1 = index(linebuffer080,'OPTMOD')
        if(substr_pos1.eq.1)then
          ! found an optional module
          nmods = nmods + 1
          !  Parse for the keyword
          read(linebuffer080,1104)mod_name
          OPTMOD_names(nmods) = trim(adjustl(mod_name))
          write(global_info,*)"     Found optional module : ",&
           OPTMOD_names(nmods),nmods
        endif
        read(10,'(a80)',iostat=ios)linebuffer080
1104    format(7x,a20)
      enddo

!     CLOSE Ash3d_ESP.inp
      write(global_info,25) infile
      write(global_log ,25) infile
      close(10)

      ! Here is the end of the file
2010  continue
      write(global_info,*)' *****************************************'
      write(global_info,*)' Finished reading/parsing input file      '
      write(global_info,*)' Custom blocks for optional modules will  '
      write(global_info,*)' be read in corresponding input_data_OPTMOD'
      write(global_info,*)' subroutines provided by the modules.     '
      write(global_info,*)' *****************************************'

      ! Write out computational scheme used
      if(useDS)then
        write(global_info,*)"Dimension splitting will be used."
        write(global_log ,*)"Dimension splitting will be used."
      else
        ! If something else besides dimension splitting is used (CTU, Semi-Lagr.)
        !  write out a note about it here.
        !write(global_info,*)""
        !write(global_log ,*)""
      endif
      write(global_info,*)limiter," limiter is used."
      write(global_log ,*)limiter," limiter is used."
      if (useCN) then
        write(global_info,*)&
         "Diffusion is calculated via Crank-Nicolson."
        write(global_log ,*)&
         "Diffusion is calculated via Crank-Nicolson."
      else
        write(global_info,*)"Diffusion is calculated explicitly."
        write(global_log ,*)"Diffusion is calculated explicitly."
      endif
 
      ! assign initial values
      !total_time = Simtime_in_hours ! total simulated time in seconds

!     CALCULATE SIZE OF GRID (cell-centered)
      if(IsLatLon) then
        nxmax = ceiling((lonUR-lonLL)/de)         !number of x nodes
        nymax = ceiling((latUR-latLL)/dn)         !number of y nodes
      else      
        nxmax = ceiling((xUR-xLL)/dx)         !number of x nodes
        nymax = ceiling((yUR-yLL)/dy)         !number of y nodes
      endif
      ! initialize active grid domain to the full grid
      imin = 1
      imax = nxmax
      jmin = 1
      jmax = nymax
      kmin = 1
      kmax = nzmax

      write(global_info,16) nxmax, nymax, nzmax
      write(global_log ,16) nxmax, nymax, nzmax
      write(global_info,*)
      write(global_log ,*)

      deallocate(iyear)
      deallocate(imonth)
      deallocate(iday)
      deallocate(hour)

      ! Set up logging logical values
        ! First check for output requests that require evaluating every time step
      if(Write_PT_Data.or.Write_PR_Data)then
        Output_every_TS = .true.
      else
        Output_every_TS = .false.
      endif
        ! Next check if there are requests at specific write times
      if (WriteDepositTS_ASCII           .or.      &
          WriteCloudConcentration_ASCII  .or.      &
          WriteCloudHeight_ASCII         .or.      &
          WriteCloudLoad_ASCII           .or.      &
          WriteDepositTS_KML             .or.      &
          WriteCloudConcentration_KML    .or.      &
          WriteCloudHeight_KML           .or.      &
          WriteCloudLoad_KML             .or.      &
          WriteReflectivity_KML          .or.      &
          Write3dFiles) then
        Output_at_WriteTimes = .true.
      else
        Output_at_WriteTimes = .false.
      endif
        ! Finally, check if we will be logging progress at log_steps
      if (log_step.gt.0)then
        Output_at_logsteps = .true.
      else
        Output_at_logsteps = .false.
      endif

      return
       
!******************************************************************************
!     ERROR TRAPS

      !ERROR TRAPS TO STDIN
1900  write(global_info,*)  'error: cannot find input file: ',infile
      write(global_info,*)  'Program stopped'
      write(global_log ,*)  'error: cannot find input file: ',infile
      write(global_log ,*)  'Program stopped'
      stop 1

      !BLOCK 1: GRID INFO
1901  write(global_info,*)  'error reading xLL or yLL.'
      write(global_info,*)  'You entered: ',cdf_b1l3
      write(global_info,*)  'Program stopped.'
      write(global_log ,*)  'error reading xLL or yLL'
      write(global_log ,*)  'You entered: ',cdf_b1l3
      write(global_log ,*)  'Program stopped.'
      stop 1
1902  write(global_info,*)  'error reading width and height of model domain.'
      write(global_info,*)  'You entered: ', cdf_b1l4
      write(global_info,*)  'Program stopped'
      write(global_log ,*)  'error reading width and height of model domain.'
      write(global_log ,*)  'You entered: ', cdf_b1l4
      write(global_log ,*)  'Program stopped'
      stop 1
!1903  write(global_info,*)  'error reading x or y of volcano.'
!      write(global_info,*)  'You entered: ',cdf_b1l5
!      write(global_info,*)  'Program stopped'
!      write(global_log ,*)  'error reading x or y of volcano.'
!      write(global_log ,*)  'You entered: ',cdf_b1l5
!      write(global_log ,*)  'Program stopped'
!      stop 1
1904  write(global_info,*)  'error reading dx or dy.'
      write(global_info,*)  'You entered: ', cdf_b1l6
      write(global_info,*)  'Program stopped'
      write(global_log ,*)  'error reading dx or dy.'
      write(global_log ,*)  'You entered: ', cdf_b1l6
      write(global_log ,*)  'Program stopped'
      stop 1
1905  write(global_info,*)  'error reading dz.'
      write(global_info,*)  'You gave: ',cdf_b1l7
      write(global_info,*)  'Program stopped.'      
      write(global_log ,*)  'error reading dz.'
      write(global_log ,*)  'You gave: ',cdf_b1l7
      write(global_log ,*)  'Program stopped.'      
      stop 1
!1906  write(global_info,*)  'error reading diffusion coefficient or Suzuki constant or plume type.'
!      write(global_info,*)  'The first value should be a number.  The second value should be either'
!      write(global_info,*)  'a number (the Suzuki constant), or the word "line", "point",'
!      write(global_info,*)  '"profile" or "umbrella" or "umbrella_air"'
!      write(global_info,*)  'You entered: ',cdf_b1l8
!      write(global_info,*)  'Program stopped'
!      write(global_log ,*)  'error reading diffusion coefficient or Suzuki constant or plume type.'
!      write(global_log ,*)  'The first value should be a number.  The second value should be either'
!      write(global_log ,*)  'a number (the Suzuki constant), or the word "line", "point"'
!      write(global_info,*)  '"profile" or "umbrella" or "umbrella_air"'
!      write(global_log ,*)  'You entered: ',cdf_b1l8
!      write(global_log ,*)  'Program stopped'
!      stop 1
1907  write(global_info,*)  'error reading number of eruptions.'
      write(global_info,*)  'You gave: ',cdf_b1l9
      write(global_info,*)  'Program stopped.'
      write(global_log ,*)  'error reading number of eruptions.'
      write(global_log ,*)  'You gave: ',cdf_b1l9
      write(global_log ,*)  'Program stopped.'
      stop 1      

      !BLOCK 2: ERUPTION PARAMETERS
1910  write(global_info,*)  'error reading start time, duration, height or',&
                  ' volume of an eruptive pulse.  Program stopped'
      write(global_log ,*)  'error reading start time, duration, height or',&
                  ' volume of an eruptive pulse.  Program stopped'
      stop 1         
1912  write(global_info,5678)  i,e_StartTime(i),e_StartTime(i-1)+e_Duration(i-1),hour(i)
      write(global_log ,5678)  i,e_StartTime(i),e_StartTime(i-1)+e_Duration(i-1),hour(i)
5678  format(4x,'error: eruption pulses are not in chronological order.',/, &
             4x,'e_StartTime(i)<(e_StartTime(i-1)+e_Duration(i-1))',/, &
             4x,'                               i=',i3,/, &
             4x,'                  e_StartTime(i)=',e15.8,/, &
             4x,'e_StartTime(i-1)+e_Duration(i-1)=',e15.8,/, &
             4x,'                        hour(i) =',f12.4,/, &
             4x,'Program stopped')
      stop 1

      !BLOCK 3: WIND PARAMETERS       
!1920  write(global_info,*)  'error reading iwind.  Program stopped'
!      write(global_log ,*)  'error reading iwind.  Program stopped'
!      stop 1        
1921  write(global_info,*)  'error reading simulation time in hours.',&
                  '  Program stopped'        
      write(global_log ,*)  'error reading simulation time in hours.',&
                  '  Program stopped'
      stop 1        
1922  write(global_info,*)  'Error reading whether to stop simulation when'
      write(global_info,*)  '99% of erupted volume has deposited.'
      write(global_info,*)  'Answer should be yes or no.'
      write(global_info,*)  'You gave: ',linebuffer080
      write(global_info,*)  'Program stopped'
      write(global_log ,*)  'Error reading whether to stop simulation when'
      write(global_log ,*)  '99% of erupted volume has deposited.'
      write(global_log ,*)  'Answer should be yes or no.'
      write(global_log ,*)  'You gave: ',linebuffer080
      write(global_log ,*)  'Program stopped'
      stop 1
1923  write(global_info,*)  'error reading number of wind files.'
      write(global_info,*)  'Answer should be a positive integer.'
      write(global_info,*)  'You gave: ',linebuffer080
      write(global_info,*)  ' Program stopped'
      write(global_log ,*)  'error reading number of wind files.'
      write(global_log ,*)  'Answer should be a positive integer.'
      write(global_log ,*)  'You gave: ',linebuffer080
      write(global_log ,*)  ' Program stopped'
      stop 1 
!1930  write(global_info,*)'iwind must be between 1 and 4. Program stopped'
!      write(global_info,*)'  IWIND OPTIONS:'
!      write(global_info,*)'  iwind = 1 read from a 1-D wind sounding'
!      write(global_info,*)'          2 read from 3D gridded ASCII files'
!      write(global_info,*)'          3 read directly from a single NetCDF file'
!      write(global_info,*)'          4 read from multiple NetCDF files'
!      write(global_log ,*)'iwind must be between 1 and 4. Program stopped'
!      write(global_log ,*)'  IWIND OPTIONS:'
!      write(global_log ,*)'  iwind = 1 read from a 1-D wind sounding'
!      write(global_log ,*)'          2 read from 3D gridded ASCII files'
!      write(global_log ,*)'          3 read directly from a single NetCDF file'
!      write(global_log ,*)'          4 read from multiple NetCDF files'
!      stop 1
!1931  write(global_info,*)'iwindformat must be between 1 and 10 or 20-24.' 
!      write(global_info,*)'iwindformat=',iwindformat,'. Program stopped'
!      write(global_log ,*)'iwindformat must be between 1 and 10 or 20-24.'
!      write(global_log ,*)'iwindformat=',iwindformat,'. Program stopped'
!      stop 1
1932  write(global_info,*) 'Error reading iHeightHandler. iHeightHandler must be 1 or 2. You entered:'
      write(global_info,*) linebuffer080
      write(global_info,*) 'Program stopped'
      write(global_log ,*) 'Error reading iHeightHandler. iHeightHandler must be 1 or 2. You entered:'
      write(global_log ,*) linebuffer080
      write(global_log ,*) 'Program stopped'
      stop 1

      !BLOCK 4: OUTPUT FILE OPTIONS
1953  write(global_info,*)'Error reading whether to print out ESRI ASCII file of',&
                ' deposit thickness.'
      write(global_info,*)'The first characters on this line should be ''yes'' or',&
                ' ''no''.  Program stopped'      
      write(global_log ,*)'Error reading whether to print out ESRI ASCII file of',&
                ' deposit thickness.'
      write(global_log ,*)'The first characters on this line should be ''yes'' or',&
                ' ''no''.  Program stopped'      
      stop 1
1954  write(global_info,*)'Error reading whether to print out KML file of',&
                ' deposit thickness.'
      write(global_info,*)'The first characters on this line should',&
                ' be yes or no.  Program stopped'      
      write(global_log ,*)'Error reading whether to print out KML file of',&
                ' deposit thickness.'
      write(global_log ,*)'The first characters on this line should',&
                ' be yes or no.  Program stopped'      

      stop 1
1955  write(global_info,*)'Error reading whether to write out KML deposit files at',&
                ' specifiied times.'
      write(global_info,*)'The first characters on this line should be ''yes'' or',&
                ' ''no''.  Program stopped'      
      write(global_log ,*)'Error reading whether to write out KML deposit files at',&
                ' specifiied times.'
      write(global_log ,*)'The first characters on this line should be ''yes'' or',&
                ' ''no''.  Program stopped'      
      stop 1
1956  write(global_info,*)'Error reading whether to print out KML deposit',&
                ' files at specifiied times.'
      write(global_info,*)'The first characters on this line should',&
                ' be yes or no.  Program stopped'      
      write(global_log ,*)'Error reading whether to print out KML deposit',&
                ' files at specifiied times.'
      write(global_log ,*)'The first characters on this line should',&
                ' be yes or no.  Program stopped'      
      stop 1
1957  write(global_info,*)'Error reading whether to print out ASCII files of',&
                ' cloud concentration at specifiied times.'
      write(global_info,*)'The first characters on this line should be ''yes'' or',&
                ' ''no''.  Program stopped'      
      write(global_log ,*)'Error reading whether to print out ASCII files of',&
                ' cloud concentration at specifiied times.'
      write(global_log ,*)'The first characters on this line should be ''yes'' or',&
                ' ''no''.  Program stopped'      
      stop 1
1958  write(global_info,*)'Error reading whether to print out a KML file of',&
                ' cloud concentration.'
      write(global_info,*)'The first characters on this line should',&
                ' be yes or no.  Program stopped'      
      write(global_log ,*)'Error reading whether to print out a KML file of',&
                ' cloud concentration.'
      write(global_log ,*)'The first characters on this line should',&
                ' be yes or no.  Program stopped'      
      stop 1
1959  write(global_info,*)'Error reading whether to print out ASCII files of',&
                ' cloud height.'
      write(global_info,*)'The first characters on this line should',&
                ' be yes or no.  Program stopped'      
      write(global_log ,*)'Error reading whether to print out ASCII files of',&
                ' cloud height.'
      write(global_log ,*)'The first characters on this line should',&
                ' be yes or no.  Program stopped'      
      stop 1
1960  write(global_info,*)'Error reading whether to print out a KML file',&
                ' of cloud height.'
      write(global_info,*)'The first characters on this line should',&
                ' be yes or no.  Program stopped'      
      write(global_log ,*)'Error reading whether to print out a KML file',&
                ' of cloud height.'
      write(global_log ,*)'The first characters on this line should',&
                ' be yes or no.  Program stopped'      
      stop 1      
19601 write(global_info,*)'Error reading whether to print out an ASCII file',&
                ' of cloud load.'
      write(global_info,*)'The first characters on this line should',&
                ' be yes or no.  Program stopped'      
      write(global_log ,*)'Error reading whether to print out an ASCII file',&
                ' of cloud load.'
      write(global_log ,*)'The first characters on this line should',&
                ' be yes or no.  Program stopped'      
      stop 1      
1961 write(global_info,*)'Error reading whether to print out a KML file ',&
                ' of cloud load.'
      write(global_info,*)'The first characters on this line should',&
                ' be yes or no.  You gave:'      
      write(global_info,*) linebuffer080
      write(global_info,*) 'program stopped'
      write(global_log ,*)'Error reading whether to print out a KML file ',&
                ' of cloud load.'
      write(global_log ,*)'The first characters on this line should',&
                ' be yes or no.  You gave:'      
      write(global_log ,*) linebuffer080
      write(global_log ,*) 'program stopped'
      stop 1      
1962 write(global_info,*)'Error reading whether to print out an ASCII file ',&
                ' of deposit arrival time.'
      write(global_info,*)'The first characters on this line should',&
                ' be yes or no.  You gave:'      
      write(global_info,*) linebuffer080
      write(global_info,*) 'program stopped'
      write(global_log ,*)'Error reading whether to print out an ASCII file ',&
                ' of deposit arrival time.'
      write(global_log ,*)'The first characters on this line should',&
                ' be yes or no.  You gave:'      
      write(global_log ,*) linebuffer080
      write(global_log ,*) 'program stopped'
      stop 1
1963 write(global_info,*)'Error reading whether to print out a KML file ',&
                ' of deposit arrival time.'
      write(global_info,*)'The first characters on this line should',&
                ' be yes or no.  You gave:'      
      write(global_info,*) linebuffer080
      write(global_info,*) 'program stopped'
      write(global_log ,*)'Error reading whether to print out a KML file ',&
                ' of deposit arrival time.'
      write(global_log ,*)'The first characters on this line should',&
                ' be yes or no.  You gave:'      
      write(global_log ,*) linebuffer080
      write(global_log ,*) 'program stopped'
      stop 1      
19621 write(global_info,*)'Error reading whether to print out an ASCII file ',&
                ' of cloud arrival time.'
      write(global_info,*)'The first characters on this line should',&
                ' be yes or no.  You gave:'      
      write(global_info,*) linebuffer080
      write(global_info,*) 'program stopped'
      write(global_log ,*)'Error reading whether to print out an ASCII file ',&
                ' of deposit arrival time.'
      write(global_log ,*)'The first characters on this line should',&
                ' be yes or no.  You gave:'      
      write(global_log ,*) linebuffer080
      write(global_log ,*) 'program stopped'
      stop 1
19631 write(global_info,*)'Error reading whether to print out a KML file ',&
                ' of cloud arrival time.'
      write(global_info,*)'The first characters on this line should',&
                ' be yes or no.  You gave:'      
      write(global_info,*) linebuffer080
      write(global_info,*) 'program stopped'
      write(global_log ,*)'Error reading whether to print out a KML file ',&
                ' of deposit arrival time.'
      write(global_log ,*)'The first characters on this line should',&
                ' be yes or no.  You gave:'      
      write(global_log ,*) linebuffer080
      write(global_log ,*) 'program stopped'
      stop 1
1964  write(global_info,*)'Error reading whether to print out 3-D ash',&
                ' concentration files at specifiied times.'
      write(global_info,*)'The first characters on this line should be ''yes'' or',&
                ' ''no''.  Program stopped'      
      write(global_info,*)'You gave:'
      write(global_info,*) linebuffer080
      write(global_log ,*)'Error reading whether to print out 3-D ash',&
                ' concentration files at specifiied times.'
      write(global_log ,*)'The first characters on this line should be ''yes'' or',&
                ' ''no''.  Program stopped'      
      write(global_log ,*)'You gave:'
      write(global_log ,*) linebuffer080
      stop 1
1965  write(global_info,*)'Error reading format of output files.'
      write(global_info,*)'The first characters on this line should',&
                ' be ''ascii'', ''binary'', or ''netcdf''.'      
      write(global_info,*)'Program stopped.'
      write(global_info,*)'You gave:'
      write(global_info,*) linebuffer080
      write(global_log ,*)'Error reading format of output files.'
      write(global_log ,*)'The first characters on this line should',&
                ' be ''ascii'', ''binary'', or ''netcdf''.'      
      write(global_log ,*)'Program stopped.'
      write(global_log ,*)'You gave:'
      write(global_log ,*) linebuffer080
      stop 1
1966  write(global_info,*)'Error reading the number of files to be written out.',&
                '  This should be a positive integer, or -1.'
      write(global_info,*)'You gave: ',nWriteTimes
      write(global_info,*) 'Program stopped.'      
      write(global_log ,*)'Error reading the number of files to be written out.',&
                '  This should be a positive integer, or -1.'
      write(global_log ,*)'You gave: ',nWriteTimes
      write(global_log ,*)'Program stopped.'      
      stop 1
1967  write(global_info,*)'Error reading the times at which output files are to',&
                ' be written out.'
      write(global_info,*)'This should be one or more real numbers.'
      write(global_info,*)'You gave: ',WriteInterval
      write(global_info,*)'  Program stopped.'
      write(global_log ,*)'Error reading the times at which output files are to',&
                ' be written out.'
      write(global_log ,*)'This should be one or more real numbers.'
      write(global_log ,*)'You gave: ',WriteInterval
      write(global_log ,*)'  Program stopped.'
      stop 1
1968  write(global_info,*)'Error: some write times are <0.  Program stopped.'
      write(global_info,*)'Error: some write times are <0.  Program stopped.'
      stop 1
1969  write(global_info,*)'Error: some write times are not in chronological',&
                ' order.  Program stopped.'      
      write(global_info,*)'Error: some write times are not in chronological',&
                ' order.  Program stopped.'
      stop 1

      !BLOCK 5: WIND FILES
1970  write(global_info,*)  'error reading wind file names.  Program stopped'
      write(global_log ,*)  'error reading wind file names.  Program stopped'
      stop 1        
!1971  write(global_info,*)  'error: cannot find file input wind file.',&
!                  '  Program stopped.'
!      write(global_log ,*)  'error: cannot find file input wind file.',&
!                  '  Program stopped.'
!      stop 1

      !BLOCK 6: AIRPORT FILE OUTPUT OPTIONS
1980  write(global_info,*) 'Error reading whether to write out ASCII airport file.'
      write(global_info,*) 'Answer must be yes or no.'
      write(global_info,*) 'You gave:',cdf_b6l1
      write(global_info,*) 'Program stopped'
      write(global_log ,*) 'Error reading whether to write out ASCII airport file.'
      write(global_log ,*) 'Answer must be yes or no.'
      write(global_log ,*) 'You gave:',cdf_b6l1
      write(global_log ,*) 'Program stopped'
      stop 1
1981  write(global_info,*) 'Error reading whether to write out grain-size distribution to ASCII airport file.'
      write(global_info,*) 'Answer must be yes or no.'
      write(global_info,*) 'You gave:',cdf_b6l2
      write(global_info,*) 'Program stopped'
      write(global_log ,*) 'Error reading whether to write out grain-size distribution to ASCII airport file.'
      write(global_log ,*) 'Answer must be yes or no.'
      write(global_log ,*) 'You gave:',cdf_b6l2
      write(global_log ,*) 'Program stopped'
      stop 1
1982  write(global_info,*) 'Error reading whether to write out KML airport file.'
      write(global_info,*) 'Answer must be yes or no.'
      write(global_info,*) 'You gave:',cdf_b6l3
      write(global_info,*) 'Program stopped'
      write(global_log ,*) 'Error reading whether to write out KML airport file.'
      write(global_log ,*) 'Answer must be yes or no.'
      write(global_log ,*) 'You gave:',cdf_b6l3
      write(global_log ,*) 'Program stopped'
      stop 1
      !If we can't open the external file
19825 write(global_info,*) 'You gave the name of the following point location file to read:'
      write(global_info,*) AirportInFile
      write(global_info,*) 'But Ash3d could not open that file.'
19827 write(global_info,*) 'Would you like Ash3d to use the interal airports database instead (y/n)?'
      read(5,'(a1)') answer
      if (answer.eq.'y') then
         ReadExtAirportFile=.false.
         AppendExtAirportFile=.false.
         AirportInFile='internal'
       else if (answer.eq.'n') then
         write(global_info,*) 'program stopped'
         stop 1
       else
         goto 19827
      endif
      goto 19828
1983  write(global_info,*) 'Error reading whether to project airport coordinates using proj4.'
      write(global_info,*) 'Answer must be yes or no.'
      write(global_info,*) 'You gave:',cdf_b6l5
      write(global_info,*) 'Program stopped'
      write(global_log ,*) 'Error reading whether to project airport coordinates using proj4.'
      write(global_log ,*) 'Answer must be yes or no.'
      write(global_log ,*) 'You gave:',cdf_b6l5
      write(global_log ,*) 'Program stopped'
      stop 1

     !BLOCK 7: GRAIN-SIZE BINS, SETTLING VELOCITY
!1990  write(global_info,*) 'Error reading number of grain-size bins to use.'
!      write(global_info,*) 'Answer should be a positive integer.'
!      write(global_info,*) 'You gave: ',linebuffer080
!      write(global_info,*) 'Program stopped'
!      write(global_log ,*) 'Error reading number of grain-size bins to use.'
!      write(global_log ,*) 'Answer should be a positive integer.'
!      write(global_log ,*) 'You gave: ',linebuffer080
!      write(global_log ,*) 'Program stopped'
!      stop 1

     !BLOCK 8: VERTICAL PROFILES
2000  write(global_info,*) 'Error reading the number of vertical profiles to use.'
      write(global_info,*) 'Answer should be an integer.'
      write(global_info,*) 'You gave: ',linebuffer080
      write(global_info,*) 'Program stopped.'
      write(global_log ,*) 'Error reading the number of vertical profiles to use.'
      write(global_log ,*) 'Answer should be an integer.'
      write(global_log ,*) 'You gave: ',linebuffer080
      write(global_log ,*) 'Program stopped.'
      stop 1
!2001  write(global_info,*) 'Error in x or y location of a vertical profile.'
!      write(global_info,*) 'Answer should be two real numbers.'
!      write(global_info,*) 'You gave: ',linebuffer080
!      write(global_info,*) 'Program stopped.'
!      write(global_log ,*) 'Error in x or y location of a vertical profile.'
!      write(global_log ,*) 'Answer should be two real numbers.'
!      write(global_log ,*) 'You gave: ',linebuffer080
!      write(global_log ,*) 'Program stopped.'
!      stop 1

      !BLOCK 9: NETCDF ANNOTATIONS

!***********************************************************************
!     format STATEMENT

!2     format(4x,'Ash3d (Rev ',a5,') run ',&
!             i4,'.',i2.2,'.',i2.2,i4,':',i2.2,' UTC')
!102   format(i4,'.',i2.2,'.',i2.2,i4,':',i2.2)
3     format(4x,'opening input file ',a130)
37    format(4x,'Volcano name:',a30,/)
4     format(4x,'MAP & GRID SETUP:',/,&
      4x,' coordinates of lower-left point:  x=',T48,f10.3,&
                 ', y=',T62,f10.3,/,&
      4x,'               extent of grid in:  x=',T48,f10.3,&
                 ', y=',T62,f10.3,/ &
      4x,'          coordinates of volcano:  x=',T48,f10.3,&
                 ', y=',T62,f10.3,/)
5     format(4x,'nodal spacing in x: ',f8.3,/,&
             4x,'              in y: ',f8.3)
43    format(4x,'              in z: ',f8.3,/)
6     format(4x,'ERUPTION SOURCE PARAMETERS:',/,&
             4x,'                diffusion coefficient (km2/hr): ',f8.4,/,&
             4x,'                               Suzuki constant: ',f8.4,/,&
             4x,'Stop calculation when 99% of ash has deposited: ',L1,/,&
             4x,'                         Simulation time (hrs): ',f8.4,/)
1438   format(4x,'ERUPTION SOURCE PARAMETERS:',/,&
             4x,'                diffusion coefficient (km2/hr): ',f8.4,/,&
             4x,'Stop calculation when 99% of ash has deposited: ',L1,/,&
             4x,'                         Simulation time (hrs): ',f8.4,/)
1439  format(4x,'                                     SourceType: ',a9)
7     format(4x,'Number of eruptions or pulses: ',i8,//, &
             4x,&
             '          plume        start time      duration   volume',/, &
             4x,&
             'number  height (km)  (yyyymmdd hh.hh)    (hrs)    (km3 DRE)')
8     format(4x,i6,f13.2,3x,i4,i2.2,i2.2,f6.2,2f10.4)
9     format(/,4x,'Number of grain-size bins:     ',i2, &
             /,4x,'               Fall Model:     ',i2)
10    format(4x,'Bins:',/,4x,&
        'mass fraction      diameter (mm)     density (kg/m3)      F     G         phi')
11    format(8x,f10.5,4x,f11.6,10x,f10.4,5x,f8.2,2x,f4.2,5x,f5.2)
2110  format(4x,'Bins:',/,4x,&
             'mass fraction      v_s (m/s)')
2111  format(8x,f5.3,4x,f11.6)
!12    format(4x,/,'WIND FILE SETUP:',/, &
!                 4x,'using 1-D wind sounding.  wind file=',a130)
13    format(4x,/,'WIND FILE SETUP:',/, &
                 4x,'Reading 4-D gridded wind data from files:')
!14    format(12x,a80)
!15    format(/,4x,'Opening Ash3d_1Dwind.inp')
16    format(4x,'Number of nodes in x:',T48,i3,/, &
             4x,'                in y:',T48,i3,/, &
             4x,'                in z:',T48,i3)
25    format(/,4x,'closing ',a80)
!26    format(4x,'oops.  iwind=2 but the current program can only', &
!                     ' use iwind=1 (1-D wind sounding).',          &
!                     '  Program stopped')
!27    format(/,4x,'1-D WIND INPUTS:',/, &
!                  4x,'time of Profile (UT):',T37,i4,/,&
!                  4x,'Number of wind levels:',T37,i4)
!28    format(4x,'Location of wind profile',T34,'x=',f7.0,T45,&
!                    'y=',f7.0,/,&
!                4x,'velocity    direction         ',&
!                   'u (E)           v (N)      Elevation',/,&
!                4x,'  m/s       deg. E of N        m/s',&
!                   '             m/s          m')
!29    format(4x,f8.2,f8.2,3(5x,f10.2))
!129   format(4x,f8.2,f8.2,5(5x,f10.2))
!30    format(4x,'MAP & GRID SETUP:',/,&
!      4x,'Lower-left point:   longitude (deg. E)=',T48,f10.3,&
!                 ', latitude (deg. N)',T62,f10.3,/,&
!      4x,'Upper-right point:  longitude (deg. E)=',T48,f10.3,&
!                 ', latitude (deg. N)',T62,f10.3)
!31     format(4x,'nodal spacing (degrees) in x: ',f8.3,/,&
!             4x,'                         in y: ',f8.3,/,&
!             4x,'                         in z: ',f8.3)
32     format(/,4x,'Warning:  some write times exceed the model',&
                   ' simulation time.',/, &
                4x,'These times will be ignored.')
33     format   (4x,'Chosen output:',/, &
                4x,'                        Output final ASCII deposit',&
                   ' file (T/F) = ',L1,/, &
                4x,'                          Output final KML deposit',&
                   ' file (T/F) = ',L1,/, &
                4x,'           Output ASCII deposit file at specified',&
                   ' times (T/F) = ',L1,/, &
                4x,'             Output KML deposit file at specified',&
                   ' times (T/F) = ',L1,/, &
                4x,'Output ASCII ash-cloud concentration at specified',&
                   ' times (T/F) = ',L1,/, &
                4x,'  Output KML ash-cloud concentration at specified',&
                   ' times (T/F) = ',L1,/, &
                4x,'       Output ASCII ash-cloud height at specified',&
                   ' times (T/F) = ',L1,/, &
                4x,'         Output KML ash-cloud height at specified',&
                   ' times (T/F) = ',L1,/, &
                4x,'             Output ASCII cloud load at specified',&
                   ' times (T/F) = ',L1,/, &
                4x,'               Output KML cloud load at specified',&
                   ' times (T/F) = ',L1,/, &
                4x,'             Output ASCII file of deposit arrival',&
                   '  time (T/F) = ',L1,/, &
                4x,'               Output KML file of deposit arrival',&
                   '  time (T/F) = ',L1,/, &
                4x,'               Output ASCII file of cloud arrival',&
                   '  time (T/F) = ',L1,/, &
                4x,'                 Output KML file of cloud arrival',&
                   '  time (T/F) = ',L1,/, &
                4x,'                     Output 3d files at specified',&
                   ' times (T/F) = ',L1,/, &
                4x,'                  Output format of 3d ash',&
                   ' concentration files = ',a6,/, &
                4x,'                                         Number',&
                   ' of time steps = ',i3,/)
34    format(/, 4x,'Files to be written out at the following hours',&
                   ' after the eruption start:')
35    format(8x,f10.3)
!38    format(/4x,'Warning: KML files cannot be written out unless the input values are given',/, &
!              4x,'in latitude and longitude.  These files will not be written out.')
!39    format( 4x,'Do you wish to continue (y/n)?')
44    format(4x,'Airport choices:',/, &
                4x,'   Read external file of locations (T/F) = ',L1,/, &
                4x,'   Append external locations to airport list (T/F) = ',L1,/, &
                4x,'   Write ASCII file of ash arrival times at',&
                   ' airports (T/F) = ',L1,/, &
                4x,'   Write deposit grain-size distribution to',&
                   ' airports file (T/F) = ',L1,/, &
                4x,'     Write KML file of ash arrival times at',&
                   ' airports (T/F) = ',L1,/, &
                4x,'     Calculate projected airport locations using',&
                   ' Ash3d (T/F) = ',L1,//, &
                4x,' Name of file containing airport locations: ',a130)
!45    format (/,4x,'Error reading whether to stop calculation when 99%',&
!                  ' of deposit  has accumulated.',/, &
!               4x,'Answer must be yes or no.  You gave:',&
!                  /,a80,//,'  Program stopped')


      end subroutine Read_Control_File
      
!******************************************************************************

      subroutine LatLonChecker(latLL,lonLL,lat_volcano,lon_volcano,gridwidth_e,gridwidth_n)
      
      !subroutine that checks the domain for errors if IsLatLon=.true.
      
      use precis_param

      use io_units

      implicit none

      real(kind=ip), intent(in)    :: latLL
      real(kind=ip), intent(in)    :: lonLL
      real(kind=ip), intent(in)    :: lat_volcano
      real(kind=ip), intent(inout) :: lon_volcano ! we might need to remap this value
      real(kind=ip), intent(in)    :: gridwidth_e
      real(kind=ip), intent(in)    :: gridwidth_n

      !MAKE SURE THAT LATITUDE IS BETWEEN -90 AND 90.
      if(abs(latLL).gt.90.0_ip.or.abs(lat_volcano).gt.90.0_ip)then
        write(global_info,*)"ERROR: latLL and lat_volcano should be",&
                  " in in range -90 - 90"
        write(global_info,*)"latLL       = ",latLL
        write(global_info,*)"lat_volcano = ",lat_volcano
        stop 1
      endif
      
      !MAKE SURE THAT LONGITUDE OF LEFT SIDE OF GRID IS BETWEEN -360 AND 360.
      if(abs(lonLL).gt.360.0_ip.or.abs(lon_volcano).gt.360.0_ip)then
        write(global_info,*) "ERROR: lonLL and lon_volcano should be",&
                   " in in range -360 - 360"
        write(global_info,*)"lonLL       = ",lonLL
        write(global_info,*)"lon_volcano = ",lon_volcano
        stop 1
      endif
      
      !MAKE SURE THAT gridwidth_e AND gridwidth_n ARE POSITIVE
      if((gridwidth_e.lt.0.0_ip).or.(gridwidth_n.lt.0.0_ip))then
        write(global_info,*)"ERROR: gridwidth_e and gridwidth_n must be positive."
        stop 1
      endif
      
      !MAKE SURE THAT THE TOP OF THE GRID doES NOT EXTEND BEYOND 90n
      if ((latLL+gridwidth_n).gt.90.0_ip) then
        write(global_info,*)"ERROR:  Latitude at the top of the grid > 90."
        stop 1        
      endif
 
      !MAKE SURE THAT lon_volcano>lonLL
      if (lonLL.gt.lon_volcano) lon_volcano=lon_volcano+360.0_ip

      !MAKE SURE THAT THE VOLCANO IS WITHIN THE MODEL REGION (LONGITUDE)
      if ((lon_volcano.lt.lonLL).or.(lon_volcano.gt.(lonLL+gridwidth_e))) then
        write(global_info,*) 'ERROR: the volcano is not within the specified longitude region'
        write(global_info,*) "lon_volcano=",lon_volcano,', lonLL=',lonLL
        write(global_info,*) 'grid_width=',gridwidth_e
        stop 1
      endif

      !MAKE SURE THE VOLCANO IS WITHIN THE MODEL REGION (LATITUDE)
      if ((lat_volcano.lt.latLL).or.(lat_volcano.gt.(latLL+gridwidth_n))) then
        write(global_info,*) 'ERROR: the volcano is not within the specified latitude region'
        write(global_info,*) "lat_volcano=",lat_volcano,', latLL=',latLL
        write(global_info,*) 'grid_height=',gridwidth_n
        stop 1
      endif

      return
      
      end subroutine LatLonChecker
      
!******************************************************************************

      subroutine xyChecker(xLL,yLL,dx,dy,x_volcano,y_volcano,gridwidth_x,gridwidth_y)

      !subroutine that checks the domain for errors if IsLatLon=.false.
      
      use precis_param

      use io_units

      use global_param,    only : &
        EPS_SMALL

      implicit none

      real(kind=ip), intent(in) :: xLL,yLL
      real(kind=ip), intent(in) :: dx,dy
      real(kind=ip), intent(in) :: x_volcano,y_volcano
      real(kind=ip), intent(in) :: gridwidth_x,gridwidth_y
      
      !MAKE SURE THAT gridwidth_x AND gridwidth_y ARE POSITIVE
      if((gridwidth_x.lt.0.0_ip).or.(gridwidth_y.lt.0.0_ip))then
        write(global_info,*) &
          "ERROR: gridwidth_x and gridwidth_y must be positive."
        stop 1
      endif
      
      !MAKE SURE THAT DX AND DY ARE POSITIVE
      if((dx.lt.0.0_ip).or.(dy.lt.0.0_ip)) then
        write(global_info,*) "ERROR: dx and dy must be positive."
        stop 1
      endif
      
      !MAKE SURE THE VOLCANO IS WITHIN THE MODEL REGION
      if (((x_volcano.lt.xLL).or.(x_volcano.gt.(xLL+gridwidth_x))).or. &
          ((y_volcano.lt.yLL).or.(y_volcano.gt.(yLL+gridwidth_y)))) then
        write(global_info,*) &
          'ERROR: the volcano is not within the model region'
        stop 1
      endif
      
      !PRINT OUT WARNING MESSAGE if DX != DY
      if (abs(dx-dy).lt.EPS_SMALL) then
          write(global_info,1)              !print out a warning message about the deposit file
          write(global_log ,1)              !if dx and dy aren't the same
      endif

      !format STATEMENTS
1     format (4x,'Warning: dx and dy are not the same.  If the',&
                ' deposit file is read by ArcMap',/, &
             4x,'they  are assumed to be the same.  The nodal',&
                ' spacing written to the ',/, &
             4x,'deposit file is dx, dy, but ArcMap will only read dx.')
               
      return               
      
      end subroutine xyChecker
      

!******************************************************************************

      subroutine vprofchecker(iprof)

!     SUBROUTINE THAT CHECKS THE LOCATIONS OF VERTICAL PROFILES SPECIFIED

      use precis_param

      use io_units

      use mesh,              only : &
        IsLatLon,de,dn,lonLL,latLL,gridwidth_e,gridwidth_n,&
        dx,dy,xLL,yLL,gridwidth_x,gridwidth_y

      use io_data,           only : &
        x_vprofile,y_vprofile,i_vprofile,j_vprofile

      implicit none

      integer, intent(in) :: iprof
      real(kind=ip) :: lon_vprof, lat_vprof

!     FIND THE I AND J VALUES OF THE NODE CONTAINING THE VERTICAL PROFILE

      if (IsLatLon) then
        lon_vprof = x_vprofile(iprof)
        lat_vprof = y_vprofile(iprof)
        if ((lonLL+gridwidth_e-lon_vprof).gt.360.0_ip) &
          lon_vprof = lon_vprof+360.0_ip
        i_vprofile(iprof) = int((lon_vprof-lonLL)/de) + 1
        j_vprofile(iprof) = int((lat_vprof-latLL)/dn) + 1
      else
        i_vprofile(iprof) = int((x_vprofile(iprof)-xLL)/dx) + 1
        j_vprofile(iprof) = int((y_vprofile(iprof)-yLL)/dy) + 1
      endif
      write(global_info,47) iprof, real(x_vprofile(iprof),kind=4),&
                                   real(y_vprofile(iprof),kind=4),&
                            i_vprofile(iprof), j_vprofile(iprof)
      write(global_log ,47) iprof, real(x_vprofile(iprof),kind=4),&
                                   real(y_vprofile(iprof),kind=4),&
                            i_vprofile(iprof), j_vprofile(iprof)
47    format(i11,2f10.3,2i6)

           !MAKE SURE THE POINT IS WITHIN THE MODEL REGION
      if (IsLatLon) then
        if (((lon_vprof.lt.lonLL)               .or. &
             (lon_vprof.gt.(lonLL+gridwidth_e))).or. &
            ((lat_vprof.lt.latLL)               .or. &
             (lat_vprof.gt.(latLL+gridwidth_n)))) then
          write(global_info,*) &
            'ERROR: this location is not within the model region'
          write(global_log ,*) &
            'ERROR: this location is not within the model region'
          stop 1
        endif
      else
        if (((x_vprofile(iprof).lt.xLL)               .or.&
             (x_vprofile(iprof).gt.(xLL+gridwidth_x))).or.&
            ((y_vprofile(iprof).lt.yLL)               .or.&
            (y_vprofile(iprof).gt.(yLL+gridwidth_y)))) then
            write(global_info,*)&
              'ERROR: this location is not within the model region'
            write(global_log ,*)&
              'ERROR: this location is not within the model region'
            stop 1
        endif
      endif

      return

      end subroutine vprofchecker

