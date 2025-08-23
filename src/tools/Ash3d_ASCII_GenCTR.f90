! Ash3d_ASCII_GenCTR

      program Ash3d_ASCII_GenCTR

      use precis_param

      use io_units

      use io_data,       only : &
         infile,concenfile,datafileIn,cdf_title,cdf_comment,VolcanoName, &
         Have_Block_NetCDF,Have_Block_ResParm,Have_Block_Topo,Have_Block_VarDiff,&
         cdf_b1l2,cdf_vardz,cdf_b3l5, &
         cdf_b4l1,cdf_b4l2,cdf_b4l3,cdf_b4l4,cdf_b4l5,cdf_b4l6,cdf_b4l7,cdf_b4l8,cdf_b4l9,cdf_b4l10,&
         cdf_b4l11,cdf_b4l12,cdf_b4l13,cdf_b4l14,cdf_b4l15,cdf_b4l17,cdf_b4l18,cdf_b5l1,  &
         cdf_b6l1,cdf_b6l2,cdf_b6l3,cdf_b6l4,cdf_b6l5, &
         nvprofiles,Site_vprofile,x_vprofile, y_vprofile

      use Ash3d_Program_Control, only : &
           Read_Control_File

      use mesh,          only : &
         lonLL,latLL,gridwidth_e,gridwidth_n,xLL,yLL,gridwidth_x,gridwidth_y, &
         IsLatLon,de,dn,dx,dy,VarDzType,dz_const

      use solution,      only : &
         StopWhenDeposited

      use time_data,     only : &
         SimStartHour,Simtime_in_hours,BaseYear,useLeap

      use Source,        only : &
         lon_volcano,lat_volcano,x_volcano,y_volcano,z_volcano,Suzuki_A,      &
         neruptions,SourceType_idx,e_StartTime,e_Duration,         &
         e_PlumeHeight,e_Volume,e_prof_dz,e_prof_nzpoints,e_prof_Volume

      use Tephra,        only : &
         n_gs_max,Tephra_gsdiam,Tephra_bin_mass,Tephra_v_s,Tephra_rho_m,FV_ID,&
         Tephra_gsF,Tephra_gsG,Tephra_gsPhi,Shape_ID,&
         LN_massfrac,LN_phi_mean,LN_phi_stddev

      use Diffusion,     only : &
         diffusivity_horz

      use Topography,    only : &
         var_User_charlines_Topo

      use Diffusivity_Variable, only : &
         var_User_charlines_VarDiff

      use Ash3d_Program_Control, only : &
         Write_input_block_header,      &
         SetWrite_input_block_01,       &
         SetWrite_input_block_02,       &
         SetWrite_input_block_03,       &
         SetWrite_input_block_04,       &
         SetWrite_input_block_05,       &
         SetWrite_input_block_06,       &
         SetWrite_input_block_07,       &
         SetWrite_input_block_08,       &
         SetWrite_input_block_09,       &
         SetWrite_input_block_ResetParam,&
         SetWrite_input_block_Topo,      &
         SetWrite_input_block_VarDiff

      use MetReader,     only : &
         MR_iWind,MR_iWindFormat,MR_iGridCode,MR_iDataFormat,MR_iHeightHandler,&
         MR_WindFiles,MR_VERB

      implicit none

      integer            :: nargs
      integer            :: stat
      character(len= 80) :: linebuffer080
      integer            :: RunID

      logical       :: IsThere
      real(kind=ip) :: LLx
      real(kind=ip) :: LLy
      real(kind=ip) :: widthx
      real(kind=ip) :: widthy
      real(kind=ip) :: x_in
      real(kind=ip) :: y_in
      real(kind=ip) :: z_in
      real(kind=ip) :: dx_in
      real(kind=ip) :: dy_in
      integer       :: dz_type

      ! Block 3 variabls
      integer          :: nwindfiles

      ! Block 4 variables
      integer          :: iostatus
      character(len=50):: iomessage
      character(len=3) :: answer
      character(len=1) :: WriteDepositFinal_ASCII_c        ! n/y Write out ESRI ASCII file of final deposit thickness?
      character(len=1) :: WriteDepositFinal_KML_c          ! n/y Write out        KML file of final deposit thickness?
      character(len=1) :: WriteDepositTS_ASCII_c           ! Write out ESRI ASCII deposit files at specified times?
      character(len=1) :: WriteDepositTS_KML_c             ! Write out        KML deposit files at specified times?
      character(len=1) :: WriteCloudConcentration_ASCII_c  ! Write out ESRI ASCII files of ash-cloud concentration?
      character(len=1) :: WriteCloudConcentration_KML_c    ! Write out        KML files of ash-cloud concentration?
      character(len=1) :: WriteCloudHeight_ASCII_c         ! Write out ESRI ASCII files of ash-cloud height?
      character(len=1) :: WriteCloudHeight_KML_c           ! Write out        KML files of ash-cloud height?
      character(len=1) :: WriteCloudLoad_ASCII_c           ! Write out ESRI ASCII files of ash-cloud load (T/km2) at specified times? 
      character(len=1) :: WriteCloudLoad_KML_c             ! Write out        KML files of ash-cloud load (T/km2) at specified times?
      character(len=1) :: WriteDepositTime_ASCII_c         ! Write out ESRI ASCII file of deposit arrival times?
      character(len=1) :: WriteDepositTime_KML_c           ! Write out        KML file of deposit arrival times?
      character(len=1) :: WriteCloudTime_ASCII_c           ! Write out ESRI ASCII file of cloud arrival times
      character(len=1) :: WriteCloudTime_KML_c             ! Write out        KML file of cloud arrival times?
      character(len=1) :: Write3dFiles_c                   ! Write out 3-D ash concentration at specified times?
      integer          :: ifm                              ! output code: 1=2d+concen,2=2d only]
      integer          :: ofm                              ! format of ash concentration files (1=ascii, 2=binary, or 3=netcdf)
      integer          :: nwt                              ! nWriteTimes
      logical          :: interval_flag                    ! indicates if nWriteTimes is actually specifying WriteTimes is interval-based
      real(kind=dp),dimension(:),allocatable :: wts        ! WriteTimes(1:nWriteTimes)

      ! Block 6 variables
      character(len=1) :: WriteAirportFile_ASCII_c
      character(len=1) :: WriteGSD_c
      character(len=1) :: WriteAirportFile_KML_c
      character(len=1) :: ProjectAirportLocations_c

      integer           :: iyear
      integer           :: idx
      character(len=4)  :: yearstr
      character(len=50) :: tmpstr

      INTERFACE
        subroutine Read_RunParam_Table(ir)
          integer,intent(in) :: ir
        end subroutine Read_RunParam_Table
        integer function HS_YearOfEvent(HoursSince,byear,useLeaps)
          real(kind=8),intent(in) :: HoursSince
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_YearOfEvent
      END INTERFACE

      ! Reset verbosity so we only are using stdout (no log file)
      VB = (/3,10/)
      MR_VERB = VB(1)

      nargs = command_argument_count()
      if (nargs.eq.0) then
          ! If no command-line arguments are given, then prompt user
          ! interactively for the input netcdf file
        if(VB(1).ge.verbosity_silent)then
          write(errlog(1),*)"Stdout is suppressed via VERB=9, but interactive input is expected."
          write(errlog(1),*)"Either recompile with VERB<9 or provide the correct command-line arguments."
          stop 1
        else
          do io=1,nio;if(VB(io).le.verbosity_production)then
            write(outlog(io),*)'Enter name of the template Ash3d control file'
          endif;enddo
        endif
        read(input_unit,*) infile

        do io=1,nio;if(VB(io).le.verbosity_production)then
          write(outlog(io),*)'Enter name of the run table file'
        endif;enddo
        read(input_unit,*)datafileIn

        do io=1,nio;if(VB(io).le.verbosity_production)then
          write(outlog(io),*)'Enter name of the run number.'
        endif;enddo
        read(input_unit,*)RunID

      elseif (nargs.gt.3) then
        do io=1,nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)'ERROR: Too many command-line arguments.'
          write(errlog(io),*)'  Usage: Ash3d_ASCII_GenCTR template.inp table.dat runID'
        endif;enddo
        stop 1
      else
        call get_command_argument(1, linebuffer080, status=stat)
        if(stat.gt.0)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Could not parse argument 1'
          endif;enddo
          stop 1
        elseif (stat.lt.0)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Argument 1 has been truncated.'
            write(errlog(io),*)'       File name length is limited to 80 char.'
          endif;enddo
          stop 1
        endif
        infile=trim(adjustl(linebuffer080))
        inquire( file=adjustl(trim(infile)), exist=IsThere )
        if (.not.IsThere)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Input file 1 could not be found'
          endif;enddo
          stop 1
        endif

        call get_command_argument(2, linebuffer080, status=stat)
        if(stat.gt.0)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Could not parse argument 1'
          endif;enddo
          stop 1
        elseif (stat.lt.0)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Argument 2 has been truncated.'
            write(errlog(io),*)'       File name length is limited to 80 char.'
          endif;enddo
          stop 1
        endif
        datafileIn=trim(adjustl(linebuffer080))
        inquire( file=adjustl(trim(datafileIn)), exist=IsThere )
        if (.not.IsThere)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Input file 2 could not be found'
          endif;enddo
          stop 1
        endif

        call get_command_argument(3, linebuffer080, status=stat)
        if(stat.gt.0)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Could not parse argument 1'
          endif;enddo
          stop 1
        endif
        read(linebuffer080,*)RunID

      endif

      call Read_Control_File

      ! Some error-checking
      if(neruptions.ne.1)then
        do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: neruptions>1'
            write(errlog(io),*)'       For now, this tool only applies to control files with a single eruptive pulse.'
        endif;enddo
        stop 1
      endif

      ! Template control file has been loaded into memory
      ! Note: any blocks for optional modules have not been read at this point

      ! Next step is to read and parse the table file. This will indicate how the
      ! output control file will be modified from the base configuration. This
      ! will mean replacing variables such as the start time, plume height, volume,
      ! vent location, etc.

      call Read_RunParam_Table(RunID)

      ! Now that we have loaded the data from the template file, start writing
      ! to a new control file
      infile = "temp.inp"
      inquire( file=infile, exist=IsThere )
      if(IsThere)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(errlog(io),*)"WARNING: Input file already exists"
        endif;enddo
      endif
      open(unit=fid_ctrlfile,file=infile,action='write')

      if(IsLatLon)then
        LLx    = lonLL
        LLy    = latLL
        widthx = gridwidth_e
        widthy = gridwidth_n
        x_in   = lon_volcano
        y_in   = lat_volcano
        z_in   = z_volcano
        dx_in  = de
        dy_in  = dn
      else
        LLx    = xLL
        LLy    = yLL
        widthx = gridwidth_x
        widthy = gridwidth_y
        x_in   = x_volcano
        y_in   = y_volcano
        z_in   = z_volcano
        dx_in  = dx
        dy_in  = dy
      endif
      if(VarDzType.eq.'dz_cons')then
        dz_type = 1
      elseif(VarDzType.eq.'dz_plin')then
        dz_type = 2
      elseif(VarDzType.eq.'dz_clog')then
        dz_type = 3
      elseif(VarDzType.eq.'dz_cust')then
        dz_type = 4
      else
        dz_type = 1
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  BLOCK 1: GRID INFO
      call Write_input_block_header(fid_ctrlfile,1)
      call SetWrite_input_block_01(fid_ctrlfile                    ,&  ! output stream ID
                                   VolcanoName                     ,&  ! volcano name
                                   cdf_b1l2                        ,&  ! projection line
                                   LLx,LLy                         ,&  ! x, y of LL corner
                                   widthx,widthy                   ,&  ! x and y width of domain
                                   x_in,y_in,z_in                  ,&  ! vent x,y,z
                                   dx_in,dy_in                     ,&  ! dx, dy
                                   dz_const                        ,&  ! dz
                                   diffusivity_horz,Suzuki_A       ,&  ! Kx, Suz param
                                   neruptions                      ,&  ! # of eruptions
                                   dz_type,cdf_vardz               ,&  ! dz_type + bonus line
                                   SourceType_idx)                     ! source_type 1=Suz,2=point,3=line,4=profile,5=umb,6=umb_air

      ! BLOCK 2: ERUPTION PARAMETERS
      call Write_input_block_header(fid_ctrlfile,2)
      call SetWrite_input_block_02(fid_ctrlfile                    ,&  ! output stream ID
                                   neruptions                      ,&  ! # of eruptions
                                   SourceType_idx                  ,&  ! source_type 1=Suz,2=point,3=line,4=profile,5=umb,6=umb_air
                                   e_StartTime+SimStartHour        ,&
                                   e_Duration                      ,&
                                   e_PlumeHeight                   ,&
                                   e_Volume                        ,&
                                   e_prof_dz                       ,&
                                   e_prof_nzpoints                 ,&
                                   e_prof_Volume)

      ! BLOCK 3: WIND PARAMETERS
      call Write_input_block_header(fid_ctrlfile,3)
      ! for iwind=5 cases, the number is windfiles is modifies, so read from cdf_b3l5
      read(cdf_b3l5,*,iostat=iostatus,iomsg=iomessage) nwindfiles
      call SetWrite_input_block_03(fid_ctrlfile                    ,&  ! output stream ID
                                   MR_iWind                        ,&
                                   MR_iWindFormat                  ,&
                                   MR_iGridCode                    ,&
                                   MR_iDataFormat                  ,&
                                   MR_iHeightHandler               ,&
                                   Simtime_in_hours                ,&
                                   StopWhenDeposited               ,&
                                   nwindfiles)

      ! BLOCK 4: OUTPUT OPTIONS
      WriteDepositFinal_ASCII_c = 'n'
      read(cdf_b4l1,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)).eq.'yes') WriteDepositFinal_ASCII_c = 'y'
      WriteDepositFinal_KML_c = 'n'
      read(cdf_b4l2,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)).eq.'yes') WriteDepositFinal_KML_c = 'y'
      WriteDepositTS_ASCII_c = 'n'
      read(cdf_b4l3,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)).eq.'yes') WriteDepositTS_ASCII_c = 'y'
      WriteDepositTS_KML_c = 'n'
      read(cdf_b4l4,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)).eq.'yes') WriteDepositTS_KML_c = 'y'
      WriteCloudConcentration_ASCII_c= 'n'
      read(cdf_b4l5,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)).eq.'yes') WriteCloudConcentration_ASCII_c = 'y'
      WriteCloudConcentration_KML_c= 'n'
      read(cdf_b4l6,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)).eq.'yes') WriteCloudConcentration_KML_c = 'y'
      WriteCloudHeight_ASCII_c = 'n'
      read(cdf_b4l7,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)).eq.'yes') WriteCloudHeight_ASCII_c = 'y'
      WriteCloudHeight_KML_c = 'n'
      read(cdf_b4l8,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)).eq.'yes') WriteCloudHeight_KML_c = 'y'
      WriteCloudLoad_ASCII_c = 'n'
      read(cdf_b4l9,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)).eq.'yes') WriteCloudLoad_ASCII_c = 'y'
      WriteCloudLoad_KML_c = 'n'
      read(cdf_b4l10,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)).eq.'yes') WriteCloudLoad_KML_c = 'y'
      WriteDepositTime_ASCII_c= 'n'
      read(cdf_b4l11,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)).eq.'yes') WriteDepositTime_ASCII_c = 'y'
      WriteDepositTime_KML_c= 'n'
      read(cdf_b4l12,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)).eq.'yes') WriteDepositTime_KML_c = 'y'
      WriteCloudTime_ASCII_c = 'n'
      read(cdf_b4l13,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)).eq.'yes') WriteCloudTime_ASCII_c = 'y'
      WriteCloudTime_KML_c = 'n'
      read(cdf_b4l14,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)).eq.'yes') WriteCloudTime_KML_c = 'y'
      Write3dFiles_c = 'y' ! This is obviously true if we are reading the netcdf file
      ifm = 2
      read(cdf_b4l15,*,iostat=iostatus,iomsg=iomessage) answer, ifm
      if(iostatus.ne.0)then ! if read fails, then make sure we set these
        Write3dFiles_c = 'y'
        ifm = 1
      endif
      ofm = 3 ! This should be 3 since we are reading a netcdf file
      read(cdf_b4l17,*)nwt
      if(nwt.gt.0)then
        interval_flag = .false.
        allocate(wts(nwt))
        read(cdf_b4l18,*)wts(1:nwt)
      else
        nwt = 1
        interval_flag = .true.
        allocate(wts(nwt))
        read(cdf_b4l18,*)wts(1)
      endif

      call Write_input_block_header(fid_ctrlfile,4)
      call SetWrite_input_block_04(fid_ctrlfile                    ,&  ! output stream ID
                                   WriteDepositFinal_ASCII_c       ,&  ! B4L1 n/y Write out ESRI ASCII file of final deposit thickness?
                                   WriteDepositFinal_KML_c         ,&  ! B4L2 n/y Write out        KML file of final deposit thickness?
                                   WriteDepositTS_ASCII_c          ,&  ! B4L3 Write out ESRI ASCII deposit files at specified times?
                                   WriteDepositTS_KML_c            ,&  ! B4L4 Write out        KML deposit files at specified times?
                                   WriteCloudConcentration_ASCII_c ,&  ! B4L5 Write out ESRI ASCII files of ash-cloud concentration?
                                   WriteCloudConcentration_KML_c   ,&  ! B4L6 Write out        KML files of ash-cloud concentration?
                                   WriteCloudHeight_ASCII_c        ,&  ! B4L7 Write out ESRI ASCII files of ash-cloud height?
                                   WriteCloudHeight_KML_c          ,&  ! B4L8 Write out        KML files of ash-cloud height?
                                   WriteCloudLoad_ASCII_c          ,&  ! B4L9 Write out ESRI ASCII files of ash-cloud load (T/km2) at specified times? 
                                   WriteCloudLoad_KML_c            ,&  ! B4L10 Write out        KML files of ash-cloud load (T/km2) at specified times?
                                   WriteDepositTime_ASCII_c        ,&  ! B4L11 Write out ESRI ASCII file of deposit arrival times?
                                   WriteDepositTime_KML_c          ,&  ! B4L12 Write out        KML file of deposit arrival times?
                                   WriteCloudTime_ASCII_c          ,&  ! B4L13 Write out ESRI ASCII file of cloud arrival times?
                                   WriteCloudTime_KML_c            ,&  ! B4L14 Write out        KML file of cloud arrival times?
                                   Write3dFiles_c                  ,&  ! B4L15 Write out 3-D ash concentration at specified times?
                                   ifm                             ,&  ! B4L15+ output code: 1=2d+concen,2=2d only]
                                   ofm                             ,&  ! B4L16 format of ash concentration files (1=ascii, 2=binary, or 3=netcdf)
                                   nwt                             ,&  ! B4L17 nWriteTimes
                                   interval_flag                   ,&  ! true if nwt sould really be -1 (e.g. write time = interval)
                                   wts)                                ! B4L18 WriteTimes(1:nWriteTimes)

      ! BLOCK 5: INPUT WIND FILES
      if(MR_iWind.eq.5)then
        nwindfiles = 1
        if(len(adjustl(trim(cdf_b5l1))).eq.0)then
          ! NetCDF file didn't have cdf_b5l1. Need to get windfile directory from MR_windfiles
          iyear = HS_YearOfEvent(SimStartHour,BaseYear,useLeap)
          write(yearstr,'(i4)')iyear
          idx=index(MR_windfiles(1),yearstr)
          write(tmpstr,*)MR_windfiles(1)(1:idx-2)
          write(MR_windfiles(1),*)adjustl(trim(tmpstr))
        else
          read(cdf_b5l1,*)MR_WindFiles(1)
        endif
      endif
      call Write_input_block_header(fid_ctrlfile,5)
      call SetWrite_input_block_05(fid_ctrlfile                    ,&  ! output stream ID
                                   nwindfiles                      ,&
                                   MR_WindFiles(1:nwindfiles))

      ! BLOCK 6: AIRPORT FILE
      WriteAirportFile_ASCII_c = 'n'
      read(cdf_b6l1,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)).eq.'yes') WriteAirportFile_ASCII_c = 'y'
      WriteGSD_c = 'n'
      read(cdf_b6l2,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)).eq.'yes') WriteGSD_c = 'y'
      WriteAirportFile_KML_c = 'n'
      read(cdf_b6l3,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)).eq.'yes') WriteAirportFile_KML_c = 'y'
      !cdf_b6l4
      ProjectAirportLocations_c = 'n'
      read(cdf_b6l5,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)).eq.'yes') ProjectAirportLocations_c = 'y'

      call Write_input_block_header(fid_ctrlfile,6)
      call SetWrite_input_block_06(fid_ctrlfile                    ,&  ! output stream ID
                                   WriteAirportFile_ASCII_c        ,&  ! Write out ash arrival times at airports to ASCII FILE?
                                   WriteGSD_c                      ,&  ! Write out grain-size distribution to ASCII airport file?
                                   WriteAirportFile_KML_c          ,&  ! Write out ash arrival times to kml file?
                                   cdf_b6l4(1:80)                  ,&  ! Name of file containing airport locations
                                   ProjectAirportLocations_c)          ! Defer to Lon/Lat coordinates? ("no" defers to projected)

      ! BLOCK 7: GRAIN-SIZE BINS, SETTLING VELOCITY
      call Write_input_block_header(fid_ctrlfile,7)
      call SetWrite_input_block_07(fid_ctrlfile                    ,&  ! output stream ID
                                   n_gs_max                        ,&  ! number of actual grain-size bins 
                                   FV_ID                           ,&  !
                                   Shape_ID                        ,&
                                   LN_massfrac                     ,&
                                   LN_phi_mean                     ,&
                                   LN_phi_stddev                   ,&
                                   Tephra_gsdiam                   ,&
                                   Tephra_bin_mass                 ,&
                                   Tephra_v_s                      ,&                   
                                   Tephra_rho_m                    ,&
                                   Tephra_gsF                      ,&
                                   Tephra_gsG                      ,&
                                   Tephra_gsPhi)

      ! BLOCK 8: VERTICAL PROFILES
      call Write_input_block_header(fid_ctrlfile,8)
      call SetWrite_input_block_08(fid_ctrlfile                    ,&  ! output stream ID
                                   nvprofiles                      ,&  ! number of profiles
                                   x_vprofile                      ,&  ! x coordinate of profile point
                                   y_vprofile                      ,&  ! y coordinate of profile point
                                   Site_vprofile)                      ! name of site


      ! BLOCK 9 (Optional): NETCDF ANNOTATIONS
      if(Have_Block_NetCDF)then
        call Write_input_block_header(fid_ctrlfile,9)
        call SetWrite_input_block_09(fid_ctrlfile                    ,&  ! output stream ID
                                     concenfile                      ,&  ! name of netcdf file
                                     cdf_title                       ,&  ! title of simulation
                                     cdf_comment)                        ! comment
      endif

      ! BLOCK 10+: OPTIONAL MODULES (RESETPARAMS)
      if(Have_Block_ResParm)then
        call Write_input_block_header(fid_ctrlfile,10)
        call SetWrite_input_block_ResetParam(fid_ctrlfile   ) !           ,&  ! output stream ID
      endif

      ! BLOCK 10+: OPTIONAL MODULES (TOPO)
      if(Have_Block_Topo)then
        call Write_input_block_header(fid_ctrlfile,11)
        call SetWrite_input_block_Topo(fid_ctrlfile                  ,&  ! output stream ID
                                       var_User_charlines_Topo(1)    ,& 
                                       var_User_charlines_Topo(2)    ,&
                                       var_User_charlines_Topo(3)    ,&
                                       var_User_charlines_Topo(4))
      endif

      ! BLOCK 10+: OPTIONAL MODULES (VARDIFF)
      if(Have_Block_VarDiff)then
        call Write_input_block_header(fid_ctrlfile,12)
        call SetWrite_input_block_VarDiff(fid_ctrlfile                  ,&  ! output stream ID
                                       var_User_charlines_VarDiff(1)    ,&
                                       var_User_charlines_VarDiff(2)    ,&
                                       var_User_charlines_VarDiff(3)    ,&
                                       var_User_charlines_VarDiff(4)    ,&
                                       var_User_charlines_VarDiff(5)    ,&
                                       var_User_charlines_VarDiff(6)    ,&
                                       var_User_charlines_VarDiff(7)    ,&
                                       var_User_charlines_VarDiff(8)    ,&
                                       var_User_charlines_VarDiff(9))
      endif


      close(fid_ctrlfile)

      end program Ash3d_ASCII_GenCTR

!##############################################################################

      subroutine Read_RunParam_Table(ir)

      use precis_param

      use io_units

      use io_data,       only : &
         datafileIn,VolcanoName

      use time_data,     only : &
         SimStartHour,BaseYear,useLeap

      use mesh,          only : &
         lonLL,latLL,gridwidth_e,gridwidth_n,de,dn,dz_const

      use Source,        only : &
         lon_volcano,lat_volcano,      &
         e_StartTime,e_Duration,e_PlumeHeight,e_Volume

      implicit none

      integer,intent(in) :: ir

      character(len=50) :: linebuffer050
      character(len=80) :: linebuffer080
      character(len=130):: linebuffer130
      character(len=130):: headerline
      character(len=130):: dataline
      integer           :: iostatus
      character(len=120):: iomessage

      integer,parameter :: MAX_COLVARS = 20
      integer,dimension(MAX_COLVARS) :: ivar_pos               = -1  ! position of variable in string
      integer,dimension(MAX_COLVARS) :: ivar_col               = -1  ! column number of variable in string
      integer,dimension(MAX_COLVARS) :: icol_var               = -1  ! variable ID for a given column
      integer           :: Ncols
      character(len=12),dimension(MAX_COLVARS) :: ivar_name1
      character(len=12),dimension(MAX_COLVARS) :: ivar_name2
      character(len=12),dimension(MAX_COLVARS) :: ivar_name3
      integer,dimension(MAX_COLVARS) :: ivar_Nnames
      integer,dimension(3)           :: itmp
      integer :: i,iv,iiv,ic,iline,irun
      integer :: itmp1,itmp2
      integer :: pos_cur,pos_diff
      logical :: HaveRunID,HaveST,HaveLoc
      real(kind=ip),dimension(MAX_COLVARS) :: values
      character(len=21) :: tmp_str1  ! Used to read start_time in 06-Sep-1996 15:37:36 format
      character(len=30) :: tmp_str2  ! Used to read location name (30 chars since VolcanoName has that length)

      integer :: iyear,imonth,iday,ihour,imin,isec
      character(len=3) :: monstr
      real(kind=8) :: hour

      INTERFACE
        real(kind=8) function HS_hours_since_baseyear(iyear,imonth,iday,hours,byear,useLeaps)
          integer     ,intent(in) :: iyear
          integer     ,intent(in) :: imonth
          integer     ,intent(in) :: iday
          real(kind=8),intent(in) :: hours
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_hours_since_baseyear
        real(kind=8)  function HS_HourOfDay(HoursSince,byear,useLeaps)
          real(kind=8),intent(in) :: HoursSince
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_HourOfDay
        integer function HS_YearOfEvent(HoursSince,byear,useLeaps)
          real(kind=8),intent(in) :: HoursSince
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_YearOfEvent
        integer function HS_MonthOfEvent(HoursSince,byear,useLeaps)
          real(kind=8),intent(in) :: HoursSince
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_MonthOfEvent
        integer function HS_DayOfEvent(HoursSince,byear,useLeaps)
          real(kind=8),intent(in) :: HoursSince
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_DayOfEvent
      END INTERFACE

      ! List the colome header variable names and synonyms
      ivar_Nnames( 1) = 2; ivar_name1( 1) = "run";          ivar_name2( 1) = "Run"
      ivar_Nnames( 2) = 3; ivar_name1( 2) = "year";         ivar_name2( 2) = "Year"; ivar_name3( 2) = "YYYY"
      ivar_Nnames( 3) = 2; ivar_name1( 3) = "month";        ivar_name2( 3) = "MM"
      ivar_Nnames( 4) = 2; ivar_name1( 4) = "day";          ivar_name2( 4) = "DD"
      ivar_Nnames( 5) = 2; ivar_name1( 5) = "hour";         ivar_name2( 5) = "HH.H"
      ivar_Nnames( 6) = 1; ivar_name1( 6) = "start time"
      ivar_Nnames( 7) = 1; ivar_name1( 7) = "Location"
      ivar_Nnames( 8) = 1; ivar_name1( 8) = "lonLL"
      ivar_Nnames( 9) = 2; ivar_name1( 9) = "latLL"
      ivar_Nnames(10) = 1; ivar_name1(10) = "dxy"
      ivar_Nnames(11) = 1; ivar_name1(11) = "dz"
      ivar_Nnames(12) = 3; ivar_name1(12) = "longitude";    ivar_name2(12) = "srcx"; ivar_name3(12) = "lon_volcano"
      ivar_Nnames(13) = 3; ivar_name1(13) = "latitude";     ivar_name2(13) = "srcy"; ivar_name3(13) = "lat_volcano"
      ivar_Nnames(14) = 2; ivar_name1(14) = "duration";     ivar_name2(14) = "EDur"
      ivar_Nnames(15) = 2; ivar_name1(15) = "plume height"; ivar_name2(15) = "EPlmH"
      ivar_Nnames(16) = 3; ivar_name1(16) = "volume";       ivar_name2(16) = "DRE";  ivar_name3(16) = "EVol"
      ivar_Nnames(17) = 2; ivar_name1(17) = "width";        ivar_name2(17) = "gridwidth_e"
      ivar_Nnames(18) = 2; ivar_name1(18) = "height";       ivar_name2(18) = "gridwidth_n"  ! Note: this must be listed after 'plume height'
      ivar_Nnames(19) = 1; ivar_name1(19) = "m_fines"
      ivar_Nnames(20) = 1; ivar_name1(20) = "mu_agg"


      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"*******************************************"
        write(outlog(io),*)"Now reading input table:"
      endif;enddo

      open(unit=fid_misc,file=datafileIn,status='old',action='read',err=9001)

      ! Reading the first header line. Should be something link: SUMMARY OF INPUT VALUES ...
      read(fid_misc,'(a130)',iostat=iostatus,iomsg=iomessage)linebuffer130
      linebuffer050 = "Reading table file for Line 1: header line"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer130(1:80),iomessage)

      ! Now read line two, which contains the column headers
      read(fid_misc,'(a130)',iostat=iostatus,iomsg=iomessage)linebuffer130
      linebuffer050 = "Reading table file for Line 2: column headers"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer130(1:80),iomessage)
      headerline = linebuffer130

      ! Now read line three, which contains the column units
      read(fid_misc,'(a130)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading table file for Line 3: column units"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080(1:80),iomessage)

      ! Parse the column header line
      Ncols = 0
      do iv = 1,MAX_COLVARS
        itmp(:) = -1
        itmp(1) = index(linebuffer130,trim(adjustl(ivar_name1(iv))))
        if(ivar_Nnames(iv).ge.2) itmp(2) = index(linebuffer130,trim(adjustl(ivar_name2(iv))))
        if(ivar_Nnames(iv).ge.3) itmp(3) = index(linebuffer130,trim(adjustl(ivar_name3(iv))))
        if(iv.eq.18)then ! Need to do an extra check if var=height since 'plume height' will catch it
          itmp1 = index(linebuffer130,'plume height')
          if (itmp(1).gt.0) then ! 'height' was found, but double-check
            if(itmp(1)-itmp1.eq.6)then ! the 'height' we found is a part of 'plume height'; look again
              itmp2 = index(linebuffer130(itmp(1)+1:),'height')
              if (itmp2.gt.0)ivar_pos(18) = itmp2
            else
              ivar_pos(18) = itmp(1)
            endif
          elseif (itmp(2).gt.0) then
            ivar_pos(18) = itmp(2) ! column found (gridwidth_n)
          endif
          if (ivar_pos(iv).gt.0) Ncols=Ncols+1
          cycle
        endif
        if (itmp(1).gt.0) then
          ivar_pos(iv) = itmp(1)
        elseif (itmp(2).gt.0) then
          ivar_pos(iv) = itmp(2)
        elseif (itmp(3).gt.0) then
          ivar_pos(iv) = itmp(3)
        endif
        if(ivar_pos(iv).gt.0)then
          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"Found column at position: ", ivar_pos(iv), trim(adjustl(ivar_name1(iv)))
          endif;enddo
        else
          if(iv.eq.1)then ! Run ID
            do io=1,2;if(VB(io).le.verbosity_debug1)then
              write(outlog(io),*)"     Column not found with: ",ivar_name1(iv)
              write(outlog(io),*)"     Using line number for run ID"
            endif;enddo
          endif
        endif
        if (ivar_pos(iv).gt.0) Ncols=Ncols+1
      enddo
      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"Found n columns",Ncols
      endif;enddo

      ! We've sorted out which variables are provided by the input table and the start position of
      ! variable name in the header. From this, we need to sort on this position and determine the
      ! order of variables. Since we only have maybe a dozen columns and only 20 variables to test,
      ! we'll use an N^2 algorithm.  Sorting aficionados should avert their eyes.

      ivar_col(:) = -1
      pos_cur  = 0   ! current position
      pos_diff = 130 ! set large initial difference
      ! Loop through variables and find the one with the closest next position
      do ic=1,Ncols
        do iiv=1,MAX_COLVARS
          if (ivar_pos(iiv).gt.pos_cur.and.&  ! if this is a var found in table and not yet flagged
              ivar_pos(iiv)-pos_cur.lt.pos_diff)then
            ivar_col(iiv) = ic
            icol_var(ic)  = iiv
            pos_diff = ivar_pos(iiv)-pos_cur
          endif
        enddo
        ! reset the start position and position difference
        pos_cur = ivar_pos(icol_var(ic))+1
        pos_diff = 130
      enddo

      do io=1,2;if(VB(io).le.verbosity_info)then
        do iv=1,Ncols
          write(outlog(io),*)iv,icol_var(iv),ivar_name1(icol_var(iv))
        enddo
      endif;enddo

      ! Error-checking of columns
      ! We have requirements:
      !  1 ) if RunID is present, it must be in the first column
      !  2 ) if start_time is present, it must be either after RunID (col 2) or
      !      in col 1 using up to position 29 of the string. Remaining cols with be in 30+
      !  3 ) if Location is present, it must be the last column of the file
      ! Check if RunID exists and make sure it is in the first column
      HaveRunID = .false.
      if(ivar_pos(1).gt.0)then
        HaveRunID=.true.
        if(icol_var(1).eq.1)then
          ! RunID is found and in the first column
          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"Requested run number will be identified from column 1 of table."
          endif;enddo
        else
          ! RunID is found but not in the first column
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"Run number is not in the first column. This will cause problems in reading"
            write(errlog(io),*)"the column, unfortunately. Please move or remove the 'run #' column."
            stop 1
          endif;enddo
        endif
      else
        ! RunID is not found. Use the line number as a proxy for run #
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Requested run number will be identified from the line number of the table."
        endif;enddo
      endif
      ! Check if start_time is provided
      HaveST = .false.
      if(ivar_pos(6).gt.0)then
        if(HaveRunID)then
          if(ivar_col(6).ne.2)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: start_time is present but not in col 2 after RunID."
            endif;enddo
            stop 1
          else
            HaveST = .true.
          endif
        else
          HaveST = .true.
          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"start_time is present and in col 1."
          endif;enddo
        endif
      endif
      ! Check if Location is provided
      HaveLoc = .false.
      if(ivar_pos(7).gt.0)then
        HaveLoc = .true.
        if(ivar_col(7).ne.Ncols)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: Location is present but not in last column."
            write(errlog(io),*)"       This will break the current line parsing."
          endif;enddo
          stop 1
        endif
      endif

      ! Start reading the 
      read(fid_misc,'(a130)',iostat=iostatus,iomsg=iomessage)linebuffer130
      linebuffer050 = "Reading table file for Line 1: header line"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer130(1:80),iomessage)
      iline = 0
      irun  = 0
      do while(iostatus.eq.0)
        iline = iline + 1
        if(ivar_pos(1).gt.0)then
          ! If the RunID from the table is what we want, exit do loop
          read(linebuffer130,*)irun
          if(irun.eq.ir) exit
        else
          ! If the line number (as a proxy for RunID) is what we want, exit do loop
          if(iline.eq.ir) exit
        endif
        read(fid_misc,'(a130)',iostat=iostatus,iomsg=iomessage)linebuffer130
      enddo
      if(ir.gt.iline.and.ir.gt.irun)then
        ! Too large of a RunID was requested
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"Could not find the requested run ID"
        endif;enddo
        stop 1
      endif
      dataline = linebuffer130

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)headerline
        write(outlog(io),*)dataline
      endif;enddo
      if(HaveRunID)then
        if(HaveST)then
          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"Reading irun start_time + other columns"
          endif;enddo
          ! First read the run # and then the start time in character format. We will
          ! have to interpret that later.
          read(dataline,'(i8,a21)')irun,tmp_str1
          ! For the rest of the values, we will read as reals starting from position 30
          if(HaveLoc)then
            read(dataline(30:),*)values(3:Ncols-1)
          else
            read(dataline(30:),*)values(3:Ncols)
          endif
        else
          ! No start_time string so we can read everything as reals except a possible site label
          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"Reading irun + other columns",2,Ncols-1
          endif;enddo
          if(HaveLoc)then
            read(dataline,*)irun,values(2:Ncols-1)
          else
            read(dataline,*)irun,values(2:Ncols)
          endif
        endif
      else
        irun = iline
        if(HaveST)then
          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"Reading start_time + other columns"
          endif;enddo
          ! First read the run # and then the start time in character format. We will
          ! have to interpret that later.
          read(dataline,'(a21)')tmp_str1
          ! For the rest of the values, we will read as reals starting from position 30
          if(HaveLoc)then
            read(dataline(30:),*)values(2:Ncols-1)
          else
            read(dataline(30:),*)values(2:Ncols)
          endif
        else
          ! No start_time string so we can read everything as reals except a possible site label
          do io=1,2;if(VB(io).le.verbosity_debug1)then
            write(outlog(io),*)"Reading other columns",2,Ncols-1
          endif;enddo
          if(HaveLoc)then
            read(dataline,*)values(1:Ncols-1)
          else
            read(dataline,*)values(1:Ncols)
          endif
        endif
      endif

      ! Finally, read the Site label if present
      if(HaveLoc)then
        ! We need the start position of the data column, but we only have the start position of the label.
        ! The strategy here is the test the label position in the data string for a space. If so, start reading
        ! from that point. Otherwise, march backwards up to five positions until a space is found.
        itmp1=ivar_pos(7)
        if(dataline(itmp1:itmp1).eq.' ')then
          read(dataline(itmp1:),*)tmp_str2
        else
          do i=1,5
            itmp1 = itmp1-1
            if(dataline(itmp1:itmp1).eq.' ')exit
          enddo
          read(dataline(itmp1:),*)tmp_str2
        endif
      endif

      ! Get start year as specified in the template file
      iyear = HS_YearOfEvent(e_StartTime(1)+SimStartHour,BaseYear,useLeap)
      imonth= HS_MonthOfEvent(e_StartTime(1)+SimStartHour,BaseYear,useLeap)
      iday  = HS_DayOfEvent(e_StartTime(1)+SimStartHour,BaseYear,useLeap)
      hour  = HS_HourOfDay(e_StartTime(1)+SimStartHour,BaseYear,useLeap)

      if(ivar_pos( 1).gt.0) write(*,*)"irun         : ",ivar_name1( 1),irun
      if(ivar_pos( 2).gt.0)then
        iyear = int(values(ivar_col( 2)))
        write(*,*)"year         : ",ivar_name1( 2),values(ivar_col( 2))
      endif
      if(ivar_pos( 3).gt.0)then
        imonth = int(values(ivar_col( 3)))
        write(*,*)"month        : ",ivar_name1( 3),values(ivar_col( 3))
      endif
      if(ivar_pos( 4).gt.0)then
        iday = int(values(ivar_col( 4)))
        write(*,*)"day          : ",ivar_name1( 4),values(ivar_col( 4))
      endif
      if(ivar_pos( 5).gt.0)then
        hour = values(ivar_col( 5))
        write(*,*)"hour         : ",ivar_name1( 5),values(ivar_col( 5))
      endif
      if(ivar_pos( 6).gt.0)then
        write(*,*)"start time   : ",ivar_name1( 6),trim(adjustl(tmp_str1))
        tmp_str1=trim(adjustl(tmp_str1))
        read(tmp_str1,51)iday,monstr,iyear,ihour,imin,isec
 51     format(i2,1x,a3,1x,i4,1x,i2,1x,i2,1x,i2)
        if(monstr.eq.'Jan')then
          imonth = 1
        elseif(monstr.eq.'Feb')then
          imonth = 2
        elseif(monstr.eq.'Mar')then
          imonth = 3
        elseif(monstr.eq.'Apr')then
          imonth = 4
        elseif(monstr.eq.'May')then
          imonth = 5
        elseif(monstr.eq.'Jun')then
          imonth = 6
        elseif(monstr.eq.'Jul')then
          imonth = 7
        elseif(monstr.eq.'Aug')then
          imonth = 8
        elseif(monstr.eq.'Sep')then
          imonth = 9
        elseif(monstr.eq.'Oct')then
          imonth = 10
        elseif(monstr.eq.'Nov')then
          imonth = 11
        elseif(monstr.eq.'Dec')then
          imonth = 12
        else
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: Cannot parse month of start_time"
          endif;enddo
          stop 1
        endif
        hour = real(ihour,kind=8) + real(imin,kind=8)/60.0_8 + real(isec,kind=8)/3600.0_8
      endif
      ! Check if we need to reset the start time by checking if any time value was reset
      if(ivar_pos( 2).gt.0.or.&
         ivar_pos( 3).gt.0.or.&
         ivar_pos( 4).gt.0.or.&
         ivar_pos( 5).gt.0.or.&
         ivar_pos( 6).gt.0)then
        SimStartHour   = HS_hours_since_baseyear(iyear,imonth,iday,hour,BaseYear,useLeap)
        e_StartTime(1) = 0.0_8
      endif
      if(ivar_pos( 7).gt.0)then
        write(*,*)"Location     : ",ivar_name1( 7),trim(adjustl(tmp_str2))
        VolcanoName = trim(adjustl(tmp_str2))
      endif
      if(ivar_pos( 8).gt.0)then
        lonLL = values(ivar_col( 8))
        write(*,*)"lonLL        : ",ivar_name1( 8),values(ivar_col( 8))
      endif
      if(ivar_pos( 9).gt.0)then
        latLL = values(ivar_col( 9))
        write(*,*)"latLL        : ",ivar_name1( 9),values(ivar_col( 9))
      endif
      if(ivar_pos(10).gt.0)then
        de = values(ivar_col(10))
        dn = values(ivar_col(10))
        write(*,*)"dxy          : ",ivar_name1(10),values(ivar_col(10))
      endif
      if(ivar_pos(11).gt.0)then
        dz_const = values(ivar_col(11))
        write(*,*)"dz           : ",ivar_name1(11),values(ivar_col(11))
      endif
      if(ivar_pos(12).gt.0)then
        lon_volcano = values(ivar_col(12))
        write(*,*)"longitude    : ",ivar_name1(12),values(ivar_col(12))
      endif
      if(ivar_pos(13).gt.0)then
        lat_volcano = values(ivar_col(13))
        write(*,*)"latitude     : ",ivar_name1(13),values(ivar_col(13))
      endif
      if(ivar_pos(14).gt.0)then
        e_Duration(1) = values(ivar_col(14))
        write(*,*)"duration     : ",ivar_name1(14),values(ivar_col(14))
      endif
      if(ivar_pos(15).gt.0)then
        e_PlumeHeight(1) = values(ivar_col(15))
        write(*,*)"plume height : ",ivar_name1(15),values(ivar_col(15))
      endif
      if(ivar_pos(16).gt.0)then
        e_Volume(1) = values(ivar_col(16))
        write(*,*)"volume       : ",ivar_name1(16),values(ivar_col(16))
      endif
      if(ivar_pos(17).gt.0)then
        gridwidth_e = values(ivar_col(17))
        write(*,*)"width        : ",ivar_name1(17),values(ivar_col(17))
      endif
      if(ivar_pos(18).gt.0)then
        gridwidth_n = values(ivar_col(18))
        write(*,*)"height       : ",ivar_name1(18),values(ivar_col(18))
      endif
      if(ivar_pos(19).gt.0)then
        write(*,*)"m_fines      : ",ivar_name1(19),values(ivar_col(19))
      endif
      if(ivar_pos(20).gt.0)then
        write(*,*)"mu_agg       : ",ivar_name1(20),values(ivar_col(20))
      endif

      close(fid_misc)

      return

9001  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error: cannot open table file: ',datafileIn
        write(errlog(io),*)  'Program stopped'
      endif;enddo
      stop 1

      end subroutine Read_RunParam_Table

!##############################################################################

