! Ash3d_NetCDF_GenCTR
! This tool will recover the control file that was used (or could be used) from
! an Ash3d NetCDF output file. This tool relies on some recent additions to
! the NetCDF output file and my not work for older files.

      program Ash3d_NetCDF_GenCTR

      use precis_param

      use io_units

      use io_data,       only : &
         infile,concenfile,cdf_title,cdf_comment,VolcanoName, &
         Have_Block_NetCDF,Have_Block_ResParm,Have_Block_Topo,Have_Block_VarDiff,&
         cdf_b1l2,cdf_vardz,cdf_b3l5, &
         cdf_b4l1,cdf_b4l2,cdf_b4l3,cdf_b4l4,cdf_b4l5,cdf_b4l6,cdf_b4l7,cdf_b4l8,cdf_b4l9,cdf_b4l10,&
         cdf_b4l11,cdf_b4l12,cdf_b4l13,cdf_b4l14,cdf_b4l15,cdf_b4l17,cdf_b4l18,cdf_b5l1,  &
         cdf_b6l1,cdf_b6l2,cdf_b6l3,cdf_b6l4,cdf_b6l5, &
         nvprofiles,Site_vprofile,x_vprofile, y_vprofile

      use mesh,          only : &
         lonLL,latLL,gridwidth_e,gridwidth_n,xLL,yLL,gridwidth_x,gridwidth_y, &
         IsLatLon,de,dn,dx,dy,VarDzType,dz_const

      use solution,      only : &
         StopWhenDeposited

      use time_data,     only : &
         SimStartHour,Simtime_in_hours,BaseYear,useLeap

      use Source,        only : &
         lon_volcano,lat_volcano,x_volcano,y_volcano,z_volcano,Suzuki_A,      &
         neruptions,SourceType_idx,e_StartTime,e_Duration,e_delta,         &
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

#ifdef USENETCDF
      use Ash3d_Netcdf_IO,  only : &
           NC_Read_Output_Products
#endif

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
         MR_iWind,MR_iWindFormat,MR_iGridCode,MR_iDataFormat,MR_iHeightHandler,MR_WindFiles

      implicit none
      !implicit none (type, external)

      integer            :: nargs
      integer            :: stat
      character(len= 80) :: linebuffer080

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
      integer             :: iostatus
      character (len= 50) :: iomessage
      character (len=  3) :: answer
      character (len=  1) :: WriteDepositFinal_ASCII_c        ! n/y Write out ESRI ASCII file of final deposit thickness?
      character (len=  1) :: WriteDepositFinal_KML_c          ! n/y Write out        KML file of final deposit thickness?
      character (len=  1) :: WriteDepositTS_ASCII_c           ! Write out ESRI ASCII deposit files at specified times?
      character (len=  1) :: WriteDepositTS_KML_c             ! Write out        KML deposit files at specified times?
      character (len=  1) :: WriteCloudConcentration_ASCII_c  ! Write out ESRI ASCII files of ash-cloud concentration?
      character (len=  1) :: WriteCloudConcentration_KML_c    ! Write out        KML files of ash-cloud concentration?
      character (len=  1) :: WriteCloudHeight_ASCII_c         ! Write out ESRI ASCII files of ash-cloud height?
      character (len=  1) :: WriteCloudHeight_KML_c           ! Write out        KML files of ash-cloud height?
      character (len=  1) :: WriteCloudLoad_ASCII_c           ! Write out ESRI ASCII files of ash-cloud load (T/km2) at spec times?
      character (len=  1) :: WriteCloudLoad_KML_c             ! Write out        KML files of ash-cloud load (T/km2) at spec times?
      character (len=  1) :: WriteDepositTime_ASCII_c         ! Write out ESRI ASCII file of deposit arrival times?
      character (len=  1) :: WriteDepositTime_KML_c           ! Write out        KML file of deposit arrival times?
      character (len=  1) :: WriteCloudTime_ASCII_c           ! Write out ESRI ASCII file of cloud arrival times
      character (len=  1) :: WriteCloudTime_KML_c             ! Write out        KML file of cloud arrival times?
      character (len=  1) :: Write3dFiles_c                   ! Write out 3-D ash concentration at specified times?
      integer             :: ifm                              ! output code: 1=2d+concen,2=2d only]
      integer             :: ofm                              ! format of ash concentration files (1=ascii, 2=binary, or 3=netcdf)
      integer             :: nwt                              ! nWriteTimes
      logical             :: interval_flag                    ! indicates if nWriteTimes is specifying WriteTimes is interval-based
      real(kind=dp),dimension(:),allocatable :: wts           ! WriteTimes(1:nWriteTimes)

      ! Block 6 variables
      character (len=  1)  :: WriteAirportFile_ASCII_c
      character (len=  1)  :: WriteGSD_c
      character (len=  1)  :: WriteAirportFile_KML_c
      character (len=  1)  :: ProjectAirportLocations_c

      integer              :: iyear
      integer              :: idx
      character (len=  4)  :: yearstr
      character (len= 50)  :: tmpstr

      integer :: i
      character (len=130) :: linebuffer130

      INTERFACE
        integer function HS_YearOfEvent(HoursSince,byear,useLeaps)
          implicit none
          !implicit none (type, external)
          integer        ,parameter   :: dp        = 8  ! double precision
          real(kind=dp)  ,intent(in)  :: HoursSince
          integer        ,intent(in)  :: byear
          logical        ,intent(in)  :: useLeaps
        end function HS_YearOfEvent
      END INTERFACE

      ! Reset verbosity so we only are using stdout (no log file)
      VB = [3,10]

      nargs = command_argument_count()
      if (nargs == 0) then
          ! If no command-line arguments are given, then prompt user
          ! interactively for the input netcdf file
        if(VB(1) >= verbosity_silent)then
          write(errlog(1),*)"Stdout is suppressed via VERB=9, but interactive input is expected."
          write(errlog(1),*)"Either recompile with VERB<9 or provide the correct command-line arguments."
          stop 1
        else
          do io=1,nio;if(VB(io) <= verbosity_production)then
            write(outlog(io),*)'Ash3d_NetCDF_GenCTR is a tool used to read an Ash3d output NetCDF'
            write(outlog(io),*)'and recreate the ASCII control file used in generating that output'
            write(outlog(io),*)'file. This tool take a single argument and writes out temp.inp.'
            write(outlog(io),*)'  Usage: Ash3d_NetCDF_GenCTR 3d_tephra_fall.nc'
            write(outlog(io),*)' '
            write(outlog(io),*)'No command-line arguments were provided, now entering interactive mode.'
            write(outlog(io),*)' '
            write(outlog(io),*)'Enter name of the Ash3d netcdf output file'
          endif;enddo
        endif
        read(input_unit,*) concenfile
      elseif (nargs > 1) then
        do io=1,nio;if(VB(io) <= verbosity_error)then
          write(errlog(io),*)'ERROR: Too many command-line arguments.'
          write(errlog(io),*)'  Usage: Ash3d_NetCDF_GenCTR output_file'
        endif;enddo
        stop 1
      else
        call get_command_argument(1, linebuffer080, status=stat)
        if(stat > 0)then
          do io=1,nio;if(VB(io) <= verbosity_error)then
            write(errlog(io),*)'ERROR: Could not parse argument 1'
          endif;enddo
          stop 1
        elseif (stat < 0)then
          do io=1,nio;if(VB(io) <= verbosity_error)then
            write(errlog(io),*)'ERROR: Argument 1 has been truncated.'
            write(errlog(io),*)'       File name length is limited to 80 char.'
          endif;enddo
          stop 1
        endif
        concenfile=trim(adjustl(linebuffer080))
        inquire( file=adjustl(trim(concenfile)), exist=IsThere )
        if (.not.IsThere)then
          do io=1,nio;if(VB(io) <= verbosity_error)then
            write(errlog(io),*)'ERROR: Input file 1 could not be found'
          endif;enddo
          stop 1
        endif
      endif

#ifdef USENETCDF
      ! Just read step 1. This brings in all the header info needed for the control file
      call NC_Read_Output_Products(1)
#else
     do io=1,nio;if(VB(io) <= verbosity_info)then
       write(errlog(io),*)'ERROR: NetCDF libraries not linked.'
       write(errlog(io),*)'       Please recompile linking NetCDF libraries.'
     endif;enddo
     stop 1
#endif

      ! Now that we have loaded the data from the NetCDF file, start writing
      ! to a control file
      infile = "temp.inp"
      inquire( file=infile, exist=IsThere )
      if(IsThere)then
        do io=1,2;if(VB(io) <= verbosity_info)then
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
      if(VarDzType == 'dz_cons')then
        dz_type = 1
      elseif(VarDzType == 'dz_plin')then
        dz_type = 2
      elseif(VarDzType == 'dz_clog')then
        dz_type = 3
      elseif(VarDzType == 'dz_cust')then
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
                                   e_delta                         ,&
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
      if(adjustl(trim(answer)) == 'yes') WriteDepositFinal_ASCII_c = 'y'
      WriteDepositFinal_KML_c = 'n'
      read(cdf_b4l2,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)) == 'yes') WriteDepositFinal_KML_c = 'y'
      WriteDepositTS_ASCII_c = 'n'
      read(cdf_b4l3,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)) == 'yes') WriteDepositTS_ASCII_c = 'y'
      WriteDepositTS_KML_c = 'n'
      read(cdf_b4l4,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)) == 'yes') WriteDepositTS_KML_c = 'y'
      WriteCloudConcentration_ASCII_c= 'n'
      read(cdf_b4l5,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)) == 'yes') WriteCloudConcentration_ASCII_c = 'y'
      WriteCloudConcentration_KML_c= 'n'
      read(cdf_b4l6,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)) == 'yes') WriteCloudConcentration_KML_c = 'y'
      WriteCloudHeight_ASCII_c = 'n'
      read(cdf_b4l7,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)) == 'yes') WriteCloudHeight_ASCII_c = 'y'
      WriteCloudHeight_KML_c = 'n'
      read(cdf_b4l8,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)) == 'yes') WriteCloudHeight_KML_c = 'y'
      WriteCloudLoad_ASCII_c = 'n'
      read(cdf_b4l9,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)) == 'yes') WriteCloudLoad_ASCII_c = 'y'
      WriteCloudLoad_KML_c = 'n'
      read(cdf_b4l10,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)) == 'yes') WriteCloudLoad_KML_c = 'y'
      WriteDepositTime_ASCII_c= 'n'
      read(cdf_b4l11,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)) == 'yes') WriteDepositTime_ASCII_c = 'y'
      WriteDepositTime_KML_c= 'n'
      read(cdf_b4l12,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)) == 'yes') WriteDepositTime_KML_c = 'y'
      WriteCloudTime_ASCII_c = 'n'
      read(cdf_b4l13,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)) == 'yes') WriteCloudTime_ASCII_c = 'y'
      WriteCloudTime_KML_c = 'n'
      read(cdf_b4l14,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)) == 'yes') WriteCloudTime_KML_c = 'y'
      Write3dFiles_c = 'y'   ! This is obviously true if we are reading the netcdf file
      ifm = 2
      read(cdf_b4l15,*,iostat=iostatus,iomsg=iomessage) answer, ifm
      if(iostatus /= 0)then  ! if read fails, then make sure we set these
        Write3dFiles_c = 'y'
        ifm = 1
      endif
      ofm = 3  ! This should be 3 since we are reading a netcdf file
      read(cdf_b4l17,*)nwt
      if(nwt > 0)then
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
      if(MR_iWind == 5)then
        nwindfiles = 1
        if(len(adjustl(trim(cdf_b5l1))) == 0)then
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
      do i=1,nwindfiles
        linebuffer130 = MR_WindFiles(i)
      enddo
      call SetWrite_input_block_05(fid_ctrlfile                    ,&  ! output stream ID
                                   nwindfiles                      ,&
                                   MR_WindFiles(1:nwindfiles))

      ! BLOCK 6: AIRPORT FILE
      WriteAirportFile_ASCII_c = 'n'
      read(cdf_b6l1,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)) == 'yes') WriteAirportFile_ASCII_c = 'y'
      WriteGSD_c = 'n'
      read(cdf_b6l2,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)) == 'yes') WriteGSD_c = 'y'
      WriteAirportFile_KML_c = 'n'
      read(cdf_b6l3,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)) == 'yes') WriteAirportFile_KML_c = 'y'
      !cdf_b6l4
      ProjectAirportLocations_c = 'n'
      read(cdf_b6l5,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)) == 'yes') ProjectAirportLocations_c = 'y'

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
                                   FV_ID                           ,&
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
        call SetWrite_input_block_ResetParam(fid_ctrlfile   )        !,&  ! output stream ID
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

      end program Ash3d_NetCDF_GenCTR

