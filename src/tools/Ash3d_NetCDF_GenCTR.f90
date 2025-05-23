! Ash3d_NetCDF_GenCTR

      program Ash3d_NetCDF_GenCTR

      use precis_param

      use io_units

      use io_data,       only : &
         infile,concenfile,VolcanoName, &
         nWriteTimes,WriteTimes,cdf_b1l2,cdf_vardz,cdf_b3l5, &
         cdf_b4l1,cdf_b4l2,cdf_b4l3,cdf_b4l4,cdf_b4l5,cdf_b4l6,cdf_b4l7,cdf_b4l8,cdf_b4l9,cdf_b4l10,&
         cdf_b4l11,cdf_b4l12,cdf_b4l13,cdf_b4l14,cdf_b4l15,cdf_b4l16,cdf_b4l17,cdf_b4l18,cdf_b5l1,  &
         cdf_b6l1,cdf_b6l2,cdf_b6l3,cdf_b6l4,cdf_b6l5

      use mesh,          only : &
         lonLL,latLL,gridwidth_e,gridwidth_n,xLL,yLL,gridwidth_x,gridwidth_y, &
         IsLatLon,de,dn,dx,dy,VarDzType,dz_const

      use solution,      only : &
         StopWhenDeposited

      use time_data,     only : &
         SimStartHour,Simtime_in_hours

      use Source,        only : &
         lon_volcano,lat_volcano,x_volcano,y_volcano,z_volcano,Suzuki_A,      &
         neruptions,SourceType,SourceType_idx,e_StartTime,e_Duration,         &
         e_PlumeHeight,e_Volume,e_prof_dz,e_prof_nzpoints,e_prof_Volume

      use Diffusion,     only : &
         diffusivity_horz

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
         MR_iWind,MR_iWindFormat,MR_iGridCode,MR_iDataFormat,MR_iHeightHandler,MR_iWindFiles,&
         MR_WindFiles

      implicit none

      integer           :: nargs
      integer           :: stat
      character(len=80) :: linebuffer080

      logical       :: IsThere
      logical       :: WriteBlock = .true.
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
      real(kind=dp),dimension(:),allocatable :: wts        ! WriteTimes(1:nWriteTimes)

      ! Block 6 variables
      character(len=1) :: WriteAirportFile_ASCII_c
      character(len=1) :: WriteGSD_c
      character(len=1) :: WriteAirportFile_KML_c
      character(len=1) :: ProjectAirportLocations_c

      write(*,*)"In Ash3d_NetCDF_GenCTR"

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
            write(outlog(io),*)'Enter name of the Ash3d netcdf output file'
          endif;enddo
        endif
        read(input_unit,*) concenfile
      elseif (nargs.gt.1) then
        do io=1,nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)'ERROR: Too many command-line arguments.'
          write(errlog(io),*)'  Usage: Ash3d_NetCDF_GenCTR output_file'
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
        concenfile=trim(adjustl(linebuffer080))
        inquire( file=adjustl(trim(concenfile)), exist=IsThere )
        if (.not.IsThere)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Input file 1 could not be found'
          endif;enddo
          stop 1
        endif
      endif

#ifdef USENETCDF
      !call NC_Read_Output_Products(-1)
      call NC_Read_Output_Products(1)
#endif  

      ! Now that we have loaded the data from the NetCDF file, start writing
      ! to a control file
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

      call Write_input_block_header(1)
      call SetWrite_input_block_01(WriteBlock                      ,&  ! indicates that write to stdout as well as set vars
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

      call Write_input_block_header(2)
      call SetWrite_input_block_02(WriteBlock                      ,&  ! indicates that write to stdout as well as set vars
                                   neruptions                      ,&  ! # of eruptions
                                   SourceType_idx                  ,&  ! source_type 1=Suz,2=point,3=line,4=profile,5=umb,6=umb_air
                                   e_StartTime+SimStartHour        ,&
                                   e_Duration                      ,&
                                   e_PlumeHeight                   ,&
                                   e_Volume                        ,&
                                   e_prof_dz                       ,&
                                   e_prof_nzpoints                 ,&
                                   e_prof_Volume)

      call Write_input_block_header(3)

      ! for iwind=5 cases, the number is windfiles is modifies, so read from cdf_b3l5
      read(cdf_b3l5,'(i2)',iostat=iostatus,iomsg=iomessage) nwindfiles
      call SetWrite_input_block_03(WriteBlock                      ,&  ! indicates that write to stdout as well as set vars
                                   MR_iWind                        ,&
                                   MR_iWindFormat                  ,&
                                   MR_iGridCode                    ,&
                                   MR_iDataFormat                  ,&
                                   MR_iHeightHandler               ,&
                                   Simtime_in_hours                ,&
                                   StopWhenDeposited               ,&
                                   nwindfiles)

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
      read(cdf_b4l15,'(a3)',iostat=iostatus,iomsg=iomessage) answer, ifm
      if(iostatus.ne.0)then ! if read fails, then make sure we set these
        Write3dFiles_c = 'y'
        ifm = 2
      endif
      ofm = 3 ! This should be 3 since we are reading a netcdf file
      nwt = nWriteTimes
      allocate(wts(nwt))
      wts = WriteTimes
 
      call Write_input_block_header(4)
      call SetWrite_input_block_04(WriteBlock                      ,&  ! indicates that write to stdout as well as set vars
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
                                   wts)                                ! B4L18 WriteTimes(1:nWriteTimes)

      if(MR_iWind.eq.5)then
        read(cdf_b5l1,*)MR_WindFiles(1)
      endif

      call Write_input_block_header(5)
      call SetWrite_input_block_05(WriteBlock                      ,&  ! indicates that write to stdout as well as set vars
                                   nwindfiles                      ,&
                                   MR_WindFiles(1:nwindfiles))

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

      call Write_input_block_header(6)
      call SetWrite_input_block_06(WriteBlock                      ,&  ! indicates that write to stdout as well as set vars
                                   WriteAirportFile_ASCII_c        ,&  ! Write out ash arrival times at airports to ASCII FILE?
                                   WriteGSD_c                      ,&  ! Write out grain-size distribution to ASCII airport file?
                                   WriteAirportFile_KML_c          ,&  ! Write out ash arrival times to kml file?
                                   cdf_b6l4(1:80)                  ,&  ! Name of file containing airport locations
                                   ProjectAirportLocations_c)          ! Defer to Lon/Lat coordinates? ("no" defers to projected)

!      call Write_input_block_header(7)
!      call SetWrite_input_block_07(WriteBlock           ) !           ,&  ! indicates that write to stdout as well as set vars
!
!      call Write_input_block_header(8)
!      call SetWrite_input_block_08(WriteBlock           ) !           ,&  ! indicates that write to stdout as well as set vars
!
!      call Write_input_block_header(9)
!      call SetWrite_input_block_09(WriteBlock           ) !           ,&  ! indicates that write to stdout as well as set vars
!
!      call Write_input_block_header(10)
!      call SetWrite_input_block_ResetParam(WriteBlock   ) !           ,&  ! indicates that write to stdout as well as set vars
!
!      call Write_input_block_header(11)
!      call SetWrite_input_block_Topo(WriteBlock         ) !           ,&  ! indicates that write to stdout as well as set vars
!
!      call Write_input_block_header(12)
!      call SetWrite_input_block_VarDiff(WriteBlock      ) !           ,&  ! indicates that write to stdout as well as set vars


      close(fid_ctrlfile)

      end program Ash3d_NetCDF_GenCTR

