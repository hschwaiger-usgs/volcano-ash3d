!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Ash3d_PostProc
!
!  This program is used for post-processing Ash3d output, primarily generating
!  output products from the netcdf output file.  The normal output products
!  available when running Ash3d, such as the kml files, 2-d ASCII output, 3-d
!  binary or vertical profile data, are all reproducible from the netcdf
!  output file using this program.  Additionally, this program can create kml,
!  ASCII, binary, shapefiles or contour maps (in png format) for many 2d variables
!  not available as options while running Ash3d.
!  
!  If no command-line arguments are provided, the program will run interactivley,
!  prompting the user for the name of the netcdf file to process, followed by
!  the output variable to process, the output format and the time step (if processing
!  a transient variable).
!  
!  If one command-line argument is provided, it is assumed to be a control file
!  that defines the input type and output products.  This allows the greatest flexibility.
!  
!  If more than one command-line argument is given, this program expects the following:
!    ./Ash3d_PostProc Ash3d_output.nc output_product_ID format_code [time_step]
!  where output_product_ID is one of:
!       0 custom (e.g. variable from user-defined module)
!       1 full concentration array
!       2 deposit granularity
!       3 deposit thickness (mm time-series)
!       4 deposit thickness (inches time-series)
!       5 deposit thickness (mm final)
!       6 deposit thickness (inches final)
!       7 ashfall arrival time (hours)
!       8 ashfall arrival at airports/POI (mm)
!       9 ash-cloud concentration (mg/m3)
!      10 ash-cloud height (km)
!      11 ash-cloud bottom (km)
!      12 ash-cloud load (T/km2 or )
!      13 ash-cloud radar reflectivity (dBz)
!      14 ash-cloud arrival time (hours)
!      15 topography
!      16 profile plots
!  and format_code is one of:
!       1 ASCII/ArcGIS
!       2 KML/KMZ
!       3 image/png
!       4 binary
!       5 shape file
!       6 grib2   (not yet implemented)
!       7 netcdf  (not yet implemented)
!       8 tecplot (not yet implemented)
!       9 vtk     (not yet implemented)
!  and time_step can optionally be provided for transient variables with -1 denoting
!  the final time data.  For data products that require contours (shape files or
!  contour plots), contour levels are defined in the module Output_Vars.
!  
!  For example, to produce a contour map of ash-cloud arrival time:
!   ./Ash3d_PostProc Ash3d_output.nc 14 3
!  
!  To produce a shapefile for deposit thickness (in mm) at the second time step:
!   ./Ash3d_PostProc Ash3d_output.nc 3 5 2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      program Ash3d_PostProc

      use precis_param

      use io_units

      use global_param,  only : &
         MM_2_IN,EPS_SMALL

      use io_data,       only : &
         iTimeNext,PP_infile,datafileIn,HaveInfile, &
         concenfile,nWriteTimes,WriteTimes,Write_PT_Data,Write_PR_Data,&
         iout3d,isFinal_TS,WriteGSD,WriteDepositTS_KML,WriteDepositTS_ASCII,&
         WriteDepositTime_KML,WriteDepositTime_ASCII,WriteDepositFinal_KML,&
         WriteDepositFinal_ASCII,WriteCloudTime_KML,WriteCloudTime_ASCII,WriteReflectivity_ASCII,&
         WriteCloudLoad_KML,WriteReflectivity_KML,WriteCloudLoad_ASCII,WriteCloudHeight_KML,&
         WriteCloudHeight_ASCII,WriteCloudConcentration_KML,WriteCloudConcentration_ASCII,&
         WriteAirportFile_KML,WriteAirportFile_ASCII,Write3dFiles,&
         nvprofiles,nvar_User2d_static_XY,nvar_User2d_XY

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,IsLatLon,dx,dy,xLL,yLL,ts0,ts1, &
           Allocate_mesh

      use solution,      only : &
         concen_pd,DepositGranularity

      use time_data,     only : &
         ntmax,time,time_native,BaseYear,useLeap,SimStartHour

      use Ash3d_Program_Control, only : &
           Set_OS_Env,                &
           Read_Control_File,         &
           Read_PostProc_Control_File

      use Output_Vars,   only : &
         DepositThickness,DepArrivalTime,CloudArrivalTime,ashcon_tot,&
         MaxConcentration,MaxHeight,CloudLoad,dbZCol,MinHeight,Mask_Cloud, &
         iplotpref,Extra2dVar,Extra2dVarName, &
           Gen_Output_Vars,  &
           Allocate_Output_Vars, &
           Set_OutVar_ContourLevel

      use help,          only : &
           help_postproc

      use Ash3d_ASCII_IO,  only : &
         A_nx,A_ny,A_XY,A_XYZ,A_xll,A_yll,A_dx,A_dy, &
           deallocate_ASCII, &
           write_2D_ASCII,   &
           read_2D_ASCII,    &
           write_3D_ASCII,   &
           read_3D_ASCII,    &
           vprofileopener,   &
           vprofilewriter,   &
           vprofilecloser,   &
           Write_PointData_Airports_ASCII

      use Ash3d_Binary_IO, only : &
         B_XY,B_XYZ,          &
           deallocate_Binary, &
           write_2D_Binary,   &
           read_2D_Binary,    &
           write_3D_Binary,   &
           read_3D_Binary

#ifdef USENETCDF
      use Ash3d_Netcdf_IO
#endif

      use Ash3d_KML_IO

#ifdef USEDISLIN
      use Ash3d_PostProc_dislin
#endif
#ifdef USEPLPLOT
      use Ash3d_PostProc_plplot
#endif
      use Ash3d_PostProc_gnuplot

      use Ash3d_PostProc_GMT

      implicit none

      integer             :: nargs
      integer             :: istat
      integer             :: iostatus
      integer             :: arglen
      character(len=120)  :: iomessage
      character(len=50)   :: linebuffer050
      character(len=80)   :: linebuffer080
      character(len=100)  :: arg
      integer             :: informat     = 3 ! input format (default = 3 for netcdf)
      integer             :: outformat        ! output format
      integer             :: iiprod           ! code for input dataset
      integer             :: iprod            ! code for output product
      integer             :: ndims           ! dimensions of the input data file
      integer             :: ivar
      integer             :: TS_Flag
      integer             :: height_flag
      integer             :: itime = -1      ! initialize time step to the last step
      integer             :: i,j,ii
      integer             :: tmp_int
      integer             :: icase
      real(kind=ip),dimension(:,:),allocatable :: OutVar
      logical      ,dimension(:,:),allocatable :: mask
      real(kind=ip),dimension(:,:),allocatable :: Topography
      real(kind=ip)       :: OutFillValue
      logical             :: IsThere
      character(len=6)    :: Fill_Value
      character(len=20)   :: filename_root
      character(len=70)   :: comd
      character(len=13)   :: cio
      character           :: testkey,testkey2
      logical             :: writeContours = .false.
      logical             :: CleanScripts  = .false.
        ! plotting library availibility and preferences
        !  1 = dislin
        !  2 = plplot
        !  3 = gnuplot
        !  4 = GMT
      integer, parameter  :: Nplot_libs = 4
      logical,dimension(Nplot_libs) :: plotlib_avail
                                                     !   -- First preference code 
                                                     !   | - Second
                                                     !   | | - Third
                                                     !   | | | - Fourth
                                                     !   V V V V
#ifdef WINDOWS
      ! For Windows systems, dislin is working; others not yet.
      integer,dimension(Nplot_libs) :: plot_pref_map = (/1,2,3,4/) ! plot preference for maps
      integer,dimension(Nplot_libs) :: plot_pref_shp = (/1,2,3,4/) ! plot preference for contours
      integer,dimension(Nplot_libs) :: plot_pref_vpr = (/1,2,3,4/) ! plot preference for vert profs.
      integer,dimension(Nplot_libs) :: plot_pref_aTS = (/1,2,3,4/) ! plot preference for Airport TS
#else
      integer,dimension(Nplot_libs) :: plot_pref_map = (/2,1,3,4/) ! plot preference for maps
      integer,dimension(Nplot_libs) :: plot_pref_shp = (/3,1,2,4/) ! plot preference for contours
      integer,dimension(Nplot_libs) :: plot_pref_vpr = (/1,2,3,4/) ! plot preference for vert profs.
      integer,dimension(Nplot_libs) :: plot_pref_aTS = (/2,3,1,4/) ! plot preference for Airport TS
#endif

      INTERFACE
        subroutine alloc_arrays
        end subroutine alloc_arrays
        subroutine calc_mesh_params
        end subroutine calc_mesh_params
        subroutine output_results
        end subroutine output_results
        subroutine write_ShapeFile_Polyline(iprod,itime)
          integer,intent(in) :: iprod
          integer,intent(in) :: itime
        end subroutine write_ShapeFile_Polyline
        character (len=13) function HS_yyyymmddhh_since(HoursSince,byear,useLeaps)
          real(kind=8)               ::  HoursSince
          integer                    ::  byear
          logical                    ::  useLeaps
        end function HS_yyyymmddhh_since
        character (len=20) function HS_xmltime(HoursSince,byear,useLeaps)
          real(kind=8)              :: HoursSince
          integer                   :: byear
          logical                   :: useLeaps
        end function HS_xmltime
        subroutine dealloc_arrays
        end subroutine dealloc_arrays
      END INTERFACE

      open(unit=fid_logfile,file='Ash3d_pp.log',status='replace',action='write')
      !   We want to call this subroutine silently, so reset the verbosity
      tmp_int = VB(1)
      VB(1)   = verbosity_dark
      VB(2)   = verbosity_dark
      call Set_OS_Env
      if(VB(1).eq.verbosity_dark)then
        VB(1)   = tmp_int
      endif

      if(VB(1).le.verbosity_debug1)then
        CleanScripts = .false.
      else
        CleanScripts = .true.
      endif
      CleanScripts_gnuplot = CleanScripts
      CleanScripts_GMT     = CleanScripts

      ! Checking to see which plotting packages we have
#ifdef USEDISLIN
      plotlib_avail(1) = .true.
#else
      plotlib_avail(1) = .false.
#endif
#ifdef USEPLPLOT
      plotlib_avail(2) = .true.
#else
      plotlib_avail(2) = .false.
#endif
      ! Test for gnuplot
#ifdef LINUX
        ! On a linux system, just try to execute gnuplot
      istat = 0
      call execute_command_line("echo 'exit' | gnuplot",exitstat=istat)
#endif
#ifdef MACOS
      ! On a MacOS system, not sure how to test yet
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Cannot test for gnuplot on MacOS for now."
        write(outlog(io),*)"Disabling gnuplot."
      endif;enddo
      istat = 1
#endif
#ifdef WINDOWS
      ! On a Windows system, not sure how to test yet
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Cannot test for gnuplot on Windows for now."
        write(outlog(io),*)"Disabling gnuplot."
      endif;enddo
      istat = 1
#endif
      if (istat.eq.0)then
        CleanScripts_gnuplot = CleanScripts
        plotlib_avail(3) = .true.
      else
        plotlib_avail(3) = .false.
      endif
      ! Test for GMT
#ifdef LINUX
        ! On a linux system, just try to execute gmt
      istat = 0
      call execute_command_line("gmt --version > /dev/null",exitstat=istat)
#endif
#ifdef MACOS
        ! On a MacOS system, not sure how to test yet
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Cannot test for gmt on MacOS for now."
          write(outlog(io),*)"Disabling gmt."
        endif;enddo
        istat = 1
#endif
#ifdef WINDOWS
        ! On a Windows system, not sure how to test yet
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Cannot test for gmt on Windows for now."
          write(outlog(io),*)"Disabling gmt."
        endif;enddo
        istat = 1
#endif
      if (istat.eq.0)then
        CleanScripts_GMT = CleanScripts
        plotlib_avail(4) = .true.
      else
        plotlib_avail(4) = .false.
      endif

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Dislin  ",plotlib_avail(1)
        write(outlog(io),*)"Plplot  ",plotlib_avail(2)
        write(outlog(io),*)"Gnuplot ",plotlib_avail(3)
        write(outlog(io),*)"GMT     ",plotlib_avail(4)
      endif;enddo

      ! Initialize all output logicals to false
      Write_PR_Data                 = .false.
      Write_PT_Data                 = .false.
      Write3dFiles                  = .false.
      WriteGSD                      = .false.
      WriteDepositTS_KML            = .false.
      WriteDepositTS_ASCII          = .false.
      WriteDepositTime_KML          = .false.
      WriteDepositTime_ASCII        = .false.
      WriteDepositFinal_KML         = .false.
      WriteDepositFinal_ASCII       = .false.
      WriteCloudTime_KML            = .false.
      WriteCloudTime_ASCII          = .false.
      WriteCloudLoad_KML            = .false.
      WriteReflectivity_ASCII       = .false.
      WriteReflectivity_KML         = .false.
      WriteCloudLoad_ASCII          = .false.
      WriteCloudHeight_KML          = .false.
      WriteCloudHeight_ASCII        = .false.
      WriteCloudConcentration_KML   = .false.
      WriteCloudConcentration_ASCII = .false.
      WriteAirportFile_KML          = .false.
      WriteAirportFile_ASCII        = .false.
      Write3dFiles                  = .false.

      ! Completed evaluating status of environment; now checking how to proceed:
      !  No command-line arguments (100) -> interactively prompt user
      !  One command-line argument (110) -> check for '-h' otherwise assume control file
      !  3+ command-line arguemnts (120) -> command-line only driven run

      ! Test read command-line arguments
      nargs = command_argument_count()
!100
      if (nargs.eq.0) then
          ! If no command-line arguments are given, then prompt user
          ! interactively for the command file name and possible a 
          ! restart file
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)'No command-line arguments detected'
          write(outlog(io),*)''
        endif;enddo

        if(VB(1).ge.verbosity_silent)then
          do io=1,2
            write(errlog(io),*)"Stdout is suppressed via verbosity=9,10, but interactive input is expected."
            write(errlog(io),*)"Either recompile with verbosity<9 or provide the correct command-line arguments."
          enddo
          stop 1
        endif
#ifndef USENETCDF
        ! If we are here, then we expect to read the netcdf output file.  If netcdf
        ! not linked, give an error and exit
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)'Expecting to prompt for a netcdf file, but the netcdf'
          write(errlog(io),*)'library is not linked.  Please recompile, linking to'
          write(errlog(io),*)'netcdf or run Ash3d_PostProc with a control file and'
          write(errlog(io),*)'ASCII/binary data.'
        endif;enddo
        stop 1
#else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)'Enter name of netcdf output file:'
        endif;enddo
        read(input_unit,*,iostat=iostatus,iomsg=iomessage) concenfile
        linebuffer080 = "concenfile"
        linebuffer050 = "Reading concenfile from stdin"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
#endif
        inquire( file=concenfile, exist=IsThere )
        if(.not.IsThere)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: Cannot find input file"
            write(errlog(io),*)"     ",concenfile
          endif;enddo
          stop 1
        endif

        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)'Select output variable (not all may be available):'
          write(outlog(io),*)' 1 full concentration array             2 deposit granularity'
          write(outlog(io),*)' 3 deposit thickness (mm time-series)   4 deposit thickness (inches time-series)'
          write(outlog(io),*)' 5 deposit thickness (mm final)         6 deposit thickness (inches final)'
          write(outlog(io),*)' 7 ashfall arrival time (hours)         8 ashfall arrival at airports/POI (mm)'
          write(outlog(io),*)' 9 ash-cloud concentration (mg/m3)     10 ash-cloud height (km)'
          write(outlog(io),*)'11 ash-cloud bottom (km)               12 ash-cloud load (T/km2)'
          write(outlog(io),*)'13 ash-cloud radar reflectivity (dBz)  14 ash-cloud arrival time (hours)'
          write(outlog(io),*)'15 topography                          16 profile plots'
          write(outlog(io),*)'  or enter 0 to be prompted for a variable name'
          write(outlog(io),*)''
          write(outlog(io),*)'Enter code for output product:'
        endif;enddo
        read(input_unit,*,iostat=iostatus,iomsg=iomessage)iprod
        linebuffer080 = "iprod"
        linebuffer050 = "Reading iprod from stdin"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        if(iprod.eq.0)then
          ! Get name of user variable
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)'Enter name of 2d variable in netcdf file:'
          endif;enddo
          read(input_unit,*,iostat=iostatus,iomsg=iomessage)Extra2dVarName
          linebuffer080 = "Extra2dVarName"
          linebuffer050 = "Reading Extra2dVarName from stdin"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)

          ! Get static vs TS
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)'Is this variable static or time-series?'
            write(outlog(io),*)' 1 static'
            write(outlog(io),*)' 2 time-series'
            write(outlog(io),*)''
            write(outlog(io),*)'Enter 1 or 2:'
          endif;enddo
          read(input_unit,*,iostat=iostatus,iomsg=iomessage)tmp_int
          linebuffer080 = "TSorNo"
          linebuffer050 = "Reading TSorNo from stdin"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
          if(tmp_int.eq.1)then
            nvar_User2d_static_XY = 1
          elseif(tmp_int.eq.2)then
            nvar_User2d_XY = 1
          else
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: Cannot find input file"
              write(errlog(io),*)"     ",concenfile
            endif;enddo
            stop 1
          endif
        elseif(iprod.eq.15)then
          nvar_User2d_static_XY = 1
          Extra2dVarName = "Topography"
        endif

        ! Before we do anything, call routine to read the netcdf file, populate
        ! the dimensions so we can see what we are dealing with.
        ! This call reads the 2d output products for the specified time
        ! If itime=-1 (for the final step), then
        !   We want to call this subroutine silently, so reset the verbosity
        tmp_int = VB(1)
        VB(1)   = verbosity_silent
#ifdef USENETCDF
        call NC_Read_Output_Products(-1)
#endif  
        VB(1)   = tmp_int

        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)'Select output format'
          write(outlog(io),*)' 1 ASCII/ArcGIS           2 KML/KMZ'
          write(outlog(io),*)' 3 image/png              4 binary'
          write(outlog(io),*)' 5 shape file'
          write(outlog(io),*)''
          write(outlog(io),*)'Enter code for output format:'
        endif;enddo
        read(input_unit,*,iostat=iostatus,iomsg=iomessage)outformat
        linebuffer080 = "outformat"
        linebuffer050 = "Reading outformat from stdin"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)

        if(iprod.eq.5.or.iprod.eq.6)then
          ! For final deposit variables, set itime to -1
          itime = -1
        elseif(iprod.eq.1.or.&
               iprod.eq.2.or.&
               iprod.eq.3.or.&
               iprod.eq.4.or.&
               iprod.eq.9.or.&
               iprod.eq.10.or.&
               iprod.eq.11.or.&
               iprod.eq.12.or.&
               iprod.eq.13)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)'Select time step:'
            do i=1,nWriteTimes
              write(outlog(io),*)i,real(WriteTimes(i),kind=sp),&
                                 HS_xmltime(SimStartHour+WriteTimes(i),BaseYear,useLeap)
            enddo
            write(outlog(io),*)'Enter index for time step:'
          endif;enddo
          read(input_unit,*,iostat=iostatus,iomsg=iomessage)itime
          linebuffer080 = "itime"
          linebuffer050 = "Reading itime from stdin"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        else
          ! iprod = 7,8,14,15 are not time-series, but set itime to -1
          itime = -1
        endif
!110
      elseif (nargs.eq.1) then
          ! If an argument is given, first test for the '-h' indicating a help
          ! request.
        call get_command_argument(1, arg, length=arglen, status=iostatus)
        testkey  = arg(1:1)
        testkey2 = arg(2:2)
        if(testkey.eq.'-'.and.testkey2.eq.'h')then
          ! This is the branch for user-requested help
          ! command is Ash3d_PostProc -h
          call help_postproc
        else
          ! Read control file
          call get_command_argument(1, arg, length=arglen, status=iostatus)
          !read(arg,*,iostat=iostatus,iomsg=iomessage)PP_infile
          PP_infile = trim(adjustl(arg))
          linebuffer050 = "Reading control file from command-line arg."
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,arg(1:80),iomessage)
          inquire( file=PP_infile, exist=IsThere )
          if(.not.IsThere)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: Cannot find input control file"
            endif;enddo
            stop 1
          endif

          call Read_PostProc_Control_File(informat,iiprod,iprod,ndims,outformat,iplotpref,itime)
!          ! Reset plotting preference if need be
!          if(iplotpref.gt.0)then
!            if(plotlib_avail(iplotpref))then
!              plot_pref_map(1:Nplot_libs) = iplotpref ! plot preference for maps
!              plot_pref_shp(1:Nplot_libs) = iplotpref ! plot preference for contours
!              plot_pref_vpr(1:Nplot_libs) = iplotpref ! plot preference for vert profs.
!              plot_pref_aTS(1:Nplot_libs) = iplotpref ! plot preference for Airport TS
!            else
!              do io=1,2;if(VB(io).le.verbosity_error)then
!                write(errlog(io),*)"WARNING: Preferred plotting library is not available."
!              endif;enddo
!            endif
!          endif
        endif
!120
          elseif (nargs.ge.3) then
        ! If we are doing command line only, then we need at least the netcdf filename, the output
        ! product code and the format.  Optionally, we can add the timestep.  If there is an
        ! inconsistency with iprod and outformat, an error message is issued before stopping.
        call get_command_argument(1, arg, length=arglen, status=iostatus)
        ! Note: unformatted read of concenfile from arg will only read the bit
        ! up to the first '/'
        !read(arg,*,iostat=iostatus,iomsg=iomessage)concenfile
        concenfile = trim(adjustl(arg))
        !linebuffer050 = "Reading concenfile from command-line arg 1"
        !if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,arg(1:80),iomessage)
        call get_command_argument(2, arg, length=arglen, status=iostatus)
        read(arg,*,iostat=iostatus,iomsg=iomessage)iprod
        linebuffer050 = "Reading iprod from command-line arg 2"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,arg(1:80),iomessage)
        call get_command_argument(3, arg, length=arglen, status=iostatus)
        read(arg,*,iostat=iostatus,iomsg=iomessage)outformat
        linebuffer050 = "Reading outformat from command-line arg 3"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,arg(1:80),iomessage)
        if (nargs.eq.4) then
          call get_command_argument(4, arg, length=arglen, status=iostatus)
          read(arg,*,iostat=iostatus,iomsg=iomessage)itime
          linebuffer050 = "Reading itime from command-line arg 4"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,arg(1:80),iomessage)
        else
          itime = -1
        endif
      else
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)' Cannot parse command line'
        endif;enddo
        stop 1
      endif

      ! Now that we have read the command line, error-check and report back what
      ! we are about to do.
      !  Arg #1
!130
      if(informat.eq.3)then
        ! Test if the Ash3d netcdf file exists
        inquire( file=concenfile, exist=IsThere )
        if(.not.IsThere)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: Cannot find input file"
            write(errlog(io),*)"     ",concenfile
          endif;enddo
          stop 1
        else
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Ash3d Output file: ",concenfile
          endif;enddo
        endif
      else
        ! input file is either ASCII or Binary
        ! Test if it can be found
        inquire( file=datafileIn, exist=IsThere )
        if(.not.IsThere)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: Cannot find input file"
            write(errlog(io),*)"     ",datafileIn
          endif;enddo
          stop 1
        else
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Ash3d Output file: ",datafileIn
          endif;enddo
        endif
      endif

      !  Arg #2
      if(iprod.lt.0.or.iprod.gt.16)then    ! Activate this line when custom variables are coded
      !if(iprod.lt.1.or.iprod.gt.16)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: output product requested is not in range 0-16."
          write(errlog(io),*)"       Run Ash3d_PostProc with no command-line"
          write(errlog(io),*)"       arguments to set usage information"
        endif;enddo
        stop 1
      else
        if(iprod.eq.1)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)'output variable = 1 full concentration array'
            write(outlog(io),*)' Currently, no output formats available for full'
            write(outlog(io),*)' granularity.  Binary output will give total ash'
            write(outlog(io),*)' concentration.'
          endif;enddo
        elseif(iprod.eq.2)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'output variable = 2 deposit granularity'
            write(errlog(io),*)' Currently, no output formats available for iprod=2'
            write(errlog(io),*)' '
          endif;enddo
          stop 1
        elseif(iprod.eq.3)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)'output variable = 3 deposit thickness; TS or step (mm)'
          endif;enddo
          ivar = 7  ! Note that kml writes out all netcdf time steps followed by the final
          TS_flag = 1      ! 1 = time series
          height_flag = 0  ! All the cells should be pinned to z=0
        elseif(iprod.eq.4)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)'output variable = 4 deposit thickness: TS or step (inches)'
          endif;enddo
          ivar = 8  ! Note that kml writes out all netcdf time steps followed by the final
          TS_flag = 1      ! 1 = time series
          height_flag = 0  ! All the cells should be pinned to z=0
        elseif(iprod.eq.5)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)'output variable = 5 deposit thickness: final (mm)'
          endif;enddo
          ivar = 7  ! Note that kml writes out all netcdf time steps followed by the final
          TS_flag = 1      ! 1 = time series (really, this final value is not a time series, but the kml writer requires this)
          height_flag = 0  ! All the cells should be pinned to z=0
        elseif(iprod.eq.6)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)'output variable = 6 deposit thickness: final (inches)'
          endif;enddo
          ivar = 8  ! Note that kml writes out all netcdf time steps followed by the final
          TS_flag = 1      ! 1 = time series (really, this final value is not a time series, but the kml writer requires this)
          height_flag = 0  ! All the cells should be pinned to z=0
        elseif(iprod.eq.7)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)'output variable = 7 ashfall arrival time (hours)'
          endif;enddo
          ivar = 9
          TS_flag = 0      ! 1 = not a time series
          height_flag = 0  ! All the cells should be pinned to z=0
        elseif(iprod.eq.8)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)'output variable = 8 ashfall arrival at airports/POI (mm)'
          endif;enddo
          !ivar = NaN There is a special KML writer for this variable
        elseif(iprod.eq.9)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)'output variable = 9 ash-cloud concentration (mg/m3)'
          endif;enddo
          ivar = 1
          TS_flag = 1      ! 1 = time series
          height_flag = 1  ! All the cells should be at cloud height
        elseif(iprod.eq.10)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)'output variable =10 ash-cloud height (km)'
          endif;enddo
          ivar = 2
          TS_flag = 1      ! 1 = time series
          height_flag = 1  ! All the cells should be at cloud height
        elseif(iprod.eq.11)then
          do io=1,2;if(VB(io).le.verbosity_info)then          
            write(outlog(io),*)'output variable =11 ash-cloud bottom (km)'
          endif;enddo
          ivar = 3
          TS_flag = 1      ! 1 = time series
          height_flag = 1  ! All the cells should be at cloud height
        elseif(iprod.eq.12)then
          do io=1,2;if(VB(io).le.verbosity_info)then          
            write(outlog(io),*)'output variable =12 ash-cloud load (T/km2)'
          endif;enddo
          ivar = 4
          TS_flag = 1      ! 1 = time series
          height_flag = 1  ! All the cells should be at cloud height
        elseif(iprod.eq.13)then
          do io=1,2;if(VB(io).le.verbosity_info)then          
            write(outlog(io),*)'output variable =13 ash-cloud radar reflectivity (dBz)'
          endif;enddo
          ivar = 6
          TS_flag = 1      ! 1 = time series
          height_flag = 1  ! All the cells should be at cloud height
        elseif(iprod.eq.14)then
          do io=1,2;if(VB(io).le.verbosity_info)then          
            write(outlog(io),*)'output variable =14 ash-cloud arrival time (hours)'
          endif;enddo
          ivar = 5
          TS_flag = 0      ! 1 = not a time series
          height_flag = 0  ! All the cells should be pinned to z=0
        elseif(iprod.eq.15)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'output variable =15 Topography (km)'
          endif;enddo
          nvar_User2d_static_XY = 1
          Extra2dVarName = "Topography"
          ivar = 10
          TS_flag = 0      ! 1 = not a time series
          height_flag = 0  ! All the cells should be pinned to z=0
        elseif(iprod.eq.16)then
          do io=1,2;if(VB(io).le.verbosity_info)then          
            write(outlog(io),*)'output variable =16 vertical concentration profile'
          endif;enddo
        endif
      endif
      !  Arg #3
      if(outformat.lt.1.or.outformat.gt.7)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: output format requested is not in range 1-7."
          write(errlog(io),*)"       Run Ash3d_PostProc with no command-line"
          write(errlog(io),*)"       arguments to set usage information"
        endif;enddo
        stop 1
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          if(outformat.eq.1)then
            write(outlog(io),*)'output format = 1 ASCII/ArcGIS'
          elseif(outformat.eq.2)then
            write(outlog(io),*)'output format = 2 KMZ'
          elseif(outformat.eq.3)then
            write(outlog(io),*)'output format = 3 image/png'
          elseif(outformat.eq.4)then
            write(outlog(io),*)'output format = 4 binary'
          elseif(outformat.eq.5)then
            write(outlog(io),*)'output format = 5 shape file'
          elseif(outformat.eq.6)then
            write(outlog(io),*)'output format = 6 grib2'
          elseif(outformat.eq.7)then
            write(outlog(io),*)'output format = 7 tecplot'
          endif
        endif;enddo
      endif
      !  Arg #4
      if(itime.eq.-1)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)'itime = -1 (Final time step)'
          write(outlog(io),*)'  This signifiies one of two conditions:'
          write(outlog(io),*)'   1. No time step provided (e.g. variable is not at time-series)'
          write(outlog(io),*)'   2. The final time step should be used'
        endif;enddo
      elseif(itime.eq.0)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: itime = 0.  Invalid time step."
        endif;enddo
        stop 1
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)'itime = ',itime
          write(outlog(io),*)'  We do not yet know the maximum number of steps available.'
        endif;enddo
      endif
      if(iprod.eq.2)then   ! deposit granularity
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)'We do not yet have processing depocon implimented.'
        endif;enddo
        stop 1
      endif
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)'Finished reading inputs.'
        write(outlog(io),*)' '
      endif;enddo

      ! Reset plotting preference if need be
      if(iplotpref.gt.0)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"  Plotting library reset"
          if(iplotpref.eq.1)then
            write(outlog(io),*)"           Using dislin if available."
          elseif(iplotpref.eq.2)then
            write(outlog(io),*)"           Using plplot if available."
          elseif(iplotpref.eq.3)then
            write(outlog(io),*)"           Using gnuplot if available."
          elseif(iplotpref.eq.4)then
            write(outlog(io),*)"           Using GMT if available."
          endif
        endif;enddo
        if(plotlib_avail(iplotpref))then
          plot_pref_map(1:Nplot_libs) = iplotpref ! plot preference for maps
          plot_pref_shp(1:Nplot_libs) = iplotpref ! plot preference for contours
          plot_pref_vpr(1:Nplot_libs) = iplotpref ! plot preference for vert profs.
          plot_pref_aTS(1:Nplot_libs) = iplotpref ! plot preference for Airport TS
        else
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"WARNING: Preferred plotting library is not available."
          endif;enddo
        endif
      endif

      ! Before we do anything, we need to set up as much as we can of the grid
      ! and populate auxilary variable.
      if(informat.eq.3)then ! netcdf
        ! call routine to read the netcdf file, populate
        ! the dimensions so we can see what we are dealing with.
        ! This call reads the 2d output products for the specified time
        ! If itime=-1 (for the final step), then
#ifndef USENETCDF
        ! If we are here, then we expect to read the netcdf output file.  If netcdf
        ! not linked, give an error and exit
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)'A netcdf file was provided in the control file, but'
          write(errlog(io),*)'the netcdf library is not linked.  Please recompile,'
          write(errlog(io),*)'linking to netcdf or run Ash3d_PostProc with a control'
          write(errlog(io),*)'file and ASCII data.'
        endif;enddo
        stop 1
#else
        call NC_Read_Output_Products(itime)
#endif
      else
        if(HaveInfile)then
          ! If we have the Ash3d control file, read it to set grid and populate
          ! volcano name, etc.
          call Read_Control_File
        endif
        if(informat.eq.1)then ! ASCII
          if(ndims.eq.2)then
            call read_2D_ASCII(datafileIn)
            if(nxmax.ne.A_nx.or.  &
               nymax.ne.A_ny)then
              do io=1,2;if(VB(io).le.verbosity_error)then
                write(errlog(io),*)"WARNING: nx,ny in ASCII file are not the same as in control file."
                write(errlog(io),*)"         Resetting nx,ny to ASCII values."
                write(errlog(io),*)"    nx = ",nxmax
                write(errlog(io),*)"    ny = ",nymax
                write(errlog(io),*)"  A_nx = ",A_nx
                write(errlog(io),*)"  A_ny = ",A_ny
              endif;enddo
            endif
            if(abs(xLL-A_xll).gt.EPS_SMALL.or.  &
               abs(yLL-A_yll).gt.EPS_SMALL)then
              do io=1,2;if(VB(io).le.verbosity_error)then
                write(errlog(io),*)"WARNING: xLL,yLL in ASCII file are not the same as in control file."
                write(errlog(io),*)"         Resetting xLL,yLL to ASCII values."
                write(errlog(io),*)"    xLL = ",xLL
                write(errlog(io),*)"    yLL = ",yLL
                write(errlog(io),*)"  A_xLL = ",A_xLL
                write(errlog(io),*)"  A_yLL = ",A_yLL
              endif;enddo
            endif
            if(abs(dx-A_dx).gt.EPS_SMALL.or.  &
               abs(dy-A_dy).gt.EPS_SMALL)then
              do io=1,2;if(VB(io).le.verbosity_error)then
                write(errlog(io),*)"WARNING: dx,dy in ASCII file are not the same as in control file."
                write(errlog(io),*)"         Resetting dx,dy to ASCII values."
                write(errlog(io),*)"    dx = ",dx
                write(errlog(io),*)"    dy = ",dy
                write(errlog(io),*)"  A_dx = ",A_dx
                write(errlog(io),*)"  A_dy = ",A_dy
              endif;enddo
            endif
          elseif(ndims.eq.3)then
            call read_3D_ASCII(datafileIn)
          endif
        elseif(informat.eq.2)then ! BINARY
          if(ndims.eq.2)then
            ! We didn't error-check nxmax and nymax on input, so do it now
            if(nxmax.le.0.or.nymax.le.0)then
              do io=1,2;if(VB(io).le.verbosity_error)then
                write(errlog(io),*)"ERROR: Either nx or ny is not a positive integer."
                stop 1
              endif;enddo
            endif
            call read_2D_Binary(nxmax,nymax,datafileIn)
          elseif(ndims.eq.3)then
            call read_3D_BINARY(nxmax,nymax,nzmax,datafileIn)
          endif
        endif ! informat = 2 or 3

        ! Allocate grid
        call Allocate_mesh
        call calc_mesh_params
        ! Load contour levels
        call Set_OutVar_ContourLevel
      endif  ! informat = 1 or not

      if(.not.IsLatLon)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(outlog(io),*)'Mapping/shapefiles of projected grids is currently only supported'
          write(outlog(io),*)'using the GMT plotting option.'
        endif;enddo
        if(plotlib_avail(4))then
          plot_pref_map = (/4,1,2,3/)
        else
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(outlog(io),*)'Mapping/shapefiles of projected grids requested, but GMT is not available.'
          endif;enddo
          stop 1
        endif
      endif

      if (itime.eq.-1) then
        iout3d = nWriteTimes
        iTimeNext = nWriteTimes
        isFinal_TS = .true.
      else
        iout3d = itime
        iTimeNext = 0
        isFinal_TS = .false.
      endif

      cio = HS_yyyymmddhh_since(SimStartHour+time,BaseYear,useLeap)

      if(    iprod.eq.1 )then ! full concentration array
        Write3dFiles = .true.
      elseif(iprod.eq.2 )then ! deposit granularity
        if(    outformat.eq.1)then
          WriteGSD                      = .true.
        elseif(outformat.eq.2)then
          WriteGSD                      = .true.
        endif
      elseif(iprod.eq.3 )then ! deposit at specified times (mm)
        if(    outformat.eq.1)then
          WriteDepositTS_ASCII          = .true.
        elseif(outformat.eq.2)then
          WriteDepositTS_KML            = .true.
        endif
      elseif(iprod.eq.4 )then ! deposit at final time (mm)
        if(    outformat.eq.1)then
          WriteDepositFinal_ASCII       = .true.
        elseif(outformat.eq.2)then
          WriteDepositFinal_KML         = .true.
        endif
      elseif(iprod.eq.5 )then ! deposit at specified times (inches)
        if(    outformat.eq.1)then
          WriteDepositTS_ASCII          = .true.
        elseif(outformat.eq.2)then
          WriteDepositTS_KML            = .true.
        endif
      elseif(iprod.eq.6 )then ! deposit at final time (inches)
        if(    outformat.eq.1)then
          WriteDepositFinal_ASCII       = .true.
        elseif(outformat.eq.2)then
          WriteDepositFinal_KML         = .true.
        endif
      elseif(iprod.eq.7 )then ! ashfall arrival time
        if(    outformat.eq.1)then
          WriteDepositTime_ASCII        = .true.
        elseif(outformat.eq.2)then
          WriteDepositTime_KML          = .true.
        endif
      elseif(iprod.eq.8 )then ! ashfall at airports/POI
        Write_PR_Data                   = .true.
        if(    outformat.eq.1)then
          WriteAirportFile_ASCII        = .true.
        elseif(outformat.eq.2)then
          WriteAirportFile_KML          = .true.
        endif
      elseif(iprod.eq.9 )then ! ash-cloud concentration
        if(    outformat.eq.1)then
          WriteCloudConcentration_ASCII = .true.
        elseif(outformat.eq.2)then
          WriteCloudConcentration_KML   = .true.
        endif
      elseif(iprod.eq.10)then ! ash-cloud height
        if(    outformat.eq.1)then
          WriteCloudHeight_ASCII        = .true.
        elseif(outformat.eq.2)then
          WriteCloudHeight_KML          = .true.
        endif
      elseif(iprod.eq.11)then ! ash-cloud bottom
        if(    outformat.eq.1)then
          WriteCloudHeight_ASCII        = .true.
        elseif(outformat.eq.2)then
          WriteCloudHeight_KML          = .true.
        endif
      elseif(iprod.eq.12)then ! ash-cloud load
        if(    outformat.eq.1)then
          WriteCloudLoad_ASCII          = .true.
        elseif(outformat.eq.2)then
          WriteCloudLoad_KML            = .true.
        endif
      elseif(iprod.eq.13)then ! ash-cloud radar reflectivity
        if(    outformat.eq.1)then
          WriteReflectivity_ASCII       = .true.
        elseif(outformat.eq.2)then
          WriteReflectivity_KML         = .true.
        endif
      elseif(iprod.eq.14)then ! ash-cloud arrival time
        if(    outformat.eq.1)then
          WriteCloudTime_ASCII          = .true.
        elseif(outformat.eq.2)then
          WriteCloudTime_KML            = .true.
        endif
      elseif(iprod.eq.15)then ! topography
        nvar_User2d_static_XY = 1
        Extra2dVarName = "Topography"
      elseif(iprod.eq.16)then ! vertical profiles of concentration
        if(    outformat.eq.1.or.outformat.eq.3)then
          ! ASCII or png
          Write_PR_Data                 = .true.
        else
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"Vertical ash concentration profiles can only be exported",&
                      " as ASCII files or png images"
          endif;enddo
          stop 1
        endif
      endif

      allocate(OutVar(nxmax,nymax))
      allocate(mask(nxmax,nymax))
      mask = .true.
      ! Load the variable OutVar from the ASCII or Binary arrays
      if(informat.eq.1)then
        if(ndims.eq.2)then
          OutVar(1:nxmax,1:nymax) = A_XY(1:nxmax,1:nymax)
        elseif(ndims.eq.3)then
          allocate(concen_pd(-1:nxmax+2,-1:nymax+2,-1:nzmax+2,1:nsmax,ts0:ts1)); concen_pd = 0.0_ip
          concen_pd(1:nxmax,1:nymax,1:nzmax,1,ts1) = A_XYZ(1:nxmax,1:nymax,1:nzmax)
          allocate(DepositGranularity(nxmax,nymax,nsmax)); DepositGranularity = 0.0_ip
          call Allocate_Output_Vars
          call Gen_Output_Vars
        endif
      elseif(informat.eq.2)then
        if(ndims.eq.2)then
          OutVar(1:nxmax,1:nymax) = B_XY(1:nxmax,1:nymax)
        elseif(ndims.eq.3)then
          allocate(concen_pd(-1:nxmax+2,-1:nymax+2,-1:nzmax+2,1:nsmax,ts0:ts1)); concen_pd = 0.0_ip
          concen_pd(1:nxmax,1:nymax,1:nzmax,1,ts1) = B_XYZ(1:nxmax,1:nymax,1:nzmax)
          allocate(DepositGranularity(nxmax,nymax,nsmax)); DepositGranularity = 0.0_ip
          call Allocate_Output_Vars
          call Gen_Output_Vars
        endif
      endif

      ! Now depending on the output product ID, copy OutVar to the named array
      if(informat.eq.1.or.informat.eq.2)then
        if(ndims.eq.2)then
          ! If we have read a 2d array, copy it to the proper named array so it can be
          ! processed correctly
          if(iprod.eq.3.or.iprod.eq.4.or.iprod.eq.5.or.iprod.eq.6)then
#ifdef USEPOINTERS
            if(.not.associated(DepositThickness))then
#else
            if(.not.allocated(DepositThickness))then
#endif
              allocate(DepositThickness(nxmax,nymax))
            endif
            DepositThickness(1:nxmax,1:nymax) = OutVar(1:nxmax,1:nymax)
          endif

          if(iprod.eq. 7)then
#ifdef USEPOINTERS
            if(.not.associated(DepArrivalTime))then
#else
            if(.not.allocated(DepArrivalTime))then
#endif
              allocate(DepArrivalTime(nxmax,nymax))
            endif
            DepArrivalTime(1:nxmax,1:nymax)   = OutVar(1:nxmax,1:nymax)
          endif
          if(iprod.eq. 8)then
            stop 1
          endif
          if(iprod.eq. 9)then
#ifdef USEPOINTERS
            if(.not.associated(MaxConcentration))then
#else
            if(.not.allocated(MaxConcentration))then
#endif
              allocate(MaxConcentration(nxmax,nymax))
            endif
            MaxConcentration(1:nxmax,1:nymax) = OutVar(1:nxmax,1:nymax)
          endif
          if(iprod.eq.10)then
#ifdef USEPOINTERS
            if(.not.associated(MaxHeight))then
#else
            if(.not.allocated(MaxHeight))then
#endif
              allocate(MaxHeight(nxmax,nymax))
            endif
            MaxHeight(1:nxmax,1:nymax)        = OutVar(1:nxmax,1:nymax)
          endif
          if(iprod.eq.11)then
#ifdef USEPOINTERS
            if(.not.associated(MinHeight))then
#else
            if(.not.allocated(MinHeight))then
#endif
              allocate(MinHeight(nxmax,nymax))
            endif
            MinHeight(1:nxmax,1:nymax)        = OutVar(1:nxmax,1:nymax)
          endif
          if(iprod.eq.12)then
#ifdef USEPOINTERS
            if(.not.associated(CloudLoad))then
#else
            if(.not.allocated(CloudLoad))then
#endif
              allocate(CloudLoad(nxmax,nymax))
            endif
            CloudLoad(1:nxmax,1:nymax)        = OutVar(1:nxmax,1:nymax)
          endif
          if(iprod.eq.13)then
#ifdef USEPOINTERS
            if(.not.associated(dbZCol))then
#else
            if(.not.allocated(dbZCol))then
#endif
              allocate(dbZCol(nxmax,nymax))
            endif
            dbZCol(1:nxmax,1:nymax)           = OutVar(1:nxmax,1:nymax)
          endif
          if(iprod.eq.14)then
#ifdef USEPOINTERS
            if(.not.associated(CloudArrivalTime))then
#else
            if(.not.allocated(CloudArrivalTime))then
#endif
              allocate(CloudArrivalTime(nxmax,nymax))
            endif
            CloudArrivalTime(1:nxmax,1:nymax) = OutVar(1:nxmax,1:nymax)
          endif
          if(iprod.eq.15)then
            write(*,*)"Allocating topo"
            if(.not.allocated(Topography)) allocate(Topography(nxmax,nymax))
            Topography(1:nxmax,1:nymax)       = OutVar(1:nxmax,1:nymax)
          endif
        elseif(ndims.eq.3)then
          ! For 3d input data, we only have the total 3d ash concentration, so we need
          ! to stop if the requested output product cannot be calculated from this.
          if(iprod.eq.3.or.iprod.eq.4.or.iprod.eq.5.or.iprod.eq.6.or.  & ! any of the deposit products
             iprod.eq.7.or.   & ! ashfall arrival time
             iprod.eq.8.or.   & ! ashfall arrival at airports
             iprod.eq.13.or.  & ! ash-cloud radar reflectivity
             iprod.eq.14.or.  & ! ash-cloud arrival time
             iprod.eq.16)then   ! profile plots
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: Requested output is not available from input data."
            endif;enddo
            stop 1
          endif
        endif
      endif !informat.eq.1.or.informat.eq.2

      ! The main differences in output products will be time-series, vs
      ! time-step output.  Currently, the time-series output will only be for
      ! the KML/KMZ files.
      if(outformat.eq.2)then
        ! KML output, separate into static vs time-series options
        !  First the static options
        if(iprod.eq.7.or.  &  ! Ashfall arrival time
           iprod.eq.14.or. &  ! Ash-cloud arrival time
           iprod.eq.15)then   ! Topography
          ! We really should have the final deposit output be in this catagory
          ! Static output
          ! We have already called the netcdf reader so we just need to write
          ! the KML output.
          ! Set up KML output files and prelim variables
          call Set_OutVar_Specs
          call OpenFile_KML(ivar)
          if(iprod.eq.7)then
            OutVar = real(DepArrivalTime(1:nxmax,1:nymax),kind=ip)
          elseif(iprod.eq.14)then
            OutVar = real(CloudArrivalTime(1:nxmax,1:nymax),kind=ip)
          elseif(iprod.eq.15)then
            OutVar = Topography(1:nxmax,1:nymax)
          endif
          call Write_2D_KML(ivar,OutVar,height_flag,TS_flag)
          call Close_KML(ivar,TS_flag)
          do io=1,2;if(VB(io).le.verbosity_info)then          
            write(outlog(io),*)"Zipping KML file."
          endif;enddo
          write(comd,*)"zip ",trim(adjustl(KMZ_filename(ivar))),' ',&
                              trim(adjustl(KML_filename(ivar)))
          call execute_command_line (comd, exitstat=istat)
        elseif(iprod.eq.3.or.iprod.eq.5.or. &  ! Deposit thickness mm (kml versions is TS + final)
               iprod.eq.4.or.iprod.eq.6.or. &  ! Deposit thickness inches (kml versions is TS + final)
               iprod.eq.9.or.               &  ! ash-cloud concentration
               iprod.eq.10.or.              &  ! ash-cloud height
               iprod.eq.11.or.              &  ! ash-cloud bottom
               iprod.eq.12.or.              &  ! ash-cloud load
               iprod.eq.13)then                ! radar reflectivity
          ! Time-series output
          ! Set up KML output files and prelim variables

          ! Currently deposit and other kml files include all of the TS so reset
          ! this flag
          isFinal_TS = .false.
          iTimeNext = 1
          time = WriteTimes(1)

          ! We need to loop over all times
          do i = 1,nWriteTimes
#ifdef USENETCDF
            call NC_Read_Output_Products(i)
#endif
            time = WriteTimes(i)
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"time = ",time
            endif;enddo
            call output_results
          enddo
          isFinal_TS = .true.
          ! For deposit output, load the final deposit variable and write to file
#ifdef USENETCDF
          call NC_Read_Output_Products(-1)
#endif
          call output_results
        elseif(iprod.eq.8)then  ! ashfall at airports
          call Write_PointData_Airports_KML
        elseif(iprod.eq.16)then  ! vertical profile plots
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: KML versions of vertical profiles not implemented."
          endif;enddo
          stop 1
        else
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: Requested iprod not implemented for KML."
          endif;enddo
          stop 1
        endif
      endif ! outformat.eq.2 (KML)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! This is the non-KML section
      ! Set some variable parameters for an ASCII output file
      OutFillValue = 0.0_ip
      if(iprod.eq.3.or.iprod.eq.5)then
        OutVar = DepositThickness
        Fill_Value = '-9999.'
        OutFillValue = 0.0_ip
        filename_root = 'DepositFile_        '
      elseif(iprod.eq.4.or.iprod.eq.6)then
        OutVar = DepositThickness*MM_2_IN
        Fill_Value = '-9999.'
        OutFillValue = 0.0_ip
        filename_root = 'DepositFile_        '
      elseif(iprod.eq.7)then
        OutVar = DepArrivalTime
        !OutVar = DepArrivalTime * merge(1.0_ip,0.0_ip,Mask_Deposit)
        Fill_Value = '-9999.'
        OutFillValue = -1.0_ip
        filename_root = 'DepositArrivalTime  '
      elseif(iprod.eq.8)then
         ! ashfall at airports/POI
         ! None of OutVar,Fill_Value,OutFillValue,filename_root need to be set
      elseif(iprod.eq.9)then
        OutVar = MaxConcentration
        !OutVar = MaxConcentration * merge(1.0_ip,0.0_ip,Mask_Cloud)
        Fill_Value = '-9999.'
        OutFillValue = 0.0_ip
        filename_root = 'CloudConcentration_ '
      elseif(iprod.eq.10)then
        OutVar = MaxHeight
        !OutVar = MaxHeight * merge(1.0_ip,0.0_ip,Mask_Cloud)
        Fill_Value = '-9999.'
        OutFillValue = 0.0_ip
        filename_root = 'CloudHeight_        '
      elseif(iprod.eq.11)then
        OutVar = MinHeight
        !OutVar = MinHeight * merge(1.0_ip,0.0_ip,Mask_Cloud)
        Fill_Value = '-9999.'
        OutFillValue = 0.0_ip
        filename_root = 'CloudHeightBot_     '
      elseif(iprod.eq.12)then
        OutVar = CloudLoad
        !OutVar = CloudLoad * merge(1.0_ip,0.0_ip,Mask_Cloud)
        Fill_Value = '-9999.'
        OutFillValue = 0.0_ip
        filename_root = 'CloudLoad_          '
      elseif(iprod.eq.13)then
        OutVar = dbZCol
        !OutVar = dbZCol * merge(1.0_ip,0.0_ip,Mask_Cloud)
        Fill_Value = '-9999.'
        OutFillValue = 0.0_ip
        filename_root = 'ClouddbZC_          '
      elseif(iprod.eq.14)then
        OutVar = real(CloudArrivalTime,kind=ip)
        !OutVar = CloudArrivalTime * merge(1.0_ip,0.0_ip,Mask_Cloud)
        Fill_Value = '-9999.'
        OutFillValue = -1.0_ip
        filename_root = 'CloudArrivalTime    '
      elseif(iprod.eq.15)then
        OutVar = real(Extra2dVar,kind=ip)
        Fill_Value = '-9999.'
        OutFillValue = -1.0_ip
        filename_root = 'Topography          '
      endif
      ! Now mask out non-cloud values
      if(iprod.eq.10.or.&  ! CloudHeight
         iprod.eq.11)then  ! CloudHeightBot
        mask = Mask_Cloud
      elseif(iprod.eq.14)then
        ! cloud mask based on cloud load does not work in this case the cloud load mask
        ! is a function of time
        mask(1:nxmax,1:nymax) = .true.
        do i=1,nxmax
          do j=1,nymax
            if(CloudArrivalTime(i,j).lt.0.0_ip)mask(i,j) = .false.
          enddo
        enddo
      endif

      if(outformat.eq.1)then  ! ASCII
        ! First check for the special cases
        if(iprod.eq.8)then
          ! Point data
          do io=1,2;if(VB(io).le.verbosity_info)then          
            write(outlog(io),*)"Calling Write_PointData_Airports_ASCII"
          endif;enddo
          call Write_PointData_Airports_ASCII
        elseif(iprod.eq.16)then
          ! Vertical profile data
          call vprofileopener
          do itime=1,ntmax
            time = time_native(itime)
            call vprofilewriter(itime)
          enddo
          call vprofilecloser
        else
          ! All other ESRI/ASCII 2d grids
          call write_2D_ASCII(nxmax,nymax,OutVar,mask,Fill_Value,filename_root)
        endif
      elseif(outformat.eq.2)then ! KML
        ! All the KML routines were called above
      elseif(outformat.eq.3)then ! image/png
        if(iprod.eq.8)then
          ! Point data
          do io=1,2;if(VB(io).le.verbosity_info)then          
            write(outlog(io),*)"No PNG output for point data output"
          endif;enddo
        elseif(iprod.eq.16)then
          ! Vertical profile data
          do i=1,nvprofiles
            icase = 0
            do ii=1,Nplot_libs
              ! Check each preference in series and see if the library is available
              if(plotlib_avail(plot_pref_vpr(ii)))then
                icase = plot_pref_vpr(ii)
                exit
              endif
            enddo

            select case (icase)
            case(1)
#ifdef USEDISLIN
              call write_2Dprof_PNG_dislin(i)
#endif
            case(2)
#ifdef USEPLPLOT
              call write_2Dprof_PNG_plplot(i)
#endif
            case(3)
              call write_2Dprof_PNG_gnuplot(i)
            case(4)
              call write_2Dprof_PNG_GMT(i)
            case default
              do io=1,2;if(VB(io).le.verbosity_error)then
                write(errlog(io),*)"ERROR: Plots requested but no plotting package is installed"
              endif;enddo
              stop 1
            end select
          enddo
        else
          ! If not point data (8) or vertical profiles (16), then do a map plot

        icase = 0
        do ii=1,Nplot_libs
          ! Check each preference in series and see if the library is available
          if(plotlib_avail(plot_pref_map(ii)))then
            icase = plot_pref_map(ii)
            exit
          endif
        enddo
        select case (icase)
        case(1)
#ifdef USEDISLIN
          call write_2Dmap_PNG_dislin(nxmax,nymax,iprod,iout3d,OutVar,writeContours)
#endif
        case(2)
#ifdef USEPLPLOT
          call write_2Dmap_PNG_plplot(nxmax,nymax,iprod,iout3d,OutVar,writeContours)
#endif
        case(3)
          call write_2Dmap_PNG_gnuplot(nxmax,nymax,iprod,iout3d,OutVar,writeContours)
        case(4)
          call write_2Dmap_PNG_GMT(nxmax,nymax,iprod,iout3d,OutVar,writeContours)
        case default
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: Plots requested but no plotting package is installed"
          endif;enddo
          stop 1
        end select

        endif
      elseif(outformat.eq.4)then ! Binary
        if(iprod.eq.1)then
          ! full concentration array but here we only output the total
          call write_3D_Binary(cio,nxmax,nymax,nzmax,ashcon_tot)
        elseif(iprod.eq.2)then
          ! deposit granularity
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: No binary output products for deposit granularity'
          endif;enddo
          stop 1
        elseif(iprod.eq.8)then
          ! ashfall arrival at airports/POI (mm)
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: No binary output products for POI ashfall arrival'
          endif;enddo
          stop 1
        elseif(iprod.eq.16)then
          ! profile plots
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: No binary output products for vertical profile plots'
          endif;enddo
          stop 1
        else
          call write_2D_Binary(nxmax,nymax,OutVar,mask,Fill_Value,filename_root)
        endif
      elseif(outformat.eq.5)then ! Shapefile
        ! For 2d contours exported from dislin, gnuplot, gmt
        ! First call plotting routine, but only get the contours
        writeContours = .true.

        icase = 0
        do ii=1,Nplot_libs
          ! Check each preference in series and see if the library is available
          if(plotlib_avail(plot_pref_shp(ii)))then
            icase = plot_pref_shp(ii)
            exit
          endif
        enddo

        select case (icase)
        case(1)
#ifdef USEDISLIN
          call write_2Dmap_PNG_dislin(nxmax,nymax,iprod,iout3d,OutVar,writeContours)
#endif
        case(2)
#ifdef USEPLPLOT
          call write_2Dmap_PNG_plplot(nxmax,nymax,iprod,iout3d,OutVar,writeContours)
#endif
        case(3)
          call write_2Dmap_PNG_gnuplot(nxmax,nymax,iprod,iout3d,OutVar,writeContours)
        case(4)
          write(*,*)'calling write_2Dmap_PNG_GMT'
          stop 66
          call write_2Dmap_PNG_GMT(nxmax,nymax,iprod,iout3d,OutVar,writeContours)
        case default
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: Plots requested but no plotting package is installed"
          endif;enddo
          stop 1
        end select

        call write_ShapeFile_Polyline(iprod,iout3d)
        !  For contours that follow topography, use
        !call write_ShapeFile_PolylineZ  (this is a place-holder, not yet implemented)
      elseif(outformat.eq.6)then
        !call write_2D_grib2
      elseif(outformat.eq.7)then
        !call write_2D_netcdf
      elseif(outformat.eq.8)then
        !call write_2D_tecplot
      elseif(outformat.eq.9)then
        !call write_2D_vtk
      endif

      ! clean up memory
      if(allocated(OutVar)) deallocate(OutVar)
      if(allocated(mask))   deallocate(mask)

      call dealloc_arrays
      call deallocate_ASCII
      call deallocate_Binary

      ! Close log file
      close(fid_logfile)

      end program Ash3d_PostProc

