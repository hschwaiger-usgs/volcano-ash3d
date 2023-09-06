!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!  that defines the input type and output products.  This is not yet implemented.
!  
!  If more than one command-line argument is given, this program expects the following:
!    ./Ash3d_PostProc Ash3d_output.nc output_product_ID format_code [time_step]
!  where output_product_ID is one of:
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
!  the final time data.  For data pproducts that require contours (shape files or
!  contour plots), contour levels are defined in the module Output_Vars.
!  
!  For example, to product a contour map of ash-cloud arrival time:
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
         MM_2_IN

      use Ash3d_Program_Control, only : &
           Set_OS_Env

      use mesh,          only : &
         nxmax,nymax,nzmax

      use time_data,     only : &
         time,time_native,BaseYear,useLeap,SimStartHour

      use io_data,       only : &
         iTimeNext, &
         concenfile,nWriteTimes,WriteTimes,Write_PT_Data,Write_PR_Data,&
         iout3d,isFinal_TS,WriteGSD,WriteDepositTS_KML,WriteDepositTS_ASCII,&
         WriteDepositTime_KML,WriteDepositTime_ASCII,WriteDepositFinal_KML,&
         WriteDepositFinal_ASCII,WriteCloudTime_KML,WriteCloudTime_ASCII,WriteReflectivity_ASCII,&
         WriteCloudLoad_KML,WriteReflectivity_KML,WriteCloudLoad_ASCII,WriteCloudHeight_KML,&
         WriteCloudHeight_ASCII,WriteCloudConcentration_KML,WriteCloudConcentration_ASCII,&
         WriteAirportFile_KML,WriteAirportFile_ASCII,Write3dFiles,&
         nvprofiles

      use Output_Vars,   only : &
         DepositThickness,DepArrivalTime,CloudArrivalTime,ashcon_tot,&
         MaxConcentration,MaxHeight,CloudLoad,dbZCol,MinHeight,Mask_Cloud

      use help,          only : &
           help_postproc

      use Ash3d_ASCII_IO,  only : &
           write_2D_ASCII, &
           write_3D_ASCII, &
           vprofileopener, &
           vprofilewriter, &
           vprofilecloser, &
           Write_PointData_Airports_ASCII

      use Ash3d_Binary_IO, only : &
           write_2D_Binary, &
           write_3D_Binary

      use Ash3d_Netcdf_IO

      use Ash3d_KML_IO

#ifdef USEDISLIN
      use Ash3d_PostProc_dislin
#endif
#ifdef USEPLPLOT
      use Ash3d_PostProc_plplot
#endif
      use Ash3d_PostProc_gnuplot
#ifdef USEGMT
      use Ash3d_PostProc_GMT
#endif

      implicit none

      integer             :: nargs
      integer             :: status
      integer             :: istat
      character (len=100) :: arg
      integer             :: iformat
      integer             :: iprod
      integer             :: ivar,TS_Flag,height_flag
      integer             :: itime = -1      ! initialize time step to the last step
      integer             :: i,ii
      integer             :: icase
      real(kind=ip),dimension(:,:),allocatable :: OutVar
      logical      ,dimension(:,:),allocatable :: mask
      real(kind=ip)       :: OutFillValue
      logical             :: IsThere
      character(len=6)    :: Fill_Value
      character(len=20)   :: filename_root
      character(len=70)   :: comd
      character(len=13)   :: cio
      logical             :: writeContours = .false.
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
      integer,dimension(Nplot_libs) :: plot_pref_map = (/2,1,3,4/) ! plot preference for maps
      integer,dimension(Nplot_libs) :: plot_pref_shp = (/3,1,2,4/) ! plot preference for contours
      integer,dimension(Nplot_libs) :: plot_pref_vpr = (/1,2,3,4/) ! plot preference for vert profs.
      integer,dimension(Nplot_libs) :: plot_pref_aTS = (/2,3,1,4/) ! plot preference for Airport TS

      INTERFACE
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
      END INTERFACE

      open(unit=global_log,file='Ash3d_pp.log',status='unknown')
      call Set_OS_Env

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
      call execute_command_line("echo 'exit' | gnuplot",exitstat=istat)
      if (istat.eq.0)then
        plotlib_avail(3) = .true.
      else
        plotlib_avail(3) = .false.
      endif
      ! Test for GMT
      plotlib_avail(4) = .false.
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Gnuplot ",plotlib_avail(1)
        write(outlog(io),*)"Plplot  ",plotlib_avail(2)
        write(outlog(io),*)"Dislin  ",plotlib_avail(3)
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

      ! Test read command-line arguments
      nargs = command_argument_count()
      if (nargs.eq.0) then
          ! If no command-line arguments are given, then prompt user
          ! interactively for the command file name and possible a 
          ! restart file
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)'No command-line arguments detected'
          write(outlog(io),*)'  '
        endif;enddo
        call help_postproc

        if(VB(1).ge.verbosity_silent)then
          do io=1,2
            write(errlog(io),*)"Stdout is suppressed via verbosity=9,10, but interactive input is expected."
            write(errlog(io),*)"Either recompile with verbosity<9 or provide the correct command-line arguments."
          enddo
          stop 1
        else
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)'Enter name of netcdf output file:'
          endif;enddo
        endif
        read(5,*)concenfile

        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)'Enter code for output product:'
        endif;enddo
        read(5,*)iprod

        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)'Enter code for output format:'
        endif;enddo
        read(5,*)iformat

        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)'Enter index for time step:'
        endif;enddo
        read(5,*)itime

      elseif (nargs.eq.1) then
        ! Read control file
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)'Reading control file not yet implemented'
        endif;enddo
        stop 1
      elseif (nargs.ge.3) then
        call get_command_argument(1, arg, status)
        read(arg,*)concenfile
        call get_command_argument(2, arg, status)
        read(arg,*)iprod
        call get_command_argument(3, arg, status)
        read(arg,*)iformat
        if (nargs.eq.4) then
          call get_command_argument(4, arg, status)
          read(arg,*)itime
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
      !  Arg #2
      if(iprod.lt.1.or.iprod.gt.16)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: output product requested is not in range 1-16."
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
            write(errlog(io),*)' Currently, no output formats available for iprod=15'
          endif;enddo
          stop 1
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
      if(iformat.lt.1.or.iformat.gt.7)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: output format requested is not in range 1-7."
          write(errlog(io),*)"       Run Ash3d_PostProc with no command-line"
          write(errlog(io),*)"       arguments to set usage information"
        endif;enddo
        stop 1
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          if(iformat.eq.1)then
            write(outlog(io),*)'output format = 1 ASCII/ArcGIS'
          elseif(iformat.eq.2)then
            write(outlog(io),*)'output format = 2 KMZ'
          elseif(iformat.eq.3)then
            write(outlog(io),*)'output format = 3 image/png'
          elseif(iformat.eq.4)then
            write(outlog(io),*)'output format = 4 binary'
          elseif(iformat.eq.5)then
            write(outlog(io),*)'output format = 5 shape file'
          elseif(iformat.eq.6)then
            write(outlog(io),*)'output format = 6 grib2'
          elseif(iformat.eq.7)then
            write(outlog(io),*)'output format = 7 tecplot'
          endif
        endif;enddo
      endif
      !  Arg #4
      if(itime.eq.-1)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)'itime = -1 (Final time step)'
          write(outlog(io),*)'Either no time step was provided or the last time step is requested.'
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
        write(outlog(io),*)'Finished reading command line'
      endif;enddo

      ! Before we do anything, call routine to read the netcdf file, populate
      ! the dimensions so we can see what we are dealing with.
      ! This call reads the 2d output products for the specified time
      ! If itime=-1 (for the final step), then
      call NC_Read_Output_Products(itime)

      if(.not.IsLatLon)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)'Only post-processing of Lon/Lat grids is currently supported.'
        endif;enddo
        stop 1
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
        if(    iformat.eq.1)then
          WriteGSD                      = .true.
        elseif(iformat.eq.2)then
          WriteGSD                      = .true.
        endif
      elseif(iprod.eq.3 )then ! deposit at specified times (mm)
        if(    iformat.eq.1)then
          WriteDepositTS_ASCII          = .true.
        elseif(iformat.eq.2)then
          WriteDepositTS_KML            = .true.
        endif
      elseif(iprod.eq.4 )then ! deposit at final time (mm)
        if(    iformat.eq.1)then
          WriteDepositFinal_ASCII       = .true.
        elseif(iformat.eq.2)then
          WriteDepositFinal_KML         = .true.
        endif
      elseif(iprod.eq.5 )then ! deposit at specified times (inches)
        if(    iformat.eq.1)then
          WriteDepositTS_ASCII          = .true.
        elseif(iformat.eq.2)then
          WriteDepositTS_KML            = .true.
        endif
      elseif(iprod.eq.6 )then ! deposit at final time (inches)
        if(    iformat.eq.1)then
          WriteDepositFinal_ASCII       = .true.
        elseif(iformat.eq.2)then
          WriteDepositFinal_KML         = .true.
        endif
      elseif(iprod.eq.7 )then ! ashfall arrival time
        if(    iformat.eq.1)then
          WriteDepositTime_ASCII        = .true.
        elseif(iformat.eq.2)then
          WriteDepositTime_KML          = .true.
        endif
      elseif(iprod.eq.8 )then ! ashfall at airports/POI
        Write_PR_Data                   = .true.
        if(    iformat.eq.1)then
          WriteAirportFile_ASCII        = .true.
        elseif(iformat.eq.2)then
          WriteAirportFile_KML          = .true.
        endif
      elseif(iprod.eq.9 )then ! ash-cloud concentration
        if(    iformat.eq.1)then
          WriteCloudConcentration_ASCII = .true.
        elseif(iformat.eq.2)then
          WriteCloudConcentration_KML   = .true.
        endif
      elseif(iprod.eq.10)then ! ash-cloud height
        if(    iformat.eq.1)then
          WriteCloudHeight_ASCII        = .true.
        elseif(iformat.eq.2)then
          WriteCloudHeight_KML          = .true.
        endif
      elseif(iprod.eq.11)then ! ash-cloud bottom
        if(    iformat.eq.1)then
          WriteCloudHeight_ASCII        = .true.
        elseif(iformat.eq.2)then
          WriteCloudHeight_KML          = .true.
        endif
      elseif(iprod.eq.12)then ! ash-cloud load
        if(    iformat.eq.1)then
          WriteCloudLoad_ASCII          = .true.
        elseif(iformat.eq.2)then
          WriteCloudLoad_KML            = .true.
        endif
      elseif(iprod.eq.13)then ! ash-cloud radar reflectivity
        if(    iformat.eq.1)then
          WriteReflectivity_ASCII       = .true.
        elseif(iformat.eq.2)then
          WriteReflectivity_KML         = .true.
        endif
      elseif(iprod.eq.14)then ! ash-cloud arrival time
        if(    iformat.eq.1)then
          WriteCloudTime_ASCII          = .true.
        elseif(iformat.eq.2)then
          WriteCloudTime_KML            = .true.
        endif
      !elseif(iprod.eq.15)then ! topography
      !  if(    iformat.eq.1)then
      !  elseif(iformat.eq.2)then
      !  endif
      elseif(iprod.eq.16)then ! vertical profiles of concentration
        if(    iformat.eq.1.or.iformat.eq.3)then
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

      ! The main differences in output products will be time-series, vs
      ! time-step output.  Currently, the time-series output will only be for
      ! the KML/KMZ files.
      if(iformat.eq.2)then
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
          !elseif(iprod.eq.15)then
          !  OutVar = Topography(1:nxmax,1:nymax)
          endif
          call Write_2D_KML(ivar,OutVar,height_flag,TS_flag)
          call Close_KML(ivar,TS_flag)
          do io=1,2;if(VB(io).le.verbosity_info)then          
            write(outlog(io),*)"Zipping KML file."
          endif;enddo
          write(comd,*)"zip ",trim(adjustl(KMZ_filename(ivar))),' ',&
                              trim(adjustl(KML_filename(ivar)))
          call execute_command_line (comd, exitstat=status)
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
            call NC_Read_Output_Products(i)
            time = WriteTimes(i)
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"time = ",time
            endif;enddo
            call output_results
          enddo
          isFinal_TS = .true.
          ! For deposit output, load the final deposit variable and write to file
          call NC_Read_Output_Products(-1)
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
      endif ! iformat.eq.2 (KML)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! This is the non-KML section
      ! Set some variable parameters for an ASCII output file
      OutFillValue = 0.0_ip
      if(iprod.eq.3.or.iprod.eq.5)then
        OutVar = DepositThickness
        Fill_Value = ' 0.000'
        OutFillValue = 0.0_ip
        filename_root = 'DepositFile_        '
      elseif(iprod.eq.4.or.iprod.eq.6)then
        OutVar = DepositThickness*MM_2_IN
        Fill_Value = ' 0.000'
        OutFillValue = 0.0_ip
        filename_root = 'DepositFile_        '
      elseif(iprod.eq.7)then
        OutVar = real(DepArrivalTime,kind=ip)
        Fill_Value = '-1.000'
        OutFillValue = -1.0_ip
        filename_root = 'DepositArrivalTime  '
      elseif(iprod.eq.8)then
         ! ashfall at airports/POI
         ! None of OutVar,Fill_Value,OutFillValue,filename_root need to be set
      elseif(iprod.eq.9)then
        OutVar = MaxConcentration
        Fill_Value = ' 0.000'
        OutFillValue = 0.0_ip
        filename_root = 'CloudConcentration_ '
      elseif(iprod.eq.10)then
        OutVar = MaxHeight
        Fill_Value = ' 0.000'
        OutFillValue = 0.0_ip
        filename_root = 'CloudHeight_        '
      elseif(iprod.eq.11)then
        OutVar = MinHeight
        Fill_Value = ' 0.000'
        OutFillValue = 0.0_ip
        filename_root = 'CloudHeightBot_     '
      elseif(iprod.eq.12)then
        OutVar = CloudLoad
        Fill_Value = ' 0.000'
        OutFillValue = 0.0_ip
        filename_root = 'CloudLoad_          '
      elseif(iprod.eq.13)then
        OutVar = dbZCol
        Fill_Value = ' 0.000'
        OutFillValue = 0.0_ip
        filename_root = 'ClouddbZC_          '
      elseif(iprod.eq.14)then
        OutVar = real(CloudArrivalTime,kind=ip)
        Fill_Value = '-1.000'
        OutFillValue = -1.0_ip
        filename_root = 'CloudArrivalTime    '
      endif
      ! Now mask out non-cloud values
      if(iprod.eq.9 .or.&
         iprod.eq.10.or.&
         iprod.eq.11.or.&
         iprod.eq.12.or.&
         iprod.eq.13.or.&
         iprod.eq.14)then
        mask = Mask_Cloud
      endif

      ! This is the ASCII section
      if(iformat.eq.1)then
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
          do itime=1,tn_len
            time = time_native(itime)
            call vprofilewriter(itime)
          enddo
          call vprofilecloser
        else
          ! All other ESRI/ASCII 2d grids
          call write_2D_ASCII(nxmax,nymax,OutVar,mask,Fill_Value,filename_root)
        endif
      elseif(iformat.eq.2)then
        ! All the KML routines were called above
      elseif(iformat.eq.3)then ! image/png
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
#ifdef USEGMT
              !call write_2Dprof_PNG_GMT(i)
#endif
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
#ifdef USEGMT
          call write_2Dmap_PNG_GMT(nxmax,nymax,iprod,iout3d,OutVar,writeContours)
#endif
        case default
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: Plots requested but no plotting package is installed"
          endif;enddo
          stop 1
        end select

        endif
      elseif(iformat.eq.4)then
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
      elseif(iformat.eq.5)then
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
#ifdef USEGMT
          call write_2Dmap_PNG_GMT(nxmax,nymax,iprod,iout3d,OutVar,writeContours)
#endif
        case default
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: Plots requested but no plotting package is installed"
          endif;enddo
          stop 1
        end select

        call write_ShapeFile_Polyline(iprod,iout3d)
        !  For contours that follow topography, use
        !call write_ShapeFile_PolylineZ  (this is a place-holder, not yet implemented)
      elseif(iformat.eq.6)then
        !call write_2D_grib2
      elseif(iformat.eq.7)then
        !call write_2D_netcdf
      elseif(iformat.eq.8)then
        !call write_2D_tecplot
      elseif(iformat.eq.9)then
        !call write_2D_vtk
      endif

      end program Ash3d_PostProc

