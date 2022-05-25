      program Ash3d_PostProc

      use precis_param

      use io_units

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,lon_cc_pd,lat_cc_pd

      use io_data,       only : &
         concenfile,nWriteTimes,WriteTimes,cdf_b3l1

      use Output_Vars,   only : &
         DepositThickness,DepArrivalTime,CloudArrivalTime,&
         MaxConcentration,MaxHeight,CloudLoad,dbZCol,MinHeight,Mask_Cloud,Mask_Deposit

      use Airports,      only : &
         Airport_Thickness_TS,Airport_Name

      use Ash3d_Netcdf

      use plplot
      use dislin

      implicit none

      integer             :: nargs
      integer             :: status
      character (len=100) :: arg
      integer             :: iformat
      integer             :: iprod
      integer             :: itime = -1      ! initialize time step to the last step
      real(kind=ip),dimension(:,:),allocatable :: OutVar
      logical             :: IsThere


!      character(len=14) :: dp_outfile
!      character(len=14) :: dp_gnufile
!      character(len=14) :: dp_pngfile1
!      character(len=14) :: dp_pngfile2
!
      integer :: plt_indx
!      character(len=12) :: Airport_Name
!      real(kind=8) :: ymaxpl
!
!
!      ! PLPLOT stuff
!      !  http://plplot.sourceforge.net/
!      real(kind=plflt) :: xmin
!      real(kind=plflt) :: xmax
!      real(kind=plflt) :: ymin
!      real(kind=plflt) :: ymax
!      real(kind=plflt), dimension(:), allocatable :: x, y, x0, y0
!      integer(kind=4) :: r1 = 0
!      integer(kind=4) :: g1 = 0
!      integer(kind=4) :: b1 = 0
!      integer(kind=4) :: r2 = 136
!      integer(kind=4) :: g2 = 136
!      integer(kind=4) :: b2 = 136
!
!
      ! TEST READ COMMAND LINE ARGUMENTS
      nargs = command_argument_count()
      if (nargs.eq.0) then
          ! If no command-line arguments are given, then prompt user
          ! interactively for the command file name and possible a 
          ! restart file
        write(global_info,*)'No command-line arguments detected'
        write(global_info,*)'Usage: Ash3d_PostProc control_file [t_index]'
        write(global_info,*)'           or'
        write(global_info,*)'       Ash3d_PostProc infile output_product format'
        write(global_info,*)'  where: infile   = the netcdf file written by Ash3d'
        write(global_info,*)'   output_product = 1 full concentration array'
        write(global_info,*)'                    2 deposit granularity'
        write(global_info,*)'                    3 final deposit thickness'
        write(global_info,*)'                    4 deposit at specified times'
        write(global_info,*)'                    5 ash-cloud concentration'
        write(global_info,*)'                    6 ash-cloud height'
        write(global_info,*)'                    7 ash-cloud bottom'
        write(global_info,*)'                    8 ash-cloud load'
        write(global_info,*)'                    9 deposit arrival time'
        write(global_info,*)'                   10 ashfall arrival time'
        write(global_info,*)'                   11 ash-cloud arrival time'
        write(global_info,*)'                   12 radar reflectivity'
        write(global_info,*)'                   13 ash arrival at airports/POI'
        write(global_info,*)'                   14 profile plots'
        write(global_info,*)'           format = 1 ASCII/ArcGIS'
        write(global_info,*)'                    2 KMZ'
        write(global_info,*)'                    3 image/png'
        write(global_info,*)'                    4 binary'
        write(global_info,*)'                    5 shape file'
        write(global_info,*)'                    6 grib2'
        write(global_info,*)'                    7 tecplot'
        write(global_info,*)'         [t_index] = index of time slice to plot (optional)'
        write(global_info,*)'  '

        write(global_info,*)'Enter name of ESP input file:'
        read(5,*)concenfile

        write(global_info,*)'Enter code for output product:'
        read(5,*)iprod

        write(global_info,*)'Enter code for output format:'
        read(5,*)iformat

        write(global_info,*)'Enter index for time step:'
        read(5,*)iformat

      elseif (nargs.eq.1) then
        ! Read control file
        write(*,*)'Reading control file not yet implemented'
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
        endif
      else
        write(global_info,*)' Cannot parse command line'
        stop 1
      endif

      ! Now that we have read the command line, error-check and report back what
      ! we are about to do.
      !  Arg #1
      inquire( file=concenfile, exist=IsThere )
      if(.not.IsThere)then
        write(global_error,*)"ERROR: Cannot file input file"
        stop 1
      else
        write(global_info,*)"Ash3d Output file: ",concenfile
      endif
      !  Arg #2
      if(iprod.lt.1.or.iprod.gt.14)then
        write(global_error,*)"ERROR: output product requested is not in range 1-14."
        write(global_error,*)"       Run Ash3d_PostProc with no command-line"
        write(global_error,*)"       arguments to set usage information"
        stop 1
      else
        if(iprod.eq.1)then
          write(global_info,*)'output variable = 1 full concentration array'
          write(global_info,*)' Currently, no output formats available for iprod=1'
          stop 1
        elseif(iprod.eq.2)then
          write(global_info,*)'output variable = 2 deposit granularity'
          write(global_info,*)' Currently, no output formats available for iprod=2'
          stop 1
        elseif(iprod.eq.3)then
          write(global_info,*)'output variable = 3 final deposit thickness'
        elseif(iprod.eq.4)then
          write(global_info,*)'output variable = 4 deposit at specified times'
        elseif(iprod.eq.5)then
          write(global_info,*)'output variable = 5 ash-cloud concentration'
        elseif(iprod.eq.6)then
          write(global_info,*)'output variable = 6 ash-cloud height'
        elseif(iprod.eq.7)then
          write(global_info,*)'output variable = 7 ash-cloud bottom'
        elseif(iprod.eq.8)then
          write(global_info,*)'output variable = 8 ash-cloud load'
        elseif(iprod.eq.9)then
          write(global_info,*)'output variable = 9 deposit arrival time'
        elseif(iprod.eq.10)then
          write(global_info,*)'output variable =10 ashfall arrival time'
        elseif(iprod.eq.11)then
          write(global_info,*)'output variable =11 ash-cloud arrival time'
        elseif(iprod.eq.12)then
          write(global_info,*)'output variable =12 radar reflectivity'
        elseif(iprod.eq.13)then
          write(global_info,*)'output variable =13 ash arrival at airports/POI'
        elseif(iprod.eq.14)then
          write(global_info,*)'output variable =14 profile plots'
          write(global_info,*)' Currently, no output formats available for iprod=14'
          stop 1
        endif
      endif
      !  Arg #3
      if(iformat.lt.1.or.iformat.gt.7)then
        write(global_error,*)"ERROR: output format requested is not in range 1-7."
        write(global_error,*)"       Run Ash3d_PostProc with no command-line"
        write(global_error,*)"       arguments to set usage information"
        stop 1
      else
        if(iformat.eq.1)then
          write(global_info,*)'output format = 1 ASCII/ArcGIS'
        elseif(iformat.eq.2)then
          write(global_info,*)'output format = 2 KMZ'
        elseif(iformat.eq.3)then
          write(global_info,*)'output format = 3 image/png'
        elseif(iformat.eq.4)then
          write(global_info,*)'output format = 4 binary'
        elseif(iformat.eq.5)then
          write(global_info,*)'output format = 5 shape file'
        elseif(iformat.eq.6)then
          write(global_info,*)'output format = 6 grib2'
        elseif(iformat.eq.7)then
          write(global_info,*)'output format = 7 tecplot'
        endif
      endif
      !  Arg #4
      if(itime.eq.-1)then
        write(global_info,*)'itime = -1 (Final time step)'
        write(global_info,*)'Either no time step was provided or the last time step is requested.'
      elseif(itime.eq.0)then
        write(global_error,*)"ERROR: itime = 0.  Invalid time step."
        stop 1
      else
        write(global_info,*)'itime = ',itime
        write(global_info,*)'  We do not yet know the maximum number of steps available.'
      endif

      if(iprod.eq.1.or. &  ! full concentration array
         iprod.eq.2)then   ! deposit granularity
        write(global_info,*)'We do not yet have a plan for processing ashcon or depocon'
        write(global_info,*)'This is probably where we would implement binary, vtk, or tecplot'
        stop 1
      endif

      ! Before we do anything, call routine to read the netcdf file, populate
      ! the dimentions so we can see what we are dealing with.
      ! This call reads the 2d output products for the specified time
      ! If itime=-1 (for the final step), then
      call NC_Read_Output_Products(itime)

      if (itime.eq.-1)itime = nWriteTimes
      allocate(OutVar(nxmax,nymax))

      ! The main differences in output products with be time-series, vs
      ! time-step output.  Currently, the time-series output will only be for
      ! the KML/KMZ files.  
      if(iformat.eq.2)then
        ! KML output, separate into static vs time-series options
        ! Note: for KML, the following variable ID's are used:
        !        ivar = 1 :: cloud concentration
        !        ivar = 2 :: cloud height (top)
        !        ivar = 3 :: cloud height (bot)
        !        ivar = 4 :: cloud load
        !        ivar = 5 :: cloud arrival time
        !        ivar = 6 :: cloud reflectivity
        !        ivar = 7 :: deposit
        !        ivar = 8 :: deposit (NWS)
        !        ivar = 9 :: deposit time
        !        ivar =10 :: topography
        if(iprod.eq.3.or.  &  ! Final Deposit Thickness
           iprod.eq.9.or.  &  ! Deposit arrival time
           iprod.eq.10.or. &  ! Ashfall arrival time
           iprod.eq.11)then   ! Ash-cloud arrival time
          ! Static output
          ! We have already called the netcdf reader so we just need to write
          ! the KML output.
          ! Set up KML output files and prelim variables
          !call Set_OutVar_Specs
          !call OpenFile_KML(ivar)
          !call Write_2D_KML(ivar,OutVar,height_flag,TS_flag)
          !call Close_KML(ivar,TS_flag)

        elseif(iprod.eq.4.or. &  ! deposit at specified times
               iprod.eq.5.or. &  ! ash-cloud concentration
               iprod.eq.6.or. &  ! ash-cloud height
               iprod.eq.7.or. &  ! ash-cloud bottom
               iprod.eq.8.or. &  ! ash-cloud load
               iprod.eq.12.or.&  ! radar reflectivity
               iprod.eq.13)then  ! ash arrival at airports/POI
          ! Time-series output
          ! Set up KML output files and prelim variables
!          call Set_OutVar_Specs
!          call OpenFile_KML(ivar)
!          ! We need to loop over all times
!          do i = 1,nWriteTimes
!            call NC_Read_Output_Products(itime)
!            call Write_2D_KML(ivar,OutVar,height_flag,TS_flag)
!          enddo
!          if(iprod.eq.13)then
!            call Write_PointData_Airports_KML
!          endif
!          call Close_KML(ivar,TS_flag)

        elseif(iprod.eq.14)then  ! vertical profile plots
          write(global_error,*)"ERROR: KML versions of vertical profiles not implemented."
          stop 1
        else
          write(global_error,*)"ERROR: Requested iprod not implemented for KML."
          stop 1
        endif
      endif ! iformat.eq.2 (KML)

      ! 
      if(iprod.eq.3.or. &  ! deposit at final time
         iprod.eq.4.or. &  ! deposit at specified times
         iprod.eq.5.or. &  ! ash-cloud concentration
         iprod.eq.6.or. &  ! ash-cloud height
         iprod.eq.7.or. &  ! ash-cloud bottom
         iprod.eq.8.or. &  ! ash-cloud load
         iprod.eq.12)then   ! radar reflectivity

        if(iprod.eq.3.or.iprod.eq.4)then
          OutVar = DepositThickness
         elseif(iprod.eq.5)then
          OutVar = MaxConcentration
         elseif(iprod.eq.6)then
          OutVar = MaxHeight
         elseif(iprod.eq.7)then
          OutVar = MinHeight
         elseif(iprod.eq.8)then
          OutVar = CloudLoad
         elseif(iprod.eq.12)then
          OutVar = dbZCol
         endif

        if(iformat.eq.1)then
          !call write_2D_ASCII(nx,ny,OutVar,Fill_Value,filename_root)
        elseif(iformat.eq.2)then
          !call Write_2D_KML(ivar,OutVar,height_flag,TS_flag)
        elseif(iformat.eq.3)then
          call write_2D_PNG_dislin(iprod,itime,OutVar)
        elseif(iformat.eq.4)then
          !call write_2D_Binary
        elseif(iformat.eq.5)then
          !call write_2D_ShapeFile
        elseif(iformat.eq.6)then
          !call write_2D_grib2
        elseif(iformat.eq.7)then
          !call write_2D_tecplot
        endif
      endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !  Pretend we have the airport data loaded and set up to call the
       !  different plotting routines.

      !allocate(Airport_Thickness_TS(nWriteTimes))

      do plt_indx = 1,nWriteTimes
        call write_DepPOI_TS_PNG_gnuplot(plt_indx)
      enddo


!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      ! PLPLOT block
!      write(dp_pngfile2,55) plt_indx,".png"
! 55   format('plplt_',i4.4,a4)
!
!      xmin=real(0,kind=plflt)
!      xmax=real(ceiling(Simtime_in_hours),kind=plflt)
!      ymin=real(0,kind=plflt)
!      ymax=real(ymaxpl,kind=plflt)
!
!      allocate(x(nWriteTimes))
!      allocate(y(nWriteTimes))
!      allocate(x0(nWriteTimes+1))
!      allocate(y0(nWriteTimes+1))
!
!      x = real(WriteTimes,kind=plflt)
!      y = real(Airport_Thickness_TS,kind=plflt)
!      x0(1:nWriteTimes)=x(1:nWriteTimes)
!      x0(nWriteTimes+1)=x(nWriteTimes)
!      y0(1:nWriteTimes)=y(1:nWriteTimes)
!      y0(nWriteTimes+1)=0.0_plflt
!
!      ! Set up for plplot
!      call plsdev("pngcairo")      ! Set output device (png, pdf, etc.)
!      call plsfnam ( dp_pngfile2 ) ! Set output filename
!
!      call plsetopt("geometry","400x300, 400x300")  ! Set image size
!      call plsetopt("bg","FFFFFF")                  ! Set background color to white
!      ! Initialize plplot
!      call plinit()
!
!      ! Create a labelled box to hold the plot.
!      call plscolbg(r1,g1,b1) ! Set color back to black
!      call plcol0(0)
!      call plschr(0.0_plflt,1.7_plflt)  ! Chage font scale
!      call plenv( xmin, xmax, ymin, ymax, 0, 0 )
!      call pllab( "Time (hours after eruption)", "Deposit Thickeness (mm)", Airport_Name )
!
!      call plscolbg(r2,g2,b2) ! Set pen color to grey
!      call plcol0(0)
!      call plline( x , y ) ! Simple line plot
!      call plfill( x0(1:nWriteTimes+1), y0(1:nWriteTimes+1) )
!
!      ! Close PLplot library
!      call plend
!      ! end PLPLOT block
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      ! DISLIN block
!      ! https://www.dislin.de/
!      ! wget https://ftp.gwdg.de/pub/grafik/dislin/linux/i586_64/dislin-11.4.linux.i586_64.tar.gz
!      !  Dislin Level 0:  before initialization or after termination
!      call metafl(CFMT)   ! set output driver/file-format (PNG); this is a 4-char string
!      call setpag('USAL') ! Set pagesize to US A Landscape (2790 x 2160)
!      call setfil('temp.png') ! Set output filename
!      call scrmod('REVERSE')  ! Default background is black; reverse to white
!
!      !  Dislin Level 1:  after initialization or a call to ENDGRF
!      call disini()       ! initialize plot (set to level 1)
!        ! setting of plot parameters
!      call pagera()       ! plot a border around the page
!      call triplx()  ! set font to triple stroke
!      call axspos(450,1800)  ! determine the position of the axis system
!      call axslen(2200,1200) ! defines the size of the axis system
!      call name('Time (hours after eruption)','X') ! Set x-axis title
!      call name('Deposit Thickeness (mm)','Y') ! Set y-axis title
!      call labdig(-1,'X') ! set number of decimal places for x label (-1 means no decimal)
!      call ticks(10,'XY') ! set number of ticks between labels
!      call titlin(Airport_Name,4)  ! Set the title to the airport name (4 is the bottom line)
!
!      !  Dislin Level 2: after a call to GRAF, GRAFP or GRAFMP
!        ! Now create graph and set to level 2
!      call graf(real(xmin,kind=4), real(xmax,kind=4), 0.0_4, 5.0_4, &
!                real(ymin,kind=4), real(ymax,kind=4), 0.0_4, 1.0_4)
!      call title() ! Actually write the title to the file
!      call setrgb(0.5_4, 0.5_4, 0.5_4)
!      !call curve(real(x,kind=4),real(y,kind=4),nWriteTimes)  ! This draws the line
!      call shdpat(16)  ! set shading pattern 
!      call shdcrv(real(x,kind=4),real(y,kind=4),nWriteTimes,& ! This fills below curve
!                  real(x,kind=4),real(0.0*y,kind=4),nWriteTimes)
!      call color('FORE') ! Reset color to defaul foreground color
!
!      !  Dislin Level 0:  before initialization or after termination
!      call disfin()
!      endif

!     (1)    setting of page format, file format and filename
!     (2)    initialization
!     (3)    setting of plot parameters
!     (4)    plotting of the axis system
!     (5)    plotting the title
!     (6)    plotting data points
!     (7)    termination.

      end program Ash3d_PostProc

!##############################################################################
!
!    write_2D_PNG_dislin
!
!    if timestep = -1, then use the last step in file
!##############################################################################

      subroutine write_2D_PNG_dislin(iprod,itime,OutVar)

      use precis_param

      use mesh,          only : &
         nxmax,nymax,x_cc_pd,y_cc_pd,lon_cc_pd,lat_cc_pd, &
         IsLatLon !,dx,dy,dz_vec_pd,nzmax,nsmax,z_cc_pd

      use Output_Vars,   only : &
         DepositThickness,DepArrivalTime,CloudArrivalTime,&
         MaxConcentration,MaxHeight,CloudLoad,dbZ,MinHeight,Mask_Cloud,Mask_Deposit

      use io_data,       only : &
         nWriteTimes,WriteTimes,cdf_b3l1,VolcanoName

      use Source,        only : &
         neruptions,e_Volume,e_Duration,e_StartTime,e_PlumeHeight

      use time_data,     only : &
         cdf_time_log,BaseYear,useLeap

      use dislin

      implicit none

      integer :: iprod
      integer :: itime
      real(kind=ip) :: OutVar(nxmax,nymax)

      integer :: i

      ! dislin stuff
      ! https://www.dislin.de/
      CHARACTER (LEN=4) :: CFMT = "PNG "
      character(len=40) :: outfile_name
      character (len=9) :: cio
      character (len=4) :: outfile_ext = '.png'

      INTEGER, PARAMETER :: N=3
      INTEGER, PARAMETER :: MAXPTS=1000
      INTEGER, PARAMETER :: MAXCRV=10
      REAL(kind=4), DIMENSION (N) :: &
         XC = (/-122.3167_4,-122.6417_4,-122.96310_4/), &
         YC = (/  47.5886_4,  45.4421_4,  49.2743_4/)
      CHARACTER (LEN=12), DIMENSION (N) :: &
         CSTR = (/'Seattle     ', 'Portland    ', 'Vancouver   '/)
      CHARACTER(len=80) :: CBUF
      integer :: NMAXLN
      character(len=7) :: zlevlab
      real(kind=4) :: XPTS(MAXPTS),YPTS(MAXPTS)
      INTEGER :: IRAY(MAXCRV)
      INTEGER :: NXP,NYP,NCLR,NCURVS
      REAL(kind=4)    :: XP,YP
      real(kind=4) :: xminDIS, xmaxDIS, yminDIS, ymaxDIS
      real(kind=4) :: dx_map, dy_map, xgrid_1, ygrid_1
      real(kind=4), dimension(11) :: depthic_zlev = &
         (/0.01_4, 0.03_4, 0.1_4, 0.3_4, 1.0_4, 3.0_4, &
           10.0_4, 30.0_4, 100.0_4, 300.0_4, 1000.0_4/)

      integer :: nzlev
      real(kind=4), dimension(:),allocatable :: zlev 

      character(len=40) :: title_plot
      character(len=15) :: title_legend

      character(len=30) :: cstr_volcname
      character(len=30) :: cstr_run_date
      character(len=30) :: cstr_windfile
      character(len=40) :: cstr_ErStartT
      character(len=27) :: cstr_ErHeight
      character(len=30) :: cstr_ErDuratn
      character(len=38) :: cstr_ErVomune
      character(len=45) :: cstr_note
      integer :: y_footer
      integer :: ioerr,iw,iwf

      INTERFACE
        character (len=20) function HS_xmltime(HoursSince,byear,useLeaps)
          real(kind=8)              :: HoursSince
          integer                   :: byear
          logical                   :: useLeaps
        end function HS_xmltime
      END INTERFACE

      if(iprod.eq.3)then
        cio='____final'
      else
        if (WriteTimes(itime).lt.10.0_ip) then
          write(cio,1) WriteTimes(itime)
1         format('00',f4.2,'hrs')
        elseif (WriteTimes(itime).lt.100.0_ip) then
          write(cio,2) WriteTimes(itime)
2         format('0',f5.2,'hrs')
        else
          write(cio,3) WriteTimes(itime)
3         format(f6.2,'hrs')
        endif
      endif

      if(iprod.eq.3)then       ! deposit at final time
        write(outfile_name,'(a13,a9,a4)')'Ash3d_Deposit',cio,outfile_ext
        title_plot = 'Final Deposit Thickness'
        title_legend = 'Dep.Thick.(mm)'
        nzlev = 11
        allocate(zlev(nzlev))
        zlev = (/0.01_4, 0.03_4, 0.1_4, 0.3_4, 1.0_4, 3.0_4, &
                10.0_4, 30.0_4, 100.0_4, 300.0_4, 1000.0_4/)
      elseif(iprod.eq.4)then   ! deposit at specified times
        write(outfile_name,'(a15,a9,a4)')'Ash3d_Deposit_t',cio,outfile_ext
        write(title_plot,'(a20,f5.2,a6)')'Deposit Thickness t=',WriteTimes(itime),' hours'
        title_legend = 'Dep.Thick.(mm)'
        nzlev = 11
        allocate(zlev(nzlev))
        zlev = (/0.01_4, 0.03_4, 0.1_4, 0.3_4, 1.0_4, 3.0_4, &
                10.0_4, 30.0_4, 100.0_4, 300.0_4, 1000.0_4/)
      elseif(iprod.eq.5)then   ! ash-cloud concentration
        write(outfile_name,'(a16,a9,a4)')'Ash3d_CloudCon_t',cio,outfile_ext
        write(title_plot,'(a26,f5.2,a6)')'Ash-cloud concentration t=',WriteTimes(itime),' hours'
        title_legend = 'Max.Con.(mg/m3)'
        nzlev = 8
        allocate(zlev(nzlev))
        zlev = (/0.1_4, 0.3_4, 1.0_4, 3.0_4, &
                10.0_4, 30.0_4, 100.0_4, 300.0_4/)
      elseif(iprod.eq.6)then   ! ash-cloud height
        write(outfile_name,'(a19,a9,a4)')'Ash3d_CloudHeight_t',cio,outfile_ext
        write(title_plot,'(a19,f5.2,a6)')'Ash-cloud height t=',WriteTimes(itime),' hours'
        title_legend = 'Cld.Height(km)'
        nzlev = 9
        allocate(zlev(nzlev))
        zlev = (/0.24_4, 3.0_4, 6.0_4, 10.0_4, 13.0_4, 16.0_4, &
                20.0_4, 25.0_4, 30.0_4/)
      elseif(iprod.eq.7)then   ! ash-cloud bottom
        write(outfile_name,'(a16,a9,a4)')'Ash3d_CloudBot_t',cio,outfile_ext
        write(title_plot,'(a19,f5.2,a6)')'Ash-cloud bottom t=',WriteTimes(itime),' hours'
        title_legend = 'Cld.Bot.(km)'
        nzlev = 9
        allocate(zlev(nzlev))
        zlev = (/0.24_4, 3.0_4, 6.0_4, 10.0_4, 13.0_4, 16.0_4, &
                20.0_4, 25.0_4, 30.0_4/)
      elseif(iprod.eq.8)then   ! ash-cloud load
        write(outfile_name,'(a17,a9,a4)')'Ash3d_CloudLoad_t',cio,outfile_ext
        write(title_plot,'(a17,f5.2,a6)')'Ash-cloud load t=',WriteTimes(itime),' hours'
        title_legend = 'Cld.Load(T/km2)'
        nzlev = 9
        allocate(zlev(nzlev))
        zlev = (/0.2_4, 1.0_4, 2.0_4, 5.0_4, 10.0_4, 30.0_4, &
                100.0_4, 300.0_4, 1000.0_4/)
      elseif(iprod.eq.12)then  ! radar reflectivity
        write(outfile_name,'(a20,a9,a4)')'Ash3d_CloudRadRefl_t',cio,outfile_ext
        write(title_plot,'(a24,f5.2,a6)')'Ash-cloud radar refl. t=',WriteTimes(itime),' hours'
        title_legend = 'Cld.Refl.(dBz)'
        nzlev = 9
        allocate(zlev(nzlev))
        zlev = (/-20._4, -10.0_4, 0.0_4, 10.0_4, 20.0_4, 30.0_4, &
                40.0_4, 50.0_4, 60.0_4/)
      endif

      ! Now map plots
      dx_map = 10.0_4
      dy_map = 5.0_4
      lon_cc_pd = lon_cc_pd -360.0_8
      xminDIS = real(minval(lon_cc_pd(1:nxmax)),kind=4)
      xmaxDIS = real(maxval(lon_cc_pd(1:nxmax)),kind=4)
      yminDIS = real(minval(lat_cc_pd(1:nymax)),kind=4)
      ymaxDIS = real(maxval(lat_cc_pd(1:nymax)),kind=4)
      xgrid_1 = real(ceiling(xminDIS/dx_map) * dx_map,kind=4)
      ygrid_1 = real(ceiling(yminDIS/dy_map) * dy_map,kind=4)

      !!!!!!!!!!!!!!!!!!!!!!!
      !  Dislin Level 0:  before initialization or after termination
      call metafl(CFMT)   ! set output driver/file-format (PNG); this is a 4-char string
      call setpag('USAL') ! Set pagesize to US A Landscape (2790 x 2160)
      call setfil(adjustl(trim(outfile_name))) ! Set output filename
      !call setvlt('SPEC')
      call scrmod('REVERSE')  ! Default background is black; reverse to white

      !  Dislin Level 1:  after initialization or a call to ENDGRF
      call disini()       ! initialize plot (set to level 1)
      call paghdr('Ash3d Simulation plotted on ','---',4,0)
      y_footer = 1900
      call filbox(2250,y_footer,130,49)
      call incfil('/opt/USGS/Ash3d/share/post_proc/USGSvid.png')
        ! setting of plot parameters
      call pagera()       ! plot a border around the page
      call triplx()  ! set font to triple stroke
      call axspos(500,1650)  ! determine the position of the axis system
      call axslen(2200,1400) ! defines the size of the axis system
      call name('Longitude','X') ! Set x-axis title
      call name('Latitude','Y') ! Set y-axis title
      call labdig(-1,'X') ! set number of decimal places for x label (-1 means no decimal)
      call ticks(1,'xy')  ! set number of ticks between labels
      call titlin(title_plot,4)  ! Set the title
      call incmrk(-1) ! selects line or symbol mode for CURVE

      !CALL LABELS('MAP','xy')
      !call projct('STER') ! defines projection
      call projct('LAMB') ! defines projection
      !call projct('MERC') ! defines projection
      call frame(3) ! bump up frame line thickness
       !  The routine GRAFMP plots a geographical axis system.
      call grafmp(xminDIS,xmaxDIS,xgrid_1,dx_map, &
                  yminDIS,ymaxDIS,ygrid_1,dy_map)

       ! set color of coastlines
      !call color('GREEN')
       ! plots coastlines and lakes or political borders
      call world()

      ! Add cities
      CALL CURVMP(XC,YC,N)
      DO I=1,N
      !    These are the points
        CALL POS2PT(XC(I),YC(I),XP,YP)
      !    These are the city lables, offset in x
        NXP=NINT(XP+30)
        NYP=NINT(YP)
        CALL MESSAG(CSTR(I),NXP,NYP)
      END DO

     !call myvlt(xr,xg,xb,nrgb)
     call shdmod('UPPER', 'CELL') ! This suppresses colors in regions above/below the zlevels pro
     call conshd(real(lon_cc_pd,kind=4),nxmax,&
                 real(lat_cc_pd,kind=4),nymax,&
                 real(OutVar,kind=4),zlev,nzlev)

       ! set color of grid lines
!      call setrgb(0.5_4, 0.5_4, 0.5_4)
       ! overlays an axis system with a longitude and latitude grid
      call gridmp(1,1)
      call height(50) ! Set character height for title
      call title() ! Actually write the title to the file

      ! Now write the legend
      call height(25) ! Reset character height to something smaller
      nmaxln = 6 ! number of characters in the longest line of text
      call legini(cbuf,nzlev,nmaxln) ! Initialize legend
      call legtit(title_legend)      ! Set legend title
      call legbgd(0)                 ! sets background color
      do i=1,nzlev
        if(zlev(i).lt.1.0_ip)then
          write(zlevlab,'(f6.2)')zlev(i)
        else
          write(zlevlab,'(f6.1)')zlev(i)
        endif
        call leglin(cbuf,zlevlab,i)
      enddo
      call legend(cbuf,6)  ! write buffer to legend and give position (6=LR)

      ! Add box below plot with run info
      ! Volcano:     Erup.start:
      ! Run date:    Plm Height:
      ! Windfile:    Duration:
      !              Volume:

      write(cstr_volcname,'(a10,a20)')'Volcano:  ' ,VolcanoName
      write(cstr_run_date,'(a10,a20)')'Run Date: ',cdf_time_log
      read(cdf_b3l1,*,iostat=ioerr) iw,iwf
      write(cstr_windfile,'(a10,i5)')'Windfile: ',iwf

      !e_StartTime,e_PlumeHeight,e_Duration,e_Volume
      write(cstr_ErStartT,'(a20,a20)')'Erup. Start Time:   ',HS_xmltime(e_StartTime(1),BaseYear,useLeap)
      write(cstr_ErHeight,'(a20,f4.1,a3)')'Erup. Plume Height: ',e_PlumeHeight(1),' km'
      write(cstr_ErDuratn,'(a20,f4.1,a6)')'Erup. Duration:     ',e_Duration(1),' hours'
      write(cstr_ErVomune,'(a20,f8.5,a10)')'Erup. Volume:       ',e_Volume(1),' km3 (DRE)'

      call messag(cstr_volcname,400 ,y_footer)
      call messag(cstr_run_date,400 ,y_footer+40)
      call messag(cstr_windfile,400 ,y_footer+80)
      if(neruptions.gt.1)then
        write(cstr_note,'(a45)')'WARNING: Multiple eruptions, only first given'
      call messag(cstr_note,400 ,y_footer+120)
      endif
      call messag(cstr_ErStartT,1200,y_footer)
      call messag(cstr_ErHeight,1200,y_footer+40)
      call messag(cstr_ErDuratn,1200,y_footer+80)
      call messag(cstr_ErVomune,1200,y_footer+120)

      !  Dislin Level 0:  before initialization or after termination
      call disfin()


      end subroutine write_2D_PNG_dislin

!##############################################################################
!
!    write_DepPOI_TS_PNG_gnuplot
!
!    if timestep = -1, then use the last step in file
!##############################################################################

      subroutine write_DepPOI_TS_PNG_gnuplot(plt_indx)

      use precis_param

      use Airports,      only : &
         nairports,Airport_Code,Airport_Name,Airport_x,Airport_y,&
         Airport_Latitude,Airport_Longitude,Airport_Thickness_TS

      use io_data,       only : &
         nWriteTimes,WriteTimes,VolcanoName

      use time_data,     only : &
         Simtime_in_hours

      implicit none

      integer :: plt_indx,i

      real(kind=8) :: ymaxpl
      character(len=14) :: dp_gnufile
      character(len=14) :: dp_outfile
      character(len=14) :: dp_pngfile
      character(len=25) :: gnucom

      write(dp_outfile,53) plt_indx,".dat"
      write(dp_gnufile,53) plt_indx,".gnu"
      write(dp_pngfile,54) plt_indx,".png"
 53   format('depTS_',i4.4,a4)
 54   format('gnupl_',i4.4,a4)

      open(54,file=dp_outfile,status='replace')
      do i = 1,nWriteTimes
        write(54,*)WriteTimes(i),Airport_Thickness_TS(plt_indx,i)
      enddo
      close(54)

      if(Airport_Thickness_TS(plt_indx,nWriteTimes).lt.0.01)then
        ymaxpl = 1.0
      elseif(Airport_Thickness_TS(plt_indx,nWriteTimes).lt.1.0)then
        ymaxpl = 1.0
      elseif(Airport_Thickness_TS(plt_indx,nWriteTimes).lt.5.0)then
        ymaxpl = 5.0
      elseif(Airport_Thickness_TS(plt_indx,nWriteTimes).lt.25.0)then
        ymaxpl = 25.0
      else
        ymaxpl = 100.0
      endif

      ! Set up to plot via gnuplot script
      open(55,file=dp_gnufile,status='replace')
      write(55,*)"set terminal png size 400,300"
      write(55,*)"set key bmargin left horizontal Right noreverse enhanced ",&
                 "autotitles box linetype -1 linewidth 1.000"
      write(55,*)"set border 31 lw 2.0 lc rgb '#000000'"
      write(55,*)"set style line 1 linecolor rgbcolor '#888888' linewidth 2.0 pt 7"
      write(55,*)"set ylabel 'Deposit Thickeness (mm)'"
      write(55,*)"set xlabel 'Time (hours after eruption)'"
      write(55,*)"set nokey"
      write(55,*)"set output '",dp_pngfile,"'"
      write(55,*)"set title '",Airport_Name,"'"
      write(55,*)"plot [0:",ceiling(Simtime_in_hours),"][0:",&
                 nint(ymaxpl),"] '",dp_outfile,"' with filledcurve x1 ls 1"
      close(55)

      write(gnucom,'(a11,a14)')'gnuplot -p ',dp_gnufile
      call execute_command_line(gnucom)

      end subroutine write_DepPOI_TS_PNG_gnuplot


!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      ! PLPLOT block
!      write(dp_pngfile2,55) plt_indx,".png"
! 55   format('plplt_',i4.4,a4)
!
!      xmin=real(0,kind=plflt)
!      xmax=real(ceiling(Simtime_in_hours),kind=plflt)
!      ymin=real(0,kind=plflt)
!      ymax=real(ymaxpl,kind=plflt)
!
!      allocate(x(nWriteTimes))
!      allocate(y(nWriteTimes))
!      allocate(x0(nWriteTimes+1))
!      allocate(y0(nWriteTimes+1))
!
!      x = real(WriteTimes,kind=plflt)
!      y = real(Airport_Thickness_TS,kind=plflt)
!      x0(1:nWriteTimes)=x(1:nWriteTimes)
!      x0(nWriteTimes+1)=x(nWriteTimes)
!      y0(1:nWriteTimes)=y(1:nWriteTimes)
!      y0(nWriteTimes+1)=0.0_plflt
!
!      ! Set up for plplot
!      call plsdev("pngcairo")      ! Set output device (png, pdf, etc.)
!      call plsfnam ( dp_pngfile2 ) ! Set output filename
!
!      call plsetopt("geometry","400x300, 400x300")  ! Set image size
!      call plsetopt("bg","FFFFFF")                  ! Set background color to white
!      ! Initialize plplot
!      call plinit()
!
!      ! Create a labelled box to hold the plot.
!      call plscolbg(r1,g1,b1) ! Set color back to black
!      call plcol0(0)
!      call plschr(0.0_plflt,1.7_plflt)  ! Chage font scale
!      call plenv( xmin, xmax, ymin, ymax, 0, 0 )
!      call pllab( "Time (hours after eruption)", "Deposit Thickeness (mm)", Airport_Name )
!
!      call plscolbg(r2,g2,b2) ! Set pen color to grey
!      call plcol0(0)
!      call plline( x , y ) ! Simple line plot
!      call plfill( x0(1:nWriteTimes+1), y0(1:nWriteTimes+1) )
!
!      ! Close PLplot library
!      call plend
!      ! end PLPLOT block
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      ! DISLIN block
!      ! https://www.dislin.de/
!      ! wget
!      ! https://ftp.gwdg.de/pub/grafik/dislin/linux/i586_64/dislin-11.4.linux.i586_64.tar.gz
!      !  Dislin Level 0:  before initialization or after termination
!      call metafl(CFMT)   ! set output driver/file-format (PNG); this is a 4-char string
!      call setpag('USAL') ! Set pagesize to US A Landscape (2790 x 2160)
!      call setfil('temp.png') ! Set output filename
!      call scrmod('REVERSE')  ! Default background is black; reverse to white
!
!      !  Dislin Level 1:  after initialization or a call to ENDGRF
!      call disini()       ! initialize plot (set to level 1)
!        ! setting of plot parameters
!      call pagera()       ! plot a border around the page
!      call triplx()  ! set font to triple stroke
!      call axspos(450,1800)  ! determine the position of the axis system
!      call axslen(2200,1200) ! defines the size of the axis system
!      call name('Time (hours after eruption)','X') ! Set x-axis title
!      call name('Deposit Thickeness (mm)','Y') ! Set y-axis title
!      call labdig(-1,'X') ! set number of decimal places for x label (-1 means no decimal)
!      call ticks(10,'XY') ! set number of ticks between labels
!      call titlin(Airport_Name,4)  ! Set the title to the airport name (4 is the bottom line)
!
!      !  Dislin Level 2: after a call to GRAF, GRAFP or GRAFMP
!        ! Now create graph and set to level 2
!      call graf(real(xmin,kind=4), real(xmax,kind=4), 0.0_4, 5.0_4, &
!                real(ymin,kind=4), real(ymax,kind=4), 0.0_4, 1.0_4)
!      call title() ! Actually write the title to the file
!      call setrgb(0.5_4, 0.5_4, 0.5_4)
!      !call curve(real(x,kind=4),real(y,kind=4),nWriteTimes)  ! This draws the line
!      call shdpat(16)  ! set shading pattern 
!      call shdcrv(real(x,kind=4),real(y,kind=4),nWriteTimes,& ! This fills below curve
!                  real(x,kind=4),real(0.0*y,kind=4),nWriteTimes)
!      call color('FORE') ! Reset color to defaul foreground color
!
!      !  Dislin Level 0:  before initialization or after termination
!      call disfin()


!      end subroutine write_DepPOI_TS_PNG

!##############################################################################
!
!    write_2D_ShapeFile
!
!    if timestep = -1, then use the last step in file
!##############################################################################

!      subroutine write_2D_ShapeFile

       !The call is:    CALL CONPTS (XRAY, N, YRAY, M, ZMAT, ZLEV, XPTRAY, YPTRAY,
       !MAXPTS, IRAY, MAXCRV, NCURVS)   level 0, 1, 2, 3
       !or:     void conpts (const float *xray, int n, const float *yray, int m, const
       !float *zmat, float zlev, float *xptray, float *yptray, int maxpts, int *iray,
       !int maxcrv, int *ncurvs);
       !
       !XRAY    is an array containing X-coordinates.
       !N       is the dimension of XRAY.
       !YRAY    is an array containing Y-coordinates.
       !M       is the dimension of YRAY.
       !ZMAT    is a matrix of the dimension (N, M) containing function values.
       !ZLEV    is a function value that defines the contour line to be calculated.
       !XPTRAY, YPTRAY  are returned arrays containing the calculated contour. The
       !arrays can contain several curves.
       !MAXPTS  is the maximal number of points that can be passed to XPTRAY and
       !YPTRAY.
       !IRAY    is a returned integer array that contains the number of points for each
       !generated contour curve.
       !MAXCRV  is the maximal number of entries that can be passed to IRAY.
       !NCURVS  is the returned number of generated curves.

!       call conpts(real(lon_cc_pd,kind=4),nxmax,&
!                 real(lat_cc_pd,kind=4),nymax,&
!                 real(DepositThickness,kind=4),0.1_4,&
!                 XPTS,YPTS,MAXPTS, &
!                 IRAY,MAXCRV,NCURVS)
       !do i=1,IRAY(1)
       !  write(*,*)XPTS(i),YPTS(i)
       !enddo

!      end subroutine write_2D_ShapeFile


