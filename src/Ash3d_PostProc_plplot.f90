!##############################################################################
!
! Ash3d_PostProc_plplot module
!
! This module provides the subroutines that use plplot for creating 2d maps,
! 2d vertical profiles, and the little deposit accumulation plots linked to
! the airport arrival kml (ash_arrivaltimes_airports.kml). The plplot
! library is linked at compile-time.
!
!      subroutine write_2Dmap_PNG_plplot
!      subroutine write_2Dprof_PNG_plplot
!      subroutine write_DepPOI_TS_PNG_plplot
!
!##############################################################################

      module Ash3d_PostProc_plplot

      use precis_param

      use io_units

!      use global_param,  only : &
!         DirDelim
!
!      use io_data,       only : &
!         Ash3dHome

      use plplot
      use iso_c_binding, only: c_ptr, c_loc, c_f_pointer

      implicit none

        ! Set everything to private by default
      private

        ! Publicly available subroutines/functions
      public write_2Dmap_PNG_plplot,    &
             write_2Dprof_PNG_plplot

        ! Publicly available variables

      integer :: lib_ver_major = 5
      !integer :: lib_ver_minor = 10
      integer :: lib_ver_minor = 14

!      character(100) :: Instit_IconFile

      contains
      !------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  write_2Dmap_PNG_plplot
!
!  Called from: Ash3d_PostProc.f90
!  Arguments:
!    nx            = x length of output array OutVar
!    ny            = y length of output array OutVar
!    iprod         = product ID
!    itime         = time index from netcdf data file
!    OutVar        = 2-d array to be written to ASCII file
!    writeContours = logical
!
!  This subroutine creates a png map of the variable in OutVar using the plplot
!  graphics package.  Annotations and contour levels are indicated via the
!  product ID (iprod).  If writeContours is set to true, then this subroutine
!  is only used for generating and storing the contours (for plplot, this
!  throws an error since contour data is not available) with no png written.
!  If timestep = -1, then use the last step in file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine write_2Dmap_PNG_plplot(nx,ny,iprod,itime,OutVar,writeContours)

      use mesh,          only : &
         x_cc_pd,y_cc_pd,lon_cc_pd,lat_cc_pd, &
         IsLatLon

      use Output_Vars,   only : &
         ContourFilled,Con_Cust,Con_Cust_N,Con_Cust_RGB,Con_Cust_Lev,&
         Con_DepThick_mm_N,Con_DepThick_mm_Lev,Con_DepThick_mm_RGB, &
         Con_DepThick_in_N,Con_DepThick_in_Lev,Con_DepThick_in_RGB, &
         Con_DepTime_N,Con_DepTime_Lev,Con_DepTime_RGB, &
         Con_CloudCon_N,Con_CloudCon_Lev,Con_CloudCon_RGB, &
         Con_CloudTop_N,Con_CloudTop_RGB,Con_CloudTop_Lev, &
         Con_CloudBot_N,Con_CloudBot_RGB,Con_CloudBot_Lev, &
         Con_CloudLoad_N,Con_CloudLoad_RGB,Con_CloudLoad_Lev, &
         Con_CloudRef_N,Con_CloudRef_RGB,Con_CloudRef_Lev, &
         Con_CloudTime_N,Con_CloudTime_RGB,Con_CloudTime_Lev, &
         ContourLev,nConLev
       ! These are needed if we can figure out how to extract contour data
!         ContourDataX,ContourDataY,ContourDataNcurves,ContourDataNpoints,&
!         CONTOUR_MAXCURVES,CONTOUR_MAXPOINTS

      use io_data,       only : &
         WriteTimes,cdf_b3l1,VolcanoName

      use Source,        only : &
         e_Volume,e_Duration,e_StartTime,e_PlumeHeight,lon_volcano,lat_volcano

      use time_data,     only : &
         os_time_log,SimStartHour,BaseYear,useLeap

      use citywriter

      integer      ,intent(in) :: nx
      integer      ,intent(in) :: ny
      integer      ,intent(in) :: iprod
      integer      ,intent(in) :: itime
      real(kind=ip),intent(in) :: OutVar(nx,ny)
      logical      ,intent(in) :: writeContours

      integer :: i
      integer      ,dimension(:,:),allocatable :: zrgb
      character(len=40) :: title_plot
      character(len=15) :: title_legend
      character(len=40) :: outfile_name
      character(len= 9) :: cio
      character(len= 4) :: outfile_ext = '.png'
      character(len=10) :: units
      character(len=80) :: outstring
      integer :: ioerr,iw,iwf

      ! PLPLOT variables
      real(kind=plflt)  :: xmin
      real(kind=plflt)  :: xmax
      real(kind=plflt)  :: ymin
      real(kind=plflt)  :: ymax
      real(kind=plflt)  :: vmin
      real(kind=plflt)  :: vmax
      real(kind=plflt), dimension(:),   allocatable :: x, y
      real(kind=plflt), dimension(:,:), allocatable :: var
      real(kind=plflt)   :: tr(6)
      real(kind=plflt)   :: clevel(1)
      integer(kind=4):: opt

      !integer, parameter :: MAX_NLEGEND = 11       ! max number of legend entries
      !integer(kind=4)    :: opt_array(MAX_NLEGEND)
      !integer(kind=4)    :: text_colors(MAX_NLEGEND)
      !integer(kind=4)    :: box_colors(MAX_NLEGEND)
      !integer(kind=4)    :: box_patterns(MAX_NLEGEND)
      !real(kind=plflt)   :: box_scales(MAX_NLEGEND)
      !real(kind=plflt)   :: box_line_widths(MAX_NLEGEND)
      !integer(kind=4)    :: line_colors(MAX_NLEGEND)
      !integer(kind=4)    :: line_styles(MAX_NLEGEND)
      !real(kind=plflt)   :: line_widths(MAX_NLEGEND)
      !integer(kind=4)    :: symbol_numbers(MAX_NLEGEND)
      !integer(kind=4)    :: symbol_colors(MAX_NLEGEND)
      !real(kind=plflt)   :: symbol_scales(MAX_NLEGEND)
      !character(len=200) :: text(MAX_NLEGEND)
      !character(len=3)   :: symbols(MAX_NLEGEND)

      integer(kind=4)   ,dimension(:),allocatable :: opt_array
      integer(kind=4)   ,dimension(:),allocatable :: text_colors
      integer(kind=4)   ,dimension(:),allocatable :: box_colors
      integer(kind=4)   ,dimension(:),allocatable :: box_patterns
      integer(kind=4)   ,dimension(:),allocatable :: line_colors
      integer(kind=4)   ,dimension(:),allocatable :: line_styles
      integer(kind=4)   ,dimension(:),allocatable :: symbol_numbers
      integer(kind=4)   ,dimension(:),allocatable :: symbol_colors
      integer(kind=4)   ,dimension(:),allocatable :: red
      integer(kind=4)   ,dimension(:),allocatable :: green
      integer(kind=4)   ,dimension(:),allocatable :: blue
      character(len=200),dimension(:),allocatable :: text
      character(len=3)  ,dimension(:),allocatable :: symbols
      real(kind=plflt)  ,dimension(:),allocatable :: box_scales
      real(kind=plflt)  ,dimension(:),allocatable :: box_line_widths
      real(kind=plflt)  ,dimension(:),allocatable :: line_widths
      real(kind=plflt)  ,dimension(:),allocatable :: symbol_scales
      real(kind=plflt)  ,dimension(:),allocatable :: alpha

      real(kind=plflt)   :: legend_width, legend_height
      real(kind=plflt)   :: x_offset,y_offset,plot_width
      real(kind=plflt)   :: text_offset
      real(kind=plflt)   :: text_scale
      real(kind=plflt)   :: text_spacing
      real(kind=plflt)   :: text_justification

      integer(kind=4):: nrow, ncolumn
      integer(kind=4):: bg_color,bb_color,bb_style
      integer(kind=4):: pos_opt

      real(kind=plflt)  :: dy_newline
      integer :: plsetopt_rc

      integer :: ncities
      real(kind=ip),dimension(:),allocatable     :: lon_cities
      real(kind=ip),dimension(:),allocatable     :: lat_cities
      character(len=26),dimension(:),allocatable :: name_cities

      INTERFACE
        character (len=20) function HS_xmltime(HoursSince,byear,useLeaps)
          real(kind=8)              :: HoursSince
          integer                   :: byear
          logical                   :: useLeaps
        end function HS_xmltime
      END INTERFACE

      if(writeContours)then
        write(errlog(io),*)"Running plplot to calculate contours lines"
        write(errlog(io),*)"Not sure yet how to save contour data with plplot"
        write(errlog(io),*)"If you want shapefiles, recompile without plplot or"
        write(errlog(io),*)" reset the plot_pref_shp variable."
        write(errlog(io),*)"Exiting"
        stop 1
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Running plplot to generate contour plot"
        endif;enddo
      endif

      ncities = 20
      allocate(lon_cities(ncities))
      allocate(lat_cities(ncities))
      allocate(name_cities(ncities))

      if(iprod.eq.5.or.iprod.eq.6)then
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

      if(Con_Cust)then
        nConLev = Con_Cust_N
        allocate(zrgb(nConLev,3))
        allocate(ContourLev(nConLev))
        ContourLev(1:nConLev) = Con_Cust_Lev(1:nConLev)
        zrgb(1:nConLev,1:3) = Con_Cust_RGB(1:nConLev,1:3)
      endif

      if(iprod.eq.3)then       ! deposit at specified times (mm)
        write(outfile_name,'(a15,a9,a4)')'Ash3d_Deposit_t',cio,outfile_ext
        write(title_plot,'(a20,f5.2,a6)')'Deposit Thickness t=',WriteTimes(itime),' hours'
        title_legend = 'Dep.Thick.(mm)'
        units = " (mm)"
        if(.not.Con_Cust)then
          nConLev = Con_DepThick_mm_N
          allocate(zrgb(nConLev,3))
          allocate(ContourLev(nConLev))
          ContourLev(1:nConLev) = Con_DepThick_mm_Lev(1:nConLev)
          zrgb(1:nConLev,1:3) = Con_DepThick_mm_RGB(1:nConLev,1:3)
        endif
      elseif(iprod.eq.4)then   ! deposit at specified times (inches)
        write(outfile_name,'(a15,a9,a4)')'Ash3d_Deposit_t',cio,outfile_ext
        write(title_plot,'(a20,f5.2,a6)')'Deposit Thickness t=',WriteTimes(itime),' hours'
        title_legend = 'Dep.Thick.(in)'
        units = " (in)"
        if(.not.Con_Cust)then
          nConLev = Con_DepThick_in_N
          allocate(zrgb(nConLev,3))
          allocate(ContourLev(nConLev))
          ContourLev(1:nConLev) = Con_DepThick_in_Lev(1:nConLev)
          zrgb(1:nConLev,1:3) = Con_DepThick_in_RGB(1:nConLev,1:3)
        endif
      elseif(iprod.eq.5)then       ! deposit at final time (mm)
        write(outfile_name,'(a13,a9,a4)')'Ash3d_Deposit',cio,outfile_ext
        title_plot = 'Final Deposit Thickness'
        title_legend = 'Dep.Thick.(mm)'
        units = " (mm)"
        if(.not.Con_Cust)then
          nConLev = Con_DepThick_mm_N
          allocate(zrgb(nConLev,3))
          allocate(ContourLev(nConLev))
          ContourLev(1:nConLev) = Con_DepThick_mm_Lev(1:nConLev)
          zrgb(1:nConLev,1:3) = Con_DepThick_mm_RGB(1:nConLev,1:3)
        endif
      elseif(iprod.eq.6)then   ! deposit at final time (inches)
        write(outfile_name,'(a13,a9,a4)')'Ash3d_Deposit',cio,outfile_ext
        title_plot = 'Final Deposit Thickness'
        title_legend = 'Dep.Thick.(in)'
        units = " (in)"
        if(.not.Con_Cust)then
          nConLev = Con_DepThick_in_N
          allocate(zrgb(nConLev,3))
          allocate(ContourLev(nConLev))
          ContourLev(1:nConLev) = Con_DepThick_in_Lev(1:nConLev)
          zrgb(1:nConLev,1:3) = Con_DepThick_in_RGB(1:nConLev,1:3)
        endif
      elseif(iprod.eq.7)then   ! ashfall arrival time (hours)
        write(outfile_name,'(a22)')'DepositArrivalTime.png'
        write(title_plot,'(a20)')'Ashfall arrival time'
        title_legend = 'Time (hours)'
        units = " (hours)"
        if(.not.Con_Cust)then
          nConLev = Con_DepTime_N
          allocate(zrgb(nConLev,3))
          allocate(ContourLev(nConLev))
          ContourLev(1:nConLev) = Con_DepTime_Lev(1:nConLev)
          zrgb(1:nConLev,1:3) = Con_DepTime_RGB(1:nConLev,1:3)
        endif
      elseif(iprod.eq.8)then   ! ashfall arrival at airports/POI (mm)
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: No map PNG output option for airport arrival time data."
          write(errlog(io),*)"       Should not be in write_2Dmap_PNG_dislin"
        endif;enddo
        stop 1
      elseif(iprod.eq.9)then   ! ash-cloud concentration
        write(outfile_name,'(a16,a9,a4)')'Ash3d_CloudCon_t',cio,outfile_ext
        write(title_plot,'(a26,f5.2,a6)')'Ash-cloud concentration t=',WriteTimes(itime),' hours'
        title_legend = 'Max.Con.(mg/m3)'
        units = " (mg/m3)"
        if(.not.Con_Cust)then
          nConLev = Con_CloudCon_N
          allocate(zrgb(nConLev,3))
          allocate(ContourLev(nConLev))
          ContourLev(1:nConLev) = Con_CloudCon_Lev(1:nConLev)
          zrgb(1:nConLev,1:3) = Con_CloudCon_RGB(1:nConLev,1:3)
        endif
      elseif(iprod.eq.10)then   ! ash-cloud height
        write(outfile_name,'(a19,a9,a4)')'Ash3d_CloudHeight_t',cio,outfile_ext
        write(title_plot,'(a19,f5.2,a6)')'Ash-cloud height t=',WriteTimes(itime),' hours'
        title_legend = 'Cld.Height(km)'
        units = " (km)"
        if(.not.Con_Cust)then
          nConLev = Con_CloudTop_N
          allocate(zrgb(nConLev,3))
          allocate(ContourLev(nConLev))
          ContourLev(1:nConLev) = Con_CloudTop_Lev(1:nConLev)
          zrgb(1:nConLev,1:3) = Con_CloudTop_RGB(1:nConLev,1:3)
        endif
      elseif(iprod.eq.11)then   ! ash-cloud bottom
        write(outfile_name,'(a16,a9,a4)')'Ash3d_CloudBot_t',cio,outfile_ext
        write(title_plot,'(a19,f5.2,a6)')'Ash-cloud bottom t=',WriteTimes(itime),' hours'
        title_legend = 'Cld.Bot.(km)'
        units = " (km)"
        if(.not.Con_Cust)then
          nConLev = Con_CloudBot_N
          allocate(zrgb(nConLev,3))
          allocate(ContourLev(nConLev))
          ContourLev(1:nConLev) = Con_CloudBot_Lev(1:nConLev)
          zrgb(1:nConLev,1:3) = Con_CloudBot_RGB(1:nConLev,1:3)
        endif
      elseif(iprod.eq.12)then   ! ash-cloud load
        write(outfile_name,'(a17,a9,a4)')'Ash3d_CloudLoad_t',cio,outfile_ext
        write(title_plot,'(a17,f5.2,a6)')'Ash-cloud load t=',WriteTimes(itime),' hours'
        title_legend = 'Cld.Load(T/km2)'
        units = " (T/km2)"
        if(.not.Con_Cust)then
          nConLev = Con_CloudLoad_N
          allocate(zrgb(nConLev,3))
          allocate(ContourLev(nConLev))
          ContourLev(1:nConLev) = Con_CloudLoad_Lev(1:nConLev)
          zrgb(1:nConLev,1:3) = Con_CloudLoad_RGB(1:nConLev,1:3)
        endif
      elseif(iprod.eq.13)then  ! radar reflectivity
        write(outfile_name,'(a20,a9,a4)')'Ash3d_CloudRadRefl_t',cio,outfile_ext
        write(title_plot,'(a24,f5.2,a6)')'Ash-cloud radar refl. t=',WriteTimes(itime),' hours'
        title_legend = 'Cld.Refl.(dBz)'
        units = " (dBz)"
        if(.not.Con_Cust)then
          nConLev = Con_CloudRef_N
          allocate(zrgb(nConLev,3))
          allocate(ContourLev(nConLev))
          ContourLev(1:nConLev) = Con_CloudRef_Lev(1:nConLev)
          zrgb(1:nConLev,1:3) = Con_CloudRef_RGB(1:nConLev,1:3)
        endif
      elseif(iprod.eq.14)then   ! ashcloud arrival time (hours)
        write(outfile_name,'(a20)')'CloudArrivalTime.png'
        write(title_plot,'(a22)')'Ash-cloud arrival time'
        title_legend = 'Time (hours)'
        units = " (hours)"
        if(.not.Con_Cust)then
          nConLev = Con_CloudTime_N
          allocate(zrgb(nConLev,3))
          allocate(ContourLev(nConLev))
          ContourLev(1:nConLev) = Con_CloudTime_Lev(1:nConLev)
          zrgb(1:nConLev,1:3) = Con_CloudTime_RGB(1:nConLev,1:3)
        endif
      elseif(iprod.eq.15)then   ! topography
        write(outfile_name,'(a14)')'Topography.png'
        write(title_plot,'(a10)')'Topography'
        title_legend = 'Elevation (km)'
        units = " (hours)"
        if(.not.Con_Cust)then
          nConLev = 8
          allocate(zrgb(nConLev,3))
          allocate(ContourLev(nConLev))
          ContourLev = (/1.0_ip, 2.0_ip, 3.0_ip, 4.0_ip, &
                  5.0_ip, 6.0_ip, 7.0_ip, 8.0_ip/)
        endif
      elseif(iprod.eq.16)then   ! profile plots
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: No map PNG output option for vertical profile data."
          write(errlog(io),*)"       Should not be in write_2Dmap_PNG_dislin"
        endif;enddo
        stop 1
      else
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: unexpected variable"
        endif;enddo
        stop 1
      endif

      if(IsLatLon)then
        xmin = real(minval(lon_cc_pd(1:nx)),kind=plflt)
        xmax = real(maxval(lon_cc_pd(1:nx)),kind=plflt)
        ymin = real(minval(lat_cc_pd(1:ny)),kind=plflt)
        ymax = real(maxval(lat_cc_pd(1:ny)),kind=plflt)
      else
        xmin = real(minval(x_cc_pd(1:nx)),kind=plflt)
        xmax = real(maxval(x_cc_pd(1:nx)),kind=plflt)
        ymin = real(minval(y_cc_pd(1:ny)),kind=plflt)
        ymax = real(maxval(y_cc_pd(1:ny)),kind=plflt)
      endif
      call citylist(0,real(xmin,kind=ip),real(xmax,kind=ip),&
                      real(ymin,kind=ip),real(ymax,kind=ip),&
                      ncities,                            &
                      lon_cities, &
                      lat_cities, &
                      name_cities)
      allocate(x(nx))
      allocate(y(ny))
      allocate(var(nx,ny))
      x(1:nx) = real(lon_cc_pd(1:nx),kind=plflt)
      y(1:ny) = real(lat_cc_pd(1:ny),kind=plflt)
      var(1:nx,1:ny) = real(OutVar(1:nx,1:ny),kind=plflt)
      vmin=real(minval(var(:,:)),kind=plflt)
      vmax=real(maxval(var(:,:)),kind=plflt)

      tr = (/ (xmax-xmin)/real(nx-1,kind=plflt), 0.0_plflt, xmin, &
              0.0_plflt, (ymax-ymin)/real(ny-1,kind=plflt), ymin /)

      ! Set up for plplot
      call plsdev("pngcairo")      ! Set output device (png, pdf, etc.)
      call plsfnam (outfile_name)  ! Set output filename

      ! set image size via command-line options tool plsetopt
      if(lib_ver_minor.lt.12)then
        ! prior to v5.12, this was a subroutine call
!        call plsetopt("geometry","854x603, 854x603")
      else
        ! after v5.12, this is a function call returing an error code plsetopt_rc
        plsetopt_rc = plsetopt("geometry","854x603, 854x603")
      endif

      !-----------------------------------
      ! sets cmap0 palette via pal file for discrete elements
      call plspal0('cmap0_black_on_white.pal')
      ! sets cmap1 palette via pal file for continuous elements
      if(lib_ver_minor.lt.12)then
        ! prior to v5.12, the second argument was an integer
!        call plspal1('cmap1_blue_yellow.pal',1)
      else
        ! after v5.12, this needs to be a logical
        call plspal1('cmap1_blue_yellow.pal',.true.)
      endif

      call plscmap0n(16)           ! sets number of colors in cmap0

      ! Initialize plplot
      call plinit()

        ! pladv: Advance the (sub-)page
      call pladv(0)
        ! plvpor: Specify viewport using normalized subpage coordinates
      call plvpor(0.1_plflt, 0.8_plflt, 0.3_plflt, 0.8_plflt)
        ! plwind: Specify window 
      call plwind(xmin,xmax,ymin,ymax)
      call plbox('bcnst', 0.0_plflt, 0, 'bcnstv', 0.0_plflt, 0)

      call plcol0(1)
      call plmap('usaglobe', xmin, xmax, ymin, ymax)

      ! Add cities
      do i=1,ncities
        if(lon_cities(i).lt.xmin)lon_cities(i)=lon_cities(i)+360.0_ip
          ! plssym: Set symbol size : default, scale
        call plssym( 0.0_plflt, 2.0_plflt )
          ! plpoin: Plot a glyph at the specified points 
        call plpoin(real(lon_cities(i:i),kind=plflt),&
                      real(lat_cities(i:i),kind=plflt),&
                      17) ! code 17 is a black dot
          ! plssym: Set symbol size : default, scale
        call plssym( 0.0_plflt, 0.5_plflt )
        call plschr( 0.0_plflt, 0.7_plflt )

          ! plptex : Write text inside the viewport (x,y,dx,dy,just,strin)
        call plptex( real(lon_cities(i)+1.0,kind=plflt),real(lat_cities(i),kind=plflt), &
             0.0_plflt, 0.0_plflt, 0.0_plflt, &
             adjustl(trim(name_cities(i))))
      enddo
      call plschr( 0.0_plflt, 1.0_plflt )

      call plssym( 0.0_plflt, 1.0_plflt )
        ! plpoin: Plot a glyph at the specified points 
      if(lon_volcano.lt.xmin)then
        lon_cities(1)=lon_volcano+360.0_ip
      else
        lon_cities(1)=lon_volcano
      endif
      lat_cities(1)=lat_volcano
      i=1
      call plpoin(real(lon_cities(i:i),kind=plflt),&
                    real(lat_cities(i:i),kind=plflt),&
                    7) ! code 7 is a triangle
                       ! (https://plplot.sourceforge.net/examples.php?demo=06&lbind=Fortran)

      do i=1,nConLev
        call plcol1(real(dble(i)/dble(nConLev),kind=plflt))
        clevel(1) = real(ContourLev(i),kind=plflt)
        call plcont(var,1,nx,1,ny,clevel, tr)
      enddo
      call pllab("Longitude", "Latitude", title_plot)
      if(lib_ver_minor.lt.12)then
        ! prior to v5.12, the second argument was an integer
!        call plstransform( 0 )
      else
        call plstransform
      endif

      !---------------------------------------------------------

      !call plcol0(2)
      ! Set the color we will use for the legend background (index 15)
      call plscol0a( 15, 255, 255, 255, 1.0_plflt )
      allocate(opt_array(nConLev))
      allocate(text_colors(nConLev))
      allocate(box_colors(nConLev))
      allocate(box_patterns(nConLev))
      allocate(line_colors(nConLev))
      allocate(line_styles(nConLev))
      allocate(symbol_numbers(nConLev))
      allocate(symbol_colors(nConLev))
      allocate(text(nConLev))
      allocate(symbols(nConLev))
      allocate(box_scales(nConLev))
      allocate(box_line_widths(nConLev))
      allocate(line_widths(nConLev))
      allocate(symbol_scales(nConLev))
      allocate(red(nConLev))
      allocate(green(nConLev))
      allocate(blue(nConLev))
      allocate(alpha(nConLev))
      do i=1,nConLev
        pos_opt = PL_POSITION_RIGHT + PL_POSITION_OUTSIDE
        opt = PL_LEGEND_BACKGROUND + PL_LEGEND_BOUNDING_BOX
        text_colors(i)   = 1 + mod( i-1, nConLev )
        line_colors(i)   = 1 + mod( i-1, nConLev )
        !call plcol1(real(dble(i)/dble(nConLev),kind=plflt))
        red(i)   = i
        green(i) = i
        blue(i)  = i
        alpha(i) = 1.0_plflt

        line_styles(i)   = 1
        line_widths(i)   = 1
        symbol_colors(i) = 1 + mod( i-1, nConLev )
        box_colors(i)     = 2
        box_patterns(i)   = 3
        box_scales(i)     = 0.8_plflt
        box_line_widths(i)= 1
        if(abs(ContourLev(i)).lt.0.01_ip.or.abs(ContourLev(i)).ge.1000.0_ip)then
          write( text(i), '(e7.2)' ) real(ContourLev(i),kind=4)
        else
          write( text(i), '(f7.2)' ) real(ContourLev(i),kind=4)
        endif
        x_offset       = 0.05_plflt
        y_offset       = 0.0_plflt
        plot_width     = 0.05_plflt
        bg_color       = 15
        bb_color       = 1
        bb_style       = 1
        nrow           = nConLev
        ncolumn        = 1    ! Note: nlegend=nrow * ncolumn
        opt_array(i)   = PL_LEGEND_LINE
        text_offset        = 1.0_plflt
        text_scale         = 0.75_plflt
        text_spacing       = 1.5_plflt
        text_justification = 0.0_plflt
        symbols(i)        = '*'
      enddo
      ! Now set the RGB values from above tothe cmap0
      !call plscmap0a(red, green, blue, alpha)
        ! pladv: Advance the (sub-)page
!      call pladv(1)
!        ! plvpor: Specify viewport using normalized subpage coordinates
!      call plvpor(0.75_plflt, 1.0_plflt, 0.6_plflt, 0.7_plflt)
!        ! plwind: Specify window 
!      call plwind(0.0_plflt, 1.0_plflt, 0.0_plflt, 1.0_plflt )
!      call plschr( 0.0_plflt, 0.7_plflt )

      call pllegend(    &
          legend_width, &  ! these are output vars
          legend_height,&  ! these are output vars
          opt,         & ! int: controls overall legend
          pos_opt,     & ! int: controls legend position
          x_offset,    & ! flt: legend offset
          y_offset,    & ! flt: legend offset
          plot_width,  & ! flt: horz width
          bg_color,    & ! int: background color from cmap0
          bb_color,    & ! int: bounding box color from cmap0
          bb_style,    & ! int: bounding box line style
          nrow,ncolumn,& ! int: rows and columns of legend
          opt_array(1:nConLev), & ! int vec: 
          text_offset, & ! flt: Offset of the text area from the plot
          text_scale,  & ! flt: Character height scale
          text_spacing, &!  flt: Vertical spacing in units of the character height
          text_justification, & ! flt: 0., 0.5, or 1.  for L, C, R
          text_colors(1:nConLev), &
          text(1:nConLev),                           &
          box_colors(1:nConLev), &
          box_patterns(1:nConLev), &
          box_scales(1:nConLev), &
          box_line_widths(1:nConLev), &
          line_colors(1:nConLev), &
          line_styles(1:nConLev), &
          line_widths(1:nConLev),                 &
          symbol_colors(1:nConLev), &
          symbol_scales(1:nConLev), &
          symbol_numbers(1:nConLev),& 
          symbols  )

      ! Now add the annotation box
        ! First the title of the legend
        ! plvpor: Specify viewport using normalized subpage coordinates
      call plvpor(0.8_plflt, 1.0_plflt, 0.7_plflt, 0.8_plflt)
        ! plwind: Specify window 
      call plwind(0.0_plflt, 1.0_plflt, 0.0_plflt, 1.0_plflt )
      write(outstring,*)trim(adjustl(title_legend))
      call plschr( 0.0_plflt, 0.8_plflt )
      call plptex(0.1_plflt, 0.5_plflt,& ! x,y
                  1.0_plflt, 0.0_plflt, & ! dx,dy
                  0.0_plflt,            & ! just
                  outstring )             ! text
      ! And the boxes below
        ! pladv: Advance the (sub-)page
      call pladv(1)
        ! plvpor: Specify viewport using normalized subpage coordinates
      call plvpor(0.05_plflt, 0.35_plflt, 0.05_plflt, 0.2_plflt)
        ! plwind: Specify window 
      call plwind(0.0_plflt, 1.0_plflt, 0.0_plflt, 1.0_plflt )
        ! plbox: Draw a box with axes, etc
        !       (xopt, xtick, nxsub, yopt, ytick, nysub)
      !call plbox('bc', 0.0_plflt, 0, 'bc', 0.0_plflt, 0 )
        !  plschr: Set character size
        !          (def, scale)
      call plschr( 0.0_plflt, 0.7_plflt )
      dy_newline = 0.13_plflt
        ! plptex: Write text inside the viewport
      write(outstring,*)"Volcano: ",trim(adjustl(VolcanoName))
      call plptex(0.02_plflt, 1.0_plflt-1.0_plflt*dy_newline, 1.0_plflt, 0.0_plflt, &
                  0.0_plflt, outstring )
      write(outstring,*)"Run Date: ",trim(adjustl(os_time_log))
      call plptex(0.02_plflt, 1.0_plflt-2.0_plflt*dy_newline, 1.0_plflt, 0.0_plflt, &
                  0.0_plflt, outstring )
      read(cdf_b3l1,*,iostat=ioerr) iw,iwf
      write(outstring,*)"Windfile: ",iwf
      call plptex(0.02_plflt, 1.0_plflt-3.0_plflt*dy_newline, 1.0_plflt, 0.0_plflt, &
                  0.0_plflt, outstring )

        ! pladv: Advance the (sub-)page
      !call pladv(2)
        ! plvpor: Specify viewport using normalized subpage coordinates
      call plvpor(0.4_plflt, 0.75_plflt, 0.05_plflt, 0.2_plflt)
        ! plwind: Specify window 
      call plwind(0.0_plflt, 1.0_plflt, 0.0_plflt, 1.0_plflt )
        ! plbox: Draw a box with axes, etc
        !       (xopt, xtick, nxsub, yopt, ytick, nysub)
      !call plbox('bc', 0.0_plflt, 0, 'bc', 0.0_plflt, 0 )

      write(outstring,*)"Erup. Start Time: ",HS_xmltime(SimStartHour+e_StartTime(1),BaseYear,useLeap)
      call plptex(0.02_plflt, 1.0_plflt-1.0_plflt*dy_newline, 1.0_plflt, 0.0_plflt, &
                  0.0_plflt, outstring )
      write(outstring,111)e_PlumeHeight(1)
 111  format(' Erup. Plume Height: ',f7.2,' km')
      call plptex(0.02_plflt, 1.0_plflt-2.0_plflt*dy_newline, 1.0_plflt, 0.0_plflt, &
                  0.0_plflt, outstring )
      write(outstring,121)e_Duration(1)
 121  format(' Erup. Duration: ',f7.2,' hours')
      call plptex(0.02_plflt, 1.0_plflt-3.0_plflt*dy_newline, 1.0_plflt, 0.0_plflt, &
                  0.0_plflt, outstring )
      write(outstring,131)e_Volume(1)
 131  format(' Erup. Volume: ',f10.5,' km3 (DRE)')
      call plptex(0.02_plflt, 1.0_plflt-4.0_plflt*dy_newline, 1.0_plflt, 0.0_plflt, &
                  0.0_plflt, outstring )

      call plend()

      ! clean up memory
      if(allocated(lon_cities))      deallocate(lon_cities)
      if(allocated(lat_cities))      deallocate(lat_cities)
      if(allocated(name_cities))     deallocate(name_cities)
      if(allocated(zrgb))            deallocate(zrgb)
      if(allocated(x))               deallocate(x)
      if(allocated(y))               deallocate(y)
      if(allocated(var))             deallocate(var)
      if(allocated(opt_array))       deallocate(opt_array)
      if(allocated(text_colors))     deallocate(text_colors)
      if(allocated(box_colors))      deallocate(box_colors)
      if(allocated(box_patterns))    deallocate(box_patterns)
      if(allocated(line_colors))     deallocate(line_colors)
      if(allocated(line_styles))     deallocate(line_styles)
      if(allocated(symbol_numbers))  deallocate(symbol_numbers)
      if(allocated(symbol_colors))   deallocate(symbol_colors)
      if(allocated(text))            deallocate(text)
      if(allocated(symbols))         deallocate(symbols)
      if(allocated(box_scales))      deallocate(box_scales)
      if(allocated(box_line_widths)) deallocate(box_line_widths)
      if(allocated(line_widths))     deallocate(line_widths)
      if(allocated(symbol_scales))   deallocate(symbol_scales)
      if(allocated(red))             deallocate(red)
      if(allocated(green))           deallocate(green)
      if(allocated(blue))            deallocate(blue)
      if(allocated(alpha))           deallocate(alpha)

      end subroutine write_2Dmap_PNG_plplot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  write_2Dprof_PNG_plplot
!
!  Called from: Ash3d_PostProc.f90
!  Arguments:
!    vprof_ID        = ID of the profile (number in list from Ash3d control file)
!
!  This subroutine creates a png plot of the transient vertical profile of ash
!  concentration above the profile point (vprof_ID) using the plplot graphics
!  package.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine write_2Dprof_PNG_plplot(vprof_ID)

      use global_param,  only : &
         KG_2_MG,KM3_2_M3

      use mesh,          only : &
         nzmax,z_cc_pd

      use Output_Vars,   only : &
         pr_ash,CLOUDCON_THRESH,CONTOUR_MAXCURVES

      use time_data,     only : &
         ntmax,time_native

      use io_data,       only : &
         Site_vprofile,x_vprofile,y_vprofile,cdf_b3l1,VolcanoName

      use Source,        only : &
         e_Volume,e_Duration,e_StartTime,e_PlumeHeight

      use time_data,     only : &
         os_time_log,SimStartHour,BaseYear,useLeap

      integer,intent(in) :: vprof_ID

      !integer,parameter :: NLEVEL     = 11
      !integer,parameter :: NLEVEL     = 30
      integer,parameter :: NUM_AXES   = 1
      integer,parameter :: NUM_LABELS = 1

      character(len=14) :: dp_pngfile
      character(len=26) :: coord_str
      character(len=76) :: title_str
      character(len=80) :: outstring
      integer :: k,i
      integer :: ioerr,iw,iwf

      ! PLPLOT variables
      real(kind=plflt)  :: tmin
      real(kind=plflt)  :: tmax
      real(kind=plflt)  :: zmin
      real(kind=plflt)  :: zmax
      real(kind=plflt)  :: cmin
      real(kind=plflt)  :: cmax
      real(kind=plflt), dimension(:),   allocatable :: t, z
      real(kind=plflt), dimension(:,:), allocatable :: conc
      real(kind=plflt), dimension(:),   allocatable :: shedge
      integer           :: cont_color
      real(kind=plflt)  :: fill_width, cont_width
      real(kind=plflt)  :: colorbar_width, colorbar_height
      real(kind=plflt)  :: dy_newline
      !integer           :: NUM_AXES, NUM_LABELS
      !parameter(NUM_AXES=1, NUM_LABELS=1)
      character(len=20) :: axis_opts(NUM_AXES)
      integer           :: num_values(NUM_AXES)
      !real(kind=plflt)  :: values(NUM_AXES,NLEVEL+1)
      real(kind=plflt)  :: values(NUM_AXES,CONTOUR_MAXCURVES+1)
      real(kind=plflt)  :: axis_ticks(NUM_AXES)
      integer           :: axis_subticks(NUM_AXES)
      character(len=100):: labels(NUM_LABELS)
      integer           :: label_opts(NUM_LABELS)
      real(kind=plflt)  :: cloudcon_thresh_mgm3
      integer :: plsetopt_rc

      real(kind=plflt)   :: tr(6)
      real(kind=plflt)   :: clevel(1)
      character(len=1)  :: defined     ! if lib_ver_minor.lt.12

      INTERFACE
        character (len=20) function HS_xmltime(HoursSince,byear,useLeaps)
          real(kind=8)              :: HoursSince
          integer                   :: byear
          logical                   :: useLeaps
        end function HS_xmltime
      END INTERFACE

      cloudcon_thresh_mgm3 = CLOUDCON_THRESH * KG_2_MG / KM3_2_M3 !convert from kg/km3 to mg/m3

      clevel(1) = cloudcon_thresh_mgm3

      write(dp_pngfile,54) vprof_ID,".png"
 54   format('plplt_',i4.4,a4)

      tmin=real(0,kind=plflt)
      tmax=real(ceiling(time_native(ntmax)),kind=plflt)
      zmin=real(0,kind=plflt)
      zmax=real(z_cc_pd(nzmax),kind=plflt)
      cmin=real(0,kind=plflt)
      cmax=real(maxval(pr_ash(:,:,vprof_ID)),kind=plflt)    ! Get the max value for this profile
      cmax=real(max(cmax,cloudcon_thresh_mgm3),kind=plflt)  ! Do not let cmax drop below the threshold

      tr = (/ tmax/real(ntmax-1,kind=plflt), 0.0_plflt, 0.0_plflt, &
              0.0_plflt, zmax/real(nzmax-1,kind=plflt), 0.0_plflt /)

      allocate(t(ntmax))
      allocate(z(nzmax))
      allocate(conc(ntmax,nzmax))

      t = real(time_native(1:ntmax),kind=plflt)
      z = real(z_cc_pd(1:nzmax),kind=plflt)
      do i=1,ntmax
        do k=1,nzmax
          conc(i,k) = pr_ash(k,i,vprof_ID)
        enddo
      enddo

!      allocate(shedge(NLEVEL+1))
      allocate(shedge(CONTOUR_MAXCURVES+1))
      ! Here we linearly interpolate color levels to the min/max of the data
!      do i = 1, NLEVEL+1
!        shedge(i) = cmin + (cmax - cmin) * real(i-1,kind=plflt) / real(NLEVEL,kind=plflt)
!      enddo
      do i = 1, CONTOUR_MAXCURVES+1
        shedge(i) = cmin + (cmax - cmin) * real(i-1,kind=plflt) / &
                            real(CONTOUR_MAXCURVES,kind=plflt)
      enddo

      ! Here we hard-wire the color levels to the same as the kml plot
      !shedge(:) = (/ 0.0_plflt,  0.1_plflt,   0.3_plflt,   1.0_plflt, 2.0_plflt, &
      !              10.0_plflt, 30.0_plflt, 100.0_plflt, 300.0_plflt, & 
      !              1000.0_plflt, 3000.0_plflt, 10000.0_plflt /)
      fill_width = 2
      cont_color = 0
      cont_width = 0
      axis_opts(1) = 'bcvtm'
      axis_ticks(1) = 0.0_plflt
      axis_subticks(1) = 0
      label_opts(1) = PL_COLORBAR_LABEL_RIGHT
      labels(1) = 'Ash conc. mg/m3'

      write(coord_str,101)x_vprofile(vprof_ID),y_vprofile(vprof_ID)
 101  format(' (lon=',f7.2,', lat=',f6.2,')')
      write(title_str,*)trim(adjustl(Site_vprofile(vprof_ID))),coord_str

      ! Set up for plplot
      call plsdev("pngcairo")      ! Set output device (png, pdf, etc.)
      call plsfnam ( dp_pngfile )  ! Set output filename

      ! set image size and background color via command-line options tool plsetopt
      if(lib_ver_minor.lt.12)then
        ! prior to v5.12, this was a subroutine call
!        call plsetopt("geometry","854x603, 854x603")
!        call plsetopt("bg","FFFFFF")                  ! Set background color to white
      else
        ! after v5.12, this is a function call returing an error code plsetopt_rc
        plsetopt_rc = plsetopt("geometry","854x603, 854x603")
        plsetopt_rc = plsetopt("bg","FFFFFF")
      endif

      ! sets cmap0 palette via pal file for discrete elements
      call plspal0('cmap0_black_on_white.pal')
      ! sets cmap1 palette via pal file for continuous elements
      if(lib_ver_minor.lt.12)then
        ! prior to v5.12, the second argument was an integer
!        call plspal1('cmap1_blue_yellow.pal',1)
      else
        ! after v5.12, this needs to be a logical
        call plspal1('cmap1_blue_yellow.pal',.true.)
      endif

      call plscmap0n(3)  ! Set number of colors in cmap0

      ! Initialize plplot
      call plinit()

        ! pladv: Advance the (sub-)page
      call pladv(0)
        ! plvpor: Specify viewport using normalized subpage coordinates
      call plvpor(0.12_plflt, 0.75_plflt, 0.3_plflt, 0.8_plflt)
        ! plwind: Specify window 
      call plwind(tmin,tmax,zmin,zmax)
        ! plpsty: Select area fill pattern (0 is for solid)
      call plpsty(0)

      ! Now plot the data
      if(lib_ver_minor.lt.12)then
        ! prior to v5.12, plshades assumed coordinate trans. would deform
        ! rectangles
!        call plshades(conc(:ntmax,:nzmax), defined, &
!          tmin,tmax,zmin,zmax, &
!          shedge, fill_width, &
!          cont_color, cont_width )
      else
        ! after v5.12, a rectangular boolean is required
        call plshades(conc(:ntmax,:nzmax), &
          tmin,tmax,zmin,zmax, &
          shedge, fill_width, &
          cont_color, cont_width , .true.)
      endif

      ! Smaller text, scale by second argument
      call plschr( 0.0_plflt, 0.5_plflt )
      ! Small ticks on the vertical axis
      call plsmaj( 0.0_plflt, 0.5_plflt )
      call plsmin( 0.0_plflt, 0.5_plflt )
      call plcol0(1)
      call plbox('bcnst', 0.0_plflt, 0, 'bcnstv', 0.0_plflt, 0)
      call plcol0(2)
      call plschr(0.0_plflt,1.0_plflt)  ! Change font scale
      call pllab("Time (hours after eruption)", "Height (km)", trim(adjustl(title_str)))

!      num_values(1) = NLEVEL + 1;
      num_values(1) = CONTOUR_MAXCURVES + 1;

      values(1,:)   = shedge;
      call plcolorbar( colorbar_width, colorbar_height, &  ! these are output values
            PL_COLORBAR_SHADE ,  &
            !ior(PL_COLORBAR_SHADE, PL_COLORBAR_SHADE_LABEL), & ! sets options
            0, &                                               ! sets position
            0.015_plflt, 0.1_plflt, 0.0375_plflt, 0.8_plflt, &  ! x,y, dimesions
            0, 1, 1, &                                         ! bg_color,bb_color,bb_style
            0.0_plflt, shedge(CONTOUR_MAXCURVES+1), & ! low/high 
!            0.0_plflt, shedge(NLEVEL+1), & ! low/high 
            cont_color, cont_width, &
            label_opts, labels, &
            axis_opts, &
            axis_ticks, axis_subticks, &
            num_values, values )

      ! This plots a contour line at the threshold level
      !call plcont(conc,1,ntmax,1,nzmax,clevel, tr)

      ! Now add the annotation box
        ! pladv: Advance the (sub-)page
      call pladv(1)
        ! plvpor: Specify viewport using normalized subpage coordinates
      call plvpor(0.05_plflt, 0.35_plflt, 0.05_plflt, 0.2_plflt)
        ! plwind: Specify window 
      call plwind(0.0_plflt, 1.0_plflt, 0.0_plflt, 1.0_plflt )
        ! plbox: Draw a box with axes, etc
        !       (xopt, xtick, nxsub, yopt, ytick, nysub)
      !call plbox('bc', 0.0_plflt, 0, 'bc', 0.0_plflt, 0 )
        !  plschr: Set character size
        !          (def, scale)
      call plschr( 0.0_plflt, 0.7_plflt )
      dy_newline = 0.13_plflt
        ! plptex: Write text inside the viewport
      write(outstring,*)"Volcano: ",trim(adjustl(VolcanoName))
      call plptex(0.02_plflt, 1.0_plflt-1.0_plflt*dy_newline, 1.0_plflt, 0.0_plflt, &
                  0.0_plflt, outstring )
      write(outstring,*)"Run Date: ",trim(adjustl(os_time_log))
      call plptex(0.02_plflt, 1.0_plflt-2.0_plflt*dy_newline, 1.0_plflt, 0.0_plflt, &
                  0.0_plflt, outstring )
      read(cdf_b3l1,*,iostat=ioerr) iw,iwf
      write(outstring,*)"Windfile: ",iwf
      call plptex(0.02_plflt, 1.0_plflt-3.0_plflt*dy_newline, 1.0_plflt, 0.0_plflt, &
                  0.0_plflt, outstring )

        ! pladv: Advance the (sub-)page
      call pladv(2)
        ! plvpor: Specify viewport using normalized subpage coordinates
      call plvpor(0.4_plflt, 0.75_plflt, 0.05_plflt, 0.2_plflt)
        ! plwind: Specify window 
      call plwind(0.0_plflt, 1.0_plflt, 0.0_plflt, 1.0_plflt )
        ! plbox: Draw a box with axes, etc
        !       (xopt, xtick, nxsub, yopt, ytick, nysub)
      !call plbox('bc', 0.0_plflt, 0, 'bc', 0.0_plflt, 0 )

      write(outstring,*)"Erup. Start Time: ",HS_xmltime(SimStartHour+e_StartTime(1),BaseYear,useLeap)
      call plptex(0.02_plflt, 1.0_plflt-1.0_plflt*dy_newline, 1.0_plflt, 0.0_plflt, &
                  0.0_plflt, outstring )
      write(outstring,111)e_PlumeHeight(1)
 111  format(' Erup. Plume Height: ',f7.2,' km')
      call plptex(0.02_plflt, 1.0_plflt-2.0_plflt*dy_newline, 1.0_plflt, 0.0_plflt, &
                  0.0_plflt, outstring )
      write(outstring,121)e_Duration(1)
 121  format(' Erup. Duration: ',f7.2,' hours')
      call plptex(0.02_plflt, 1.0_plflt-3.0_plflt*dy_newline, 1.0_plflt, 0.0_plflt, &
                  0.0_plflt, outstring )
      write(outstring,131)e_Volume(1)
 131  format(' Erup. Volume: ',f10.5,' km3 (DRE)')
      call plptex(0.02_plflt, 1.0_plflt-4.0_plflt*dy_newline, 1.0_plflt, 0.0_plflt, &
                  0.0_plflt, outstring )


      ! Close PLplot library
      call plend

      ! clean up memory
      if(allocated(t))      deallocate(t)
      if(allocated(z))      deallocate(z)
      if(allocated(conc))   deallocate(conc)
      if(allocated(shedge)) deallocate(shedge)

      end subroutine write_2Dprof_PNG_plplot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  write_DepPOI_TS_PNG_plplot
!
!  Called from: Ash3d_PostProc.f90
!  Arguments:
!    pt_indx       = index of point in Ash3d netcdf output file
!
!  This subroutine creates a png plot of the transient deposit accumulation at
!  the airport/POI given by pt_index using the plplot graphics package.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine write_DepPOI_TS_PNG_plplot(pt_indx)

      use Airports,      only : &
         Airport_Name,Airport_Thickness_TS

      use io_data,       only : &
         nWriteTimes,WriteTimes

      use time_data,     only : &
         Simtime_in_hours

      integer,intent(in) :: pt_indx

      real(kind=dp)     :: ymaxpl
      character(len=14) :: dp_pngfile
      integer,save      :: plot_index = 0
      integer           :: plsetopt_rc

      real(kind=plflt) :: xmin
      real(kind=plflt) :: xmax
      real(kind=plflt) :: ymin
      real(kind=plflt) :: ymax
      real(kind=plflt), dimension(:), allocatable :: x, y, x0, y0
      integer(kind=4) :: r1 = 0
      integer(kind=4) :: g1 = 0
      integer(kind=4) :: b1 = 0
      integer(kind=4) :: r2 = 136
      integer(kind=4) :: g2 = 136
      integer(kind=4) :: b2 = 136

      if(Airport_Thickness_TS(pt_indx,nWriteTimes).lt.0.01_ip)then
        return
      else
        plot_index = plot_index + 1
      endif

      write(dp_pngfile,55) plot_index,".png"
 55   format('plplt_',i4.4,a4)

      if(Airport_Thickness_TS(plot_index,nWriteTimes).lt.0.01)then
        ymaxpl = 1.0
      elseif(Airport_Thickness_TS(plot_index,nWriteTimes).lt.1.0)then
        ymaxpl = 1.0
      elseif(Airport_Thickness_TS(plot_index,nWriteTimes).lt.5.0)then
        ymaxpl = 5.0
      elseif(Airport_Thickness_TS(plot_index,nWriteTimes).lt.25.0)then
        ymaxpl = 25.0
      else
        ymaxpl = 100.0
      endif

      xmin=real(0,kind=plflt)
      xmax=real(ceiling(Simtime_in_hours),kind=plflt)
      ymin=real(0,kind=plflt)
      ymax=real(ymaxpl,kind=plflt)

      allocate(x(nWriteTimes))
      allocate(y(nWriteTimes))
      allocate(x0(nWriteTimes+1))
      allocate(y0(nWriteTimes+1))

      x = real(WriteTimes,kind=plflt)
      y = real(Airport_Thickness_TS(pt_indx,1:nWriteTimes),kind=plflt)
      x0(1:nWriteTimes)=x(1:nWriteTimes)
      x0(nWriteTimes+1)=x(nWriteTimes)
      y0(1:nWriteTimes)=y(1:nWriteTimes)
      y0(nWriteTimes+1)=0.0_plflt

      ! Set up for plplot
      call plsdev("pngcairo")      ! Set output device (png, pdf, etc.)
      call plsfnam ( dp_pngfile )  ! Set output filename

      ! set image size and background colog via command-line options tool plsetopt
      if(lib_ver_minor.lt.12)then
        ! prior to v5.12, this was a subroutine call
!        call plsetopt("geometry","400x300, 400x300")
!        call plsetopt("bg","FFFFFF")
      elseif(lib_ver_major.le.5.and.lib_ver_minor.eq.15)then
        ! after v5.12, this is a function call returing an error code plsetopt_rc
        plsetopt_rc = plsetopt("geometry","400x300, 400x300")  ! Set image size
        plsetopt_rc = plsetopt("bg","FFFFFF")                  ! Set background color to white
      endif


      ! Initialize plplot
      call plinit()

      ! Create a labelled box to hold the plot.
      call plscolbg(r1,g1,b1) ! Set color back to black
      call plcol0(0)
      call plschr(0.0_plflt,1.7_plflt)  ! Chage font scale
      call plenv( xmin, xmax, ymin, ymax, 0, 0 )
      call pllab( "Time (hours after eruption)", "Deposit Thickeness (mm)", Airport_Name(pt_indx))

      call plscolbg(r2,g2,b2) ! Set pen color to grey
      call plcol0(0)
      call plline( x , y ) ! Simple line plot
      call plfill( x0(1:nWriteTimes+1), y0(1:nWriteTimes+1) )

      ! Close PLplot library
      call plend

      ! clean up memory
      if(allocated(x))  deallocate(x)
      if(allocated(y))  deallocate(y)
      if(allocated(x0)) deallocate(x0)
      if(allocated(y0)) deallocate(y0)

      end subroutine write_DepPOI_TS_PNG_plplot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module Ash3d_PostProc_plplot

!##############################################################################
