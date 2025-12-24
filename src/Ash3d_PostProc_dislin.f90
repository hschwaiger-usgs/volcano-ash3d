!##############################################################################
!
! Ash3d_PostProc_dislin module
!
! This module provides the subroutines that use DISLIN for creating 2d maps,
! 2d vertical profiles, and the little deposit accumulation plots linked to
! the airport arrival kml (ash_arrivaltimes_airports.kml).
! The dislin library is linked at compile-time. DISLIN is available from
!   https://www.dislin.de
!
!      subroutine write_2Dmap_PNG_dislin
!      subroutine write_2Dprof_PNG_dislin
!      subroutine write_DepPOI_TS_PNG_dislin
!
!##############################################################################

      module Ash3d_PostProc_dislin

      use precis_param

      use io_units

!      use global_param,  only : &
!         DirDelim

      use io_data,       only : &
         Instit_IconFile

      use dislin

      implicit none

        ! Set everything to private by default
      private

        ! Publicly available subroutines/functions
      public write_2Dmap_PNG_dislin,    &
             write_2Dprof_PNG_dislin,   &
             write_DepPOI_TS_PNG_dislin

        ! Publicly available variables

      integer,parameter :: DS = 8
      
      contains
      !------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  write_2Dmap_PNG_dislin
!
!  Called from: Ash3d_PostProc.f90
!  Arguments:
!    nx            = x length of output array OutVar
!    ny            = y length of output array OutVar
!    iprod         = product ID
!    itime         = time index from netcdf data file
!    OutVar        = 2-d array to be written to ASCII file
!    Fill_Value    = NaN value
!    writeContours = logical
!
!  This subroutine creates a png map of the variable in OutVar using the dislin
!  graphics package.  Annotations and contour levels are indicated via the
!  product ID (iprod).  If writeContours is set to true, then this subroutine
!  is only used for generating and storing the contours (for dislin, this is
!  directly in memory as a variable) with no png written.
!  If timestep = -1, then use the last step in file.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine write_2Dmap_PNG_dislin(nx,ny,iprod,itime,OutVar,Fill_Value,writeContours)

      use mesh,          only : &
         IsLatLon,lon_cc_pd,lat_cc_pd,de,dn, &
         dx,dy,                              &
         latLL,lonLL,latUR,lonUR,            &
         xLL,yLL,xUR,yUR

      use Output_Vars,   only : &
         ContourFilled, &
         Con_Cust,Con_Cust_N,Con_Cust_RGB,Con_Cust_Lev,&
         Con_DepThick_mm_N,Con_DepThick_mm_Lev,Con_DepThick_mm_RGB, &
         Con_DepThick_in_N,Con_DepThick_in_Lev,Con_DepThick_in_RGB, &
         Con_DepTime_N,Con_DepTime_Lev,Con_DepTime_RGB, &
         Con_CloudCon_N,Con_CloudCon_Lev,Con_CloudCon_RGB, &
         Con_CloudTop_N,Con_CloudTop_RGB,Con_CloudTop_Lev, &
         Con_CloudBot_N,Con_CloudBot_RGB,Con_CloudBot_Lev, &
         Con_CloudLoad_N,Con_CloudLoad_RGB,Con_CloudLoad_Lev, &
         Con_CloudRef_N,Con_CloudRef_RGB,Con_CloudRef_Lev, &
         Con_CloudTime_N,Con_CloudTime_RGB,Con_CloudTime_Lev, &
         ContourLev,nConLev, &
         ContourDataX,ContourDataY,ContourDataNcurves,ContourDataNpoints,&
         CONTOUR_MAXCURVES,CONTOUR_MAXPOINTS

      use io_data,       only : &
         WriteTimes,cdf_b3l1,VolcanoName

      use Source,        only : &
         neruptions,e_Volume,e_Duration,e_StartTime,e_PlumeHeight,lon_volcano,lat_volcano


      use time_data,     only : &
         os_time_log,SimStartHour,BaseYear,useLeap

      use citywriter

      integer      ,intent(in) :: nx
      integer      ,intent(in) :: ny
      integer      ,intent(in) :: iprod
      integer      ,intent(in) :: itime
      real(kind=ip),intent(in) :: OutVar(nx,ny)
      real(kind=ip),intent(in) :: Fill_Value
      logical      ,intent(in) :: writeContours

      integer :: tmp_int
      integer,dimension(:,:),allocatable :: zrgb
      character(len=40) :: title_plot
      character(len=30) :: cstr_xlabel = 'Longitude'
      character(len=30) :: cstr_ylabel = 'Latitude'
      character(len=30) :: cstr_zlabel
      character(len=30) :: cstr_volcname
      character(len=30) :: cstr_run_date
      character(len=30) :: cstr_windfile
      character(len=40) :: cstr_ErStartT
      character(len=27) :: cstr_ErHeight
      character(len=30) :: cstr_ErDuratn
      character(len=38) :: cstr_ErVolume
      character(len=45) :: cstr_note
      character(len=20) :: varname
      character(len= 9) :: cio
      character(len= 4) :: outfile_ext = '.png'
      character(len=10) :: units
      integer           :: iostatus
      character(len=120):: iomessage
      integer           :: iw,iwf
      logical           :: HaveIconFile
      character(len=80) :: linebuffer080

      ! Plot dimensions
      real(kind=ip)  :: xmin
      real(kind=ip)  :: xmax
      real(kind=ip)  :: ymin
      real(kind=ip)  :: ymax
      logical        :: IsRegGrid

      ! Aux. File names
      character(len= 8) :: filename_root
      character(len=40) :: filename_png
      ! Citywriter variables
      integer :: icty
      integer :: ncities
      integer :: cityname_offset_px = 30
      real(kind=ip),dimension(:),allocatable     :: lon_cities
      real(kind=ip),dimension(:),allocatable     :: lat_cities
      character(len=26),dimension(:),allocatable :: name_cities

      ! Contour variables
      integer           :: ilev        ! number of contour levels
      !integer           :: substr_pos1
      integer           :: icurve      ! number of curves for level ilev
      integer           :: npts        ! number of points in curve ilev,icurve

      ! Plotting variables
      integer :: plotx_px       = 2200
      integer :: ploty_px       = 1400
      integer :: plotposx_px    = 500
      integer :: plotposy_px    = 1650
      integer :: y_footer_px    = 1900
      integer :: x_leg1_px      = 400
      integer :: x_leg2_px      = 1200
      integer :: x_leg3_px      = 2250
      integer :: dy_newline_px  = 40

      ! DISLIN variables
      real(kind=DS) :: xminDS
      real(kind=DS) :: xmaxDS
      real(kind=DS) :: yminDS
      real(kind=DS) :: ymaxDS
      real(kind=DS) :: dx_map, dy_map    ! lon and lat grid line increment
      real(kind=DS) :: xgrid_1, ygrid_1  ! gridlines

      character(len=4) :: cfmt = "PNG " ! output driver/file-format (PNG); this is a 4-char string
      character(len=4) :: cfsz = "USAL" ! US A Landscape (2790 x 2160)
      character(len=3) :: cmode = "ON " ! Plotting mode; can be off if just calc. contours
      integer          :: nclr          ! color index
      integer          :: nmaxln        ! number of characters in the longest line of text
      character(len=7) :: zlevlab       ! legend level label

      integer       :: NCURVS                   ! number of curves for each level: ContourDataNcurves(ilev)
      integer       :: iray(CONTOUR_MAXCURVES)  ! number of point on curve ContourDataNpoints(ilev,icurve)
      real(kind=DS) :: xpts(CONTOUR_MAXPOINTS)  ! coordinates of the points on the curve
      real(kind=DS) :: ypts(CONTOUR_MAXPOINTS)

      real(kind=DS) :: xp,yp          ! point coordinate (in plot system) as real
      integer       :: nxp,nyp        ! point coordinate (in plot system) as integer
      real(kind=DS) :: xr,xg,xb       ! RGB colors

      INTERFACE
        character (len=20) function HS_xmltime(HoursSince,byear,useLeaps)
          real(kind=8),intent(in) :: HoursSince
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_xmltime
      END INTERFACE

      ! Test for icon file
      inquire( file=trim(adjustl(Instit_IconFile)), exist=HaveIconFile)

      ncities = 20
      allocate(lon_cities(ncities))
      allocate(lat_cities(ncities))
      allocate(name_cities(ncities))
      filename_root        = "outvar"

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
        varname = "depothick"
        write(filename_png,'(a15,a9,a4)')'Ash3d_Deposit_t',cio,outfile_ext
        write(title_plot,'(a20,f5.2,a6)')'Deposit Thickness t=',WriteTimes(itime),' hours'
        cstr_zlabel = 'Dep.Thick.(mm)'
        units = " (mm)"
        if(.not.Con_Cust)then
          nConLev = Con_DepThick_mm_N
          allocate(zrgb(nConLev,3))
          allocate(ContourLev(nConLev))
          ContourLev(1:nConLev) = Con_DepThick_mm_Lev(1:nConLev)
          zrgb(1:nConLev,1:3) = Con_DepThick_mm_RGB(1:nConLev,1:3)
        endif
      elseif(iprod.eq.4)then   ! deposit at specified times (inches)
        varname = "depothick"
        write(filename_png,'(a15,a9,a4)')'Ash3d_Deposit_t',cio,outfile_ext
        write(title_plot,'(a20,f5.2,a6)')'Deposit Thickness t=',WriteTimes(itime),' hours'
        cstr_zlabel = 'Dep.Thick.(in)'
        units = " (in)"
        if(.not.Con_Cust)then
          nConLev = Con_DepThick_in_N
          allocate(zrgb(nConLev,3))
          allocate(ContourLev(nConLev))
          ContourLev(1:nConLev) = Con_DepThick_in_Lev(1:nConLev)
          zrgb(1:nConLev,1:3) = Con_DepThick_in_RGB(1:nConLev,1:3)
        endif
      elseif(iprod.eq.5)then       ! deposit at final time (mm)
        varname = "depothickFin"
        write(filename_png,'(a13,a9,a4)')'Ash3d_Deposit',cio,outfile_ext
        title_plot = 'Final Deposit Thickness'
        cstr_zlabel = 'Dep.Thick.(mm)'
        units = " (mm)"
        if(.not.Con_Cust)then
          nConLev = Con_DepThick_mm_N
          allocate(zrgb(nConLev,3))
          allocate(ContourLev(nConLev))
          ContourLev(1:nConLev) = Con_DepThick_mm_Lev(1:nConLev)
          zrgb(1:nConLev,1:3) = Con_DepThick_mm_RGB(1:nConLev,1:3)
        endif
      elseif(iprod.eq.6)then   ! deposit at final time (inches)
        varname = "depothickFin"
        write(filename_png,'(a13,a9,a4)')'Ash3d_Deposit',cio,outfile_ext
        title_plot = 'Final Deposit Thickness'
        cstr_zlabel = 'Dep.Thick.(in)'
        units = " (in)"
        if(.not.Con_Cust)then
          nConLev = Con_DepThick_in_N
          allocate(zrgb(nConLev,3))
          allocate(ContourLev(nConLev))
          ContourLev(1:nConLev) = Con_DepThick_in_Lev(1:nConLev)
          zrgb(1:nConLev,1:3) = Con_DepThick_in_RGB(1:nConLev,1:3)
        endif
      elseif(iprod.eq.7)then   ! ashfall arrival time (hours)
        varname = "depotime"
        write(filename_png,'(a22)')'DepositArrivalTime.png'
        write(title_plot,'(a20)')'Ashfall arrival time'
        cstr_zlabel = 'Time (hours)'
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
        varname = "ashcon_max"
        write(filename_png,'(a16,a9,a4)')'Ash3d_CloudCon_t',cio,outfile_ext
        write(title_plot,'(a26,f5.2,a6)')'Ash-cloud concentration t=',WriteTimes(itime),' hours'
        cstr_zlabel = 'Max.Con.(mg/m3)'
        units = " (mg/m3)"
        if(.not.Con_Cust)then
          nConLev = Con_CloudCon_N
          allocate(zrgb(nConLev,3))
          allocate(ContourLev(nConLev))
          ContourLev(1:nConLev) = Con_CloudCon_Lev(1:nConLev)
          zrgb(1:nConLev,1:3) = Con_CloudCon_RGB(1:nConLev,1:3)
        endif
      elseif(iprod.eq.10)then   ! ash-cloud height
        varname = "cloud_height"
        write(filename_png,'(a19,a9,a4)')'Ash3d_CloudHeight_t',cio,outfile_ext
        write(title_plot,'(a19,f5.2,a6)')'Ash-cloud height t=',WriteTimes(itime),' hours'
        cstr_zlabel = 'Cld.Height(km)'
        units = " (km)"
        if(.not.Con_Cust)then
          nConLev = Con_CloudTop_N
          allocate(zrgb(nConLev,3))
          allocate(ContourLev(nConLev))
          ContourLev(1:nConLev) = Con_CloudTop_Lev(1:nConLev)
          zrgb(1:nConLev,1:3) = Con_CloudTop_RGB(1:nConLev,1:3)
        endif
      elseif(iprod.eq.11)then   ! ash-cloud bottom
        varname = "cloud_bottom"
        write(filename_png,'(a16,a9,a4)')'Ash3d_CloudBot_t',cio,outfile_ext
        write(title_plot,'(a19,f5.2,a6)')'Ash-cloud bottom t=',WriteTimes(itime),' hours'
        cstr_zlabel = 'Cld.Bot.(km)'
        units = " (km)"
        if(.not.Con_Cust)then
          nConLev = Con_CloudBot_N
          allocate(zrgb(nConLev,3))
          allocate(ContourLev(nConLev))
          ContourLev(1:nConLev) = Con_CloudBot_Lev(1:nConLev)
          zrgb(1:nConLev,1:3) = Con_CloudBot_RGB(1:nConLev,1:3)
        endif
      elseif(iprod.eq.12)then   ! ash-cloud load
        varname = "cloud_load"
        write(filename_png,'(a17,a9,a4)')'Ash3d_CloudLoad_t',cio,outfile_ext
        write(title_plot,'(a17,f5.2,a6)')'Ash-cloud load t=',WriteTimes(itime),' hours'
        cstr_zlabel = 'Cld.Load(T/km2)'
        units = " (T/km2)"
        if(.not.Con_Cust)then
          nConLev = Con_CloudLoad_N
          allocate(zrgb(nConLev,3))
          allocate(ContourLev(nConLev))
          ContourLev(1:nConLev) = Con_CloudLoad_Lev(1:nConLev)
          zrgb(1:nConLev,1:3) = Con_CloudLoad_RGB(1:nConLev,1:3)
        endif
      elseif(iprod.eq.13)then  ! radar reflectivity
        varname = "radar_reflectivity"
        write(filename_png,'(a20,a9,a4)')'Ash3d_CloudRadRefl_t',cio,outfile_ext
        write(title_plot,'(a24,f5.2,a6)')'Ash-cloud radar refl. t=',WriteTimes(itime),' hours'
        cstr_zlabel = 'Cld.Refl.(dBz)'
        units = " (dBz)"
        if(.not.Con_Cust)then
          nConLev = Con_CloudRef_N
          allocate(zrgb(nConLev,3))
          allocate(ContourLev(nConLev))
          ContourLev(1:nConLev) = Con_CloudRef_Lev(1:nConLev)
          zrgb(1:nConLev,1:3) = Con_CloudRef_RGB(1:nConLev,1:3)
        endif
      elseif(iprod.eq.14)then   ! ashcloud arrival time (hours)
        varname = "ash_arrival_time"
        write(filename_png,'(a20)')'CloudArrivalTime.png'
        write(title_plot,'(a22)')'Ash-cloud arrival time'
        cstr_zlabel = 'Time (hours)'
        units = " (hours)"
        if(.not.Con_Cust)then
          nConLev = Con_CloudTime_N
          allocate(zrgb(nConLev,3))
          allocate(ContourLev(nConLev))
          ContourLev(1:nConLev) = Con_CloudTime_Lev(1:nConLev)
          zrgb(1:nConLev,1:3) = Con_CloudTime_RGB(1:nConLev,1:3)
        endif
      elseif(iprod.eq.15)then   ! topography
        varname = "topography"
        write(filename_png,'(a14)')'Topography.png'
        write(title_plot,'(a10)')'Topography'
        cstr_zlabel = 'Elevation (km)'
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
          write(errlog(io),*)"         iprod = ",iprod
          write(errlog(io),*)"       Cannot map this variable."
        endif;enddo
        stop 1
      endif
      ! Now have string vars (varname,cstr_zlabel, etc.) and contour info (nConLev,zrgb,ContourLev)

      ! Dislin allows us to directly generate contours from OutVar, whereas the
      ! scripted graphics packages require building a script, 'plotting', then
      ! reading in the contour lines. Here, we can just contour then exit the subroutine.
      !call nancrv(cmode)
      if(writeContours)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Running dislin to calculate contours lines"
        endif;enddo
        allocate(ContourDataNcurves(nConLev))
        allocate(ContourDataNpoints(nConLev,CONTOUR_MAXCURVES))
        allocate(ContourDataX(nConLev,CONTOUR_MAXCURVES,CONTOUR_MAXPOINTS))
        allocate(ContourDataY(nConLev,CONTOUR_MAXCURVES,CONTOUR_MAXPOINTS))
        ContourDataNcurves(:)   = 0
        ContourDataNpoints(:,:) = 0
        ContourDataX(:,:,:)     = 0.0_ip
        ContourDataY(:,:,:)     = 0.0_ip
        do ilev=1,nConLev
          ! This part calculates the contours
          xpts(1:CONTOUR_MAXPOINTS) = 0.0_DS
          ypts(1:CONTOUR_MAXPOINTS) = 0.0_DS
          iray(1:CONTOUR_MAXCURVES) = 0
          call conpts(real(lon_cc_pd(1:nx),kind=DS),nx,&  ! x coord and size
                      real(lat_cc_pd(1:ny),kind=DS),ny,&  ! y coord and size
                      real(OutVar(1:nx,1:ny),kind=DS), &  ! matrix with function values
                      real(ContourLev(ilev),kind=DS),  &  ! level to contour
                      xpts(1:CONTOUR_MAXPOINTS),   &  ! x of contour (may have mul. curvex)
                      ypts(1:CONTOUR_MAXPOINTS),   &  ! y of contour (may have mul. curves)
                      CONTOUR_MAXPOINTS,           &  ! max # of points for contour arrays
                      iray(1:CONTOUR_MAXCURVES),   &  ! num of points for each contour
                      CONTOUR_MAXCURVES,           &  ! max number of curves
                      NCURVS)                         ! actual number of curves

          if(NCURVS.gt.CONTOUR_MAXCURVES)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: Maximum number of curves for this level exceeded by conpts"
              write(errlog(io),*)"       Current maximum set to CONTOUR_MAXCURVES = ",CONTOUR_MAXCURVES
              write(errlog(io),*)"       Please increase CONTOUR_MAXCURVES and recompile."
              write(errlog(io),*)"  Output_Vars.f90:CONTOUR_MAXCURVES"
            endif;enddo
            stop 1
          endif
          ContourDataNcurves(ilev)=NCURVS
          do icurve=1,NCURVS
            ContourDataNpoints(ilev,icurve)=iray(icurve)
            if(ContourDataNpoints(ilev,icurve).gt.CONTOUR_MAXPOINTS)then
              do io=1,2;if(VB(io).le.verbosity_error)then
                write(errlog(io),*)"ERROR: Maximum number of points for this curve exceeded by conpts"
                write(errlog(io),*)"       Current maximum set to CONTOUR_MAXPOINTS = ",CONTOUR_MAXPOINTS
                write(errlog(io),*)"       Please increase CONTOUR_MAXPOINTS and recompile."
                write(errlog(io),*)"  Output_Vars.f90:CONTOUR_MAXPOINTS"
              endif;enddo
              stop 1 
            endif
            npts = iray(icurve)
            !do k=1,iray(icurve)
            ContourDataX(ilev,icurve,1:npts) = real(xpts(1:npts),kind=ip)
            ContourDataY(ilev,icurve,1:npts) = real(ypts(1:npts),kind=ip)
            !enddo
            ! These data could be plotted to produce the same plot as from contur
            !call curvmp(xpts,ypts,iray(icurve))
          enddo
        enddo
        ! Once we've loaded contours, we are all done here
        return
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Running dislin to generate contour plot"
        endif;enddo
      endif

      ! This is the section where we actually start plotting the map
      ! Evaluate grid
      if(IsLatLon)then
        xmin = lonLL
        xmax = lonUR
        ymin = latLL
        ymax = latUR
        if(abs(dn-de).lt.1.0e-4_ip)then
          IsRegGrid = .true.
        else
          IsRegGrid = .false.
        endif
      else
        xmin = xLL
        xmax = xUR
        ymin = yLL
        ymax = yUR
        if(abs(dx-dy).lt.1.0e-4_ip)then
          IsRegGrid = .true.
        else
          IsRegGrid = .false.
        endif
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Currenntly, plotting with dislin only enabled for lon/lat grids."
          write(errlog(io),*)"       Please use GMT to plot projected maps."
          write(errlog(io),*)"       ./ASH3DPLOT=4 ./Ash3d_PostProc ...."
        endif;enddo
        stop 1
      endif

      ! No prep needed for the data in a form that dislin can read since it is internal

      call citylist(0,                        &  ! 0 is for internal list only (no file)
                    xmin,xmax,ymin,ymax,      &
                    ncities,                  &
                    lon_cities,               &
                    lat_cities,               &
                    name_cities)

      ! Build strings with run info for legend
      ! Volcano:     Erup.start:
      ! Run date:    Plm Height:
      ! Windfile:    Duration:
      !              Volume:
      write(cstr_volcname,'(a10,a20)')'Volcano:  ' ,VolcanoName
      write(cstr_run_date,'(a10,a20)')'Run Date: ',os_time_log
      read(cdf_b3l1,*,iostat=iostatus,iomsg=iomessage) iw,iwf
      write(cstr_windfile,'(a10,i5)')'Windfile: ',iwf
      if(neruptions.gt.1)then
        write(cstr_note,'(a45)')'WARNING: Multiple eruptions, only first given'
      endif

      !e_StartTime,e_PlumeHeight,e_Duration,e_Volume
      write(cstr_ErStartT,'(a20,a20)')'Erup. Start Time:   ',&
            HS_xmltime(SimStartHour+e_StartTime(1),BaseYear,useLeap)
      write(cstr_ErHeight,'(a20,f4.1,a3)')'Erup. Plume Height: ',e_PlumeHeight(1),' km'
      write(cstr_ErDuratn,'(a20,f4.1,a6)')'Erup. Duration:     ',e_Duration(1),' hours'
      write(cstr_ErVolume,'(a20,f8.5,a10)')'Erup. Volume:       ',e_Volume(1),' km3 (DRE)'


      ! Set gridlines
      dx_map = 10.0_DS
      dy_map = 5.0_DS
      xminDS = real(xmin- 360.0_ip,kind=DS)
      xmaxDS = real(xmax- 360.0_ip,kind=DS)
      yminDS = real(ymin,kind=DS)
      ymaxDS = real(ymax,kind=DS)
      if(xmaxDS-xminDS.gt.40.0_DS)then
        dx_map = 10.0_DS
      elseif(xmaxDS-xminDS.gt.10.0_DS)then
        dx_map = 5.0_DS
      elseif(xmaxDS-xminDS.gt.5.0_DS)then
        dx_map = 2.0_DS
      else
        dx_map = 1.0_DS
      endif
      if(ymaxDS-yminDS.gt.40.0_DS)then
        dy_map = 10.0_DS
      elseif(ymaxDS-yminDS.gt.10.0_DS)then
        dy_map = 5.0_DS
      elseif(ymaxDS-yminDS.gt.5.0_DS)then
        dy_map = 2.0_DS
      else
        dy_map = 1.0_DS
      endif
      xgrid_1 = real(ceiling(xminDS/dx_map) * dx_map,kind=DS)
      ygrid_1 = real(ceiling(yminDS/dy_map) * dy_map,kind=DS)

      !!!!!!!!!!!!!!!!!!!!!!!
      !  Dislin Level 0:  before initialization or after termination
      call metafl(cfmt) ! set output driver/file-format (PNG); this is a 4-char string
      call setpag(cfsz) ! Set pagesize to US A Landscape (2790 x 2160)
      call setfil(trim(adjustl(filename_png))) ! Set output filename
      call scrmod('REVERSE')  ! Default background is black; reverse to white

      !  Dislin Level 1:  after initialization or a call to ENDGRF
      call disini()       ! initialize plot (set to level 1)
        ! Set the color table : SPEC,RAIN,GREY,TEMP
      call setvlt('RAIN')
      call paghdr('Ash3d Simulation plotted on ','---',4,0)

      if(HaveIconFile)then
        call filbox(x_leg3_px,y_footer_px,130,49)
        call incfil(trim(adjustl(Instit_IconFile)))
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Institution Logo not found."
          write(outlog(io),*)"If you would like your logo on these maps, copy a"
          write(outlog(io),*)"small (130x50) to ASH3DHOME/share/post_proc/logo.png"
        endif;enddo
      endif

       ! setting of plot parameters
      call bmpfnt('HELVE')

      call axspos(plotposx_px,plotposy_px)  ! determine the position of the axis system
      call axslen(plotx_px,ploty_px)        ! defines the size of the axis system
      call name(cstr_xlabel,'X') ! Set x-axis title
      call name(cstr_ylabel,'Y') ! Set y-axis title
      call labdig(1,'X') ! set number of decimal places for x label (-1 means no decimal)
      call labdig(1,'y') 
      call ticks(1,'xy')  ! set number of ticks between labels
      call titlin(title_plot,4)  ! Set the title
      call incmrk(0) ! selects line (0) or symbol (-1) mode for CURVE

      !call LABELS('MAP','xy')
      ! set projection : STER,LAMB,CYLI,MERC
      call projct('CYLI') ! defines projection
      call frame(3)       ! bump up frame line thickness

      !  Dislin Level 2: after a call to GRAF, GRAFP or GRAFMP
        ! Now create graph and set to level 2
       !  The routine GRAFMP plots a geographical axis system.
      call grafmp(xminDS,xmaxDS,xgrid_1,dx_map, &
                  yminDS,ymaxDS,ygrid_1,dy_map)
      call getlev(tmp_int)
      if(tmp_int.ne.2)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"grafmp is supposed to return plot level 2, but we have ",tmp_int
          write(errlog(io),*)"Exiting"
        endif;enddo
        stop 1
      endif
       ! set color of coastlines
       ! plots coastlines and lakes or political borders
      call world()
      ! Add cities
      do icty=1,ncities
      !    These are the points
        call pos2pt(real(lon_cities(icty),kind=DS),real(lat_cities(icty),kind=DS),&
                    xp,yp)
        nxp=nint(xp)
        nyp=nint(yp)
        call hsymbl(35) ! Set size of city marker
        call symbol(21,nxp,nyp) ! Symbol #21 is a filled circle
        !    These are the city labels, offset in x
        call messag(adjustl(trim(name_cities(icty))),nxp+cityname_offset_px,nyp)
      enddo

      ! Add volcano
      call pos2pt(real(lon_volcano,kind=DS),real(lat_volcano,kind=DS),&
                  xp,yp)
      nxp=nint(xp)
      nyp=nint(yp)
      call hsymbl(70) ! Set size of volcano marker
      call color('red')
      call symbol(18,nxp,nyp) ! Symbol #18 is a filled triangle

      lon_cc_pd(:) = lon_cc_pd(:) - 360.0_ip
      if(ContourFilled)then
        call shdmod('UPPER', 'CELL') ! This suppresses colors in regions above/below the zlevels pro
        call conshd(real(lon_cc_pd(1:nx),kind=DS),nx,&
                    real(lat_cc_pd(1:ny),kind=DS),ny,&
                    real(OutVar,kind=DS),real(ContourLev(1:nConLev),kind=DS),nConLev)
      else
        do ilev=1,nConLev
          xr   = real(zrgb(ilev,1),kind=DS)/real(255,kind=DS)
          xg   = real(zrgb(ilev,2),kind=DS)/real(255,kind=DS)
          xb   = real(zrgb(ilev,3),kind=DS)/real(255,kind=DS)
          nclr = intrgb(xr,xg,xb)

          call setclr(nclr)
          call contur(real(lon_cc_pd(1:nx),kind=DS),nx,&
                      real(lat_cc_pd(1:ny),kind=DS),ny,&
                      real(OutVar(1:nx,1:ny),kind=DS), &
                      real(ContourLev(ilev),kind=DS))
        enddo
      endif

       ! set color of grid lines
      call setrgb(0.0_DS, 0.0_DS, 0.0_DS)
       ! overlays an axis system with a longitude and latitude grid
      call gridmp(1,1)
      call height(50) ! Set character height for title
      call title() ! Actually write the title to the file

      ! Now write the legend
      call height(25) ! Reset character height to something smaller
      nmaxln = 6 ! number of characters in the longest line of text
      call legini(linebuffer080,nConLev,nmaxln) ! Initialize legend
      call legtit(cstr_zlabel)       ! Set legend title
      call legbgd(0)                 ! sets background color
      do ilev=1,nConLev
        if(ContourLev(ilev).lt.1.0_ip)then
          write(zlevlab,'(f6.2)')real(ContourLev(ilev),kind=4)
        else
          write(zlevlab,'(f6.1)')real(ContourLev(ilev),kind=4)
        endif
        call leglin(linebuffer080,zlevlab,ilev)
      enddo
      call legend(linebuffer080,6)  ! write buffer to legend and give position (6=LR)

      ! Add boxes below plot with run info
      call messag(cstr_volcname,x_leg1_px,y_footer_px+0*dy_newline_px)
      call messag(cstr_run_date,x_leg1_px,y_footer_px+1*dy_newline_px)
      call messag(cstr_windfile,x_leg1_px,y_footer_px+2*dy_newline_px)
      if(neruptions.gt.1)then
        call messag(cstr_note,x_leg1_px,y_footer_px+3*dy_newline_px)
      endif
      call messag(cstr_ErStartT,x_leg2_px,y_footer_px+0*dy_newline_px)
      call messag(cstr_ErHeight,x_leg2_px,y_footer_px+1*dy_newline_px)
      call messag(cstr_ErDuratn,x_leg2_px,y_footer_px+2*dy_newline_px)
      call messag(cstr_ErVolume,x_leg2_px,y_footer_px+3*dy_newline_px)

      !  Dislin Level 0:  before initialization or after termination
      call disfin()

      ! clean up memory
      if(allocated(lon_cities))         deallocate(lon_cities)
      if(allocated(lat_cities))         deallocate(lat_cities)
      if(allocated(name_cities))        deallocate(name_cities)
      if(allocated(zrgb))               deallocate(zrgb)

      end subroutine write_2Dmap_PNG_dislin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  write_2Dprof_PNG_dislin
!
!  Called from: Ash3d_PostProc.f90
!  Arguments:
!    vprof_ID        = ID of the profile (number in list from Ash3d control file)
!
!  This subroutine creates a png plot of the transient vertical profile of ash
!  concentration above the profile point (vprof_ID) using the dislin graphics
!  package.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine write_2Dprof_PNG_dislin(vprof_ID)

      use global_param,  only : &
         KG_2_MG,KM3_2_M3

      use mesh,          only : &
         IsLatLon,nzmax,z_cc_pd

      use Output_Vars,   only : &
         pr_ash,CLOUDCON_THRESH

      use time_data,     only : &
         ntmax,time_native

      use io_data,       only : &
         Site_vprofile,x_vprofile,y_vprofile,cdf_b3l1,VolcanoName

      use Source,        only : &
         neruptions,e_Volume,e_Duration,e_StartTime,e_PlumeHeight

      use time_data,     only : &
         os_time_log,SimStartHour,BaseYear,useLeap

      integer,intent(in) :: vprof_ID

      logical           :: HaveIconFile
      character(len=76) :: title_plot
      character(len=30) :: cstr_xlabel = 'Time (hours after eruption)'
      character(len=30) :: cstr_ylabel = 'Height (km)'
      character(len=30) :: cstr_zlabel = 'Ash conc. mg/m3'
      character(len=30) :: cstr_volcname
      character(len=30) :: cstr_run_date
      character(len=30) :: cstr_windfile
      character(len=40) :: cstr_ErStartT
      character(len=27) :: cstr_ErHeight
      character(len=30) :: cstr_ErDuratn
      character(len=38) :: cstr_ErVolume
      character(len=45) :: cstr_note

      character(len=10) :: filename_root
      character(len=14) :: filename_png
      character(len=26) :: coord_str
      integer           :: i,k
      integer           :: iostatus
      character(len=120):: iomessage
      integer           :: iw,iwf

      ! Plotting variables
      integer :: plotx_px       = 1800
      integer :: ploty_px       = 1200
      integer :: plotz_px       = 1200
      integer :: plotposx_px    = 500
      integer :: plotposy_px    = 1650
      integer :: y_footer_px    = 1900
      integer :: x_leg1_px      = 400
      integer :: x_leg2_px      = 1200
      integer :: x_leg3_px      = 2250
      integer :: dy_newline_px  = 40

      ! DISLIN variables
      character(len=4) :: cfmt = "PNG " ! output driver/file-format (PNG); this is a 4-char string
      character(len=4) :: cfsz = "USAL" ! US A Landscape (2790 x 2160)
      real(kind=DS) :: tmin    , zmin    , cmin     ! graph minima
      real(kind=DS) :: tmax    , zmax    , cmax     ! graph maxima
      real(kind=DS) :: tlab1   , zlab1   , clab1    ! graph first label
      real(kind=DS) :: tlabstep, zlabstep, clabstep ! graph label increment
      real(kind=DS) :: cloudcon_thresh_mgm3
      real(kind=DS), dimension(:),   allocatable :: t, z
      real(kind=DS), dimension(:,:), allocatable :: conc

      INTERFACE
        character (len=20) function HS_xmltime(HoursSince,byear,useLeaps)
          real(kind=8),intent(in) :: HoursSince
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_xmltime
      END INTERFACE

      ! Test for icon file
      inquire( file=trim(adjustl(Instit_IconFile)), exist=HaveIconFile)

      ! Build strings with run info for legend
      ! Volcano:     Erup.start:
      ! Run date:    Plm Height:
      ! Windfile:    Duration:
      !              Volume:
      write(cstr_volcname,'(a10,a20)')'Volcano:  ' ,VolcanoName
      write(cstr_run_date,'(a10,a20)')'Run Date: ',os_time_log
      read(cdf_b3l1,*,iostat=iostatus,iomsg=iomessage) iw,iwf
      write(cstr_windfile,'(a10,i5)')'Windfile: ',iwf
      if(neruptions.gt.1)then
        write(cstr_note,'(a45)')'WARNING: Multiple eruptions, only first given'
      endif

      !e_StartTime,e_PlumeHeight,e_Duration,e_Volume
      write(cstr_ErStartT,'(a20,a20)')'Erup. Start Time:   ',&
            HS_xmltime(SimStartHour+e_StartTime(1),BaseYear,useLeap)
      write(cstr_ErHeight,'(a20,f4.1,a3)')'Erup. Plume Height: ',e_PlumeHeight(1),' km'
      write(cstr_ErDuratn,'(a20,f4.1,a6)')'Erup. Duration:     ',e_Duration(1),' hours'
      write(cstr_ErVolume,'(a20,f8.5,a10)')'Erup. Volume:       ',e_Volume(1),' km3 (DRE)'

      write(filename_root,52)vprof_ID
 52   format('vprof_',i4.4)
      filename_png     = trim(adjustl(filename_root)) // ".png"

      ! Get min/max and label interval for all three axies.
      tmin=real(0,kind=DS)
      tmax=real(ceiling(time_native(ntmax)),kind=DS)
      tlab1    = 0.0_DS
      if(tmax.gt.240.0_DS)then
        tlabstep = 48.0_DS
      elseif(tmax.gt.120.0_DS)then
        tlabstep = 24.0_DS
      elseif(tmax.gt.30.0_DS)then
        tlabstep = 10.0_DS
      elseif(tmax.gt.15.0_DS)then
        tlabstep = 5.0_DS
      elseif(tmax.gt.6.0_DS)then
        tlabstep = 2.0_DS
      else
        tlabstep = 1.0_DS
      endif

      zmin=real(0,kind=DS)
      zmax=real(z_cc_pd(nzmax),kind=DS)
      zlab1    = 0.0_DS
      if(zmax.gt.30.0_DS)then
        zlabstep = 10.0_DS
      elseif(zmax.gt.15.0_DS)then
        zlabstep = 5.0_DS
      elseif(zmax.gt.6.0_DS)then
        zlabstep = 2.0_DS
      else
        zlabstep = 1.0_DS
      endif

      cloudcon_thresh_mgm3 = real(CLOUDCON_THRESH * KG_2_MG / KM3_2_M3,kind=DS) !convert from kg/km3 to mg/m3
      cmin=real(0,kind=DS)
      cmax=real(maxval(pr_ash(:,:,vprof_ID)),kind=DS)    ! Get the max value for this profile
      cmin=real(min(cmin,cloudcon_thresh_mgm3),kind=DS)  ! Do not let cmax drop below the threshold
      if    (cmax.gt.4.0e4_DS)then
          clabstep = 5.0e3_DS
      elseif(cmax.gt.1.0e4_DS)then
          clabstep = 2.0e3_DS
      elseif(cmax.gt.4.0e3_DS)then
          clabstep = 5.0e2_DS
      elseif(cmax.gt.1.0e3_DS)then
          clabstep = 2.0e2_DS
      elseif(cmax.gt.4.0e2_DS)then
          clabstep = 5.0e1_DS
      elseif(cmax.gt.1.0e2_DS)then
          clabstep = 2.0e1_DS
      elseif(cmax.gt.4.0e1_DS)then
          clabstep = 5.0e0_DS
      elseif(cmax.gt.1.0e1_DS)then
          clabstep = 2.0e0_DS
      elseif(cmax.gt.1.0e0_DS)then
          clabstep = 5.0e-1_DS
      else
          clabstep = 1.0e-1_DS
      endif
      clab1    = 0.0_DS

      ! Prep data: dislin stores a 2d array internally
      allocate(t(ntmax))
      allocate(z(nzmax))
      allocate(conc(ntmax,nzmax))
      t = real(time_native(1:ntmax),kind=DS)
      z = real(z_cc_pd(1:nzmax),kind=DS)
      do i=1,ntmax
        do k=1,nzmax
          conc(i,k) = real(pr_ash(k,i,vprof_ID),kind=DS)
        enddo
      enddo

      ! Build the plot title
      if(IsLatLon)then
        write(coord_str,101)x_vprofile(vprof_ID),y_vprofile(vprof_ID)
      else
        write(coord_str,102)x_vprofile(vprof_ID),y_vprofile(vprof_ID)
      endif
 101  format(' (lon=',f7.2,', lat=',f6.2,')')
 102  format(' (x=',f9.3,', y=',f9.3,')')
      write(title_plot,*)trim(adjustl(Site_vprofile(vprof_ID))),coord_str

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! DISLIN block
      ! https://www.dislin.de/
      ! wget https://www.dislin.de/downloads/linux/i586_64/dislin-11.5.linux.i586_64.tar.gz
      !(1)    setting of page format, file format and filename
      !(2)    initialization
      !(3)    setting of plot parameters
      !(4)    plotting of the axis system
      !(5)    plotting the title
      !(6)    plotting data points
      !(7)    termination.

      !  Dislin Level 0:  before initialization or after termination
      call metafl(cfmt) ! set output driver/file-format (PNG); this is a 4-char string
      call setpag(cfsz) ! Set pagesize to US A Landscape (2790 x 2160)
      call setfil(trim(adjustl(filename_png))) ! Set output filename
      call scrmod('REVERSE')  ! Default background is black; reverse to white

      !  Dislin Level 1:  after initialization or a call to ENDGRF
      call disini()       ! initialize plot (set to level 1)
        ! Set the color table : SPEC,RAIN,GREY,TEMP
      call setvlt('RAIN')
      call paghdr('Ash3d Simulation plotted on ','---',4,0)

      if(HaveIconFile)then
        call filbox(x_leg3_px,y_footer_px,130,49)
        call incfil(trim(adjustl(Instit_IconFile)))
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Institution Logo not found."
          write(outlog(io),*)"If you would like your logo on these maps, copy a"
          write(outlog(io),*)"small image (130x50) to ASH3DHOME/share/post_proc/logo.png"
        endif;enddo
      endif

        ! setting of plot parameters
      !call pagera()       ! plot a border around the page
      !call triplx()  ! set font to triple stroke
      call bmpfnt('HELVE')

      !call helves()
      call titlin(title_plot,2)
      call name(cstr_xlabel,'X')
      call name(cstr_ylabel,'Y')
      call name(cstr_zlabel,'Z')
      call rvynam() ! reverse the axis labels

      call intax() ! With the routine INTAX, all axes will be labeled with integers.
      call autres(ntmax,nzmax) !The size of coloured rectangles will be automatically calculated by GRAF3 or CRVMAT
      call axspos(plotposx_px,plotposy_px)    ! determine the position of the axis system
      call ax3len(plotx_px,ploty_px,plotz_px) ! defines the size of the axis system : NXL, NYL, NZL

      call graf3(tmin,tmax,tlab1,tlabstep,& ! plots a 3-D axis system where the Z-axis
                 zmin,zmax,zlab1,zlabstep,& ! is plotted as a colour bar.
                 cmin,cmax,clab1,clabstep)
      call crvmat(conc,ntmax,nzmax,1,1) ! Interpolated data onto grid with spec. interp. points

      call height(50) ! Set character height for title
      call title() ! Actually write the title to the file

      call height(30) ! Reset character height to something smaller

      ! Add boxes below plot with run info
      call messag(cstr_volcname,x_leg1_px,y_footer_px+0*dy_newline_px)
      call messag(cstr_run_date,x_leg1_px,y_footer_px+1*dy_newline_px)
      call messag(cstr_windfile,x_leg1_px,y_footer_px+2*dy_newline_px)
      if(neruptions.gt.1)then
        call messag(cstr_note,x_leg1_px,y_footer_px+3*dy_newline_px)
      endif
      call messag(cstr_ErStartT,x_leg2_px,y_footer_px+0*dy_newline_px)
      call messag(cstr_ErHeight,x_leg2_px,y_footer_px+1*dy_newline_px)
      call messag(cstr_ErDuratn,x_leg2_px,y_footer_px+2*dy_newline_px)
      call messag(cstr_ErVolume,x_leg2_px,y_footer_px+3*dy_newline_px)

      !  Dislin Level 0:  before initialization or after termination
      call disfin()

      ! clean up memory
      if(allocated(t))      deallocate(t)
      if(allocated(z))      deallocate(z)
      if(allocated(conc))   deallocate(conc)

      end subroutine write_2Dprof_PNG_dislin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  write_DepPOI_TS_PNG_dislin
!
!  Called from: Ash3d_PostProc.f90
!  Arguments:
!    pt_indx       = index of point in Ash3d netcdf output file
!
!  This subroutine creates a png plot of the transient deposit accumulation at
!  the airport/POI given by pt_index using the dislin graphics package.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine write_DepPOI_TS_PNG_dislin(pt_indx)

      use Airports,      only : &
         Airport_Name,Airport_Thickness_TS

      use io_data,       only : &
         nWriteTimes,WriteTimes

      use time_data,     only : &
         Simtime_in_hours

      integer,intent(in) :: pt_indx

      character(len=14) :: filename_png
      integer,save      :: plot_index = 0

      real(kind=DS) :: xmin
      real(kind=DS) :: xmax
      real(kind=DS) :: ymin
      real(kind=DS) :: ymax
      real(kind=DS), dimension(:), allocatable :: x, y

      ! dislin stuff
      ! https://www.dislin.de/
      character(len=4) :: cfmt = "PNG "

      if(Airport_Thickness_TS(pt_indx,nWriteTimes).lt.0.01_ip)then
        return
      else
        plot_index = plot_index + 1
      endif

      write(filename_png,55) plot_index,".png"
 55   format('depTS_',i4.4,a4)

      if(Airport_Thickness_TS(plot_index,nWriteTimes).lt.0.01_ip)then
        ymax = 1.0_DS
      elseif(Airport_Thickness_TS(plot_index,nWriteTimes).lt.1.0_ip)then
        ymax = 1.0_DS
      elseif(Airport_Thickness_TS(plot_index,nWriteTimes).lt.5.0_ip)then
        ymax = 5.0_DS
      elseif(Airport_Thickness_TS(plot_index,nWriteTimes).lt.25.0_ip)then
        ymax = 25.0_DS
      else
        ymax = 100.0_DS
      endif

      xmin=0.0_DS
      xmax=real(ceiling(Simtime_in_hours),kind=DS)
      ymin=0.0_DS

      allocate(x(nWriteTimes))
      allocate(y(nWriteTimes))
      x = real(WriteTimes,kind=DS)
      y = real(Airport_Thickness_TS(pt_indx,1:nWriteTimes),kind=DS)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! DISLIN block
      ! https://www.dislin.de/
      ! wget
      ! https://ftp.gwdg.de/pub/grafik/dislin/linux/i586_64/dislin-11.4.linux.i586_64.tar.gz
      !(1)    setting of page format, file format and filename
      !(2)    initialization
      !(3)    setting of plot parameters
      !(4)    plotting of the axis system
      !(5)    plotting the title
      !(6)    plotting data points
      !(7)    termination.

      !  Dislin Level 0:  before initialization or after termination
      call metafl(cfmt)   ! set output driver/file-format (PNG); this is a 4-char string
      call winsiz(400, 300) ! Set window size to 400 x 300
      call sclmod("FULL")
      call setfil(trim(adjustl(filename_png))) ! Set output filename
      call scrmod('REVERSE')  ! Default background is black; reverse to white

      !  Dislin Level 1:  after initialization or a call to ENDGRF
      call disini()       ! initialize plot (set to level 1)
      call bmpfnt('HELVE')

        ! setting of plot parameters
      call name('Time (hours after eruption)','X') ! Set x-axis title
      call name('Deposit Thickeness (mm)','Y') ! Set y-axis title
      call labdig(-1,'X') ! set number of decimal places for x label (-1 means no decimal)
      call ticks(10,'XY') ! set number of ticks between labels
      call titlin(Airport_Name(pt_indx),4)  ! Set the title to the airport name (4 is the bottom line)

      !  Dislin Level 2: after a call to GRAF, GRAFP or GRAFMP
        ! Now create graph and set to level 2
      call graf(xmin, xmax, 0.0_DS, 5.0_DS, &
                ymin, ymax, 0.0_DS, 1.0_DS)
      call title() ! Actually write the title to the file
      call setrgb(0.5_DS, 0.5_DS, 0.5_DS)
      !call curve(real(x,kind=4),real(y,kind=4),nWriteTimes)  ! This draws the line
      call shdpat(16)  ! set shading pattern 
      call shdcrv(x,y,nWriteTimes,x,0.0_DS*y,nWriteTimes) ! This fills below curve
      call color('FORE') ! Reset color to defaul foreground color

      !  Dislin Level 0:  before initialization or after termination
      call disfin()

      ! clean up memory
      if(allocated(x))      deallocate(x)
      if(allocated(y))      deallocate(y)

      end subroutine write_DepPOI_TS_PNG_dislin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module Ash3d_PostProc_dislin

!##############################################################################
