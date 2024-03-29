!##############################################################################
!
! Ash3d_PostProc_dislin module
!
! This module provides the subroutines that use DISLIN for creating 2d maps,
! 2d vertical profiles, and the little deposit accumulation plots linked to
! the airport arrival kml (ash_arrivaltimes_airports.kml).  The dislin
! library is linked at compile-time.
!
!      subroutine write_2Dmap_PNG_dislin
!      subroutine write_2Dprof_PNG_dislin
!      subroutine write_DepPOI_TS_PNG_dislin
!
!##############################################################################

      module Ash3d_PostProc_dislin

      use precis_param

      use io_units

      use global_param,  only : &
         DirDelim

      use io_data,       only : &
         Ash3dHome

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
      character(100)    :: USGSIconFile

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
!    writeContours = logical
!
!  This subroutine creates a png map of the variable in OutVar using the dislin
!  graphics package.  Annotations and contour levels are indicated via the
!  product ID (iprod).  If writeContours is set to true, then this subroutine
!  is only used for generating and storing the contours (for dislin, this is
!  directly in memory as a variable) with no png written. If timestep = -1,
!  then use the last step in file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine write_2Dmap_PNG_dislin(nx,ny,iprod,itime,OutVar,writeContours)

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
         ContourDataX,ContourDataY,ContourDataNcurves,ContourDataNpoints,&
         Contour_MaxCurves,Contour_MaxPoints,ContourLev,nConLev

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
      logical      ,intent(in) :: writeContours

      integer :: i,j,k
      integer     , dimension(:,:),allocatable :: zrgb
      character(len=40) :: title_plot
      character(len=15) :: title_legend
      character(len=40) :: outfile_name
      character(len= 9) :: cio
      character(len= 4) :: outfile_ext = '.png'
      character(len=10) :: units
      integer :: ioerr,iw,iwf
      integer :: tmp_int

      real(kind=ip)  :: xmin
      real(kind=ip)  :: xmax
      real(kind=ip)  :: ymin
      real(kind=ip)  :: ymax
      real(kind=ip)  :: dx,dy

      ! dislin stuff
      ! https://www.dislin.de/
      character(len= 4) :: cfmt = "PNG "
      character(len=80) :: cbuf
      integer :: nmaxln  ! number of characters in the longest line of text
      character(len=7) :: zlevlab
      real(kind=DS) :: xpts(Contour_MaxPoints)
      real(kind=DS) :: ypts(Contour_MaxPoints)

      integer :: iray(Contour_MaxCurves)
      integer :: nxp,nyp,nclr,NCURVS

      real(kind=DS) :: xp,yp
      real(kind=DS) :: xminDIS, xmaxDIS, yminDIS, ymaxDIS
      real(kind=DS) :: dx_map, dy_map, xgrid_1, ygrid_1
      real(kind=DS) :: xr,xg,xb

      character(len=30) :: cstr_volcname
      character(len=30) :: cstr_run_date
      character(len=30) :: cstr_windfile
      character(len=40) :: cstr_ErStartT
      character(len=27) :: cstr_ErHeight
      character(len=30) :: cstr_ErDuratn
      character(len=38) :: cstr_ErVolume
      character(len=45) :: cstr_note
      integer :: y_footer
      logical :: UseShadedContours

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

      ncities = 20
      allocate(lon_cities(ncities))
      allocate(lat_cities(ncities))
      allocate(name_cities(ncities))

      ! logical switch to change from contour lines to shaded contours
      UseShadedContours = .false.
      !UseShadedContours = .true.

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
          ContourLev = (/0.1_ip, 0.3_ip, 1.0_ip, 3.0_ip, &
                  10.0_ip, 30.0_ip, 100.0_ip, 300.0_ip/)
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

      if(writeContours)then
        allocate(ContourDataNcurves(nConLev))
        allocate(ContourDataNpoints(nConLev,Contour_MaxCurves))
        allocate(ContourDataX(nConLev,Contour_MaxCurves,Contour_MaxPoints))
        allocate(ContourDataY(nConLev,Contour_MaxCurves,Contour_MaxPoints))
        ContourDataNcurves(:)   = 0
        ContourDataNpoints(:,:) = 0
        ContourDataX(:,:,:)     = 0.0_ip
        ContourDataY(:,:,:)     = 0.0_ip
        do i=1,nConLev
          ! This part calculates the contours
          xpts(1:Contour_MaxPoints) = 0.0_DS
          ypts(1:Contour_MaxPoints) = 0.0_DS
          iray(1:Contour_MaxCurves) = 0
          call conpts(real(lon_cc_pd(1:nx),kind=DS),nx,&  ! x coord and size
                      real(lat_cc_pd(1:ny),kind=DS),ny,&  ! y coord and size
                      real(OutVar(1:nx,1:ny),kind=DS), &  ! matrix with function values
                      real(ContourLev(i),kind=DS),     &  ! level to contour
                      xpts(1:Contour_MaxPoints),   &  ! x of contour (may have mul. curvex)
                      ypts(1:Contour_MaxPoints),   &  ! y of contour (may have mul. curves)
                      Contour_MaxPoints,           &  ! max # of points for contour arrays
                      iray(1:Contour_MaxCurves),   &  ! num of points for each contour
                      Contour_MaxCurves,           &  ! max number of curves
                      NCURVS)                         ! actual number of curves

          ContourDataNcurves(i)=NCURVS
          do j=1,NCURVS
            ContourDataNpoints(i,j)=iray(j)
            do k=1,iray(j)
              ContourDataX(i,j,k) = real(xpts(k),kind=ip)
              ContourDataY(i,j,k) = real(ypts(k),kind=ip)
            enddo
            ! These data could be plotted to produce the same plot as from contur
            !call curvmp(xpts,ypts,iray(j))
          enddo
        enddo
        ! Once we've loaded contours, we are all done here
        return
      endif

      ! This is the section where we actually start plotting the map
      if(IsLatLon)then
        xmin = minval(lon_cc_pd(1:nx))
        xmax = maxval(lon_cc_pd(1:nx))
        ymin = minval(lat_cc_pd(1:ny))
        ymax = maxval(lat_cc_pd(1:ny))
        dx = lon_cc_pd(2) - lon_cc_pd(1)
        dy = lat_cc_pd(2) - lat_cc_pd(1)
      else
        xmin = minval(x_cc_pd(1:nx))
        xmax = maxval(x_cc_pd(1:nx))
        ymin = minval(y_cc_pd(1:ny))
        ymax = maxval(y_cc_pd(1:ny))
        dx = x_cc_pd(2) - x_cc_pd(1)
        dy = y_cc_pd(2) - y_cc_pd(1)
      endif
      call citylist(0,xmin,xmax,ymin,ymax,    &
                    ncities,                  &
                    lon_cities, &
                    lat_cities, &
                    name_cities)

      dx_map = 10.0_DS
      dy_map = 5.0_DS
      lon_cc_pd(:) = lon_cc_pd(:) - 360.0_ip
      xminDIS = real(xmin- 360.0_ip-0.5_ip*dx,kind=DS)
      xmaxDIS = real(xmax- 360.0_ip+0.5_ip*dx,kind=DS)
      yminDIS = real(ymin-0.5_ip*dy,kind=DS)
      ymaxDIS = real(ymax+0.5_ip*dy,kind=DS)
      xgrid_1 = real(ceiling(xminDIS/dx_map) * dx_map,kind=DS)
      ygrid_1 = real(ceiling(yminDIS/dy_map) * dy_map,kind=DS)

      !!!!!!!!!!!!!!!!!!!!!!!
      !  Dislin Level 0:  before initialization or after termination
      call metafl(cfmt)   ! set output driver/file-format (PNG); this is a 4-char string
      call setpag('USAL') ! Set pagesize to US A Landscape (2790 x 2160)
      call setfil(trim(adjustl(outfile_name))) ! Set output filename
      call scrmod('REVERSE')  ! Default background is black; reverse to white

      !  Dislin Level 1:  after initialization or a call to ENDGRF
      call disini()       ! initialize plot (set to level 1)

        ! Set the color table : SPEC,RAIN,GREY,TEMP
      call setvlt('RAIN')

      call paghdr('Ash3d Simulation plotted on ','---',4,0)
      y_footer = 1900
      call filbox(2250,y_footer,130,49)
      USGSIconFile = trim(Ash3dHome) // &
                        DirDelim // 'share' // &
                        DirDelim // 'post_proc' // &
                        DirDelim // 'USGSvid.png'
      call incfil(USGSIconFile)

       ! setting of plot parameters
      call triplx()  ! set font to triple stroke
      call axspos(500,1650)  ! determine the position of the axis system
      call axslen(2200,1400) ! defines the size of the axis system
      call name('Longitude','X') ! Set x-axis title
      call name('Latitude','Y') ! Set y-axis title
      call labdig(1,'X') ! set number of decimal places for x label (-1 means no decimal)
      call labdig(1,'y') 
      call ticks(1,'xy')  ! set number of ticks between labels
      call titlin(title_plot,4)  ! Set the title
      call incmrk(0) ! selects line (0) or symbol (-1) mode for CURVE

      !call LABELS('MAP','xy')
      ! set projection : STER,LAMB,CYLI,MERC
      call projct('CYLI') ! defines projection
      call frame(3) ! bump up frame line thickness

      !  Dislin Level 2: after a call to GRAF, GRAFP or GRAFMP
        ! Now create graph and set to level 2
       !  The routine GRAFMP plots a geographical axis system.
      call grafmp(xminDIS,xmaxDIS,xgrid_1,dx_map, &
                  yminDIS,ymaxDIS,ygrid_1,dy_map)
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
      do i=1,ncities
      !    These are the points
        call pos2pt(real(lon_cities(i),kind=DS),real(lat_cities(i),kind=DS),&
                    xp,yp)
        nxp=nint(xp)
        nyp=nint(yp)
        call symbol(21,nxp,nyp)
        !    These are the city labels, offset in x
        call messag(adjustl(trim(name_cities(i))),nxp+30,nyp)
      enddo

      ! Add volcano
      call pos2pt(real(lon_volcano,kind=DS),real(lat_volcano,kind=DS),&
                  xp,yp)
      nxp=nint(xp)
      nyp=nint(yp)
      call symbol(18,nxp,nyp)

      if(UseShadedContours)then
        call shdmod('UPPER', 'CELL') ! This suppresses colors in regions above/below the zlevels pro
        call conshd(real(lon_cc_pd(1:nx),kind=DS),nx,&
                    real(lat_cc_pd(1:ny),kind=DS),ny,&
                    real(OutVar,kind=DS),real(ContourLev(i:nConLev),kind=DS),nConLev)
      else
        do i=1,nConLev
          xr   = real(zrgb(i,1),kind=DS)/real(255,kind=DS)
          xg   = real(zrgb(i,2),kind=DS)/real(255,kind=DS)
          xb   = real(zrgb(i,3),kind=DS)/real(255,kind=DS)
          nclr = intrgb(xr,xg,xb)

          call setclr(nclr)
          call contur(real(lon_cc_pd(1:nx),kind=DS),nx,&
                      real(lat_cc_pd(1:ny),kind=DS),ny,&
                      real(OutVar(1:nx,1:ny),kind=DS), &
                      real(ContourLev(i),kind=DS))
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
      call legini(cbuf,nConLev,nmaxln) ! Initialize legend
      call legtit(title_legend)      ! Set legend title
      call legbgd(0)                 ! sets background color
      do i=1,nConLev
        if(ContourLev(i).lt.1.0_ip)then
          write(zlevlab,'(f6.2)')real(ContourLev(i),kind=4)
        else
          write(zlevlab,'(f6.1)')real(ContourLev(i),kind=4)
        endif
        call leglin(cbuf,zlevlab,i)
      enddo
      call legend(cbuf,6)  ! write buffer to legend and give position (6=LR)

      ! Add boxes below plot with run info
      ! Volcano:     Erup.start:
      ! Run date:    Plm Height:
      ! Windfile:    Duration:
      !              Volume:

      write(cstr_volcname,'(a10,a20)')'Volcano:  ' ,VolcanoName
      write(cstr_run_date,'(a10,a20)')'Run Date: ',os_time_log
      read(cdf_b3l1,*,iostat=ioerr) iw,iwf
      write(cstr_windfile,'(a10,i5)')'Windfile: ',iwf

      !e_StartTime,e_PlumeHeight,e_Duration,e_Volume
      write(cstr_ErStartT,'(a20,a20)')'Erup. Start Time:   ',&
            HS_xmltime(SimStartHour+e_StartTime(1),BaseYear,useLeap)
      write(cstr_ErHeight,'(a20,f4.1,a3)')'Erup. Plume Height: ',e_PlumeHeight(1),' km'
      write(cstr_ErDuratn,'(a20,f4.1,a6)')'Erup. Duration:     ',e_Duration(1),' hours'
      write(cstr_ErVolume,'(a20,f8.5,a10)')'Erup. Volume:       ',e_Volume(1),' km3 (DRE)'

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
      call messag(cstr_ErVolume,1200,y_footer+120)

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
         nzmax,z_cc_pd

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

      integer, intent (in) :: vprof_ID

      character(len=14) :: dp_pngfile
      character(len=26) :: coord_str
      character(len=76) :: title_str
      integer :: k,i
      integer :: ioerr,iw,iwf

      real(kind=DS) :: tmin
      real(kind=DS) :: tmax
      real(kind=DS) :: tlab1
      real(kind=DS) :: tlabstep
      real(kind=DS) :: zmin
      real(kind=DS) :: zmax
      real(kind=DS) :: zlab1
      real(kind=DS) :: zlabstep
      real(kind=DS) :: cloudcon_thresh_mgm3
      real(kind=DS) :: cmin
      real(kind=DS) :: cmax
      real(kind=DS) :: clab1
      real(kind=DS) :: clabstep
      real(kind=DS), dimension(:),   allocatable :: t, z
      real(kind=DS), dimension(:,:), allocatable :: conc

      ! dislin stuff
      ! https://www.dislin.de/
      character(len=4)  :: cfmt = "PNG "
      character(len=30) :: cstr_volcname
      character(len=30) :: cstr_run_date
      character(len=30) :: cstr_windfile
      character(len=40) :: cstr_ErStartT
      character(len=27) :: cstr_ErHeight
      character(len=30) :: cstr_ErDuratn
      character(len=38) :: cstr_ErVolume
      character(len=45) :: cstr_note
      integer :: y_footer

      INTERFACE
        character (len=20) function HS_xmltime(HoursSince,byear,useLeaps)
          real(kind=8)              :: HoursSince
          integer                   :: byear
          logical                   :: useLeaps
        end function HS_xmltime
      END INTERFACE

      write(dp_pngfile,54) vprof_ID,".png"
 54   format('dslin_',i4.4,a4)

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
      cmax=real(max(cmax,cloudcon_thresh_mgm3),kind=DS)  ! Do not let cmax drop below the threshold

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
      else
          clabstep = 1.0e-1_DS
      endif
      clab1    = clabstep

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

      write(coord_str,101)x_vprofile(vprof_ID),y_vprofile(vprof_ID)
 101  format(' (lon=',f7.2,', lat=',f6.2,')')
      write(title_str,*)trim(adjustl(Site_vprofile(vprof_ID))),coord_str

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
      call setpag('USAL') ! Set pagesize to US A Landscape (2790 x 2160)
      call setfil(trim(adjustl(dp_pngfile))) ! Set output filename
      call scrmod('REVERSE')  ! Default background is black; reverse to white

      !  Dislin Level 1:  after initialization or a call to ENDGRF
      call disini()       ! initialize plot (set to level 1)
      y_footer = 1900
      call filbox(2250,y_footer,130,49)
      USGSIconFile = trim(Ash3dHome) // &
                        DirDelim // 'share' // &
                        DirDelim // 'post_proc' // &
                        DirDelim // 'USGSvid.png'
      call incfil(USGSIconFile)

        ! setting of plot parameters
      !call pagera()       ! plot a border around the page
      call helves()
      call titlin(title_str,2)
      call name('Time (hours after eruption)','X')
      call name('Height (km)','Y')
      call name('Ash conc. mg/m3','Z')
      call rvynam() ! reverse the axis labels

      call intax() ! With the routine INTAX, all axes will be labeled with integers.
      call autres(ntmax,nzmax) !The size of coloured rectangles will be automatically calculated by GRAF3 or CRVMAT
      call axspos(500,1650)  ! determine the position of the axis system
      call ax3len(1800,1200,1200) ! defines the size of the axis system : NXL, NYL, NZL

      call graf3(tmin,tmax,tlab1,tlabstep,& !plots a 3-D axis system where the Z-axis
                 zmin,zmax,zlab1,zlabstep,& ! is plotted as a colour bar.
                 cmin,cmax,clab1,clabstep)
      call crvmat(conc,ntmax,nzmax,1,1) ! Interpolated data onto grid with spec. interp. points

      call height(50) ! Set character height for title
      call title() ! Actually write the title to the file

      call height(30) ! Reset character height to something smaller
      write(cstr_volcname,'(a10,a20)')'Volcano:  ' ,VolcanoName
      write(cstr_run_date,'(a10,a20)')'Run Date: ',os_time_log
      read(cdf_b3l1,*,iostat=ioerr) iw,iwf
      write(cstr_windfile,'(a10,i5)')'Windfile: ',iwf

      !e_StartTime,e_PlumeHeight,e_Duration,e_Volume
      write(cstr_ErStartT,'(a20,a20)')'Erup. Start Time:   ',HS_xmltime(SimStartHour+e_StartTime(1),BaseYear,useLeap)
      write(cstr_ErHeight,'(a20,f4.1,a3)')'Erup. Plume Height: ',e_PlumeHeight(1),' km'
      write(cstr_ErDuratn,'(a20,f4.1,a6)')'Erup. Duration:     ',e_Duration(1),' hours'
      write(cstr_ErVolume,'(a20,f8.5,a10)')'Erup. Volume:       ',e_Volume(1),' km3 (DRE)'

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
      call messag(cstr_ErVolume,1200,y_footer+120)

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

      character(len=14) :: dp_pngfile
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

      write(dp_pngfile,55) plot_index,".png"
 55   format('dslin_',i4.4,a4)

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
      call setpag('USAL') ! Set pagesize to US A Landscape (2790 x 2160)
      call setfil(trim(adjustl(dp_pngfile))) ! Set output filename
      call scrmod('REVERSE')  ! Default background is black; reverse to white

      !  Dislin Level 1:  after initialization or a call to ENDGRF
      call disini()       ! initialize plot (set to level 1)
        ! setting of plot parameters
      call pagera()       ! plot a border around the page
      call triplx()  ! set font to triple stroke
      call axspos(450,1800)  ! determine the position of the axis system
      call axslen(2200,1200) ! defines the size of the axis system
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
