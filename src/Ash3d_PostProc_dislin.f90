!##############################################################################
!
!    write_2Dmap_PNG_dislin
!
!    if timestep = -1, then use the last step in file
!##############################################################################

      subroutine write_2Dmap_PNG_dislin(iprod,itime,OutVar,writeContours)

      use precis_param

      use mesh,          only : &
         nxmax,nymax,x_cc_pd,y_cc_pd,lon_cc_pd,lat_cc_pd, &
         IsLatLon

      use Output_Vars,   only : &
         DepositThickness,DepArrivalTime,CloudArrivalTime,&
         MaxConcentration,MaxHeight,CloudLoad,dbZ,MinHeight,Mask_Cloud,Mask_Deposit,&
         Con_DepThick_mm_N,Con_DepThick_mm_Lev,Con_DepThick_mm_RGB, &
         Con_DepThick_in_N,Con_DepThick_in_Lev,Con_DepThick_in_RGB, &
         Con_DepTime_N,Con_DepTime_Lev,Con_DepTime_RGB, &
         Con_CloudCon_N,Con_CloudCon_Lev,Con_CloudCon_RGB, &
         Con_CloudTop_N,Con_CloudTop_RGB,Con_CloudTop_Lev, &
         Con_CloudBot_N,Con_CloudBot_RGB,Con_CloudBot_Lev, &
         Con_CloudLoad_N,Con_CloudLoad_RGB,Con_CloudLoad_Lev, &
         Con_CloudRef_N,Con_CloudRef_RGB,Con_CloudRef_Lev, &
         Con_CloudTime_N,Con_CloudTime_RGB,Con_CloudTime_Lev,&
         ContourDataX,ContourDataY,ContourDataNcurves,ContourDataNpoints,&
         Contour_MaxCurves,Contour_MaxPoints,ContourLev

      use io_data,       only : &
         nWriteTimes,WriteTimes,cdf_b3l1,VolcanoName

      use Source,        only : &
         neruptions,e_Volume,e_Duration,e_StartTime,e_PlumeHeight,lon_volcano,lat_volcano

      use time_data,     only : &
         os_time_log,BaseYear,useLeap

      use dislin

      implicit none

      integer :: iprod
      integer :: itime
      real(kind=ip) :: OutVar(nxmax,nymax)
      logical :: writeContours

      integer :: i,j,k
      integer :: nzlev
      real(kind=4), dimension(:)  ,allocatable :: zlev
      integer     , dimension(:,:),allocatable :: zrgb
      character(len=40) :: title_plot
      character(len=15) :: title_legend
      character(len=40) :: outfile_name
      character (len=9) :: cio
      character (len=4) :: outfile_ext = '.png'
      integer :: ioerr,iw,iwf

      real(kind=8)  :: xmin
      real(kind=8)  :: xmax
      real(kind=8)  :: ymin
      real(kind=8)  :: ymax

      ! dislin stuff
      ! https://www.dislin.de/
      CHARACTER (LEN=4) :: CFMT = "PNG "
      CHARACTER(len=80) :: CBUF
      integer :: nmaxln  ! number of characters in the longest line of text
      character(len=7) :: zlevlab
      real(kind=4) :: XPTS(Contour_MaxPoints),YPTS(Contour_MaxPoints)

      INTEGER :: IRAY(Contour_MaxCurves)
      INTEGER :: NXP,NYP,nclr,NCURVS

      REAL(kind=4) :: XP,YP
      real(kind=4) :: xminDIS, xmaxDIS, yminDIS, ymaxDIS
      real(kind=4) :: dx_map, dy_map, xgrid_1, ygrid_1
      ! Hot colormap RGB's and breakpoints
      real(kind=4), dimension(4) :: cpt_break_hot = (/ 0.0_4, 0.38_4, 0.76_4, 1.0_4 /)
      real(kind=4), dimension(4) :: cpt_r_hot     = (/ 0.0_4, 1.0_4, 1.0_4, 1.0_4 /)
      real(kind=4), dimension(4) :: cpt_g_hot     = (/ 0.0_4, 0.0_4, 1.0_4, 1.0_4 /)
      real(kind=4), dimension(4) :: cpt_b_hot     = (/ 0.0_4, 0.0_4, 0.0_4, 0.9_4 /) ! hot actually ends in 1.0
      real(kind=4) :: i_flt,idel,cdel,xr,xg,xb

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
      real(kind=8),dimension(:),allocatable     :: lon_cities
      real(kind=8),dimension(:),allocatable     :: lat_cities
      character(len=26),dimension(:),allocatable :: name_cities

      INTERFACE
        character (len=20) function HS_xmltime(HoursSince,byear,useLeaps)
          real(kind=8)              :: HoursSince
          integer                   :: byear
          logical                   :: useLeaps
        end function HS_xmltime
        subroutine citylist(outCode,lonLL,lonUR,latLL,latUR,ncities, &
                            CityLon_out,CityLat_out,CityName_out)
          integer      :: outCode
          real(kind=8) :: lonLL
          real(kind=8) :: lonUR
          real(kind=8) :: latLL
          real(kind=8) :: latUR
          integer      :: ncities

          real(kind=8),dimension(ncities) :: CityLon_out
          real(kind=8),dimension(ncities) :: CityLat_out
          character(len=26),dimension(ncities) :: CityName_out
        end subroutine citylist
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

      allocate(ContourLev(Contour_MaxCurves))
      if(iprod.eq.3)then       ! deposit at specified times (mm)
        write(outfile_name,'(a15,a9,a4)')'Ash3d_Deposit_t',cio,outfile_ext
        write(title_plot,'(a20,f5.2,a6)')'Deposit Thickness t=',WriteTimes(itime),' hours'
        title_legend = 'Dep.Thick.(mm)'
        nzlev = Con_DepThick_mm_N
        allocate(zlev(nzlev))
        allocate(zrgb(nzlev,3))
        ContourLev(1:nzlev) = Con_DepThick_mm_Lev(1:nzlev)
        zlev(1:nzlev) = real(Con_DepThick_mm_Lev(1:nzlev),kind=4)
        zrgb(1:nzlev,1:3) = Con_DepThick_mm_RGB(1:nzlev,1:3)
      elseif(iprod.eq.4)then   ! deposit at specified times (inches)
        write(outfile_name,'(a15,a9,a4)')'Ash3d_Deposit_t',cio,outfile_ext
        write(title_plot,'(a20,f5.2,a6)')'Deposit Thickness t=',WriteTimes(itime),' hours'
        title_legend = 'Dep.Thick.(in)'
        nzlev = Con_DepThick_in_N
        allocate(zlev(nzlev))
        allocate(zrgb(nzlev,3))
        ContourLev(1:nzlev) = Con_DepThick_in_Lev(1:nzlev)
        zlev(1:nzlev) = real(Con_DepThick_in_Lev(1:nzlev),kind=4)
        zrgb(1:nzlev,1:3) = Con_DepThick_in_RGB(1:nzlev,1:3)
      elseif(iprod.eq.5)then       ! deposit at final time (mm)
        write(outfile_name,'(a13,a9,a4)')'Ash3d_Deposit',cio,outfile_ext
        title_plot = 'Final Deposit Thickness'
        title_legend = 'Dep.Thick.(mm)'
        nzlev = Con_DepThick_mm_N
        allocate(zlev(nzlev))
        allocate(zrgb(nzlev,3))
        ContourLev(1:nzlev) = Con_DepThick_mm_Lev(1:nzlev)
        zlev(1:nzlev) = real(Con_DepThick_mm_Lev(1:nzlev),kind=4)
        zrgb(1:nzlev,1:3) = Con_DepThick_mm_RGB(1:nzlev,1:3)
      elseif(iprod.eq.6)then   ! deposit at final time (inches)
        write(outfile_name,'(a13,a9,a4)')'Ash3d_Deposit',cio,outfile_ext
        title_plot = 'Final Deposit Thickness'
        title_legend = 'Dep.Thick.(in)'
        nzlev = Con_DepThick_in_N
        allocate(zlev(nzlev))
        allocate(zrgb(nzlev,3))
        ContourLev(1:nzlev) = Con_DepThick_in_Lev(1:nzlev)
        zlev(1:nzlev) = real(Con_DepThick_in_Lev(1:nzlev),kind=4)
        zrgb(1:nzlev,1:3) = Con_DepThick_in_RGB(1:nzlev,1:3)
      elseif(iprod.eq.7)then   ! ashfall arrival time (hours)
        write(outfile_name,'(a22)')'DepositArrivalTime.png'
        write(title_plot,'(a20)')'Ashfall arrival time'
        title_legend = 'Time (hours)'
        nzlev = Con_DepTime_N
        allocate(zlev(nzlev))
        allocate(zrgb(nzlev,3))
        ContourLev(1:nzlev) = Con_DepTime_Lev(1:nzlev)
        zlev(1:nzlev) = real(Con_DepTime_Lev(1:nzlev),kind=4)
        zrgb(1:nzlev,1:3) = Con_DepTime_RGB(1:nzlev,1:3)
      elseif(iprod.eq.8)then   ! ashfall arrival at airports/POI (mm)
        write(*,*)"ERROR: No map PNG output option for airport arrival time data."
        write(*,*)"       Should not be in write_2Dmap_PNG_dislin"
        stop 1
      elseif(iprod.eq.9)then   ! ash-cloud concentration
        write(outfile_name,'(a16,a9,a4)')'Ash3d_CloudCon_t',cio,outfile_ext
        write(title_plot,'(a26,f5.2,a6)')'Ash-cloud concentration t=',WriteTimes(itime),' hours'
        title_legend = 'Max.Con.(mg/m3)'
        nzlev = Con_CloudCon_N
        allocate(zlev(nzlev))  
        allocate(zrgb(nzlev,3))
        ContourLev(1:nzlev) = Con_CloudCon_Lev(1:nzlev)
        zlev(1:nzlev) = real(Con_CloudCon_Lev(1:nzlev),kind=4)
        zrgb(1:nzlev,1:3) = Con_CloudCon_RGB(1:nzlev,1:3)
      elseif(iprod.eq.10)then   ! ash-cloud height
        write(outfile_name,'(a19,a9,a4)')'Ash3d_CloudHeight_t',cio,outfile_ext
        write(title_plot,'(a19,f5.2,a6)')'Ash-cloud height t=',WriteTimes(itime),' hours'
        title_legend = 'Cld.Height(km)'
        nzlev = Con_CloudTop_N
        allocate(zlev(nzlev))
        allocate(zrgb(nzlev,3))
        ContourLev(1:nzlev) = Con_CloudTop_Lev(1:nzlev)
        zlev(1:nzlev) = real(Con_CloudTop_Lev(1:nzlev),kind=4)
        zrgb(1:nzlev,1:3) = Con_CloudTop_RGB(1:nzlev,1:3)
      elseif(iprod.eq.11)then   ! ash-cloud bottom
        write(outfile_name,'(a16,a9,a4)')'Ash3d_CloudBot_t',cio,outfile_ext
        write(title_plot,'(a19,f5.2,a6)')'Ash-cloud bottom t=',WriteTimes(itime),' hours'
        title_legend = 'Cld.Bot.(km)'
        nzlev = Con_CloudBot_N
        allocate(zlev(nzlev))
        allocate(zrgb(nzlev,3))
        ContourLev(1:nzlev) = Con_CloudBot_Lev(1:nzlev)
        zlev(1:nzlev) = real(Con_CloudBot_Lev(1:nzlev),kind=4)
        zrgb(1:nzlev,1:3) = Con_CloudBot_RGB(1:nzlev,1:3)
      elseif(iprod.eq.12)then   ! ash-cloud load
        write(outfile_name,'(a17,a9,a4)')'Ash3d_CloudLoad_t',cio,outfile_ext
        write(title_plot,'(a17,f5.2,a6)')'Ash-cloud load t=',WriteTimes(itime),' hours'
        title_legend = 'Cld.Load(T/km2)'
        nzlev = Con_CloudLoad_N
        allocate(zlev(nzlev))
        allocate(zrgb(nzlev,3))
        ContourLev(1:nzlev) = Con_CloudLoad_Lev(1:nzlev)
        zlev(1:nzlev) = real(Con_CloudLoad_Lev(1:nzlev),kind=4)
        zrgb(1:nzlev,1:3) = Con_CloudLoad_RGB(1:nzlev,1:3)
      elseif(iprod.eq.13)then  ! radar reflectivity
        write(outfile_name,'(a20,a9,a4)')'Ash3d_CloudRadRefl_t',cio,outfile_ext
        write(title_plot,'(a24,f5.2,a6)')'Ash-cloud radar refl. t=',WriteTimes(itime),' hours'
        title_legend = 'Cld.Refl.(dBz)'
        nzlev = Con_CloudRef_N
        allocate(zlev(nzlev))
        allocate(zrgb(nzlev,3))
        ContourLev(1:nzlev) = Con_CloudRef_Lev(1:nzlev)
        zlev(1:nzlev) = real(Con_CloudRef_Lev(1:nzlev),kind=4)
        zrgb(1:nzlev,1:3) = Con_CloudRef_RGB(1:nzlev,1:3)
      elseif(iprod.eq.14)then   ! ashcloud arrival time (hours)
        write(outfile_name,'(a20)')'CloudArrivalTime.png'
        write(title_plot,'(a22)')'Ash-cloud arrival time'
        title_legend = 'Time (hours)'
        nzlev = Con_CloudTime_N
        allocate(zlev(nzlev))
        allocate(zrgb(nzlev,3))
        ContourLev(1:nzlev) = Con_CloudTime_Lev(1:nzlev)
        zlev(1:nzlev) = real(Con_CloudTime_Lev(1:nzlev),kind=4)
        zrgb(1:nzlev,1:3) = Con_CloudTime_RGB(1:nzlev,1:3)
      elseif(iprod.eq.15)then   ! topography
        write(outfile_name,'(a14)')'Topography.png'
        write(title_plot,'(a10)')'Topography'
        title_legend = 'Elevation (km)'
        nzlev = 8
        zlev = (/0.1_4, 0.3_4, 1.0_4, 3.0_4, &
                10.0_4, 30.0_4, 100.0_4, 300.0_4/)
      elseif(iprod.eq.16)then   ! profile plots
        write(*,*)"ERROR: No map PNG output option for vertical profile data."
        write(*,*)"       Should not be in write_2Dmap_PNG_dislin"
        stop 1
      else
        write(*,*)"ERROR: unexpected variable"
        stop 1
      endif

      if(writeContours)then
        allocate(ContourDataNcurves(Contour_MaxCurves))
        allocate(ContourDataNpoints(Contour_MaxCurves,Contour_MaxPoints))
        allocate(ContourDataX(nzlev,Contour_MaxCurves,Contour_MaxPoints))
        allocate(ContourDataY(nzlev,Contour_MaxCurves,Contour_MaxPoints))
        ContourDataNcurves(:)   = 0
        ContourDataNpoints(:,:) = 0
        ContourDataX(:,:,:)     = 0.0_8
        ContourDataY(:,:,:)     = 0.0_8
        do i=1,nzlev
          ! This part calculated the contours
          XPTS(:) = 0.0_4
          YPTS(:) = 0.0_4
          IRAY(:) = 0
          call conpts(real(lon_cc_pd(1:nxmax),kind=4),nxmax,&  ! x coord and size
                      real(lat_cc_pd(1:nymax),kind=4),nymax,&  ! y coord and size
                      real(OutVar(1:nxmax,1:nymax),kind=4),&   ! matrix with function values
                      zlev(i),            &           ! level to contour
                      XPTS,YPTS,          &           ! x,y of contour (may have mul. curves)
                      Contour_MaxPoints,  &           ! max # of points for contour arrays
                      IRAY,               &           ! num of points for each contour
                      Contour_MaxCurves,  &           ! max number of curves
                      NCURVS)                         ! actual number of curves
          ContourDataNcurves(i)=NCURVS
          do j=1,NCURVS
            ContourDataNpoints(i,j)=IRAY(j)
            do k=1,IRAY(j)
              ContourDataX(i,j,k) = real(XPTS(k),kind=8)
              ContourDataY(i,j,k) = real(YPTS(k),kind=8)
              !write(*,*) i,j,k,XPTS(k),YPTS(k),zlev(i)
            enddo

            !call curvmp(XPTS,YPTS,IRAY(j))
          enddo
        enddo
        return
      endif

      ! Now map plots
      xmin = real(minval(lon_cc_pd(1:nxmax)),kind=8)
      xmax = real(maxval(lon_cc_pd(1:nxmax)),kind=8)
      ymin = real(minval(lat_cc_pd(1:nymax)),kind=8)
      ymax = real(maxval(lat_cc_pd(1:nymax)),kind=8)

      dx_map = 10.0_4
      dy_map = 5.0_4
      lon_cc_pd(:) = lon_cc_pd(:) - 360.0_8
      xminDIS = real(minval(lon_cc_pd(1:nxmax)),kind=4)
      xmaxDIS = real(maxval(lon_cc_pd(1:nxmax)),kind=4)
      yminDIS = real(minval(lat_cc_pd(1:nymax)),kind=4)
      ymaxDIS = real(maxval(lat_cc_pd(1:nymax)),kind=4)
      xgrid_1 = real(ceiling(xminDIS/dx_map) * dx_map,kind=4)
      ygrid_1 = real(ceiling(yminDIS/dy_map) * dy_map,kind=4)

      call citylist(0,xmin,xmax,ymin,ymax, &
                    ncities,                        &
                    lon_cities,lat_cities,          &
                    name_cities)

      !!!!!!!!!!!!!!!!!!!!!!!
      !  Dislin Level 0:  before initialization or after termination
      call metafl(CFMT)   ! set output driver/file-format (PNG); this is a 4-char string
      call setpag('USAL') ! Set pagesize to US A Landscape (2790 x 2160)
      call setfil(trim(adjustl(outfile_name))) ! Set output filename
      call scrmod('REVERSE')  ! Default background is black; reverse to white

      !  Dislin Level 1:  after initialization or a call to ENDGRF
      call disini()       ! initialize plot (set to level 1)


        ! Set the color table
      !call setvlt('SPEC')
      call setvlt('RAIN')
      !call setvlt('GREY')
      !call setvlt('TEMP')

      call paghdr('Ash3d Simulation plotted on ','---',4,0)
      y_footer = 1900
      call filbox(2250,y_footer,130,49)
      call incfil('/opt/USGS/Ash3d/share/post_proc/USGSvid.png')
        ! setting of plot parameters
      !call pagera()       ! plot a border around the page
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
      !call projct('STER') ! defines projection
      !call projct('LAMB') ! defines projection
      call projct('CYLI') ! defines projection
      !call projct('MERC') ! defines projection
      call frame(3) ! bump up frame line thickness
       !  The routine GRAFMP plots a geographical axis system.
      call grafmp(xminDIS,xmaxDIS,xgrid_1,dx_map, &
                  yminDIS,ymaxDIS,ygrid_1,dy_map)

       ! set color of coastlines
      !call color('GREEN')
       ! plots coastlines and lakes or political borders
      call world()
      !call shdmap('GSHH')

      ! Add cities
      do i=1,ncities
      !    These are the points
        call pos2pt(real(lon_cities(i),kind=4),real(lat_cities(i),kind=4),&
                    XP,YP)
        NXP=NINT(XP)
        NYP=NINT(YP)
        call symbol(21,NXP,NYP)
      !    These are the city labels, offset in x
        call messag(adjustl(trim(name_cities(i))),NXP+30,NYP)
      enddo

      ! Add volcano
      call pos2pt(real(lon_volcano,kind=4),real(lat_volcano,kind=4),&
                  XP,YP)
      NXP=NINT(XP)
      NYP=NINT(YP)
      call symbol(18,NXP,NYP)

      if(UseShadedContours)then
        !call myvlt(xr,xg,xb,nrgb)
        call shdmod('UPPER', 'CELL') ! This suppresses colors in regions above/below the zlevels pro
        call conshd(real(lon_cc_pd(1:nxmax),kind=4),nxmax,&
                    real(lat_cc_pd(1:nymax),kind=4),nymax,&
                    real(OutVar,kind=4),zlev,nzlev)
      else
        do i=1,nzlev
        !do i=1,1
          xr   = real(zrgb(i,1),kind=4)/real(255,kind=4)
          xg   = real(zrgb(i,2),kind=4)/real(255,kind=4)
          xb   = real(zrgb(i,3),kind=4)/real(255,kind=4)
          nclr = intrgb(xr,xg,xb)

          call setclr(nclr)
          !if(writeContours)then
          !  ! This part calculated the contours
          !  XPTS(:) = 0.0_4
          !  YPTS(:) = 0.0_4
          !  IRAY(:) = 0
          !  call conpts(real(lon_cc_pd(1:nxmax),kind=4),nxmax,&  ! x coord and size
          !              real(lat_cc_pd(1:nymax),kind=4),nymax,&  ! y coord and size
          !              real(OutVar(1:nxmax,1:nymax),kind=4),&   ! matrix with function values
          !              zlev(i),            &           ! level to contour
          !              XPTS,YPTS,          &           ! x,y of contour (may have mul. curves)
          !              Contour_MaxPoints,  &           ! max # of points for contour arrays
          !              IRAY,               &           ! num of points for each contour
          !              Contour_MaxCurves,  &           ! max number of curves
          !              NCURVS)                         ! actual number of curves
          !  ContourDataNcurves(i)=NCURVS
          !  do j=1,NCURVS
          !    ContourDataNpoints(i,j)=IRAY(j)
          !    do k=1,IRAY(j)
          !      ContourDataX(i,j,k) = XPTS(k)
          !      ContourDataY(i,j,k) = YPTS(k)
          !      write(*,*) i,j,k,XPTS(k),YPTS(k),zlev(i)
          !    enddo

          !    !call curvmp(XPTS,YPTS,IRAY(j))
          !  enddo
          !else
            ! If no contour data are needed, just call contur
            call contur(real(lon_cc_pd(1:nxmax),kind=4),nxmax,&
                        real(lat_cc_pd(1:nymax),kind=4),nymax,&
                        real(OutVar(1:nxmax,1:nymax),kind=4),zlev(i))
          !endif
        enddo
      endif

       ! set color of grid lines
      call setrgb(0.0_4, 0.0_4, 0.0_4)
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
      write(cstr_ErStartT,'(a20,a20)')'Erup. Start Time:   ',HS_xmltime(e_StartTime(1),BaseYear,useLeap)
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

      end subroutine write_2Dmap_PNG_dislin

!##############################################################################
!
!    write_2Dprof_PNG_dislin
!
!##############################################################################

      subroutine write_2Dprof_PNG_dislin(vprof_ID)

      use precis_param

      use mesh,          only : &
         nzmax,z_cc_pd

      use Output_Vars,   only : &
         pr_ash

      use time_data,     only : &
         ntmax,time_native

      use io_data,       only : &
         Site_vprofile,x_vprofile,y_vprofile,cdf_b3l1,VolcanoName

      use Source,        only : &
         neruptions,e_Volume,e_Duration,e_StartTime,e_PlumeHeight

      use time_data,     only : &
         os_time_log,BaseYear,useLeap

      use dislin

      implicit none

      integer, intent (in) :: vprof_ID

      character(len=14) :: dp_pngfile
      character(len=26) :: coord_str
      character(len=76) :: title_str
      character(len=80) :: outstring
      integer :: k,i
      integer :: ioerr,iw,iwf

      real(kind=4) :: tmin
      real(kind=4) :: tmax
      real(kind=4) :: tlab1
      real(kind=4) :: tlabstep
      real(kind=4) :: zmin
      real(kind=4) :: zmax
      real(kind=4) :: zlab1
      real(kind=4) :: zlabstep
      real(kind=4) :: cmin
      real(kind=4) :: cmax
      real(kind=4) :: clab1
      real(kind=4) :: clabstep
      real(kind=4), dimension(:),   allocatable :: t, z
      real(kind=4), dimension(:,:), allocatable :: conc

      ! dislin stuff
      ! https://www.dislin.de/
      CHARACTER (LEN=4) :: CFMT = "PNG "
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
      tmin=real(0,kind=4)
      tmax=real(ceiling(time_native(ntmax)),kind=4)
      tlab1    = 0.0_4
      if(tmax.gt.240.0_4)then
        tlabstep = 48.0_4
      elseif(tmax.gt.120.0_4)then
        tlabstep = 24.0_4
      elseif(tmax.gt.30.0_4)then
        tlabstep = 10.0_4
      elseif(tmax.gt.15.0_4)then
        tlabstep = 5.0_4
      elseif(tmax.gt.6.0_4)then
        tlabstep = 2.0_4
      else
        tlabstep = 1.0_4
      endif

      zmin=real(0,kind=4)
      zmax=real(z_cc_pd(nzmax),kind=4)
      zlab1    = 0.0_4
      if(zmax.gt.30.0_4)then
        zlabstep = 10.0_4
      elseif(zmax.gt.15.0_4)then
        zlabstep = 5.0_4
      elseif(zmax.gt.6.0_4)then
        zlabstep = 2.0_4
      else
        zlabstep = 1.0_4
      endif

      cmin=real(0,kind=4)
      cmax=real(maxval(pr_ash(:,:,vprof_ID)),kind=4)
      if    (cmax.gt.4.0e4_4)then
          clabstep = 5.0e3_4
      elseif(cmax.gt.1.0e4_4)then
          clabstep = 2.0e3_4
      elseif(cmax.gt.4.0e3_4)then
          clabstep = 5.0e2_4
      elseif(cmax.gt.1.0e3_4)then
          clabstep = 2.0e2_4
      elseif(cmax.gt.4.0e2_4)then
          clabstep = 5.0e1_4
      elseif(cmax.gt.1.0e2_4)then
          clabstep = 2.0e1_4
      elseif(cmax.gt.4.0e1_4)then
          clabstep = 5.0e0_4
      elseif(cmax.gt.1.0e1_4)then
          clabstep = 2.0e0_4
      else
          clabstep = 1.0e-1_4
      endif
      clab1    = clabstep

      allocate(t(ntmax))
      allocate(z(nzmax))
      allocate(conc(ntmax,nzmax))
      t = real(time_native(1:ntmax),kind=4)
      z = real(z_cc_pd(1:nzmax),kind=4)
      do i=1,ntmax
        do k=1,nzmax
          conc(i,k) = real(pr_ash(k,i,vprof_ID),kind=4)
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
!     (1)    setting of page format, file format and filename
!     (2)    initialization
!     (3)    setting of plot parameters
!     (4)    plotting of the axis system
!     (5)    plotting the title
!     (6)    plotting data points
!     (7)    termination.

      !  Dislin Level 0:  before initialization or after termination
      call metafl(CFMT)   ! set output driver/file-format (PNG); this is a 4-char string
      call setpag('USAL') ! Set pagesize to US A Landscape (2790 x 2160)
      call setfil(trim(adjustl(dp_pngfile))) ! Set output filename
      call scrmod('REVERSE')  ! Default background is black; reverse to white

      !  Dislin Level 1:  after initialization or a call to ENDGRF
      call disini()       ! initialize plot (set to level 1)
      y_footer = 1900
      call filbox(2250,y_footer,130,49)
      call incfil('/opt/USGS/Ash3d/share/post_proc/USGSvid.png')

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
      write(cstr_ErStartT,'(a20,a20)')'Erup. Start Time:   ',HS_xmltime(e_StartTime(1),BaseYear,useLeap)
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

      end subroutine write_2Dprof_PNG_dislin

!##############################################################################
!
!    write_DepPOI_TS_PNG_dislin
!
!##############################################################################

      subroutine write_DepPOI_TS_PNG_dislin(pt_indx)

      use precis_param

      use Airports,      only : &
         nairports,Airport_Code,Airport_Name,Airport_x,Airport_y,&
         Airport_Latitude,Airport_Longitude,Airport_Thickness_TS

      use io_data,       only : &
         nWriteTimes,WriteTimes,VolcanoName

      use time_data,     only : &
         Simtime_in_hours

      use dislin

      implicit none

      integer :: pt_indx,i

      real(kind=4) :: ymaxpl
      character(len=14) :: dp_pngfile
      character(len=25) :: gnucom
      integer,save      :: plot_index = 0

      real(kind=4) :: xmin
      real(kind=4) :: xmax
      real(kind=4) :: ymin
      real(kind=4) :: ymax
      real(kind=4), dimension(:), allocatable :: x, y

      ! dislin stuff
      ! https://www.dislin.de/
      CHARACTER (LEN=4) :: CFMT = "PNG "

      if(Airport_Thickness_TS(pt_indx,nWriteTimes).lt.0.01_ip)then
        !write(*,*)"No deposit at ",pt_indx,Airport_Name(pt_indx)
        return
      else
        plot_index = plot_index + 1
        !write(*,*)"Processing ",pt_indx,plot_index,Airport_Name(pt_indx),'--'
      endif

      write(dp_pngfile,55) plot_index,".png"
 55   format('dslin_',i4.4,a4)

      if(Airport_Thickness_TS(plot_index,nWriteTimes).lt.0.01_4)then
        ymaxpl = 1.0_4
      elseif(Airport_Thickness_TS(plot_index,nWriteTimes).lt.1.0_4)then
        ymaxpl = 1.0_4
      elseif(Airport_Thickness_TS(plot_index,nWriteTimes).lt.5.0_4)then
        ymaxpl = 5.0_4
      elseif(Airport_Thickness_TS(plot_index,nWriteTimes).lt.25.0_4)then
        ymaxpl = 25.0_4
      else
        ymaxpl = 100.0_4
      endif

      xmin=0.0_4
      xmax=real(ceiling(Simtime_in_hours),kind=4)
      ymin=0.0_4
      ymax=ymaxpl

      allocate(x(nWriteTimes))
      allocate(y(nWriteTimes))
      x = real(WriteTimes,kind=4)
      y = real(Airport_Thickness_TS(pt_indx,1:nWriteTimes),kind=4)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! DISLIN block
      ! https://www.dislin.de/
      ! wget
      ! https://ftp.gwdg.de/pub/grafik/dislin/linux/i586_64/dislin-11.4.linux.i586_64.tar.gz
!     (1)    setting of page format, file format and filename
!     (2)    initialization
!     (3)    setting of plot parameters
!     (4)    plotting of the axis system
!     (5)    plotting the title
!     (6)    plotting data points
!     (7)    termination.

      !  Dislin Level 0:  before initialization or after termination
      call metafl(CFMT)   ! set output driver/file-format (PNG); this is a 4-char string
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
      call graf(real(xmin,kind=4), real(xmax,kind=4), 0.0_4, 5.0_4, &
                real(ymin,kind=4), real(ymax,kind=4), 0.0_4, 1.0_4)
      call title() ! Actually write the title to the file
      call setrgb(0.5_4, 0.5_4, 0.5_4)
      !call curve(real(x,kind=4),real(y,kind=4),nWriteTimes)  ! This draws the line
      call shdpat(16)  ! set shading pattern 
      call shdcrv(x,y,nWriteTimes,x,0.0_4*y,nWriteTimes) ! This fills below curve
      call color('FORE') ! Reset color to defaul foreground color

      !  Dislin Level 0:  before initialization or after termination
      call disfin()

      end subroutine write_DepPOI_TS_PNG_dislin
