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

      real(kind=8) :: ymaxpl
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
      call setfil(adjustl(trim(dp_pngfile))) ! Set output filename
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
