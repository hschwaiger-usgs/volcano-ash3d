!##############################################################################
!
!    write_2Dmap_PNG_gnuplot
!
!    if timestep = -1, then use the last step in file
!##############################################################################

      subroutine write_2Dmap_PNG_gnuplot(iprod,itime,OutVar,writeContours)

      use precis_param

      use global_param,  only : &
         EPS_SMALL

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
         Con_CloudTime_N,Con_CloudTime_RGB,Con_CloudTime_Lev, &
         ContourDataX,ContourDataY,ContourDataNcurves,ContourDataNpoints,&
         Contour_MaxCurves,Contour_MaxPoints,ContourLev,Con_N

      use time_data,     only : &
         os_time_log,BaseYear,useLeap

      use io_data,       only : &
         nWriteTimes,WriteTimes,cdf_b3l1,VolcanoName

      use Source,        only : &
         neruptions,e_Volume,e_Duration,e_StartTime,e_PlumeHeight,lon_volcano,lat_volcano

      implicit none

      integer :: iprod
      integer :: itime
      real(kind=ip) :: OutVar(nxmax,nymax)
      logical :: writeContours

      integer :: i,j,k,ii
      integer :: nzlev
      real(kind=4), dimension(:)  ,allocatable :: zlev
      integer     , dimension(:,:),allocatable :: zrgb
      character(len=40) :: title_plot
      character(len=15) :: title_legend
      character(len=40) :: outfile_name
      character (len=9) :: cio
      character (len=4) :: outfile_ext = '.png'

      real(kind=8)  :: xmin
      real(kind=8)  :: xmax
      real(kind=8)  :: ymin
      real(kind=8)  :: ymax

      character(len=10) :: dp_gnufile
      character(len=10) :: dp_outfile
      character(len=10) :: dp_confile
      character(len=26) :: coord_str
      character(len=25) :: gnucom
      integer :: ioerr,ioerr2,iw,iwf,istat

      integer :: ncities
      real(kind=8),dimension(:),allocatable     :: lon_cities
      real(kind=8),dimension(:),allocatable     :: lat_cities
      character(len=26),dimension(:),allocatable :: name_cities
      logical           :: IsThere1,IsThere2
      character(len=80) :: linebuffer080
      character         :: testkey
      integer           :: ilev,ignulev
      integer           :: lev_i,substr_pos1,substr_pos2,substr_pos3
      real(kind=4)      :: lev_r4
      integer           :: icurve,ipt

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

      inquire(file="world_50m.txt",exist=IsThere1)
      if(.not.IsThere1)then
        inquire(file="/opt/USGS/Ash3d/share/post_proc/world_50m.txt",exist=IsThere2)
        if(.not.IsThere2)then
          write(*,*)"Could not find required file world_50m.txt"
          write(*,*)"This file is available at:"
          write(*,*)"  http://www.gnuplotting.org_data_world_50m.txt"
          write(*,*)"Please download this file to the current working directory or"
          write(*,*)"copy to /opt/USGS/Ash3d/share/post_proc/"
        endif
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

      allocate(ContourLev(nzlev))
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
        zlev = (/0.1_op, 0.3_op, 1.0_op, 3.0_op, &
                10.0_op, 30.0_op, 100.0_op, 300.0_op/)
      elseif(iprod.eq.16)then   ! profile plots
        write(*,*)"ERROR: No map PNG output option for vertical profile data."
        write(*,*)"       Should not be in write_2Dmap_PNG_dislin"
        stop 1
      else
        write(*,*)"ERROR: unexpected variable"
        stop 1
      endif
      Con_N = nzlev

      if(writeContours)then
        write(*,*)"Running Gnuplot to calculate contours lines"
        allocate(ContourDataNcurves(nzlev))
        allocate(ContourDataNpoints(nzlev,Contour_MaxCurves))
        allocate(ContourDataX(nzlev,Contour_MaxCurves,Contour_MaxPoints))
        allocate(ContourDataY(nzlev,Contour_MaxCurves,Contour_MaxPoints))
        ContourDataNcurves(:)   = 0
        ContourDataNpoints(:,:) = 0
        ContourDataX(:,:,:)     = 0.0_8
        ContourDataY(:,:,:)     = 0.0_8
      else
        write(*,*)"Running Gnuplot to generate contour plot"
      endif

      write(dp_outfile,53) "outvar.dat"
      write(dp_confile,53) "outvar.con"
      write(dp_gnufile,53) "outvar.gpi"
 53   format(a10)

      xmin = real(minval(lon_cc_pd(1:nxmax)),kind=8)
      ! Make sure xmin is in the range -180->180
      if (xmin.gt.180.0_8)then
        xmin = real(minval(lon_cc_pd(1:nxmax)),kind=8)-360.0
        xmax = real(maxval(lon_cc_pd(1:nxmax)),kind=8)-360.0
      else
        xmax = real(maxval(lon_cc_pd(1:nxmax)),kind=8)
      endif
      ymin = real(minval(lat_cc_pd(1:nymax)),kind=8)
      ymax = real(maxval(lat_cc_pd(1:nymax)),kind=8)

      call citylist(2,xmin,xmax,ymin,ymax, &
                    ncities,                        &
                    lon_cities,lat_cities,          &
                    name_cities)

      if(lon_volcano.gt.xmax)lon_volcano=lon_volcano-360.0 

      ! write out the data in a form that gnuplot can read
      open(54,file=dp_outfile,status='replace')
      do i = 1,nxmax
        do j = 1,nymax
          write(54,*)lon_cc_pd(i)-360.0_8,lat_cc_pd(j),OutVar(i,j)
        enddo
        write(54,*)" "
      enddo
      close(54)
      open(54,file="volc.dat",status='replace')
      write(54,*)real(lon_volcano,kind=4),real(lat_volcano,kind=4),'""'
      close(54)

      ! Set up to plot via gnuplot script
      open(55,file=dp_gnufile,status='replace')
      write(55,*)"set terminal pngcairo font 'sans,12' size 854,603"   ! Set the image size
      write(55,*)"set origin 0.05, .20"
      write(55,*)"set size 0.85, 0.8"              ! Set x and y scale for plot
      write(55,*)"set ylabel 'Latitude'"
      write(55,*)"set xlabel 'Longitude'"
      write(55,*)"set output '",trim(adjustl(outfile_name)),"'"
      write(55,*)"set title '",title_plot,"'"
      write(55,*)"XMIN = ",real(xmin,kind=4)
      write(55,*)"YMIN = ",real(ymin,kind=4)
      write(55,*)"XMAX = ",real(xmax,kind=4)
      write(55,*)"YMAX = ",real(ymax,kind=4)

      write(55,*)"set xrange [XMIN:XMAX]"
      write(55,*)"set yrange [YMIN:YMAX]"
      write(55,*)"set contour base"
      write(55,*)"set cntrparam bspline"
      write(55,*)"set cntrparam levels discrete \"
      do i=1,nzlev-1
        write(55,*)zlev(i),', \'
      enddo
      write(55,*)zlev(nzlev)
      write(55,*)"unset surface"
      ! Now write out the contours to a datafile
      write(55,*)"set table 'outvar.con'"
      write(55,*)"splot 'outvar.dat' using 1:2:3"
      write(55,*)"unset table"
      write(55,*)"set style line 2 lc rgb '#808080' lt 0 lw 1"
      write(55,*)"set grid front ls 2"
      write(55,*)"unset key"

      write(55,*)"XVAL = XMIN-(XMAX-XMIN)*0.1"
      write(55,*)"YVAL = YMIN-(YMAX-YMIN)*0.25"

      write(55,*)"set label 'Volcano: " ,VolcanoName,&
                  "' at XVAL, YVAL font 'sans,9'"
      write(55,*)"set label 'Run Date: ",os_time_log,&
                  "' at XVAL, YVAL font 'sans,9' offset character 0,-1"
      read(cdf_b3l1,*,iostat=ioerr) iw,iwf
      write(55,*)"set label 'Windfile: ",iwf,&
                  "' at XVAL, YVAL font 'sans,9' offset character 0,-2"

      write(55,*)"XVAL = XMIN+(XMAX-XMIN)*0.4"
      write(55,*)"set label 'Erup. Start Time: ",HS_xmltime(e_StartTime(1),BaseYear,useLeap),&
                  "' at XVAL, YVAL font 'sans,9'"
      write(55,*)"set label 'Erup. Plume Height: ",real(e_PlumeHeight(1),kind=4),&
                  " km' at XVAL, YVAL font 'sans,9' offset character 0,-1"
      write(55,*)"set label 'Erup. Duration: ",real(e_Duration(1),kind=4),&
                  " hours' at XVAL, YVAL font 'sans,9' offset character 0,-2"
      write(55,*)"set label 'Erup. Volume: ",real(e_Volume(1),kind=4),&
                  " km3 (DRE)' at XVAL, YVAL font 'sans,9' offset character 0,-3"

      write(55,*)" plot 'world_50m.txt' with filledcurves linetype rgb '#dddddd' , \"
      write(55,*)"   'outvar.con' using 1:2 with l lc rgb '#888888' , \"
      write(55,*)"   '' every 1000 with labels font ',6' , \"
      write(55,*)"   'cities.xy' using 1:2 , \"
      write(55,*)"   '' using 1:2:3 with labels font ',10' point pointtype 7 offset char 1,1, \"
      write(55,*)"   'volc.dat' using 1:2 , \"
      write(55,*)"   '' using 1:2:3 with labels point pointtype 22 pointsize 2 lt rgb 'red'"

      close(55)

      write(gnucom,'(a11,a14)')'gnuplot -p ',dp_gnufile
      call execute_command_line(gnucom,exitstat=istat)

      if(writeContours)then

        ! Read outvar.con
        write(*,*)"Now reading outvar.con and loading contour data."
        open(54,file=dp_confile,status='old')
        ! In the gnuplot contour file, all contours of a certain level have a header
        ! in the following format:
        !# Contour 0, label:      300
        ! Each curve for that level is separated by a blank line
        ! Gnuplot writes the contour data to file starting with the highest level
        ! so we will need to check with zlev(:) to make sure we populate ContourDat
        ! correctly.
        ilev = -1
        ignulev = -1
        read(54,'(a80)',iostat=ioerr)linebuffer080
        do while(ioerr.ge.0)
          ! Check if this is a header line
          read(linebuffer080,*,iostat=ioerr2)testkey
          if(ioerr2.lt.0)then
            ! if there was an error trying to read a character, then this is a blank
            ! line; check if we are in a contour block or still in the file header
            if(ignulev.eq.-1)then
              ! Still in file header
              ! Read the next line and cycle
              read(54,'(a80)',iostat=ioerr)linebuffer080
              cycle
            else
              ! Blank line in a contour block means we are starting another curve
              ! Increment the number of curves for this level
              ContourDataNcurves(ilev) = ContourDataNcurves(ilev) + 1
              ! This is an easier index to used
              icurve = ContourDataNcurves(ilev)
            endif
          elseif(testkey.eq.'#')then
            ! This is a header line
            ! There are two possibilities:
            !  (1) a line of the file header
            !  (2) the start of a contour block
            !   if (2), then it will have this format:# Contour 0, label:      300
            substr_pos1 = index(linebuffer080,'Contour')

            if(substr_pos1.eq.0)then
              ! This is a file header line
              read(54,'(a80)',iostat=ioerr)linebuffer080
              cycle
            else
              read(linebuffer080,4,iostat=ioerr2)ignulev
4             format(9x,i2)
              ! Now read the level. Look for the ':' to isolate the last bit
              substr_pos1 = index(linebuffer080,':')
              ! Also look for a '.' or an 'e' to see if the level value is written as a real or int
              substr_pos2 = index(linebuffer080,'.')
              substr_pos3 = index(linebuffer080(20:),'e')
              if(substr_pos2.gt.0.or.substr_pos3.gt.0)then
                ! level is written as real
                read(linebuffer080(substr_pos1+1:28),*)lev_r4
              else
                read(linebuffer080(substr_pos1+1:28),*)lev_i
                lev_r4 = real(lev_i,kind=4)
              endif
              ! The value for this new level ignulev is lev_r4, but we need to find which
              ! zlev this corresponds to
              do ii = 1,nzlev
                if(abs(lev_r4-zlev(ii)).lt.EPS_SMALL)then
                  ilev=ii
                endif
              enddo
              ! When we have a new level, initialize the curve index for this level to 1
              ContourDataNcurves(ilev) = 1
              icurve = ContourDataNcurves(ilev)
            endif
          else
            ! This is the data section
            ! Increment the number of points
            ContourDataNpoints(ilev,icurve) = ContourDataNpoints(ilev,icurve) + 1
            ipt = ContourDataNpoints(ilev,icurve) 
            read(linebuffer080,*)ContourDataX(ilev,icurve,ipt),ContourDataY(ilev,icurve,ipt)
          endif

          ! Try to read the next line
          read(54,'(a80)',iostat=ioerr)linebuffer080
        enddo
        close(54)

        ! Loop through all the levels and curves and trim any curves with zero length
        do i=1,nzlev
          icurve = Contour_MaxCurves + 1
          do ii=Contour_MaxCurves,1,-1
            if(ContourDataNpoints(i,ii).le.0)then
              ! log each curve number with no points
              icurve = ii
            endif
          enddo
          ContourDataNcurves(i) = max(0,icurve-1)
        enddo

      endif

      end subroutine write_2Dmap_PNG_gnuplot


!##############################################################################
!
!    write_2Dprof_PNG_gnuplot
!
!##############################################################################

      subroutine write_2Dprof_PNG_gnuplot(vprof_ID)

      use precis_param

      use mesh,          only : &
         nzmax,z_cc_pd

      use Output_Vars,   only : &
         pr_ash

      use io_data,       only : &
         Site_vprofile,x_vprofile,y_vprofile,cdf_b3l1,VolcanoName

      use Source,        only : &
         neruptions,e_Volume,e_Duration,e_StartTime,e_PlumeHeight

      use time_data,     only : &
         os_time_log,BaseYear,useLeap,ntmax,time_native

      implicit none

      integer, intent (in) :: vprof_ID

      character(len=14) :: dp_gnufile
      character(len=14) :: dp_outfile
      character(len=14) :: dp_pngfile
      character(len=26) :: coord_str
      character(len=25) :: gnucom
      integer :: k,i
      integer :: ioerr,iw,iwf

      INTERFACE
        character (len=20) function HS_xmltime(HoursSince,byear,useLeaps)
          real(kind=8)              :: HoursSince
          integer                   :: byear
          logical                   :: useLeaps
        end function HS_xmltime
      END INTERFACE

      write(dp_outfile,53) vprof_ID,".dat"
      write(dp_gnufile,53) vprof_ID,".gpi"
      write(dp_pngfile,54) vprof_ID,".png"
 53   format('vprof_',i4.4,a4)
 54   format('gnupl_',i4.4,a4)

      open(54,file=dp_outfile,status='replace')
      do i = 1,ntmax
        do k = 1,nzmax
          write(54,*)time_native(i),z_cc_pd(k),pr_ash(k,i,vprof_ID)
        enddo
        write(54,*)" "
      enddo
      close(54)

      write(coord_str,101)x_vprofile(vprof_ID),y_vprofile(vprof_ID)
 101  format(' (lon=',f7.2,', lat=',f6.2,')')
      ! Set up to plot via gnuplot script
      open(55,file=dp_gnufile,status='replace')
      write(55,*)"set terminal pngcairo font 'sans,12' size 854,603"   ! Set the image size
      write(55,*)"set origin 0, .10"
      write(55,*)"set size 0.85, 0.9"              ! Set x and y scale for plot
      write(55,*)"set ylabel 'Height (km)'"
      write(55,*)"set xlabel 'Time (hours after eruption)'"
      write(55,*)"set output '",dp_pngfile,"'"
      write(55,*)"set title '",&
                  trim(adjustl(Site_vprofile(vprof_ID))),&
                  coord_str,"'"
      write(55,*)"set isosamples 50"
      write(55,*)"set pm3d"
      write(55,*)"set palette cubehelix negative"
      write(55,*)"unset surface"
      write(55,*)"set view map"
      write(55,*)"set key off"

      write(55,*)"XMIN = 0.0"
      write(55,*)"YMIN = 0.0"
      write(55,*)"XMAX = ",time_native(ntmax)
      write(55,*)"YMAX = ",z_cc_pd(nzmax)
      write(55,*)"XVAL = -XMAX*0.1"
      write(55,*)"YVAL = -YMAX*0.25"
      
      write(55,*)"set label 'Volcano: " ,VolcanoName,&
                  "' at XVAL, YVAL font 'sans,9'"
      write(55,*)"set label 'Run Date: ",os_time_log,&
                  "' at XVAL, YVAL font 'sans,9' offset character 0,-1"
      read(cdf_b3l1,*,iostat=ioerr) iw,iwf
      write(55,*)"set label 'Windfile: ",iwf,&
                  "' at XVAL, YVAL font 'sans,9' offset character 0,-2"

      write(55,*)"XVAL = XMAX*0.4"
      write(55,*)"set label 'Erup. Start Time: ",HS_xmltime(e_StartTime(1),BaseYear,useLeap),&
                  "' at XVAL, YVAL font 'sans,9'"
      write(55,*)"set label 'Erup. Plume Height: ",real(e_PlumeHeight(1),kind=4),&
                  " km' at XVAL, YVAL font 'sans,9' offset character 0,-1"
      write(55,*)"set label 'Erup. Duration: ",real(e_Duration(1),kind=4),&
                  " hours' at XVAL, YVAL font 'sans,9' offset character 0,-2"
      write(55,*)"set label 'Erup. Volume: ",real(e_Volume(1),kind=4),&
                  " km3 (DRE)' at XVAL, YVAL font 'sans,9' offset character 0,-3"

      write(55,*)"set cblabel 'Ash con. in mg/m3'"
      write(55,*)"splot '",dp_outfile,"'"

      close(55)

      write(gnucom,'(a11,a14)')'gnuplot -p ',dp_gnufile
      call execute_command_line(gnucom)

      end subroutine write_2Dprof_PNG_gnuplot

!##############################################################################

!##############################################################################
!
!    write_DepPOI_TS_PNG_gnuplot
!
!##############################################################################

      subroutine write_DepPOI_TS_PNG_gnuplot(pt_indx)

      use precis_param

      use Output_Vars,   only : &
         THICKNESS_THRESH

      use Airports,      only : &
         nairports,Airport_Code,Airport_Name,Airport_x,Airport_y,&
         Airport_Latitude,Airport_Longitude,Airport_Thickness_TS

      use io_data,       only : &
         nWriteTimes,WriteTimes,VolcanoName

      use time_data,     only : &
         Simtime_in_hours

      implicit none

      integer :: pt_indx,i

      real(kind=8) :: ymaxpl
      character(len=14) :: dp_gnufile
      character(len=14) :: dp_outfile
      character(len=14) :: dp_pngfile
      character(len=25) :: gnucom
      integer,save      :: plot_index = 0

      if(Airport_Thickness_TS(pt_indx,nWriteTimes).lt.THICKNESS_THRESH)then
        return
      else
        plot_index = plot_index + 1
      endif

      write(dp_outfile,53) plot_index,".dat"
      write(dp_gnufile,53) plot_index,".gpi"
      write(dp_pngfile,54) plot_index,".png"
 53   format('depTS_',i4.4,a4)
 54   format('gnupl_',i4.4,a4)

      open(54,file=dp_outfile,status='replace')
      do i = 1,nWriteTimes
        write(54,*)WriteTimes(i),Airport_Thickness_TS(pt_indx,i)
      enddo
      close(54)

      if(Airport_Thickness_TS(plot_index,nWriteTimes).lt.THICKNESS_THRESH)then
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
      write(55,*)"set title '",Airport_Name(pt_indx),"'"
      write(55,*)"plot [0:",ceiling(Simtime_in_hours),"][0:",&
                 nint(ymaxpl),"] '",dp_outfile,"' with filledcurve x1 ls 1"
      close(55)

      write(gnucom,'(a11,a14)')'gnuplot -p ',dp_gnufile
      call execute_command_line(gnucom)

      end subroutine write_DepPOI_TS_PNG_gnuplot
