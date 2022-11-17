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

      use time_data,     only : &
         ntmax,time_native

      use io_data,       only : &
         Site_vprofile,x_vprofile,y_vprofile,cdf_b3l1,VolcanoName

      use Source,        only : &
         neruptions,e_Volume,e_Duration,e_StartTime,e_PlumeHeight

      use time_data,     only : &
         os_time_log,BaseYear,useLeap

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
      write(dp_gnufile,53) vprof_ID,".gnu"
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
      write(dp_gnufile,53) plot_index,".gnu"
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
