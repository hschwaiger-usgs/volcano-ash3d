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
