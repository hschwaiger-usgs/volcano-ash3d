!##############################################################################
!
!    write_DepPOI_TS_PNG_plplot
!
!##############################################################################

      subroutine write_DepPOI_TS_PNG_plplot(pt_indx)

      use precis_param

      use Airports,      only : &
         nairports,Airport_Code,Airport_Name,Airport_x,Airport_y,&
         Airport_Latitude,Airport_Longitude,Airport_Thickness_TS

      use io_data,       only : &
         nWriteTimes,WriteTimes,VolcanoName

      use time_data,     only : &
         Simtime_in_hours

      use plplot

      implicit none

      integer :: pt_indx,i

      real(kind=8) :: ymaxpl
      character(len=14) :: dp_pngfile
      integer,save      :: plot_index = 0

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
        !write(*,*)"No deposit at ",pt_indx,Airport_Name(pt_indx)
        return
      else
        plot_index = plot_index + 1
        !write(*,*)"Processing ",pt_indx,plot_index,Airport_Name(pt_indx),'--'
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

      call plsetopt("geometry","400x300, 400x300")  ! Set image size
      call plsetopt("bg","FFFFFF")                  ! Set background color to white
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

      end subroutine write_DepPOI_TS_PNG_plplot
