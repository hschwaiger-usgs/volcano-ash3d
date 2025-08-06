! Ash3d_ASCII_DepThin calculates the deposit thickness as a function of vent distance

      program Ash3d_ASCII_DepThin

      use precis_param

      use io_units

      use Ash3d_ASCII_IO,  only : &
        A_XY,A_nx,A_ny,A_xll,A_yll,A_dx,A_dy,A_Fill, &
          read_2D_ASCII

      implicit none

      real(kind=8), parameter :: DEG2RAD   = 1.74532925e-2
      real(kind=8), parameter :: Re        = 6371.229

      integer           :: nargs
      integer           :: stat
      character(len=80) :: linebuffer080

      character(len=80) :: file1
      logical :: IsThere1

      integer :: nx_1
      integer :: ny_1
      real(kind=8) :: xll_1
      real(kind=8) :: yll_1
      real(kind=8) :: dx_1
      real(kind=8) :: dy_1
      real(kind=8),dimension(:),allocatable :: lon_1
      real(kind=8),dimension(:),allocatable :: lat_1
      real(kind=8),dimension(:,:),allocatable :: XY_1
      integer :: i,j,ii

      real(kind=8) :: srcx,srcy
      real(kind=8) :: lon1,lat1     ! Input coords in rad for point 1
      real(kind=8) :: dlon,dlat
      real(kind=8) :: Rng,Rngmax
      real(kind=8) :: a,c
      real(kind=8) :: dx,dy,dxbin
      real(kind=8),dimension(:),allocatable :: Rngvec
      real(kind=8),dimension(:),allocatable :: Thickvec

!      character(len=14) :: filename_script
!      character(len=14) :: filename_outdata
!      integer           :: fid_outdata  = 54
!      integer           :: fid_script  = 55
      character(len=25) :: gnucom

      nio = 1  ! Turn off logging by setting output streams to stdout/stderr only

      nargs = command_argument_count()
      if (nargs.eq.0) then
          ! If no command-line arguments are given, then prompt user
          ! interactively for the ASCII file name and source coordinate
        write(output_unit,*)'Enter name of the ESRI ASCII deposit file:'
        read(input_unit,*) file1
        do io=1,nio;if(VB(io).le.verbosity_production)then
          write(outlog(io),*)'Enter lon (or x) of vent:'
        endif;enddo
        read(input_unit,*) srcx
        do io=1,nio;if(VB(io).le.verbosity_production)then
          write(outlog(io),*)'Enter lat (or y) of vent:'
        endif;enddo
        read(input_unit,*) srcy
      elseif (nargs.eq.1.or.nargs.gt.3) then
        do io=1,nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)'ERROR: Too few command-line arguments.'
          write(errlog(io),*)'  Usage: Ash3d_ASCII_DepThin file1 srcx srcy'
        endif;enddo
        stop 1
      else
        call get_command_argument(1, linebuffer080, status=stat)
        if(stat.gt.0)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Could not parse argument 1'
          endif;enddo
          stop 1
        elseif (stat.lt.0)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Argument 1 has been truncated.'
            write(errlog(io),*)'       File name length is limited to 80 char.'
          endif;enddo
          stop 1
        endif
        file1=trim(adjustl(linebuffer080))
        inquire( file=adjustl(trim(file1)), exist=IsThere1 )
        if (.not.IsThere1)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Input file 1 could not be found'
          endif;enddo
          stop 1
        endif

        call get_command_argument(2, linebuffer080, status=stat)
        if(stat.gt.0)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Could not parse argument 2'
          endif;enddo
          stop 1
        endif
        read(linebuffer080,*)srcx

        call get_command_argument(3, linebuffer080, status=stat)
        if(stat.gt.0)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Could not parse argument 3'
          endif;enddo
          stop 1
        endif
        read(linebuffer080,*)srcy
      endif

      srcx = srcx*DEG2RAD
      srcy = srcy*DEG2RAD

      call read_2D_ASCII(file1)
      nx_1  = A_nx
      ny_1  = A_ny
      xll_1 = A_xll
      yll_1 = A_yll
      dx_1  = A_dx
      dy_1  = A_dy
      allocate(lon_1(nx_1))
      allocate(lat_1(ny_1))
      allocate(XY_1(nx_1,ny_1))
      do i=1,nx_1
        lon_1(i) = xll_1 + 0.5_8*dx_1 + (i-1)*dx_1
      enddo
      do j=1,ny_1
        lat_1(j) = yll_1 + 0.5_8*dy_1 + (j-1)*dy_1
      enddo
      ! zero out any NaN so that it doesn't throw off the error
      do i=1,nx_1
        do j=1,ny_1
          if(abs(A_XY(i,j)-A_Fill).lt.0.01_8)A_XY(i,j)=0.0_8
        enddo
      enddo
      XY_1  = A_XY
      deallocate(A_XY)

      ! get dx,dy in km
      dx = Re*dx_1*DEG2RAD*cos(DEG2RAD*minval(abs(lat_1(:))))
      dy = Re*dy_1*DEG2RAD
      dxbin = ceiling(max(dx,dy)*1.4_8)

      allocate(Rngvec(nx_1+ny_1))
      allocate(Thickvec(nx_1+ny_1))
      do i=1,nx_1+ny_1
        Rngvec(i) = i*dxbin
      enddo
      Thickvec(:) = 0.0_8

      Rngmax=0.0_8
      do i=1,nx_1
        do j=1,ny_1
          if (XY_1(i,j).gt.1.0e-3_8)then
            lon1 = DEG2RAD * lon_1(i)
            lat1 = DEG2RAD * lat_1(j)
            dlon = srcx - lon1
            dlat = srcy - lat1
            a = sin(0.5_8*dlat)**2.0_8 + cos(lat1)*cos(srcy)*sin(0.5_8*dlon)**2.0_8
            c = 2.0_8*atan2(sqrt(a),sqrt(1.0_8-a))
            Rng = Re * c
            if(Rng.gt.Rngmax)Rngmax=Rng
            ii=ceiling(Rng/dxbin)
            if(XY_1(i,j).gt.Thickvec(ii))then
              Thickvec(ii) = XY_1(i,j)
            endif
          endif
        enddo
      enddo

      open(fid_outdata,file="Tvd.dat",status='replace')
      do i=1,nx_1+ny_1
        if(Rngvec(i).gt.Rngmax)exit
        write(fid_outdata,*)Rngvec(i)-0.5_8*dxbin,Thickvec(i)
      enddo
      close(fid_outdata)

      !filename_script = "Tvd.gnu"
      open(fid_script,file="Tvd.gnu",status='replace')
      write(fid_script,*)"set terminal pngcairo font 'sans,12' size 854,603"
      write(fid_script,*)"set origin 0, .10"
      write(fid_script,*)"set size 0.85, 0.9"
      write(fid_script,*)"set ylabel 'Deposit Thickness (mm)'"
      write(fid_script,*)"set xlabel 'Distance from vent (km)'"
      write(fid_script,*)"set output 'DepoThick_vs_distance.png'"
      write(fid_script,*)"set key off"
      write(fid_script,*)"plot 'Tvd.dat' using 1:2 with lines"
      close(fid_script)

      write(gnucom,'(a18)')'gnuplot -p Tvd.gnu'
      call execute_command_line(gnucom)


      end program Ash3d_ASCII_DepThin

