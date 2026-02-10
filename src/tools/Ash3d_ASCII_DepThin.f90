! Ash3d_ASCII_DepThin calculates the deposit thickness as a function of vent distance
! In calculating the Thick v Dist function, the distance from each cell-center
! to the vent is calculated via the Haversine formula. The distance axis is binned
! with bin-width of 1.4*dx. If a cell has a greater thickness than previously tested,
! the T(d) value is replace so the final table represents the maximum thickness for
! a given distance bin. The mas thickness vs. distance data are written out to the file
! Tvd.dat and a gnuplot script written to Tvd.gnu. If gnuplot is found on the system,
! a plot is automatically created to display the results (DepoThick_vs_distance.png).


      program Ash3d_ASCII_DepThin

      use precis_param

      use io_units

      use Ash3d_ASCII_IO,  only : &
        A_XY,A_nx,A_ny,A_xll,A_yll,A_dx,A_dy,A_Fill, &
          read_2D_ASCII

      implicit none
      !implicit none (type, external)

      real(kind=dp), parameter :: DEG2RAD   = 1.74532925e-2_dp
      real(kind=dp), parameter :: Re        = 6371.229_dp

      integer           :: nargs
      integer           :: stat
      character (len= 80) :: linebuffer080

      character (len= 80) :: file1
      logical :: IsThere1

      integer :: nx_1
      integer :: ny_1
      real(kind=dp) :: xll_1
      real(kind=dp) :: yll_1
      real(kind=dp) :: dx_1
      real(kind=dp) :: dy_1
      real(kind=dp),dimension(:),allocatable :: lon_1
      real(kind=dp),dimension(:),allocatable :: lat_1
      real(kind=dp),dimension(:,:),allocatable :: XY_1
      integer :: i,j,ii

      real(kind=dp) :: srcx,srcy
      real(kind=dp) :: lon1,lat1     ! Input coords in rad for point 1
      real(kind=dp) :: dlon,dlat
      real(kind=dp) :: Rng,Rngmax
      real(kind=dp) :: a,c
      real(kind=dp) :: dx,dy,dxbin
      real(kind=dp),dimension(:),allocatable :: Rngvec
      real(kind=dp),dimension(:),allocatable :: Thickvec

      character(len=25) :: gnucom

      nio = 1  ! Turn off logging by setting output streams to stdout/stderr only

      nargs = command_argument_count()
      if (nargs == 0) then
          ! If no command-line arguments are given, then prompt user
          ! interactively for the ASCII file name and source coordinate
        do io=1,2;if(VB(io) <= verbosity_production)then
          write(outlog(io),*)'Ash3d_ASCII_DepThin calculates the deposit thickness as a function'
          write(outlog(io),*)'of distance from the vent. The expected usage is non-interactive via'
          write(outlog(io),*)'command-line arguments.'
          write(errlog(io),*)'  Usage: Ash3d_ASCII_DepThin file1 srcx srcy'
          write(outlog(io),*)'where file1 is a ESRI ASCII deposit file and srcx, srcy are the'
          write(outlog(io),*)'longitude and latitude of the vent. If gnuplot is available on the system'
          write(outlog(io),*)'the data are plotted to file DepoThick_vs_distance.png'
          write(outlog(io),*)' '
          write(outlog(io),*)'No command-line arguments were provided, now entering interactive mode.'
          write(outlog(io),*)' '
        endif;enddo
        write(output_unit,*)'Enter name of the ESRI ASCII deposit file:'
        read(input_unit,*) file1
        do io=1,nio;if(VB(io) <= verbosity_production)then
          write(outlog(io),*)'Enter lon (or x) of vent:'
        endif;enddo
        read(input_unit,*) srcx
        do io=1,nio;if(VB(io) <= verbosity_production)then
          write(outlog(io),*)'Enter lat (or y) of vent:'
        endif;enddo
        read(input_unit,*) srcy
      elseif (nargs == 1.or.nargs > 3) then
        do io=1,nio;if(VB(io) <= verbosity_error)then
          write(errlog(io),*)'ERROR: Too few command-line arguments.'
          write(errlog(io),*)'  Usage: Ash3d_ASCII_DepThin file1 srcx srcy'
        endif;enddo
        stop 1
      else
        call get_command_argument(1, linebuffer080, status=stat)
        if(stat > 0)then
          do io=1,nio;if(VB(io) <= verbosity_error)then
            write(errlog(io),*)'ERROR: Could not parse argument 1'
          endif;enddo
          stop 1
        elseif (stat < 0)then
          do io=1,nio;if(VB(io) <= verbosity_error)then
            write(errlog(io),*)'ERROR: Argument 1 has been truncated.'
            write(errlog(io),*)'       File name length is limited to 80 char.'
          endif;enddo
          stop 1
        endif
        file1=trim(adjustl(linebuffer080))
        inquire( file=adjustl(trim(file1)), exist=IsThere1 )
        if (.not.IsThere1)then
          do io=1,nio;if(VB(io) <= verbosity_error)then
            write(errlog(io),*)'ERROR: Input file 1 could not be found'
          endif;enddo
          stop 1
        endif

        call get_command_argument(2, linebuffer080, status=stat)
        if(stat > 0)then
          do io=1,nio;if(VB(io) <= verbosity_error)then
            write(errlog(io),*)'ERROR: Could not parse argument 2'
          endif;enddo
          stop 1
        endif
        read(linebuffer080,*)srcx

        call get_command_argument(3, linebuffer080, status=stat)
        if(stat > 0)then
          do io=1,nio;if(VB(io) <= verbosity_error)then
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
        lon_1(i) = xll_1 + 0.5_dp*dx_1 + (i-1)*dx_1
      enddo
      do j=1,ny_1
        lat_1(j) = yll_1 + 0.5_dp*dy_1 + (j-1)*dy_1
      enddo
      ! zero out any NaN so that it doesn't throw off the error
      do i=1,nx_1
        do j=1,ny_1
          if(abs(A_XY(i,j)-A_Fill) < 0.01_dp)A_XY(i,j)=0.0_dp
        enddo
      enddo
      XY_1  = A_XY
      deallocate(A_XY)

      ! get dx,dy in km
      dx = Re*dx_1*DEG2RAD*cos(DEG2RAD*minval(abs(lat_1(:))))
      dy = Re*dy_1*DEG2RAD
      dxbin = ceiling(max(dx,dy)*1.4_dp)

      allocate(Rngvec(nx_1+ny_1))
      allocate(Thickvec(nx_1+ny_1))
      do i=1,nx_1+ny_1
        Rngvec(i) = i*dxbin
      enddo
      Thickvec(:) = 0.0_dp

      Rngmax=0.0_dp
      do i=1,nx_1
        do j=1,ny_1
          if (XY_1(i,j) > 1.0e-3_dp)then
            lon1 = DEG2RAD * lon_1(i)
            lat1 = DEG2RAD * lat_1(j)
            dlon = srcx - lon1
            dlat = srcy - lat1
            ! Use Haversine formula
            a = sin(0.5_dp*dlat)**2.0_dp + cos(lat1)*cos(srcy)*sin(0.5_dp*dlon)**2.0_dp
            c = 2.0_dp*atan2(sqrt(a),sqrt(1.0_dp-a))
            Rng = Re * c
            if(Rng > Rngmax)Rngmax=Rng
            ii=ceiling(Rng/dxbin)
            if(XY_1(i,j) > Thickvec(ii))then
              Thickvec(ii) = XY_1(i,j)
            endif
          endif
        enddo
      enddo

      open(fid_outdata,file="Tvd.dat",status='replace')
      do i=1,nx_1+ny_1
        if(Rngvec(i) > Rngmax)exit
        write(fid_outdata,*)Rngvec(i)-0.5_dp*dxbin,Thickvec(i)
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

