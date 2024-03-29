! Ash3d_ASCII_check calculates the difference between two ascii grid files

      program Ash3d_ASCII_check

      use precis_param

      use io_units

      use global_param,   only : &
        EPS_TINY

      use Ash3d_ASCII_IO,  only : &
        A_XY,A_nx,A_ny,A_xll,A_yll,A_dx,A_dy, &
          read_2D_ASCII

      implicit none

      integer           :: nargs
      integer           :: stat
      character(len=80) :: linebuffer080

      character(len=80) :: file1,file2
      logical :: IsThere1, IsThere2

      integer :: nx_1,nx_2
      integer :: ny_1,ny_2
      real(kind=ip) :: xll_1,xll_2
      real(kind=ip) :: yll_1,yll_2
      real(kind=ip) :: dx_1,dx_2
      real(kind=ip) :: dy_1,dy_2
      real(kind=ip),dimension(:,:),allocatable :: XY_1, XY_2
      integer :: i,j
      real(kind=ip) :: err

      real(kind=ip) :: tmp_ip
      real(kind=ip) :: L2_toterror
      real(kind=ip) :: L2_tol = 1.0e-3

      nio = 1  ! Turn off logging by setting output streams to stdout/stderr only

      nargs = command_argument_count()
      if (nargs.eq.0) then
          ! If no command-line arguments are given, then prompt user
          ! interactively for the two file names and the L2 tolerance
        if(VB(1).ge.verbosity_silent)then
          write(errlog(1),*)"Stdout is suppressed via VERB=9, but interactive input is expected."
          write(errlog(1),*)"Either recompile with VERB<9 or provide the correct command-line arguments."
          stop 1
        else
          do io=1,nio;if(VB(io).le.verbosity_production)then
            write(outlog(io),*)'Enter name of the first ESRI ASCII file:'
          endif;enddo
        endif
        read(input_unit,*) file1
        do io=1,nio;if(VB(io).le.verbosity_production)then
          write(outlog(io),*)'Enter name of the second ESRI ASCII file:'
        endif;enddo
        read(input_unit,*) file2
        do io=1,nio;if(VB(io).le.verbosity_production)then
          write(outlog(io),*)'Enter the L2 tolerance (default is 1.0e-7):'
        endif;enddo
        read(input_unit,*) L2_tol
      elseif (nargs.eq.1.or.nargs.gt.3) then
        do io=1,nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)'ERROR: Too few command-line arguments.'
          write(errlog(io),*)'  Usage: Ash3d_ASCII_check file1 file2 (tol.)'
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

        call get_command_argument(2, linebuffer080, status=stat)
        if(stat.gt.0)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Could not parse argument 2'
          endif;enddo
          stop 1
        elseif (stat.lt.0)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Argument 2 has been truncated.'
            write(errlog(io),*)'       File name length is limited to 80 char.'
          endif;enddo
          stop 1
        endif
        file2=trim(adjustl(linebuffer080))
        inquire( file=adjustl(trim(file2)), exist=IsThere2 )

        if (.not.IsThere1.and..not.IsThere2)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Neither input files could be found'
          endif;enddo
          stop 1
        elseif (.not.IsThere1)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Input file 1 could not be found'
          endif;enddo
          stop 1
        elseif (.not.IsThere2)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Input file 2 could not be found'
          endif;enddo
          stop 1
        endif

        if (nargs.eq.3)then
          call get_command_argument(3, linebuffer080, status=stat)
          read(linebuffer080,*)tmp_ip
          if(stat.eq.0)then
            L2_tol = tmp_ip
          else
            do io=1,nio;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)'WARNING: argument 3 not read correctly.'
              write(errlog(io),*)'         Using default value.'
            endif;enddo
            stop 1
          endif
        endif
      endif

      !do io=1,nio;if(VB(io).le.verbosity_error)then
      !  write(outlog(io),*)"Reading ASCII file 1"
      !endif;enddo
      call read_2D_ASCII(file1)
      nx_1  = A_nx
      ny_1  = A_ny
      xll_1 = A_xll
      yll_1 = A_yll
      dx_1  = A_dx
      dy_1  = A_dy
      allocate(XY_1(nx_1,ny_1))
      XY_1  = A_XY
      deallocate(A_XY)

      !do io=1,nio;if(VB(io).le.verbosity_error)then
      !  write(outlog(io),*)"Reading ASCII file 2"
      !endif;enddo
      call read_2D_ASCII(file2)
      nx_2  = A_nx
      ny_2  = A_ny
      xll_2 = A_xll
      yll_2 = A_yll
      dx_2  = A_dx
      dy_2  = A_dy
      allocate(XY_2(nx_2,ny_2))
      XY_2  = A_XY
      deallocate(A_XY)

      if(nx_1.ne.nx_2)then
        do io=1,nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"FAIL : nx differs"
        endif;enddo
        stop 1
      endif
      if(ny_1.ne.ny_2)then
        do io=1,nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"FAIL : ny differs"
        endif;enddo
        stop 1
      endif
      if(abs(dx_1-dx_2).gt.EPS_TINY)then
        do io=1,nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"FAIL : dx differs"
        endif;enddo
        stop 1
      endif
      if(abs(dy_1-dy_2).gt.EPS_TINY)then
        do io=1,nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"FAIL : dy differs"
        endif;enddo
        stop 1
      endif
      if(abs(xll_1-xll_2).gt.EPS_TINY)then
        do io=1,nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"FAIL : xll differs"
        endif;enddo
        stop 1
      endif
      if(abs(yll_1-yll_2).gt.EPS_TINY)then
        do io=1,nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"FAIL : yll differs"
        endif;enddo
        stop 1
      endif

      L2_toterror = 0.0_ip
      do i=1,nx_1
        do j=1,ny_1
          err=abs(XY_1(i,j)-XY_2(i,j))
          L2_toterror = L2_toterror + err*err
        enddo
      enddo
      L2_toterror = sqrt(L2_toterror)
      L2_toterror = L2_toterror/nx_1/ny_1

      do io=1,nio;if(VB(io).le.verbosity_error)then
        if(abs(L2_toterror).lt.L2_tol)then
          write(outlog(io),*)'PASS : ',L2_toterror
        else
          write(outlog(io),*)'FAIL : ',L2_toterror
        endif
      endif;enddo

      end program Ash3d_ASCII_check

