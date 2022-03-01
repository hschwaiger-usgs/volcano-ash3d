! Ash3d_ASCII_check calculates the difference between two ascii grid files

      program Ash3d_ASCII_check

      use precis_param

      use io_units

      use global_param,   only : &
        EPS_SMALL,EPS_TINY

      use Output_Vars

      implicit none

      integer             :: nargs
      integer             :: stat
      character(len=80):: linebuffer50

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

      INTERFACE
        subroutine read_2D_ASCII(filename)
          character(len=50),intent(in) :: filename
        end subroutine
      END INTERFACE

      nargs = command_argument_count()
      if (nargs.eq.0) then
          ! If no command-line arguments are given, then prompt user
          ! interactively for the two file names and the L2 tolerance
        write(global_info,*)'Enter name of the first ESRI ASCII file:'
        read(5,*) file1
        write(global_info,*)'Enter name of the second ESRI ASCII file:'
        read(5,*) file2
        write(global_info,*)'Enter the L2 tolerance (default is 1.0e-7):'
        read(5,*) L2_tol
      elseif (nargs.eq.1.or.nargs.gt.3) then
        write(global_error,*)'ERROR: Too few command-line arguments.'
        write(global_error,*)'  Usage: Ash3d_ASCII_check file1 file2 (tol.)'
        stop 1
      else
        call get_command_argument(1, linebuffer50, status=stat)
        if(stat.gt.0)then
          write(global_error,*)'ERROR: Could not parse argument 1'
          stop 1
        elseif (stat.lt.0)then
          write(global_error,*)'ERROR: Argument 1 has been truncated.'
          write(global_error,*)'       File name length is limited to 80 char.'
          stop 1
        endif
        file1=trim(adjustl(linebuffer50))
        inquire( file=adjustl(trim(file1)), exist=IsThere1 )

        call get_command_argument(2, linebuffer50, status=stat)
        if(stat.gt.0)then
          write(global_error,*)'ERROR: Could not parse argument 2'
          stop 1
        elseif (stat.lt.0)then
          write(global_error,*)'ERROR: Argument 2 has been truncated.'
          write(global_error,*)'       File name length is limited to 80 char.'
          stop 1
        endif
        file2=trim(adjustl(linebuffer50))
        inquire( file=adjustl(trim(file1)), exist=IsThere2 )

        if (.not.IsThere1.and..not.IsThere2)then
          write(global_error,*)'ERROR: Neither input files could be found'
          stop 1
        elseif (.not.IsThere1)then
          write(global_error,*)'ERROR: Input file 1 could not be found'
          stop 1
        elseif (.not.IsThere2)then
          write(global_error,*)'ERROR: Input file 2 could not be found'
          stop 1
        endif

        if (nargs.eq.3)then
          call get_command_argument(3, linebuffer50, status=stat)
          read(linebuffer50,*)tmp_ip
          if(stat.eq.0)then
            L2_tol = tmp_ip
          else
            write(global_error,*)'WARNING: argument 3 not read correctly.'
            write(global_error,*)'         Using default value.'
            stop 1
          endif
        endif
      endif

      call read_2D_ASCII(file1)
      nx_1  = R_nx
      ny_1  = R_ny
      xll_1 = R_xll
      yll_1 = R_yll
      dx_1  = R_dx
      dy_1  = R_dy
      allocate(XY_1(nx_1,ny_1))
      XY_1  = R_XY
      deallocate(R_XY)

      call read_2D_ASCII(file2)
      nx_2  = R_nx
      ny_2  = R_ny
      xll_2 = R_xll
      yll_2 = R_yll
      dx_2  = R_dx
      dy_2  = R_dy
      allocate(XY_2(nx_2,ny_2))
      XY_2  = R_XY
      deallocate(R_XY)

      if(nx_1.ne.nx_2)then
        write(*,*)"FAIL : nx differs"
        stop
      endif
      if(ny_1.ne.ny_2)then
        write(*,*)"FAIL : ny differs"
        stop
      endif
      if(abs(dx_1-dx_2).gt.EPS_TINY)then
        write(*,*)"FAIL : dx differs"
        stop
      endif
      if(abs(dy_1-dy_2).gt.EPS_TINY)then
        write(*,*)"FAIL : dy differs"
        stop
      endif
      if(abs(xll_1-xll_2).gt.EPS_TINY)then
        write(*,*)"FAIL : xll differs"
        stop
      endif
      if(abs(yll_1-yll_2).gt.EPS_TINY)then
        write(*,*)"FAIL : yll differs"
        stop
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

      if(abs(L2_toterror).lt.L2_tol)then
        write(global_info,*)'PASS : ',L2_toterror
      else
        write(global_info,*)'FAIL : ',L2_toterror
      endif

      end program Ash3d_ASCII_check

