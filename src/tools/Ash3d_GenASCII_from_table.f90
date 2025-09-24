! Ash3d_GenASCII_from_table reads a table of coordinates + values and reformats data
! to an ESRI ASCII format
! Expects as input: table.txt ncols nrows XLLcorner YLLcorner dx dy

      program Ash3d_GenASCII_from_table

      use precis_param

      use io_units

      use global_param,   only : &
        EPS_TINY

      use Ash3d_ASCII_IO,  only : &
          write_2D_ASCII_flt

      implicit none

      integer           :: nargs
      integer           :: stat
      character(len=80) :: linebuffer080

      character(len=80) :: file1,file2
      logical :: IsThere1, IsThere2

      integer :: nx
      integer :: ny
      logical :: IsCC
      logical :: IsLL
      real(kind=sp) :: xll
      real(kind=sp) :: yll
      real(kind=sp) :: dx
      real(kind=sp) :: dy
      real(kind=sp),dimension(:)  ,allocatable :: X_dim
      real(kind=sp),dimension(:)  ,allocatable :: Y_dim
      real(kind=sp),dimension(:,:),allocatable :: OutVar
      character(len=20) :: filename =       'coord_tab_values.dat'
      integer :: i,j
      real(kind=sp)    :: inlon,inlat,value1
      integer          :: iostatus
      character(len=50):: iomessage
      character(len=6) :: Fill_Value = '-9999.'

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
            write(outlog(io),*)'Ash3d_GenASCII_from_table creates an ESRI ASCII file from an input table of'
            write(outlog(io),*)'coordinates + values.'
            write(outlog(io),*)'command-line arguments.'
            write(outlog(io),*)'  Usage: Ash3d_GenASCII_from_table table.txt ncols nrows XLLcorner YLLcorner dx dy'
            write(outlog(io),*)' '
            write(outlog(io),*)'No command-line arguments were provided, now entering interactive mode.'
            write(outlog(io),*)' '
            write(outlog(io),*)'Enter name of the file with tabular data:'
          endif;enddo
        endif
        read(input_unit,*) file1
        do io=1,nio;if(VB(io).le.verbosity_production)then
          write(outlog(io),*)'Enter the number of columns'
        endif;enddo
        read(input_unit,*) nx
        do io=1,nio;if(VB(io).le.verbosity_production)then
          write(outlog(io),*)'Enter the number of rows'
        endif;enddo
        read(input_unit,*) ny
        do io=1,nio;if(VB(io).le.verbosity_production)then
          write(outlog(io),*)'Enter the longitude of the lower-left corner'
        endif;enddo
        read(input_unit,*) xLL
        do io=1,nio;if(VB(io).le.verbosity_production)then
          write(outlog(io),*)'Enter the latitude of the lower-left corner'
        endif;enddo
        read(input_unit,*) yLL
        do io=1,nio;if(VB(io).le.verbosity_production)then
          write(outlog(io),*)'Enter the dx of the cells.'
        endif;enddo
        read(input_unit,*) dx
        do io=1,nio;if(VB(io).le.verbosity_production)then
          write(outlog(io),*)'Enter the dy of the cells.'
        endif;enddo
        read(input_unit,*) dy

      elseif (nargs.ne.7) then
        do io=1,nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)'ERROR: Too few command-line arguments.'
          write(errlog(io),*)'  Usage: Ash3d_GenASCII_from_table table.txt ncols nrows XLLcorner YLLcorner dx dy'
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
        if(stat.ne.0)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Could not parse argument 2'
          endif;enddo
          stop 1
        endif
        read(linebuffer080,*)nx

        call get_command_argument(3, linebuffer080, status=stat)
        if(stat.ne.0)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Could not parse argument 3'
          endif;enddo
          stop 1
        endif
        read(linebuffer080,*)ny

        call get_command_argument(4, linebuffer080, status=stat)
        if(stat.ne.0)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Could not parse argument 4'
          endif;enddo
          stop 1
        endif
        read(linebuffer080,*)xLL

        call get_command_argument(5, linebuffer080, status=stat)
        if(stat.ne.0)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Could not parse argument 5'
          endif;enddo
          stop 1
        endif
        read(linebuffer080,*)yLL

        call get_command_argument(6, linebuffer080, status=stat)
        if(stat.ne.0)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Could not parse argument 6'
          endif;enddo
          stop 1
        endif
        read(linebuffer080,*)dx

        call get_command_argument(7, linebuffer080, status=stat)
        if(stat.ne.0)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Could not parse argument 7'
          endif;enddo
          stop 1
        endif
        read(linebuffer080,*)dy

        if (.not.IsThere1)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Neither input file of tablular data could be found'
          endif;enddo
          stop 1
        endif
      endif

      allocate(X_dim(nx))
      do i=1,nx
        X_dim(i) = xLL + (i-1)*dx
      enddo
      allocate(Y_dim(ny))
      do i=1,ny
        Y_dim(i) = yLL + (i-1)*dy
      enddo
      allocate(OutVar(nx,ny))
      OutVar(:,:) = -9999.0_sp

      open(unit=fid_ctrlfile,file=file1,action='read')
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      do while (iostatus.eq.0)
        read(linebuffer080,*)inlon,inlat,value1
        i=nint((inlon-xLL)/dx) + 1
        j=nint((inlat-yLL)/dy) + 1
        OutVar(i,j) = value1
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      enddo
      close(fid_ctrlfile)

      IsLL = .true.
      IsCC = .false.
      call write_2D_ASCII_flt(nx,ny,IsLL,xLL,yLL,IsCC,dx,dy,Fill_Value,OutVar,filename)


      end program Ash3d_GenASCII_from_table

