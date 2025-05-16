! Ash3d_NetCDF_GenCTR

      program Ash3d_NetCDF_GenCTR

      use precis_param

      use io_units

      use io_data,       only : &
         concenfile

#ifdef USENETCDF
      use Ash3d_Netcdf_IO,  only : &
           NC_Read_Output_Products
#endif

      implicit none

      integer           :: nargs
      integer           :: stat
      character(len=80) :: linebuffer080

      logical :: IsThere1

      write(*,*)"In Ash3d_NetCDF_GenCTR"

      nargs = command_argument_count()
      if (nargs.eq.0) then
          ! If no command-line arguments are given, then prompt user
          ! interactively for the input netcdf file
        if(VB(1).ge.verbosity_silent)then
          write(errlog(1),*)"Stdout is suppressed via VERB=9, but interactive input is expected."
          write(errlog(1),*)"Either recompile with VERB<9 or provide the correct command-line arguments."
          stop 1
        else
          do io=1,nio;if(VB(io).le.verbosity_production)then
            write(outlog(io),*)'Enter name of the Ash3d netcdf output file'
          endif;enddo
        endif
        read(input_unit,*) concenfile
      elseif (nargs.gt.1) then
        do io=1,nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)'ERROR: Too many command-line arguments.'
          write(errlog(io),*)'  Usage: Ash3d_NetCDF_GenCTR output_file'
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
        concenfile=trim(adjustl(linebuffer080))
        inquire( file=adjustl(trim(concenfile)), exist=IsThere1 )
        if (.not.IsThere1)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Input file 1 could not be found'
          endif;enddo
          stop 1
        endif
      endif

#ifdef USENETCDF
        !call NC_Read_Output_Products(-1)
        call NC_Read_Output_Products(1)
#endif  

      end program Ash3d_NetCDF_GenCTR

