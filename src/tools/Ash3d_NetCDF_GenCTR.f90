! Ash3d_NetCDF_GenCTR

      program Ash3d_NetCDF_GenCTR

      use precis_param

      use io_units

      use io_data,       only : &
         infile,concenfile,VolcanoName, &
         cdf_b1l2,cdf_vardz

      use mesh,          only : &
         lonLL,latLL,gridwidth_e,gridwidth_n,xLL,yLL,gridwidth_x,gridwidth_y, &
         IsLatLon,de,dn,dx,dy,VarDzType,dz_const

      use Source,        only : &
         lon_volcano,lat_volcano,x_volcano,y_volcano,z_volcano,Suzuki_A,      &
         neruptions,SourceType

      use Diffusion,     only : &
         diffusivity_horz

#ifdef USENETCDF
      use Ash3d_Netcdf_IO,  only : &
           NC_Read_Output_Products
#endif

      use Ash3d_Program_Control, only : &
         Write_input_block_header,      &
         SetWrite_input_block_1

      implicit none

      integer           :: nargs
      integer           :: stat
      character(len=80) :: linebuffer080

      logical       :: IsThere
      logical       :: WriteBlock = .true.
      real(kind=ip) :: LLx
      real(kind=ip) :: LLy
      real(kind=ip) :: widthx
      real(kind=ip) :: widthy
      real(kind=ip) :: x_in
      real(kind=ip) :: y_in
      real(kind=ip) :: z_in
      real(kind=ip) :: dx_in
      real(kind=ip) :: dy_in

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
        inquire( file=adjustl(trim(concenfile)), exist=IsThere )
        if (.not.IsThere)then
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

      ! Now that we have loaded the data from the NetCDF file, start writing
      ! to a control file
      infile = "temp.inp"
      inquire( file=infile, exist=IsThere )
      if(IsThere)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(errlog(io),*)"WARNING: Input file already exists"
        endif;enddo
      endif
      open(unit=fid_ctrlfile,file=infile,status='new',action='write')

      if(IsLatLon)then
        LLx    = lonLL
        LLy    = latLL
        widthx = gridwidth_e
        widthy = gridwidth_n
        x_in   = lon_volcano
        y_in   = lat_volcano
        z_in   = z_volcano
        dx_in  = de
        dy_in  = dn
      else
        LLx    = xLL
        LLy    = yLL
        widthx = gridwidth_x
        widthy = gridwidth_y
        x_in   = x_volcano
        y_in   = y_volcano
        z_in   = z_volcano
        dx_in  = dx
        dy_in  = dy
      endif

      call SetWrite_input_block_1(WriteBlock                       ,&  ! indicates that write to stdout as well as set vars
                                   VolcanoName                     ,&  ! volcano name
                                   cdf_b1l2                        ,&  ! projection line
                                   LLx,LLy                         ,&  ! x, y of LL corner
                                   widthx,widthy                   ,&  ! x and y width of domain
                                   x_in,y_in,z_in                  ,&  ! vent x,y,z
                                   dx_in,dy_in                     ,&  ! dx, dy
                                   dz_const                        ,&  ! dz
                                  0.0_ip,4.0_ip                    ,&  ! Kx, Suz param
                                  9                                ,&  ! # of eruptions
                                  1,cdf_vardz                      ,&  ! dz_type + bonus line
                                  1)                                   ! source_type



      close(fid_ctrlfile)

      end program Ash3d_NetCDF_GenCTR

