!##############################################################################
!
!  This file contains all the subroutines for provided user-requested help
!  or usage information.
!
!      subroutine help_general
!      subroutine help_make
!      subroutine help_run
!      subroutine help_input
!      subroutine help_postproc
!      subroutine help_inputfile(blockID)
!
!##############################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  help_general
!
!  Called from: Parse_Command_Line
!  Arguments:
!    none
!
!  This subroutine is called if general Ash3d interactive help is requested
!  via 'Ash3d -h'.  Further information on different help options are printed
!  to stdout along with a call to help_run, giving information on how to
!  run Ash3d on the command line.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine help_general

      use io_units

      implicit none
      !implicit none (type, external)

      INTERFACE
        subroutine help_run
          implicit none
          !implicit none (type, external)
        end subroutine help_run
      END INTERFACE

#ifdef WINDOWS
      integer           :: iostatus
      character(len=120):: iomessage
      character(len=1)  :: key
#endif

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine help_general"
      endif;enddo

      io = 1

      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)' Ash3d help                                                                     '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)' You can get more specific help by typing:                                      '
      write(outlog(io),1)'   Ash3d -h make     : Information on build options, if you want to rebuild     '
      write(outlog(io),1)'   Ash3d -h run      : Information on running Ash3d from the command line       '
      write(outlog(io),1)'   Ash3d -h input    : Information on the structure of the input file           '
      write(outlog(io),1)'   Ash3d -h postproc : Information on the post-processing output results        '
      write(outlog(io),1)'   Ash3d -h info     : Information on this executable                           '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)' Writing run help information:                                                  '
      write(outlog(io),1)'                                                                                '

      call help_run

      ! For windows systems, allow the console to remain on the screen before exiting.
#ifdef WINDOWS
      write(outlog(io),*)'Press any key to exit.'
      read(input_unit,'(a1)',iostat=iostatus,iomsg=iomessage) key
#endif
      stop 1

 1    format(a80)

      end subroutine help_general

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  help_make
!
!  Called from: Parse_Command_Line
!  Arguments:
!    none
!
!  This subroutine writes to stdout the options that can be set by the user in
!  the makefile.  Lastly, controll is returned to Parse_Command_Line so the
!  subroutine Set_OS_Env can be called, printing out the settings used in this
!  current executable.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine help_make

      use io_units

      implicit none
      !implicit none (type, external)

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine help_make"
      endif;enddo

      io = 1

      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'To build Ash3d, we currently use a user-edited makefile.  All the main          '
      write(outlog(io),1)'variables to edit are in the top block of the makefile up to the lines:         '
      write(outlog(io),1)'############################################################################### '
      write(outlog(io),1)'#####  END OF USER SPECIFIED FLAGS  ########################################### '
      write(outlog(io),1)'############################################################################### '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'The following are the variables available to edit:                              '
      write(outlog(io),1)'SYSTEM      = [gfortran], ifort, aocc, nchpc                                    '
      write(outlog(io),1)'              This controls which compiler to use and which libraries to link.  '
      write(outlog(io),1)'              If you use a different compiler than these options, you will have '
      write(outlog(io),1)'              to edit the variables in the SYSTEM blocks to suit your system.   '
      write(outlog(io),1)'RUN         = DEBUG, PROF, [OPT], or OMPOPT                                     '
      write(outlog(io),1)'              This variable sets which set of compiler flags will be used for   '
      write(outlog(io),1)'              the executable, either with debugging, profiling, optimization,   '
      write(outlog(io),1)'              or with OMP enabled.                                              '
      write(outlog(io),1)'OS          = [LINUX], MACOS, WINDOWS                                           '
      write(outlog(io),1)'              This variable is used to indicate how Ash3d should build file     '
      write(outlog(io),1)'              paths when reading or writing files that require full paths.      '
      write(outlog(io),1)'USGSROOT    = [/opt/USGS]                                                       '
      write(outlog(io),1)'              This is the location of the USGS libraries and include files      '
      write(outlog(io),1)'              needed from volcano-ash3d-hourssince, volcano-ash3d-projection,   '
      write(outlog(io),1)'              and volcano-ash3d-metreader                                       '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'ASH3DCCSRC  = [./]                                                              '
      write(outlog(io),1)'              Location of the Ash3d core-code source files. This is normally    '
      write(outlog(io),1)'              just the cwd unless you are building an optional module and       '
      write(outlog(io),1)'              need to link to the code Ash3d files.                             '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'INSTALLDIR  = [/opt/USGS/Ash3d]                                                 '
      write(outlog(io),1)'              Location of the installation directory                            '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'USENETCDF   = [T] or F                                                          '
      write(outlog(io),1)'              Whether or not to build netCDF capabilities. If you only have     '
      write(outlog(io),1)'              netCDF v.3 available, then you will need to edit ncFPPFLAG        '
      write(outlog(io),1)'              to include -DNC3 so as to disable v.4 specific features.          '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'USEGRIB     = [T] or F                                                          '
      write(outlog(io),1)'              Whether or not to build GRIB2 capabilities via eccodes.           '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'USEPOINTERS = T or [F]                                                          '
      write(outlog(io),1)'              The default is for arrays to be allocatable, but if you are       '
      write(outlog(io),1)'              building Ash3d to be called from C++ codes, such as forestclaw,   '
      write(outlog(io),1)'              then some arrays need to be defined as pointer. This does not     '
      write(outlog(io),1)'              work with older versions of gfortran.                             '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'USEEXTDATA  = [T] or F                                                          '
      write(outlog(io),1)'              Ash3d can be built with lists of Airports and global volcanoes    '
      write(outlog(io),1)'              as data variables set at compile-time.  Some low-memory systems   '
      write(outlog(io),1)'              might not be able to compile with these data variables (e.g. on   '
      write(outlog(io),1)'              Raspberry Pi systems).  Alternatively, these lists can be read    '
      write(outlog(io),1)'              at run-time from files installed on the system.                   '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'FASTFPPFLAG = [] -DFAST_DT -DFAST_SUBGRID                                       '
      write(outlog(io),1)'              These are additional pre-processor flags that can be included     '
      write(outlog(io),1)'              to speed up the computation.                                      '
      write(outlog(io),1)'              -DFAST_DT limits DT evaluations to just on the steps of the       '
      write(outlog(io),1)'               windfiles with linear interpolation of dt between steps. This    '
      write(outlog(io),1)'               cannot be used with sources that effect the winds (umbrella).    '
      write(outlog(io),1)'              -DFAST_SUBGRID adjusts the computational domain to be just the    '
      write(outlog(io),1)'               min/max in x,y,z of the region where ash concentration exceeds   '
      write(outlog(io),1)'               some threshhold.                                                 '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'USEZIP      = [T] or F                                                          '
      write(outlog(io),1)'              If point data are requested in kml, Ash3d used zip to bundle the  '
      write(outlog(io),1)'              time-series accumulation plots into a kmz file. If zip in not     '
      write(outlog(io),1)'              available on the system, set this to F.                           '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'USEDISLIN   = T or [F]                                                          '
      write(outlog(io),1)'              If Ash3d_PostProc is to be built, this variable indicates if      '
      write(outlog(io),1)'              the dislin library is available on the system for plotting maps   '
      write(outlog(io),1)'              and/or vertical profiles.                                         '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'USEPLPLOT   = T or [F]                                                          '
      write(outlog(io),1)'              If Ash3d_PostProc is to be built, this variable indicates if      '
      write(outlog(io),1)'              the plplot library is available on the system for plotting maps   '
      write(outlog(io),1)'              and/or vertical profiles.                                         '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'USEGNUPLOT  = [T] or F                                                          '
      write(outlog(io),1)'              If point data are requested in kml format, Ash3d uses gnuplot, or '
      write(outlog(io),1)'              one of the other plotting packages, to generate accumulation      '
      write(outlog(io),1)'              time-series plots. If gnuplot is not available on the system, set '
      write(outlog(io),1)'              this to F.                                                        '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'USEGMT      = [F]                                                               '
      write(outlog(io),1)'              If Ash3d_PostProc is to be built, this variable indicates if      '
      write(outlog(io),1)'              the GMT library is available on the system for plotting maps      '
      write(outlog(io),1)'              and/or vertical profiles. Currently, this is not enabled.         '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'MFILE       = [makefile]                                                        '
      write(outlog(io),1)'              The name of the makefile.                                         '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'LIMITER     = LIM_NONE, LIM_LAXWEN, LIM_BW, LIM_FROMM, LIM_MINMOD,              '
      write(outlog(io),1)'              [LIM_SUPERBEE], LIM_MC                                            '
      write(outlog(io),1)'              This variable specifies the limiter to use in the advection       '
      write(outlog(io),1)'              routines.                                                         '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'DIFFMETH    = EXPLDIFF, [CRANKNIC]                                              '
      write(outlog(io),1)'              This variable specifies the diffusion scheme to be used.          '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'PII         = ON, [OFF]                                                         '
      write(outlog(io),1)'              This variable allows information about a run to be included in    '
      write(outlog(io),1)'              the output netCDF file. This is set to OFF by default as the      '
      write(outlog(io),1)'              full path of the run directory, the hostname and the username     '
      write(outlog(io),1)'              might be considered sensitive information, but this information   '
      write(outlog(io),1)'              may be desired for logging of runs.                               '
      write(outlog(io),1)'                                                                                '

 1    format(a80)

      end subroutine help_make

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  help_run
!
!  Called from: Parse_Command_Line and help_general
!  Arguments:
!    none
!
!  This subroutine gives a bried description of what is required to run Ash3d
!  from the command line.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine help_run

      use io_units

      implicit none
      !implicit none (type, external)

#ifdef WINDOWS
      integer           :: iostatus
      character(len=120):: iomessage
      character(len=1)  :: key
#endif

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine help_run"
      endif;enddo

      io = 1

      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'Ash3d requires a control file to run.  Example control files are given in this  '
      write(outlog(io),1)'repository in the examples folder.  To run Ash3d with a control file, simply    '
      write(outlog(io),1)'pass the control file name as a command-line argument to the executable:        '
      write(outlog(io),1)'   ./Ash3d Spurr_081992_ESP.inp                                                 '
      write(outlog(io),1)'Alternatively, Ash3d can be run interactively, where the user will be prompted  '
      write(outlog(io),1)'for the control file, followed by a prompt for restarting a previous run, if    '
      write(outlog(io),1)'desired.  If a restart case is requested, the control file and the netcdf file  '
      write(outlog(io),1)'must be consistent.                                                             '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'There are several environment variables that can be used to modify the run      '
      write(outlog(io),1)'conditions of Ash3d, if desired.                                                '
      write(outlog(io),1)'  ASH3DVERB = 1-10                                                              '
      write(outlog(io),1)'   This overrides the default level of 3 for output.                            '
      write(outlog(io),1)'    1 = debug2     : Additional debugging information only written to stdout    '
      write(outlog(io),1)'    2 = debug1     : Debugging information only written to stdout               '
      write(outlog(io),1)'    3 = log        : Time step information (limit for writing to logfile)       '
      write(outlog(io),1)'    4 = info       : Additional information on run set up and shutdown          '
      write(outlog(io),1)'    5 = statistics : Details on health of run (timing, mass conservation)       '
      write(outlog(io),1)'    6 = production : Major program flow info only                               '
      write(outlog(io),1)'    7 = essential  : Only start up and shutdown messages                        '
      write(outlog(io),1)'    8 = error      : No logging to stdout, only stderr (and logfile)            '
      write(outlog(io),1)'    9 = silent     : No logging to stdout,stderr. Logfile written as normal     '
      write(outlog(io),1)'   10 = dark       : No logging to stdout,stderr or logfile                     '
      write(outlog(io),1)'To run Ash3d with no output to the terminal, but with a logfile, for example:   '
      write(outlog(io),1)'       ASH3DVERB=9 ./Ash3d control.inp                                          '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'  ASH3DHOME = path to installation location                                     '
      write(outlog(io),1)'   If the executable is compiled such that it needs access to external files in '
      write(outlog(io),1)'   the installation directory, the path can be changed at run-time, if needed.  '
      write(outlog(io),1)'       ASH3DHOME=/opt/USGS/Ash3d ./Ash3d control.inp                            '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'  ASH3DCFL = 0.0-1.0                                                            '
      write(outlog(io),1)'   The Courant-Friedrichs-Lewy (CFL) factor can be changed at run-time from the '
      write(outlog(io),1)'   default of 0.8.                                                              '
      write(outlog(io),1)'       ASH3DCFL=0.5 ./Ash3d control.inp                                         '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'If the executable was compiled with OpenMP, then the number of threads can be   '
      write(outlog(io),1)'specified at run-time.                                                          '
      write(outlog(io),1)'       OMP_NUM_THREADS=4 ./Ash3d_opt control.in                                 '
      write(outlog(io),1)'                                                                                '

      ! For windows systems, allow the console to remain on the screen before exiting.
#ifdef WINDOWS
      write(outlog(io),*)'Press any key to exit.'
      read(input_unit,'(a1)',iostat=iostatus,iomsg=iomessage) key
#endif
      stop 1

 1    format(a80)

      end subroutine help_run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  help_input
!
!  Called from: Parse_Command_Line
!  Arguments:
!    none
!
!  This subroutine first writes to stdout a general description of the structure
!  of the Ash3d control file.  Then the user is prompted to request detailed
!  information about each of the blocks of the control file via call to
!  help_inputfile(blockID).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine help_input

      use io_units

      implicit none
      !implicit none (type, external)

      INTERFACE
        subroutine help_inputfile(blockID)
          implicit none
          !implicit none (type, external)
          integer,intent(in) :: blockID
        end subroutine help_inputfile
      END INTERFACE

      character(len=3)  :: answer
      character(len=50) :: linebuffer050
      character(len=80) :: linebuffer080
      integer           :: blockID
      integer           :: iostatus
      character(len=120):: iomessage
      integer           :: i

#ifdef WINDOWS
      character(len=1)  :: key
#endif

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine help_input"
      endif;enddo

      io = 1

      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'For information of the structure of the blocks of the input file, please enter  '
      write(outlog(io),1)'the block number 1-10, or press return to exit.                                 '
      write(outlog(io),1)'The control file is structures with at least 8 blocks:                          '
      write(outlog(io),1)' BLOCK 1: GRID INFO                                                             '
      write(outlog(io),1)' BLOCK 2: ERUPTION PARAMETERS                                                   '
      write(outlog(io),1)' BLOCK 3: WIND PARAMETERS                                                       '
      write(outlog(io),1)' BLOCK 4: OUTPUT OPTIONS                                                        '
      write(outlog(io),1)' BLOCK 5: INPUT WIND FILES                                                      '
      write(outlog(io),1)' BLOCK 6: AIRPORT FILE                                                          '
      write(outlog(io),1)' BLOCK 7: GRAIN-SIZE BINS, SETTLING VELOCITY                                    '
      write(outlog(io),1)' BLOCK 8: VERTICAL PROFILES                                                     '
      write(outlog(io),1)' BLOCK 9 (Optional): NETCDF ANNOTATIONS                                         '
      write(outlog(io),1)' BLOCK 10+: OPTIONAL MODULES                                                    '
      write(outlog(io),1)'   RESETPARAMS                                                                  '
      write(outlog(io),1)'   TOPO                                                                         '
      write(outlog(io),1)'   VARDIFF                                                                      '
      write(outlog(io),1)'Would you like information on the structure of these blocks? (y or n)           '

      read(input_unit,'(a3)',iostat=iostatus,iomsg=iomessage) answer
      linebuffer080 = answer
      linebuffer050 = "Reading from stdin, answer"
      if(iostatus /= 0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      do while (adjustl(trim(answer)) == 'y'.or.adjustl(trim(answer)) == 'yes')
        write(outlog(io),*)'  Please enter the block ID (1-10) or 0 for all:'
        read(input_unit,*,err=10,iostat=iostatus,iomsg=iomessage)blockID
        linebuffer080 = "blockID"
        linebuffer050 = "Reading from stdin, blockID"
        if(iostatus /= 0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        if(blockID < 0.or.blockID > 10)then
          write(errlog(io),*)'  Invalid range for blockID; should be between 0 and 10.'
        elseif(blockID == 0)then
          do i=1,10
            call help_inputfile(i)
          enddo
          answer = 'no '
        else
          call help_inputfile(blockID)
          write(outlog(io),1)'                                                                                '
          write(outlog(io),1)'Would you like information on another block? (y or n)                           '
          read(input_unit,'(a3)',iostat=iostatus,iomsg=iomessage) answer
          linebuffer080 = answer
          linebuffer050 = "Reading from stdin, answer"
          if(iostatus /= 0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        endif
      enddo

      ! For windows systems, allow the console to remain on the screen before exiting.
#ifdef WINDOWS
      write(outlog(io),*)'Press any key to exit.'
      read(input_unit,'(a1)',iostat=iostatus,iomsg=iomessage) key
#endif
 10   stop 1

 1    format(a80)

      end subroutine help_input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  help_inputfile(blockID)
!
!  Called from: Parse_Command_Line and help_input
!  Arguments:
!    blockID = block number of the control file to print
!
!  This subroutine writes an example of the requested block of the control file
!  with descriptive comments above the block.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine help_inputfile(blockID)

      use precis_param

      use io_units

      use Ash3d_Program_Control, only : &
         Write_input_block_header,      &
         SetWrite_input_block_01,       &
         SetWrite_input_block_02,       &
         SetWrite_input_block_03,       &
         SetWrite_input_block_04,       &
         SetWrite_input_block_05,       &
         SetWrite_input_block_06,       &
         SetWrite_input_block_07,       &
         SetWrite_input_block_08,       &
         SetWrite_input_block_09,       &
         SetWrite_input_block_ResetParam, &
         SetWrite_input_block_Topo,     &
         SetWrite_input_block_VarDiff

      implicit none
      !implicit none (type, external)

      integer,intent(in) :: blockID

      !logical            :: WriteBlock = .true.
      character(len= 30) :: vname
      character(len=130) :: dz_line
      character(len= 80) :: projline

      dz_line = '                                                                               ' // &
                '                                                 '

      ! The idea with the blockID is that help for only a particular block is
      ! called if there is an error reading something in the input file or if
      ! the user requests it.

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine help_inputfile"
      endif;enddo

      select case (blockID)
        case(1)  ! BLOCK 1: GRID INFO
          call Write_input_block_header(output_unit,blockID)
          vname =  "Eyjafjallajokull             "
          projline = '1 1 0.0 90.0 0.933 6367.470                                                    '
          !call SetWrite_input_block_01(WriteBlock                       ,&  ! indicates that write to stdout as well as set vars
          !                             vname                           ,&  ! volcano name
          !                             projline                        ,&  ! projection line
          !                            -2500.0_ip,-4500.0_ip            ,&  ! x, y of LL corner
          !                             55.0_ip,25.0_ip                 ,&  ! x and y width of domain
          !                            -19.62_ip,63.63_ip,0.0_ip        ,&  ! vent x,y,z
          !                            3.0_ip,3.0_ip                    ,&  ! dx, dy
          !                            1.0_ip                           ,&  ! dz
          !                            0.0_ip,4.0_ip                    ,&  ! Kx, Suz param
          !                            9                                ,&  ! # of eruptions
          !                            1,dz_line                        ,&  ! dz_type + bonus line
          !                            1)                                   ! source_type
      write(output_unit,1)'******************* BLOCK 1 ***************************************************                        '
      write(output_unit,1)'Eyjafjallajokull               # Volcano name (character*30)                                           '
      write(output_unit,1)'1 1 0.0 90.0 0.933 6367.470    # Proj flags and params; first term (LLflag) is 1, so all else ignored  '
      write(output_unit,1)'-25.0    45.0                  # x, y of LL corner of grid (km, or deg. if latlongflag=1)              '
      write(output_unit,1)'55.0     25.0                  # grid width and height (km, or deg. if latlonflag=1)                   '
      write(output_unit,1)'-19.62   63.63                 # vent location         (km, or deg. if latlonflag=1)                   '
      write(output_unit,1)'3.0      3.0                   # DX, DY of grid cells  (km, or deg. if latlonflag=1)                   '
      write(output_unit,1)'1.0                            # DZ of grid cells      (always km)                                     '
      write(output_unit,1)'0.     4.                      # diffusion coefficient (m2/s), Suzuki constant                         '
      write(output_unit,1)'9                              # neruptions, number of eruptions or pulses                             '
      write(output_unit,1)'*******************************************************************************                        '
        case(2)  ! BLOCK 2: ERUPTION PARAMETERS
          call Write_input_block_header(output_unit,blockID)
          !call SetWrite_input_block_02(WriteBlock)
      write(output_unit,1)'******************* BLOCK 2 ***************************************************                        '
      write(output_unit,1)'2010 4 14  9.0  3.0 7.4 8.46E-004                                                                      '
      write(output_unit,1)'2010 4 14 12.0  3.0 8.4 1.38E-003                                                                      '
      write(output_unit,1)'2010 4 14 15.0  3.0 6.8 6.11E-004                                                                      '
      write(output_unit,1)'2010 4 14 18.0  3.0 5.6 2.89E-004                                                                      '
      write(output_unit,1)'2010 4 14 21.0  3.0 5.1 2.01E-004                                                                      '
      write(output_unit,1)'2010 4 15  0.0  3.0 5.3 2.33E-004                                                                      '
      write(output_unit,1)'2010 4 15  3.0  3.0 5.2 2.17E-004                                                                      '
      write(output_unit,1)'2010 4 15  6.0  3.0 5.3 2.33E-004                                                                      '
      write(output_unit,1)'2010 4 15  9.0  3.0 5.7 3.09E-004                                                                      '
        case(3)  ! BLOCK 3: WIND PARAMETERS
          call Write_input_block_header(output_unit,blockID)
          !call SetWrite_input_block_03(WriteBlock)
      write(output_unit,1)'******************* BLOCK 3 ***************************************************                        '
      write(output_unit,1)'4  20               # iwind, iwindformat, [igrid, idata]                                               '
      write(output_unit,1)'2                   # iHeightHandler                                                                   '
      write(output_unit,1)'60                  # Simulation time in hours                                                         '
      write(output_unit,1)'no                  # stop computation when 99% of erupted mass has deposited?                         '
      write(output_unit,1)'16                  # nWindFiles, number of gridded wind files (used if iwind>1)                       '
        case(4)  ! BLOCK 4: OUTPUT OPTIONS
          call Write_input_block_header(output_unit,blockID)
          !call SetWrite_input_block_04(WriteBlock)
      write(output_unit,1)'******************* BLOCK 4 ***************************************************                        '
      write(output_unit,1)'no      # Write out ESRI ASCII file of final deposit thickness?                                        '
      write(output_unit,1)'yes     # Write out        KML file of final deposit thickness?                                        '
      write(output_unit,1)'no      # Write out ESRI ASCII deposit files at specified times?                                       '
      write(output_unit,1)'no      # Write out        KML deposit files at specified times?                                       '
      write(output_unit,1)'no      # Write out ESRI ASCII files of ash-cloud concentration?                                       '
      write(output_unit,1)'no      # Write out        KML files of ash-cloud concentration?                                      '
      write(output_unit,1)'no      # Write out ESRI ASCII files of ash-cloud height?                                              '
      write(output_unit,1)'no      # Write out        KML files of ash-cloud height?                                              '
      write(output_unit,1)'yes     # Write out ESRI ASCII files of ash-cloud load (T/km2) at specified times?                     '
      write(output_unit,1)'yes     # Write out        KML files of ash-cloud load (T/km2) at specified times?                     '
      write(output_unit,1)'yes     # Write out ESRI ASCII file of deposit arrival times?                                          '
      write(output_unit,1)'yes     # Write out        KML file of deposit arrival times?                                          '
      write(output_unit,1)'yes     # write out ESRI ASCII file of cloud arrival times?                                            '
      write(output_unit,1)'yes     # Write out        KML file of cloud arrival times?                                            '
      write(output_unit,1)'yes 1   # Write out 3-D ash concentration at specified times? / [output code: 1=2d+concen,2=2d only]   '
      write(output_unit,1)'netcdf  # format of ash concentration files     (ascii, binary, or netcdf)                             '
      write(output_unit,1)'-1      # nWriteTimes                                                                                  '
      write(output_unit,1)'1       # WriteTimes (hours since eruption start)                                                      '
        case(5)  ! BLOCK 5: INPUT WIND FILES
          call Write_input_block_header(output_unit,blockID)
          !call SetWrite_input_block_05(WriteBlock)
      write(output_unit,1)'******************* BLOCK 5 ***************************************************                        '
      write(output_unit,1)'Wind_nc/gfs/gfs.2010041400/2010041400.f000.nc                                                          '
      write(output_unit,1)'Wind_nc/gfs/gfs.2010041400/2010041400.f003.nc                                                          '
      write(output_unit,1)'Wind_nc/gfs/gfs.2010041400/2010041400.f006.nc                                                          '
      write(output_unit,1)'Wind_nc/gfs/gfs.2010041400/2010041400.f009.nc                                                          '
      write(output_unit,1)'Wind_nc/gfs/gfs.2010041400/2010041400.f012.nc                                                          '
      write(output_unit,1)'Wind_nc/gfs/gfs.2010041400/2010041400.f015.nc                                                          '
      write(output_unit,1)'Wind_nc/gfs/gfs.2010041400/2010041400.f018.nc                                                          '
      write(output_unit,1)'Wind_nc/gfs/gfs.2010041400/2010041400.f021.nc                                                          '
      write(output_unit,1)'Wind_nc/gfs/gfs.2010041400/2010041400.f024.nc                                                          '
      write(output_unit,1)'Wind_nc/gfs/gfs.2010041400/2010041400.f027.nc                                                          '
      write(output_unit,1)'Wind_nc/gfs/gfs.2010041400/2010041400.f030.nc                                                          '
      write(output_unit,1)'Wind_nc/gfs/gfs.2010041400/2010041400.f033.nc                                                          '
      write(output_unit,1)'Wind_nc/gfs/gfs.2010041400/2010041400.f036.nc                                                          '
      write(output_unit,1)'Wind_nc/gfs/gfs.2010041400/2010041400.f039.nc                                                          '
      write(output_unit,1)'Wind_nc/gfs/gfs.2010041400/2010041400.f042.nc                                                          '
      write(output_unit,1)'Wind_nc/gfs/gfs.2010041400/2010041400.f045.nc                                                          '
        case(6)  ! BLOCK 6: AIRPORT FILE
          call Write_input_block_header(output_unit,blockID)
          !call SetWrite_input_block_06(WriteBlock)
      write(output_unit,1)'******************* BLOCK 6 ***************************************************                        '
      write(output_unit,1)'no                            # Write out ash arrival times at airports to ASCII FILE?                 '
      write(output_unit,1)'no                            # Write out grain-size distribution to ASCII airports file?              '
      write(output_unit,1)'no                            # Write out ash arrival times to kml file?                               '
      write(output_unit,1)'GlobalAirports.txt            # Name of file containing airport locations                              '
      write(output_unit,1)'yes                           # Defer to Lon/Lat coordinates? ("no" defers to projected)               '
        case(7)  ! BLOCK 7: GRAIN-SIZE BINS, SETTLING VELOCITY
          call Write_input_block_header(output_unit,blockID)
          !call SetWrite_input_block_07(WriteBlock)
      write(output_unit,1)'******************* BLOCK 7 ***************************************************                        '
      write(output_unit,1)'15                                 # Number of grain-size bins. FV_ID not given; defaults to 1         '
      write(output_unit,1)'2.000      0.0208  2003.   0.44    # grain size (mm), mass fraction, density (kg/m3), [F=0.44]         '
      write(output_unit,1)'1.414      0.0084  2350.   0.44                                                                        '
      write(output_unit,1)'1.000      0.0141  2005.   0.44                                                                        '
      write(output_unit,1)'0.707      0.0214  2248.   0.44                                                                        '
      write(output_unit,1)'0.500      0.0459  2624.   0.44                                                                        '
      write(output_unit,1)'0.354      0.0723  2644.   0.44                                                                        '
      write(output_unit,1)'0.250      0.0532  2639.   0.44                                                                        '
      write(output_unit,1)'0.177      0.0219  2690.   0.44                                                                        '
      write(output_unit,1)'0.125      0.0165  2691.   0.44                                                                        '
      write(output_unit,1)'0.088      0.0115  2730.   0.44                                                                        '
      write(output_unit,1)'0.2176     0.0714  600.    1.00    # Note that these are bigger again, low density and round           '
      write(output_unit,1)'0.2031     0.1428  600.    1.00    # These bottom five bins represent an aggregate distribution        '
      write(output_unit,1)'0.1895     0.2856  600.    1.00                                                                        '
      write(output_unit,1)'0.1768     0.1428  600.    1.00                                                                        '
      write(output_unit,1)'0.1649     0.0714  600.    1.00                                                                        '
        case(8)  ! BLOCK 8: VERTICAL PROFILES
          call Write_input_block_header(output_unit,blockID)
          !call SetWrite_input_block_08(WriteBlock)
      write(output_unit,1)'******************* BLOCK 8 ***************************************************                        '
      write(output_unit,1)'4                             # number of locations for vertical profiles (nlocs)                      '
      write(output_unit,1)'12.4  51.4  Leipzig           # x,y (or lon/lat) [Site name]                                           '
      write(output_unit,1)'11.3  48.2  Munich            # Munich (Maisach)                                                       '
      write(output_unit,1)'11.0  47.4  Schneefernerhaus  # Schneefernerhaus (Zugspitze)                                           '
      write(output_unit,1)'11.0  47.8  Hohenpeissenberg  # Hohenpeissenberg                                                       '
        case(9)  ! BLOCK 9: (Optional): NETCDF ANNOTATIONS
          call Write_input_block_header(output_unit,blockID)
          !call SetWrite_input_block_09(WriteBlock)
      write(output_unit,1)'******************* BLOCK 9 ***************************************************                        '
      write(output_unit,1)'3d_tephra_fall.nc             # Name of output file                                                    '
      write(output_unit,1)'Eyjafjallajokull              # Title of simulation                                                    '
      write(output_unit,1)'no comment                    # Comment                                                                '
        case(10)  ! BLOCK 10 (OPTMOD): Optional module blocks
          call Write_input_block_header(output_unit,blockID)
          !call SetWrite_input_block_ResetParam(WriteBlock)
                  !   First RESETPARAMS
      write(output_unit,1)'******************* BLOCK 10+ *************************************************                        '
      write(output_unit,1)'OPTMOD=RESETPARAMS                                                                                     '
      write(output_unit,1)' MagmaDensity         = 2500.0                                                                         '
      write(output_unit,1)' DepositDensity       = 1000.0                                                                         '
      write(output_unit,1)' LAM_GS_THRESH        = 250.0                                                                          '
      write(output_unit,1)' AIRBORNE_THRESH      = 1.0e-3                                                                         '
      write(output_unit,1)' GRAV                 = 9.81                                                                           '
      write(output_unit,1)' RAD_EARTH            = 6371.229                                                                       '
      write(output_unit,1)' CFL                  = 0.80                                                                           '
      write(output_unit,1)' DT_MIN               = 1.0e-5                                                                         '
      write(output_unit,1)' DT_MAX               = 1.0                                                                            '
      write(output_unit,1)' ZPADDING             = 1.3                                                                            '
      write(output_unit,1)' DEPO_THRESH          = 1.0e-1                                                                         '
      write(output_unit,1)' DEPRATE_THRESH       = 1.0e-2                                                                         '
      write(output_unit,1)' CLOUDCON_THRESH      = 1.0e-3                                                                         '
      write(output_unit,1)' CLOUDCON_GRID_THRESH = 1.0e-7                                                                         '
      write(output_unit,1)' CLOUDLOAD_THRESH     = 1.0e-2                                                                         '
      write(output_unit,1)' THICKNESS_THRESH     = 1.0e-3                                                                         '
      write(output_unit,1)' StopValue_FracAshDep = 0.99                                                                           '
      write(output_unit,1)' DBZ_THRESH           = -2.0e+1                                                                        '
      write(output_unit,1)' Imp_fac              = 0.5                                                                            '
      write(output_unit,1)' Imp_DT_fac           = 4.0                                                                            '
      write(output_unit,1)' VelMod_umb           = 1                                                                              '
      write(output_unit,1)' lambda_umb           = 0.2                                                                            '
      write(output_unit,1)' N_BV_umb             = 0.02                                                                           '
      write(output_unit,1)' k_entrainment_umb    = 0.1                                                                            '
      write(output_unit,1)' SuzK_umb             = 12.0                                                                           '
      write(output_unit,1)' useMoistureVars      = F                                                                              '
      write(output_unit,1)' useWindVars          = F                                                                              '
      write(output_unit,1)' useOutprodVars       = T                                                                              '
      write(output_unit,1)' useRestartVars       = T                                                                              '
      write(output_unit,1)' useVz_rhoG           = T                                                                              '
      write(output_unit,1)' cdf_institution      = USGS                                                                           '
      write(output_unit,1)' cdf_run_class        = Analysis                                                                       '
      write(output_unit,1)' cdf_url              = https://vsc-ash.wr.usgs.gov/ash3d-gui                                          '
                 !   TOPO
          call Write_input_block_header(output_unit,blockID+1)
          !call SetWrite_input_block_Topo(WriteBlock)
      write(output_unit,1)'******************* BLOCK 10+ *************************************************                        '
      write(output_unit,1)'OPTMOD=TOPO                                                                                            '
      write(output_unit,1)'yes 2                         # use topography?; z-mod (0=none,1=shift,2=sigma)                        '
      write(output_unit,1)'1 1.0                         # Topofile format, smoothing radius                                      '
      write(output_unit,1)'GEBCO_2023.nc                 # topofile name                                                          '

                 !   VARDIFF
          call Write_input_block_header(output_unit,blockID+2)
          !call SetWrite_input_block_VarDiff(WriteBlock)
      write(output_unit,1)'******************* BLOCK 10+ *************************************************                        '
      write(output_unit,1)'OPTMOD=VARDIFF                                                                                         '
      write(output_unit,1)'yes 2 0.2                   # use horizontal variable diffusivity                                      '
      write(output_unit,1)'yes                         # use vertical variable diffusivity                                        '
      write(output_unit,1)'4                           # boundary layer model                                                     '
      write(output_unit,1)'3                           # free-air model                                                           '
      write(output_unit,1)'0.4                         # vonKarman                                                                '
      write(output_unit,1)'150.0                       # LambdaC                                                                  '
      write(output_unit,1)'0.25                        # RI_CRIT                                                                  '
      write(output_unit,1)'*******************************************************************************                        '

        case default
      do io=1,2;if(VB(io) <= verbosity_error)then
        write(errlog(io),*) &
          "WARNING: Block Header Number Not Specified"
      endif;enddo

      end select

 1    format(a103)

      end subroutine help_inputfile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  help_postproc
!
!  Called from: Parse_Command_Line or Ash3d_PostProc.f90
!  Arguments:
!    none
!
!  This subroutine writes the usage of the tool Ash3d_PostProc to the screen.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine help_postproc

      use io_units

      implicit none
      !implicit none (type, external)

#ifdef WINDOWS
      integer           :: iostatus
      character(len=120):: iomessage
      character(len=1)  :: key
#endif

      do io=1,2;if(VB(io) <= verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine help_postproc"
      endif;enddo

      io = 1

      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)' Ash3d post-processing tool: Ash3d_PostProc                                     '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)'Usage: Ash3d_PostProc control_file [t_index]                                    '
      write(outlog(io),1)'           or                                                                   '
      write(outlog(io),1)'       Ash3d_PostProc infile output_product format                              '
      write(outlog(io),1)'  where: infile   = the netcdf file written by Ash3d                            '
      write(outlog(io),1)'   output_product = 1 full concentration array                                  '
      write(outlog(io),1)'                    2 deposit granularity                                       '
      write(outlog(io),1)'                    3 deposit thickness (mm time-series)                        '
      write(outlog(io),1)'                    4 deposit thickness (inches time-series)                    '
      write(outlog(io),1)'                    5 deposit thickness (mm final)                              '
      write(outlog(io),1)'                    6 deposit thickness (inches final)                          '
      write(outlog(io),1)'                    7 ashfall arrival time (hours)                              '
      write(outlog(io),1)'                    8 ashfall arrival at airports/POI (mm)                      '
      write(outlog(io),1)'                    9 ash-cloud concentration (mg/m3)                           '
      write(outlog(io),1)'                   10 ash-cloud height (km)                                     '
      write(outlog(io),1)'                   11 ash-cloud bottom (km)                                     '
      write(outlog(io),1)'                   12 ash-cloud load (T/km2 or )                                '
      write(outlog(io),1)'                   13 ash-cloud radar reflectivity (dBz)                        '
      write(outlog(io),1)'                   14 ash-cloud arrival time (hours)                            '
      write(outlog(io),1)'                   15 topography                                                '
      write(outlog(io),1)'                   16 profile plots                                             '
      write(outlog(io),1)'           format = 1 ASCII/ArcGIS                                              '
      write(outlog(io),1)'                    2 KML/KMZ                                                   '
      write(outlog(io),1)'                    3 image/png                                                 '
      write(outlog(io),1)'                    4 binary                                                    '
      write(outlog(io),1)'                    5 shape file                                                '
       ! Place-holders for planned data formats
      !write(outlog(io),1)'                    6 grib2                                                     '
      !write(outlog(io),1)'                    7 netcdf                                                    '
      !write(outlog(io),1)'                    8 tecplot                                                   '
      !write(outlog(io),1)'                    9 vtk                                                       '
      write(outlog(io),1)'         [t_index] = index of time slice to plot; -1 for final (optional)       '

      ! For windows systems, allow the console to remain on the screen before exiting.
#ifdef WINDOWS
      write(outlog(io),*)'Press any key to exit.'
      read(input_unit,'(a1)',iostat=iostatus,iomsg=iomessage) key
#endif
      stop 1

 1    format(a80)

      end subroutine help_postproc

!##############################################################################
