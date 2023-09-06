!##############################################################################
!
!  Help module
!
!  This module contains all the subroutines for provided user-requested help
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

      module help

      use io_units

      implicit none

        ! Set everything to private by default
      private

       ! Publicly available subroutines/functions
      public help_general,   &
             help_make,      &
             help_run,       &
             help_input,     &
             help_inputfile, &
             help_postproc

      contains
      !------------------------------------------------------------------------

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

      io = 1

      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)' Ash3d help                                                                     '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)' You can get more specific help by typing:                                      '
      write(outlog(io),1)'   Ash3d -h make     : Information on build options, if you want to rebuild     '
      write(outlog(io),1)'   Ash3d -h run      : Information on running Ash3d from the command line       '
      write(outlog(io),1)'   Ash3d -h input    : Information on the structure of the input file           '
      write(outlog(io),1)'   Ash3d -h postproc : Information on the post-processing output results        '
      write(outlog(io),1)'                                                                                '
      write(outlog(io),1)' Writing run help information:                                                       '
      write(outlog(io),1)'                                                                                '

      call help_run

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

      io = 1

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

      character(len=3)  :: answer
      integer           :: blockID

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
      write(outlog(io),1)' BLOCK 10: OPTIONAL MODULES                                                     '
      write(outlog(io),1)'Would you like information on the structure of these blocks? (y or n)           '

      read(input_unit,'(a3)') answer
      do while (answer.eq.'y'.or.answer.eq.'yes')
        write(outlog(io),*)'  Please enter the block ID (1-10):'
        read(input_unit,*,err=10)blockID
        if(blockID.lt.1.or.blockID.gt.10)then
          write(errlog(io),*)'  Invalid range for blockID; should be between 1 and 10.'
        else
         call help_inputfile(blockID)
        endif
        write(outlog(io),1)'                                                                                '
        write(outlog(io),1)'Would you like information on another block? (y or n)                           '
        read(input_unit,'(a3)') answer
      enddo 

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

      integer,intent(in) :: blockID

      ! The idea with the blockID is that help for only a particular block is
      ! called if there is an error reading something in the input file or if
      ! the user requests it.

      io = 1

      select case (blockID)
        case(1) ! BLOCK 1: GRID INFO
      write(outlog(io),1)'# The following is an input file to the model Ash3d, v.1.0                                             '
      write(outlog(io),1)'# Created by L.G. Mastin, R.P. Denlinger and H.F. Schwaiger, U.S. Geological Survey, 2009.             '
      write(outlog(io),1)'#                                                                                                      '
      write(outlog(io),1)'# GENERAL SOURCE PARAMETERS. DO NOT DELETE ANY NON-COMMENT LINES                                       '
      write(outlog(io),1)'#  The first line of this block identifies the volcano by name.                                        '
      write(outlog(io),1)'#  If the volcano name begins with either 0 or 1, then the volcano                                     '
      write(outlog(io),1)'#  is assumed to be in the Smithonian database and default values for                                  '
      write(outlog(io),1)'#  Plume Height, Duration, Mass Flux Rate, Volume, and mass fraction of                                '
      write(outlog(io),1)'#  fines are loaded.  These can be over-written by entering non-negative                               '
      write(outlog(io),1)'#  values in the appropriate locations in this input file.                                             '
      write(outlog(io),1)'#                                                                                                      '
      write(outlog(io),1)'#  The second line of this block identifies the projection used and the form of                        '
      write(outlog(io),1)'#  the input coordinates and is of the following format:                                               '
      write(outlog(io),1)'#    latlonflag, projflag,  followed by a variable list of projection parameters                       '
      write(outlog(io),1)'#  projflag describes the projection used for the Ash3d run. Windfiles can have a                      '
      write(outlog(io),1)'#  different projection.                                                                               '
      write(outlog(io),1)'#  For a particular projflag, additional values are read defining the projection.                      '
      write(outlog(io),1)'#    latlonflag = 0 if computational grid is projected                                                 '
      write(outlog(io),1)'#               = 1 if copmutational grid is lat/lon (all subsequent projection parameters ignored.)   '
      write(outlog(io),1)'#    projflag   = 1 -- polar stereographic projection                                                  '
      write(outlog(io),1)'#           lambda0 -- longitude of projection point                                                   '
      write(outlog(io),1)'#           phi0    -- latitude of projection point                                                    '
      write(outlog(io),1)'#           k0      -- scale factor at projection point                                                '
      write(outlog(io),1)'#           radius  -- earth radius for spherical earth                                                '
      write(outlog(io),1)'#     e.g. for NAM 104,198, 216: 0 1 -105.0 90.0 0.933 6371.229                                        '
      write(outlog(io),1)'#               = 2 -- Alberts Equal Area ( not yet verified !!)                                       '
      write(outlog(io),1)'#               = 3 -- UTM ( not yet verified !!)                                                      '
      write(outlog(io),1)'#               = 4 -- Lambert conformal conic                                                         '
      write(outlog(io),1)'#           lambda0 -- longitude of origin                                                             '
      write(outlog(io),1)'#              phi0 -- latitude of origin                                                              '
      write(outlog(io),1)'#              phi1 -- latitude of secant1                                                             '
      write(outlog(io),1)'#              phi2 -- latitude of secant2                                                             '
      write(outlog(io),1)'#            radius -- earth radius for a spherical earth                                              '
      write(outlog(io),1)'#     e.g. for NAM 212: 0 4 265.0 25.0 25.0 25.0 6371.22                                               '
      write(outlog(io),1)'#               = 5 -- Mercator                                                                        '
      write(outlog(io),1)'#           lambda0 -- longitude of origin                                                             '
      write(outlog(io),1)'#              phi0 -- latitude of origin                                                              '
      write(outlog(io),1)'#            radius -- earth radius for a spherical earth                                              '
      write(outlog(io),1)'#     e.g. for NAM 196: 0 5 198.475 20.0 6371.229                                                      '
      write(outlog(io),1)'#                                                                                                      '
      write(outlog(io),1)'# On line 3, the vent coordinates can optionally include a third value for elevation in km.            '
      write(outlog(io),1)'# If the vent elevation is not given, 0 is used if topography is turned off.                           '
      write(outlog(io),1)'#                                                                                                      '
      write(outlog(io),1)'# Line 4 is the width and height of the computational grid in km (if projected) or degrees.            '
      write(outlog(io),1)'# Line 5 is the vent x,y (or lon, lat) coordinates.                                                    '
      write(outlog(io),1)'# Line 6, DX and DY resolution in km or degrees (for projected or lon/lat grid, respectively)          '
      write(outlog(io),1)'# Line 7, DZ can be given as a real number, indicating the vertical spacing in km.                     '
      write(outlog(io),1)'# Alternatively, it can be given as dz_plin (piece-wise linear), dz_clog (constant-                    '
      write(outlog(io),1)'# logarithmic), or dz_cust (custom specification)                                                      '
      write(outlog(io),1)'# If dz_plin, then a second line is read containing:                                                   '
      write(outlog(io),1)'#   number of line segments (N) followed by the steps and step-size of each segment                    '
      write(outlog(io),1)'#   e.g. 4 6 0.25 5 0.5 5 1.0 10 2.0                                                                   '
      write(outlog(io),1)'#         This corresponds to 4 line segments with 6 cells of 0.25, then 5 cells of 0.5,               '
      write(outlog(io),1)'#         5 cells of 1.0, and finally 10 cells of 2.0                                                  '
      write(outlog(io),1)'# If dz_clog, then a second line is read containing:                                                   '
      write(outlog(io),1)'#   maximum z and number of steps of constant dlogz                                                    '
      write(outlog(io),1)'#   e.g. 30.0 30                                                                                       '
      write(outlog(io),1)'#         This corresponds to 30 steps from 0-30km with constant log-spacing                           '
      write(outlog(io),1)'# If dz_cust, then a second line is read containing:                                                   '
      write(outlog(io),1)'#   the number of dz values to read (ndz), followed by dz(1:ndz)                                       '
      write(outlog(io),1)'#   e.g. 20 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 5.5            '
      write(outlog(io),1)'#         This corresponds to 10 steps of 0.5, 9 steps of 1.5, followed by 1 step of 5.5               '
      write(outlog(io),1)'#                                                                                                      '
      write(outlog(io),1)'#                                                                                                      '
      write(outlog(io),1)'# Line 8 is the the diffusivity (m2/s) followed by the eruption specifier.  The                        '
      write(outlog(io),1)'# eruption specifier can be a real number, in which case it is assumed to be the                       '
      write(outlog(io),1)'# positive constant specifying the Suzuki distribution.  Alternatively, it can be                      '
      write(outlog(io),1)'#  umbrella     : Suzuki (const. = 12) with radial spreading of the plume                              '
      write(outlog(io),1)'#  umbrella_air : Suzuki (const. = 12) with radial spreading of the plume scaled to 5% of vol.         '
      write(outlog(io),1)'#  point        : all mass inserted in cell containing PlmH                                            '
      write(outlog(io),1)'#  linear       : mass uniformly distributed from z-vent to PlmH                                       '
      write(outlog(io),1)'# Line 9 : number of pulses to be read in BLOCK 2                                                      '
      write(outlog(io),1)'#'
      write(outlog(io),1)'******************* BLOCK 1 ***************************************************                        '
      write(outlog(io),1)'Eyjafjallajokull               # Volcano name (character*30)                                           '
      write(outlog(io),1)'1 1 0.0 90.0 0.933 6367.470    # Proj flags and params; first term (LLflag) is 1, so all else ignored  '
      write(outlog(io),1)'-25.0    45.0                  # x, y of LL corner of grid (km, or deg. if latlongflag=1)              '
      write(outlog(io),1)'55.0     25.0                  # grid width and height (km, or deg. if latlonflag=1)                   '
      write(outlog(io),1)'-19.62   63.63                 # vent location         (km, or deg. if latlonflag=1)                   '
      write(outlog(io),1)'3.0      3.0                   # DX, DY of grid cells  (km, or deg. if latlonflag=1)                   '
      write(outlog(io),1)'1.0                            # DZ of grid cells      (always km)                                     '
      write(outlog(io),1)'0.     4.                      # diffusion coefficient (m2/s), Suzuki constant                         '
      write(outlog(io),1)'9                              # neruptions, number of eruptions or pulses                             '
      write(outlog(io),1)'*******************************************************************************                        '
        case(2) ! BLOCK 2: ERUPTION PARAMETERS
      write(outlog(io),1)'# ERUPTION LINES (number = neruptions)                                                                 '
      write(outlog(io),1)'# In the following line, each line represents one eruptive pulse.                                      '
      write(outlog(io),1)'# Parameters are (1-4) start time (yyyy mm dd h.hh (UT)); (5) duration (hrs);                          '
      write(outlog(io),1)'#                  (6) plume height;                      (7) eruped volume (km3 DRE)                  '
      write(outlog(io),1)'# If neruptions=1 and the year is 0, then the model run in forecast mode where mm dd h.hh are          '
      write(outlog(io),1)'# interpreted as the time after the start of the windfile.  In this case, duration, plume              '
      write(outlog(io),1)'# height and erupted volume are replaced with ESP if the values are negative.                          '
      write(outlog(io),1)'# This applies to source types: suzuki, point, line, umbrella and umbrella_air.                        '
      write(outlog(io),1)'# For profile sources, an additional two values are read: dz and nz                                    '
      write(outlog(io),1)'# 2010 04 14   0.00   1.0     18.0  0.16 1.0 18                                                        '
      write(outlog(io),1)'# 0.01 0.02 0.03 0.03 0.04 0.04 0.05 0.06 0.06 0.070 0.08 0.08 0.09 0.09 0.09 0.08 0.06 0.02           '
      write(outlog(io),1)'******************* BLOCK 2 ***************************************************                        '
      write(outlog(io),1)'2010 4 14  9.0  3.0 7.4 8.46E-004                                                                      '
      write(outlog(io),1)'2010 4 14 12.0  3.0 8.4 1.38E-003                                                                      '
      write(outlog(io),1)'2010 4 14 15.0  3.0 6.8 6.11E-004                                                                      '
      write(outlog(io),1)'2010 4 14 18.0  3.0 5.6 2.89E-004                                                                      '
      write(outlog(io),1)'2010 4 14 21.0  3.0 5.1 2.01E-004                                                                      '
      write(outlog(io),1)'2010 4 15  0.0  3.0 5.3 2.33E-004                                                                      '
      write(outlog(io),1)'2010 4 15  3.0  3.0 5.2 2.17E-004                                                                      '
      write(outlog(io),1)'2010 4 15  6.0  3.0 5.3 2.33E-004                                                                      '
      write(outlog(io),1)'2010 4 15  9.0  3.0 5.7 3.09E-004                                                                      '
        case(3) ! BLOCK 3: WIND PARAMETERS
      write(outlog(io),1)'*******************************************************************************                        '
      write(outlog(io),1)'# WIND OPTIONS                                                                                         '
      write(outlog(io),1)'# Ash3d will read from either a single 1-D wind sounding, or gridded, time-                            '
      write(outlog(io),1)'# dependent 3-D wind data, depending on the value of the parameter iwind.                              '
      write(outlog(io),1)'# For iwind = 1, read from a 1-D wind sounding                                                         '
      write(outlog(io),1)'#             2, read from 3D gridded ASCII files                                                      '
      write(outlog(io),1)'#             3/4, read directly from a single or multiple NetCDF files.                               '
      write(outlog(io),1)'#             5, read directly from multiple multi-timestep NetCDF files.                              '
      write(outlog(io),1)'# The parameter iwindformat specifies the format of the wind files, as follows:                        '
      write(outlog(io),1)'#  iwindformat =  0: User-defined via template                                                         '
      write(outlog(io),1)'#                 1: User-specified ASCII files                                                        '
      write(outlog(io),1)'#                 2: Global radiosonde data                                                            '
      write(outlog(io),1)'#                 3: NARR 221 Reanalysis (32 km)                                                       '
      write(outlog(io),1)'#                 4: NAM Regional North America 221 Forecast (32 km)                                   '
      write(outlog(io),1)'#                 5: NAM 216 Regional Alaska Forecast (45 km)                                          '
      write(outlog(io),1)'#                 6: NAM 104 Northern Hemisphere Forecast (90 km)                                      '
      write(outlog(io),1)'#                 7: NAM 212 40km Cont. US Forecast (40 km)                                            '
      write(outlog(io),1)'#                 8: NAM 218 12km Cont. US Forecast (12 km)                                            '
      write(outlog(io),1)'#                 9: NAM 227 Cont. US Forecast (5.08 km)                                               '
      write(outlog(io),1)'#                10: NAM 242 11km Regional Alaska Forecast (11.25 km)                                  '
      write(outlog(io),1)'#                11: NAM 196 Regional Hawaii Forecast (2.5 km)                                         '
      write(outlog(io),1)'#                12: NAM 198 Regional Alaska Forecast (5.953 km)                                       '
      write(outlog(io),1)'#                13: NAM 91 Regional Alaska Forecast (2.976 km)                                        '
      write(outlog(io),1)'#                14: NAM Regional Cont. US Forecast (3.0 km)                                           '
      write(outlog(io),1)'#                20: GFS 0.5 degree files Forecast                                                     '
      write(outlog(io),1)'#                21: GFS 1.0 degree files Forecast                                                     '
      write(outlog(io),1)'#                22: GFS 0.25 degree files Forecast                                                    '
      write(outlog(io),1)'#                23: NCEP DOE Reanalysis 2.5 degree                                                    '
      write(outlog(io),1)'#                24: NASA MERRA-2 Reanalysis                                                           '
      write(outlog(io),1)'#                25: NCEP1 2.5 global Reanalysis (1948-pres)                                           '
      write(outlog(io),1)'#                      Note: use nWindFiles=1 for iwindformat=25                                       '
      write(outlog(io),1)'#                26: JRA-55 Reanalysis                                                                 '
      write(outlog(io),1)'#                27: NOAA-CIRES II 2-deg global Reanalysis (1870-2010)                                 '
      write(outlog(io),1)'#                28: ECMWF ERA-Interim Reanalysis                                                      '
      write(outlog(io),1)'#                29: ECMWA ERA-5 Reanalysis                                                            '
      write(outlog(io),1)'#                30: ECMWA ERA-20C Reanalysis                                                          '
      write(outlog(io),1)'#                32: Air Force Weather Agency                                                          '
      write(outlog(io),1)'#                33: CCSM 3.0 Community Atmospheric Model                                              '
      write(outlog(io),1)'#                40: NASA GEOS-5 Cp                                                                    '
      write(outlog(io),1)'#                41: NASA GEOS-5 Np                                                                    '
      write(outlog(io),1)'#                50: Weather Research and Forecast (WRF) output                                        '
      write(outlog(io),1)'#                                                                                                      '
      write(outlog(io),1)'# igrid (optional, defaults to that associated with iwindformat) is the NCEP grid ID,                  '
      write(outlog(io),1)'# if a NWP product is used, or the number of stations of sonde data, if iwind = 1.                     '
      write(outlog(io),1)'# idata (optional, defaults to 2) is a flag for data type (1=ASCII, 2=netcdf, 3=grib).                 '
      write(outlog(io),1)'#                                                                                                      '
      write(outlog(io),1)'# Many plumes extend higher than the maximum height of mesoscale models.                               '
      write(outlog(io),1)'# Ash3d handles this as determined by the parameter iHeightHandler, as follows:                        '
      write(outlog(io),1)'# for iHeightHandler = 1, stop the program if the plume height exceeds mesoscale height                '
      write(outlog(io),1)'#                      2, wind velocity at levels above the highest node                               '
      write(outlog(io),1)'#                         equal that of the highest node.  Temperatures in the                         '
      write(outlog(io),1)'#                         upper nodes do not change between 11 and 20 km; above                        '
      write(outlog(io),1)'#                         20 km they increase by 2 C/km, as in the Standard                            '
      write(outlog(io),1)'#                         atmosphere.  A warning is written to the log file.                           '
      write(outlog(io),1)'# Simulation time in hours is the maximal length of the simulation.                                    '
      write(outlog(io),1)'# Ash3d can end the simulation early if desired, once 99% of the ash has deposited.                    '
      write(outlog(io),1)'# The last line of this block is the number of windfiles listed in block 5 below.  If                  '
      write(outlog(io),1)'# iwind=5 and one of the NWP products is used that require a special file structure,                   '
      write(outlog(io),1)'# then nWindFiles should be set to 1 and only the root folder of the windfiles listed.                 '
      write(outlog(io),1)'******************* BLOCK 3 ***************************************************                        '
      write(outlog(io),1)'4  20               # iwind, iwindformat, [igrid, idata]                                               '
      write(outlog(io),1)'2                   # iHeightHandler                                                                   '
      write(outlog(io),1)'60                  # Simulation time in hours                                                         '
      write(outlog(io),1)'no                  # stop computation when 99% of erupted mass has deposited?                         '
      write(outlog(io),1)'16                  # nWindFiles, number of gridded wind files (used if iwind>1)                       '
        case(4) ! BLOCK 4: OUTPUT OPTIONS
      write(outlog(io),1)'*******************************************************************************                        '
      write(outlog(io),1)'# OUTPUT OPTIONS:                                                                                      '
      write(outlog(io),1)'# The list below allows users to specify the output options                                            '
      write(outlog(io),1)'# All but the final deposit file can be written out at specified                                       '
      write(outlog(io),1)'# times using the following parameters:                                                                '
      write(outlog(io),1)'# Line 15 asks for 3d output (yes/no) followed by an optional output format code;                      '
      write(outlog(io),1)'#   1 = (default) output all the normal 2d products to the output file as well as the 3d concentrations'
      write(outlog(io),1)'#   2 = only output the 2d products                                                                    '
      write(outlog(io),1)'# nWriteTimes   = if >0,  number of times output are to be written. The following                      '
      write(outlog(io),1)'#                  line contains nWriteTimes numbers specifying the times of output                    '
      write(outlog(io),1)'#                 if =-1, it specifies that the following line gives a constant time                   '
      write(outlog(io),1)'#                  interval in hours between write times.                                              '
      write(outlog(io),1)'# WriteTimes    = Hours between output (if nWritetimes=-1), or                                         '
      write(outlog(io),1)'#                 Times (hours since start of first eruption) for each output                          '
      write(outlog(io),1)'#                (if nWriteTimes >1)                                                                   '
      write(outlog(io),1)'******************* BLOCK 4 ***************************************************                        '
      write(outlog(io),1)'no      # Write out ESRI ASCII file of final deposit thickness?                                        '
      write(outlog(io),1)'yes     # Write out        KML file of final deposit thickness?                                        '
      write(outlog(io),1)'no      # Write out ESRI ASCII deposit files at specified times?                                       '
      write(outlog(io),1)'no      # Write out        KML deposit files at specified times?                                       '
      write(outlog(io),1)'no      # Write out ESRI ASCII files of ash-cloud concentration?                                       '
      write(outlog(io),1)'no      # Write out        KML files of ash-cloud concentration ?                                      '
      write(outlog(io),1)'no      # Write out ESRI ASCII files of ash-cloud height?                                              '
      write(outlog(io),1)'no      # Write out        KML files of ash-cloud height?                                              '
      write(outlog(io),1)'yes     # Write out ESRI ASCII files of ash-cloud load (T/km2) at specified times?                     '
      write(outlog(io),1)'yes     # Write out        KML files of ash-cloud load (T/km2) at specified times?                     '
      write(outlog(io),1)'yes     # Write out ESRI ASCII file of deposit arrival times?                                          '
      write(outlog(io),1)'yes     # Write out        KML file of deposit arrival times?                                          '
      write(outlog(io),1)'yes     # write out ESRI ASCII file of cloud arrival times?                                            '
      write(outlog(io),1)'yes     # Write out        KML file of cloud arrival times?                                            '
      write(outlog(io),1)'yes 1   # Write out 3-D ash concentration at specified times? / [output code: 1=2d+concen,2=2d only]   '
      write(outlog(io),1)'netcdf  # format of ash concentration files     (ascii, binary, or netcdf)                             '
      write(outlog(io),1)'-1      # nWriteTimes                                                                                  '
      write(outlog(io),1)'1       # WriteTimes (hours since eruption start)                                                      '
        case(5) ! BLOCK 5: INPUT WIND FILES
      write(outlog(io),1)'*******************************************************************************                        '
      write(outlog(io),1)'# WIND INPUT FILES                                                                                     '
      write(outlog(io),1)'# The following block of data contains names of wind files.                                            '
      write(outlog(io),1)'# If we are reading from a 1-D wind sounding (i.e. iwind=1) then there should                          '
      write(outlog(io),1)'# be only one wind file.                                                                               '
      write(outlog(io),1)'# If we are reading gridded data there should be iWinNum wind files, each having                       '
      write(outlog(io),1)'# the format volcano_name_yyyymmddhh_FHhh.win                                                          '
      write(outlog(io),1)'# For iwindformat=25 (continuous wind data from 1948-pres), just give                                  '
      write(outlog(io),1)'# the directory with the windfiles (e.g. Wind_nc)                                                      '
      write(outlog(io),1)'******************* BLOCK 5 ***************************************************                        '
      write(outlog(io),1)'Wind_nc/gfs/gfs.2010041400/2010041400.f000.nc                                                          '
      write(outlog(io),1)'Wind_nc/gfs/gfs.2010041400/2010041400.f003.nc                                                          '
      write(outlog(io),1)'Wind_nc/gfs/gfs.2010041400/2010041400.f006.nc                                                          '
      write(outlog(io),1)'Wind_nc/gfs/gfs.2010041400/2010041400.f009.nc                                                          '
      write(outlog(io),1)'Wind_nc/gfs/gfs.2010041400/2010041400.f012.nc                                                          '
      write(outlog(io),1)'Wind_nc/gfs/gfs.2010041400/2010041400.f015.nc                                                          '
      write(outlog(io),1)'Wind_nc/gfs/gfs.2010041400/2010041400.f018.nc                                                          '
      write(outlog(io),1)'Wind_nc/gfs/gfs.2010041400/2010041400.f021.nc                                                          '
      write(outlog(io),1)'Wind_nc/gfs/gfs.2010041400/2010041400.f024.nc                                                          '
      write(outlog(io),1)'Wind_nc/gfs/gfs.2010041400/2010041400.f027.nc                                                          '
      write(outlog(io),1)'Wind_nc/gfs/gfs.2010041400/2010041400.f030.nc                                                          '
      write(outlog(io),1)'Wind_nc/gfs/gfs.2010041400/2010041400.f033.nc                                                          '
      write(outlog(io),1)'Wind_nc/gfs/gfs.2010041400/2010041400.f036.nc                                                          '
      write(outlog(io),1)'Wind_nc/gfs/gfs.2010041400/2010041400.f039.nc                                                          '
      write(outlog(io),1)'Wind_nc/gfs/gfs.2010041400/2010041400.f042.nc                                                          '
      write(outlog(io),1)'Wind_nc/gfs/gfs.2010041400/2010041400.f045.nc                                                          '
        case(6) ! BLOCK 6: AIRPORT FILE
      write(outlog(io),1)'*******************************************************************************                        '
      write(outlog(io),1)'# AIRPORT LOCATION FILE                                                                                '
      write(outlog(io),1)'# The following lines allow the user to specify whether times of ash arrival                           '
      write(outlog(io),1)'# at airports and other locations will be written out, and which file                                  '
      write(outlog(io),1)'# to read for a list of airport locations.                                                             '
      write(outlog(io),1)'# PLEASE NOTE:  Each line in the airport location file should contain the                              '
      write(outlog(io),1)'#               airport latitude, longitude, projected x and y coordinates,                            '
      write(outlog(io),1)'#               and airport name.  If you are using a projected grid,                                  '
      write(outlog(io),1)'#               THE X AND Y MUST BE IN THE SAME PROJECTION as the computational grid.                  '
      write(outlog(io),1)'#               Alternatively, if coordinates can be projected via libprojection                       '
      write(outlog(io),1)'#               by typing "yes" to the last parameter                                                  '
      write(outlog(io),1)'******************* BLOCK 6 ***************************************************                        '
      write(outlog(io),1)'no                            # Write out ash arrival times at airports to ASCII FILE?                 '
      write(outlog(io),1)'no                            # Write out grain-size distribution to ASCII airports file?              '
      write(outlog(io),1)'no                            # Write out ash arrival times to kml file?                               '
      write(outlog(io),1)'GlobalAirports.txt            # Name of file containing aiport locations                               '
      write(outlog(io),1)'no                            # Have libprojection calculate projected coordinates?                    '
        case(7) ! BLOCK 7: GRAIN-SIZE BINS, SETTLING VELOCITY
      write(outlog(io),1)'*******************************************************************************                        '
      write(outlog(io),1)'#GRAIN SIZE GROUPS                                                                                     '
      write(outlog(io),1)'******************* BLOCK 7 ***************************************************                        '
      write(outlog(io),1)'# The first line must contain the number of settling velocity groups, but                              '
      write(outlog(io),1)'# can optionally also include a flag for the fall velocity model to be used.                           '
      write(outlog(io),1)'#    FV_ID = 1, Wilson and Huang                                                                       '
      write(outlog(io),1)'#          = 2, Wilson and Huang + Cunningham slip                                                     '
      write(outlog(io),1)'#          = 3, Wilson and Huang + Mod by Pfeiffer Et al.                                              '
      write(outlog(io),1)'#          = 4, Ganser (assuming prolate ellipsoids)                                                   '
      write(outlog(io),1)'#          = 5, Stokes flow for spherical particles + slip                                             '
      write(outlog(io),1)'# If no fall model is specified, FV_ID = 1, by default                                                 '
      write(outlog(io),1)'# The grain size bins can be enters with 2, 3, or 4 parameters.                                        '
      write(outlog(io),1)'# If TWO are given, they are read as:   FallVel (in m/s), mass fraction                                '
      write(outlog(io),1)'# If THREE are given, they are read as: diameter (mm), mass fraction, density (kg/m3)                  '
      write(outlog(io),1)'# If FOUR are given, they are read as:  diameter (mm), mass fraction, density (kg/m3), Shape F         '
      write(outlog(io),1)'# The shape factor is given as in Wilson and Huang: F=(b+c)/(2*a), but converted                       '
      write(outlog(io),1)'# to sphericity (assuming b=c) for the Ganser model.                                                   '
      write(outlog(io),1)'# If a shape factor is not given, a default value of F=0.4 is used.                                    '
      write(outlog(io),1)'# If FIVE are given, they are read as:  diameter (mm), mass fraction, density (kg/m3), Shape F, G      '
      write(outlog(io),1)'#  where G is an additional Ganser shape factor equal to c/b                                           '
      write(outlog(io),1)'#                                                                                                      '
      write(outlog(io),1)'# If the last grain size bin has a negative diameter, then the remaining mass fraction                 '
      write(outlog(io),1)'# will be distributed over the previous bins via a log-normal distribution in phi.                     '
      write(outlog(io),1)'# The last bin would be interpreted as:                                                                '
      write(outlog(io),1)'# diam (neg value) , phi_mean, phi_stddev                                                              '
      write(outlog(io),1)'#                              # nsmax FV_ID                                                           '
      write(outlog(io),1)'*******************************************************************************                        '
      write(outlog(io),1)'15                                 # Number of grain-size bins. FV_ID not given; defaults to 1         '
      write(outlog(io),1)'2.000      0.0208  2003.   0.44    # grain size (mm), mass fraction, density (kg/m3), F=0.44           '
      write(outlog(io),1)'1.414      0.0084  2350.   0.44                                                                        '
      write(outlog(io),1)'1.000      0.0141  2005.   0.44                                                                        '
      write(outlog(io),1)'0.707      0.0214  2248.   0.44                                                                        '
      write(outlog(io),1)'0.500      0.0459  2624.   0.44                                                                        '
      write(outlog(io),1)'0.354      0.0723  2644.   0.44                                                                        '
      write(outlog(io),1)'0.250      0.0532  2639.   0.44                                                                        '
      write(outlog(io),1)'0.177      0.0219  2690.   0.44                                                                        '
      write(outlog(io),1)'0.125      0.0165  2691.   0.44                                                                        '
      write(outlog(io),1)'0.088      0.0115  2730.   0.44                                                                        '
      write(outlog(io),1)'0.2176     0.0714  600.    1.00    # Note that these are bigger again, low density and round           '
      write(outlog(io),1)'0.2031     0.1428  600.    1.00    # These bottom five bins represent an aggregate distribution        '
      write(outlog(io),1)'0.1895     0.2856  600.    1.00                                                                        '
      write(outlog(io),1)'0.1768     0.1428  600.    1.00                                                                        '
      write(outlog(io),1)'0.1649     0.0714  600.    1.00                                                                        '
        case(8) ! BLOCK 8: VERTICAL PROFILES
      write(outlog(io),1)'*******************************************************************************                        '
      write(outlog(io),1)'# Options for writing vertical profiles                                                                '
      write(outlog(io),1)'# The first line below gives the number of locations (nlocs) where vertical                            '
      write(outlog(io),1)'# profiles are to be written.  That is followed by nlocs lines, each of which                          '
      write(outlog(io),1)'# contain the location, in the same coordinate as the computational grid.                              '
      write(outlog(io),1)'# Optionally, a site name can be provided in after the location.                                       '
      write(outlog(io),1)'******************* BLOCK 8 ***************************************************                        '
      write(outlog(io),1)'4                             # number of locations for vertical profiles (nlocs)                      '
      write(outlog(io),1)'12.4  51.4  Leipzig           # x,y (or lon/lat) [Site name]                                           '
      write(outlog(io),1)'11.3  48.2  Munich            # Munich (Maisach)                                                       '
      write(outlog(io),1)'11.0  47.4  Schneefernerhaus  # Schneefernerhaus (Zugspitze)                                           '
      write(outlog(io),1)'11.0  47.8  Hohenpeissenberg  # Hohenpeissenberg                                                       '
        case(9) ! BLOCK 9: (Optional): NETCDF ANNOTATIONS
      write(outlog(io),1)'*******************************************************************************                        '
      write(outlog(io),1)'# netCDF output options                                                                                '
      write(outlog(io),1)'# This last block is optional.                                                                         '
      write(outlog(io),1)'# The output file name can be give, but will default to 3d_tephra_fall.nc if absent                    '
      write(outlog(io),1)'# The title and comment lines are passed through to the netcdf header of the                           '
      write(outlog(io),1)'# output file.                                                                                         '
      write(outlog(io),1)'******************* BLOCK 9 ***************************************************                        '
      write(outlog(io),1)'3d_tephra_fall.nc             # Name of output file                                                    '
      write(outlog(io),1)'Eyjafjallajokull              # Title of simulation                                                    '
      write(outlog(io),1)'no comment                    # Comment                                                                '
        case(10) ! BLOCK 10 (OPTMOD): Optional module blocks
      write(outlog(io),1)'******************* BLOCK 10+ *************************************************                        '
      write(outlog(io),1)'# Optional Modules are identified by the text string at the top of the block                           '
      write(outlog(io),1)'# OPTMOD=[module name]                                                                                 '
      write(outlog(io),1)'# There will need to be a custom block reader in the module to read this section                       '
      write(outlog(io),1)'# section of the input file.  Below is the built-in example for resetting parameters                   '
      write(outlog(io),1)'*******************************************************************************                        '
      write(outlog(io),1)'OPTMOD=RESETPARAMS                                                                                     '
      write(outlog(io),1)' MagmaDensity         = 3500.0                                                                         '
      write(outlog(io),1)' DepositDensity       = 1300.0                                                                         '
      write(outlog(io),1)' LAM_GS_THRESH        = 250.0                                                                          '
      write(outlog(io),1)' AIRBORNE_THRESH      = 1.0e-3                                                                         '
      write(outlog(io),1)' GRAV                 = 9.81                                                                           '
      write(outlog(io),1)' RAD_EARTH            = 6371.229                                                                       '
      write(outlog(io),1)' CFL                  = 0.80                                                                           '
      write(outlog(io),1)' DT_MIN               = 1.0e-5                                                                         '
      write(outlog(io),1)' DT_MAX               = 1.0                                                                            '
      write(outlog(io),1)' ZPADDING             = 1.3                                                                            '
      write(outlog(io),1)' DEPO_THRESH          = 1.0e-1                                                                         '
      write(outlog(io),1)' DEPRATE_THRESH       = 1.0e-2                                                                         '
      write(outlog(io),1)' CLOUDCON_THRESH      = 2.0e-1                                                                         '
      write(outlog(io),1)' CLOUDCON_GRID_THRESH = 2.0e-1                                                                         '
      write(outlog(io),1)' CLOUDLOAD_THRESH     = 1.0e-2                                                                         '
      write(outlog(io),1)' THICKNESS_THRESH     = 1.0e-3                                                                         '
      write(outlog(io),1)' DBZ_THRESH           = -2.0e+1                                                                        '
      write(outlog(io),1)' VelMod_umb           = 1                                                                              '
      write(outlog(io),1)' lambda_umb           = 0.2                                                                            '
      write(outlog(io),1)' N_BV_umb             = 0.02                                                                           '
      write(outlog(io),1)' k_entrainment_umb    = 0.1                                                                            '
      write(outlog(io),1)' SuzK_umb             = 12.0                                                                           '
      write(outlog(io),1)' useMoistureVars      = F                                                                              '
      write(outlog(io),1)' useVz_rhoG           = T                                                                              '
      write(outlog(io),1)' cdf_institution      = USGS                                                                           '
      write(outlog(io),1)' cdf_run_class        = Analysis                                                                       '
      write(outlog(io),1)' cdf_url              = https://vsc-ash.wr.usgs.gov/ash3d-gui                                          '

!        case default
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
      write(outlog(io),1)'                    6 grib2                                                     '
      write(outlog(io),1)'                    7 netcdf                                                    '
      write(outlog(io),1)'                    8 tecplot                                                   '
      write(outlog(io),1)'                    9 vtk                                                       '
      write(outlog(io),1)'         [t_index] = index of time slice to plot; -1 for final (optional)       '

      stop 1

 1    format(a80)

      end subroutine help_postproc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module help

!##############################################################################

