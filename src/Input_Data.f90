!##############################################################################
!
!      Ash3d_Program_Control module
!
!  This module manages parsing of the command line, reading and vetting the
!  control file.
!
!      subroutine Parse_Command_Line
!      subroutine Set_OS_Env
!      subroutine check_endian
!      subroutine Read_Control_File
!      subroutine LatLonChecker
!      subroutine xyChecker
!      subroutine vprofchecker
!      subroutine Read_PostProc_Control_File
!      subroutine Write_input_block_header
!      subroutine SetWrite_input_block_01
!      subroutine SetWrite_input_block_02
!      subroutine SetWrite_input_block_03
!      subroutine SetWrite_input_block_04
!      subroutine SetWrite_input_block_05
!      subroutine SetWrite_input_block_06
!      subroutine SetWrite_input_block_07
!      subroutine SetWrite_input_block_08
!      subroutine SetWrite_input_block_09
!      subroutine SetWrite_input_block_ResetParam
!      subroutine SetWrite_input_block_Topo
!      subroutine SetWrite_input_block_VarDiff
!
!##############################################################################

      module Ash3d_Program_Control

      use precis_param

      use io_units

      implicit none

        ! Set everything to private by default
      private

        ! Publicly available subroutines/functions
      public Parse_Command_Line,         &
             Set_OS_Env,                 &
             Read_Control_File,          &
             Read_PostProc_Control_File, &
             Write_input_block_header,   &
             SetWrite_input_block_01,    &
             SetWrite_input_block_02,    &
             SetWrite_input_block_03,    &
             SetWrite_input_block_04,    &
             SetWrite_input_block_05,    &
             SetWrite_input_block_06,    &
             SetWrite_input_block_07,    &
             SetWrite_input_block_08,    &
             SetWrite_input_block_09,    &
             SetWrite_input_block_ResetParam,    &
             SetWrite_input_block_Topo,  &
             SetWrite_input_block_VarDiff

      contains
      !------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Subroutine Parse_Command_Line
!
!  Called from: Ash3d.F90
!  Arguments:
!    none
!
! This subroutine parses the Ash3d command line to determine if the user needs
! help, is running Ash3d interactively, or is running Ash3d via a control file.
! Returns variables: LoadConcen, infile
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Parse_Command_Line

      ! This module requires Fortran 2003 or later
      use iso_c_binding

      ! This module requires Fortran 2003 or later
      use iso_fortran_env, only : &
         input_unit

      use io_data,       only : &
         concenfile,infile,LoadConcen

#ifdef USENETCDF
      use Ash3d_Netcdf_IO
#endif

      integer           :: nargs          ! number of command-line arguments
      character(len=3)  :: answer
      character(len=130):: linebuffer130
      character         :: testkey,testkey2
      integer           :: arglen
      integer           :: iostatus
      character(len=120):: iomessage
      integer           :: blockID

      !character(kind=c_char), dimension(1:130) :: fc_inputfile
      character, dimension(1:130) :: fc_inputfile
      ! Size matches length of infile (specified in module io_data)
      integer fc_len

      INTERFACE
        subroutine help_general
        end subroutine help_general
        subroutine help_make
        end subroutine help_make
        subroutine help_run
        end subroutine help_run
        subroutine help_input
        end subroutine help_input
        subroutine help_postproc
        end subroutine help_postproc
      END INTERFACE

      ! Since we haven't opened a logfile yet, only write out to stdout if not a
      ! control file case.
      io = 1

      ! Test read command-line arguments
      nargs = command_argument_count()
      if(nargs.eq.0) then
          ! If no command-line arguments are given, then prompt user
          ! interactively for the command file name and possible a 
          ! restart file
        if(VB(1).ge.verbosity_silent)then
          write(errlog(io),*)"Stdout is suppressed via ASH3DVERB=9, but interactive input is expected."
          write(errlog(io),*)"Either recompile with ASH3DVERB<9, over-ride with the environment variable"
          write(errlog(io),*)"(ASH3DVERB) or provide the correct command-line arguments."
          stop 1
        else
          write(outlog(io),*)'Enter name of ESP input file:'
        endif
        read(input_unit,*,iostat=iostatus,iomsg=iomessage) infile 
        write(outlog(io),*)'Load concentration file?'
        read(input_unit,'(a3)',iostat=iostatus,iomsg=iomessage) answer
        if(adjustl(trim(answer)).eq.'y'.or.adjustl(trim(answer)).eq.'yes') then
          LoadConcen = .true.
        elseif(adjustl(trim(answer)).eq.'n'.or.adjustl(trim(answer)).eq.'no') then
          LoadConcen = .false.
        else
          write(errlog(io),*) 'Sorry, I cannot understand your answer.'
          write(errlog(io),*) "Expected either 'yes' or 'no', but you provided:",answer
          stop 1
        endif
        if(LoadConcen)then
          ! We are initializing the concentration and time from an output file
          ! Currently, Ash3d assumes the concentration file is compatible with
          ! the computational grid and grainsize distribution
          write(outlog(io),*)'Enter name of concentration file'
          read(input_unit,*,iostat=iostatus,iomsg=iomessage) concenfile
#ifdef USENETCDF
          call NC_RestartFile_ReadTimes
#else
          write(errlog(io),*)"ERROR: ",&
           "Loading concentration files requires previous netcdf"
          write(errlog(io),*)&
           "       output.  This Ash3d executable was not compiled with"
          write(errlog(io),*)&
           "       netcdf support.  Please recompile Ash3d with"
          write(errlog(io),*)&
           "       USENETCDF=T, or select another source."
          stop 1
#endif
        endif
      elseif(nargs.ge.1) then
          ! If an argument is given, first test for the '-h' indicating a help
          ! request.
        call get_command_argument(1, linebuffer130, length=arglen, status=iostatus)
        testkey  = linebuffer130(1:1)
        testkey2 = linebuffer130(2:2)
        if(testkey.eq.'-')then
          if(testkey2.eq.'h')then
            ! This is the branch for user-requested help
            if(nargs.eq.1) then
              ! command is Ash3d -h
              call help_general
            else
              ! command is Ash3d -h [help topic]
              call get_command_argument(2, linebuffer130, length=arglen, status=iostatus)
              if(trim(adjustl(linebuffer130)).eq.'make')then
                call help_make
                write(outlog(io),*) ' --------------------------------------------'
                write(outlog(io),*) ' Here are the build parameters and run state:'
                write(outlog(io),*) ' --------------------------------------------'
                call Set_OS_Env
                stop 1
              elseif(trim(adjustl(linebuffer130)).eq.'run')then
                call help_run
              elseif(trim(adjustl(linebuffer130)).eq.'input')then
                if(nargs.ge.3) then
                  ! So far, we have Ash3d -h input
                  ! Check if there is an additional command-line parameter specifying
                  ! the block number

                  call get_command_argument(3, linebuffer130, length=arglen, status=iostatus)

                  read(linebuffer130,*,err=1600,iostat=iostatus,iomsg=iomessage)blockID
                  if(blockID.lt.1.or.blockID.gt.10)then
                    write(outlog(io),*) 'Input file block IS out of range 1-10'
                    stop 1
                  else
                    call help_inputfile(blockID)
                  endif
                else
                  call help_input
                endif
              elseif(trim(adjustl(linebuffer130)).eq.'postproc')then
                call help_postproc
              elseif(trim(adjustl(linebuffer130)).eq.'info')then
                call Set_OS_Env
                stop 1
              else
                write(outlog(io),*) 'Unknown help option'
                call help_general
              endif
            endif
          else
            !  This branch is reserved for when we can't figure out
            !  what help the user wants
            write(outlog(io),*) 'Unknown command-line options'
            call help_general
          endif
          stop 1
        else
          ! If the first argument does not begin with '-', then
          ! assume it is the input file name
          read(linebuffer130,*,iostat=iostatus,iomsg=iomessage)infile
          if(iostatus.lt.0)then
            write(errlog(io),*)'ERROR:  EOR encountered'
            stop 1
          elseif(iostatus.gt.0)then
            write(errlog(io),*)'ERROR:  Error reading control file name'
            write(errlog(io),*)'        From the command-line argument'
            write(errlog(io),*)linebuffer130
            write(errlog(io),*)'System Message: ',trim(adjustl(iomessage))
            stop 1
          endif
        endif
      elseif(nargs.lt.0) then
        ! When code called from ForestClaw, nargs is -1
        fc_len = 0
        do
          if(fc_inputfile(fc_len+1).eq.C_NULL_CHAR) exit
          fc_len = fc_len + 1
          infile(fc_len:fc_len) = fc_inputfile(fc_len)
        end do
        write(outlog(io),*) 'Reading input file ''',&
                           infile,''' from ForestClaw'
      endif

      return

1600  write(errlog(io),*)'ERROR: Unknown third command-line argument'
      write(errlog(io),*)'System Message: ',trim(adjustl(iomessage))
      stop 1

      end subroutine Parse_Command_Line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Subroutine Set_OS_Env 
!
!  This subroutine evaluates the state of all aspect of the run.
!  It is called from Ash3d.F90, or the help subroutine help_make
!  First, runtime environment variables are checked, including:
!    ASH3DVERB : (1-10), allows the verbosity level from the default (3) at runtime
!      = 1 : debug2        Additional debugging information only written to stdout
!      = 2 : debug1        Debugging information only written to stdout
!      = 3 : log           Time step information (this is the limit for writing to logfile)
!      = 4 : info          Additional information on run set up and shutdown
!      = 5 : statistics    Details on health of run (timing, mass conservation)
!      = 6 : production    Major program flow info
!      = 7 : essential     Only start up and shutdown messages
!      = 8 : error         No logging to stdout, only stderr (and logfile)
!      = 9 : silent        No logging to stdout,stderr. Logfile written as normal
!      =10 : dark          No logging to stdout,stderr or logfile
!    ASH3DHOME : The default installation directory is /opt/USGS/Ash3d, but may have been
!                customized in the makefile when compiled.  It can also be changed at runtime
!                via this environment variable.  This effects the location of data files Ash3d
!                might read at runtime.
!    ASH3DCFL  : (0.0<CFL<=1.0) The default CFL factor is 0.8, but this can be overwritten by
!                this environment variable.  Note; this could be subsequently overwritten
!                if there is a RESETPARAMS block in the input file
!    ASH3DPLOT : This environment variable allows overriding the default graphics package.
!      = 1 : DISLIN
!      = 2 : PLPLOT
!      = 3 : GNUPLOT
!      = 4 : GMT
!      = 5 : matlab/octave
!  Next, details of the system state are logged, including OS type (Linux, Mac, Windows),
!  endian flavor of hardware, fortran compiler version and flags, command-line arguments,
!  and date/time of the run.  Additionally, if PII=ON was set in the makefile when this
!  executable was compiled, then the user, system hostname and run directory are also logged.
!
!  Finally, all pre-processor flags are checked here with logging to stdout of which flags
!  are invoked. Pre-processor flags checked:
!   LINUX, MACOS, WINDOWS     : OS declaration
!   GFORTRAN,IFORT,AOCC,NVHPC : fortran compiler (for turning on/off non-standard subroutines)
!   FAST_DT, FAST_SUBGRID     : Speed-up tools for dt and grid calculations
!   EXPLDIFF, CRANKNIC        : Explicit vs implicit diffusion algorithm
!   LIM_NONE,LIM_LAXWEN,LIM_BW,LIM_FROMM,LIM_MINMOD,LIM_SUPERBEE,LIM_MC  : Limiter for advection
!   USENETCDF, USEGRIB        : Invokes netcdf and/or grib functionality
!   USEPOINTERS               : determines if variables are pointers or allocatable arrays
!   USEEXTDATA                : determines if Ash3d will rely on external lists for airport and
!                               volcano data
!   USEZIP                    : specifies that 'zip' should be used to compress/bundle kml files and pngs
!   USEDISLIN,USEPLPLOT,USEGNUPLOT,USEGMT : plotting packages linked or declared as available
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Set_OS_Env

      use global_param,  only : &
        DirPrefix,DirDelim,IsLitEnd,IsLinux,IsWindows,IsMacOS, &
        version,version_major,version_minor,version_patch,&
        CFL,OS_TYPE,OS_Flavor,os_full_command_line,os_cwd,os_host,os_user,&
        Comp_Code,Comp_Flavor,useFastDt,FastDt_suppress, &
        usezip,zippath,usegnuplot,gnuplotpath

      use io_data,       only : &
        Ash3dHome,cdf_source

      use time_data,     only : &
        BaseYear,useLeap,os_time_log, &
        RunStartDay,RunStartHour_ch,RunStartHr,RunStartMinute,RunStartMonth,RunStartYear

      use Output_Vars,   only : &
         iplotpref

      use MetReader,     only : &
         MR_OS_TYPE,MR_DirPrefix,MR_DirDelim,MR_VERB,MR_nio

      ! This module requires Fortran 2003 or later
      use iso_fortran_env, only : &
           compiler_version,&
           compiler_options

      integer            :: iostatus
      integer            :: cstat
      character(len=120) :: iomessage

      character(len=130) :: tmp_str
        ! variables to hold results of date_and_time
      character(len=8)  :: date
      character(len=10) :: time2
      character(len=5)  :: zone
      integer           :: values(8)
      integer           :: timezone
      real(kind=dp)     :: StartHour
      real(kind=dp)     :: RunStartHour    ! Start time of model run, in hours since BaseYear
      character(len=100):: CompVer
      character(len=608):: CompOpt
      logical           :: IsThere

      INTERFACE
        real(kind=8) function HS_hours_since_baseyear(iyear,imonth,iday,hours,byear,useLeaps)
          integer     ,intent(in) :: iyear
          integer     ,intent(in) :: imonth
          integer     ,intent(in) :: iday
          real(kind=8),intent(in) :: hours
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_hours_since_baseyear
        character (len=13) function HS_yyyymmddhhmm_since(HoursSince,byear,useLeaps)
          real(kind=8),intent(in) ::  HoursSince
          integer     ,intent(in) ::  byear
          logical     ,intent(in) ::  useLeaps
        end function HS_yyyymmddhhmm_since
        character (len=20) function HS_xmltime(HoursSince,byear,useLeaps)
          real(kind=8),intent(in) :: HoursSince
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_xmltime
      END INTERFACE

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Set_OS_Env"
      endif;enddo

      ! Reset OS varibles based on PP flags set in the makefile
#ifdef LINUX
      OS_TYPE = 1
      OS_Flavor = 'LINUX'
      DirPrefix = ''
      DirDelim  = '/'
      IsLinux   = .true.
      IsWindows = .false.
      IsMacOS   = .false.
#endif
#ifdef MACOS
      OS_TYPE = 2
      OS_Flavor = 'MACOS'
      DirPrefix = ''
      DirDelim = '/'
      IsLinux   = .false.
      IsWindows = .false.
      IsMacOS   = .true.
#endif
#ifdef WINDOWS
      OS_TYPE = 3
      OS_Flavor = 'WINDOWS'
      DirPrefix = 'C:'
      DirDelim  = '\'
      IsLinux   = .false.
      IsWindows = .true.
      IsMacOS   = .false.
#endif

#ifdef USEZIP
      !zippath = '/usr/bin/zip'
#include "zippath.h"
#endif

!#ifdef USEDISLIN
!  For now, no Dislin-specific varibles are set
!#endif
!#ifdef USEPLPLOT
!  For now, no plplot-specific varibles are set
!#endif
#ifdef USEGNUPLOT
      !gnuplotpath = '/usr/bin/gnuplot'
#include "gnuplotpath.h"
#endif
!#ifdef USEGMT
!  For now, nor GMT-specific varibles are set
!#endif

#ifdef GFORTRAN
      Comp_Code   = 1
      Comp_Flavor = 'gfortran'
#endif
#ifdef IFORT
      Comp_Code   = 2
      Comp_Flavor = 'ifort'
#endif
#ifdef AOCC
      Comp_Code   = 3
      Comp_Flavor = 'flang'
#endif
#ifdef NVHPC
      Comp_Code   = 4
      Comp_Flavor = 'nvfortran'
#endif

#ifdef FAST_DT
      ! With Fast_Dt on, time step criteria is only checked at the meso steps
      ! We can always suppress this if we have an active process that could
      ! effect the time step (e.g. umbrella spreading or scrubbing)
      ! Set the corresponding logical variables
      useFastDt       = .true.
      FastDt_suppress = .false.
#else
      useFastDt       = .false.
      FastDt_suppress = .true.
#endif

      ! Environment variables
      ! First order of business is to get the environment variable (if present) for verbosity
      call get_environment_variable(name="ASH3DVERB",value=tmp_str,status=iostatus)
      if(iostatus.eq.0)then
        ! read the value of the ASH3DVERB environment variable to the local variable VB(1)
        read(tmp_str,*,iostat=iostatus,iomsg=iomessage)VB(1)
        if(iostatus.ne.0)then
          write(outlog(1),*)"WARNING: ASH3DVERB found, but expecting an integer value"
          write(outlog(1),*)"         Instead, env. variable set to: ",tmp_str
          write(outlog(1),*)"         Resetting to ASH3DVERB=3"
          VB(1) = 3
        endif
        if(VB(1).le.verbosity_info)then
          write(outlog(1),*)"Checking for run-time environment variable: ASH3DVERB"
          write(outlog(1),*)"  Verbosity reset by environment variable"
        endif
      else
        if(VB(1).le.verbosity_info)then
          write(outlog(1),*)"Checking for run-time environment variable: ASH3DVERB"
          write(outlog(1),*)"  ASH3DVERB environment variable not found."
        endif
      endif
      if(VB(1).eq.1)then
        vlevel = "debug2"
      elseif(VB(1).eq.2)then
        vlevel = "debug1"
      elseif(VB(1).eq.3)then
        vlevel = "log"
      elseif(VB(1).eq.4)then
        vlevel = "info"
      elseif(VB(1).eq.5)then
        vlevel = "statistics"
      elseif(VB(1).eq.6)then
        vlevel = "production"
      elseif(VB(1).eq.7)then
        vlevel = "essential"
      elseif(VB(1).eq.8)then
        vlevel = "error"
      elseif(VB(1).eq.9)then
        vlevel = "silent"
      elseif(VB(1).eq.10)then
        vlevel = "dark"
      else
        write(errlog(1),*)"WARNING: Verbosity level not recognized. Value should be between 1 and 10."
        write(errlog(1),*)"         verbosity level : ",VB(1)
        write(outlog(1),*)"         Resetting to ASH3DVERB=3"
        VB(1) = 3
      endif
      if(VB(1).lt.9)then
        write(outlog(1),*)"    verbosity level : ",VB(1), vlevel
      endif

      if(VB(1).eq.verbosity_dark)then
        ! If the output verbosity is for nothing at all, reset the log verbosity to that too
        VB(2) = verbosity_dark
      else
         ! Before we do anything else, start a log file
        open(unit=fid_logfile,file=logfile,status='replace',action='write')
        ! Mirror the above output to the logfile
        if(iostatus.eq.0)then
          if(VB(2).le.verbosity_info)then
            write(outlog(2),*)" "
            write(outlog(2),*)"Checking for run-time environment variable: ASH3DVERB"
            write(outlog(2),*)"  Verbosity reset by environment variable to: ",VB(1)
          endif
        else
          if(VB(2).le.verbosity_info)then
            write(outlog(2),*)" "
            write(outlog(2),*)"Checking for run-time environment variable: ASH3DVERB"
            write(outlog(2),*)"  ASH3DVERB environment variable not found."
            write(outlog(2),*)"    verbosity level : ",VB(1)
          endif
        endif
      endif
      ! Harmonizing verbosity levels with MetReader
      MR_VERB = VB(1)
      MR_nio  = 2    ! Ash3d uses a logfile so set the output streams to 2 for stdin/stderr + logfile

      ! Next, check for environment variables ASH3DHOME
      ! Set the default installation path
      ! This is only needed if shared data files with fixed paths are read
      ! in such as the global airport and volcano ESP files.
      Ash3dHome = trim(adjustl(DirPrefix)) // DirDelim // &
                  "opt" // DirDelim // "USGS" // DirDelim // "Ash3d"

      ! Here it is over-written by compile-time path, if available
#ifndef WINDOWS
#include "installpath.h"
#endif
      ! This can be over-written if an environment variable is set
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(2),*)" "
        write(outlog(io),*)"Checking for run-time environment variable: ASH3DHOME"
      endif;enddo

      call get_environment_variable(name="ASH3DHOME",value=tmp_str,status=iostatus)
      if(iostatus.eq.0)then
        ! Environment variable ASH3DHOME found, now reading it
        Ash3dHome = tmp_str
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"  Install path reset by environment variable to: ",trim(adjustl(Ash3dHome))
        endif;enddo
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"  ASH3DHOME environment variable not found."
          write(outlog(io),*)"   Install path set at compilation is: ",trim(adjustl(Ash3dHome))
        endif;enddo
      endif

      ! Testing for the presence of a directory is compiler-specific
#ifdef IFORT
      ! With ifort, we need to test for a directory as follows
      inquire(directory=trim(adjustl(Ash3dHome)),exist=IsThere)
#else
      ! For testing the existance of a directory with gfortran or aocc, append a delimiter
      ! and . to make a file
      tmp_str = trim(adjustl(Ash3dHome)) // DirDelim // '.'
      inquire(file=trim(adjustl(tmp_str)),exist=IsThere)
#endif

      if(IsThere)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"  Path to ASH3DHOME found on system.  Good."
        endif;enddo
      else
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"WARNING: Cannot find ASH3DHOME=",trim(adjustl(Ash3dHome))
          write(errlog(io),*)"         If this run requires shared volcano or airport files,"
          write(errlog(io),*)"         it will fail. This directory could be missing if you have"
          write(errlog(io),*)"         not yet executed 'make install', or it could be that 'make install'"
          write(errlog(io),*)"         failed if permissions were not adaquate."
          write(errlog(io),*)"         Please either recompile with the install directory set,"
          write(errlog(io),*)"         or set the environment variable ASH3DHOME"
        endif;enddo
      endif

      ! Now, check for environment variables ASH3DCFL
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)" "
        write(outlog(io),*)"Checking for run-time environment variable: ASH3DCFL"
      endif;enddo
      call get_environment_variable(name="ASH3DCFL",value=tmp_str,status=iostatus)
      if(iostatus.eq.0)then
        ! Environment variable ASH3DCFL found, now reading it
        read(tmp_str,*,iostat=iostatus,iomsg=iomessage)CFL
        if(iostatus.ne.0)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: ASH3DCFL found, but expecting a floating point value"
            write(errlog(io),*)"       Instead, env. variable set to: ",tmp_str
            write(errlog(io),*)'System Message: ',trim(adjustl(iomessage))
          endif;enddo
          stop 1
        endif
        if(CFL.le.0.0_ip)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),'(a22)')     "ERROR: CFL must be > 0"
            write(errlog(io),'(a24,f8.2)')"       Currently set to ",CFL
          endif;enddo
          stop 1
        elseif(CFL.ge.1.0_ip)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),'(a22)')     "ERROR: CFL must be < 1"
            write(errlog(io),'(a24,f8.2)')"       Currently set to ",CFL
          endif;enddo
          stop 1
        endif
        ! CFL variable seems valid, proceeding.
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),'(a50,f8.2)')"  CFL condition reset by environment variable to: ",CFL
          write(outlog(io),*)"   Note: It is possible this may be subsequently reset via the"
          write(outlog(io),*)"         control file in a RESETPARAM block. Check the log file"
          write(outlog(io),*)"         or the netcdf output file for final CFL used."
        endif;enddo
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),'(a40,f8.2)')"  ASH3DCFL not found.  CFL condition : ",CFL
        endif;enddo
      endif

      ! Now, check for environment variables ASH3DPLOT
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)" "
        write(outlog(io),*)"Checking for run-time environment variable: ASH3DPLOT"
      endif;enddo
      call get_environment_variable(name="ASH3DPLOT",value=tmp_str,status=iostatus)
      if(iostatus.eq.0)then
        ! Environment variable ASH3DPLOT found, now reading it
        read(tmp_str,*,iostat=iostatus,iomsg=iomessage)iplotpref
        if(iostatus.ne.0)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: ASH3DPLOT found, but expecting an integer value"
            write(errlog(io),*)"       Instead, env. variable set to: ",tmp_str
            write(errlog(io),*)'System Message: ',trim(adjustl(iomessage))
          endif;enddo
          stop 1
        endif
        if(iplotpref.le.0.or.iplotpref.gt.5)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: ASH3DPLOT must be any of:"
            write(errlog(io),*)"         1 = dislin"
            write(errlog(io),*)"         2 = plplot"
            write(errlog(io),*)"         3 = gnuplot"
            write(errlog(io),*)"         4 = GMT"
            write(errlog(io),*)"         5 = matlab/octave"
            ! Placeholders for other post-processing graphics packages
            !write(errlog(io),*)"         6 = cartopy"
            !write(errlog(io),*)"         7 = R"
            write(errlog(io),*)"       Currently set to ",iplotpref
          endif;enddo
          stop 1
        else
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"  Plotting library reset from environment variable ASH3DPLOT=",iplotpref
            if(iplotpref.eq.1)then
              write(outlog(io),*)"           Using dislin if available."
            elseif(iplotpref.eq.2)then
              write(outlog(io),*)"           Using plplot if available."
            elseif(iplotpref.eq.3)then
              write(outlog(io),*)"           Using gnuplot if available."
            elseif(iplotpref.eq.4)then
              write(outlog(io),*)"           Using GMT if available."
            elseif(iplotpref.eq.5)then
              write(outlog(io),*)"           Using matlab/octave if available."
            endif
          endif;enddo
        endif
      endif

      ! Operating System details
      MR_OS_TYPE   = OS_TYPE
      MR_DirPrefix = DirPrefix
      MR_DirDelim  = DirDelim
      ! Find out if we are running on a little-endian or big-endian system
      call check_endian(IsLitEnd)

      ! Get info on the compiler version and options
      CompVer = compiler_version()
      CompOpt = compiler_options()

      ! Get some run-specific and system-specific information
      call get_command(os_full_command_line)

      ! Fill these with N/A in case we don't have the fuctions to fill them
      os_user = 'N/A'
      os_host = 'N/A'
      os_cwd  = 'N/A'
#ifdef USEPII
      call get_environment_variable(name="USER",value=tmp_str,status=iostatus)
      os_user = adjustl(trim(tmp_str))
      call get_environment_variable(name="HOSTNAME",value=tmp_str,status=iostatus)
      os_host = adjustl(trim(tmp_str))
      call get_environment_variable(name="PWD",value=tmp_str,status=iostatus)
      os_cwd = adjustl(trim(tmp_str))
#endif
      call date_and_time(date,time2,zone,values)
        ! date  = ccyymmdd
        ! time2 = hhmmss.sss
        ! zone  = (+-)hhmm, representing the difference with respect to Coordinated Universal Time (UTC)
        ! VALUE(1) = The year
        ! VALUE(2) = The month
        ! VALUE(3) = The day of the month
        ! VALUE(4) = Time difference with UTC in minutes
        ! VALUE(5) = The hour of the day
        ! VALUE(6) = The minutes of the hour
        ! VALUE(7) = The seconds of the minute
        ! VALUE(8) = The milliseconds of the second

      do io=1,2;if(VB(io).le.verbosity_essential)then
        write(outlog(io),*)""
        write(outlog(io),*)"Running Ash3d with command line: ",&
                    trim(adjustl(os_full_command_line))
      endif;enddo

      ! Determining the run start time
      !  We skip error-handling here since we are reading from the string 'zone' from date_and_time
      read(zone,'(i3)',iostat=iostatus,iomsg=iomessage) timezone
      ! Find time in UTC
      StartHour = real(values(5)-timezone,kind=ip) + real(values(6)/60.0,kind=ip)    ! add offset to UTC
        ! find time in HoursSinceBaseYear
        !  Note: This will be relative to the BaseYear in time_data (default is 1900). 
        !        That BaseYear might be changed if the eruption start time is
        !        before BaseYear
      RunStartHour    = HS_hours_since_baseyear(values(1),values(2),values(3), &
                                                StartHour,BaseYear,useLeap)
      RunStartHour_ch = HS_yyyymmddhhmm_since(RunStartHour,BaseYear,useLeap)
      read(RunStartHour_ch,'(i4)',iostat=iostatus,iomsg=iomessage) RunStartYear
      read(RunStartHour_ch,'(4x,i2)',iostat=iostatus,iomsg=iomessage) RunStartMonth
      read(RunStartHour_ch,'(6x,i2)',iostat=iostatus,iomsg=iomessage) RunStartDay
      read(RunStartHour_ch,'(8x,i2)',iostat=iostatus,iomsg=iomessage) RunStartHr
      read(RunStartHour_ch,'(11x,i2)',iostat=iostatus,iomsg=iomessage) RunStartMinute

        ! Prepare a note to include in the netcdf output file
      os_time_log = HS_xmltime(RunStartHour,BaseYear,useLeap)

      do io=1,2;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)" System Information"
        if(IsLitEnd)then
          write(outlog(io),*)"   host: ",trim(adjustl(os_host)), &
                              ' (',trim(adjustl(OS_Flavor)),' little-endian)'
        else
          write(outlog(io),*)"   host: ",trim(adjustl(os_host)), &
                              ' (',trim(adjustl(OS_Flavor)),' little-endian)'
        endif
        write(outlog(io),*)"    cwd: ",trim(adjustl(os_cwd))
        write(outlog(io),*)"   user: ",trim(adjustl(os_user))

        write(outlog(io),*)
        write(outlog(io),*)"This executable was compiled with the following compiler and options:"
        write(outlog(io),*)"    ",trim(adjustl(CompVer))
        write(outlog(io),*)"    ",trim(adjustl(CompOpt))
        write(outlog(io),*)"and with the following pre-processor flags:"
#if defined GFORTRAN
        write(outlog(io),*)"      GFORTRAN: Compiler specified in makefile"
#elif defined IFORT
        write(outlog(io),*)"         IFORT: Compiler specified in makefile"
#elif defined AOCC
        write(outlog(io),*)"          AOCC: Compiler specified in makefile"
#elif defined NVHPC
        write(outlog(io),*)"         NVHPC: Compiler specified in makefile"
#endif

#if defined LINUX
        write(outlog(io),*)"         LINUX: System specified as linux"
#elif defined MACOS
        write(outlog(io),*)"         MACOS: System specified as MacOS"
#elif defined WINDOWS
        write(outlog(io),*)"       WINDOWS: System specified as MS Windows"
#endif

#ifdef USEZIP
        if(IsWindows)then
          write(outlog(io),*)"      USEZIP: zip set to T, but this is a Windows system"
          write(outlog(io),*)"              zip currently not integrated with Ash3d on Windows."
          write(outlog(io),*)"       Deactivating zip"
          usezip = .false.
        else
          usezip = .true.
          write(outlog(io),*)"        USEZIP: zip will be used to bundle kmz files using path: ",&
                             trim(adjustl(zippath))
          ! First check if zip is actually in the path
          write(outlog(io),*)"                  Checking for default path for zip"
          call execute_command_line('which zip',wait=.true.,exitstat=iostatus)
          if(iostatus.ne.0)then
            write(outlog(io),*)"Error: 'which zip' failed. No zip executable in default path"
            write(outlog(io),*)"       Deactivating zip"
            usezip = .false.
          endif
          if(usezip)then
            ! Next check if zippath exists
            inquire( file=trim(adjustl(zippath)), exist=IsThere)
            if(.not.IsThere)then
              write(outlog(io),*)"Error: user-specified path does not exist and is"
              write(outlog(io),*)"       inconsistent with 'which zip'."
              write(outlog(io),*)"       Please correct zippath in makefile."
              write(outlog(io),*)"       Deactivating zip"
              usezip = .false.
            endif
          endif
          if(usezip)then
            ! Finally, do a test run of the zip executable
            write(outlog(io),*)"                  Checking if zip executes."
            call execute_command_line("echo 'exit' | zip --version > /dev/null",&
                                      WAIT=.true., EXITSTAT=iostatus, CMDSTAT=cstat, CMDMSG=iomessage)
            if(iostatus.eq.0)then
              write(outlog(io),*)"                  Success"
            else
              write(outlog(io),*)"Error: Something is wrong with the zip executable."
              write(outlog(io),*)"         zip is returing an error code = ",iostatus
              write(outlog(io),*)"       execute_command_line command status = ",cstat
              write(outlog(io),*)"       execute_command_line error message: ",trim(adjustl(iomessage))
              write(outlog(io),*)"       Deactivating zip"
              usezip = .false.
            endif
          endif
        endif
#endif

#ifdef USEDISLIN
        write(outlog(io),*)"     USEDISLIN: API to Dislin plotting package is enabled"
#endif
#ifdef USEPLPLOT
        write(outlog(io),*)"     USEPLPLOT: API to plplot plotting package is enabled"
#endif
#ifdef USEGNUPLOT
        if(IsWindows)then
          write(outlog(io),*)"    USEGNUPLOT: gnuplot set to T, but this is a Windows system"
          write(outlog(io),*)"                gnuplot currently not integrated with Ash3d on Windows."
          write(outlog(io),*)"       Deactivating gnuplot"
          usegnuplot = .false.
        else
          usegnuplot = .true.
          write(outlog(io),*)"    USEGNUPLOT: gnuplot plotting package is installed"
          ! First check if gnuplot is actually in the path
          write(outlog(io),*)"                  Checking for default path for gnuplot"
          call execute_command_line('which gnuplot',wait=.true.,exitstat=iostatus)
          if(iostatus.ne.0)then
            write(outlog(io),*)"Error: 'which gnuplot' failed. No gnuplot executable in default path"
            write(outlog(io),*)"       Deactivating gnuplot"
            usegnuplot = .false.
          endif
          if(usegnuplot)then
            ! Next check if gnuplotpath exists
            inquire( file=trim(adjustl(gnuplotpath)), exist=IsThere)
            if(.not.IsThere)then
              write(outlog(io),*)"Error: user-specified path does not exist and is"
              write(outlog(io),*)"       inconsistent with 'which gnuplot'."
              write(outlog(io),*)"       Please correct gnuplotpath in makefile."
              write(outlog(io),*)"       Deactivating gnuplot"
              usegnuplot = .false.
            endif
          endif
          if(usegnuplot)then
            ! Finally, do a test run of the gnuplot executable
            write(outlog(io),*)"                  Checking if gnuplot executes."
            call execute_command_line("echo 'exit' | gnuplot",&
                                      WAIT=.true., EXITSTAT=iostatus, CMDSTAT=cstat, CMDMSG=iomessage)
            if(iostatus.eq.0)then
              write(outlog(io),*)"                  Success"
            else
              write(outlog(io),*)"Error: Something is wrong with the gnuplot executable."
              write(outlog(io),*)"       gnuplot is returing an error code",iostatus
              write(outlog(io),*)"       execute_command_line command status = ",cstat
              write(outlog(io),*)"       execute_command_line error message: ",trim(adjustl(iomessage))
              write(outlog(io),*)"       Deactivating gnuplot"
              usegnuplot = .false.
            endif
          endif
        endif
#endif
#ifdef USEGMT
        write(outlog(io),*)"        USEGMT: API to GMT plotting package is enabled"
#endif

#ifdef FAST_DT
        write(outlog(io),*)"       FAST_DT: ON"
        write(outlog(io),*)"                dt will only be evaluated on the time steps"
        write(outlog(io),*)"                in the wind files.  If there are processes"
        write(outlog(io),*)"                that affect the wind speeds (e.g. umbrella"
        write(outlog(io),*)"                spreading), this can cause job failure."
#else
        write(outlog(io),*)"       FAST_DT: OFF"
        write(outlog(io),*)"                dt will be evaluated each time step"
#endif
#ifdef FAST_SUBGRID
        write(outlog(io),*)"  FAST_SUBGRID: ON"
        write(outlog(io),*)"                Advection and diffusion routines will only"
        write(outlog(io),*)"                be calculated in the region where the cloud"
        write(outlog(io),*)"                concentration exceeds a threshold"
#else
        write(outlog(io),*)"  FAST_SUBGRID: OFF"
        write(outlog(io),*)"                Advection and diffusion routines will be"
        write(outlog(io),*)"                applied on all cells, including those with"
        write(outlog(io),*)"                negligible concentrations."
#endif
#ifdef EXPLDIFF
        write(outlog(io),*)"      EXPLDIFF: Diffusion will be calculated via the explicit solver."
#endif
#ifdef CRANKNIC
        write(outlog(io),*)"      CRANKNIC: Diffusion will be calculated implicitly (via Crank-Nicolson)"
#endif
#ifdef LIM_NONE
        write(outlog(io),*)"      LIM_NONE: Advection routines use no limiters"
#endif
#ifdef LIM_LAXWEN
        write(outlog(io),*)"    LIM_LAXWEN: Advection routines use a Lax-Wendroff limiter"
#endif
#ifdef LIM_BW
        write(outlog(io),*)"        LIM_BW: Advection routines use a Beam-Warming limiter"
#endif
#ifdef LIM_FROMM
        write(outlog(io),*)"     LIM_FROMM: Advection routines use a Fromm limiter"
#endif
#ifdef LIM_MINMOD
        write(outlog(io),*)"    LIM_MINMOD: Advection routines use a minmod limiter"
#endif
#ifdef LIM_SUPERBEE
        write(outlog(io),*)"  LIM_SUPERBEE: Advection routines use a superbee limiter"
#endif
#ifdef LIM_MC
        write(outlog(io),*)"        LIM_MC: Advection routines use a MC limiter"
#endif
#ifdef USENETCDF
        write(outlog(io),*)"     USENETCDF: ON"
        write(outlog(io),*)"                NetCDF functionality is included"
#else
        write(outlog(io),*)"     USENETCDF: OFF"
        write(outlog(io),*)"                NetCDF functionality is not included"
#endif
#ifdef USEGRIB
        write(outlog(io),*)"       USEGRIB: ON"
        write(outlog(io),*)"                Grib functionality is included"
#else
        write(outlog(io),*)"       USEGRIB: OFF"
        write(outlog(io),*)"                Grib functionality is not included"
#endif
#ifdef USEPOINTERS
        write(outlog(io),*)"   USEPOINTERS: ON"
        write(outlog(io),*)"                Arrays are defined as pointers"
        write(outlog(io),*)"                This helps Ash3d subroutines to be called via C++"
#else
        write(outlog(io),*)"   USEPOINTERS: OFF"
        write(outlog(io),*)"                All arrays are allocatable."
#endif
#ifdef USEEXTDATA
        write(outlog(io),*)"    USEEXTDATA: ON"
        write(outlog(io),*)"                Data files for airports and volcanoes are"
        write(outlog(io),*)"                read at run-time"
#else
        write(outlog(io),*)"    USEEXTDATA: OFF"
        write(outlog(io),*)"                Data arrays for airports and volcanoes are"
        write(outlog(io),*)"                included when compiled."
#endif

        ! Write out start time in UTC
        write(outlog(io),*)
        write(version,'(i0,a1,i0,a1,i0)')version_major,'.',version_minor,'.',version_patch
        write(outlog(io),2) version,RunStartYear,RunStartMonth,RunStartDay,RunStartHr,RunStartMinute
        write(outlog(io),*)
      endif;enddo

      write(cdf_source,'(a7,i0,a1,i0,a1,i0)')'Ash3d v',version_major,'.',version_minor,'.',version_patch
      ! Format statements
2     format(4x,'Ash3d ( v. ',a8,') run ',&
             i4,'.',i2.2,'.',i2.2,i4,':',i2.2,' UTC')

      end subroutine Set_OS_Env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Subroutine CHECK_ENDIAN checks if the local system uses Big-Endian
!  or Little-Endian byte ordering.  Returns the logical value
!  IsLitEnd = .true. if the system is Little-Endian, .false. otherwise.
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine check_endian(IsLitEnd)

      logical,intent(inout)  :: IsLitEnd

      integer(kind=2)  :: s = 1

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine check_endian"
      endif;enddo

      if(btest(transfer(int((/1,0/),kind=1),s),0)) then
        ! System is Little-Endian
        IsLitEnd = .true.
      else
        ! System is Big-Endian
        IsLitEnd = .false.
      endif

      end subroutine check_endian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Subroutine Read_Control_File()
!
! This subroutine sets up the parameters for the Ash3d run.
! 
! The control file is opened and read block-by-block.
!       ! BLOCK 1: GRID INFO
!       ! BLOCK 2: ERUPTION PARAMETERS
!       ! BLOCK 3: WIND PARAMETERS
!   Sets source term variables
!       ! BLOCK 4: OUTPUT OPTIONS
!       ! BLOCK 5: INPUT WIND FILES
!   Sets variable in MetReader
!       ! BLOCK 6: AIRPORT FILE
!       ! BLOCK 7: GRAIN-SIZE BINS, SETTLING VELOCITY
!   Sets parameters in Tephra module
!       ! BLOCK 8: VERTICAL PROFILES
!       ! BLOCK 9 (Optional): NETCDF ANNOTATIONS
! After these standard block of the input file are read, the remaining part of the
! file is searched for blocks that are part of user-built optional modules.  These
! are indentified by the keywork OPTMOD and are read by customized subroutines in
! the optional module.
!
! Lastly, some notes are written to stdout and the logfile specifying some aspects
! of the run.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !subroutine Read_Control_File(fc_inputfile)
      subroutine Read_Control_File

      ! This module requires Fortran 2003 or later
      use iso_c_binding

      ! This module requires Fortran 2003 or later
      use iso_fortran_env, only : &
         input_unit

      use global_param,  only : &
         EPS_SMALL,EPS_TINY,nmods,OPTMOD_names,limiter,&
         useDS,useTemperature,useCalcFallVel,useLogNormGSbins,&
         useDiffusion,useCN,useVz_rhoG,M2PS_2_KM2PHR,MAXNUM_OPTMODS

      use io_data,       only : &
         cdf_b1l1,cdf_b1l2,cdf_b1l3,cdf_b1l4,cdf_b1l5,cdf_b1l6,cdf_b1l7,cdf_vardz,cdf_b1l8,cdf_b1l9,&
         cdf_b3l1,cdf_b3l2,cdf_b3l3,cdf_b3l4,cdf_b3l5,cdf_b4l1,cdf_b4l2,cdf_b4l3,cdf_b4l4,&
         cdf_b4l5,cdf_b4l6,cdf_b4l7,cdf_b4l8,cdf_b4l9,cdf_b4l10,cdf_b4l11,cdf_b4l12,cdf_b4l13,&
         cdf_b4l14,cdf_b4l15,cdf_b4l16,cdf_b4l17,cdf_b4l18,cdf_b5l1,cdf_b6l1,cdf_b6l2,cdf_b6l3,cdf_b6l4,&
         cdf_b6l5,Have_Block_NetCDF,cdf_comment,cdf_title,cdf_institution,cdf_source_url,cdf_history,cdf_references,&
         concenfile,VolcanoName,WriteTimes,nWriteTimes,cdf_conventions,cdf_run_class,cdf_url,&
         x_vprofile,y_vprofile,i_vprofile,j_vprofile,Site_vprofile,&
         infile,ioutputFormat,LoadConcen,log_step,NextWriteTime,Have_Block_ResParm,&
         AppendExtAirportFile,WriteInterval,WriteGSD,WriteDepositTS_KML,WriteDepositTS_ASCII,&
         WriteDepositTime_KML,WriteDepositTime_ASCII,WriteDepositFinal_KML,&
         WriteDepositFinal_ASCII,WriteCloudTime_KML,WriteCloudTime_ASCII,&
         WriteCloudLoad_KML,WriteReflectivity_KML,WriteCloudLoad_ASCII,WriteCloudHeight_KML,&
         WriteCloudHeight_ASCII,WriteCloudConcentration_KML,WriteCloudConcentration_ASCII,&
         WriteAirportFile_KML,WriteAirportFile_ASCII,Write3dFiles,ReadExtAirportFile,&
         Output_every_TS,Output_at_WriteTimes,Output_at_logsteps,nvprofiles,iTimeNext,&
         Write_PT_Data,Write_PR_Data,isFinal_TS

      use mesh,          only : &
         de,dn,dx,dy,z_vec_init,dz_const,nxmax,nymax,nzmax,nsmax,VarDzType,ivent,jvent,kvent,&
         gridwidth_e,gridwidth_n,gridwidth_x,gridwidth_y,&
         lonLL,latLL,lonUR,latUR,xLL,yLL,xUR,yUR,&
         A3d_iprojflag,A3d_k0,A3d_phi0,A3d_lam0,A3d_lam1,A3d_phi1,A3d_lam2,&
         A3d_phi2,A3d_Re,IsLatLon,IsPeriodic,ZPADDING,Ztop

      use solution,      only : &
         StopWhenDeposited,StopValue_FracAshDep,StopValue_FracAshDep_Default,imin,imax,jmin,jmax,kmin,kmax

      use time_data,     only : &
         BaseYear,useLeap,time,SimStartHour,Simtime_in_hours,xmlSimStartTime

      use Source,        only : &
         neruptions,e_Duration,e_Volume,e_PlumeHeight,e_prof_Volume,e_prof_dz,&
         e_prof_maxpoints,e_prof_nzpoints,e_StartTime,MAX_ER_PROFPOINTS,&
         ESP_duration,ESP_height,ESP_Vol,&
         lat_volcano,lon_volcano,x_volcano,y_volcano,z_volcano,Suzuki_A,&
         IsCustom_SourceType,SourceType,SourceType_idx,&
           Allocate_Source_eruption

      use Source_Umbrella, only : &
         SuzK_umb

      use Tephra,        only : &
         DepositDensity,Tephra_Ncols,Tephra_v_s,Tephra_gsdiam,Tephra_bin_mass,Tephra_rho_m,Tephra_gsPhi,&
         Tephra_gsF,Tephra_gsG,FV_ID,Shape_ID,LN_phi_mean,LN_phi_stddev,LN_massfrac,n_gs_max,n_gs_aloft,&
           Calculate_Tephra_Shape,&
           Allocate_Tephra, &
           Sort_Tephra_Size

      use Output_Vars,   only : &
         useRestartVars

      use Airports,      only : &
         AirportInFile,ProjectAirportLocations

      use VotW_ESP,      only : &
           get_ESP

      use Diffusion,     only : &
         diffusivity_horz,diffusivity_vert,Imp_fac,Imp_DT_fac, &
           Allocate_Diff

      use projection,    only : &
         PJ_iprojflag,PJ_k0,PJ_lam0,PJ_lam1,PJ_lam2,PJ_phi0,PJ_phi1,PJ_phi2,PJ_Re,&
           PJ_Set_Proj_Params

      use MetReader,     only : &
         MR_iWindFiles,MR_WindFiles,MR_BaseYear,MR_useLeap,MR_Comp_StartHour,&
         MR_WindFiles_GRIB_index,MR_WindFiles_Have_GRIB_index,MR_Comp_Time_in_hours,&
         MR_WindFile_starthour,MR_windfile_stephour,MR_iHeightHandler,&
         MR_iwf_template,MR_iwind,MR_Comp_StartYear,MR_Comp_StartMonth,MR_ztop,&
           MR_Allocate_FullMetFileList, &
           MR_Read_Met_DimVars

#ifdef USENETCDF
      use Ash3d_Netcdf_IO
#endif

      INTERFACE
        subroutine help_input
        end subroutine help_input
      END INTERFACE

      integer           :: i,k,ii,isize

      integer,       allocatable, dimension(:) :: e_iyear  ! time data read from files
      integer,       allocatable, dimension(:) :: e_imonth
      integer,       allocatable, dimension(:) :: e_iday
      real(kind=dp), allocatable, dimension(:) :: e_hour   ! Start time of eruption in
                                                           !  hour (UT)
      character(len=50) :: linebuffer050 
      character(len=80) :: linebuffer080
      character(len=130):: linebuffer130
      character(len=400):: linebuffer400 ! Used for reading line lists of values (write times, etc)
      character(len=3)  :: answer
      character(len=6)  :: formatanswer
      character(len=20) :: mod_name
      character(len=20) :: dumstr20

      integer           :: iw,iwf,igrid,idf,iwfiles
      integer           :: ivalue1, ivalue2, ivalue3, ivalue4
      integer           :: loc
      integer           :: iform

      character(len=80) :: Comp_projection_line
      integer           :: ilatlonflag
      character         :: testkey
      integer           :: iostatus
      character(len=120):: iomessage
      integer           :: iendstr,init_n_gs_max
      real(kind=ip)     :: value1, value2, value3, value4, value5
      real(kind=ip),allocatable,dimension(:) :: values
      real(kind=ip)     :: tmp_ip
      real(kind=dp)     :: tmp_dp
      real(kind=ip)     :: sum_bins
      character(len=8)  :: volc_code
      real(kind=ip),allocatable,dimension(:) :: dum_prof
      integer           :: clog_nsteps,cust_nsteps
      real(kind=ip)     :: clog_zmax


      real(kind=ip),allocatable,dimension(:) :: temp_v_s,temp_gsdiam
      real(kind=ip),allocatable,dimension(:) :: temp_bin_mass,temp_rho_m
      real(kind=ip),allocatable,dimension(:) :: temp_gsF,temp_gsG,temp_phi
      real(kind=ip)     :: fracfine = 0.0_ip
      !real(kind=ip)     :: CompGrid_height
      real(kind=ip)     :: last_z
      integer           :: nz_init,nsegments
      integer      ,allocatable,dimension(:) :: nz_plin_segments
      real(kind=ip),allocatable,dimension(:) :: dz_plin_segments
      integer           :: substr_pos1
      integer           :: substr_pos2
      logical           :: IsThere
      logical           :: runAsForecast       = .false.  ! This will be changed if year=0
      real(kind=dp)     :: FC_Offset = 0.0_dp
      real(kind=ip)     :: Davg,Aaxis,Baxis,Caxis
      logical           :: IsOpen
      logical           :: IsSat
      logical           :: IsComment

      INTERFACE
        subroutine input_data_ResetParams
        end subroutine input_data_ResetParams
        real(kind=8) function HS_hours_since_baseyear(iyear,imonth,iday,hours,byear,useLeaps)
          integer     ,intent(in) :: iyear
          integer     ,intent(in) :: imonth
          integer     ,intent(in) :: iday
          real(kind=8),intent(in) :: hours
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_hours_since_baseyear
        character (len=13) function HS_yyyymmddhhmm_since(HoursSince,byear,useLeaps)
          real(kind=8),intent(in) ::  HoursSince
          integer     ,intent(in) ::  byear
          logical     ,intent(in) ::  useLeaps
        end function HS_yyyymmddhhmm_since
        character (len=20) function HS_xmltime(HoursSince,byear,useLeaps)
          real(kind=8),intent(in) :: HoursSince
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_xmltime
        subroutine MR_Set_Gen_Index_GRIB(grib_file)
          character(len=130),intent(in)  :: grib_file
        end subroutine MR_Set_Gen_Index_GRIB
        subroutine help_inputfile(blockID)
          integer,intent(in) :: blockID
        end subroutine help_inputfile
      END INTERFACE

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"--------------------------------------------------"
        write(outlog(io),*)"---------- READ_CONTROL_FILE ---------------------"
        write(outlog(io),*)"--------------------------------------------------"
      endif;enddo

      ! Initialize output
      formatanswer  = 'null'
      isFinal_TS    = .false.
      nWriteTimes   = 0                 ! number of output files to write (default=0)
      NextWriteTime = 1.0_ip/EPS_TINY   ! Time to write the next file (default = never)

      ! Open and read control file
      do io=1,2;if(VB(io).le.verbosity_production)then
        write(outlog(io),3) infile
      endif;enddo

      inquire( file=infile, exist=IsThere )
      if(.not.IsThere)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Cannot find input file"
        endif;enddo
        stop 1
      endif
      open(unit=fid_ctrlfile,file=infile,status='old',action='read',err=9001)

      !************************************************************************
      ! Searching for optional blocks labled by OPTMOD
      !  These are expected to be at the end of the control file, but we need to
      !  know now if we are calling input_data_ResetParams as this can effect the
      !  grid
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)' *****************************************'
        write(outlog(io),*)' Reading control file for optional modules   '
        write(outlog(io),*)' *****************************************'
      endif;enddo

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Searching for blocks with OPTMOD"
      endif;enddo
      nmods = 0
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      ! if there are no further blocks, then we will skip over this while loop
      do while(iostatus.eq.0)
        substr_pos1 = index(linebuffer080,'OPTMOD')
        if(substr_pos1.eq.1)then
          ! found an optional module
          nmods = nmods + 1
          if(nmods.gt.MAXNUM_OPTMODS)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: Maximum number of optional modules exceeded"
              write(errlog(io),*)"       Current maximum set to MAXNUM_OPTMODS = ",MAXNUM_OPTMODS
              write(errlog(io),*)"       Please increase MAXNUM_OPTMODS and recompile."
              write(errlog(io),*)"  Ash3d_VariableModules.f90:global_param:MAXNUM_OPTMODS"
            endif;enddo
            stop 1
          endif
          !  Parse for the keyword
          read(linebuffer080,1104,iostat=iostatus,iomsg=iomessage)mod_name
          linebuffer050 = "Reading control file blk9+ (OPTMOD)"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
          OPTMOD_names(nmods) = trim(adjustl(mod_name))
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"     Found optional module : ",&
                                OPTMOD_names(nmods),nmods
          endif;enddo
        endif
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.lt.0)then
          ! end of file reached; exit do loop
          exit
        elseif(iostatus.gt.0)then
          ! Some non-EOF error
          linebuffer050 = "Reading control file blk9+ (OPTMOD)"
          call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        endif
1104    format(7x,a20)
      enddo
      if(nmods.eq.0)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"No OPTMOD blocks found."
        endif;enddo
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Number of OPTMOD blocks found = ",nmods
        endif;enddo
      endif

      ! Now checking if we need to reset any parameters
      do i=1,nmods
        do io=1,2;if(VB(io).le.verbosity_essential)then
          write(outlog(io),*)"Testing for ",OPTMOD_names(i),i
        endif;enddo
        if(OPTMOD_names(i).eq.'RESETPARAMS')then
          do io=1,2;if(VB(io).le.verbosity_essential)then
            write(outlog(io),*)"  Reading input block for RESETPARAMS"
          endif;enddo
          Have_Block_ResParm = .true.
          call input_data_ResetParams
        endif
      enddo
      ! Reopen or rewind to begining so we can start parsing block 1
      !   check to make sure the control file is open
      inquire(unit=fid_ctrlfile,opened=IsOpen)
      if(IsOpen)then
        rewind(fid_ctrlfile)
      else
        open(unit=fid_ctrlfile,file=infile,status='old',action='read',err=9001)
      endif

      !************************************************************************
      ! BLOCK 1: GRID INFO
      ! Start reading the input file assuming there is a variable length
      ! header with each header line flagged by a '#' or '*' in the first position
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)' *******************************************'
        write(outlog(io),*)' Reading Block 1: Volcano/grid specification'
        write(outlog(io),*)' *******************************************'
      endif;enddo
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading first line of control file."
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,'(a1)',iostat=iostatus,iomsg=iomessage)testkey
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      do while (IsComment)
         ! Line is a comment, read next line
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        linebuffer050 = "Reading comment line of control file (until Blk1)"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
        if(iostatus.ne.0)then
          do io=1,nio;if(VB(io).le.verbosity_error)then
            if(iostatus.lt.0)then
              write(errlog(io),*)'ERROR:  EOR encountered'
            else
              write(errlog(io),*)'ERROR:  Error reading character from string'
              write(errlog(io),*)'           From the following line from the file: '
              write(errlog(io),*)linebuffer080
              write(errlog(io),*)'System Message: ',trim(adjustl(iomessage))
            endif
          endif;enddo
          stop 1
        else
          call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
        endif
      enddo

      ! Block 1 Line 1
      ! Read volcano name
      cdf_b1l1 = linebuffer080
      iendstr = scan(linebuffer080, "#")
      if(iendstr.eq.1)then
        ! End-of-string marker is in spot 1
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
           "Volcano name cannot start with #"
        endif;enddo
        stop 1
      elseif(iendstr.eq.0)then
        ! End-of-string marker is not found
        iendstr = scan(linebuffer080, " ",.true.) ! rescan for space with back=.true.
      endif
      iendstr = min(iendstr,30)  ! limit string to 30 characters
      VolcanoName = trim(adjustl(linebuffer080(1:iendstr-1)))
      ! Check if the volcano name is a text name or a Smithsonian
      ! database ID
      read(VolcanoName,*,iostat=iostatus,iomsg=iomessage)testkey
      if(testkey.eq.'0'.or.testkey.eq.'1')then
        ! the 'name' is the CAVW Smithsonian ID
        ! get the source parameters for this volcano
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Calling get_ESP"
        endif;enddo
        volc_code = VolcanoName(1:8)
        call get_ESP(volc_code)
      endif

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)
        write(outlog(io),37) VolcanoName
      endif;enddo

      ! Block 1 Line 2
      ! Read projection parameters
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage) linebuffer080
      linebuffer050 = "Reading control file, Block 1, Line 2."
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      cdf_b1l2 = linebuffer080
      Comp_projection_line = linebuffer080
      read(Comp_projection_line,*,iostat=iostatus,iomsg=iomessage)ilatlonflag
      if(iostatus.ne.0)then
        do io=1,nio;if(VB(io).le.verbosity_error)then
          if(iostatus.lt.0)then
            write(errlog(io),*)'ERROR:  EOR encountered'
          else
            write(errlog(io),*)'ERROR:  Error reading ilatlonflag from string'
            write(errlog(io),*)'           From the following line from the file: '
            write(errlog(io),*)linebuffer080
            write(errlog(io),*)'System Message: ',trim(adjustl(iomessage))
          endif
        endif;enddo
        stop 1
      endif

      if(ilatlonflag.eq.0) then
        ! expecting input variables to be in the same projection as
        ! specified by iprojflag and parameters
        IsLatLon          = .false.
      else
        ! expecting input variables to be in lat/lon
       IsLatLon          = .true.
      endif

      ! Set Projection Parameters
      if(IsLatLon.eqv..false.)then
        call PJ_Set_Proj_Params(Comp_projection_line)
        A3d_iprojflag  = PJ_iprojflag
        A3d_k0         = PJ_k0
        A3d_Re         = PJ_Re
        A3d_lam0       = PJ_lam0
        A3d_lam1       = PJ_lam1
        A3d_lam2       = PJ_lam2
        A3d_phi0       = PJ_phi0
        A3d_phi1       = PJ_phi1
        A3d_phi2       = PJ_phi2
      endif

      ! Read boundaries of model domain
      if(IsLatLon)then
        ! If input coordinates are in lat/lon, interpret lines as follows
        ! Block 1 Line 3
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        linebuffer050 = "Reading control file, Block 1, Line 3."
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        cdf_b1l3 = linebuffer080
        read(linebuffer080,*,err=9103,iostat=iostatus,iomsg=iomessage) lonLL, latLL            ! lat/lon of LL corner
        ! Block 1 Line 4
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        linebuffer050 = "Reading control file, Block 1, Line 4."
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        cdf_b1l4 = linebuffer080
        read(linebuffer080,*,err=9104,iostat=iostatus,iomsg=iomessage) gridwidth_e, gridwidth_n   ! Dimensions (in degrees) of the grid
        ! Block 1 Line 5
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage) linebuffer080
        linebuffer050 = "Reading control file, Block 1, Line 5."
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        cdf_b1l5 = linebuffer080
        read(linebuffer080,*,err=9105,iostat=iostatus,iomsg=iomessage) value1, value2   ! First read two values and flag
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) value1, value2, value3 ! Try for three
        if(iostatus.eq.0)then
          ! Successfully read 3 values; third is interpreted as elevation (in km)
          lon_volcano = value1
          lat_volcano = value2
          z_volcano   = value3
        else
          ! third value unsuccessful, assign vent elevation to 0
          lon_volcano = value1
          lat_volcano = value2
          z_volcano   = 0.0_ip
        endif
        ! Block 1 Line 6
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        linebuffer050 = "Reading control file, Block 1, Line 6."
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        cdf_b1l6 = linebuffer080
        read(linebuffer080,*,err=9106,iostat=iostatus,iomsg=iomessage) de, dn                 ! cell size in degrees 

        !Make sure longitudes are between 0 and 360 degrees
        if(lonLL.lt.-360.0_ip) then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)&
                  "Please give longitude values between -360 and 360."
          endif;enddo
          stop 1
        endif
        if(lonLL.lt.  0.0_ip) lonLL=lonLL+360.0_ip
        if(lonLL.ge.360.0_ip) lonLL=mod(lonLL,360.0_ip)
        if(lon_volcano.lt.-360.0) then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)&
                  "Please give longitude values between -360 and 360."
          endif;enddo
          stop 1
        endif
        if(lon_volcano.lt.  0.0_ip) lon_volcano = lon_volcano+360.0_ip
        if(lon_volcano.ge.360.0_ip) lon_volcano = mod(lon_volcano,360.0_ip)

        if(IsLatLon.and.&
           (gridwidth_e.ge.360.0_ip.or.&
            abs(gridwidth_e-360.0_ip).lt.EPS_TINY))then
          IsPeriodic  = .true.
          lonLL       = 0.0_ip
          gridwidth_e = 360.0_ip
        endif
        lonUR = lonLL + gridwidth_e
        latUR = latLL + gridwidth_n

        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),'(a13,f10.4)')'lonLL      = ',lonLL
          write(outlog(io),'(a13,f10.4)')'lonUR      = ',lonUR
          write(outlog(io),'(a13,f10.4)')'latLL      = ',latLL
          write(outlog(io),'(a13,f10.4)')'latUR      = ',latUR
          write(outlog(io),'(a13,f10.4)')'lon_volcano= ',lon_volcano
          write(outlog(io),'(a13,f10.4)')'lat_volcano= ',lat_volcano
          write(outlog(io),4) lonLL, latLL, gridwidth_e, gridwidth_n, &
                              lon_volcano, lat_volcano
          write(outlog(io),'(a13,f10.4,a3)')'z_volcano  = ',z_volcano,' km'
          write(outlog(io),5) de, dn
       endif;enddo

       !check for errors in input
        call LatLonChecker(latLL,lonLL,lat_volcano,lon_volcano,gridwidth_e,gridwidth_n)
      else  ! IsLatLon
        ! Block 1 Line 3
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        linebuffer050 = "Reading control file, Block 1, Line 3."
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        cdf_b1l3 = linebuffer080
        read(linebuffer080,*,err=9103,iostat=iostatus,iomsg=iomessage) xLL, yLL                ! LL corner in km 
        ! Block 1 Line 4
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        linebuffer050 = "Reading control file, Block 1, Line 4."
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        cdf_b1l4 = linebuffer080
        read(linebuffer080,*,err=9104,iostat=iostatus,iomsg=iomessage) gridwidth_x, gridwidth_y ! width and height of simulation area in km
        xUR = xLL + gridwidth_x
        yUR = yLL + gridwidth_y

        ! Block 1 Line 5
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        linebuffer050 = "Reading control file, Block 1, Line 5."
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        cdf_b1l5 = linebuffer080
        read(linebuffer080,*,err=9105,iostat=iostatus,iomsg=iomessage) value1, value2   ! First read two values and flag
                                                   ! an error if unable
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) value1, value2, value3 ! Try for three
        if(iostatus.eq.0)then
          x_volcano = value1
          y_volcano = value2
          z_volcano = value3
        else
          x_volcano = value1
          y_volcano = value2
          z_volcano = 0.0_ip
        endif
        ! Block 1 Line 6
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        linebuffer050 = "Reading control file, Block 1, Line 6."
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        cdf_b1l6 = linebuffer080
        read(linebuffer080,*,err=9106,iostat=iostatus,iomsg=iomessage) dx, dy                 ! cell size in horizontal, vertical, in km
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),'(a12,f10.4)')'xLL       = ',xLL
          write(outlog(io),'(a12,f10.4)')'xUR       = ',xUR
          write(outlog(io),'(a12,f10.4)')'yLL       = ',yLL
          write(outlog(io),'(a12,f10.4)')'yUR       = ',yUR
          write(outlog(io),'(a12,f10.4)')'x_volcano = ',x_volcano
          write(outlog(io),'(a12,f10.4)')'y_volcano = ',y_volcano
          write(outlog(io),4) xLL, yLL, gridwidth_x, gridwidth_y, &
                              x_volcano,y_volcano
          write(outlog(io),'(a12,f10.4,a3)')'z_volcano = ',z_volcano,' km'
          write(outlog(io),5) dx, dy
        endif;enddo
        call xyChecker(xLL,yLL,dx,dy,x_volcano,y_volcano,gridwidth_x,gridwidth_y)
      endif

      ! Block 1 Line 7
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading control file, Block 1, Line 7."
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      cdf_b1l7 = linebuffer080
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) dz_const    ! nodal spacing in z (always km)
      if(iostatus.eq.0)then
        ! numeric value read for nodal spacing in z; assume constant
        VarDzType = "dz_cons"
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),43) dz_const
        endif;enddo
        ! Set up initial z_vector up to 50km or so.  This is to match the variable
        ! dz cases in which the max height is specified.  The computational grid
        ! height will be truncated below to just that needed to cover the plume
        nz_init  = ceiling(100.0_ip/dz_const)+1
        allocate(z_vec_init(0:nz_init))
        z_vec_init = 0.0_ip
        do k=1,nz_init
          z_vec_init(k)=dz_const*(k) ! This the top of cell-boundaries (lower bound at 0)
        enddo
      else
        ! dz unsucessfully read, try to VarDzType stringg
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)&
                     "Could not read dz. Trying to reinterpret as alternate z-spacing"
          write(outlog(io),*)linebuffer080
        endif;enddo
        read(linebuffer080,*,err=9107,iostat=iostatus,iomsg=iomessage) VarDzType
        if(VarDzType.eq.'dz_plin')then
          ! Piece-wise linear
          !  Read another line with: n-segments, nz1, dz1, nz2, dz2, ...
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"z is piecewise linear:  Now reading the segments."
          endif;enddo
          ! Block 1 Line 7+1 (Reading the next line into cdf_vardz)
          read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer130
          linebuffer050 = "Reading control file, Block 1, Line 7+ (dz_plin)."
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer130(1:80),iomessage)
          cdf_vardz = linebuffer130
          read(linebuffer130,*,err=9107,iostat=iostatus,iomsg=iomessage) nsegments
          if(nsegments.lt.1)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: ",&
                           "nsegments must be positive integer"
              write(errlog(io),*)&
                           "       nsegments = ",nsegments
            endif;enddo
            stop 1
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Found n-segments: ",nsegments
            write(outlog(io),*)"      segment :    nz  :     dz "
          endif;enddo
          allocate(nz_plin_segments(nsegments))
          allocate(dz_plin_segments(nsegments))
          allocate(values(1+2*nsegments))
          read(linebuffer130,*,err=9107,iostat=iostatus,iomsg=iomessage)values(1:1+2*nsegments)
          do i=1,nsegments
            nz_plin_segments(i) = nint(values(1+(i-1)*2 + 1))
            dz_plin_segments(i) = values(1+(i-1)*2 + 2)
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),'(7x,i5,5x,i5,6x,f5.2)')i,nz_plin_segments(i),dz_plin_segments(i)
            endif;enddo
          enddo
          deallocate(values)
          nz_init = sum(nz_plin_segments(:))
          allocate(z_vec_init(0:nz_init))
          z_vec_init = 0.0_ip
          k = 0
          do i=1,nsegments
            do ii = 1,nz_plin_segments(i)
              if(k.eq.0)then
                last_z = 0.0_ip
              else
                last_z = z_vec_init(k)
              endif
              k=k+1
              z_vec_init(k)=last_z+dz_plin_segments(i)
            enddo
          enddo
        elseif(VarDzType.eq.'dz_clog')then
          ! constant log steps
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Logrithmic z (constant steps of dlog(z))"
            write(outlog(io),*)"Now reading the clog_zmax and number of steps."
          endif;enddo
          ! Block 1 Line 7+1 (Reading the next line into cdf_vardz)
          read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer130
          linebuffer050 = "Reading control file, Block 1, Line 7+ (dz_clog)."
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer130(1:80),iomessage)
          cdf_vardz = linebuffer130
          read(linebuffer130,*,err=9107,iostat=iostatus,iomsg=iomessage) clog_zmax, clog_nsteps
          if(clog_zmax.le.0.0)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: ",&
                           "z-level must be strictly positive"
              write(errlog(io),*)&
                           "       z-max = ",clog_zmax
            endif;enddo
            stop 1
          endif
          if(clog_nsteps.le.1)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: ",&
                           "Please choose more that 1 segments in log-z"
              write(errlog(io),*)&
                           "       n-segments = ",clog_nsteps
            endif;enddo
            stop 1
          endif
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Generating log profile with ",clog_nsteps," steps up to ",clog_zmax
          endif;enddo
          nz_init = clog_nsteps+1
          allocate(z_vec_init(0:nz_init))
          z_vec_init = 0.0_ip
          do k=1,nz_init
            tmp_ip = (       real(k,kind=ip)/real(clog_nsteps)) * log10(clog_zmax+1.0_ip)
            z_vec_init(k)=10.0**(tmp_ip) - 1.0_ip
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)k,z_vec_init(k)
            endif;enddo
          enddo
        elseif(VarDzType.eq.'dz_cust')then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Custom dz"
            write(outlog(io),*)"Now reading number of steps (ndz) followed by values(1:ndz)"
          endif;enddo
          ! Block 1 Line 7+1 (Reading the next line into cdf_vardz)
          read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer130
          linebuffer050 = "Reading control file, Block 1, Line 7+ (cust)."
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer130(1:80),iomessage)
          cdf_vardz = linebuffer130
          read(linebuffer130,*,err=9107,iostat=iostatus,iomsg=iomessage) cust_nsteps
          if(cust_nsteps.le.1)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: ",&
                           "Must have a postive number of steps in z"
              write(errlog(io),*)&
                           "       n steps = ",cust_nsteps
            endif;enddo
            stop 1
          endif
          allocate(values(cust_nsteps))
          read(linebuffer130,*,err=9107,iostat=iostatus,iomsg=iomessage) cust_nsteps, values(1:cust_nsteps)
          nz_init = cust_nsteps
          allocate(z_vec_init(0:nz_init))
          z_vec_init = 0.0_ip
          do k=1,nz_init
            z_vec_init(k)=z_vec_init(k-1)+values(k)
          enddo
        else
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)&
                  "dz type must be either a number (in km) for constant dz, or"
            write(errlog(io),*)&
                  "dz_plin, dz_clog, or dz_cust for variable dz"
            write(errlog(io),*)&
                  "You entered: ",cdf_b1l7
            write(errlog(io),*)&
                  "Interpreted as: ",VarDzType
          endif;enddo
          stop 1
        endif
      endif ! dz vs VarDzType

      ! Block 1 Line 8
      ! Read this line looking for diffusion coefficient and either a Suzuki constant, 
      ! or a plume type ('line', 'point', 'profile', 'umbrella', or 'umbrella_air')
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading control file, Block 1, Line 8."
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      cdf_b1l8 = linebuffer080
      ! First read the diffusivity
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) diffusivity_horz
      if(iostatus.ne.0)then
        ! Cannot read a valid diffusivity value
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)'ERROR:  Error reading diffusivity'
          write(errlog(io),*)'           From the following line from the file: '
          write(errlog(io),*)linebuffer080
          write(errlog(io),*)'System Message: ',trim(adjustl(iomessage))
        endif;enddo
        stop 1
      endif
      ! Now try both diffusivity and a Suzuki coefficient
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) diffusivity_horz, Suzuki_A
      if(iostatus.eq.0)then
        SourceType     = 'suzuki'
        SourceType_idx = 1
      else
        ! if the second item is not a number, read SourceType
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)&
            "Source type is not suzuki. Trying to read another standard type"
        endif;enddo
        read(linebuffer080,*,err=9108,iostat=iostatus,iomsg=iomessage) diffusivity_horz, SourceType
        if((SourceType.eq.'point').or. &
            (SourceType.eq.'Point').or. &
            (SourceType.eq.'POINT')) then
            SourceType     = 'point'
            SourceType_idx = 2
        elseif((SourceType.eq.'line').or. &
                   (SourceType.eq.'Line').or. &
                   (SourceType.eq.'LINE')) then
            SourceType     = 'line'
            SourceType_idx = 3
        elseif((SourceType.eq.'profile').or. &
                   (SourceType.eq.'Profile').or. &
                   (SourceType.eq.'PROFILE')) then
            SourceType     = 'profile'
            SourceType_idx = 4
        elseif((SourceType.eq.'umbrella').or. &
                   (SourceType.eq.'Umbrella').or. &
                   (SourceType.eq.'UMBRELLA')) then
            SourceType     = 'umbrella'
            SourceType_idx = 5
            Suzuki_A = SuzK_umb
        elseif((SourceType.eq.'umbrella_air').or. &
                   (SourceType.eq.'Umbrella_air').or. &
                   (SourceType.eq.'UMBRELLA_AIR')) then
            ! umbrella_air is the same as 'umbrella'
            ! but it is assumed to be an airborne run.
            ! Thus if gsbins=1, the MER is multiplied by 20
            ! to obtain the right rate of umbrella growth.
            SourceType     = 'umbrella_air'
            SourceType_idx = 6
            Suzuki_A = SuzK_umb
        else
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)&
             "SourceType is not point, line, profile, umbrella or umbrella_air."
            write(outlog(io),*)&
             "Assuming this is a custom source type."
            write(outlog(io),*)&
            "For now, just read eruptions start time, duration, and height."
          endif;enddo
          IsCustom_SourceType = .true.
        endif
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"  SourceType = ",SourceType
        endif;enddo
      endif
      ! convert diffusion coefficient from m2/s to km2/hr
      diffusivity_horz = diffusivity_horz*M2PS_2_KM2PHR
      diffusivity_vert = diffusivity_horz

      if(abs(diffusivity_horz).lt.EPS_SMALL)then
        useDiffusion = .false.
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Not using turbulent diffusivity."
        endif;enddo
      elseif(diffusivity_horz.lt.0.0)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                     "Diffusivity must be non-negative."
        endif;enddo
        stop 1
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),'(a39,f10.3,a5)')"Using constant turbulent diffusivity:  ",&
                  diffusivity_horz/M2PS_2_KM2PHR," m2/s"
        endif;enddo
        useDiffusion = .true.
      endif

      ! Block 1 Line 9
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading control file, Block 1, Line 9."
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      cdf_b1l9 = linebuffer080
      read(linebuffer080,*,err=9109,iostat=iostatus,iomsg=iomessage) neruptions  ! read in number of eruptions or pulses
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*) 'Expecting to read ',neruptions,&
                           ' eruptions lines in Block 2.'
      endif;enddo
      if(((SourceType.eq.'umbrella').or.(SourceType.eq.'umbrella_air')) &
           .and.(neruptions.gt.1)) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                'when SourceType=umbrella, neruptions must equal 1'
          write(errlog(io),*)&
                'You gave neruptions=',neruptions
          write(errlog(io),*)&
                'Program stopped'
        endif;enddo
        stop 1
      endif
      if(SourceType.eq.'suzuki')then
        if(Suzuki_A.le.0.0_ip)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: ",&
                  "Suzuki_A must be positive, not ",Suzuki_A
          endif;enddo
          stop 1
        endif
      else
        ! Source type is not Suzuki, so Suzuki_A may have been mangled trying
        ! to read something into it.  Reinitialize.
        Suzuki_A = SuzK_umb
      endif
      if(neruptions.le.0)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "neruptions must be positive, not ",neruptions
        endif;enddo
        stop 1
      endif
      ! END OF BLOCK 1
      !************************************************************************

      ! Allocate arrays of eruptive properties
      call Allocate_Source_eruption

      allocate (e_iyear(neruptions))
      allocate (e_imonth(neruptions))
      allocate (e_iday(neruptions))
      allocate (e_hour(neruptions))
      !************************************************************************
      ! BLOCK 2: ERUPTION PARAMETERS
      ! Again, assuming there is a variable length
      ! header with each header line flagged by a '#' or '*' in the first position
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading control file, past Blk1, looking for Blk2"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (.not.IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                'Expecting a comment line separating blocks.'
          write(errlog(io),*)&
                '       Check that Block 1 is correct.'
        endif;enddo
        stop 1
      endif
      do while (IsComment)
         ! Line is a comment, read next line
        read(fid_ctrlfile,'(a130)',iostat=iostatus,iomsg=iomessage)linebuffer130
        linebuffer050 = "Reading control file, past Blk1, looking for Blk2"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer130(1:80),iomessage)
        read(linebuffer130,*,iostat=iostatus,iomsg=iomessage)testkey
        linebuffer050 = "Reading testkey from linebuffer"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer130(1:80),iomessage)
        call FileIO_Check_testkey(testkey,linebuffer130(1:80),IsComment)
      enddo
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)' *******************************************'
        write(outlog(io),*)' Reading Block 2: Eruption parameters'
        write(outlog(io),*)' *******************************************'
      endif;enddo
      ! Begin reading times of eruptive pulses
      do i=1,neruptions  
        ! Always check if we have overshot the block
        read(linebuffer130,*,iostat=iostatus,iomsg=iomessage)testkey
        linebuffer050 = "Reading control file, Blk2, Line 1+."
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer130(1:80),iomessage)
        call FileIO_Check_testkey(testkey,linebuffer130(1:80),IsComment)
        if (IsComment) then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: ",&
                  'Trying to read Blk2 and detecting comment line'
            write(errlog(io),*)&
                  '  Eruption ',i,'of',neruptions
            write(errlog(io),*)&
                  '  Offending line: ',linebuffer130
          endif;enddo
          stop 1
        endif
        if(i.eq.1)then
          read(linebuffer130,*,err=9201,iostat=iostatus,iomsg=iomessage) e_iyear(i),e_imonth(i)
          if(e_iyear(i).ne.0.and.e_iyear(i).lt.BaseYear.or.e_iyear(i)-BaseYear.gt.100)then
            ! Reset BaseYear to the start of the century containing the eruption year
            MR_Comp_StartYear  = e_iyear(i)
            MR_Comp_StartMonth = e_imonth(i)
            BaseYear = e_iyear(i) - mod(e_iyear(i),100)
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"WARNING: Resetting BaseYear to ",BaseYear
            endif;enddo
          endif
          if(e_iyear(i).eq.0)then
            runAsForecast = .true.
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"Year = 0; Running as forecast."
            endif;enddo
          endif
        endif
        if(SourceType.eq.'suzuki'      .or. &
           SourceType.eq.'point'       .or. &
           SourceType.eq.'line'        .or. &
           SourceType.eq.'umbrella'    .or. &
           SourceType.eq.'umbrella_air')then
          ! read start time, duration, plume height, volume of each pulse
          read(linebuffer130,*,err=9201,iostat=iostatus,iomsg=iomessage) &
                                e_iyear(i),e_imonth(i),e_iday(i),e_hour(i), &
                                e_Duration(i), e_PlumeHeight(i), e_Volume(i)
        elseif(SourceType.eq.'profile')then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Start reading eruption profile number ",i
          endif;enddo
          ! read start time, duration, plume height, volume of each pulse
          read(linebuffer130,*,err=9201,iostat=iostatus,iomsg=iomessage) &
                                e_iyear(i),e_imonth(i),e_iday(i),e_hour(i), &
                                e_Duration(i), e_PlumeHeight(i), e_Volume(i),&
                                e_prof_dz(i),e_prof_nzpoints(i)
          if(e_prof_nzpoints(i).gt.MAX_ER_PROFPOINTS)then

          endif
          allocate(dum_prof(e_prof_nzpoints(i)))
          read(fid_ctrlfile,*,iostat=iostatus,iomsg=iomessage)dum_prof(1:e_prof_nzpoints(i))
          linebuffer080="Direct read of profile points; no line buffer"
          linebuffer050 = "Reading control file, Blk2, Line 1+ (profile)."
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
          ! Check to make sure the sum of the percentages add to 1.0
          if(abs(sum(dum_prof(1:e_prof_nzpoints(i)))-1.0_ip).gt.EPS_SMALL)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)&
               'ERROR:  The profile fractions do not sum to 1.0 for eruptive pulse ',i
              write(errlog(io),*)'             i          z (km)           %'
              do ii=1,e_prof_nzpoints(i)
                write(errlog(io),204)ii,e_prof_dz(i)*ii,dum_prof(ii)
 204            format(10x,i5,f15.3,f15.3)
              enddo
              write(errlog(io),'(a12,f5.3)')'     Sum =  ',sum(dum_prof(1:e_prof_nzpoints(i)))
            endif;enddo
            stop 1
          endif
          ! Check consistency between e_PlumeHeight(i) and the profile
          tmp_ip = 0.0_ip
          do ii=1,e_prof_nzpoints(i)
            if(dum_prof(ii).gt.EPS_SMALL)then
              tmp_ip = max(e_PlumeHeight(i),e_prof_dz(i)*ii)
            endif
          enddo
          if(tmp_ip.gt.e_PlumeHeight(i))then
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"Warning: Eruption pulse profile is higher than reported height."
              write(outlog(io),*)"         Reported height = ",e_PlumeHeight(i)
              write(outlog(io),*)"         Resetting height to ",tmp_ip
            endif;enddo
            e_PlumeHeight(i) = tmp_ip
          endif
          e_prof_Volume(i,1:e_prof_nzpoints(i))=dum_prof(1:e_prof_nzpoints(i))*e_Volume(i)
          if(e_prof_nzpoints(i).gt.e_prof_maxpoints)e_prof_maxpoints=e_prof_nzpoints(i)
          deallocate(dum_prof)
        else
          ! This is the custom source.  A special call to a source reader
          ! will need to made from Ash3d_??.F90.  For now, just read the
          ! start time, duration, and plume height
          read(linebuffer130,*,err=9201,iostat=iostatus,iomsg=iomessage)&
                                       e_iyear(i),e_imonth(i),e_iday(i),e_hour(i),&
                                       e_Duration(i), e_PlumeHeight(i)
          e_Volume(i)    = 0.0_ip
          if(neruptions.gt.1)then
            ! For more than one custom source, the next iteration might cause problems
            ! since a custom source might require multiple input lines per source.
            ! For now, copy slot 1 to all the others and break out of the do loop.
            ! The full source list must be populated by the user-provided custom source
            ! readers.
            e_iyear(2:neruptions)       = e_iyear(1)
            e_imonth(2:neruptions)      = e_imonth(1)
            e_iday(2:neruptions)        = e_iday(1)
            e_hour(2:neruptions)        = e_hour(1)
            e_Duration(2:neruptions)    = e_Duration(1)
            e_PlumeHeight(2:neruptions) = e_PlumeHeight(1)
            e_Volume(2:neruptions)      = e_Volume(1)
            exit
          endif

        endif ! SourceType

        ! if we're doing a forecast run, we'll need to add the windfile
        ! reference time to SimStartHour and e_StartTime(i) so
        ! temporarily set year and month to Jan, 1900.  If iday is
        ! considered 'days after start of wind file', then we need to
        ! add 1 so that the hours are calculated properly.
        if(runAsForecast)then
          e_iyear(i)  = BaseYear
          e_imonth(i) = 1
          e_iday(i)   = e_iday(i) + 1
          FC_Offset   = real(e_hour(1),sp)
        endif
        if(e_Duration(i).lt.0.0_ip)    e_Duration(i)    = ESP_duration
        if(e_PlumeHeight(i).lt.0.0_ip) e_PlumeHeight(i) = ESP_height
        if(e_Volume(i).lt.0.0_ip)      e_Volume(i)      = ESP_Vol

        ! read next line of input file
        read(fid_ctrlfile,'(a130)',iostat=iostatus,iomsg=iomessage)linebuffer130
        linebuffer050 = "Reading control file, next line of Blk2"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer130(1:80),iomessage)
      enddo

      !Error trap if more pulses are entered than are specified
      if(IsCustom_SourceType)then
        ! If we are using custom sources, suppress the error-checking
        ! since we don't know the number of lines of input we need and
        ! we will need to loop through here again anyway.
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)&
            'Skipping error-checking of eruption source lines since '
          write(outlog(io),*)&
            'the source is a custom type'
        endif;enddo
        ! For the custom source, we will need to read to the end of block 2
        call FileIO_Check_testkey(testkey,linebuffer130(1:80),IsComment)
        do while (IsComment)
           ! Line is a comment, read next line
          read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
          linebuffer050 = "Reading ctr file, past Blk 2, looking for Blk 3."
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
          read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
          linebuffer050 = "Reading testkey from linebuffer"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
          call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
        enddo
      else
        if(linebuffer130(1:5).ne.'*****') then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)&
             'The beginning of the line following the list of', &
             ' eruptive pulses did not'
            write(errlog(io),*)&
             'start with ''*****''.  Did you enter the correct', &
             '  number of eruptive pulses?'
            write(errlog(io),*) 'Program stopped.'
          endif;enddo
          stop 1
        endif
      endif
      ! END OF BLOCK 2
      !************************************************************************

      !************************************************************************
      ! BLOCK 3: WIND PARAMETERS
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading control file, Blk3 Line 1"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk3)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if(.not.IsCustom_SourceType.and..not.IsComment) then !&  ! only perform this check for standard src
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                'Expecting a comment line separating blks.'
          write(errlog(io),*)&
                '       Check that Block 2 is correct.'
        endif;enddo
        stop 1
      endif
      do while (IsComment)
         ! Line is a comment, read next line
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        linebuffer050 = "Reading control file, Blk3 Line 1"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
        linebuffer050 = "Reading testkey from linebuffer (Blk3)"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      enddo
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)' *******************************************'
        write(outlog(io),*)' Reading Block 3: Windfile parameters'
        write(outlog(io),*)' *******************************************'
      endif;enddo
      ! Block 3 Line 1
      cdf_b3l1 = linebuffer080
      ! Try to read at least two integers (iwind and iwindformat) or throw and error
      read(linebuffer080,*,err=9301,iostat=iostatus,iomsg=iomessage) iw,iwf
      ! Note: the validity of iw and iwf will be checked in the call to MR_Allocate_FullMetFileList
      idf = 0
      ! Succeeded in reading the two required values, try for three
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) iw, iwf, ivalue3
      if(iostatus.eq.0)then
        ! Success reading three values, try for four
        igrid = ivalue3
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) iw, iwf, ivalue3, ivalue4
        if(iostatus.eq.0)then
          ! Success!, set data format (ascii, netcdf, grib)
          idf = ivalue4
        endif
      else
        igrid = 0
      endif
      if(idf.lt.1)then
        ! Data format is not given, assume netcdf unless ascii specified
        if(iw.eq.1.or.iw.eq.2)then
          idf = 1 ! ASCII
        else
          idf = 2 ! Netcdf
        endif
      endif

      if(iwf.eq.0)then
        ! If iwindformat = 0, then the input file is a not a known format
        ! Read an extra line given the name of a template file.
        if(idf.ne.2)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)&
               "Currently only netcdf reader implemented for templates.",&
               "  Resetting idf to 2"
          endif;enddo
          idf = 2
        endif
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        linebuffer050 = "Reading control file, Blk3 Line 1+ (template)"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
        linebuffer050 = "Reading testkey from linebuffer (Blk3)"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
        if (IsComment) then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: ",&
                  "Trying to read template name and detecting comment line"
          endif;enddo
          stop 1
        endif
        read(linebuffer080,'(a80)',err=93011,iostat=iostatus,iomsg=iomessage) MR_iwf_template
      else
        ! iwf is a known format.
        if(iwf.eq.33)then
            ! The only known format that does not use leap years is the
            ! paleoclimate CAM files.
          useLeap=.false.
        else
          useLeap=.true.
        endif
      endif

      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading from control file, Bloc 3 line 2, iHeight"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk3)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Trying to read MR_iHeightHandler and detecting comment line"
        endif;enddo
        stop 1
      endif
      ! Block 3 Line 2
      cdf_b3l2 = linebuffer080
      read(linebuffer080,*,err=9302,iostat=iostatus,iomsg=iomessage)&
                      MR_iHeightHandler ! parameter that determines what to do if the
                                        ! plume height exceeds the wind sounding max. height

      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading from control file, Bloc 3 line 3, SimTime"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk3)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Trying to read Simtime_in_hours and detecting comment line"
        endif;enddo
        stop 1
      endif
      ! Block 3 Line 3
      cdf_b3l3 = linebuffer080
      read(linebuffer080,*,err=9303,iostat=iostatus,iomsg=iomessage)&
                                     Simtime_in_hours   ! simulated transport time
                                                        ! for ash cloud, in hours

      ! Read whether to stop calculation when percent_accumulated>0.99
      ! Block 3 Line 4
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage) linebuffer080
      linebuffer050 = "Reading from control file, Bloc 3 line 4"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk3)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Trying to read StopWhenDeposited and detecting comment line"
        endif;enddo
        stop 1
      endif
      ! Block 3 Line 4
      cdf_b3l4 = linebuffer080
      read(linebuffer080,'(a3)',err=9304,iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)).eq.'yes') then
        StopWhenDeposited = .true.
        StopValue_FracAshDep = StopValue_FracAshDep_Default
       else if(adjustl(trim(answer(1:2))).eq.'no') then
        StopWhenDeposited = .false.
        StopValue_FracAshDep = 1.0e2_ip
       else
        goto 9304
      endif

      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading from control file, Bloc 3 line 5"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk3)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment) 
      if (IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Trying to read nWindFiles and detecting comment line"
        endif;enddo
        stop 1
      endif
      ! Block 3 Line 5
      cdf_b3l5 = linebuffer080
      read(linebuffer080,*,err=9305,iostat=iostatus,iomsg=iomessage) iwfiles      ! number of wind files to read

      ! Now that we know which calendar we are using (BaseYear, useLeap), now we
      ! can set the HoursSince time for the source terms
      do i=1,neruptions
        !set start time of simulation
        ! Note: This will be updated for forecast runs once we know the start
        !       time of the windfiles
        if(i.eq.1) then
          SimStartHour = HS_hours_since_baseyear(e_iyear(i),e_imonth(i),  &
                        e_iday(i),e_hour(i),BaseYear,useLeap)
          xmlSimStartTime = HS_xmltime(SimStartHour,BaseYear,useLeap)
          MR_Comp_StartHour     = SimStartHour
          MR_Comp_Time_in_hours = Simtime_in_hours
        endif
        e_StartTime(i) = HS_hours_since_baseyear(e_iyear(i),e_imonth(i),  &
                                e_iday(i),e_hour(i),BaseYear,useLeap) - SimStartHour
        ! error trap if eruptions are not in chronological order
        if(.not.IsCustom_SourceType)then
          ! relax the chronological requirement for custom sources
          if(i.ge.2)then
            !(add 0.001 hours to make sure that rounding error does not cause
            !the program to stop)
            if((e_StartTime(i)+0.001_ip).lt.(e_StartTime(i-1)+e_Duration(i-1))) goto 9202
          endif
        endif
      enddo

      ! write out eruption information
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),7) neruptions
      endif;enddo
      do i=1,neruptions
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),8) i, e_PlumeHeight(i), e_iyear(i), e_imonth(i), &
                     e_iday(i), e_hour(i), e_Duration(i), e_Volume(i)
        endif;enddo
        if(SourceType.eq.'profile')then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)'             i          z (km)           km3'
            do ii=1,e_prof_nzpoints(i)
              write(outlog(io),205)ii,e_prof_dz(i)*ii,e_prof_Volume(i,ii)
            enddo
 205            format(10x,i5,f15.3,f15.7)
            write(outlog(io),'(a39,g12.5,a7)')"         Total Volume for this pulse = ",&
                                sum(e_prof_Volume(i,:))," km3 DRE"
          endif;enddo
        endif
      enddo
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),'(a32,g12.5,a8)')"Total volume of all eruptions = ",&
                            sum(e_volume)," km3 DRE"
      endif;enddo

      ! Now that we know the requested dz profile and the plume heights, we can
      ! set up the z-grid for computation
      Ztop = ZPADDING*maxval(e_PlumeHeight(1:neruptions))
      MR_ztop         = real(Ztop,kind=sp)   ! Set the MetReader copy in case we scale the grid
      nzmax = 0
      do k = 1,nz_init-1
        if(z_vec_init(k+1).gt.Ztop.and. &
           z_vec_init(k).le.Ztop)then
          nzmax = k
        endif
      enddo
      if(nzmax.eq.0)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                   "Specified z-grid does not extend high enough"
          write(errlog(io),*)"        for given plume heights."
          write(errlog(io),*)"    e_PlumeHeight = ",e_PlumeHeight(1:neruptions)
          write(errlog(io),*)"         ZPADDING = ",ZPADDING
          write(errlog(io),*)"  CompGrid_height = ",Ztop
          write(errlog(io),*)"       z_vec_init = "
          do k = 1,nz_init-1
            if(maxval(e_PlumeHeight(1:neruptions)).gt.z_vec_init(k).and.&
               maxval(e_PlumeHeight(1:neruptions)).lt.z_vec_init(k+1))then
              kvent=k
            endif
            write(errlog(io),*)"                ",k,real(z_vec_init(k),kind=4)
          enddo
        endif;enddo
        stop 1
      endif

      ! Error-check initial wind info and allocate windfile arrays
      MR_BaseYear = BaseYear
      MR_useLeap  = useLeap
      call MR_Allocate_FullMetFileList(iw,iwf,igrid,idf,iwfiles)
      if(MR_useLeap.neqv.useLeap)then
        useLeap  = MR_useLeap
        BaseYear = MR_BaseYear
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Change in calandar; resetting e_StartTime"
        endif;enddo
        tmp_dp = HS_hours_since_baseyear(e_iyear(1),e_imonth(1),  &
                         e_iday(1),e_hour(1),BaseYear,useLeap)
        tmp_dp = tmp_dp - SimStartHour   ! Recast tmp_dp as the difference in calandars
        SimStartHour = SimStartHour + tmp_dp
        xmlSimStartTime = HS_xmltime(SimStartHour,BaseYear,useLeap)
      endif
      if(SourceType.eq.'suzuki') then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),6) diffusivity_horz, Suzuki_A, &
                               StopWhenDeposited, Simtime_in_hours
        endif;enddo
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),1438) diffusivity_horz, &
                                  StopWhenDeposited, Simtime_in_hours
        endif;enddo
      endif
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),1439) SourceType
      endif;enddo

      ! Find the i,j index values of the node containing the source volcano
      if(IsLatLon) then
        ivent = int((lon_volcano-lonLL)/de) + 1
        jvent = int((lat_volcano-latLL)/dn) + 1
      else
        ivent = int((x_volcano-xLL)/dx) + 1
        jvent = int((y_volcano-yLL)/dy) + 1
      endif

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),104) ivent,jvent
      endif;enddo
104   format(4x,'i and j coordinates of volcano:',/, &
             4x,'i=',i4,/, &
             4x,'j=',i4,/)

      ! Now back to reading the input file
      !************************************************************************
      ! BLOCK 4: OUTPUT OPTIONS
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading ctr file, past Blk 3, looking for Blk 4"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk4)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (.not.IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Expecting a comment line separating blocks."
          write(errlog(io),*)'       Check that Block 3 is correct.'
        endif;enddo
        stop 1
      endif
      do while (IsComment)
         ! Line is a comment, read next line
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        linebuffer050 = "Reading test line from control file, blk4"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
        linebuffer050 = "Reading testkey from linebuffer (Blk4)"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      enddo
      ! Block 4
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)' *******************************************'
        write(outlog(io),*)' Reading Block 4: Output options '
        write(outlog(io),*)' *******************************************'
      endif;enddo
      WriteDepositFinal_ASCII       = .false.
      WriteDepositFinal_KML         = .false.
      WriteDepositTS_ASCII          = .false.
      WriteDepositTS_KML            = .false.
      WriteCloudConcentration_ASCII = .false.
      WriteCloudConcentration_KML   = .false.
      WriteCloudHeight_ASCII        = .false.
      WriteCloudHeight_KML          = .false.
      WriteCloudLoad_ASCII          = .false.
      WriteCloudLoad_KML            = .false.
      WriteDepositTime_ASCII        = .false.
      WriteDepositTime_KML          = .false.
      WriteCloudTime_ASCII          = .false.
      WriteCloudTime_KML            = .false.
      WriteAirportFile_ASCII        = .false.
      WriteAirportFile_KML          = .false.
      Write_PT_Data                 = .false.
      Write_PR_Data                 = .false.
      cdf_b4l1 = linebuffer080
      ! Block 4 Line 1
      ! Read whether to write out final ESRI ASCII deposit file
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk4)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Trying to read WriteDepositFinal_ASCII and detecting",& 
                " comment line"
        endif;enddo
        stop 1
      endif
      read(linebuffer080,'(a3)',err=9401,iostat=iostatus,iomsg=iomessage) answer
      linebuffer050 = "Reading control file Blk 4, line 1"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      if(adjustl(trim(answer)).eq.'yes') then
        WriteDepositFinal_ASCII = .true.
       else if(adjustl(trim(answer(1:2))).eq.'no') then
        WriteDepositFinal_ASCII = .false.
       else
        goto 9401
      endif

      ! Block 4 Line 2
      ! Read whether to write out final KML deposit file
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading control file Blk 4, line 2"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk4)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Trying to read WriteDepositFinal_KML and detecting", &
                " comment line"
        endif;enddo
        stop 1
      endif
      read(linebuffer080,'(a3)',err=9402,iostat=iostatus,iomsg=iomessage) answer
      linebuffer050 = "Reading control file Blk 4, line 2"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      cdf_b4l2 = linebuffer080
      if(adjustl(trim(answer)).eq.'yes') then
        WriteDepositFinal_KML = .true.
       else if(adjustl(trim(answer(1:2))).eq.'no') then
        WriteDepositFinal_KML = .false.
       else
        goto 9402
      endif

      ! Block 4 Line 3
      ! Read whether to write out ESRI ASCII deposit files at specified times
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading control file Blk 4, line 3"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk4)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Trying to read WriteDepositTS_ASCII and detecting",&
                " comment line"
        endif;enddo
        stop 1
      endif
      cdf_b4l3 = linebuffer080
      read(linebuffer080,'(a3)',err=9403,iostat=iostatus,iomsg=iomessage) answer
      linebuffer050 = "Reading control file Blk 4, line 3"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      if(adjustl(trim(answer)).eq.'yes') then
        WriteDepositTS_ASCII = .true.
       else if(adjustl(trim(answer(1:2))).eq.'no') then
        WriteDepositTS_ASCII = .false.
       else
        goto 9403
      endif

      ! Block 4 Line 4
      ! Read whether to write out KML deposit files at specified times
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading control file Blk 4, line 4"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk4)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Trying to read WriteDepositTS_KML and detecting comment line"
        endif;enddo
        stop 1
      endif
      cdf_b4l4 = linebuffer080
      read(linebuffer080,'(a3)',err=9404,iostat=iostatus,iomsg=iomessage) answer
      linebuffer050 = "Reading control file Blk 4, line 4"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      if(adjustl(trim(answer)).eq.'yes') then
        WriteDepositTS_KML = .true.
       else if(adjustl(trim(answer(1:2))).eq.'no') then
        WriteDepositTS_KML = .false.
       else
        goto 9404
      endif

      ! Block 4 Line 5
      ! Read whether to write out ESRI ASCII files of cloud concentration
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading control file Blk 4, line 5"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk4)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Trying to read WriteCloudConcentration_ASCII and detecting",&
                " comment line"
        endif;enddo
        stop 1
      endif
      cdf_b4l5 = linebuffer080
      read(linebuffer080,'(a3)',err=9405,iostat=iostatus,iomsg=iomessage) answer
      linebuffer050 = "Reading control file Blk 4, line 5"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      if(adjustl(trim(answer)).eq.'yes') then
        WriteCloudConcentration_ASCII = .true.
       else if(adjustl(trim(answer(1:2))).eq.'no') then
        WriteCloudConcentration_ASCII = .false.
       else
        goto 9405
      endif

      ! Block 4 Line 6
      ! Read whether to write out KML files of cloud concentration
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading control file Blk 4, line 6"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk4)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Trying to read WriteCloudConcentration_KML and detecting",&
                " comment line"
        endif;enddo
        stop 1
      endif
      cdf_b4l6 = linebuffer080
      read(linebuffer080,'(a3)',err=9406,iostat=iostatus,iomsg=iomessage) answer
      if(adjustl(trim(answer)).eq.'yes') then
        WriteCloudConcentration_KML = .true.
       else if(adjustl(trim(answer(1:2))).eq.'no') then
        WriteCloudConcentration_KML = .false.
       else
        goto 9406
      endif

      ! Block 4 Line 7
      ! Read whether to write out ESRI ASCII files of cloud height
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading control file Blk 4, line 7"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk4)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Trying to read WriteCloudHeight_ASCII and detecting",&
                " comment line"
        endif;enddo
        stop 1
      endif
      cdf_b4l7 = linebuffer080
      read(linebuffer080,'(a3)',err=9407,iostat=iostatus,iomsg=iomessage) answer
      linebuffer050 = "Reading control file Blk 4, line 7"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      if(adjustl(trim(answer)).eq.'yes') then
        WriteCloudHeight_ASCII = .true.
       else if(adjustl(trim(answer(1:2))).eq.'no') then
        WriteCloudHeight_ASCII = .false.
       else
        goto 9407
      endif

      ! Block 4 Line 8
      ! Read whether to write out KML files of cloud height
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading control file Blk 4, line 8"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk4)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Trying to read WriteCloudHeight_KML and detecting",&
                " comment line"
        endif;enddo
        stop 1
      endif
      cdf_b4l8 = linebuffer080
      read(linebuffer080,'(a3)',err=9408,iostat=iostatus,iomsg=iomessage) answer
      linebuffer050 = "Reading control file Blk 4, line 8"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      if(adjustl(trim(answer)).eq.'yes') then
        WriteCloudHeight_KML = .true.
       else if(adjustl(trim(answer(1:2))).eq.'no') then
        WriteCloudHeight_KML = .false.
       else
        goto 9408
      endif

      ! Block 4 Line 9
      ! Read whether to write out ESRI ASCII files of ashcloud load
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading control file Blk 4, line 9"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk4)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Trying to read WriteCloudLoad_ASCII and detecting",&
                " comment line"
        endif;enddo
        stop 1
      endif
      cdf_b4l9 = linebuffer080
      read(linebuffer080,'(a3)',err=9409,iostat=iostatus,iomsg=iomessage) answer
      linebuffer050 = "Reading control file Blk 4, line 9"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      if(adjustl(trim(answer)).eq.'yes') then
        WriteCloudLoad_ASCII = .true.
       else if(adjustl(trim(answer(1:2))).eq.'no') then
        WriteCloudLoad_ASCII = .false.
       else
        goto 9409
      endif

      ! Block 4 Line 10
      ! Read whether to write out KML files of ashcloud load
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading control file Blk 4, line 10"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk4)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Trying to read WriteCloudLoad_KML and detecting comment line"
        endif;enddo
        stop 1
      endif
      cdf_b4l10 = linebuffer080
      read(linebuffer080,'(a3)',err=9410,iostat=iostatus,iomsg=iomessage) answer
      linebuffer050 = "Reading control file Blk 4, line 10"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      if(adjustl(trim(answer)).eq.'yes') then
        WriteCloudLoad_KML = .true.
       else if(adjustl(trim(answer(1:2))).eq.'no') then
        WriteCloudLoad_KML = .false.
       else
        goto 9410
      endif

      ! Block 4 Line 11
      ! Read whether to write out ASCII file of deposit arrival time
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading control file Blk 4, line 11"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk4)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Trying to read WriteDepositTime_ASCII and detecting",&
                " comment line"
        endif;enddo
        stop 1
      endif
      cdf_b4l11 = linebuffer080
      read(linebuffer080,'(a3)',err=9411,iostat=iostatus,iomsg=iomessage) answer
      linebuffer050 = "Reading control file Blk 4, line 11"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      if(adjustl(trim(answer)).eq.'yes') then
        WriteDepositTime_ASCII = .true.
       else if(adjustl(trim(answer(1:2))).eq.'no') then
        WriteDepositTime_ASCII = .false.
       else
        goto 9411
      endif

      ! Block 4 Line 12
      ! Read whether to write out KML files of deposit arrival time
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading control file Blk 4, line 12"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk4)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Trying to read WriteDepositTime_KML and detecting",&
                " comment line"
        endif;enddo
        stop 1
      endif
      cdf_b4l12 = linebuffer080
      read(linebuffer080,'(a3)',err=9412,iostat=iostatus,iomsg=iomessage) answer
      linebuffer050 = "Reading control file Blk 4, line 12"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      if(adjustl(trim(answer)).eq.'yes') then
        WriteDepositTime_KML = .true.
       else if(adjustl(trim(answer(1:2))).eq.'no') then
        WriteDepositTime_KML = .false.
       else
        goto 9412
      endif

      ! Block 4 Line 13
      ! Read whether to write out ESRI ASCII file of cloud arrival time
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading control file Blk 4, line 13"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk4)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Trying to read WriteCloudTime_ASCII and detecting",&
                " comment line"
        endif;enddo
        stop 1
      endif
      cdf_b4l13 = linebuffer080
      read(linebuffer080,'(a3)',err=9413,iostat=iostatus,iomsg=iomessage) answer
      linebuffer050 = "Reading control file Blk 4, line 13"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      if(adjustl(trim(answer)).eq.'yes') then
        WriteCloudTime_ASCII = .true.
       else if(adjustl(trim(answer(1:2))).eq.'no') then
        WriteCloudTime_ASCII = .false.
       else
        goto 9413
      endif

      ! Block 4 Line 14
      ! Read whether to write out KML files of cloud arrival time
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading control file Blk 4, line 14"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk4)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Trying to read WriteCloudTime_KML and detecting comment line"
        endif;enddo
        stop 1
      endif
      cdf_b4l14 = linebuffer080
      read(linebuffer080,'(a3)',err=9414,iostat=iostatus,iomsg=iomessage) answer
      linebuffer050 = "Reading control file Blk 4, line 14"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      if(adjustl(trim(answer)).eq.'yes') then
        WriteCloudTime_KML = .true.
       else if(adjustl(trim(answer(1:2))).eq.'no') then
        WriteCloudTime_KML = .false.
       else
        goto 9414
      endif

      ! Block 4 Line 15
      ! Read whether to write out 3D files of ash concentration
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading control file Blk 4, line 15"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk4)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Trying to read Write3dFiles and detecting comment line"
        endif;enddo
        stop 1
      endif
      cdf_b4l15 = linebuffer080
      read(linebuffer080,'(a3)',err=9415,iostat=iostatus,iomsg=iomessage) answer
      linebuffer050 = "Reading control file Blk 4, line 15"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      useRestartVars = .true.
      if(adjustl(trim(answer)).eq.'yes') then
        Write3dFiles = .true.
        ! if a consolidated output file will be written, assume both standard
        ! variables and the 3d ash concentrations will be written

        ! Try to read an output code
        loc = index(linebuffer080,'yes')
        dumstr20 = linebuffer080(loc+4:loc+24)
        read(dumstr20,*,iostat=iostatus,iomsg=iomessage) iform
        if(iostatus.eq.0)then
          ! Succeeded in reading the format code
          if(iform.eq.1)then
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)" Successfully read format code=1"
              write(outlog(io),*)"  Both output products and ash concentrations will be written"
            endif;enddo
            useRestartVars = .true.
          elseif(iform.eq.2)then
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)" Successfully read format code=2"
              write(outlog(io),*)"  Only output products will be written"
            endif;enddo
            useRestartVars = .false.
          else
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)" Could not read format code"
              write(outlog(io),*)"  Assuming both output products and ash concentration will be written"
            endif;enddo
            useRestartVars = .true.
          endif
        endif
       else if(adjustl(trim(answer(1:2))).eq.'no') then
        Write3dFiles = .false.
       else
        goto 9415
      endif

      ! Block 4 Line 16
      ! Read output file format
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading control file Blk 4, line 16"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk4)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Trying to read output format and detecting comment line"
        endif;enddo
        stop 1
      endif
      cdf_b4l16 = linebuffer080
      linebuffer080=adjustl(linebuffer080)
      if(Write3dFiles) then
        read(linebuffer080,'(a6)',err=9416,iostat=iostatus,iomsg=iomessage) formatanswer
        linebuffer050 = "Reading control file Blk 4, line 16"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        if(formatanswer(1:5).eq.'ascii') then
          ioutputFormat = 1
        else if(formatanswer(1:6).eq.'binary') then
          ioutputFormat = 2
        else if(formatanswer(1:6).eq.'netcdf') then
          ioutputFormat = 3
        else
          goto 9416
        endif
      endif
      
      ! Block 4 Line 17
      ! Read number of output steps to write out
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage) linebuffer080
      linebuffer050 = "Reading control file Blk 4, line 17"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk4)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Trying to read nWriteTimes and detecting comment line"
        endif;enddo
        stop 1
      endif
      cdf_b4l17 = linebuffer080

      ! Block 4 Line 18
      read(fid_ctrlfile,'(a400)',iostat=iostatus,iomsg=iomessage)linebuffer400
      linebuffer050 = "Reading control file Blk 4, line 18"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer400(1:80),iomessage)
      read(linebuffer400,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk4)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer400(1:80),iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Trying to read WriteTimes or WriteInterval and detecting",&
                " comment line"
        endif;enddo
        stop 1
      endif
      cdf_b4l18 = linebuffer400(1:80)
      if(WriteDepositFinal_ASCII      .or. &
          WriteDepositFinal_KML        .or. &
          WriteDepositTS_ASCII         .or. &
          WriteDepositTS_KML           .or. &
          WriteCloudConcentration_ASCII.or. &
          WriteCloudConcentration_KML  .or. &
          Write3dFiles                 .or. &
          WriteCloudHeight_ASCII       .or. &
          WriteCloudHeight_KML         .or. &
          WriteCloudLoad_KML) then
        read(linebuffer080,*,err=9417,iostat=iostatus,iomsg=iomessage) nWriteTimes
          ! Check how to interpret nWriteTimes
        if(nWriteTimes.gt.0) then
          ! If a positive number, then we're reading an array of times
          allocate(WriteTimes(nWriteTimes))
          read(linebuffer400,*,err=9418,iostat=iostatus,iomsg=iomessage) WriteTimes(1:nWriteTimes)
        elseif(nWriteTimes.eq.0) then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"nWriteTimes = 0: Running without output"
          endif;enddo
        elseif(nWriteTimes.ne.-1) then
          ! If not a positive number, then it should be -1
          ! Report error otherwise
          goto 9417
        else
          ! If nWriteTimes=-1, then read a single WriteTimes and interpret it as a time interval
          read(cdf_b4l18,*,err=94181,iostat=iostatus,iomsg=iomessage) WriteInterval
          ! Redefine nWriteTimes since it was read in as -1
          nWriteTimes = int(Simtime_in_hours/WriteInterval)+1
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),'(a21,f10.3)')"  WriteInterval    = ",WriteInterval
            write(outlog(io),'(a21,f10.3)')"  Simtime_in_hours = ",Simtime_in_hours
            write(outlog(io),'(a21,i3)')"  nWriteTimes      = ",nWriteTimes
          endif;enddo
          allocate(WriteTimes(nWriteTimes))

          forall (i=1:nWriteTimes)  WriteTimes(i) = (i-1)*WriteInterval
          do i=1,nWriteTimes
            WriteTimes(i) = (i-1)*WriteInterval
          enddo
          do i=1,nWriteTimes     !check writetimes for errors
            if(WriteTimes(i).lt.0.0_ip) then  !if the time <0
                ! Abort the program
              goto 94182
            elseif(i.gt.1)then
              if(WriteTimes(i).lt.Writetimes(i-1))then !if times are not in chronological order
                  ! Abort the program
                goto 94183
              endif
            elseif(WriteTimes(i).gt.Simtime_in_hours) then   !if some times exceed the simulation time
              do io=1,2;if(VB(io).le.verbosity_info)then
                write(outlog(io),32)
              endif;enddo
              nWriteTimes = i-1
              exit
            endif
          enddo
        endif

        if(LoadConcen)then
          ! Find the output time index that is next
          if(time.le.WriteTimes(1))then
            iTimeNext = 1
          else
            do i = 2, nWriteTimes
              if(time.ge.WriteTimes(i-1).and. &
                 time.lt.WriteTimes(i))then
                iTimeNext = i
              endif
            enddo
          endif
        else
          iTimeNext = 1
        endif
        NextWriteTime = WriteTimes(iTimeNext)
      endif

      ! Write out the types of output to be written
      ! output options
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),33) WriteDepositFinal_ASCII,       &
                             WriteDepositFinal_KML,         &
                             WriteDepositTS_ASCII,          &
                             WriteDepositTS_KML,            &
                             WriteCloudConcentration_ASCII, &
                             WriteCloudConcentration_KML,   &
                             WriteCloudHeight_ASCII,        &
                             WriteCloudHeight_KML,          &
                             WriteCloudLoad_ASCII,          &
                             WriteCloudLoad_KML,            &
                             WriteDepositTime_ASCII,        &
                             WriteDepositTime_KML,          &
                             WriteCloudTime_ASCII,          &
                             WriteCloudTime_KML,            &
                             Write3dFiles,                  &
                             formatanswer,                  &
                             nWriteTimes
      endif;enddo
      if(WriteDepositTS_ASCII          .or. &
          WriteDepositTS_KML            .or. &
          WriteCloudConcentration_ASCII .or. &
          WriteCloudConcentration_KML   .or. &
          WriteCloudHeight_ASCII        .or. &
          WriteCloudHeight_KML          .or. &
          WriteCloudLoad_ASCII          .or. &
          WriteCloudLoad_KML            .or. &
          WriteCloudTime_ASCII          .or. &
          WriteCloudTime_KML            .or. &
          Write3dFiles) then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),34)                                !write out the types of output specified
          do i=1,nWriteTimes
            write(outlog(io),35) WriteTimes(i)                !write out the times when  files will be written
          enddo
        endif;enddo
      endif
      ! END OF BLOCK 4
      !************************************************************************

      !************************************************************************
      ! BLOCK 5: INPUT WIND FILES
      if(MR_iWindFiles.gt.0)then
        read(fid_ctrlfile,'(a130)',iostat=iostatus,iomsg=iomessage)linebuffer130
        linebuffer050 = "Reading ctr file, past Blk 4, looking for Blk 5."
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer130(1:80),iomessage)
        read(linebuffer130,*,iostat=iostatus,iomsg=iomessage)testkey
        linebuffer050 = "Reading testkey from linebuffer (Blk5)"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer130(1:80),iomessage)
        call FileIO_Check_testkey(testkey,linebuffer130(1:80),IsComment)
        if (.not.IsComment) then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: ",&
                  "Expecting a comment line separating blocks."
            write(errlog(io),*)'       Check that Block 4 is correct.'
          endif;enddo
          stop 1
        endif      
        do while (IsComment)
           ! Line is a comment, read next line
          !iostatus=0
          !iomessage='M'
          read(fid_ctrlfile,'(a130)',iostat=iostatus,iomsg=iomessage)linebuffer130
          linebuffer050 = "Reading test line of blk5"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer130(1:80),iomessage)
          read(linebuffer130,*,iostat=iostatus,iomsg=iomessage)testkey
          linebuffer050 = "Reading testkey from linebuffer (Blk5)"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer130(1:80),iomessage)
          call FileIO_Check_testkey(testkey,linebuffer130(1:80),IsComment)
        enddo
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)' *****************************************'
          write(outlog(io),*)' Reading Block 5: Windfile names'
          write(outlog(io),*)' *****************************************'
          write(outlog(io),13)
        endif;enddo
          ! Read list of windfiles.
        if(MR_iwind.eq.5)then
          ! For NCEP 2.5 degree (25), NOAA product (27), ERA5 (29), or ERA-20C (30)
          ! just read the path to the files
          read(linebuffer130,'(a130)',err=9501,iostat=iostatus,iomsg=iomessage) MR_WindFiles(1)
          cdf_b5l1(1:80) = linebuffer130(1:80)
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),1034) 1,trim(adjustl(MR_WindFiles(1)))
          endif;enddo
          read(fid_ctrlfile,'(a130)',iostat=iostatus,iomsg=iomessage)linebuffer130
          linebuffer050 = "Reading control file Blk 5, line 1"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer130(1:80),iomessage)
        else
          ! For all other iwf (MR_iwindformats), read the full list
          do i=1,iwfiles
            ! Always check if we have overshot the block
            testkey=linebuffer130(1:1)
            call FileIO_Check_testkey(testkey,linebuffer130(1:80),IsComment)
            if (IsComment) then
              do io=1,2;if(VB(io).le.verbosity_error)then
                write(errlog(io),*)"ERROR: ",&
                     "Trying to read Block 5 and detecting comment line"
                write(errlog(io),*)'  Windfile ',i,'of',iwfiles
                write(errlog(io),*)'  Offending line:',linebuffer130
              endif;enddo
              stop 1
            endif
            read(linebuffer130,'(a130)',err=9501,iostat=iostatus,iomsg=iomessage) MR_WindFiles(i)
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),1034) i,trim(adjustl(MR_WindFiles(i)))
            endif;enddo
            if(idf.eq.3)then
              ! If we are reading grib files, check that the index file has been
              ! generated
#ifndef USEGRIB
             do io=1,2;if(VB(io).le.verbosity_error)then
                write(errlog(io),*)"ERROR: ",&
                      "This input file specifies that the Met files are"
                write(errlog(io),*)&
                      "       in grib format, but Ash3d has not been compiled"
                write(errlog(io),*)&
                      "       with grib support.  Please recompile with grib"
                write(errlog(io),*)&
                      "       enabled.  MetReader must also support grib."
              endif;enddo
              stop 1
#else
              MR_WindFiles_GRIB_index(i) = trim(adjustl(MR_WindFiles(i))) // ".index"
              inquire( file=MR_WindFiles_GRIB_index(i), exist=IsThere )
              MR_WindFiles_Have_GRIB_index(i) = IsThere
              if(.not.IsThere)then
                ! Grib index file is not there, Try to generate it.
                ! Note, we might have some permission problems here and we should
                ! set up a fail-safe to the cwd or something.
                do io=1,2;if(VB(io).le.verbosity_info)then
                  write(outlog(io),*)" Grib index file not found; attempting to create it."
                endif;enddo
                call MR_Set_Gen_Index_GRIB(MR_WindFiles(i))
              endif
#endif
            endif
            read(fid_ctrlfile,'(a130)',iostat=iostatus,iomsg=iomessage)linebuffer130
            linebuffer050 = "Reading next windfile of blk5"
            if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer130(1:80),iomessage)
          enddo
        endif
1034    format(' i=',i3,'  MR_WindFiles(i) = ',a)
      else
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "MR_iWindFiles = 0"
          write(errlog(io),*)&
                "       Either the number of windfiles specified = 0, or"
          write(errlog(io),*)&
                "       MR_Allocate_FullMetFileList has not been called."
        endif;enddo
        stop 1
      endif

      ! Error trap if more windfiles are entered than are specified
      if(linebuffer130(1:5).ne.'*****') then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)&
                'The beginning of the line following the list of', &
                ' windfiles did not'
          write(errlog(io),*)&
                'start with ''*****''.  Did you enter the correct', &
                '  number of windfiles?'
          write(errlog(io),*) 'Program stopped.'
        endif;enddo
        stop 1
      endif
      ! END OF BLOCK 5
      !************************************************************************

        ! Check for existance and compatibility with simulation time requirements
      call MR_Read_Met_DimVars(e_iyear(1))
      if(MR_BaseYear.ne.BaseYear)then
        ! Base year was reset, probably because a windfile had an old base year
        useLeap  = MR_useLeap
        BaseYear = MR_BaseYear
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Change in calandar; resetting e_StartTime"
        endif;enddo
        tmp_dp = HS_hours_since_baseyear(e_iyear(1),e_imonth(1),  &
                         e_iday(1),e_hour(1),BaseYear,useLeap)
        tmp_dp = tmp_dp - SimStartHour   ! Recast tmp_dp as the difference in calandars
        SimStartHour = SimStartHour + tmp_dp
        xmlSimStartTime = HS_xmltime(SimStartHour,BaseYear,useLeap)
      endif

        ! Now that we have the actual times available from the Met files, we can reset
        ! the Simulation Start times for forecast runs
      if(runAsForecast)then
        MR_Comp_StartHour = MR_windfile_starthour(1) + MR_windfile_stephour(1,1) + FC_Offset
        SimStartHour      = MR_Comp_StartHour
        xmlSimStartTime   = HS_xmltime(SimStartHour,BaseYear,useLeap)
      endif

      !************************************************************************
      ! BLOCK 6: AIRPORT/POI FILE
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading ctr file, past Blk 5, looking for Blk 6."
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk6)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (.not.IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Expecting a comment line separating blocks."
          write(errlog(io),*)'       Check that Block 5 is correct.'
        endif;enddo
        stop 1
      endif
      do while (IsComment)
         ! Line is a comment, read next line
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        linebuffer050 = "Reading test line from blk6"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
        linebuffer050 = "Reading testkey from linebuffer (Blk6)"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      enddo
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)' *******************************************'
        write(outlog(io),*)' Reading Block 6: Airport/POI location output'
        write(outlog(io),*)' *******************************************'
      endif;enddo
      ! Block 6 Line 1
      ! Read whether to write out ASCII airport file
      cdf_b6l1 = linebuffer080
      read(linebuffer080,'(a3)',err=9601,iostat=iostatus,iomsg=iomessage) answer
      linebuffer050 = "Reading control file, blk6 line 1"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      if(adjustl(trim(answer)).eq.'yes') then
        WriteAirportFile_ASCII = .true.
        Write_PT_Data          = .true.
      else if(adjustl(trim(answer(1:2))).eq.'no') then
        WriteAirportFile_ASCII = .false.
      else
        goto 9601
      endif

      ! Block 6 Line 2
      ! Read whether to write out grain-size distribution to airport file
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage) linebuffer080
      linebuffer050 = "Reading control file, blk6 line 2"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      ! Block 6 Line 2
      cdf_b6l2 = linebuffer080
      read(linebuffer080,'(a3)',err=9602,iostat=iostatus,iomsg=iomessage) answer
      linebuffer050 = "Reading answer from linebuffer"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      if(adjustl(trim(answer)).eq.'yes') then
        WriteGSD = .true.
      else if(adjustl(trim(answer(1:2))).eq.'no') then
        WriteGSD = .false.
      else
        goto 9602
      endif

      ! Block 6 Line 3
      ! Read whether to write out kml airport file
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage) linebuffer080
      linebuffer050 = "Reading control file, blk6 line 3"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      cdf_b6l3 = linebuffer080
      read(linebuffer080,'(a3)',err=9603,iostat=iostatus,iomsg=iomessage) answer
      linebuffer050 = "Reading answer from linebuffer"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      if(adjustl(trim(answer)).eq.'yes') then
        WriteAirportFile_KML = .true.
        Write_PT_Data        = .true.
      else if(adjustl(trim(answer(1:2))).eq.'no') then
        WriteAirportFile_KML = .false.
      else
        goto 9603
      endif
            
      ! Block 6 Line 4
      ! Read name of input file containing airport locations
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage) linebuffer080
      linebuffer050 = "Reading control file, blk6 line 4"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      cdf_b6l4 = linebuffer080
      AirportInFile = cdf_b6l4(1:scan(cdf_b6l4,' ')-1)     !Read to the first blank space

      ! See if we need to read an external airport file
      ReadExtAirportFile   = .false.
      AppendExtAirportFile = .false.
      if((AirportInFile.ne.'internal').and. &
          (AirportInFile.ne.'')) then
        ReadExtAirportFile=.true.              ! read external data
        if(AirportInFile(1:1).eq.'+') then
          AppendExtAirportFile=.true.          ! read and append external data to master list
          AirportInFile = AirportInFile(2:)    ! strip off the "plus" at the beginning
        else
          AppendExtAirportFile=.false.         ! read external data; do not append to master list
        endif
        ! Make sure the external file exists and can be opened.
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)'Making sure the external airport file exists and can be opened'
          write(outlog(io),*)'Opening ',trim(adjustl(AirportInFile))
        endif;enddo
        ! try opening the external file
        open(unit=fid_airport,file=AirportInFile,status='old',action='read',err=9604)
        close(fid_airport)                         ! if it opens, close it back up.
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)'Success.'
        endif;enddo
      else
         ReadExtAirportFile = .false.              ! Do not read external data
                                                   ! Internal list might still be used
      endif

      ! Read whether to project airport coordinates
96042 continue

      ! Block 6 Line 5
      ! Have libprojection calculate projected coordinates?
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage) linebuffer080
      linebuffer050 = "Reading control file, blk6 line 5"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      cdf_b6l5 = linebuffer080
      read(linebuffer080,'(a3)',err=9605,iostat=iostatus,iomsg=iomessage) answer
      linebuffer050 = "Reading answer from linebuffer"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      if(adjustl(trim(answer)).eq.'yes') then
        ProjectAirportLocations = .true.
      else if(adjustl(trim(answer(1:2))).eq.'no') then
        ProjectAirportLocations = .false.
      else
        goto 9605
      endif

      ! Write out parameters
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),44) ReadExtAirportFile, AppendExtAirportFile, &
                  WriteAirportFile_ASCII, WriteGSD, WriteAirportFile_KML, &
                  ProjectAirportLocations
44      format(4x,'Airport choices:',/, &
                4x,'   Read external file of locations (T/F)                        = ',L1,/, &
                4x,'   Append external locations to airport list (T/F)              = ',L1,/, &
                4x,'   Write ASCII file of ash arrival times at airports (T/F)      = ',L1,/, &
                4x,'   Write deposit grain-size distribution to airports file (T/F) = ',L1,/, &
                4x,'   Write KML file of ash arrival times at airports (T/F)        = ',L1,/, &
                4x,'   Calculate projected airport locations using Ash3d (T/F)      = ',L1,/)
        if(ReadExtAirportFile)then
          write(outlog(io),*)'   Name of file containing airport locations: ', &
                             adjustl(trim(AirportInFile))
          if(AppendExtAirportFile)then
            write(outlog(io),*)'   Airports/POI in this file will be appended to the internal'
            write(outlog(io),*)'   global list.'
          else
            write(outlog(io),*)'   This file will replace the internal global list of airports.'
          endif
        else
          write(outlog(io),*)'   Internal list of global airport locations will be used.'
        endif
          write(outlog(io),*)' '
      endif;enddo
      if(ProjectAirportLocations)then
        if(IsLatLon)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"     The control file indicates that the projected coordinates"
            write(outlog(io),*)"     in columns 3 and 4 of the Airport/POI file should be used,"
            write(outlog(io),*)"     but the current coordinate system is lon/lat.  Projected"
            write(outlog(io),*)"     coordinates will be ignored and the lon/lat provided in"
            write(outlog(io),*)"     the file used instead."
          endif;enddo
        else
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"     The control file indicates that the projected coordinates"
            write(outlog(io),*)"     of the Airport/POI file should be used.  Please make sure"
            write(outlog(io),*)"     that the projected coordinates match the projection of"
            write(outlog(io),*)"     the computational grid."
          endif;enddo
        endif
      endif

      ! END OF BLOCK 6
      !************************************************************************

      !************************************************************************
      ! BLOCK 7: GRAIN-SIZE BINS, SETTLING VELOCITY
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading ctr file, past Blk 6, looking for Blk 7"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk7)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (.not.IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Expecting a comment line separating blocks."
          write(errlog(io),*)'       Check that Block 6 is correct.'
        endif;enddo
        stop 1
      endif
      do while (IsComment)
         ! Line is a comment, read next line
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        linebuffer050 = "Reading test line from control file (Blk7)"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
        linebuffer050 = "Reading testkey from linebuffer (Blk7)"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      enddo
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)' *******************************************'
        write(outlog(io),*)' Reading Block 7: Grain Size Groups'
        write(outlog(io),*)' *******************************************'
      endif;enddo
      ! READ GRAIN-SIZE BINS
      ! First, get the number of tephra bins to read
      ! Note: This might be one greater than what is calculated if the last bin
      !       has a negative diameter.  In this case, the remaining mass fraction
      !       neglecting the last bin is distributed over the previous bins with
      !       a gaussian distribution given by a LN_phi_mean and LN_phi_stdev
      !    e.g.  -1 4 2
      !       Also note that the number of tephra bins can be zero if the species
      !       will be defined in optional modules such as gas, aggregates, etc.
      ! Block 7 Line 1
      ! We at least need the number of tephra bins
      read(linebuffer080,*,err=9701,iostat=iostatus,iomsg=iomessage) ivalue1
      init_n_gs_max = ivalue1
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) ivalue1, ivalue2
      FV_ID    = 1
      Shape_ID = 1
      ! Assume we can read at least read one value, try for two with the second being
      ! the fall model:
      !  FV_ID = 0 -> No fall (just tracer)
      !          1 -> Wilson and Huang
      !          2 -> Wilson and Huang + Cunningham slip
      !          3 -> Wilson and Huang + Mod by Pfeiffer Et al.
      !          4 -> Ganser
      !          5 -> Ganser + Cunningham slip
      !          6 -> Stokes flow for spherical particles + slip
      if(iostatus.eq.0)then
        FV_ID = ivalue2
        if(FV_ID.ne.0.and.FV_ID.ne.1.and.FV_ID.ne.2.and.&
           FV_ID.ne.3.and.FV_ID.ne.4.and.FV_ID.ne.5.and.FV_ID.ne.6)then
          FV_ID = 1 ! Default to Wilson and Huang
        endif
        ! Try for a third value which specifies shape factory (F vs phi)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) ivalue1, ivalue2, ivalue3
        if(iostatus.eq.0)then
          Shape_ID = ivalue3
          if(Shape_ID.ne.1.and.Shape_ID.ne.2)then
            Shape_ID = 1 ! Default to Wilson and Huang
          endif
        endif
      else
        FV_ID = 1 ! Wilson and Huang
      endif

      allocate(temp_v_s(init_n_gs_max))
      allocate(temp_gsdiam(init_n_gs_max))
      allocate(temp_bin_mass(init_n_gs_max))
      allocate(temp_rho_m(init_n_gs_max))
      allocate(temp_gsF(init_n_gs_max))
      allocate(temp_gsG(init_n_gs_max))

      if(init_n_gs_max.lt.0)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*) 'ERROR: Number of grainsizes must be non-negative.'
          write(errlog(io),*) 'Program stopped'
        endif;enddo
        stop 1
      elseif(init_n_gs_max.eq.0)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(errlog(io),*) 'WARNING: Number of grainsizes is 0.'
          write(errlog(io),*) '         Assuming a custom module is being used to'
          write(errlog(io),*) '         add to the concen array.'
        endif;enddo
        n_gs_max = init_n_gs_max
      else
        ! This is the normal case with actual grain size bins specified
        LN_phi_mean   = 0.0_ip
        LN_phi_stddev = 0.0_ip
        LN_massfrac   = 0.0_ip
        do isize=1,init_n_gs_max
          value1 = -1.99_ip
          value2 = -1.99_ip
          value3 = -1.99_ip
          read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
          linebuffer050 = "Reading control file, blk7, line 2+"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
          ! Always check if we have overshot the block
          testkey = linebuffer080(1:1)
          call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
          if (IsComment) then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: ",&
                    "Error in specifying grain sizes.  You specified ",&
                            init_n_gs_max,', sizes,'
              write(errlog(io),*) 'but only ',isize-1,&
                    ' size classes were listed in the input file.'
              write(errlog(io),*)' Offending line: ',linebuffer080
              write(errlog(io),*) 'Program stopped'
            endif;enddo
            stop 1
          endif
          ! Read at least two values or throw an error
          read(linebuffer080,*,err=9702,iostat=iostatus,iomsg=iomessage) value1, value2
          if(isize.eq.1) Tephra_Ncols = 2
          ! Two values successfully read; try for three
          read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) value1, value2, value3
          if(iostatus.eq.0)then
            ! Three values were successfully read, interpret as:
            ! grain-size, mass fraction, density
            ! W&H suggest 800 kg/m3 for d>300um and 2000 for d<88um for pumice
            ! fragments
            if(isize.eq.1) Tephra_Ncols = 3
            useCalcFallVel = .true. 
            useTemperature = .true. ! When calculating Fall Vel. we need T
            temp_gsdiam(isize) = value1
            temp_bin_mass(isize) = value2
            temp_rho_m(isize) = value3
            ! Try for a forth value for shape
            read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) value1, value2, value3, value4
            if(iostatus.eq.0)then
              ! Fourth value was successfully read, interpret as W/H shape
              ! parameter
              if(isize.eq.1) Tephra_Ncols = 4
              temp_gsF(isize) = value4
              ! Try for a fifth value for ratio of minor axies of ellipsoid
              read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) value1, value2, value3, value4, value5
              if(iostatus.eq.0)then
                ! Fifth value was successfully read, interpret as c/b
                if(isize.eq.1) Tephra_Ncols = 5
                temp_gsG(isize) = value5
              else
                temp_gsG(isize) = 1.0_ip
              endif
            else
              temp_gsF(isize) = 0.44_ip
              temp_gsG(isize) = 1.0_ip
            endif
              ! Initialize this to zero
            temp_v_s(isize) = 0.0_ip
            if(temp_gsdiam(isize).lt.0.0_ip)then
              if(isize.lt.init_n_gs_max)then
                do io=1,2;if(VB(io).le.verbosity_error)then
                  write(errlog(io),*)"ERROR: ",&
                        "diameter must be positive",isize,init_n_gs_max,temp_gsdiam(isize)
                endif;enddo
                stop 1
              else
                LN_phi_mean   = value2
                LN_phi_stddev = value3
                LN_massfrac   = 1.0_ip-sum(temp_bin_mass(1:init_n_gs_max-1))
                do io=1,2;if(VB(io).le.verbosity_info)then
                  write(outlog(io),*) &
                        "Last grain-size bin will be partitioned across all previous."
                  write(outlog(io),*)"Volume fraction partitioned = ",&
                        LN_massfrac
                  write(outlog(io),*)&
                        "  Assuming remainder is Gaussian in phi"
                  write(outlog(io),*)&
                        "    LN_phi_mean   = ", LN_phi_mean
                  write(outlog(io),*)&
                        "    LN_phi_stddev = ", LN_phi_stddev
                endif;enddo
                useLogNormGSbins = .true.
              endif
            endif
          else
            ! Only two values were successfully read, interpret with
            ! old format as:
            ! FallVel, mass fraction
            useCalcFallVel = .false.
            temp_v_s(isize)     = value1
            if(temp_v_s(isize).lt.0.0_ip)then
              do io=1,2;if(VB(io).le.verbosity_info)then
                write(outlog(io),*)&
                   "WARNING: fall velocity is negative.  Grains will 'fall' upward"
              endif;enddo
            endif
            temp_bin_mass(isize)= value2
              ! Initialize these
            temp_gsdiam(isize)  = 0.1_ip
            temp_rho_m(isize)   = 2000.0_ip
            temp_gsF(isize)     = 0.44_ip
            temp_gsG(isize)     = 1.0_ip
          endif
        enddo ! isize=1,init_n_gs_max
        ! Set the number of grain-size bins
        if(useLogNormGSbins)then
          ! In this case, the last bin is the remainder to be distributed over the
          ! previous bins.  So we need to decrement the tephra bins by 1
          n_gs_max = init_n_gs_max-1
        else
          n_gs_max = init_n_gs_max
        endif
      endif ! if init_n_gs_max > 0
      ! END OF BLOCK 7
      !************************************************************************

      ! Since this subroutine is called before any optional modules, we can
      ! initialize nsmax to the number of tephra bins
      nsmax      = n_gs_max  ! Total tracked bins
      n_gs_aloft = n_gs_max  ! Number of tephra species aloft

      if(n_gs_max.gt.0)then
        call Allocate_Tephra
        allocate(temp_phi(n_gs_max))

        Tephra_v_s(1:n_gs_max)      = -1.0_ip * temp_v_s(1:n_gs_max) ! make sure 'fall velocity'
                                                                     ! is in the -z direction
        Tephra_gsdiam(1:n_gs_max)   = temp_gsdiam(1:n_gs_max)
        Tephra_bin_mass(1:n_gs_max) = temp_bin_mass(1:n_gs_max)
        Tephra_rho_m(1:n_gs_max)    = temp_rho_m(1:n_gs_max)
        if(Shape_ID.eq.1)then
          ! Interpret shape columns as F [and G]
          Tephra_gsF(1:n_gs_max)      = temp_gsF(1:n_gs_max)
          Tephra_gsG(1:n_gs_max)      = temp_gsG(1:n_gs_max)
          Tephra_gsPhi(1:n_gs_max)    = 1.0_ip
        elseif(Shape_ID.eq.2)then
          ! Interpret shape columns as Phi
          Tephra_gsF(1:n_gs_max)      = 1.0_ip
          Tephra_gsG(1:n_gs_max)      = 1.0_ip
          Tephra_gsPhi(1:n_gs_max)    = temp_gsF(1:n_gs_max)
        endif
      endif

      deallocate(temp_v_s,temp_gsdiam,temp_bin_mass,temp_rho_m,temp_gsF)

      if(n_gs_max.gt.0)then
        call Calculate_Tephra_Shape

        ! If a log-normal distribution is to be added, make sure the grainsize
        ! bins are sorted by size (smallest first)
        if(useLogNormGSbins) call Sort_Tephra_Size
        temp_phi = -log(Tephra_gsdiam)/log(2.0)

        if(useCalcFallVel) Tephra_gsdiam = Tephra_gsdiam/1000.0_ip   ! convert diameter from mm to m

        ! Find the fraction of fine (<= phi4, 63um)
        fracfine = 0.0_ip
        do isize=1,n_gs_max
          if(Tephra_gsdiam(isize).lt.6.4e-5_ip)then
            fracfine = fracfine + Tephra_bin_mass(isize)
          endif
        enddo

!       Make sure that Tephra_bin_mass>0 for  each size class
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),2531)
        endif;enddo
2531    format('Checking to make sure Tephra_bin_mass>0 for all size classes')
        do isize=1,n_gs_max
          if(Tephra_bin_mass(isize).lt.0.0_ip) then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),25311) isize
            endif;enddo
            stop 1
          endif
        enddo
25311        format('Error: mass of bin ',i2,' is less than zero.  Program stopped')

!       Send error message if sum(bin_mass(1:n_gs_max)) does not equal 1
        sum_bins=sum(Tephra_bin_mass(1:n_gs_max))
        if(abs(sum_bins-1.0_ip).gt.0.02_ip) then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),2532)
            do isize=1,n_gs_max
              write(errlog(io),2533) isize, Tephra_bin_mass(isize)
            enddo
            write(errlog(io),2534) sum_bins
          endif;enddo
2532      format('Error.  Sum of mass fractions of grain-size bins',/, &
                 'does not equal 1.',/, &
                 'bin   mass')
2533      format(i3,f7.4)
2534      format(3x,f7.4,'  total',/,'Program stopped')
          stop 1
            ! If it differs just slightly from 1, adjust automatically
        else if(abs(sum_bins-1.0_ip).gt.0.001_ip) then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),2535) sum_bins
          endif;enddo
2535      format('Warning: sum(bin_mass(1:n_gs_max))=',f10.5,/, &
                 'This differs slightly from 1.0',/, &
                 'adjusting bin masses automatically.')
          do isize=1,n_gs_max
            Tephra_bin_mass(isize) = Tephra_bin_mass(isize)*1.0_ip/sum_bins
          enddo
        endif

        ! Write out grain-size bins
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),9) n_gs_max, FV_ID
          if(FV_ID.eq.0)then
            write(outlog(io),*)"               None (tracer)"
          elseif(FV_ID.eq.1)then
            write(outlog(io),*)"               Wilson and Huang"
            write(outlog(io),*)"     Input particle sizes are interpreted to be the mean; Da=(A+B+C)/3"
          elseif(FV_ID.eq.2)then
            write(outlog(io),*)"               Wilson and Huang + Cunningham slip"
            write(outlog(io),*)"     Input particle sizes are interpreted to be the mean; Da=(A+B+C)/3"
          elseif(FV_ID.eq.3)then
            write(outlog(io),*)"               Wilson and Huang + Mod by PCM"
            write(outlog(io),*)"     Input particle sizes are interpreted to be the mean; Da=(A+B+C)/3"
          elseif(FV_ID.eq.4)then
            write(outlog(io),*)"               Ganser"
            write(outlog(io),*)"     Input particle sizes are interpreted to be the geometric mean; Dv=(ABC)^0.33"
            write(outlog(io),*)"     which is also the diameter of a volume-equivalent sphere."
          elseif(FV_ID.eq.5)then
            write(outlog(io),*)"               Ganser + Cunningham slip"
            write(outlog(io),*)"     Input particle sizes are interpreted to be the geometric mean; Dv=(ABC)^0.33"
            write(outlog(io),*)"     which is also the diameter of a volume-equivalent sphere."
          elseif(FV_ID.eq.6)then
            write(outlog(io),*)"               Stokes flow + slip"
          else
            write(outlog(io),*)"Default Fall Model = Wilson and Huang"
            write(outlog(io),*)"     Input particle sizes are interpreted to be the mean; Da=(A+B+C)/3"
          endif
          if(Shape_ID.eq.1)then
            write(outlog(io),*)"Shape Specification : axis ratios : F=(b+c)/2a [and G=c/b]"
          elseif(Shape_ID.eq.2)then
            write(outlog(io),*)"Shape Specification : sphericity : Area-of-vol.eq.sphere/Area-of-particle"
          endif
        endif;enddo
        if(useCalcFallVel)then
          ! Note: F and G are currently not calculated from sphericity
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),10)
            do isize=1,n_gs_max
              ! write out diameter in mm, not m
              write(outlog(io),11) Tephra_bin_mass(isize), Tephra_gsdiam(isize)*1000.0_ip, &
                                   Tephra_rho_m(isize),    Tephra_gsF(isize),              &
                                   Tephra_gsG(isize),      Tephra_gsPhi(isize), temp_phi(isize)
            enddo
          endif;enddo
          if(FV_ID.gt.0.and.Shape_ID.eq.1)then
            ! Now write some details of the particle axes given Da and F (and G)
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"Particle shapes were defined with F and G.  This corresponds to the"
              if(FV_ID.le.3)then
                write(outlog(io),*)"following axes lengths (Assuming diameter = arithmetic average):"
              elseif(FV_ID.le.5)then
                write(outlog(io),*)"following axes lengths (Assuming diameter = geometric average):"
              else
                write(outlog(io),*)"following axes lengths:"
              endif
              write(outlog(io),*)"   Bin # : diameter (mm)   :  A (mm)   :  B (mm)   :  C (mm)"
              do isize=1,n_gs_max
                ! write out diameter in mm, not m
                if(FV_ID.le.3)then
                  ! Wilson/Huang models uses arithmetic average
                  Davg = Tephra_gsdiam(isize)*1000.0_ip
                  Aaxis= Davg*(3.0_ip/(1.0_ip+2.0_ip*Tephra_gsF(isize)))
                  Baxis= Davg*6.0_ip*Tephra_gsF(isize)/&
                   ( (1.0_ip+Tephra_gsG(isize)) * (1.0_ip+2.0_ip*Tephra_gsF(isize)))
                  Caxis=3.0_ip*Davg - Aaxis - Baxis
                elseif(FV_ID.le.5)then
                  ! Ganser model uses geometric average
                  Davg = Tephra_gsdiam(isize)*1000.0_ip
                  Aaxis= Davg*(((1.0_ip+Tephra_gsG(isize))**2.0_ip)/&
                                (4.0_ip*Tephra_gsF(isize)*Tephra_gsF(isize)*Tephra_gsG(isize)))**(1.0_ip/3.0_ip)
                  Baxis= sqrt(Davg*Davg*Davg/Aaxis/Tephra_gsG(isize))
                  Caxis= Davg**3.0_ip/Aaxis/Baxis
                else
                  ! Stokes flow is for a spherical particle
                  Davg = Tephra_gsdiam(isize)*1000.0_ip
                  Aaxis= Davg
                  Baxis= Davg
                  Caxis= Davg
                endif
                write(outlog(io),'(4x,i5,4x,4f12.4)')isize,Davg,Aaxis,Baxis,Caxis
              enddo
            endif;enddo
          endif
        else
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),2110)
            do isize=1,n_gs_max
              write(outlog(io),2111) Tephra_bin_mass(isize), Tephra_v_s(isize)
            enddo
          endif;enddo
        endif
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)
          write(outlog(io),'(a26,f10.3,a6)')"Using deposit density of: ",DepositDensity," kg/m3"
          write(outlog(io),*)
          write(outlog(io),'(a45,f10.3)')"Using a mass-fraction of fines (< 63um) of : ",fracfine
        endif;enddo
        ! if bin masses sum up close to 1 (within 1%), adjust automatically      
        if((abs(sum(Tephra_bin_mass)-1.0_ip).gt.1.0e-5_ip).and. &  
            (abs(sum(Tephra_bin_mass)-1.0_ip).le.1.0e-02_ip)) then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),1001) sum(Tephra_bin_mass)
          endif;enddo
1001      format(4x,'The sum of the bin masses=',f5.3,&
                 '  Adjusting to equal 1',/, &
                 4x,'New bin masses:')
          Tephra_bin_mass(1:n_gs_max) = Tephra_bin_mass(1:n_gs_max)*1.0_ip/sum(Tephra_bin_mass(1:n_gs_max))

          if(useCalcFallVel)then
            do io=1,2;if(VB(io).le.verbosity_info)then
              do isize=1,n_gs_max
                ! write out diameter in mm, not m
                write(outlog(io),11) Tephra_bin_mass(isize), Tephra_gsdiam(isize)*1000.0_ip, &
                                     Tephra_rho_m(isize)
              enddo
            endif;enddo
          else
            do io=1,2;if(VB(io).le.verbosity_info)then
              do isize=1,n_gs_max
                write(outlog(io),2111) Tephra_bin_mass(isize), Tephra_v_s(isize)
              enddo
            endif;enddo
          endif
        else if(abs(sum(Tephra_bin_mass)-1.0_ip).gt.1.0e-2_ip) then
          ! if bin masses are do not sum to within 1% of unity; stop
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),1000)
          endif;enddo
1000      format(4x,'Error: Sum of the mass fraction of the grain',&
                 ' size bins differs from 1 by more than 0.01',/, &
                 4x,'Program stopped')
          stop 1
        endif
      endif

      !************************************************************************
      ! BLOCK 8: VERTICAL PROFILES
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading ctr file, past Blk 7, looking for Blk 8"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk8)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (.not.IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Expecting a comment line separating blocks."
          write(errlog(io),*)'       Check that Block 7 is correct.'
        endif;enddo
        stop 1
      endif
      do while (IsComment)
         ! Line is a comment, read next line
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        linebuffer050 = "Reading test line from control file (Blk8)"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
        linebuffer050 = "Reading testkey from linebuffer (Blk8)"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      enddo
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)' *******************************************'
        write(outlog(io),*)' Reading Block 8: Vertical profile output'
        write(outlog(io),*)' *******************************************'
      endif;enddo

      ! Block 8 Line 1
      ! Read number of vertical profiles
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*) 'Reading vertical profile information'
      endif;enddo
      read(linebuffer080,*,err=9801,iostat=iostatus,iomsg=iomessage) nvprofiles
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*) 'number of vertical profiles=',nvprofiles
      endif;enddo

      if(nvprofiles.gt.0) then
        Write_PR_Data = .true.
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Allocating profile arrays:",nvprofiles
        endif;enddo
        allocate(x_vprofile(nvprofiles))
        allocate(y_vprofile(nvprofiles))
        allocate(i_vprofile(nvprofiles))
        allocate(j_vprofile(nvprofiles))
        allocate(Site_vprofile(nvprofiles))
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),46)
        endif;enddo
46      format(/,'     vertical profile locations',/, &
                  '          #         x         y     i     j')
        do i=1,nvprofiles
          read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage) linebuffer080
          linebuffer050 = "Reading Blk8, line 1+ (vprof)"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
          ! Always check if we have overshot the block
          testkey=linebuffer080(1:1)
          call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
          if (IsComment) then
          do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: ",&
                    "Trying to read Block 8 and detecting comment line"
              write(errlog(io),*)'  Vert Prof ',i,'of',nvprofiles
              write(errlog(io),*)'  Offending line:',linebuffer080
            endif;enddo
            stop 1
          endif

          ! Block 8 Line 2+
          read(linebuffer080,*,err=9802,iostat=iostatus,iomsg=iomessage) value1, value2
          x_vprofile(i) = value1
          y_vprofile(i) = value2
          write(Site_vprofile(i),'(a14,1x,i3)')"Vertical Prof ",i
          ! Assume we can read at least read two values, try for three
          read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) value1, value2, Site_vprofile(i)
          if(iostatus.eq.0)then
            substr_pos1 = index(linebuffer080,trim(adjustl(Site_vprofile(i))))
            substr_pos2 = index(linebuffer080,'#')
            if(substr_pos2.eq.0)then
              ! comment indicator '#' not found, set end of string to length
              substr_pos2 = min(len(linebuffer080),substr_pos1+50)
            endif
            Site_vprofile(i) = trim(adjustl(linebuffer080(substr_pos1:substr_pos2)))
          endif
          call vprofchecker(i)
        enddo
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*) 'No vertical profile specified'
        endif;enddo
      endif
      ! END OF BLOCK 8
      !************************************************************************

      !************************************************************************
      ! BLOCK 9: NETCDF ANNOTATIONS
      !  This block is only optional in the sense that variables will have default
      !  values if the input file ends before this block.  However, the presence of
      !  optional_module blocks will not read correctly.
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading ctr file, past Blk 8, looking for Blk 9"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
      linebuffer050 = "Reading testkey from linebuffer (Blk9)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      if (.not.IsComment) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: ",&
                "Expecting a comment line separating blocks."
          write(errlog(io),*)'       Check that Block 8 is correct.'
        endif;enddo
        stop 1
      endif      
      do while(iostatus.eq.0.and.IsComment)
         ! Line is a comment, read next line
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        linebuffer050 = "Reading test line from control file blk9"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)testkey
        linebuffer050 = "Reading testkey from linebuffer (Blk9)"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        call FileIO_Check_testkey(testkey,linebuffer080,IsComment)
      enddo

      ! Here are the default output file name and comments if Block 9 is not given
      concenfile  = "3d_tephra_fall.nc"
      cdf_title   = infile
      cdf_comment = ""
      if(iostatus.ne.0)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)'  Setting outfile to: 3d_tephra_fall.nc'
          write(outlog(io),*)'  Setting Title to: ',infile
        endif;enddo
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)' *****************************************'
          write(outlog(io),*)' Reading Block 9: Output file / comments'
          write(outlog(io),*)' *****************************************'
        endif;enddo
        ! Start reading annotation info

        ! First line is the output file name
        Have_Block_NetCDF = .true.
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) concenfile
        linebuffer050 = "Reading control file, blk9, line 1"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        concenfile = trim(adjustl(concenfile))

        ! We only really need the name of the output concentration file. If we don't
        ! have the next two lines, error in reading will sent control to label 2010,
        ! which just allows continuation of code.

        ! Next line is the title of the job
          ! Read title line up until the first '#', then truncate
        read(fid_ctrlfile,'(a80)',err=2010,iostat=iostatus,iomsg=iomessage)linebuffer080
        iendstr = scan(linebuffer080, "#")
        if(iendstr.eq.0)then
             ! '#' not found, just copy linebuffer080 to title
          cdf_title = trim(adjustl(linebuffer080))
        else
            ! clip title at key
          cdf_title = trim(linebuffer080(1:iendstr-1))
        endif

          ! Read comment line up until the first '#', then truncate
        read(fid_ctrlfile,'(a80)',err=2010,iostat=iostatus,iomsg=iomessage)linebuffer080
        if(iostatus.ne.0)goto 2010
        iendstr = scan(linebuffer080, "#")
        if(iendstr.eq.0)then
             ! '#' not found, just copy linebuffer080 to comment
          cdf_comment = trim(adjustl(linebuffer080))
        else
            ! clip comment at key
          cdf_comment = trim(linebuffer080(1:iendstr-1))
        endif

        ! This is the end of the standard blocks
      endif
      ! END OF BLOCK 9
      !************************************************************************

      !************************************************************************
      ! Might have optional modules; these were identified at the start of this
      ! subroutine.
      !************************************************************************

     ! close input file 
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),25) infile
      endif;enddo
      close(fid_ctrlfile)
      ! Here is the end of the input file

      ! Now write out what we found
2010  continue
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)' *****************************************'
        write(outlog(io),*)' Finished reading/parsing input file      '
        write(outlog(io),*)' Custom blocks for optional modules will  '
        write(outlog(io),*)' be read in corresponding input_data_OPTMOD'
        write(outlog(io),*)' subroutines provided by the modules.     '
        write(outlog(io),*)' *****************************************'
      endif;enddo

      do io=1,2;if(VB(io).le.verbosity_info)then
        ! Write out computational scheme used
        if(useDS)then
          write(outlog(io),*)"Dimension splitting will be used."
        else
          ! If something else besides dimension splitting is used (CTU, Semi-Lagr.)
          !  write out a note about it here.
          !write(outlog(io),*)""
        endif
        write(outlog(io),*)limiter," limiter is used."
        if(useCN) then
          write(outlog(io),*)&
           "Diffusion is calculated implicitly using lapack routines for solving Ax=b."
          write(outlog(io),*)&
           "Note, Imp_fac controls the amount of the t+1 step that is used in the stencil."
          write(outlog(io),*)&
           "  Imp_fac = 1.0 (  0% t, 100% t+1) :: Backward Euler"
          write(outlog(io),*)&
           "  Imp_fac = 0.5 ( 50% t,  50% t+1) :: Crank-Nicolson"
          write(outlog(io),*)&
           "  Imp_fac = 0.0 (100% t,   0% t+1) :: Forward Euler"
          write(outlog(io),*)&
           "  In this run, we are using Imp_fac = ",real(Imp_fac,kind=4)
          IsSat = .false.
          if(Imp_fac.lt.0.0_ip)then
            ! error
            write(outlog(io),*)"    ERROR: Imp_fac must be > 0"
            stop 1
          elseif(Imp_fac.lt.EPS_SMALL)then
            ! FE
            write(outlog(io),*)"    Forward Euler"
            write(outlog(io),*)"     CFL condition requires that Imp_DT_fac<=0.5"
            if(Imp_DT_fac.le.0.5_ip)IsSat = .true.
          elseif(abs(Imp_fac-0.5_ip).lt.EPS_SMALL)then
            ! CN
            write(outlog(io),*)"    Crank-Nicolson"
            write(outlog(io),*)"     Unconditionally stable, but for accuracy, please use Imp_DT_fac<=4.0"
            if(Imp_DT_fac.le.4.0_ip)IsSat = .true.
          elseif(abs(Imp_fac-1.0_ip).lt.EPS_SMALL)then
            ! BE
            write(outlog(io),*)"    Backward Euler"
            write(outlog(io),*)"     Unconditionally stable, but for accuracy, please use Imp_DT_fac<=4.0"
            if(Imp_DT_fac.le.4.0_ip)IsSat = .true.
          else
            ! Some other blend
            write(outlog(io),*)"    Non-standard mix of t and t+1"
            write(outlog(io),*)"     Unknown stability; probably should use Imp_DT_fac<=0.5"
            if(Imp_DT_fac.le.0.5_ip)IsSat = .true.
          endif
          if(IsSat)then
             write(outlog(io),*)"     Good, it satisfies the condition. Imp_DT_fac=",real(Imp_DT_fac,kind=4)
          else
            write(outlog(io),*)"     WARNING: Condition not satisfied. Imp_DT_fac=",real(Imp_DT_fac,kind=4)
          endif
        else
          write(outlog(io),*)"Diffusion is calculated explicitly."
        endif
   
        ! Write out Vz calculation scheme used
        if(useVz_rhoG)then
          write(outlog(io),*)"useVz_rhoG=.true. : Vz calculated PVV (if avail.) and density"
        else
          write(outlog(io),*)"Vz calculated via PVV and finite-differencing dp/dz"
        endif
      endif;enddo

      ! Calculate size of grid (cell-centered)
      if(IsLatLon) then
        nxmax = ceiling((lonUR-lonLL)/de)     ! number of x nodes
        nymax = ceiling((latUR-latLL)/dn)     ! number of y nodes
      else      
        nxmax = ceiling((xUR-xLL)/dx)         ! number of x nodes
        nymax = ceiling((yUR-yLL)/dy)         ! number of y nodes
      endif
      ! initialize active grid domain to the full grid
      imin = 1
      imax = nxmax
      jmin = 1
      jmax = nymax
      kmin = 1
      kmax = nzmax

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),16) nxmax, nymax, nzmax
        write(outlog(io),*)
      endif;enddo

      deallocate(e_iyear)
      deallocate(e_imonth)
      deallocate(e_iday)
      deallocate(e_hour)

      ! Set up logging logical values
        ! First check for output requests that require evaluating every time step
      if(Write_PT_Data.or.Write_PR_Data)then
        Output_every_TS = .true.
      else
        Output_every_TS = .false.
      endif
        ! Next check if there are requests at specific write times
      if(WriteDepositTS_ASCII           .or.      &
          WriteCloudConcentration_ASCII  .or.      &
          WriteCloudHeight_ASCII         .or.      &
          WriteCloudLoad_ASCII           .or.      &
          WriteDepositTS_KML             .or.      &
          WriteCloudConcentration_KML    .or.      &
          WriteCloudHeight_KML           .or.      &
          WriteCloudLoad_KML             .or.      &
          WriteReflectivity_KML          .or.      &
          Write3dFiles) then
        Output_at_WriteTimes = .true.
      else
        Output_at_WriteTimes = .false.
      endif
        ! Finally, check if we will be logging progress at log_steps
      if(log_step.gt.0)then
        Output_at_logsteps = .true.
      else
        Output_at_logsteps = .false.
      endif

      return

!******************************************************************************
      ! Error traps (starting with 9000)
      ! For this subroutine, the 100's position refers to block # of control file

9001  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error: cannot open input file: ',infile
        write(errlog(io),*)  'Program stopped'
      endif;enddo
      stop 1

      !Block 1/Line 3: x/y or lon/lat lower-left corner of grid
9103  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error reading xLL or yLL.'
        write(errlog(io),*)  'You entered: ',cdf_b1l3
        write(errlog(io),*)  'Program stopped.'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(1)
      stop 1

      !Block 1/Line 4: grid width/height
9104  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error reading width and height of model domain.'
        write(errlog(io),*)  'You entered: ', cdf_b1l4
        write(errlog(io),*)  'Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(1)
      stop 1

      !Block 1/Line 5: vent location
9105  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error reading vent coordinates.'
        write(errlog(io),*)  'You entered: ', cdf_b1l4
        write(errlog(io),*)  'Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(1)
      stop 1

      !Block 1/Line 6: dx and dy of computational grid
9106  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error reading dx or dy.'
        write(errlog(io),*)  'You entered: ', cdf_b1l6
        write(errlog(io),*)  'Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(1)
      stop 1

      !Block 1/Line 7: dz of computational grid
9107  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error reading dz or dz type.'
        write(errlog(io),*)  'You gave: ',cdf_b1l7
        write(errlog(io),*)  'Program stopped.'      
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(1)
      stop 1

      !Block 1/Line 8: Diffusivity and source type / Suzuki paramter
9108  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error reading diffusion coefficient or Suzuki constant or plume type.'
        write(errlog(io),*)  'The first value should be a number.  The second value should be either'
        write(errlog(io),*)  'a number (the Suzuki constant), or the word "line", "point",'
        write(errlog(io),*)  '"profile" or "umbrella" or "umbrella_air"'
        write(errlog(io),*)  'You entered: ',cdf_b1l8
        write(errlog(io),*)  'Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(1)
      stop 1

      !Block 1/Line 9: Number of eruptions
9109  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error reading number of eruptions.'
        write(errlog(io),*)  'You gave: ',cdf_b1l9
        write(errlog(io),*)  'Program stopped.'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(1)
      stop 1      

      !Block 2: ERUPTION PARAMETERS
9201  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error reading start time, duration, height or',&
                    ' volume of an eruptive pulse.  Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(2)
      stop 1 


9202  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),5678)  i,e_StartTime(i),e_StartTime(i-1)+e_Duration(i-1),e_hour(i)
      endif;enddo
5678  format(4x,'error: eruption pulses are not in chronological order.',/, &
             4x,'e_StartTime(i)<(e_StartTime(i-1)+e_Duration(i-1))',/, &
             4x,'                               i=',i3,/, &
             4x,'                  e_StartTime(i)=',e15.8,/, &
             4x,'e_StartTime(i-1)+e_Duration(i-1)=',e15.8,/, &
             4x,'                        hour(i) =',f12.4,/, &
             4x,'Program stopped')
      stop 1

      !Block 3/Line 1: wind file specifications 
9301  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error reading iwind, iwindformat.  Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(3)
      stop 1

      !Block 3/Line 1.1: template file name
93011 do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error reading template file name.  Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(3)
      stop 1

      !Block 3/Line 2: iHeightHandler
9302  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*) 'Error reading iHeightHandler. iHeightHandler must be 1 or 2. You entered:'
        write(errlog(io),*) linebuffer080
        write(errlog(io),*) 'Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(3)
      stop 1

      !Block 3/Line 3: Simulation time in hours
9303  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error reading simulation time in hours.',&
                    '  Program stopped'        
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(3)
      stop 1
        
      !Block 3/Line 4: stop computation when 99% of erupted mass has deposited? (yes/no)
9304  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'Error reading whether to stop simulation when'
        write(errlog(io),*)  '99% of erupted volume has deposited.'
        write(errlog(io),*)  'Answer should be yes or no.'
        write(errlog(io),*)  'You gave: ',linebuffer080
        write(errlog(io),*)  'Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(3)
      stop 1

      !Block 3/Line 5: nWindFiles
9305  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error reading number of wind files.'
        write(errlog(io),*)  'Answer should be a positive integer.'
        write(errlog(io),*)  'You gave: ',linebuffer080
        write(errlog(io),*)  ' Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(3)
      stop 1 

      !Block 4/Line 1: WriteDepositFinal_ASCII
      !                Write out ESRI ASCII file of final deposit thickness?
9401  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)'Error reading whether to print out ESRI ASCII file of',&
                  ' final deposit thickness.'
        write(errlog(io),*)'The first characters on this line should be ''yes'' or',&
                  ' ''no''.  Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(4)
      stop 1

      !Block 4/Line 2: WriteDepositFinal_KML
      !                Write out KML file of final deposit thickness?
9402  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)'Error reading whether to print out KML file of',&
                  ' final deposit thickness.'
        write(errlog(io),*)'The first characters on this line should',&
                  ' be yes or no.  Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(4)
      stop 1

      !Block 4/Line 3: WriteDepositTS_ASCII
      !                Write out ESRI ASCII deposit files at specified times?
9403  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)'Error reading whether to write out ESRI ASCII deposit files at',&
                  ' specifiied times.'
        write(errlog(io),*)'The first characters on this line should be ''yes'' or',&
                  ' ''no''.  Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(4)
      stop 1

      !Block 4/Line 4: WriteDepositTS_KML
      !                Write out KML deposit files at specified times?
9404  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)'Error reading whether to print out KML deposit',&
                  ' files at specifiied times.'
        write(errlog(io),*)'The first characters on this line should',&
                  ' be yes or no.  Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(4)
      stop 1

      !Block 4/Line 5: WriteCloudConcentration_ASCII
      !                Write out ESRI ASCII files of ash-cloud concentration?
9405  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)'Error reading whether to print out ESRI ASCII files of',&
                  ' cloud concentration at specifiied times.'
        write(errlog(io),*)'The first characters on this line should be ''yes'' or',&
                  ' ''no''.  Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(4)
      stop 1

      !Block 4/Line 6: WriteCloudConcentration_KLM
      !                Write out KLM files of ash-cloud concentration?
9406  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)'Error reading whether to print out a KML file of',&
                 ' cloud concentration.'
        write(errlog(io),*)'The first characters on this line should',&
                  ' be yes or no.  Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(4)
      stop 1

      !Block 4/Line 7: WriteCloudHeight_ASCII
      !                Write out ESRI ASCII files of ash-cloud height?
9407  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)'Error reading whether to print out ESRI ASCII files of',&
                  ' cloud height.'
        write(errlog(io),*)'The first characters on this line should',&
                  ' be yes or no.  Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(4)
      stop 1

      !Block 4/Line 8: WriteCloudHeight_KML
      !                Write out KLM files of ash-cloud height?
9408  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)'Error reading whether to print out a KML file',&
                  ' of cloud height.'
        write(errlog(io),*)'The first characters on this line should',&
                  ' be yes or no.  Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(4)
      stop 1

      !Block 4/Line 9: WriteCloudLoad_ASCII
      !                Write out ESRI ASCII files of ash-cloud load (T/km2) at specified times?
9409  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)'Error reading whether to print out an ESRI ASCII file',&
                 ' of cloud load.'
        write(errlog(io),*)'The first characters on this line should',&
                  ' be yes or no.  Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(4)
      stop 1
   
      !Block 4/Line 10: WriteCloudLoad_KML
      !                 Write out KML files of ash-cloud load (T/km2) at specified times?
9410  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)'Error reading whether to print out a KML file ',&
                  ' of cloud load.'
        write(errlog(io),*)'The first characters on this line should',&
                  ' be yes or no.  You gave:'
        write(errlog(io),*) linebuffer080
        write(errlog(io),*) 'program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(4)
      stop 1

      !Block 4/Line 11: WriteDepositTime_ASCII
      !                 Write out ESRI ASCII file of deposit arrival times?
9411  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)'Error reading whether to print out an ESRI ASCII file ',&
                  ' of deposit arrival time.'
        write(errlog(io),*)'The first characters on this line should',&
                  ' be yes or no.  You gave:'
        write(errlog(io),*) linebuffer080
        write(errlog(io),*) 'program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(4)
      stop 1

      !Block 4/Line 12: WriteDepositTime_KML
      !                 Write out KML file of deposit arrival times?
9412  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)'Error reading whether to print out a KML file ',&
                  ' of deposit arrival time.'
        write(errlog(io),*)'The first characters on this line should',&
                  ' be yes or no.  You gave:'
        write(errlog(io),*) linebuffer080
        write(errlog(io),*) 'program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(4)
      stop 1

      !Block 4/Line 13: WriteCloudTime_ASCII
      !                 Write out ESRI ASCII file of cloud arrival times?
9413  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)'Error reading whether to print out an ESRI ASCII file ',&
                  ' of cloud arrival time.'
        write(errlog(io),*)'The first characters on this line should',&
                  ' be yes or no.  You gave:'
        write(errlog(io),*) linebuffer080
        write(errlog(io),*) 'program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(4)
      stop 1

      !Block 4/Line 14: WriteCloudTime_KML
      !                 Write out KML file of cloud arrival times?
9414  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)'Error reading whether to print out a KML file ',&
                  ' of cloud arrival time.'
        write(errlog(io),*)'The first characters on this line should',&
                  ' be yes or no.  You gave:'
        write(errlog(io),*) linebuffer080
        write(errlog(io),*) 'program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(4)
      stop 1

      !Block 4/Line 15: Write3dFiles [iform]
      !                 Write out 3-D ash concentration at specified times? / [output code: 1=2d+concen,2=2d only]
9415  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)'Error reading whether to print out 3-D ash',&
                  ' concentration files at specifiied times.'
        write(errlog(io),*)'The first characters on this line should be ''yes'' or',&
                  ' ''no''.  Program stopped'
        write(errlog(io),*)'You gave:'
        write(errlog(io),*) linebuffer080
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(4)
      stop 1

      !Block 4/Line 16: formatanswer
      !                 format of ash concentration files   ("ascii", "binary", or "netcdf")
9416  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)'Error reading format of output files.'
        write(errlog(io),*)'The first characters on this line should',&
                  ' be ''ascii'', ''binary'', or ''netcdf''.'
        write(errlog(io),*)'Program stopped.'
        write(errlog(io),*)'You gave:'
        write(errlog(io),*) linebuffer080
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(4)
      stop 1

      !Block 4/Line 17: 
      !
9417  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)'Error reading the number of files to be written out.',&
                  '  This should be a positive integer, or -1.'
        write(errlog(io),*)'You gave: ',nWriteTimes
        write(errlog(io),*) 'Program stopped.'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(4)
      stop 1

      !Block 4/Line 18: WriteTimes()
9418  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)'Error reading the times at which output files are to',&
                  ' be written out.'
        write(errlog(io),*)'This should be one or more real numbers.'
        write(errlog(io),*)'You gave: ',WriteTimes
        write(errlog(io),*)'  Program stopped.'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(4)
      stop 1

      !Block 4/Line 18: WriteInterval
94181 do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)'Error reading the times at which output files are to',&
                  ' be written out.'
        write(errlog(io),*)'This should be one or more real numbers.'
        write(errlog(io),*)'You gave: ',WriteInterval
        write(errlog(io),*)'  Program stopped.'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(4)
      stop 1

      !Block 4/Line 18: Check on negative write times
94182 do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)'Error: some write times are <0.  Program stopped.'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(4)
      stop 1

      !Block 4/Line 18: Check on non-increasing write times
94183  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)'Error: some write times are not in chronological',&
                  ' order.  Program stopped.'      
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(4)
      stop 1

      !Block 5/all lines: Wind file name
9501  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error reading wind file names.  Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(5)
      stop 1
       
      !Block 6/Line 1: WriteAirportFile_ASCII
9601  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*) 'Error reading whether to write out ASCII airport file.'
        write(errlog(io),*) 'Answer must be yes or no.'
        write(errlog(io),*) 'You gave:',cdf_b6l1
        write(errlog(io),*) 'Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(6)
      stop 1

      !Block 6/Line 2: WriteGSD
9602  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*) 'Error reading whether to write out grain-size distribution to ASCII airport file.'
        write(errlog(io),*) 'Answer must be yes or no.'
        write(errlog(io),*) 'You gave:',cdf_b6l2
        write(errlog(io),*) 'Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(6)
      stop 1

      !Block 6/Line 3: WriteAirportFile_KML
9603  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*) 'Error reading whether to write out KML airport file.'
        write(errlog(io),*) 'Answer must be yes or no.'
        write(errlog(io),*) 'You gave:',cdf_b6l3
        write(errlog(io),*) 'Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(6)
      stop 1

      !Block 6/Line 4
      !If we can't open the external file
9604  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*) 'You gave the name of the following point location file to read:'
        write(errlog(io),*) AirportInFile
        write(errlog(io),*) 'But Ash3d could not open that file.'
      endif;enddo
      ! This is not an error with a stop point, program flow continues down to prompt user
      ! for the internal list of airports.

      !Block 6/Line 4: 
      !               Name of file containing airport locations
96041 do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*) 'Would you like Ash3d to use the internal airports database instead (y/n)?'
      endif;enddo
      read(input_unit,'(a1)',iostat=iostatus,iomsg=iomessage) answer
      linebuffer080 = answer
      linebuffer050 = "Reading answer from stdin (airport)"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      if(adjustl(trim(answer)).eq.'y') then
        ReadExtAirportFile=.false.
        AppendExtAirportFile=.false.
        AirportInFile='internal'
      else if(adjustl(trim(answer)).eq.'n') then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*) 'program stopped'
          write(errlog(io),*) '------------------------------'
        endif;enddo
        call help_inputfile(6)
        stop 1
      else
        goto 96041
      endif
      ! Return to normal program flow now that we have the airport file
      goto 96042

      !Block 6/Line 5: ProjectAirportLocations
      !                Have libprojection calculate projected coordinates?
9605  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*) 'Error reading whether to project airport coordinates using libprojection.'
        write(errlog(io),*) 'Answer must be yes or no.'
        write(errlog(io),*) 'You gave:',cdf_b6l5
        write(errlog(io),*) 'Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(6)
      stop 1

     !BLOCK 7: GRAIN-SIZE BINS, SETTLING VELOCITY
     !Block 7/Line 1: n_gs_max
9701  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*) 'Error reading number of grain-size bins to use.'
        write(errlog(io),*) 'Answer should be a positive integer.'
        write(errlog(io),*) 'You gave: ',linebuffer080
        write(errlog(io),*) 'Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(7)
      stop 1

     !Block 7/Line 2: Grain-size specification
9702  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*) 'Error reading grain-size specification.'
        write(errlog(io),*) 'There should be at least two values provided (FallVel (in m/s), mass fraction).'
        write(errlog(io),*) 'You gave: ',linebuffer080
        write(errlog(io),*) 'Program stopped'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(7)
      stop 1

     !BLOCK 8: VERTICAL PROFILES
     !Block 8/Line 1: nvprofiles
9801  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*) 'Error reading the number of vertical profiles to use.'
        write(errlog(io),*) 'Answer should be an integer.'
        write(errlog(io),*) 'You gave: ',linebuffer080
        write(errlog(io),*) 'Program stopped.'
        write(errlog(io),*) '------------------------------'
      endif;enddo
      call help_inputfile(8)
      stop 1

     !Block 8/Line 2+: x,y (or lon/lat) [Site name]
9802  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*) 'Error in x or y location of a vertical profile.'
        write(errlog(io),*) 'Answer should be two real numbers [+ Site name, if desired].'
        write(errlog(io),*) 'You gave: ',linebuffer080
        write(errlog(io),*) 'Program stopped.'
      endif;enddo
      stop 1

      !BLOCK 9: NETCDF ANNOTATIONS

!***********************************************************************
!     format statement

3     format(4x,'opening input file ',a130)
37    format(4x,'Volcano name:',a30,/)
4     format(4x,'MAP & GRID SETUP:',/,&
      4x,' coordinates of lower-left point:  x=',T48,f10.3,&
                 ', y=',T62,f10.3,/,&
      4x,'               extent of grid in:  x=',T48,f10.3,&
                 ', y=',T62,f10.3,/ &
      4x,'          coordinates of volcano:  x=',T48,f10.3,&
                 ', y=',T62,f10.3,/)
5     format(4x,'nodal spacing in x: ',f8.3,/,&
             4x,'              in y: ',f8.3)
43    format(4x,'              in z: ',f8.3,/)
6     format(4x,'ERUPTION SOURCE PARAMETERS:',/,&
             4x,'                diffusion coefficient (km2/hr): ',f8.4,/,&
             4x,'                               Suzuki constant: ',f8.4,/,&
             4x,'Stop calculation when 99% of ash has deposited: ',L1,/,&
             4x,'                         Simulation time (hrs): ',f8.4,/)
1438   format(4x,'ERUPTION SOURCE PARAMETERS:',/,&
             4x,'                diffusion coefficient (km2/hr): ',f8.4,/,&
             4x,'Stop calculation when 99% of ash has deposited: ',L1,/,&
             4x,'                         Simulation time (hrs): ',f8.4,/)
1439  format(4x,'                                     SourceType: ',a9)
7     format(4x,'Number of eruptions or pulses: ',i8,//, &
             4x,&
             '          plume        start time      duration   volume',/, &
             4x,&
             'number  height (km)  (yyyymmdd hh.hh)    (hrs)    (km3 DRE)')
8     format(4x,i6,f13.2,3x,i4,i2.2,i2.2,f6.2,2f10.4)
9     format(/,4x,'Number of grain-size bins:     ',i2, &
             /,4x,'               Fall Model:     ',i2)
10    format(4x,'Bins:',/,4x,&
        'mass fraction      diameter (mm)     density (kg/m3)      F     G     Sphr       phi')
11    format(8x,f10.5,4x,f11.6,10x,f10.4,5x,f8.2,2x,f4.2,2x,f4.2,5x,f5.2)
2110  format(4x,'Bins:',/,4x,&
             'mass fraction      v_s (m/s)')
2111  format(8x,f5.3,4x,f11.6)
13    format(4x,/,'WIND FILE SETUP:',/, &
                 4x,'Reading 4-D gridded wind data from files:')
16    format(4x,'Number of nodes in x:',T48,i3,/, &
             4x,'                in y:',T48,i3,/, &
             4x,'                in z:',T48,i3)
25    format(/,4x,'closing ',a80)
32     format(/,4x,'Warning:  some write times exceed the model',&
                   ' simulation time.',/, &
                4x,'These times will be ignored.')
33     format   (4x,'Chosen output:',/, &
                4x,'                        Output final ASCII deposit',&
                   ' file (T/F) = ',L1,/, &
                4x,'                          Output final KML deposit',&
                   ' file (T/F) = ',L1,/, &
                4x,'           Output ASCII deposit file at specified',&
                   ' times (T/F) = ',L1,/, &
                4x,'             Output KML deposit file at specified',&
                   ' times (T/F) = ',L1,/, &
                4x,'Output ASCII ash-cloud concentration at specified',&
                   ' times (T/F) = ',L1,/, &
                4x,'  Output KML ash-cloud concentration at specified',&
                   ' times (T/F) = ',L1,/, &
                4x,'       Output ASCII ash-cloud height at specified',&
                   ' times (T/F) = ',L1,/, &
                4x,'         Output KML ash-cloud height at specified',&
                   ' times (T/F) = ',L1,/, &
                4x,'             Output ASCII cloud load at specified',&
                   ' times (T/F) = ',L1,/, &
                4x,'               Output KML cloud load at specified',&
                   ' times (T/F) = ',L1,/, &
                4x,'             Output ASCII file of deposit arrival',&
                   '  time (T/F) = ',L1,/, &
                4x,'               Output KML file of deposit arrival',&
                   '  time (T/F) = ',L1,/, &
                4x,'               Output ASCII file of cloud arrival',&
                   '  time (T/F) = ',L1,/, &
                4x,'                 Output KML file of cloud arrival',&
                   '  time (T/F) = ',L1,/, &
                4x,'                     Output 3d files at specified',&
                   ' times (T/F) = ',L1,/, &
                4x,'                  Output format of 3d ash',&
                   ' concentration files = ',a6,/, &
                4x,'                                         Number',&
                   ' of time steps = ',i3,/)
34    format(/, 4x,'Files to be written out at the following hours',&
                   ' after the eruption start:')
35    format(8x,f10.3)

      end subroutine Read_Control_File

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     LatLonChecker
!
!     This subroutine checks the domain for errors if IsLatLon=.true.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine LatLonChecker(latLL,lonLL,lat_volcano,lon_volcano,gridwidth_e,gridwidth_n)

      real(kind=ip), intent(in)    :: latLL
      real(kind=ip), intent(in)    :: lonLL
      real(kind=ip), intent(in)    :: lat_volcano
      real(kind=ip), intent(inout) :: lon_volcano ! we might need to remap this value
      real(kind=ip), intent(in)    :: gridwidth_e
      real(kind=ip), intent(in)    :: gridwidth_n

      ! Make sure that latitude is between -90 and 90.
      if(abs(latLL).gt.90.0_ip.or.abs(lat_volcano).gt.90.0_ip)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: latLL and lat_volcano should be",&
                    " in in range -90 - 90"
          write(errlog(io),*)"latLL       = ",latLL
          write(errlog(io),*)"lat_volcano = ",lat_volcano
        endif;enddo
        stop 1
      endif
      
      ! Make sure that longitude of left side of grid is between -360 and 360.
      if(abs(lonLL).gt.360.0_ip.or.abs(lon_volcano).gt.360.0_ip)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*) "ERROR: lonLL and lon_volcano should be",&
                     " in in range -360 - 360"
          write(errlog(io),*)"lonLL       = ",lonLL
          write(errlog(io),*)"lon_volcano = ",lon_volcano
        endif;enddo
        stop 1
      endif

      ! Make sure that gridwidth_e and gridwidth_n are positive
      if((gridwidth_e.lt.0.0_ip).or.(gridwidth_n.lt.0.0_ip))then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: gridwidth_e and gridwidth_n must be positive."
        endif;enddo
        stop 1
      endif

      ! Make sure that the top of the grid does not extend beyond 90n
      if((latLL+gridwidth_n).gt.90.0_ip) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR:  Latitude at the top of the grid > 90."
        endif;enddo
        stop 1        
      endif

      ! Make sure that lon_volcano > lonLL
      if(lonLL.gt.lon_volcano) lon_volcano=lon_volcano+360.0_ip

      ! Make sure the volcano is within the model region (longitude)
      if((lon_volcano.lt.lonLL).or.(lon_volcano.gt.(lonLL+gridwidth_e))) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*) 'ERROR: the volcano is not within the specified longitude region'
          write(errlog(io),*) "lon_volcano=",lon_volcano,', lonLL=',lonLL
          write(errlog(io),*) 'grid_width=',gridwidth_e
        endif;enddo
        stop 1
      endif

      ! Make sure the volcano is within the model region (latitude)
      if((lat_volcano.lt.latLL).or.(lat_volcano.gt.(latLL+gridwidth_n))) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*) 'ERROR: the volcano is not within the specified latitude region'
          write(errlog(io),*) "lat_volcano=",lat_volcano,', latLL=',latLL
          write(errlog(io),*) 'grid_height=',gridwidth_n
        endif;enddo
        stop 1
      endif

      return
      
      end subroutine LatLonChecker
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     xyChecker
!
!     This subroutine checks the domain for errors if IsLatLon=.false.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine xyChecker(xLL,yLL,dx,dy,x_volcano,y_volcano,gridwidth_x,gridwidth_y)

      use global_param,    only : &
        EPS_SMALL

      real(kind=ip), intent(in) :: xLL,yLL
      real(kind=ip), intent(in) :: dx,dy
      real(kind=ip), intent(in) :: x_volcano,y_volcano
      real(kind=ip), intent(in) :: gridwidth_x,gridwidth_y

      ! Make sure thae gridwidth_x and gridwidth_y are positive
      if((gridwidth_x.lt.0.0_ip).or.(gridwidth_y.lt.0.0_ip))then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*) &
            "ERROR: gridwidth_x and gridwidth_y must be positive."
        endif;enddo
        stop 1
      endif

      ! Make sure that dx and dy are positive
      if((dx.lt.0.0_ip).or.(dy.lt.0.0_ip)) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*) "ERROR: dx and dy must be positive."
        endif;enddo
        stop 1
      endif

      ! Make sure the volcano is within the model region
      if(((x_volcano.lt.xLL).or.(x_volcano.gt.(xLL+gridwidth_x))).or. &
          ((y_volcano.lt.yLL).or.(y_volcano.gt.(yLL+gridwidth_y)))) then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*) &
                'ERROR: the volcano is not within the model region'
        endif;enddo
        stop 1
      endif

      ! Print out warning message if dx != dy
      if(abs(dx-dy).gt.EPS_SMALL) then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),1)              ! print out a warning message about the deposit file
        endif;enddo
      endif

      return

      ! Format statements
1     format (4x,'Warning: dx and dy are not the same.  If the',&
                ' deposit file is read by ArcMap',/, &
             4x,'they are assumed to be the same.  The nodal',&
                ' spacing written to the ',/, &
             4x,'deposit file is dx, dy, but ArcMap will only read dx.')

      end subroutine xyChecker
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    vprofchecker
!
!    Subroutine that checks the locations of vertical profiles specified
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine vprofchecker(iprof)

      use mesh,              only : &
        IsLatLon,de,dn,lonLL,latLL,gridwidth_e,gridwidth_n,&
        dx,dy,xLL,yLL,gridwidth_x,gridwidth_y

      use io_data,           only : &
        x_vprofile,y_vprofile,i_vprofile,j_vprofile

      integer, intent(in) :: iprof
      real(kind=ip) :: lon_vprof, lat_vprof

      ! Find the i and j values of the node containing the vertical profile
      if(IsLatLon) then
        lon_vprof = x_vprofile(iprof)
        lat_vprof = y_vprofile(iprof)
        if((lonLL+gridwidth_e-lon_vprof).gt.360.0_ip) &
          lon_vprof = lon_vprof+360.0_ip
        i_vprofile(iprof) = int((lon_vprof-lonLL)/de) + 1
        j_vprofile(iprof) = int((lat_vprof-latLL)/dn) + 1
      else
        i_vprofile(iprof) = int((x_vprofile(iprof)-xLL)/dx) + 1
        j_vprofile(iprof) = int((y_vprofile(iprof)-yLL)/dy) + 1
      endif
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),47) iprof, real(x_vprofile(iprof),kind=4),&
                                     real(y_vprofile(iprof),kind=4),&
                              i_vprofile(iprof), j_vprofile(iprof)
      endif;enddo
47    format(i11,2f10.3,2i6)

      ! Make sure the point is within the model region
      if(IsLatLon) then
        if(((lon_vprof.lt.lonLL)               .or. &
             (lon_vprof.gt.(lonLL+gridwidth_e))).or. &
            ((lat_vprof.lt.latLL)               .or. &
             (lat_vprof.gt.(latLL+gridwidth_n)))) then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*) &
              'ERROR: this location is not within the model region'
          endif;enddo
          stop 1
        endif
      else
        if(((x_vprofile(iprof).lt.xLL)               .or.&
             (x_vprofile(iprof).gt.(xLL+gridwidth_x))).or.&
            ((y_vprofile(iprof).lt.yLL)               .or.&
            (y_vprofile(iprof).gt.(yLL+gridwidth_y)))) then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)&
                'ERROR: this location is not within the model region'
            endif;enddo
            stop 1
        endif
      endif

      return

      end subroutine vprofchecker

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Subroutine Read_PostProc_Control_File()
!
! This subroutine sets up the parameters for the Ash3d Post-processing task via
! a control file.
!
! 3d_tephra_fall.nc       data_filename
! 3                       input format code (1=ascii, 2=binary, 3=netcdf, 4=grib, 5=tecplot, 6=vtk)
! test.inp                Ash3d_control_file_name   (skipped if datafile is netcdf)
! 5                       invar [outvar]    input variable and output variable, if different
!  if invar=0 varname
! 2                       ndims Only needed for format 1 or 2
! 0 0                     nx ny [nz] Also only needed for format 1 or 2
! 0.0 0.0                 dx dy [dz] Needed if not a part of the data file
! 0.0 0.0                 srtx srty         Start x and y
! 3                       output format     (1=ascii, 2=KML 3=image, 4=binary, 5=shapefile 6=grib, 7=tecplot, 8=vtk)
! 3                       plot_pref         (1=dislin, 2=plplot, 3=gnuplot, 4=GMT)
! -1                      time_step         Only needed if input file is multi-timestep (eg netcdf) (0 for static, -1 for final, -2 for all)
! 0                       Filled contour flag (1 for filled, 0 for lines)
! 1 5                     custom contour flag (1 for true, 0 for false), number of contours
! 1.0 3.0 10.0 50.0 100.0 lev(ncont)  : contour levels
! 100 100 100 100 255     R(ncont)    : Red channel of RGB
! 100 150 200 150   0     G(ncont)    : Green channel of RGB
! 200 150 100  50   0     B(ncont)    : Blue channel of RGB
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Read_PostProc_Control_File(informat,iprod1,iprod2,ndims,outformat,iplotpref,itime)

      use io_data,       only : &
         PP_infile,concenfile,datafileIn,infile,HaveInfile,&
         nvar_User2d_static_XY,nvar_User2d_XY

      use mesh,          only : &
         nxmax,nymax,nzmax,dx,dy,dz_const,xLL,yLL,IsLatLon,lonLL,latLL,de,dn

      use Output_Vars,   only : &
         ContourFilled,Con_Cust,Con_Cust_N,Con_Cust_RGB,Con_Cust_Lev, &
         Extra2dVarName

      integer, intent(out) :: informat
      integer, intent(out) :: iprod1
      integer, intent(out) :: iprod2
      integer, intent(out) :: ndims
      integer, intent(out) :: outformat
      integer, intent(out) :: iplotpref
      integer, intent(out) :: itime

      character(len=50)  :: linebuffer050 
      character(len=80)  :: linebuffer080
      character(len=10)  :: tmpstr
      integer            :: iostatus
      character(len=120) :: iomessage
      integer            :: ivalue
      real(kind=ip)      :: rvalue
      integer            :: i
      integer            :: ilatlonflag

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"--------------------------------------------------"
        write(outlog(io),*)"------ READ_POSTPROC_CONTROL_FILE ----------------"
        write(outlog(io),*)"--------------------------------------------------"
      endif;enddo

      ! Open the control file. Note that we have already checked for existance
      open(unit=fid_ctrlfile,file=PP_infile,status='old',action='read',err=9001)
      ! Line 1:
      !  The first line is the data file.  This is expected to be the netcdf file
      !  with all the run information, but could be a binary or ASCII output file.
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading line 1 of post-proc control file"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) datafileIn
      linebuffer050 = "Reading datafileIn from linebuffer"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      datafileIn = adjustl(trim(datafileIn))

      ! Line 2:
      !  This is the format code of the data file from line 1.  Currently, we expect
      !  this to be the netcdf run output file, but could be binary or ASCII.  The
      !  other formats are placeholders and not yet implemented.
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading line 2 of post-proc control file"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) informat
      linebuffer050 = "Reading informat from linebuffer"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      if(informat.ne.1.and.&
         informat.ne.2.and.&
         informat.ne.3)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Invalid format code for data file."
        endif;enddo
        stop 1
      else
        if(informat.eq.1)tmpstr="ASCII"
        if(informat.eq.2)tmpstr="Binary"
        if(informat.eq.3)tmpstr="NetCDF"
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"File format of input data    = ",informat,tmpstr
        endif;enddo
      endif
      ! If this is a netcdf file, copy name to 
      if(informat.eq.3)then
        concenfile = adjustl(trim(datafileIn))
      endif

      ! Line 3:
      !  The name of the input file used for this run. If the datafile is the
      !  netcdf concentratiln file, then the contents of the input file are
      !  already available and this line will be ignored.  However, ASCII and
      !  binary files can use the additional information from the input file
      !  for the output products.
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading line 3 of post-proc control file"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) infile
      linebuffer050 = "Reading infile from linebuffer"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      infile = adjustl(trim(infile))
      if(infile(1:4).eq.'none')then
        HaveInfile = .false.
      else
        HaveInfile = .true.
      endif

      ! Line 4:
      !  The code for the variable to read and optionally an output code.
      !  Normally, the variable read wil be what is written out, but you
      !  could have a binary 3d concentration and want cloud_load or have
      !  an ASCII deposit (in mm) and want a plot of deposit in inches.
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading line 4 of post-proc control file"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) iprod1
      linebuffer050 = "Reading iprod` from linebuffer"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) iprod1, ivalue
      if(iostatus.eq.0)then
        if(iprod1.eq.0)then
          ! If we are reading in a non-standard product, then interpret the second integer
          ! to be the TS or no flag
          iprod2 = iprod1
          if(ivalue.eq.1)then
            nvar_User2d_static_XY = 1
          else
            nvar_User2d_XY = 1
          endif
        else
          ! If this is a standard product code, then second int is the outproduct
          iprod2 = ivalue
        endif
      else
        if(iprod1.eq.0)then
          ! non-standard product given, but no TS code, assume TS
          nvar_User2d_XY = 1
        else
          ! standard product 1, but no product 2; assume the same
          iprod2 = iprod1
        endif
      endif
      ! Error-checking these values
      if(iprod1.lt.0.or.iprod1.gt.16)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Invalid format code input variable type."
        endif;enddo
        stop 1
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Variable code of input data  = ",iprod1
        endif;enddo
      endif
      if(iprod2.lt.0.or.iprod2.gt.16)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Invalid format code output variable type."
        endif;enddo
        stop 1
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Variable code of output data  = ",iprod2
        endif;enddo
      endif
        ! Now the bonus line if the user requests a custom variable
      if(iprod1.eq.0)then

        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        linebuffer050 = "Reading custom variable name from control file"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) Extra2dVarName
        linebuffer050 = "Reading datafileIn from linebuffer"

        !do io=1,2;if(VB(io).le.verbosity_error)then
        !  write(errlog(io),*)"ERROR: Need to code custom variables in control file"
        !endif;enddo
        !stop 1
      endif

      ! Line 5:
      !  number of dimensions of the input data file (2 or 3)
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading line 5 of post-proc control file"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) ndims, ilatlonflag
      linebuffer050 = "Reading ndims,ilatlonfalg from linebuffer"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      ! Error-checking ndims
      if(ndims.ne.2.and.ndims.ne.3)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Invalid number of dimensions. Should be 2 or 3."
        endif;enddo
        stop 1
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Number of dimensions of data = ",ndims
        endif;enddo
      endif
      if(ilatlonflag.eq.0) then
        ! expecting input variables to be in the same projection as
        ! specified by iprojflag and parameters
        IsLatLon          = .false.
      else
        ! expecting input variables to be in lat/lon
       IsLatLon          = .true.
      endif

      ! Line 6:
      !  size of array in nx,ny,nz
      !  This is overwritten when reading netcdf or ASCII data, but is needed for binary
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading line 6 of post-proc control file"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) nxmax,nymax
      linebuffer050 = "Reading nx,ny from linebuffer"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) nxmax,nymax,ivalue
      if(iostatus.eq.0)then
        nzmax = ivalue
      else
        ! Do a hard stop if nz needs to be provided.
        if(ndims.eq.3.and.informat.eq.3)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: 3D binary data requires a nz."
          endif;enddo
          stop 1
        endif
        nzmax = 1
      endif
      ! Error-checking nx,ny,nz
      if(nxmax.lt.1.or.nymax.lt.1.or.nzmax.lt.1)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"WARNING: One or more dimension has non-positive size."
        endif;enddo
        if(informat.eq.2)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: Datafile is binary but dimesion size is invalid."
          endif;enddo
          stop 1
        endif
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Size of input data grid = ",nxmax,nymax,nzmax
        endif;enddo
      endif

      ! Line 7:
      !  cell-size dx,dy,dz
      !  This is overwritten when reading netcdf or ASCII data, but is needed for binary
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading line 7 of post-proc control file"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) dx,dy
      linebuffer050 = "Reading dx,dy from linebuffer"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) dx,dy,rvalue
      if(iostatus.eq.0)then
        dz_const = rvalue
      else
        ! Do a hard stop if dz needs to be provided.
        if(ndims.eq.3.and.informat.eq.3)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: 3D binary data requires a dz."
          endif;enddo
          stop 1
        endif
        dz_const = 1.0_ip
      endif
      ! Error-checking cell dimensions
      if(dx.le.0.0_ip.or.dy.le.0.0_ip.or.dz_const.le.0.0_ip)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"WARNING: One or more cell dimension has an invalid thickness."
        endif;enddo
        if(informat.eq.2)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: Datafile is binary but cell dimesions are invalid."
          endif;enddo
          stop 1
        endif
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Size of input data cell = ",dx,dy,dz_const
        endif;enddo
      endif

      ! Line 8:
      !  x,y of lower-left corner of grid
      !  This is overwritten when reading netcdf or ASCII data, but is needed for binary
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading line 8 of post-proc control file"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) xLL,yLL
      linebuffer050 = "Reading xLL,yLL from linebuffer"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      ! Since currently, this only works for Lon/Lat grids, copy these values
      lonLL = xLL
      latLL = yLL
      de    = dx
      dn    = dy
      ! Error-checking cell coordinate
      !if()then
      !  do io=1,2;if(VB(io).le.verbosity_error)then
      !    write(errlog(io),*)"WARNING: Lower-left cell coordinate is invalid"
      !  endif;enddo
      !  if(informat.eq.2)then
      !    do io=1,2;if(VB(io).le.verbosity_error)then
      !      write(errlog(io),*)"ERROR: Datafile is binary but cell starting coordinate is invalid."
      !    endif;enddo
      !    stop 1
      !  endif
      !else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Lower-left coordinate of computational grid = ",xLL,yLL
        endif;enddo
      !endif

      ! Line 9:
      !  This is the format code of the output file.
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading line 9 of post-proc control file"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)outformat
      linebuffer050 = "Reading outformat from linebuffer"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      if(outformat.ne.1.and.&
          outformat.ne.2.and.&
          outformat.ne.3.and.&
          outformat.ne.4.and.&
          outformat.ne.5)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Invalid format code for output product."
        endif;enddo
        stop 1
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Output product format = ",outformat
        endif;enddo
      endif

      ! Line 10:
      !  preferred plotting library
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading line 10 of post-proc control file"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) iplotpref
      linebuffer050 = "Reading iplotpref from linebuffer"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      ! Error-checking iplotpref
      if(iplotpref.lt.0.or.iplotpref.gt.4)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Invalid code for plotting library."
        endif;enddo
        stop 1
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Preferred plotting package = ",iplotpref
        endif;enddo
      endif

      ! Line 11:
      !  output time step
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading line 11 of post-proc control file"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) itime
      linebuffer050 = "Reading itime from linebuffer"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      ! Error-checking itime
      if(itime.lt.-2)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Invalid timestep requested."
        endif;enddo
        stop 1
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Time step requested = ",itime
        endif;enddo
      endif

      ! Line 12:
      !  The flag for filled contours as opposed to the default lines.
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading line 12 of post-proc control file"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) ivalue
      linebuffer050 = "Reading ivalue from linebuffer"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      if(ivalue.eq.1)then
        ContourFilled = .true.
      endif
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Draw flooded/filled contours = ",ContourFilled
      endif;enddo 

     ! Line 13+:
      !  Custom contours
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      linebuffer050 = "Reading line 13 of post-proc control file"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) ivalue
      linebuffer050 = "Reading ivalue from linebuffer"
      if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
      if(ivalue.eq.1)then
        Con_Cust = .true.
      endif
      if(Con_Cust)then
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) ivalue, Con_Cust_N
        if(iostatus.eq.0)then
          ! Success reading the number of custom levels; run error-check
          if(Con_Cust_N.le.0)then
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: Number of contour levels must be positive."
            endif;enddo
            stop 1
          else
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"Custom contours (# of levels) ",Con_Cust,Con_Cust_N
            endif;enddo
          endif
        else
          ! User wants custom contours, but we can't read the number
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: Number of contour levels must be positive."
          endif;enddo
          stop 1
        endif
        ! Now read an additional 4 lines with contour info:
        allocate(Con_Cust_RGB(Con_Cust_N,3)); Con_Cust_RGB(:,:) = 0
        allocate(Con_Cust_Lev(Con_Cust_N));   Con_Cust_Lev(:)   = 0.0_ip
        ! Line 14: custom levels
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        linebuffer050 = "Reading line 14 of post-proc control file"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)Con_Cust_Lev(1:Con_Cust_N)
        linebuffer050 = "Reading Con_Cust_Lev from linebuffer"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)

        ! Line 15: custom color (R)
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        linebuffer050 = "Reading line 15 of post-proc control file"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)Con_Cust_RGB(1:Con_Cust_N,1)
        linebuffer050 = "Reading Con_Cust_RGB from linebuffer"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        ! Line 16: custom color (G)
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        linebuffer050 = "Reading line 16 of post-proc control file"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)Con_Cust_RGB(1:Con_Cust_N,2)
        linebuffer050 = "Reading Con_Cust_RGB from linebuffer"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        ! Line 17: custom color (B)
        read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
        linebuffer050 = "Reading line 17 of post-proc control file"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)Con_Cust_RGB(1:Con_Cust_N,3)
        linebuffer050 = "Reading Con_Cust_RGB from linebuffer"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Custom contours levels"
          do i=1,Con_Cust_N
            write(outlog(io),*)i,real(Con_Cust_Lev(i),kind=4),Con_Cust_RGB(i,1:3)
          enddo
        endif;enddo

      endif

      close(fid_ctrlfile)

      return

!******************************************************************************
      ! Error traps (starting with 9000)
      ! For this subroutine, the 100's position refers to block # of control file

      !ERROR TRAPS TO STDIN
9001  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error: cannot open input file: ',PP_infile
        write(errlog(io),*)  'Program stopped'
      endif;enddo
      stop 1

      end subroutine Read_PostProc_Control_File

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Write_input_block_header(outunit,blockID)
!
!  Called from: help_inputfile
!  Arguments:
!    outunit = output unit for write stream (6 for stdout, etc.)
!    blockID = block number of the control file to print
!
!  This subroutine writes an example header of the requested block of the
!  control file.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Write_input_block_header(outunit,blockID)

      use io_units

      integer,intent(in) :: outunit
      integer,intent(in) :: blockID

      ! The idea with the blockID is that help for only a particular block is
      ! called if there is an error reading something in the input file or if
      ! the user requests it.

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine write_input_block_header"
      endif;enddo

      select case (blockID)
        case(1) ! BLOCK 1: GRID INFO
      write(outunit,1)'# The following is an input file to the model Ash3d, v1.0 https://code.usgs.gov/vsc/ash3d/volcano-ash3d'
      write(outunit,1)'# Created by L.G. Mastin, R.P. Denlinger and H.F. Schwaiger, U.S. Geological Survey, 2009.             '
      write(outunit,1)'#                                                                                                      '
      write(outunit,1)'# GENERAL SOURCE PARAMETERS. DO NOT DELETE ANY NON-COMMENT LINES                                       '
      write(outunit,1)'#  The first line of this block identifies the volcano by name.                                        '
      write(outunit,1)'#  If the volcano name begins with either 0 or 1, then the volcano                                     '
      write(outunit,1)'#  is assumed to be in the Smithsonian database and default values for                                 '
      write(outunit,1)'#  Plume Height, Duration, Mass Flux Rate, Volume, and mass fraction of                                '
      write(outunit,1)'#  fines are loaded.  These can be over-written by entering non-negative                               '
      write(outunit,1)'#  values in the appropriate locations in this input file.                                             '
      write(outunit,1)'#                                                                                                      '
      write(outunit,1)'#  The second line of this block identifies the projection used and the form of                        '
      write(outunit,1)'#  the input coordinates and is of the following format:                                               '
      write(outunit,1)'#    latlonflag, projflag,  followed by a variable list of projection parameters                       '
      write(outunit,1)'#  projflag describes the projection used for the Ash3d run. Windfiles can have a                      '
      write(outunit,1)'#  different projection.                                                                               '
      write(outunit,1)'#  For a particular projflag, additional values are read defining the projection.                      '
      write(outunit,1)'#    latlonflag = 0 if computational grid is projected                                                 '
      write(outunit,1)'#               = 1 if computational grid is lat/lon (all subsequent projection parameters ignored.)   '
      write(outunit,1)'#    projflag   = 1 -- polar stereographic projection                                                  '
      write(outunit,1)'#           lambda0 -- longitude of projection point                                                   '
      write(outunit,1)'#           phi0    -- latitude of projection point                                                    '
      write(outunit,1)'#           k0      -- scale factor at projection point                                                '
      write(outunit,1)'#           radius  -- earth radius for spherical earth                                                '
      write(outunit,1)'#     e.g. for NAM 104,198, 216: 0 1 -105.0 90.0 0.933 6371.229                                        '
      write(outunit,1)'#               = 2 -- Albers Equal Area ( not yet implemented)                                        '
      write(outunit,1)'#               = 3 -- UTM ( not yet implemented)                                                      '
      write(outunit,1)'#               = 4 -- Lambert conformal conic                                                         '
      write(outunit,1)'#           lambda0 -- longitude of origin                                                             '
      write(outunit,1)'#              phi0 -- latitude of origin                                                              '
      write(outunit,1)'#              phi1 -- latitude of secant1                                                             '
      write(outunit,1)'#              phi2 -- latitude of secant2                                                             '
      write(outunit,1)'#            radius -- earth radius for a spherical earth                                              '
      write(outunit,1)'#     e.g. for NAM 212: 0 4 265.0 25.0 25.0 25.0 6371.22                                               '
      write(outunit,1)'#               = 5 -- Mercator                                                                        '
      write(outunit,1)'#           lambda0 -- longitude of origin                                                             '
      write(outunit,1)'#              phi0 -- latitude of origin                                                              '
      write(outunit,1)'#            radius -- earth radius for a spherical earth                                              '
      write(outunit,1)'#     e.g. for NAM 196: 0 5 198.475 20.0 6371.229                                                      '
      write(outunit,1)'#                                                                                                      '
      write(outunit,1)'# On line 3, the vent coordinates can optionally include a third value for elevation in km.            '
      write(outunit,1)'# If the vent elevation is not given, 0 is used if topography is turned off.                           '
      write(outunit,1)'#                                                                                                      '
      write(outunit,1)'# Line 4 is the width and height of the computational grid in km (if projected) or degrees.            '
      write(outunit,1)'# Line 5 is the vent x,y (or lon, lat) coordinates.                                                    '
      write(outunit,1)'# Line 6, DX and DY resolution in km or degrees (for projected or lon/lat grid, respectively)          '
      write(outunit,1)'# Line 7, DZ can be given as a real number, indicating the vertical spacing in km.                     '
      write(outunit,1)'# Alternatively, it can be given as dz_plin (piece-wise linear), dz_clog (constant-                    '
      write(outunit,1)'# logarithmic), or dz_cust (custom specification)                                                      '
      write(outunit,1)'# If dz_plin, then a second line is read containing:                                                   '
      write(outunit,1)'#   number of line segments (N) followed by the steps and step-size of each segment                    '
      write(outunit,1)'#   e.g. 4 6 0.25 5 0.5 5 1.0 10 2.0                                                                   '
      write(outunit,1)'#         This corresponds to 4 line segments with 6 cells of 0.25, then 5 cells of 0.5,               '
      write(outunit,1)'#         5 cells of 1.0, and finally 10 cells of 2.0                                                  '
      write(outunit,1)'# If dz_clog, then a second line is read containing:                                                   '
      write(outunit,1)'#   maximum z and number of steps of constant dlogz                                                    '
      write(outunit,1)'#   e.g. 30.0 30                                                                                       '
      write(outunit,1)'#         This corresponds to 30 steps from 0-30km with constant log-spacing                           '
      write(outunit,1)'# If dz_cust, then a second line is read containing:                                                   '
      write(outunit,1)'#   the number of dz values to read (ndz), followed by dz(1:ndz)                                       '
      write(outunit,1)'#   e.g. 20 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 5.5            '
      write(outunit,1)'#         This corresponds to 10 steps of 0.5, 9 steps of 1.5, followed by 1 step of 5.5               '
      write(outunit,1)'#                                                                                                      '
      write(outunit,1)'#                                                                                                      '
      write(outunit,1)'# Line 8 is the the diffusivity (m2/s) followed by the eruption specifier.  The                        '
      write(outunit,1)'# eruption specifier can be a real number, in which case it is assumed to be the                       '
      write(outunit,1)'# positive constant specifying the Suzuki distribution.  Alternatively, it can be                      '
      write(outunit,1)'#  umbrella     : Suzuki (const. = 12) with radial spreading of the plume                              '
      write(outunit,1)'#  umbrella_air : Suzuki (const. = 12) with radial spreading of the plume scaled to 5% of vol.         '
      write(outunit,1)'#  point        : all mass inserted in cell containing PlmH                                            '
      write(outunit,1)'#  linear       : mass uniformly distributed from z-vent to PlmH                                       '
      write(outunit,1)'# Line 9 : number of pulses to be read in BLOCK 2                                                      '
        case(2) ! BLOCK 2: ERUPTION PARAMETERS
      write(outunit,1)'# ERUPTION LINES (number = neruptions)                                                                 '
      write(outunit,1)'# In the following line, each line represents one eruptive pulse.                                      '
      write(outunit,1)'# Parameters are (1-4) start time (yyyy mm dd h.hh (UT)); (5) duration (hrs);                          '
      write(outunit,1)'#                  (6) plume height;                      (7) erupted volume (km3 DRE)                 '
      write(outunit,1)'# If neruptions=1 and the year is 0, then the model run in forecast mode where mm dd h.hh are          '
      write(outunit,1)'# interpreted as the time after the start of the windfile.  In this case, duration, plume              '
      write(outunit,1)'# height and erupted volume are replaced with ESP if the values are negative.                          '
      write(outunit,1)'# This applies to source types: suzuki, point, line, umbrella and umbrella_air.                        '
      write(outunit,1)'# For profile sources, an additional two values are read: dz and nz                                    '
      write(outunit,1)'# 2010 04 14   0.00   1.0     18.0  0.16 1.0 18                                                        '
      write(outunit,1)'# 0.01 0.02 0.03 0.03 0.04 0.04 0.05 0.06 0.06 0.070 0.08 0.08 0.09 0.09 0.09 0.08 0.06 0.02           '
        case(3) ! BLOCK 3: WIND PARAMETERS
      write(outunit,1)'# WIND OPTIONS                                                                                         '
      write(outunit,1)'# Ash3d will read from either a single 1-D wind sounding, or gridded, time-                            '
      write(outunit,1)'# dependent 3-D wind data, depending on the value of the parameter iwind.                              '
      write(outunit,1)'# For iwind = 1, read from a 1-D wind sounding                                                         '
      write(outunit,1)'#             2, read from 3D gridded ASCII files                                                      '
      write(outunit,1)'#             3/4, read directly from a single or multiple NetCDF files.                               '
      write(outunit,1)'#             5, read directly from multiple multi-timestep NetCDF files.                              '
      write(outunit,1)'# The parameter iwindformat specifies the format of the wind files, as follows:                        '
      write(outunit,1)'#  iwindformat =  0: User-defined via template                                                         '
      write(outunit,1)'#                 1: User-specified ASCII files                                                        '
      write(outunit,1)'#                 2: Global radiosonde data                                                            '
      write(outunit,1)'#                 3: NARR 221 Reanalysis (32 km)                                                       '
      write(outunit,1)'#                 4: NAM Regional North America 221 Forecast (32 km)                                   '
      write(outunit,1)'#                 5: NAM 216 Regional Alaska Forecast (45 km)                                          '
      write(outunit,1)'#                 6: NAM 104 Northern Hemisphere Forecast (90 km)                                      '
      write(outunit,1)'#                 7: NAM 212 40km Cont. US Forecast (40 km)                                            '
      write(outunit,1)'#                 8: NAM 218 12km Cont. US Forecast (12 km)                                            '
      write(outunit,1)'#                 9: NAM 227 Cont. US Forecast (5.08 km)                                               '
      write(outunit,1)'#                10: NAM 242 11km Regional Alaska Forecast (11.25 km)                                  '
      write(outunit,1)'#                11: NAM 196 Regional Hawaii Forecast (2.5 km)                                         '
      write(outunit,1)'#                12: NAM 198 Regional Alaska Forecast (5.953 km)                                       '
      write(outunit,1)'#                13: NAM 91 Regional Alaska Forecast (2.976 km)                                        '
      write(outunit,1)'#                14: NAM Regional Cont. US Forecast (3.0 km)                                           '
      write(outunit,1)'#                20: GFS 0.5 degree files Forecast                                                     '
      write(outunit,1)'#                21: GFS 1.0 degree files Forecast                                                     '
      write(outunit,1)'#                22: GFS 0.25 degree files Forecast                                                    '
      write(outunit,1)'#                23: NCEP DOE Reanalysis 2.5 degree                                                    '
      write(outunit,1)'#                24: NASA MERRA-2 Reanalysis                                                           '
      write(outunit,1)'#                25: NCEP1 2.5 global Reanalysis (1948-pres)                                           '
      write(outunit,1)'#                      Note: use nWindFiles=1 for iwindformat=25                                       '
      write(outunit,1)'#                26: JRA-55 Reanalysis                                                                 '
      write(outunit,1)'#                27: NOAA-CIRES II 2-deg global Reanalysis (1870-2010)                                 '
      write(outunit,1)'#                28: ECMWF ERA-Interim Reanalysis                                                      '
      write(outunit,1)'#                29: ECMWA ERA-5 Reanalysis                                                            '
      write(outunit,1)'#                30: ECMWA ERA-20C Reanalysis                                                          '
      write(outunit,1)'#                32: Air Force Weather Agency                                                          '
      write(outunit,1)'#                33: CCSM 3.0 Community Atmospheric Model                                              '
      write(outunit,1)'#                34: ECMWF 0.25-degree forecast                                                        '
      write(outunit,1)'#                40: NASA GEOS-5 Cp                                                                    '
      write(outunit,1)'#                41: NASA GEOS-5 Np                                                                    '
      write(outunit,1)'#                50: Weather Research and Forecast (WRF) output                                        '
      write(outunit,1)'#                                                                                                      '
      write(outunit,1)'# igrid (optional, defaults to that associated with iwindformat) is the NCEP grid ID,                  '
      write(outunit,1)'# if a NWP product is used, or the number of stations of sonde data, if iwind = 1.                     '
      write(outunit,1)'# idata (optional, defaults to 2) is a flag for data type (1=ASCII, 2=netcdf, 3=grib).                 '
      write(outunit,1)'#                                                                                                      '
      write(outunit,1)'# Many plumes extend higher than the maximum height of mesoscale models.                               '
      write(outunit,1)'# Ash3d handles this as determined by the parameter iHeightHandler, as follows:                        '
      write(outunit,1)'# for iHeightHandler = 1, stop the program if the plume height exceeds mesoscale height                '
      write(outunit,1)'#                      2, wind velocity at levels above the highest node                               '
      write(outunit,1)'#                         equal that of the highest node.  Temperatures in the                         '
      write(outunit,1)'#                         upper nodes do not change between 11 and 20 km; above                        '
      write(outunit,1)'#                         20 km they increase by 2 C/km, as in the Standard                            '
      write(outunit,1)'#                         atmosphere.  A warning is written to the log file.                           '
      write(outunit,1)'# Simulation time in hours is the maximal length of the simulation.                                    '
      write(outunit,1)'# Ash3d can end the simulation early if desired, once 99% of the ash has deposited.                    '
      write(outunit,1)'# The last line of this block is the number of windfiles listed in block 5 below.  If                  '
      write(outunit,1)'# iwind=5 and one of the NWP products is used that require a special file structure,                   '
      write(outunit,1)'# then nWindFiles should be set to 1 and only the root folder of the windfiles listed.                 '
        case(4) ! BLOCK 4: OUTPUT OPTIONS
      write(outunit,1)'# OUTPUT OPTIONS:                                                                                      '
      write(outunit,1)'# The list below allows users to specify the output options                                            '
      write(outunit,1)'# All but the final deposit file can be written out at specified                                       '
      write(outunit,1)'# times using the following parameters:                                                                '
      write(outunit,1)'# Line 15 asks for 3d output (yes/no) followed by an optional output format code;                      '
      write(outunit,1)'#   1 = (default) output all the normal 2d products to the output file as well as the 3d concentrations'
      write(outunit,1)'#   2 = only output the 2d products                                                                    '
      write(outunit,1)'# nWriteTimes   = if >0,  number of times output are to be written. The following                      '
      write(outunit,1)'#                  line contains nWriteTimes numbers specifying the times of output                    '
      write(outunit,1)'#                 if =-1, it specifies that the following line gives a constant time                   '
      write(outunit,1)'#                  interval in hours between write times.                                              '
      write(outunit,1)'# WriteTimes    = Hours between output (if nWritetimes=-1), or                                         '
      write(outunit,1)'#                 Times (hours since start of first eruption) for each output                          '
      write(outunit,1)'#                (if nWriteTimes >1)                                                                   '
        case(5) ! BLOCK 5: INPUT WIND FILES
      write(outunit,1)'# WIND INPUT FILES                                                                                     '
      write(outunit,1)'# The following block of data contains names of wind files. There should be one line for               '
      write(outunit,1)'# each of nWindFiles (from Block 3 Line 5) windfiles. Files should be given in                         '
      write(outunit,1)'# chronological order, should have names with only letters and numbers (no spaces)                     '
      write(outunit,1)'# and should not exceed 130 characters in length.                                                      '
      write(outunit,1)'# For iwind=5 (files with hard-coded paths), just provide the directory with the                       '
      write(outunit,1)'# windfiles or the root of the dataset (if files are sorted by year).                                  '
      write(outunit,1)'# For example, iwind=5, iwindformat=25 for NCEP reanalysis, data might look like:                      '
      write(outunit,1)'# /data/WindFiles/NCEP                                                                                 '
      write(outunit,1)'# |-- 2016                                                                                             '
      write(outunit,1)'# |   |-- air.2016.nc                                                                                  '
      write(outunit,1)'# |   |-- hgt.2016.nc                                                                                  '
      write(outunit,1)'# |   |-- omega.2016.nc                                                                                '
      write(outunit,1)'# |   |-- uwnd.2016.nc                                                                                 '
      write(outunit,1)'# |   `-- vwnd.2016.nc                                                                                 '
      write(outunit,1)'# |-- 2017                                                                                             '
      write(outunit,1)'#     |-- air.2017.nc                                                                                  '
      write(outunit,1)'# In this case, Block 5 will just contain one line: /data/WindFiles/NCEP or just NCEP                  '
      write(outunit,1)'# if you have a soft link in the run directory.                                                        '
      write(outunit,1)'# For a network of radiosonde data, please see the MetReader documentation for                         '
      write(outunit,1)'# the input specification https://code.usgs.gov/vsc/ash3d/volcano-ash3d-metreader.                     '
        case(6) ! BLOCK 6: AIRPORT FILE
      write(outunit,1)'# AIRPORT LOCATION FILE                                                                                '
      write(outunit,1)'# The following lines allow the user to specify whether times of ash arrival                           '
      write(outunit,1)'# at airports and other locations will be written out, and which file                                  '
      write(outunit,1)'# to read for a list of airport locations.                                                             '
      write(outunit,1)'# PLEASE NOTE:  Each line in the airport location file should contain the                              '
      write(outunit,1)'#               airport latitude, longitude, projected x and y coordinates,                            '
      write(outunit,1)'#               and airport name.  If you are using a projected grid,                                  '
      write(outunit,1)'#               THE X AND Y MUST BE IN THE SAME PROJECTION as the computational grid.                  '
      write(outunit,1)'#               Alternatively, coordinates can be projected via libprojection                          '
      write(outunit,1)'#               by typing "yes" to the last parameter                                                  '
        case(7) ! BLOCK 7: GRAIN-SIZE BINS, SETTLING VELOCITY
      write(outunit,1)'# GRAIN SIZE GROUPS                                                                                    '
      write(outunit,1)'# The first line must contain the number of settling velocity groups, but                              '
      write(outunit,1)'# can optionally also include a flag for the fall velocity model to be used.                           '
      write(outunit,1)'#    FV_ID = 1, Wilson and Huang                                                                       '
      write(outunit,1)'#          = 2, Wilson and Huang + Cunningham slip                                                     '
      write(outunit,1)'#          = 3, Wilson and Huang + Mod by Pfeiffer Et al.                                              '
      write(outunit,1)'#          = 4, Ganser (assuming prolate ellipsoids)                                                   '
      write(outunit,1)'#          = 5, Ganser + Cunningham slip                                                               '
      write(outunit,1)'#          = 6, Stokes flow for spherical particles + slip                                             '
      write(outunit,1)'# If no fall model is specified, FV_ID = 1, by default                                                 '
      write(outunit,1)'# The first line can also optionally contain the Shape_ID code where 1= Wilson and Huang (F,G),        '
      write(outunit,1)'# and 2 = sphericity                                                                                   '
      write(outunit,1)'# The grain size bins can be enters with 2, 3, or 4 parameters.                                        '
      write(outunit,1)'# If TWO are given, they are read as:   FallVel (in m/s), mass fraction                                '
      write(outunit,1)'# If THREE are given, they are read as: diameter (mm), mass fraction, density (kg/m3)                  '
      write(outunit,1)'# For Shape_ID = 1:                                                                                    '
      write(outunit,1)'# If FOUR are given, they are read as:  diameter (mm), mass fraction, density (kg/m3), Shape F         '
      write(outunit,1)'# The shape factor is given as in Wilson and Huang: F=(b+c)/(2*a), but converted                       '
      write(outunit,1)'# to sphericity (assuming b=c) for the Ganser model.                                                   '
      write(outunit,1)'# If a shape factor is not given, a default value of F=0.4 is used.                                    '
      write(outunit,1)'# If FIVE are given, they are read as:  diameter (mm), mass fraction, density (kg/m3), Shape F, G      '
      write(outunit,1)'#  where G is an additional minor axis ratio shape factor equal to c/b                                 '
      write(outunit,1)'# For Shape_ID = 2:                                                                                    '
      write(outunit,1)'# If FOUR are given, they are read as:  diameter (mm), mass fraction, density (kg/m3), Sphericity      '
      write(outunit,1)'#                                                                                                      '
      write(outunit,1)'# If the last grain size bin has a negative diameter, then the remaining mass fraction                 '
      write(outunit,1)'# will be distributed over the previous bins via a log-normal distribution in phi.                     '
      write(outunit,1)'# The last bin would be interpreted as:                                                                '
      write(outunit,1)'# diam (neg value) , LN_phi_mean, LN_phi_stddev                                                        '
        case(8) ! BLOCK 8: VERTICAL PROFILES
      write(outunit,1)'# Options for writing vertical profiles                                                                '
      write(outunit,1)'# The first line below gives the number of locations (nlocs) where vertical                            '
      write(outunit,1)'# profiles are to be written.  That is followed by nlocs lines, each of which                          '
      write(outunit,1)'# contain the location, in the same coordinates as the computational grid.                             '
      write(outunit,1)'# Optionally, a site name can be provided in after the location.                                       '
        case(9) ! BLOCK 9: (Optional): NETCDF ANNOTATIONS
      write(outunit,1)'# netCDF output options                                                                                '
      write(outunit,1)'# This last block is optional.                                                                         '
      write(outunit,1)'# The output file name can be give, but will default to 3d_tephra_fall.nc if absent                    '
      write(outunit,1)'# The title and comment lines are passed through to the netcdf header of the                           '
      write(outunit,1)'# output file.                                                                                         '
        case(10) ! BLOCK 10 (OPTMOD): Optional module blocks
                 !   First RESETPARAMS
      write(outunit,1)'# Optional Modules are identified by the text string at the top of the block                           '
      write(outunit,1)'# OPTMOD=[module name]                                                                                 '
      write(outunit,1)'# There will need to be a custom block reader in the module to read this section                       '
      write(outunit,1)'# section of the input file.  Below is the built-in example for resetting parameters.                  '
      write(outunit,1)'# You only need to include the line(s) for the parameters you want to reset. All                       '
      write(outunit,1)'# options are listed below.                                                                            '
        case(11) ! BLOCK 10+1 (OPTMOD):
                 !   TOPO
      write(outunit,1)'# Topography                                                                                           '
      write(outunit,1)'# Line 1 indicates whether or not to use topography followed by the integer flag                       '
      write(outunit,1)'#        describing how topography will modify the vertical grid.                                      '
      write(outunit,1)'#          0 = no vertical modification; z-grid remains 0-> top throughout the domain                  '
      write(outunit,1)'#          1 = shifted; s = z-z_surf; computational grid is uniformly shifted upward                   '
      write(outunit,1)'#              everywhere by topography                                                                '
      write(outunit,1)'#          2 = sigma-altitude; s=z_top(z-z_surf)/(z_top-z_surf); topography has decaying               '
      write(outunit,1)'#              influence with height                                                                   '
      write(outunit,1)'# Line 2 indicates the topography data format followed by the smoothing radius in km                   '
      write(outunit,1)'# Topofile format must be one of                                                                       '
      write(outunit,1)'#   1 : Gridded lon/lat (netcdf): ETOPO, GEBCO                                                         '
      write(outunit,1)'#   2 : Gridded Binary: NOAA GLOBE, GTOPO30                                                            '
      write(outunit,1)'#   3 : ESRI ASCII                                                                                     '
      write(outunit,1)'#  Line 3 is the file name of the topography data.                                                     '
      write(outunit,1)'#                                                                                                      '
        case(12) ! BLOCK 10+2 (OPTMOD):
                 !   VARDIFF
      write(outunit,1)'# Variable Diffusivity                                                                                 '
      write(outunit,1)'#   Line 1 indicates whether or not to use horizontal diffusivity followed by the                      '
      write(outunit,1)'#          type ID and value with                                                                      '
      write(outunit,1)'#             1  500.0 # constant horizontal diffusivity with specified value (m2/s)                   '
      write(outunit,1)'#             2  0.2   # Smagorinsky model with coefficient C ()                                       '
      write(outunit,1)'#             3  0.2   # Pielke model with coefficient C ()                                            '
      write(outunit,1)'#   Line 2 indicates whether or not to use vertical diffusivity                                        '
      write(outunit,1)'#   Line 3 indicates the boundary layer model and value (if model requires)                            '
      write(outunit,1)'#             1 500.0         # BL model 1=const ; value (m2/s)                                        '
      write(outunit,1)'#             2               #          2=none (Use Free-air throughout)                              '
      write(outunit,1)'#             3               #          3=Troen and Mahrt                                             '
      write(outunit,1)'#             4               #          4=Ulke                                                        '
      write(outunit,1)'#             5               #          5=Shir / Businger,Ayer                                        '
      write(outunit,1)'#   Line 4 indicates the Free-Air model and value (if model requires)                                  '
      write(outunit,1)'#             1 500.0         # Free-Air 1=const ; value (m2/s)                                        '
      write(outunit,1)'#             2               #          2=F(Ri)=Louis                                                 '
      write(outunit,1)'#             3               #          3=F(Ri)=Jacobson/Stull                                        '
      write(outunit,1)'#             4               #          4=F(Ri)=Betts                                                 '
      write(outunit,1)'#             5               #          5=F(Ri)=Hong                                                  '
      write(outunit,1)'#             6               #          6=F(Ri)=Collins                                               '
      write(outunit,1)'#   Line 5 contains the von Karman constant                                                            '
      write(outunit,1)'#   Line 6 contains the free-air mixing length (m)                                                     '
      write(outunit,1)'#   Line 7 is the critical Richardson number used in calculating atmospheric stability.                '
      write(outunit,1)'#                                                                                                      '

!        case default
      end select

 1    format(a103)

      end subroutine Write_input_block_header

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  SetWrite_input_block_01
!
!  Called from: help_inputfile
!  Arguments:
!    outunit         = output stream ID or 0 for Set Vars only with no output
!    vname           = volcano name
!    projline        = projection specification
!    LLx,LLy         = coordinates of lower-left corner of computational grid
!    widthx,widthy   = width (in x,y) of computational grid
!    x_in,y_in,z_in  = source coordinates
!    dx_in,dy_in     = cell width in x,y
!    dz_in           = dz, which may be over-written
!    kdiff           = diffusivity
!    Suzk            = Susuki parameter
!    nerup           = # of eruptions
!    dz_type         = index of dz class 1=const, 2= plin, 3=clog, 4=custom
!    dz_line         = bonus line for dz specification
!    src_type        = source_type 1=Suz,2=point,3=line,4=profile,5=umb,6=umb_air
!    
!  This subroutine writes the content of block 1 (Grid/Src info) of the Ash3d
!  control file.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SetWrite_input_block_01(outunit,vname,projline,LLx,LLy,widthx,widthy,&
                                        x_in,y_in,z_in,dx_in,dy_in,dz_in,kdiff,Suzk,nerup,   &
                                        dz_type,dz_line,src_type)

      use io_units

      use io_data,       only : &
         VolcanoName

      use mesh,          only : &
         lonLL,latLL,gridwidth_e,gridwidth_n,xLL,yLL,gridwidth_x,gridwidth_y, &
         IsLatLon,de,dn,dx,dy,VarDzType,dz_const

      use Source,        only : &
         lon_volcano,lat_volcano,x_volcano,y_volcano,z_volcano,Suzuki_A,      &
         neruptions,SourceType,SourceType_idx

      use Diffusion,     only : &
         diffusivity_horz

      integer           ,intent(in) :: outunit
      character(len=30) ,intent(in) :: vname
      character(len=80) ,intent(in) :: projline
      real(kind=ip)     ,intent(in) :: LLx
      real(kind=ip)     ,intent(in) :: LLy
      real(kind=ip)     ,intent(in) :: widthx
      real(kind=ip)     ,intent(in) :: widthy
      real(kind=ip)     ,intent(in) :: x_in
      real(kind=ip)     ,intent(in) :: y_in
      real(kind=ip)     ,intent(in) :: z_in
      real(kind=ip)     ,intent(in) :: dx_in
      real(kind=ip)     ,intent(in) :: dy_in
      real(kind=ip)     ,intent(in) :: dz_in
      real(kind=ip)     ,intent(in) :: kdiff
      real(kind=ip)     ,intent(in) :: SuzK  ! These might not be needed, but must be included even with dummy values
      integer           ,intent(in) :: nerup
       ! Some extra bits for dz and source options
      integer           ,intent(in) :: dz_type  ! 1=const, 2=dz_plin, 3=dz_clog, 4=dz_cust
      character(len=130),intent(in) :: dz_line
      integer           ,intent(in) :: src_type ! 1=Suz, 2=point, 3=linear, 4=profile, 5=umb, 6=umb_a

      integer :: ilatlonflag

      VolcanoName          = adjustl(trim(vname))

      read(projline,*)ilatlonflag
      if(ilatlonflag.eq.0) then
        ! expecting input variables to be in the same projection as
        ! specified by iprojflag and parameters
        IsLatLon          = .false.
      else
        ! expecting input variables to be in lat/lon
       IsLatLon          = .true.
      endif

      if(IsLatLon)then
        lonLL       = LLx
        latLL       = LLy
        gridwidth_e = widthx
        gridwidth_n = widthy
        lon_volcano = x_in
        lat_volcano = y_in
        z_volcano   = z_in
        de          = dx_in
        dn          = dy_in
      else
        xLL         = LLx
        yLL         = LLy
        gridwidth_x = widthx
        gridwidth_y = widthy
        x_volcano   = x_in
        y_volcano   = y_in
        z_volcano   = z_in
        dx          = dx_in
        dy          = dy_in
      endif

      if(dz_type.eq.1)then
        VarDzType = "dz_cons"
        dz_const  = dz_in
      elseif(dz_type.eq.2)then
        VarDzType = 'dz_plin'
        !   number of line segments (N) followed by the steps and step-size of each segment                    
        !   e.g. 4 6 0.25 5 0.5 5 1.0 10 2.0                                                                   
        !         This corresponds to 4 line segments with 6 cells of 0.25, then 5 cells of 0.5,               
        !         5 cells of 1.0, and finally 10 cells of 2.0   
        ! dz_line should contain this      
      elseif(dz_type.eq.3)then
        VarDzType = 'dz_clog'
        !   maximum z and number of steps of constant dlogz                                                    
        !   e.g. 30.0 30                                                                                       
        !         This corresponds to 30 steps from 0-30km with constant log-spacing        
        ! dz_line should contain this      
      elseif(dz_type.eq.4)then
        VarDzType = 'dz_cust'
        !   the number of dz values to read (ndz), followed by dz(1:ndz)                                       
        !   e.g. 20 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 5.5            
        !         This corresponds to 10 steps of 0.5, 9 steps of 1.5, followed by 1 step of 5.5               
        ! dz_line should contain this      
      endif
      diffusivity_horz = kdiff
      SourceType_idx = src_type
      if(src_type.eq.1)then
        SourceType='suzuki'
        Suzuki_A = SuzK
      elseif(src_type.eq.2)then
        SourceType='point'
      elseif(src_type.eq.3)then
        SourceType='line'
      elseif(src_type.eq.4)then
        SourceType='profile'
      elseif(src_type.eq.5)then
        SourceType='umbrella'
      elseif(src_type.eq.6)then
        SourceType='umbrella_air'
      endif
      neruptions = nerup

      if(outunit.gt.0)then
        write(outunit,1)&
         '******************* BLOCK 1 ****************************************************'
        write(outunit,2)VolcanoName
        write(outunit,3)adjustl(trim(projline))
        write(outunit,4)LLx,LLy
        write(outunit,5)widthx,widthy
        write(outunit,6)x_in,y_in,z_in
        write(outunit,7)dx_in,dy_in
        if(dz_type.eq.1)then
          write(outunit,8)dz_const
        else
          write(outunit,9)VarDzType
          write(outunit,10)dz_line
        endif
        if(src_type.eq.1)then
          write(outunit,11)diffusivity_horz,Suzuki_A
        else
          write(outunit,12)diffusivity_horz,SourceType
        endif
        write(outunit,13)neruptions
        write(outunit,1)&
         '********************************************************************************'
      endif

 1    format(a80)
 2    format(a30   ,10x,' # Volcano name (character*30)')
 3    format(a40       ,' # Proj flags and params; first term (LLflag) is 1, so all else ignored')
 4    format(2f13.3,14x,' # x, y of LL corner of grid (km, or deg. if latlongflag=1)')
 5    format(2f13.3,14x,' # grid width and height (km, or deg. if latlonflag=1)')
 6    format(3f13.3, 1x,' # vent location         (km, or deg. if latlonflag=1)')
 7    format(2f13.3,14x,' # DX, DY of grid cells  (km, or deg. if latlonflag=1)')
 8    format(1f13.3,27x,' # DZ of grid cells      (always km)')
 9    format(a7    ,33x,' # DZ of grid cells      (always km)')
 10   format(a80)
 11   format(2f13.3,14x,' # diffusion coefficient (m2/s), Suzuki constant')
 12   format(1f13.3,5x,a12,' # diffusion coefficient (m2/s), Suzuki constant')
 13   format(i5    ,35x,' # neruptions, number of eruptions or pulses')

      end subroutine SetWrite_input_block_01

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  SetWrite_input_block_02
!
!  Called from: help_inputfile
!  Arguments:
!    outunit    = output stream ID or 0 for Set Vars only with no output
!    nerup      = # of eruptions
!    src_type   = source_type 1=Suz,2=point,3=line,4=profile,5=umb,6=umb_air
!    e_ST       = eruption start times given as hourssince BaseYear
!    e_Dur      = eruption durations
!    e_PmH      = eruption plume heights
!    e_Vol      = eruption volumes
!    ep_dz      = dz of eruption profile
!    ep_nz      = nz of eruption profile
!    ep_Vol     = normalized volume of eruption profile
!
!  This subroutine writes the content of block 2 (Eruption Parameters) of the
!  Ash3d control file.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SetWrite_input_block_02(outunit,nerup,src_type,     &
                                        e_ST,e_Dur,e_PmH,e_Vol,         &
                                        ep_dz,ep_nz,ep_Vol)

      use io_units

      use time_data,     only : &
          BaseYear,useLeap

      integer           ,intent(in) :: outunit
      integer           ,intent(in) :: nerup
      integer           ,intent(in) :: src_type
      real(kind=dp),dimension(nerup),intent(in) :: e_ST
      real(kind=dp),dimension(nerup),intent(in) :: e_Dur
      real(kind=ip),dimension(nerup),intent(in) :: e_PmH
      real(kind=ip),dimension(nerup),intent(in) :: e_Vol
      real(kind=ip),dimension(nerup),intent(in) :: ep_dz
      integer      ,dimension(nerup),intent(in) :: ep_nz
      real(kind=ip),dimension(nerup,50),intent(in) :: ep_Vol

      integer :: i
      integer :: iyear
      integer :: imonth
      integer :: iday
      real(kind=dp) :: hour

      INTERFACE
        real(kind=8)  function HS_HourOfDay(HoursSince,byear,useLeaps)
          real(kind=8),intent(in) :: HoursSince
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_HourOfDay
        integer function HS_YearOfEvent(HoursSince,byear,useLeaps)
          real(kind=8),intent(in) :: HoursSince
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_YearOfEvent
        integer function HS_MonthOfEvent(HoursSince,byear,useLeaps)
          real(kind=8),intent(in) :: HoursSince
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_MonthOfEvent
        integer function HS_DayOfEvent(HoursSince,byear,useLeaps)
          real(kind=8),intent(in) :: HoursSince
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_DayOfEvent
      END INTERFACE

      if(outunit.gt.0)then
        write(outunit,1)&
         '******************* BLOCK 2 ****************************************************'
        do i=1,nerup
          iyear = HS_YearOfEvent(e_ST(i),BaseYear,useLeap)
          imonth= HS_MonthOfEvent(e_ST(i),BaseYear,useLeap)
          iday  = HS_DayOfEvent(e_ST(i),BaseYear,useLeap)
          hour  = HS_HourOfDay(e_ST(i),BaseYear,useLeap)
          if(src_type.eq.1.or.&
             src_type.eq.2.or.&
             src_type.eq.3.or.&
             src_type.eq.5.or.&
             src_type.eq.6)then
            write(outunit,2)iyear,imonth,iday,hour,e_Dur(i),e_PmH(i),e_Vol(i)
          else
            ! src_type = 4 (profile)
            write(outunit,3)iyear,imonth,iday,hour,e_Dur(i),e_PmH(i),e_Vol(i),ep_dz(i),ep_nz(i)
            write(outunit,4)real(ep_Vol(i,1:ep_nz(i)),kind=4)
          endif
        enddo
        write(outunit,1)&
         '********************************************************************************'
      endif

 1    format(a80)
 2    format(3i5,1x,1f8.3,2f15.5,g15.5)
 3    format(3i5,1x,1f8.3,2f15.5,g15.5,1f8.2,i5)
 4    format(*(f7.4))

      end subroutine SetWrite_input_block_02

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  SetWrite_input_block_03
!
!  Called from: help_inputfile  
!  Arguments:
!    outunit    = output stream ID or 0 for Set Vars only with no output
!    iw         = iwind        (windfile class)
!    iwf        = iwindformat  (windfile product)
!    igrid      = NCEP grid ID
!    idf        = data format
!    iHH        = iHeightHandler
!    sim_time   = Simulation time
!    comp_stop  = y/n on whether to stop the simulation early
!    nwindfiles = number of windfiles
!
!  This subroutine writes the content of block 3 (Wind Parameters) of the Ash3d
!  control file.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SetWrite_input_block_03(outunit,iw,iwf,igrid,idf,iHH,&
                                         sim_time,comp_stop,nwindfiles)

      use io_units

      integer           ,intent(in) :: outunit
      integer           ,intent(in) :: iw
      integer           ,intent(in) :: iwf
      integer           ,intent(in) :: igrid
      integer           ,intent(in) :: idf
      integer           ,intent(in) :: iHH
      real(kind=dp)     ,intent(in) :: sim_time
      logical           ,intent(in) :: comp_stop
      integer           ,intent(in) :: nwindfiles

      if(outunit.gt.0)then
        write(outunit,1)&
         '******************* BLOCK 3 ****************************************************'
        if(igrid.eq.0)then
          write(outunit,2)iw,iwf
        else
          write(outunit,3)iw,iwf,igrid,idf
        endif
        write(outunit,4)iHH
        write(outunit,5)sim_time
        if(comp_stop)then
          write(outunit,6)
        else
          write(outunit,7)
        endif
        write(outunit,8)nwindfiles
        write(outunit,1)&
         '********************************************************************************'
      endif

 1    format(a80)
 2    format(2i5,20x,'# iwind, iwindformat, [igrid, idata]')
 3    format(4i5,10x,'# iwind, iwindformat, [igrid, idata]')
 4    format(1i5,25x,'# iHeightHandler')
 5    format(f13.3,17x,'# Simulation time in hours')
 6    format('yes                           # stop computation when 99% of erupted mass has deposited?')
 7    format('no                            # stop computation when 99% of erupted mass has deposited?')
 8    format(1i5,25x,'# nWindFiles, number of gridded wind files (used if iwind>1)')

      end subroutine SetWrite_input_block_03

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  SetWrite_input_block_04
!
!  Called from: help_inputfile  
!  Arguments:
!    outunit                         = output stream ID or 0 for Set Vars only with no output
!    WriteDepositFinal_ASCII_c       = char:  B4L1 n/y Write out ESRI ASCII file of final deposit thickness?
!    WriteDepositFinal_KML_c         = char:  B4L2 n/y Write out        KML file of final deposit thickness?
!    WriteDepositTS_ASCII_c          = char:  B4L3 Write out ESRI ASCII deposit files at specified times?
!    WriteDepositTS_KML_c            = char:  B4L4 Write out        KML deposit files at specified times?
!    WriteCloudConcentration_ASCII_c = char:  B4L5 Write out ESRI ASCII files of ash-cloud concentration?
!    WriteCloudConcentration_KML_c   = char:  B4L6 Write out        KML files of ash-cloud concentration?
!    WriteCloudHeight_ASCII_c        = char:  B4L7 Write out ESRI ASCII files of ash-cloud height?
!    WriteCloudHeight_KML_c          = char:  B4L8 Write out        KML files of ash-cloud height?
!    WriteCloudLoad_ASCII_c          = char:  B4L9 Write out ESRI ASCII files of ash-cloud load (T/km2) at specified times? 
!    WriteCloudLoad_KML_c            = char:  B4L10 Write out        KML files of ash-cloud load (T/km2) at specified times?
!    WriteDepositTime_ASCII_c        = char:  B4L11 Write out ESRI ASCII file of deposit arrival times?
!    WriteDepositTime_KML_c          = char:  B4L12 Write out        KML file of deposit arrival times?
!    WriteCloudTime_ASCII_c          = char:  B4L13 Write out ESRI ASCII file of cloud arrival times
!    WriteCloudTime_KML_c            = char:  B4L14 Write out        KML file of cloud arrival times?
!    Write3dFiles_c                  = char:  B4L15 Write out 3-D ash concentration at specified times?
!    ifm                             = integer : B4L15+ output code: 1=2d+concen,2=2d only]
!    ofm                             = integer : B4L16 format of ash concentration files (1=ascii, 2=binary, or 3=netcdf)
!    nwt                             = integer : B4L17 nWriteTimes
!    wts)                            = real(8),dim(nwt) : B4L18 WriteTimes(1:nWriteTimes)

!  This subroutine writes the content of block 4 (Output Options) of the Ash3d
!  control file.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SetWrite_input_block_04(outunit,                                                       &
                                         WriteDepositFinal_ASCII_c,WriteDepositFinal_KML_c,             &
                                         WriteDepositTS_ASCII_c,WriteDepositTS_KML_c,                   &
                                         WriteCloudConcentration_ASCII_c,WriteCloudConcentration_KML_c, &
                                         WriteCloudHeight_ASCII_c,WriteCloudHeight_KML_c,               &
                                         WriteCloudLoad_ASCII_c,WriteCloudLoad_KML_c,                   &
                                         WriteDepositTime_ASCII_c,WriteDepositTime_KML_c,               &
                                         WriteCloudTime_ASCII_c,WriteCloudTime_KML_c,                   &
                                         Write3dFiles_c,ifm,ofm,nwt,intflg,wts)

      use io_units

      integer           ,intent(in) :: outunit
      character(len=1)  ,intent(in) :: WriteDepositFinal_ASCII_c
      character(len=1)  ,intent(in) :: WriteDepositFinal_KML_c
      character(len=1)  ,intent(in) :: WriteDepositTS_ASCII_c
      character(len=1)  ,intent(in) :: WriteDepositTS_KML_c
      character(len=1)  ,intent(in) :: WriteCloudConcentration_ASCII_c
      character(len=1)  ,intent(in) :: WriteCloudConcentration_KML_c
      character(len=1)  ,intent(in) :: WriteCloudHeight_ASCII_c
      character(len=1)  ,intent(in) :: WriteCloudHeight_KML_c
      character(len=1)  ,intent(in) :: WriteCloudLoad_ASCII_c
      character(len=1)  ,intent(in) :: WriteCloudLoad_KML_c
      character(len=1)  ,intent(in) :: WriteDepositTime_ASCII_c
      character(len=1)  ,intent(in) :: WriteDepositTime_KML_c
      character(len=1)  ,intent(in) :: WriteCloudTime_ASCII_c
      character(len=1)  ,intent(in) :: WriteCloudTime_KML_c
      character(len=1)  ,intent(in) :: Write3dFiles_c
      integer           ,intent(in) :: ifm
      integer           ,intent(in) :: ofm
      integer           ,intent(in) :: nwt
      logical           ,intent(in) :: intflg
      real(kind=dp),dimension(nwt),intent(in) :: wts

      if(outunit.gt.0)then
        write(outunit,1)&
         '******************* BLOCK 4 ****************************************************'
        if(WriteDepositFinal_ASCII_c      .eq.'n')write(outunit,11)'no '
        if(WriteDepositFinal_ASCII_c      .eq.'y')write(outunit,11)'yes'
        if(WriteDepositFinal_KML_c        .eq.'n')write(outunit,21)'no '
        if(WriteDepositFinal_KML_c        .eq.'y')write(outunit,21)'yes'
        if(WriteDepositTS_ASCII_c         .eq.'n')write(outunit,31)'no '
        if(WriteDepositTS_ASCII_c         .eq.'y')write(outunit,31)'yes'
        if(WriteDepositTS_KML_c           .eq.'n')write(outunit,41)'no '
        if(WriteDepositTS_KML_c           .eq.'y')write(outunit,41)'yes'
        if(WriteCloudConcentration_ASCII_c.eq.'n')write(outunit,51)'no '
        if(WriteCloudConcentration_ASCII_c.eq.'y')write(outunit,51)'yes'
        if(WriteCloudConcentration_KML_c  .eq.'n')write(outunit,61)'no '
        if(WriteCloudConcentration_KML_c  .eq.'y')write(outunit,61)'yes'
        if(WriteCloudHeight_ASCII_c       .eq.'n')write(outunit,71)'no '
        if(WriteCloudHeight_ASCII_c       .eq.'y')write(outunit,71)'yes'
        if(WriteCloudHeight_KML_c         .eq.'n')write(outunit,81)'no '
        if(WriteCloudHeight_KML_c         .eq.'y')write(outunit,81)'yes'
        if(WriteCloudLoad_ASCII_c         .eq.'n')write(outunit,91)'no '
        if(WriteCloudLoad_ASCII_c         .eq.'y')write(outunit,91)'yes'
        if(WriteCloudLoad_KML_c           .eq.'n')write(outunit,101)'no '
        if(WriteCloudLoad_KML_c           .eq.'y')write(outunit,101)'yes'
        if(WriteDepositTime_ASCII_c       .eq.'n')write(outunit,111)'no '
        if(WriteDepositTime_ASCII_c       .eq.'y')write(outunit,111)'yes'
        if(WriteDepositTime_KML_c         .eq.'n')write(outunit,121)'no '
        if(WriteDepositTime_KML_c         .eq.'y')write(outunit,121)'yes'
        if(WriteCloudTime_ASCII_c         .eq.'n')write(outunit,131)'no '
        if(WriteCloudTime_ASCII_c         .eq.'y')write(outunit,131)'yes'
        if(WriteCloudTime_KML_c           .eq.'n')write(outunit,141)'no '
        if(WriteCloudTime_KML_c           .eq.'y')write(outunit,141)'yes'
        if(Write3dFiles_c                 .eq.'n')write(outunit,151)'no ',ifm
        if(Write3dFiles_c                 .eq.'y')write(outunit,151)'yes',ifm
        if(ofm.eq.1)then
          write(outunit,161)'ascii '
        elseif(ofm.eq.2)then
          write(outunit,161)'binary'
        elseif(ofm.eq.3)then
          write(outunit,161)'netcdf'
        endif
        if(intflg)then
          ! Output steps are interval-based
          write(outunit,171)-1
          write(outunit,181)wts(1)
        else
          write(outunit,171)nwt
          write(outunit,181)wts(1:nwt)
        endif
        write(outunit,1)&
         '********************************************************************************'
      endif

 1    format(a80)
 11   format(a3,'     # Write out ESRI ASCII file of final deposit thickness?')
 21   format(a3,'     # Write out        KML file of final deposit thickness?')
 31   format(a3,'     # Write out ESRI ASCII deposit files at specified times?')
 41   format(a3,'     # Write out        KML deposit files at specified times?')
 51   format(a3,'     # Write out ESRI ASCII files of ash-cloud concentration?')
 61   format(a3,'     # Write out        KML files of ash-cloud concentration?')
 71   format(a3,'     # Write out ESRI ASCII files of ash-cloud height?')
 81   format(a3,'     # Write out        KML files of ash-cloud height?')
 91   format(a3,'     # Write out ESRI ASCII files of ash-cloud load (T/km2) at specified times?')
 101  format(a3,'     # Write out        KML files of ash-cloud load (T/km2) at specified times?')
 111  format(a3,'     # Write out ESRI ASCII file of deposit arrival times?')
 121  format(a3,'     # Write out        KML file of deposit arrival times?')
 131  format(a3,'     # Write out ESRI ASCII file of cloud arrival times?')
 141  format(a3,'     # Write out        KML file of cloud arrival times?')
 151  format(a3,2x,i1,'  # Write out 3-D ash concentration at specified times?')
 161  format(a6,3x,'# format of ash concentration files     (ascii, binary, or netcdf)')
 171  format(i3,6x,'# nWriteTimes')
 181  format(*(f8.3))

      end subroutine SetWrite_input_block_04

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  SetWrite_input_block_05
!
!  Called from: help_inputfile  
!  Arguments:
!    outunit    = output stream ID or 0 for Set Vars only with no output
!    nwindfiles = number of windfiles
!    windfiles  = names of windfiles
!
!  This subroutine writes the content of block 5 (Wind file list) of the Ash3d
!  control file.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SetWrite_input_block_05(outunit,                                                       &
                                         nwindfiles,                                                    &
                                         windfiles)

      use io_units

      integer           ,intent(in) :: outunit
      integer           ,intent(in) :: nwindfiles
      character(len=130),dimension(nwindfiles),intent(in) :: windfiles

      integer :: i

      if(outunit.gt.0)then
        write(outunit,1)&
         '******************* BLOCK 5 ****************************************************'
        do i=1,nwindfiles
          write(outunit,2)adjustl(windfiles(i))
        enddo
        write(outunit,1)&
         '********************************************************************************'
      endif

 1    format(a80)
 2    format(a130)

      end subroutine SetWrite_input_block_05

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  SetWrite_input_block_06
!
!  Called from: help_inputfile  
!  Arguments:
!    outunit                   = output stream ID or 0 for Set Vars only with no output
!    WriteAirportFile_ASCII_c  = Write out ash arrival times at airports to ASCII FILE?
!    WriteGSD_c                = Write out grain-size distribution to ASCII airport file?
!    WriteAirportFile_KML_c    = Write out ash arrival times to kml file?
!    b6l4                      = Name of file containing airport locations
!    ProjectAirportLocations_c = Defer to Lon/Lat coordinates? ("no" defers to projected)
!    
!  This subroutine writes the content of block 6 (Airport/POI info.) of the
!  Ash3d control file.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SetWrite_input_block_06(outunit,  &
                                         WriteAirportFile_ASCII_c,WriteGSD_c, &
                                         WriteAirportFile_KML_c,b6l4, &
                                         ProjectAirportLocations_c)

      use io_units

      integer           ,intent(in) :: outunit
      character(len=1)  ,intent(in) :: WriteAirportFile_ASCII_c
      character(len=1)  ,intent(in) :: WriteGSD_c
      character(len=1)  ,intent(in) :: WriteAirportFile_KML_c
      character(len=80) ,intent(in) :: b6l4
      character(len=1)  ,intent(in) :: ProjectAirportLocations_c


      if(outunit.gt.0)then
        write(outunit,1)&
         '******************* BLOCK 6 ****************************************************'
        if(WriteAirportFile_ASCII_c.eq.'n')write(outunit,11)'no '
        if(WriteAirportFile_ASCII_c.eq.'y')write(outunit,11)'yes'

        if(WriteGSD_c.eq.'n')write(outunit,21)'no '
        if(WriteGSD_c.eq.'y')write(outunit,21)'yes'

        if(WriteAirportFile_KML_c.eq.'n')write(outunit,31)'no '
        if(WriteAirportFile_KML_c.eq.'y')write(outunit,31)'yes'

        write(outunit,41)b6l4(1:30)

        if(ProjectAirportLocations_c.eq.'n')write(outunit,51)'no '
        if(ProjectAirportLocations_c.eq.'y')write(outunit,51)'yes'

        write(outunit,1)&
         '********************************************************************************'
      endif

 1    format(a80)
 11   format(a3,'     # Write out ash arrival times at airports to ASCII FILE?')
 21   format(a3,'     # Write out grain-size distribution to ASCII airport file?')
 31   format(a3,'     # Write out ash arrival times to kml file?')
 41   format(a30,'    # Name of file containing airport locations')
 51   format(a3,'     # Defer to Lon/Lat coordinates? ("no" defers to projected)')

      end subroutine SetWrite_input_block_06

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  SetWrite_input_block_07
!    outunit         = output stream ID or 0 for Set Vars only with no output
!
!  Called from: help_inputfile  
!  Arguments:
!    
!  This subroutine writes the content of block 7 (GSD specification) of the
!  Ash3d control file.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SetWrite_input_block_07(outunit                   ,&
                                         ns                        ,&  ! number of actual grain-size bins 
                                         fv_idx                    ,&  ! 
                                         shape_idx                 ,&  ! 
                                         mf                        ,&  ! 
                                         phim                      ,&  ! 
                                         phisig                    ,&  ! 
                                         T_diam                    ,&  ! 
                                         T_mf                      ,&  ! 
                                         T_fv                      ,&  !
                                         T_rho                     ,&  ! 
                                         T_F                       ,&  ! 
                                         T_G                       ,&  ! 
                                         T_phi)                        !

      use io_units

      use global_param,  only : &
         useLogNormGSbins

      use Tephra,        only : &
         n_gs_max,Tephra_gsdiam,Tephra_v_s,Tephra_bin_mass,Tephra_rho_m,FV_ID,&
         Tephra_gsF,Tephra_gsG,Tephra_gsPhi,Shape_ID,Tephra_Ncols,&
         LN_massfrac,LN_phi_mean,LN_phi_stddev,LN_suppl_frac

      integer                    ,intent(in) :: outunit
      integer                    ,intent(in) :: ns
      integer                    ,intent(in) :: fv_idx
      integer                    ,intent(in) :: shape_idx
      real(kind=ip)              ,intent(in) :: mf
      real(kind=ip)              ,intent(in) :: phim
      real(kind=ip)              ,intent(in) :: phisig
      real(kind=ip),dimension(ns),intent(in) :: T_diam
      real(kind=ip),dimension(ns),intent(in) :: T_mf
      real(kind=ip),dimension(ns),intent(in) :: T_fv
      real(kind=ip),dimension(ns),intent(in) :: T_rho
      real(kind=ip),dimension(ns),intent(in) :: T_F
      real(kind=ip),dimension(ns),intent(in) :: T_G
      real(kind=ip),dimension(ns),intent(in) :: T_phi

      integer :: isize

      n_gs_max        = ns
      FV_ID           = fv_idx
      Shape_ID        = shape_idx
      LN_massfrac     = mf
      LN_phi_mean     = phim
      LN_phi_stddev   = phisig
#ifdef USEPOINTERS
      if(.not.associated(Tephra_gsdiam))then
#else
      if(.not.allocated(Tephra_gsdiam))then
#endif
        allocate(Tephra_gsdiam(n_gs_max))
      endif
      Tephra_gsdiam(1:ns)   = T_diam(1:ns)
#ifdef USEPOINTERS
      if(.not.associated(Tephra_bin_mass))then
#else
      if(.not.allocated(Tephra_bin_mass))then
#endif
        allocate(Tephra_bin_mass(n_gs_max))
      endif
      Tephra_bin_mass(1:ns) = T_mf(1:ns)
      if(useLogNormGSbins)then
        ! If a Log-normal supplement was used, recover the original list.
        ! Note, this may have been sorted by size and in a different order than originally
        ! given.
        do isize = 1,n_gs_max
          Tephra_bin_mass(isize) = Tephra_bin_mass(isize) - LN_suppl_frac(isize)*LN_massfrac
        enddo
      endif

#ifdef USEPOINTERS
      if(.not.associated(Tephra_rho_m))then
#else
      if(.not.allocated(Tephra_rho_m))then
#endif
        allocate(Tephra_rho_m(n_gs_max))
      endif
      Tephra_rho_m(1:ns)    = T_rho(1:ns)
#ifdef USEPOINTERS
      if(.not.associated(Tephra_gsF))then
#else
      if(.not.allocated(Tephra_gsF))then
#endif
        allocate(Tephra_gsF(n_gs_max))
      endif
      Tephra_gsF(1:ns)      = T_F(1:ns)
#ifdef USEPOINTERS
      if(.not.associated(Tephra_gsG))then
#else
      if(.not.allocated(Tephra_gsG))then
#endif
        allocate(Tephra_gsG(n_gs_max))
      endif
      Tephra_gsG(1:ns)      = T_G(1:ns)
#ifdef USEPOINTERS
      if(.not.associated(Tephra_gsPhi))then
#else
      if(.not.allocated(Tephra_gsPhi))then
#endif
        allocate(Tephra_gsPhi(n_gs_max))
      endif
      Tephra_gsPhi(1:ns)    = T_phi(1:ns)

      ! Finally, write out block
      if(outunit.gt.0)then
        write(outunit,1)&
         '******************* BLOCK 7 ****************************************************'
        if(useLogNormGSbins)then
          ! For the log-normal supplemental line, we need to bump n_gs_max by one
          write(outunit,11)n_gs_max+1,FV_ID,Shape_ID
        else
          write(outunit,11)n_gs_max,FV_ID,Shape_ID
        endif
        if(Tephra_Ncols.eq.2)then
          ! 2-columns is just: fall-velocity (m/s), mass-fraction
          do isize = 1,n_gs_max
            write(outunit,12)Tephra_v_s(isize),Tephra_bin_mass(isize)
          enddo        
        elseif(Tephra_Ncols.eq.3)then
          ! 3-columns is: diameter (mm), mass fraction, density (kg/m3)
          do isize = 1,n_gs_max
            write(outunit,13)Tephra_gsdiam(isize),Tephra_bin_mass(isize),Tephra_rho_m(isize)
          enddo
        elseif(Tephra_Ncols.eq.4)then
          ! 4-columns is: diameter (mm), mass fraction, density (kg/m3), Shape F (or Psi)
          do isize = 1,n_gs_max  
            if(Shape_ID.eq.2)then
              write(outunit,14)Tephra_gsdiam(isize),Tephra_bin_mass(isize),Tephra_rho_m(isize),&
                               Tephra_gsPhi(isize)
            else
              write(outunit,14)Tephra_gsdiam(isize),Tephra_bin_mass(isize),Tephra_rho_m(isize),&
                               Tephra_gsF(isize)
            endif
          enddo
        elseif(Tephra_Ncols.eq.5)then
          ! 5-columns is: diameter (mm), mass fraction, density (kg/m3), Shape F, G
          do isize = 1,n_gs_max
            write(outunit,15)Tephra_gsdiam(isize),Tephra_bin_mass(isize),Tephra_rho_m(isize),&
                             Tephra_gsF(isize),Tephra_gsG(isize)
          enddo
        else
          ! Number of columns of original control file not recorded in this NetCDF file, assume 4
          do isize = 1,n_gs_max
            if(Shape_ID.eq.2)then
              write(outunit,14)Tephra_gsdiam(isize),Tephra_bin_mass(isize),Tephra_rho_m(isize),&
                               Tephra_gsPhi(isize)
            else
              write(outunit,14)Tephra_gsdiam(isize),Tephra_bin_mass(isize),Tephra_rho_m(isize),&
                               Tephra_gsF(isize)
            endif
          enddo
        endif
  
        if(useLogNormGSbins)then
          write(outunit,16)-1,LN_phi_mean,LN_phi_stddev
        endif
        write(outunit,1)&
         '********************************************************************************'
      endif

 1    format(a80)
 11   format(i2,2x,i1,2x,i1,24x,'# Number of grain-size bins. FV_ID not given; defaults to 1')
 12   format(f8.5,1x,f8.5)
 13   format(f8.5,1x,f8.5,1x,f8.1)
 14   format(f8.5,1x,f8.5,1x,f8.1,1x,f8.3)
 15   format(f8.5,1x,f8.5,1x,f8.1,1x,f8.3,1x,f8.3)
 16   format(i2,1x,f8.5,1x,f8.1)

      end subroutine SetWrite_input_block_07

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  SetWrite_input_block_08
!
!  Called from: help_inputfile  
!  Arguments:
!    outunit         = output stream ID or 0 for Set Vars only with no output
!    
!  This subroutine writes the content of block 8 (Vertical profile info.) of
!  the Ash3d control file.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SetWrite_input_block_08(outunit,nprof,x_prof,y_prof,name_prof)

      use io_units

      use io_data,       only : &
         nvprofiles,Site_vprofile,x_vprofile, y_vprofile

      integer                           ,intent(in) :: outunit
      integer                           ,intent(in) :: nprof
      real(kind=ip)    ,dimension(nprof),intent(in) :: x_prof
      real(kind=ip)    ,dimension(nprof),intent(in) :: y_prof
      character(len=50),dimension(nprof),intent(in) :: name_prof

      integer :: iprof

      nvprofiles   = nprof

#ifdef USEPOINTERS
      if(.not.associated(x_vprofile))then
#else
      if(.not.allocated(x_vprofile))then
#endif
        allocate(x_vprofile(nvprofiles))
      endif
      x_vprofile(1:nvprofiles) = x_prof(1:nvprofiles)

#ifdef USEPOINTERS
      if(.not.associated(y_vprofile))then
#else
      if(.not.allocated(y_vprofile))then
#endif
        allocate(y_vprofile(nvprofiles))
      endif
      y_vprofile(1:nvprofiles) = y_prof(1:nvprofiles)

#ifdef USEPOINTERS
      if(.not.associated(Site_vprofile))then
#else
      if(.not.allocated(Site_vprofile))then
#endif
        allocate(Site_vprofile(nvprofiles))
      endif
      Site_vprofile(1:nvprofiles) = name_prof(1:nvprofiles)

      ! Finally, write out block
      if(outunit.gt.0)then
        write(outunit,1)&
         '******************* BLOCK 8 ****************************************************'
          write(outunit,11)nvprofiles
          do iprof = 1,nvprofiles
            write(outunit,12)x_vprofile(iprof),y_vprofile(iprof),Site_vprofile(iprof)
          enddo
        write(outunit,1)&
         '********************************************************************************'
      endif

 1    format(a80)
 11   format(i3,27x,'# number of locations for vertical profiles (nlocs)')
 12   format(f13.5,2x,f13.5,2x,a50)

      end subroutine SetWrite_input_block_08

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  SetWrite_input_block_09
!
!  Called from: help_inputfile  
!  Arguments:
!    outunit         = output stream ID or 0 for Set Vars only with no output
!    
!  This subroutine writes the content of block 9 (NetCDF annotations) of the
!  Ash3d control file.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SetWrite_input_block_09(outunit,ofilename,title,comment)

      use io_units

      use io_data,       only : &
         concenfile,cdf_title,cdf_comment

      integer           ,intent(in) :: outunit
      character(len=80) ,intent(in) :: ofilename
      character(len=130),intent(in) :: title
      character(len=80) ,intent(in) :: comment

      concenfile  = ofilename
      cdf_title   = title
      cdf_comment = comment

      if(outunit.gt.0)then
        write(outunit,1)&
         '******************* BLOCK 9 ****************************************************'
        write(outunit,1)adjustl(concenfile)
        write(outunit,2)adjustl(cdf_title)
        write(outunit,1)adjustl(cdf_comment)
        write(outunit,1)&
         '********************************************************************************'
      endif

 1    format(a80)
 2    format(a130)

      end subroutine SetWrite_input_block_09

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  SetWrite_input_block_ResetParam
!
!  Called from: help_inputfile  
!  Arguments:
!    outunit         = output stream ID or 0 for Set Vars only with no output
!    
!  This subroutine writes the content of block 10+ (OPTMOD=RESETPARAMS) of the
!  Ash3d control file.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SetWrite_input_block_ResetParam(outunit)

      use io_units

      use global_param,  only : &
         GRAV,GRAV_Default,RAD_EARTH,RAD_EARTH_Default,CFL,CFL_Default,DT_MIN,DT_MIN_Default,&
         DT_MAX,DT_MAX_Default,useVz_rhoG,useVz_rhoG_Default,useMoistureVars,useMoistureVars_Default

      use io_data,       only : &
         cdf_institution,cdf_institution_Default,cdf_run_class,cdf_run_class_Default,&
         cdf_url,cdf_url_Default

      use mesh,          only : &
         ZPADDING,ZPADDING_Default

      use Tephra,        only : &
         MagmaDensity,MagmaDensity_Default,DepositDensity,DepositDensity_Default,&
         LAM_GS_THRESH,LAM_GS_THRESH_Default,AIRBORNE_THRESH,AIRBORNE_THRESH_Default

      use Output_Vars,   only : &
         DEPO_THRESH,DEPO_THRESH_Default,DEPRATE_THRESH,DEPRATE_THRESH_Default,&
         CLOUDCON_THRESH,CLOUDCON_THRESH_Default,CLOUDCON_GRID_THRESH,CLOUDCON_GRID_THRESH_Default,&
         CLOUDLOAD_THRESH,CLOUDLOAD_THRESH_Default,THICKNESS_THRESH,THICKNESS_THRESH_Default,&
         DBZ_THRESH,DBZ_THRESH_Default,useWindVars,useWindVars_Default,useOutprodVars,useOutprodVars_Default,&
         useRestartVars,useRestartVars_Default

      use Diffusion,     only : &
         Imp_fac,Imp_fac_Default,Imp_DT_fac,Imp_DT_fac_Default

      use Source_Umbrella,only : &
         VelMod_umb,VelMod_umb_Default,k_entrainment_umb,k_entrainment_umb_Default,&
         lambda_umb,lambda_umb_Default,N_BV_umb,N_BV_umb_Default,&
         SuzK_umb,SuzK_umb_Default

      integer           ,intent(in) :: outunit

      real(kind=ip) :: tmp1,tmp2,rerr
      integer       :: dumint1,dumint2

      if(outunit.gt.0)then
        write(outunit,1)&
         '******************* BLOCK 10+ **************************************************'
        write(outunit,'(a18)')'OPTMOD=RESETPARAMS'
        ! Parameters from Tephra
        tmp1=MagmaDensity
        tmp2=MagmaDensity_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'MagmaDensity        ',tmp1
        tmp1=DepositDensity
        tmp2=DepositDensity_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'DepositDensity      ',tmp1
        tmp1=LAM_GS_THRESH
        tmp2=LAM_GS_THRESH_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'LAM_GS_THRESH       ',tmp1
        tmp1=AIRBORNE_THRESH
        tmp2=AIRBORNE_THRESH_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'AIRBORNE_THRESH     ',tmp1
        ! Parameters from global_param
        tmp1=GRAV
        tmp2=GRAV_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'GRAV                ',tmp1
        tmp1=CFL
        tmp2=CFL_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'CFL                 ',tmp1
        tmp1=RAD_EARTH
        tmp2=RAD_EARTH_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'RAD_EARTH           ',tmp1
        tmp1=DT_MIN
        tmp2=DT_MIN_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'DT_MIN              ',tmp1
        tmp1=DT_MAX
        tmp2=DT_MAX_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'DT_MAX               ',tmp1
        if(useVz_rhoG)then
          dumint1 = 1
        else
          dumint1 = 0
        endif
        if(useVz_rhoG.neqv.useVz_rhoG_Default) &
                          write(outunit,3)'useVz_rhoG           ',dumint1
        if(useMoistureVars)then
          dumint1 = 1
        else
          dumint1 = 0
        endif
        if(useMoistureVars.neqv.useMoistureVars_Default) &
                          write(outunit,3)'useMoistureVars      ',dumint1
        ! Parameters from mesh
        tmp1=ZPADDING
        tmp2=ZPADDING_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'ZPADDING             ',tmp1

        ! Parameters from Output_Vars
        tmp1=DEPO_THRESH
        tmp2=DEPO_THRESH_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'DEPO_THRESH         ',tmp1
        tmp1=DEPRATE_THRESH
        tmp2=DEPRATE_THRESH_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'DEPRATE_THRESH      ',tmp1
        tmp1=CLOUDCON_THRESH
        tmp2=CLOUDCON_THRESH_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'CLOUDCON_THRESH     ',tmp1
        tmp1=CLOUDCON_GRID_THRESH
        tmp2=CLOUDCON_GRID_THRESH_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'CLOUDCON_GRID_THRESH',tmp1
        tmp1=CLOUDLOAD_THRESH
        tmp2=CLOUDLOAD_THRESH_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'CLOUDLOAD_THRESH    ',tmp1
        tmp1=THICKNESS_THRESH
        tmp2=THICKNESS_THRESH_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'THICKNESS_THRESH    ',tmp1
        tmp1=DBZ_THRESH
        tmp2=DBZ_THRESH_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'DBZ_THRESH          ',tmp1
        if(useWindVars)then
          dumint1 = 1
        else
          dumint1 = 0
        endif
        if(useWindVars.neqv.useWindVars_Default) &
                          write(outunit,3)'useWindVars         ',dumint1
        if(useOutprodVars)then
          dumint1 = 1
        else
          dumint1 = 0
        endif
        if(useOutprodVars.neqv.useOutprodVars_Default) &
                          write(outunit,3)'useOutprodVars       ',dumint1
        if(useRestartVars)then
          dumint1 = 1
        else
          dumint1 = 0
        endif
        if(useRestartVars.neqv.useRestartVars_Default) &
                          write(outunit,3)'useRestartVars       ',dumint1

        ! Parameters from Diffusion
        tmp1=Imp_fac
        tmp2=Imp_fac_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'Imp_fac              ',tmp1
        tmp1=Imp_DT_fac
        tmp2=Imp_DT_fac_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'Imp_DT_fac          ',tmp1

        ! Parameters from Source_Umbrella
        if(VelMod_umb.ne.VelMod_umb_Default) &
                          write(outunit,3)'VelMod_umb           ',tmp1
        tmp1=k_entrainment_umb
        tmp2=k_entrainment_umb_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'k_entrainment_umb    ',tmp1
        tmp1=lambda_umb
        tmp2=lambda_umb_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'lambda_umb           ',tmp1
        tmp1=N_BV_umb
        tmp2=N_BV_umb_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'N_BV_umb             ',tmp1
        tmp1=SuzK_umb
        tmp2=SuzK_umb_Default
        rerr=abs(tmp1-tmp2)/tmp1
        if(rerr.gt.1.0e-3)write(outunit,2)'SuzK_umb             ',tmp1

        ! Parameters from io_data
        dumint1=index(cdf_institution,cdf_institution_Default)
        dumint2=index(cdf_institution_Default,cdf_institution)
        if(dumint1.gt.0.and.dumint2.gt.0.and.dumint1.ne.dumint2)&
                          write(outunit,4)'cdf_institution      ',trim(adjustl(cdf_institution))
        dumint1=index(cdf_run_class,cdf_run_class_Default)
        dumint2=index(cdf_run_class_Default,cdf_run_class)
        if(dumint1.gt.0.and.dumint2.gt.0.and.dumint1.ne.dumint2)&
                          write(outunit,4)'cdf_run_class        ',trim(adjustl(cdf_run_class))

        dumint1=index(cdf_url,cdf_url_Default)
        dumint2=index(cdf_url_Default,cdf_url)
        if(dumint1.gt.0.and.dumint2.gt.0.and.dumint1.ne.dumint2)&
                          write(outunit,4)'cdf_url              ',trim(adjustl(cdf_url))

        write(outunit,1)&
         '********************************************************************************'
      endif

 1    format(a80)
 2    format(a20,' = ',g10.3)
 3    format(a20,' = ',i1)
 4    format(a20,' = ',g0)

      end subroutine SetWrite_input_block_ResetParam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  SetWrite_input_block_Topo
!
!  Called from: help_inputfile  
!  Arguments:
!    outunit         = output stream ID or 0 for Set Vars only with no output
!    
!  This subroutine writes the content of block 10+ (OPTMOD=TOPO) of the Ash3d
!  control file.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SetWrite_input_block_Topo(outunit)

      use io_units

      integer           ,intent(in) :: outunit

      end subroutine SetWrite_input_block_Topo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  SetWrite_input_block_VarDiff
!
!  Called from: help_inputfile  
!  Arguments:
!    outunit         = output stream ID or 0 for Set Vars only with no output
!    
!  This subroutine writes the content of block 10+ (OPTMOD=VARDIFF) of the
!  Ash3d control file.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SetWrite_input_block_VarDiff(outunit)

      use io_units

      integer           ,intent(in) :: outunit

      end subroutine SetWrite_input_block_VarDiff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module Ash3d_Program_Control

!##############################################################################

