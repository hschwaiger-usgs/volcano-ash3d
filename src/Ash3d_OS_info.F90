!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Subroutine Set_OS_Env 
!
!  This subroutine evaluates the state of all aspect of the run.
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
!  Next, details of the system state are logged, including OS type (Linux, Mac, Windows),
!  endian flavor of hardware, fortran compiler version and flags, command-line arguments,
!  and date/time of the run.  Additionally, if PII=ON was set in the makefile when this
!  executable was compiled, then the user, system hostname and run directory are also logged.
!
!  Finally, all pre-processor flags are checked here with logging to stdout of which flags
!  are invoked. Pre-processor flags checked:
!   LINUX, MACOS, WINDOWS  : OS declaration
!   FAST_DT, FAST_SUBGRID  : Speed-up tools for dt and grid calculations
!   EXPLDIFF, CRANKNIC     : Explicit vs implicit diffusion algorithm
!   LIM_NONE,LIM_LAXWEN,LIM_BW,LIM_FROMM,LIM_MINMOD,LIM_SUPERBEE,LIM_MC  : Limiter for advection
!   USENETCDF, USEGRIB     : Invokes netcdf and/or grib functionality
!   USEPOINTERS            : determines if variables are pointers or allocatable arrays
!   USEEXTDATA             : determines if Ash3d will rely on external lists for airport and volcano data
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Set_OS_Env

      use precis_param

      use io_units

      use global_param,  only : &
        DirPrefix,DirDelim,IsLitEnd,IsLinux,IsWindows,IsMacOS,version, &
        CFL,OS_TYPE,OS_Flavor,os_full_command_line,os_cwd,os_host,os_user

      use io_data,       only : &
        Ash3dHome

      use time_data,     only : &
        BaseYear,useLeap,os_time_log, &
        RunStartDay,RunStartHour_ch,RunStartHr,RunStartMinute,RunStartMonth,RunStartYear

      use MetReader,     only : &
         MR_OS_TYPE,MR_DirPrefix,MR_DirDelim,MR_VERB,MR_nio

      ! This module requires Fortran 2003 or later
      use iso_fortran_env, only : &
           compiler_version,&
           compiler_options

      implicit none

      integer              :: iostatus,ierr
      character(len=130)   :: tmp_str
        ! variables to hold results of date_and_time
      character(len=8)  :: date
      character(len=10) :: time2
      character(len=5)  :: zone
      integer           :: values(8)
      integer           :: timezone
      real(kind=dp)     :: StartHour
      real(kind=dp)     :: RunStartHour    ! Start time of model run, in hours since BaseYear
      character(len=80) :: linebuffer080
      character(len=100):: CompVer
      character(len=508):: CompOpt
      logical           :: IsThere

      INTERFACE
        real(kind=8) function HS_hours_since_baseyear(iyear,imonth,iday,hours,byear,useLeaps)
          integer            :: iyear
          integer            :: imonth
          integer            :: iday
          real(kind=8)       :: hours
          integer            :: byear
          logical            :: useLeaps
        end function HS_hours_since_baseyear
        character (len=13) function HS_yyyymmddhhmm_since(HoursSince,byear,useLeaps)
          real(kind=8)               ::  HoursSince
          integer                    ::  byear
          logical                    ::  useLeaps
        end function HS_yyyymmddhhmm_since
        character (len=20) function HS_xmltime(HoursSince,byear,useLeaps)
          real(kind=8)              :: HoursSince
          integer                   :: byear
          logical                   :: useLeaps
        end function HS_xmltime
        subroutine check_endian(IsLitEnd)
          logical,intent(inout)  :: IsLitEnd
        end subroutine check_endian
!#ifdef USEPII
!        ! These are included in gfortran, but are not standard
!        subroutine getlog(os_user)
!          character(len=32) :: os_user
!        end subroutine getlog
!        subroutine hostnm(os_host)
!          character(*) :: os_host
!        end subroutine hostnm
!        subroutine getcwd(os_cwd)
!          character(len=255) :: os_cwd
!        end subroutine getcwd
!#endif
      END INTERFACE

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

      ! ENVIRONMENT VARIABLES
      ! First order of business is to get the environment variable (if present) for verbosity
      call get_environment_variable(name="ASH3DVERB",VALUE=tmp_str,STATUS=iostatus)
      if(iostatus.eq.0)then
        ! read the value of the ASH3DVERB environment variable to the local variable VB(1)
        read(tmp_str,*,iostat=ierr)VB(1)
        if(ierr.ne.0)then
          write(errlog(1),*)"ERROR: ASH3DVERB found, but expecting an integer value"
          write(errlog(1),*)"       Instead, env. variable set to: ",tmp_str
          stop 1
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
        write(errlog(1),*)"ERROR: Verbosity level not recognized. Value should be between 1 and 10."
        write(errlog(1),*)"    verbosity level : ",VB(1)
        stop 1
      endif
      write(outlog(1),*)"    verbosity level : ",VB(1), vlevel

      if(VB(1).eq.verbosity_dark)then
        ! If the output verbosity is for nothing at all, reset the log verbosity to that too
        VB(2) = verbosity_dark
      else
         ! Before we do anything else, start a log file
        open(unit=global_log,file=logfile,status='unknown')
        ! Mirror the above output to the logfile
        if(iostatus.eq.0)then
          if(VB(2).le.verbosity_info)then
            write(outlog(2),*)"Checking for run-time environment variable: ASH3DVERB"
            write(outlog(2),*)"  Verbosity reset by environment variable to: ",VB(1)
          endif
        else
          if(VB(2).le.verbosity_info)then
            write(outlog(2),*)"Checking for run-time environment variable: ASH3DVERB"
            write(outlog(2),*)"  ASH3DVERB environment variable not found.  verbosity level : ",VB(1)
          endif
        endif
      endif
      ! Harmonizing verbosity levels with MetReader
      MR_VERB = VB(1)
      MR_nio  = 2    ! Ash3d uses a logfile so set the output streams to stdin/stderr + logfile

      ! Next, check for environment variables ASH3DHOME

      ! Set the default installation path
      ! This is only needed if shared data files with fixed paths are read
      ! in such as the global airport and volcano ESP files.
      Ash3dHome = trim(adjustl(DirPrefix)) // DirDelim // &
                  "opt" // DirDelim // "USGS" // DirDelim // "Ash3d"

      ! Here it is over-written by compile-time path, if available
#include "installpath.h"
      ! This can be over-written if an environment variable is set
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Checking for run-time environment variable: ASH3DHOME"
      endif;enddo

      call get_environment_variable(name="ASH3DHOME",VALUE=tmp_str,STATUS=iostatus)
      if(iostatus.eq.0)then
        Ash3dHome = tmp_str
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"  Install path reset by environment variable to: ",trim(adjustl(Ash3dHome))
        endif;enddo
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"  ASH3DHOME not found. Install path set to: ",trim(adjustl(Ash3dHome))
        endif;enddo
      endif
      inquire(file=Ash3dHome,exist=IsThere)
      if(IsThere)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"  Path to ASH3DHOME found on system.  Good."
        endif;enddo
      else
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Cannot find ASH3DHOME=",Ash3dHome
        endif;enddo
        stop 1
      endif

      ! Now, check for environment variables ASH3DCFL
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Checking for run-time environment variable: ASH3DCFL"
      endif;enddo
      call get_environment_variable(name="ASH3DCFL",VALUE=tmp_str,STATUS=iostatus)
      if(iostatus.eq.0)then
        read(tmp_str,*,iostat=ierr)CFL
        if(ierr.ne.0)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: ASH3DCFL found, but expecting a floating point value"
            write(errlog(io),*)"       Instead, env. variable set to: ",tmp_str
          endif;enddo
          stop 1
        endif
        if (CFL.le.0.0_ip)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: CFL must be > 0"
            write(errlog(io),*)"       Currently set to ",CFL
          endif;enddo
          stop 1
        elseif (CFL.ge.1.0_ip)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: CFL must be < 1"
            write(errlog(io),*)"       Currently set to ",CFL
          endif;enddo
          stop 1
        endif
        ! variable seems valid, proceeding.
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"  CFL condition reset by environment variable to: ",real(CFL,kind=4)
          write(outlog(io),*)"   Note: It is possible this may be subsequently reset via"
          write(outlog(io),*)"         the control file in a RESETPARAM block. Check the log file"
          write(outlog(io),*)"         or the netcdf output file for final CFL used."
        endif;enddo
      else
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"  ASH3DCFL not found.  CFL condition : ",real(CFL,kind=4)
        endif;enddo
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

#ifdef USEPII
      call getlog(os_user)
      call hostnm(os_host)
      call getcwd(os_cwd)
#else
      os_user = 'N/A'
      os_host = 'N/A'
      os_cwd  = 'N/A'
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
        write(outlog(io),*)" "
        write(outlog(io),*)"Running Ash3d with command-line: ",&
                    trim(adjustl(os_full_command_line))
      endif;enddo

      ! Determining the run start time
      read(zone,'(i3)') timezone
      ! FIND TIME IN UTC
      StartHour = real(values(5)-timezone,kind=ip) + real(values(6)/60.0,kind=ip)    !add offset to UTC
        ! find time in HoursSinceBaseYear
        !  Note: This will be relative to the BaseYear in time_data (default is 1900). 
        !        That BaseYear might be changed if the eruption start time is
        !        before BaseYear
      RunStartHour    = HS_hours_since_baseyear(values(1),values(2),values(3), &
                                                StartHour,BaseYear,useLeap)
      RunStartHour_ch = HS_yyyymmddhhmm_since(RunStartHour,BaseYear,useLeap)
      read(RunStartHour_ch,'(i4)') RunStartYear
      read(RunStartHour_ch,'(4x,i2)') RunStartMonth
      read(RunStartHour_ch,'(6x,i2)') RunStartDay
      read(RunStartHour_ch,'(8x,i2)') RunStartHr
      read(RunStartHour_ch,'(11x,i2)') RunStartMinute

        ! Prepare a note to include in the netcdf output file
!      write(linebuffer080,102) RunStartYear,RunstartMonth,RunStartDay,RunStartHr,RunStartMinute
!      os_time_log = linebuffer080(1:17)
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
  
        write(outlog(io),*)"  "
        write(outlog(io),*)"This executable was compiled with the following compiler and options:"
        write(outlog(io),*)"    ",trim(adjustl(CompVer))
        write(outlog(io),*)"    ",trim(adjustl(CompOpt))
        write(outlog(io),*)"and with the following pre-proc flags:"
#ifdef LINUX
        write(outlog(io),*)"         LINUX: System specified as linux"
#endif
#ifdef MACOS
        write(outlog(io),*)"         MACOS: System specified as MacOS"
#endif
#ifdef WINDOWS
        write(outlog(io),*)"       WINDOWS: System specified as MS Windows"
#endif
#ifdef FAST_DT
        write(outlog(io),*)"       FAST_DT: ON"
        write(outlog(io),*)"              : dt will only be evaluated on the time steps"
        write(outlog(io),*)"                in the wind files.  If there are processes"
        write(outlog(io),*)"                that affect the wind speeds (e.g. umbrella"
        write(outlog(io),*)"                spreading), this can cause job failure."
#else
        write(outlog(io),*)"       FAST_DT: OFF"
        write(outlog(io),*)"              : dt will be evaluated each time step"
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
        write(outlog(io),*)"      CRANKNIC: Diffusion will be calculated via Crank-Nicolson"
#endif
#ifdef LIM_NONE
        write(outlog(io),*)"      LIM_NONE: Advection routines use no limiters"
#endif
#ifdef LIM_LAXWEN
        write(outlog(io),*)"    LIM_LAXWEN: Advection routines use a Lax-Wendrof limiter"
#endif
#ifdef LIM_BW
        write(outlog(io),*)"        LIM_BW: Advection routines use a Beam-Warming limiter"
#endif
#ifdef LIM_FROMM
        write(outlog(io),*)"     LIM_FROMM: Advection routines use a Fromm limiter"
#endif
#ifdef LIM_MINMOD
        write(outlog(io),*)"    LIM_MINMOD: Advection routines use a MinMod limiter"
#endif
#ifdef LIM_SUPERBEE
        write(outlog(io),*)"  LIM_SUPERBEE: Advection routines use a SuperBee limiter"
#endif
#ifdef LIM_MC
        write(outlog(io),*)"        LIM_MC: Advection routines use a MC limiter"
#endif
#ifdef USENETCDF
        write(outlog(io),*)"     USENETCDF: NetCDF functionality is included"
#else
        write(outlog(io),*)"     USENETCDF not found: NetCDF functionality is not included"
#endif
#ifdef USEGRIB
        write(outlog(io),*)"       USEGRIB: Grib functionality is included"
#else
        write(outlog(io),*)"       USEGRIB not found: Grib functionality is not included"
#endif
#ifdef USEPOINTERS
        write(outlog(io),*)"   USEPOINTERS: Arrays are defined as pointers"
        write(outlog(io),*)"                This helps Ash3d subroutines to be called via C++"
#else
        write(outlog(io),*)"   USEPOINTERS not found: All arrays are allocatable."
#endif
#ifdef USEEXTDATA
        write(outlog(io),*)"   USEEXTDATA: Data files for airports and volcanoes are"
        write(outlog(io),*)"                read at run-time"
#else
        write(outlog(io),*)"   USEEXTDATA not found: airport and volcano list included when"
        write(outlog(io),*)"                          compiled."
#endif
  
        ! WRITE OUT START TIME IN UTC
        write(outlog(io),*)
        write(outlog(io),2) version,RunStartYear,RunStartMonth,RunStartDay,RunStartHr,RunStartMinute
        write(outlog(io),*)
      endif;enddo

      ! FORMAT STATEMENTS
2     format(4x,'Ash3d (Rev ',a5,') run ',&
             i4,'.',i2.2,'.',i2.2,i4,':',i2.2,' UTC')
102   format(i4,'.',i2.2,'.',i2.2,i4,':',i2.2)

      end subroutine Set_OS_Env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Subroutine CHECK_ENDIAN checks if the local system uses Big-Endian
!  or Little-Endian byte ordering.  Returns the logical value
!  IsLitEnd = .true. if the system is Little-Endian, .false. otherwise.
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine check_endian(IsLitEnd)

      implicit none

      logical,intent(inout)  :: IsLitEnd

      integer(kind=2)  :: s = 1

      if(btest(transfer(int((/1,0/),kind=1),s),0)) then
        ! System is Little-Endian
        IsLitEnd = .true.
      else
        ! System is Big-Endian
        IsLitEnd = .false.
      endif 

      end subroutine check_endian

