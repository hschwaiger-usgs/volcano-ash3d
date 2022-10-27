      subroutine Set_OS_Env

      use precis_param

      use io_units

      use global_param,  only : &
        DirPrefix,DirDelim,IsLinux,IsWindows,IsMacOS, &
        CFL,OS_TYPE,OS_Flavor,os_full_command_line,os_cwd,os_host,os_user

      use io_data,       only : &
        Ash3dHome

      use time_data,     only : &
        BaseYear,useLeap,os_time_log, &
        RunStartDay,RunStartHour_ch,RunStartHr,RunStartMinute,RunStartMonth,RunStartYear

      implicit none

      integer              :: iostatus
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
      character(len=8)  :: version             =  ' 1.0  '

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
      END INTERFACE

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

      ! Get some run-specific and system-specific information
      call get_command(os_full_command_line)
      call getlog(os_user)
      call hostnm(os_host)
      call getcwd(os_cwd)
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

      write(global_info,*)" "
      write(global_info,*)"Running Ash3d with command-line: ",adjustl(trim(os_full_command_line))

      ! Check for environment variables ASH3DHOME and ASH3DCFL

      ! Set the default installation path
      ! This is only needed if shared data files with fixed paths are read
      ! in such as the global airport and volcano ESP files.
      Ash3dHome = adjustl(trim(DirPrefix)) // DirDelim // &
                  "opt" // DirDelim // "USGS" // DirDelim // "Ash3d"

      ! Here it is over-written by compile-time path, if available
#include "installpath.h"
      ! This can be over-written if an environment variable is set
      write(global_info,*)" "
      write(global_info,*)"Checking for run-time environment variable: ASH3DHOME"
      call GET_ENVIRONMENT_VARIABLE(NAME="ASH3DHOME",VALUE=tmp_str,STATUS=iostatus)
      if(iostatus.eq.0)then
        Ash3dHome = tmp_str
        write(global_info,*)&
          "  Install path reset by environment variable to: ",adjustl(trim(Ash3dHome))
      else
        write(global_info,*)&
          "  ASH3DHOME not found. Install path set to: ",adjustl(trim(Ash3dHome))
      endif
      write(global_info,*)"Checking for run-time environment variable: ASH3DCFL"
      call GET_ENVIRONMENT_VARIABLE(NAME="ASH3DCFL",VALUE=tmp_str,STATUS=iostatus)
      if(iostatus.eq.0)then
        read(tmp_str,*)CFL
        write(global_info,*)&
          "  CFL condition reset by environment variable to: ",real(CFL,kind=4)
      else
        write(global_info,*)&
          "  ASH3DCFL not found.  CFL condition : ",real(CFL,kind=4)
      endif

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

      write(global_info,*)" System Information"
      write(global_info,*)"   host: ",adjustl(trim(os_host)),' (',adjustl(trim(OS_Flavor)),')'
      write(global_info,*)"    cwd: ",adjustl(trim(os_cwd))
      write(global_info,*)"   user: ",adjustl(trim(os_user))

      write(global_info,*)"  "
      write(global_info,*)"This executable was compiled with the following pre-proc flags:"
#ifdef LINUX
      write(global_info,*)"         LINUX: System specified as linux"
#endif
#ifdef MACOS
      write(global_info,*)"         MACOS: System specified as MacOS"
#endif
#ifdef WINDOWS
      write(global_info,*)"       WINDOWS: System specified as MS Windows"
#endif
#ifdef FAST_DT
      write(global_info,*)"       FAST_DT: dt will only be evaluated on the time steps"
      write(global_info,*)"                in the wind files.  If there are processes"
      write(global_info,*)"                that affect the wind speeds (e.g. umbrella"
      write(global_info,*)"                spreading), this can cause job failure."
#endif
#ifdef FAST_SUBGRID
      write(global_info,*)"  FAST_SUBGRID: Advection and diffusion routines will only"
      write(global_info,*)"                be calculated in the region where the cloud"
      write(global_info,*)"                concentration exceeds a threshold"
#endif
#ifdef EXPLDIFF
      write(global_info,*)"      EXPLDIFF: Diffusion will be calculated via the explicit solver."
#endif
#ifdef CRANKNIC
      write(global_info,*)"      CRANKNIC: Diffusion will be calculated via Crank-Nicolson"
#endif
#ifdef LIM_NONE
      write(global_info,*)"      LIM_NONE: Advection routines use no limiters"
#endif
#ifdef LIM_LAXWEN
      write(global_info,*)"    LIM_LAXWEN: Advection routines use a Lax-Wendrof limiter"
#endif
#ifdef LIM_BW
      write(global_info,*)"        LIM_BW: Advection routines use a Beam-Warming limiter"
#endif
#ifdef LIM_FROMM
      write(global_info,*)"     LIM_FROMM: Advection routines use a Fromm limiter"
#endif
#ifdef LIM_MINMOD
      write(global_info,*)"    LIM_MINMOD: Advection routines use a MinMod limiter"
#endif
#ifdef LIM_SUPERBEE
      write(global_info,*)"  LIM_SUPERBEE: Advection routines use a SuperBee limiter"
#endif
#ifdef LIM_MC
      write(global_info,*)"        LIM_MC: Advection routines use a MC limiter"
#endif
#ifdef VERBOSE_L0
      write(global_info,*)"    VERBOSE_L0: Output verbosity = 0 (suppress all output to stdout)"
#endif
#ifdef VERBOSE_L1
      write(global_info,*)"    VERBOSE_L1: Output verbosity = 1"
#endif
#ifdef VERBOSE_L2
      write(global_info,*)"    VERBOSE_L2: Output verbosity = 2"
#endif
#ifdef VERBOSE_L3
      write(global_info,*)"    VERBOSE_L3: Output verbosity = 3"
#endif
#ifdef USENETCDF
      write(global_info,*)"     USENETCDF: NetCDF functionality is included"
#endif
#ifdef USEGRIB
      write(global_info,*)"       USEGRIB: Grib functionality is included"
#endif
#ifdef USEPOINTERS
      write(global_info,*)"   USEPOINTERS: Arrays are defined as pointers"
      write(global_info,*)"                This helps Ash3d work with ForestClaw"
#endif
#ifdef USEEXTDATA
      write(global_info,*)"    USEEXTDATA: Data files for airports and volcanoes are"
      write(global_info,*)"                read at run-time"
#endif

      ! WRITE OUT START TIME IN UTC
      write(global_info,*)
      write(global_log ,*)
      write(global_info,2) version,RunStartYear,RunStartMonth,RunStartDay,RunStartHr,RunStartMinute
      write(global_log ,2) version,RunStartYear,RunStartMonth,RunStartDay,RunStartHr,RunStartMinute
        ! Prepare a note to include in the netcdf output file
      write(linebuffer080,102) RunStartYear,RunstartMonth,RunStartDay,RunStartHr,RunStartMinute
      os_time_log = linebuffer080(1:17)
      write(global_info,*)
      write(global_log ,*)

2     format(4x,'Ash3d (Rev ',a5,') run ',&
             i4,'.',i2.2,'.',i2.2,i4,':',i2.2,' UTC')
102   format(i4,'.',i2.2,'.',i2.2,i4,':',i2.2)

      end subroutine Set_OS_Env

