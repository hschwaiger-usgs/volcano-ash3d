!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  TimeStepTotals(itime)
!
!  Called from: Ash3d.F90, once within the time integration loop and once at the
!               end.
!  Arguments:
!    itime = integer time step
!
!  This subroutine call subroutine to calculate the volume of ash aloft, in the
!  deposit and that exited the domain through the sides and possibly the top
!  at the specified time step itime.  This information, along with the time and
!  cloud area are writen to stdout and the logfile.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine TimeStepTotals(itime)

      use precis_param

      use io_units

      use solution,      only : &
         aloft_vol,dep_vol,outflow_vol,tot_vol,SourceCumulativeVol

      use Output_Vars,   only : &
         CloudArea,&
         Calc_AshVol_Aloft,&
         Calc_AshVol_Deposit,&
         Calc_AshVol_Outflow

      use time_data,     only : &
         time,SimStartHour,Simtime_in_hours,BaseYear,useLeap

      use io_data,       only : &
         OutputStep_Marker

      implicit none

      integer, intent(in) :: itime

      character(len=13) :: DateTime  !function that calculates date string given hours since base year (1900)

      INTERFACE
        character (len=13) function HS_yyyymmddhhmm_since(HoursSince,byear,useLeaps)
          real(kind=8)               ::  HoursSince
          integer                    ::  byear
          logical                    ::  useLeaps
        end function HS_yyyymmddhhmm_since
      END INTERFACE

      ! Sum up the total ash aloft, on the ground and that left the domain
      tot_vol     = 0.0_ip
      dep_vol     = 0.0_ip 
      aloft_vol   = 0.0_ip
      outflow_vol = 0.0_ip

      call Calc_AshVol_Aloft(aloft_vol)

      call Calc_AshVol_Deposit(dep_vol)

      call Calc_AshVol_Outflow(outflow_vol)

      tot_vol = aloft_vol + dep_vol + outflow_vol

      DateTime = HS_yyyymmddhhmm_since(time+SimStartHour,BaseYear,useLeap)

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),2) itime,OutputStep_Marker,time,DateTime,&
                         SourceCumulativeVol,dep_vol,aloft_vol,&
                         outflow_vol,tot_vol,CloudArea
      endif;enddo

      open(unit=19,file='progress.txt',status='replace')
      write(19,*)real(time/Simtime_in_hours,kind=4)
      close(19)

      OutputStep_Marker = ' ' 

2     format(6x,i6,a1,f11.3,3x,a,5f12.5,f13.1)

      return

      end subroutine TimeStepTotals
