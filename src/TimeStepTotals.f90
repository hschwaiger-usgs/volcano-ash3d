!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  TimeStepTotals(itime)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine TimeStepTotals(itime)

!  Subroutine that sums the volume erupted, deposited, and area of ash cloud at specified time steps.
      
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

      !integer            :: isize !, ix, iy, iz

!  SUM UP TOTAL ASH ON THE GROUND, IN THE AIR, AND DRIFTED OUT OF THE DOMAIN
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

      !do io=1,2;if(VB(io).le.verbosity_info)then
      !  write(outlog(io),*)sum(outflow_yz1(     1:ny,1:nz,1)*kappa(0     ,1:ny  ,1:nz  ))/MagmaDensity/KM3_2_M3,&
      !            sum(outflow_yz2(     1:ny,1:nz,1)*kappa(  nx+1,1:ny  ,1:nz  ))/MagmaDensity/KM3_2_M3,&
      !            sum(outflow_xz1(1:nx,     1:nz,1)*kappa(1:nx  ,0     ,1:nz  ))/MagmaDensity/KM3_2_M3,&
      !            sum(outflow_xz2(1:nx,     1:nz,1)*kappa(1:nx  ,  ny+1,1:nz  ))/MagmaDensity/KM3_2_M3,&
      !            sum(outflow_xy1(1:nx,1:ny,     1)*kappa(1:nx  ,1:ny  ,0     ))/MagmaDensity/KM3_2_M3,&
      !            sum(outflow_xy2(1:nx,1:ny,     1)*kappa(1:nx  ,1:ny  ,  nz+1))/MagmaDensity/KM3_2_M3
      !endif;enddo

      open(unit=19,file='progress.txt',status='replace')
      write(19,*)real(time/Simtime_in_hours,kind=4)
      close(19)

!      do io=1,2;if(VB(io).le.verbosity_info)then
!        write(outlog(io),*)"Total Mass Aloft = ",real(aloft_vol*MagmaDensity*KM3_2_M3*1.0e-9,kind=sp),&
!                " Tg"
!      endif;enddo

      OutputStep_Marker = ' ' 

2     format(6x,i6,a1,f11.3,3x,a,5f12.5,f13.1)

      return

      end subroutine TimeStepTotals
