!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  MesoInterpolater(TimeNow,Load_MesoSteps,Interval_Frac)
!
!  Called from: Main ash3d (Ash3d.F90) in two places, once before the time loop
!  and once within it
!  Arguments:
!    TimeNow        = current time, in hours since start of simulation
!    Load_MesoSteps = logical flag to indicate loading the next met time-step
!    Interval_Frac  = the fraction of the met step for the current time
!
! This purpose of this subroutine is to interpolate data from the mesoscale NWP
! files onto the current time of the simulation.  
! Initially, if the 'first_time' is set to true, then this is the call from
! Ash3d.F90 prior to the time-integration loop.  During the first call, both
! braketing time steps need to be loaded, velocities interpolated, then the
! time step calculated. The time step is calculated once all the velocities
! are determined, including the additional component from umbrella cloud
! spreading, if an umbrella source is used.  If compiled with FAST_DT set,
! then Adjust_DT is not called for each step, but only on the intevals
! of the NWP files with linear interpolation in between.
! The variables are read into memory via the subroutine Read_NextMesoStep
! After the current velocities are determined, some error-checking is
! aplied at the end of the subroutine to ensure velocities are valid.
!
! Variables set: vx_pd, vy_pd, vz_pd, vf_pd, dt
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine MesoInterpolater(TimeNow,Load_MesoSteps,Interval_Frac)

      use precis_param

      use io_units

      use global_param,    only : &
         EPS_SMALL,useFastDt,FastDt_suppress,MPS_2_KMPHR,useCalcFallVel

      use mesh,            only : &
         nxmax,nymax,nzmax

      use solution,        only : &
         vx_pd,vy_pd,vz_pd,vf_pd

      use time_data,       only : &
         Simtime_in_hours,time,dt,dt_ip,SimStartHour,dt_meso_last,dt_meso_next, &
         tw1,tw2,tw_tot

      use io_data,       only : &
         NextWriteTime

      use wind_grid,       only : &
          vx_meso_last_step_sp,vx_meso_next_step_sp,&
          vy_meso_last_step_sp,vy_meso_next_step_sp,&
          vz_meso_last_step_sp,vz_meso_next_step_sp,&
          vf_meso_last_step_sp,vf_meso_next_step_sp,&
          vx_meso_1_sp,vy_meso_1_sp,vz_meso_1_sp, &
          vx_meso_2_sp,vy_meso_2_sp,vz_meso_2_sp, &
          Meso_toggle

      use Source,          only : &
         SourceType,e_EndTime

      use Source_Umbrella,        only : &
         uvx_pd,uvy_pd,ibase,itop, &
           umbrella_winds

      use Atmosphere,      only : &
           Set_Atmosphere_Meso

      use MetReader,       only : &
         MR_iMetStep_Now,&
         MR_MetSteps_Total,&
         MR_MetStep_Hour_since_baseyear,MR_MetStep_Interval,&
           MR_Read_HGT_arrays,&
           MR_Read_3d_Met_Variable_to_CompH,&
!           MR_Rotate_UV_GR2ER_Met,&
           MR_Rotate_UV_ER2GR_Comp,&
           MR_Regrid_MetP_to_CompH,&
           MR_Read_3d_MetP_Variable
 
      implicit none

      real(kind=dp),intent(in)    :: TimeNow                ! current time, in hours since start of simulation
      real(kind=dp),intent(out)   :: Interval_Frac
      logical      ,intent(inout) :: Load_MesoSteps

      integer           :: i,j,k
      character(len=1)  :: answer
      character(len=50) :: linebuffer050 
      character(len=80) :: linebuffer080
      logical,save      :: first_time = .true.  ! There is a bit of extra work the first time
                                                ! this subroutine is called since we need to
                                                ! fill the step prior to the start time.
      integer           :: iostatus
      character(len=120):: iomessage

      real(kind=dp) :: HoursIntoInterval ! hours since the last windfile timestep
      real(kind=dp) :: TimeNow_fromRefTime

      INTERFACE
        subroutine Adjust_DT(mesostep)
          logical, intent(in), optional :: mesostep
        end subroutine Adjust_DT
        subroutine Read_NextMesoStep(Load_MesoSteps)
          logical      ,intent(inout) :: Load_MesoSteps
        end subroutine Read_NextMesoStep
      END INTERFACE

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine MesoInterpolator"
      endif;enddo

      TimeNow_fromRefTime = SimStartHour+TimeNow  ! hours since reference time (1-1-1900)
      ! MesoInterpolater is called once before the time loop in order to
      ! initilize velocities on the computational grid and to determine the
      ! start time relative to the MetSteps
      if(first_time)then
        vx_meso_last_step_sp = 0.0_sp
        vy_meso_last_step_sp = 0.0_sp
        vz_meso_last_step_sp = 0.0_sp
        Meso_toggle = 0
        ! Find the first MetStep that we need
        do i = 1,MR_MetSteps_Total-1
          if(TimeNow_fromRefTime.ge.MR_MetStep_Hour_since_baseyear(i).and.&
             TimeNow_fromRefTime.lt.MR_MetStep_Hour_since_baseyear(i+1))then
            MR_iMetStep_Now = i
            cycle
          endif
        enddo
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),'(a18,i4,2f15.3)')"MR_iMetStep_Now = ",MR_iMetStep_Now, &
                    TimeNow_fromRefTime, &
                    MR_MetStep_Hour_since_baseyear(MR_iMetStep_Now)
        endif;enddo

        call cpu_time(tw1)
        call Read_NextMesoStep(Load_MesoSteps)
        call cpu_time(tw2)
        tw_tot = tw_tot + (tw2-tw1)

        ! We only loaded one step so set load flag to True
        Load_MesoSteps = .true.
        if(useFastDt)then
          ! Here we call adjust_dt for the meso step even if it might be suppressed later
          ! since once the suppression turns off, we will need to have these available.
          call Adjust_DT(Load_MesoSteps)
          ! Since this is the first time, copy the output dt_meso_next to dt_meso_last
          dt_meso_last = dt_meso_next
        endif
      else
        ! If this is NOT the first time MesoInterpolater is called, then check
        ! if we need to load the next step
        Load_MesoSteps = .false.    ! Initialize
          ! Check if we've crossed the MetStep Boundary
        if(TimeNow_fromRefTime.gt.MR_MetStep_Hour_since_baseyear(MR_iMetStep_Now+1))then
          Load_MesoSteps = .true.
          MR_iMetStep_Now = MR_iMetStep_Now+1
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"  Need to load next step"
            write(outlog(io),'(a18,i4,2f15.3)')"MR_iMetStep_Now = ",MR_iMetStep_Now, &
                      TimeNow_fromRefTime, &
                      MR_MetStep_Hour_since_baseyear(MR_iMetStep_Now)
          endif;enddo
        endif
      endif !first_time

        ! Now load the next MetStep if needed
      if(Load_MesoSteps)then
        ! Toggle next to last
          ! Set pointer toggle
        if(Meso_toggle.eq.0)then
          Meso_toggle = 1
          vx_meso_last_step_sp = vx_meso_1_sp
          vy_meso_last_step_sp = vy_meso_1_sp
          vz_meso_last_step_sp = vz_meso_1_sp
        else
          Meso_toggle = 0
          vx_meso_last_step_sp = vx_meso_2_sp
          vy_meso_last_step_sp = vy_meso_2_sp
          vz_meso_last_step_sp = vz_meso_2_sp
        endif

        ! Copy the timestep from next to last
        dt_meso_last = dt_meso_next
        call cpu_time(tw1)
        call Read_NextMesoStep(Load_MesoSteps)
        call cpu_time(tw2)
        tw_tot = tw_tot + (tw2-tw1)

        if(useFastDt)then
          ! Now calculate the dt at the next timestep (again, even though it might
          ! be suppressed later)
          call Adjust_DT(Load_MesoSteps)
        endif

      endif   ! Load_MesoSteps

      ! Now we have the bracketing wind data for vx,vy,vz
      ! Set up to do the interpolation to the current time
      HoursIntoInterval   = TimeNow_fromRefTime - MR_MetStep_Hour_since_baseyear(MR_iMetStep_Now)
      Interval_Frac = HoursIntoInterval /  MR_MetStep_Interval(MR_iMetStep_Now)

      ! Interpolate to get wind fields at current time and in km/hr
      ! Note: v[x,y,z,f]_pd are all in km/hr and
      !       v[x,y,z,f]_meso_[last,next]step_sp are in m/s
      vx_pd(1:nxmax,1:nymax,1:nzmax) = &
            (real(vx_meso_last_step_sp(1:nxmax,1:nymax,1:nzmax),kind=ip) + &
           real((vx_meso_next_step_sp(1:nxmax,1:nymax,1:nzmax) - &
                 vx_meso_last_step_sp(1:nxmax,1:nymax,1:nzmax)),kind=ip) * &
                                     real(Interval_Frac,kind=ip))*MPS_2_KMPHR

      vy_pd(1:nxmax,1:nymax,1:nzmax) = &
            (real(vy_meso_last_step_sp(1:nxmax,1:nymax,1:nzmax),kind=ip) + &
           real((vy_meso_next_step_sp(1:nxmax,1:nymax,1:nzmax) - &
                 vy_meso_last_step_sp(1:nxmax,1:nymax,1:nzmax)),kind=ip) * &
                                     real(Interval_Frac,kind=ip))*MPS_2_KMPHR

      vz_pd(1:nxmax,1:nymax,1:nzmax) = &
            (real(vz_meso_last_step_sp(1:nxmax,1:nymax,1:nzmax),kind=ip) + &
           real((vz_meso_next_step_sp(1:nxmax,1:nymax,1:nzmax) - &
                 vz_meso_last_step_sp(1:nxmax,1:nymax,1:nzmax)),kind=ip) * &
                                     real(Interval_Frac,kind=ip))*MPS_2_KMPHR

      ! Now interpolate onto current time
      if(useCalcFallVel)then
        vf_pd(1:nxmax,1:nymax,1:nzmax,:) = (real(vf_meso_last_step_sp(:,:,:,:),kind=ip) + &
                                            real((vf_meso_next_step_sp(:,:,:,:) - &
                                                vf_meso_last_step_sp(:,:,:,:)),kind=ip) * &
                                            real(Interval_Frac,kind=ip))*MPS_2_KMPHR
      endif

      ! if we're calculating an umbrella cloud, add the winds in the cloud
      if (((SourceType.eq.'umbrella')     .or.   &
           (SourceType.eq.'umbrella_air')))then
        if(TimeNow.lt.e_EndTime(1))then
          FastDt_suppress = .true.
          call umbrella_winds(first_time)
          vx_pd(1:nxmax,1:nymax,ibase:itop) = vx_pd(1:nxmax,1:nymax,ibase:itop) + &
                                           uvx_pd(1:nxmax,1:nymax,ibase:itop)
          vy_pd(1:nxmax,1:nymax,ibase:itop) = vy_pd(1:nxmax,1:nymax,ibase:itop) + &
                                           uvy_pd(1:nxmax,1:nymax,ibase:itop)
        else
          ! Even though the umbrella winds change the windfield, we can turn Fast_DT
          ! back on since the umbrella cloud has stopped spreading.
          FastDt_suppress = .false.
        endif
      endif

      if(useFastDt.and..not.FastDt_suppress)then
        dt = dt_meso_last + (dt_meso_next-dt_meso_last)*Interval_Frac
        if (((NextWriteTime-time).gt.EPS_SMALL).and.(NextWriteTime-time.lt.dt))then
          ! Don't let time advance past the next output time step
          dt = NextWritetime-time
        endif
        if(time+dt.gt.Simtime_in_hours)then
          ! Don't let time advance past the requested stop time
          dt = Simtime_in_hours-time
        endif
        dt_ip = real(dt,kind=ip)
      else
        call Adjust_DT(.false.)
      endif

      ! If vxmax or vymax > 2000 km/hr, send a warning and see which values of vx and vy
      ! are so high.
      if (abs(maxval(vx_pd(1:nxmax,1:nymax,1:nzmax))).gt.5.0e3_ip) then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),1041) nxmax, nymax, nzmax, SimStartHour, &
                      MR_MetStep_Hour_since_baseyear(MR_iMetStep_Now), &
                      MR_MetStep_Interval(MR_iMetStep_Now)
        endif;enddo
1041    format(4x,'Warning: max(vx) in Mesointerpolater is > 5.0e3.',&
                  '  Writing out locations of all ',/, &
                  'values greater than 5e3',/,4x,' nxmax=',i4,&
                  ' nymax=',i4,' nzmax=',i3,/, &
                  'SimStartHour=',f12.2,' MetStep_Hour_since_baseyear(MR_iMetStep_Now)',&
                  f12.2,' ForecastInterval=',f5.2)
        do k=1,nzmax
          do i=1,nxmax
            do j=1,nymax
              if (abs(vx_pd(i,j,k)).gt.5.0e3_ip) then
                do io=1,2;if(VB(io).le.verbosity_info)then
                  write(outlog(io),1042) i,j,k,vx_pd(i,j,k), & 
                              vx_meso_last_step_sp(i,j,k), vx_meso_next_step_sp(i,j,k),&
                              HoursIntoInterval
                endif;enddo
1042            format(4x,'i=',i4,' iy=',i4,' iz=',i4,' vx=',e12.4, &
                          ' vx_regrid(1)=',e12.4,' vx_regrid(2)=',e12.4,e12.4)
              endif
            enddo
          enddo
        enddo
        if(VB(1).ge.verbosity_silent)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*) 'Continue (y/n)?'
          endif;enddo
          read(input_unit,'(a1)',iostat=iostatus,iomsg=iomessage) answer
          linebuffer080 = answer
          linebuffer050 = "Reading from stdin, answer"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
          if (adjustl(trim(answer)).eq.'n') stop 1
        else
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: velocities seem to be out of expected range."
          endif;enddo
          stop 1
        endif
      endif

      if (abs(maxval(vy_pd(1:nxmax,1:nymax,1:nzmax))).gt.5.0e3_ip) then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),1043) nxmax, nymax, nzmax, SimStartHour, &
                      MR_MetStep_Hour_since_baseyear(MR_iMetStep_Now), &
                      MR_MetStep_Interval(MR_iMetStep_Now)
        endif;enddo
1043    format(4x,'Warning: max(vy) in Mesointerpolater is > 5e3.',&
                  '  Writing out locations of all ',/, &
                  'values greater than 5e3',/,4x,' nxmax=',i4,&
                  ' nymax=',i4,' nzmax=',i3,/, &
                  'SimStartHour=',f12.2,' MetStep_Hour_since_baseyear(MR_iMetStep_Now)',&
                  f12.2,' ForecastInterval=',f5.2)

        do k=1,nzmax
          do i=1,nxmax
            do j=1,nymax
              if (vy_pd(i,j,k).gt.5.0e3_ip) then
                do io=1,2;if(VB(io).le.verbosity_info)then
                  write(outlog(io),1044) i,j,k,vy_pd(i,j,k),  & 
                              vy_meso_last_step_sp(i,j,k), vy_meso_next_step_sp(i,j,k)
                endif;enddo
1044            format(4x,'i=',i4,' j=',i4,' k=',i4,' vy=',e12.4,/, &
                          ' vy_regrid(1)=',e12.4,' vy_regrid(2)=',e12.4)
              endif
            enddo
          enddo
        enddo
        if(VB(1).ge.verbosity_silent)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*) 'Continue (y/n)?'
          endif;enddo
          read(input_unit,'(a1)',iostat=iostatus,iomsg=iomessage) answer
          linebuffer080 = answer
          linebuffer050 = "Reading from stdin, answer"
          if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
          if (adjustl(trim(answer)).eq.'n') stop 1
        else
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: velocities seem to be out of expected range."
          endif;enddo
          stop 1
        endif
      endif

      first_time = .false.

      return

      end subroutine MesoInterpolater

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Read_NextMesoStep(Load_MesoSteps)
!
!  Called from: MesoInterpolater
!  Arguments:
!    Load_MesoSteps = logical flag to indicate loading the next met time-step
!
! This subroutine manages the calls to MetReader to populate the variables at
! the braketing timesteps of the NWP files, reading the next step of the
! simulation time crosses over to the next time step of the NWP files.  The
! primary variables read from the NWP files are vx, vy, and vz, but if temperatures
! are needed (as is the case for most fall models), then those are read via calls
! to Set_Atmosphere_Meso.  Fall velocities are calculated via calls to Set_Vf_Meso.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Read_NextMesoStep(Load_MesoSteps)

      use precis_param

      use io_units

      use global_param,    only : &
         MPS_2_KMPHR,GRAV, &
         useTemperature,useCalcFallVel,useVz_rhoG

      use solution,        only : &
         vf_pd

      use wind_grid,       only : &
          vx_meso_next_step_sp,&
          vy_meso_next_step_sp,&
          vz_meso_next_step_sp,&
          vf_meso_last_step_sp,vf_meso_next_step_sp,&
          vx_meso_1_sp,vy_meso_1_sp,vz_meso_1_sp, &
          vx_meso_2_sp,vy_meso_2_sp,vz_meso_2_sp, &
          Meso_toggle

      use Tephra,          only : &
         n_gs_max,Tephra_v_s,   &
           Set_Vf_Meso

      use Atmosphere,      only : &
         AirDens_meso_next_step_MetP_sp, &
           Set_Atmosphere_Meso

      use MetReader,       only : &
         MR_dum3d_compH,MR_dum3d_compH_2,MR_iMetStep_Now,&
         Met_var_IsAvailable,isGridRelative,Map_Case,&
         MR_dum3d_metP,Met_var_GRIB_names,&
           MR_Read_HGT_arrays,&
           MR_Read_3d_Met_Variable_to_CompH,&
           MR_Rotate_UV_GR2ER_Met,&
           MR_Rotate_UV_ER2GR_Comp,&
           MR_Regrid_MetP_to_CompH,&
           MR_Read_3d_MetP_Variable

      implicit none

      logical      ,intent(inout) :: Load_MesoSteps

      logical,save :: first_time = .true.

      integer           :: isize
      integer           :: ivar
      integer           :: istep

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Read_NextMesoStep"
      endif;enddo

        ! Before reading state variables, we need to load the height grid
        ! which will be used in the QC checking
      if(first_time)then
        call MR_Read_HGT_arrays(MR_iMetStep_Now,first_time)
        istep = MR_iMetStep_Now
      else
        ! This subroutine sets both last and next geoH arrays so call with
        ! MR_iMetStep_Now
        call MR_Read_HGT_arrays(MR_iMetStep_Now)
        istep = MR_iMetStep_Now+1
      endif
      if(Map_Case.eq.1.or.Map_Case.eq.2)then
        ! Either both the comp and met grids are LL (Map_Case = 1)
        ! or they are both the same projection (Map_Case = 2) so
        ! we can read the velocity components individually and interpolate onto
        ! the computational grid

         ! Fill array from the step prior/equal to current time
        ivar = 2 ! U winds
        call MR_Read_3d_Met_Variable_to_CompH(ivar,istep,.true.)
        if(Meso_toggle.eq.0)then
          vx_meso_1_sp = MR_dum3d_compH
          vx_meso_next_step_sp = vx_meso_1_sp
        else
          vx_meso_2_sp = MR_dum3d_compH
          vx_meso_next_step_sp = vx_meso_2_sp
        endif

        ivar = 3 ! V winds 
        call MR_Read_3d_Met_Variable_to_CompH(ivar,istep,.true.)
        if(Meso_toggle.eq.0)then
          vy_meso_1_sp = MR_dum3d_compH
          vy_meso_next_step_sp = vy_meso_1_sp
        else
          vy_meso_2_sp = MR_dum3d_compH
          vy_meso_next_step_sp = vy_meso_2_sp
        endif
      else
        ! Grids are different, we will need to rotate vectors
        ! In all these cases, we need:
        !    MR_dum3d_compH   holding U
        !    MR_dum3d_compH_2 holding V
        if(Map_Case.eq.3)then
            ! Met grid is natively LL and Comp grid is projected
          call MR_Rotate_UV_ER2GR_Comp(istep)
        elseif(Map_Case.eq.4)then
            ! Met grid is projected and comp grid is LL
          if(isGridRelative)then
            call MR_Rotate_UV_GR2ER_Met(istep,.true.,.true.) ! optional argument returns data on compH
          else
            ! if the projected data is already Earth-relative (NARR), then just read it
            ivar = 3 ! Vy
            call MR_Read_3d_Met_Variable_to_CompH(ivar,istep)
            MR_dum3d_compH_2 = MR_dum3d_compH
            ivar = 2 ! Vx
            call MR_Read_3d_Met_Variable_to_CompH(ivar,istep)
          endif
        elseif(Map_Case.eq.5)then
          ! Both comp and met grids are projected, but with different projections
          ! First, rotate winds from met grid to earth-relative on the met nodes
          call MR_Rotate_UV_GR2ER_Met(istep)
          ! Now rotate and interpolate those earth-relative values to the projected comp grid.
          call MR_Rotate_UV_ER2GR_Comp(istep)
        endif

        if(Meso_toggle.eq.0)then
          vx_meso_1_sp = MR_dum3d_compH
          vx_meso_next_step_sp = vx_meso_1_sp
          vy_meso_1_sp = MR_dum3d_compH_2
          vy_meso_next_step_sp = vy_meso_1_sp
        else
          vx_meso_2_sp = MR_dum3d_compH
          vx_meso_next_step_sp = vx_meso_2_sp
          vy_meso_2_sp = MR_dum3d_compH_2
          vy_meso_next_step_sp = vy_meso_2_sp
        endif
      endif

      vf_meso_last_step_sp = vf_meso_next_step_sp
      ! Only bother getting Vz if it is available in the wind files
      if(useTemperature)then
        call Set_Atmosphere_Meso(Load_MesoSteps,1.0_ip,first_time)
      endif
      if(Met_var_IsAvailable(4))then
        if(useTemperature.and.useVz_rhoG)then
          ! Now that we have temperature and density, we can get a better Vz
          ivar = 7 ! Pressure Vertical Velocity
          if(Met_var_IsAvailable(ivar))then
            call MR_Read_3d_MetP_Variable(ivar,istep)
            MR_dum3d_MetP = MR_dum3d_MetP/      &
             real((-AirDens_meso_next_step_MetP_sp*GRAV),kind=sp)
            call MR_Regrid_MetP_to_CompH(istep)
          else
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"Tried to read variable, but it's not available: ",ivar,&
                        Met_var_GRIB_names(ivar)
            endif;enddo
            MR_dum3d_compH = 0.0_sp
            stop 1
          endif
          if(Meso_toggle.eq.0)then
            vz_meso_1_sp = MR_dum3d_compH
            vz_meso_next_step_sp = vz_meso_1_sp
          else
            vz_meso_2_sp = MR_dum3d_compH
            vz_meso_next_step_sp = vz_meso_2_sp
          endif
        else
          ivar = 4 ! W winds
          if(Met_var_IsAvailable(ivar))then
            call MR_Read_3d_Met_Variable_to_CompH(ivar,istep)
          else
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"Tried to read variable, but it's not available: ",ivar,&
                      Met_var_GRIB_names(ivar)
            endif;enddo
            MR_dum3d_compH = 0.0_sp
            stop 1
          endif
          if(Meso_toggle.eq.0)then
            vz_meso_1_sp = MR_dum3d_compH
            vz_meso_next_step_sp = vz_meso_1_sp
          else
            vz_meso_2_sp = MR_dum3d_compH
            vz_meso_next_step_sp = vz_meso_2_sp
          endif
        endif
      else
        vz_meso_next_step_sp = 0.0_sp
      endif

      if(useCalcFallVel)then
        ! Populate the fall velocities for the next meso step
        call Set_Vf_Meso(.true.)
      else
        if(first_time)then
          ! This is the case where the fall velocity is assigned in the input file
          do isize=1,n_gs_max
            vf_meso_next_step_sp(:,:,:,isize) = real(Tephra_v_s(isize),kind=sp)
            vf_pd(:,:,:,isize)                =      Tephra_v_s(isize)*MPS_2_KMPHR
          enddo
        endif
      endif

      first_time = .false.

      return

      end subroutine Read_NextMesoStep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
