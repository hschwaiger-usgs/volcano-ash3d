      subroutine MesoInterpolater(TimeNow,Load_MesoSteps,Interval_Frac,first_time)
      ! Fclaw subroutine MesoInterpolater(TimeNow,Load_MesoSteps,Interval_Frac,first_time,Meso_toggle)

!    Subroutine that interpolates to obtain the current wind field
!
      use precis_param

      use io_units

      use global_param,    only : &
         useTemperature,useCalcFallVel,MPS_2_KMPHR

      use solution,        only : &
         vx_pd,vy_pd,vz_pd

      use wind_grid,       only : &
          vx_meso_last_step_sp,vx_meso_next_step_sp,&
          vy_meso_last_step_sp,vy_meso_next_step_sp,&
          vz_meso_last_step_sp,vz_meso_next_step_sp,&
          vx_meso_1_sp,vy_meso_1_sp,vz_meso_1_sp, &
          vx_meso_2_sp,vy_meso_2_sp,vz_meso_2_sp, &
          Meso_toggle

      use time_data,       only : &
         SimStartHour

      use mesh,            only : &
         nxmax,nymax,nzmax

      use Source,          only : &
         uvx_pd,uvy_pd,ibase,itop,SourceType,e_EndTime

      use Tephra,          only : &
         Set_Vf_Meso

      use Atmosphere,      only : &
           Set_Atmosphere_Meso

      use MetReader,       only : &
         MR_dum3d_compH,MR_dum3d_compH_2,MR_iMetStep_Now,&
         MR_MetSteps_Total,Met_var_IsAvailable,isGridRelative,Map_Case,&
         MR_MetStep_Hour_since_baseyear,MR_MetStep_Interval,&
         MR_dum3d_compH,&
           MR_Read_HGT_arrays,&
           MR_Read_3d_Met_Variable_to_CompGrid,&
           MR_Rotate_UV_GR2ER_Met,&
           MR_Rotate_UV_ER2GR_Comp
 
      implicit none

      real(kind=dp),intent(in)  :: TimeNow                !current time, in hours since start of simulation
      real(kind=dp),intent(out) :: Interval_Frac
      logical      ,intent(out) :: Load_MesoSteps
      logical      ,intent(in)  :: first_time

      integer           :: i
      integer           :: ivar
      integer           :: ix,iy,iz
      character(len=1)  :: answer

      real(kind=dp):: HoursIntoInterval ! hours since the last windfile timestep
      real(kind=ip) :: TimeNow_fromRefTime
      !logical,save :: first_time = .true.
      ! Fclaw
      !logical, intent(in) :: first_time
      !integer, intent(inout) :: Meso_toggle

      INTERFACE
        subroutine umbrella_winds
        end subroutine umbrella_winds
      END INTERFACE

      TimeNow_fromRefTime = SimStartHour+TimeNow  ! hours since reference time (1-1-1900)

      ! MesoInterpolater is called once before the time loop in order to
      ! initilize velocities on the computational grid and to determine the
      ! start time relative to the MetSteps
      if(first_time)then
        Meso_toggle = 0
        ! Find the first MetStep that we need
        do i = 1,MR_MetSteps_Total-1
          if(TimeNow_fromRefTime.ge.MR_MetStep_Hour_since_baseyear(i).and.&
             TimeNow_fromRefTime.lt.MR_MetStep_Hour_since_baseyear(i+1))then
            MR_iMetStep_Now = i
            cycle
          endif
        enddo
        write(global_info,*)"MR_iMetStep_Now = ",MR_iMetStep_Now, &
                    TimeNow_fromRefTime, &
                    MR_MetStep_Hour_since_baseyear(MR_iMetStep_Now)

        ! Before reading state variables, we need to load the height grid
        ! which will be used in the QC checking
        call MR_Read_HGT_arrays(MR_iMetStep_Now,first_time)

        if(Map_Case.eq.1.or.Map_Case.eq.2)then
          ! Either both the comp and met1 grids are LL (Map_Case = 1)
          ! or they are both the same projection (Map_Case = 2) so
          ! we can read the velocity components individually and interpolate onto
          ! the computational grid

           ! Fill array from the step prior/equal to current time
          ivar = 2 ! U winds
          call MR_Read_3d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now,.true.)
            vx_meso_1_sp = MR_dum3d_compH
            vx_meso_next_step_sp = vx_meso_1_sp

          ivar = 3 ! V winds 
          call MR_Read_3d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now,.true.)
            vy_meso_1_sp = MR_dum3d_compH
            vy_meso_next_step_sp = vy_meso_1_sp
        else
          ! Grids are different, we will need to rotate vectors
          ! In all these cases, we need:
          !    MR_dum3d_compH   holding U
          !    MR_dum3d_compH_2 holding V
          if(Map_Case.eq.3)then
              ! Met grid is natively LL and Comp grid is projected
            call MR_Rotate_UV_ER2GR_Comp(MR_iMetStep_Now)
          elseif(Map_Case.eq.4)then
              ! Met grid is projected and comp grid is LL
            if(isGridRelative)then
              call MR_Rotate_UV_GR2ER_Met(MR_iMetStep_Now,.true.) ! optional argument returns data on compH
            else
              ! if the projected data is already Earth-relative (NARR), then just read it
              ivar = 3 ! Vy
              call MR_Read_3d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now)
              MR_dum3d_compH_2 = MR_dum3d_compH
              ivar = 2 ! Vx
              call MR_Read_3d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now)
            endif
          elseif(Map_Case.eq.5)then
            ! Both comp and met grids are projected, but with different projections
            ! First, rotate winds from met grid to earth-relative on the met nodes
            call MR_Rotate_UV_GR2ER_Met(MR_iMetStep_Now)
            ! Now rotate and interpolate those earth-relative values to the projected comp grid.
            call MR_Rotate_UV_ER2GR_Comp(MR_iMetStep_Now)
          endif

          vx_meso_1_sp = MR_dum3d_compH
          vx_meso_next_step_sp = vx_meso_1_sp
          vy_meso_1_sp = MR_dum3d_compH_2
          vy_meso_next_step_sp = vy_meso_1_sp
        endif

        ivar = 4 ! W winds
        if(Met_var_IsAvailable(ivar))then
          call MR_Read_3d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now)
        else
          MR_dum3d_compH = 0.0_sp
        endif
        vz_meso_1_sp = MR_dum3d_compH
        vz_meso_next_step_sp = vz_meso_1_sp

        ! We only loaded one step so set load flag to True
        Load_MesoSteps = .true.
        !first_time = .false.
      else
        ! If this is NOT the first time MesoInterpolater is called, then check
        ! if we need to load the next step
        Load_MesoSteps = .false.    ! Initialize
          ! Check if we've crossed the MetStep Boundary
        if(TimeNow_fromRefTime.gt.MR_MetStep_Hour_since_baseyear(MR_iMetStep_Now+1))then
          Load_MesoSteps = .true.
          MR_iMetStep_Now = MR_iMetStep_Now+1
          write(global_info,*)"  Need to load next step"
          write(global_info,*)"MR_iMetStep_Now = ",MR_iMetStep_Now, &
                      TimeNow_fromRefTime, &
                      MR_MetStep_Hour_since_baseyear(MR_iMetStep_Now)
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

        ! This subroutine sets both last and next geoH arrays so call with
        ! MR_iMetStep_Now
        call MR_Read_HGT_arrays(MR_iMetStep_Now)

        if(Map_Case.eq.1.or.Map_Case.eq.2)then
          ! Either both the comp and met1 grids are LL (Map_Case = 1)
          ! or they are both the same projection (Map_Case = 2) so
          ! we can read the velocity components individually and interpolate onto
          ! the computational grid

           ! Fill array from the step prior/equal to current time
          ivar = 2 ! U winds
          call MR_Read_3d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now+1,.true.)
          if(Meso_toggle.eq.0)then
            vx_meso_1_sp = MR_dum3d_compH
            vx_meso_next_step_sp = vx_meso_1_sp
          else
            vx_meso_2_sp = MR_dum3d_compH
            vx_meso_next_step_sp = vx_meso_2_sp
          endif

          ivar = 3 ! V winds 
          call MR_Read_3d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now+1,.true.)
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
            call MR_Rotate_UV_ER2GR_Comp(MR_iMetStep_Now+1)
          elseif(Map_Case.eq.4)then
              ! Met grid is projected and comp grid is LL
            if(isGridRelative)then
              call MR_Rotate_UV_GR2ER_Met(MR_iMetStep_Now+1,.true.) ! optional argument returns data on compH
            else
              ! if the projected data is already Earth-relative (NARR), then just read it
              ivar = 3 ! Vy
              call MR_Read_3d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now+1)
              MR_dum3d_compH_2 = MR_dum3d_compH
              ivar = 2 ! Vx
              call MR_Read_3d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now+1)
            endif
          elseif(Map_Case.eq.5)then
            ! Both comp and met grids are projected, but with different projections
            ! First, rotate winds from met grid to earth-relative on the met nodes
            call MR_Rotate_UV_GR2ER_Met(MR_iMetStep_Now+1)
            ! Now rotate and interpolate those earth-relative values to the projected comp grid.
            call MR_Rotate_UV_ER2GR_Comp(MR_iMetStep_Now+1)
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

        ivar = 4 ! W winds
        if(Met_var_IsAvailable(ivar))then
          call MR_Read_3d_Met_Variable_to_CompGrid(ivar,MR_iMetStep_Now+1)
        else
          MR_dum3d_compH = 0.0_sp
        endif
        if(Meso_toggle.eq.0)then
          vz_meso_1_sp = MR_dum3d_compH
          vz_meso_next_step_sp = vz_meso_1_sp
        else
          vz_meso_2_sp = MR_dum3d_compH
          vz_meso_next_step_sp = vz_meso_2_sp
        endif
      endif   ! Load_MesoSteps

      ! Now we have the bracketing wind data for vx,vy,vz
      ! Set up to do the interpolation to the current time
      HoursIntoInterval   = TimeNow_fromRefTime - MR_MetStep_Hour_since_baseyear(MR_iMetStep_Now)
      Interval_Frac = HoursIntoInterval /  MR_MetStep_Interval(MR_iMetStep_Now)
      if(useTemperature)then
        call Set_Atmosphere_Meso(Load_MesoSteps,Interval_Frac,first_time)
      endif
      if(useCalcFallVel)then
        call Set_Vf_Meso(Load_MesoSteps,Interval_Frac)
      endif

      ! INTERPOLATE TO GET CURRENT WIND FIELDS VX, VY 
      vx_pd = 0.0_ip
      vy_pd = 0.0_ip
      vz_pd = 0.0_ip
      vx_pd(1:nxmax,1:nymax,1:nzmax) = &
            real(vx_meso_last_step_sp(1:nxmax,1:nymax,1:nzmax),kind=ip) + &
           real((vx_meso_next_step_sp(1:nxmax,1:nymax,1:nzmax) - &
                 vx_meso_last_step_sp(1:nxmax,1:nymax,1:nzmax)),kind=ip) * &
                                     Interval_Frac

      vy_pd(1:nxmax,1:nymax,1:nzmax) = &
            real(vy_meso_last_step_sp(1:nxmax,1:nymax,1:nzmax),kind=ip) + &
           real((vy_meso_next_step_sp(1:nxmax,1:nymax,1:nzmax) - &
                 vy_meso_last_step_sp(1:nxmax,1:nymax,1:nzmax)),kind=ip) * &
                                     Interval_Frac


      vz_pd(1:nxmax,1:nymax,1:nzmax) = &
            real(vz_meso_last_step_sp(1:nxmax,1:nymax,1:nzmax),kind=ip) + &
           real((vz_meso_next_step_sp(1:nxmax,1:nymax,1:nzmax) - &
                 vz_meso_last_step_sp(1:nxmax,1:nymax,1:nzmax)),kind=ip) * &
                                     Interval_Frac

      !if we're calculating an umbrella cloud, add the winds in the cloud
      if (((SourceType.eq.'umbrella')     .or.   &
           (SourceType.eq.'umbrella_air')).and.  &   !(TimeNow.gt.0.0_ip).and. &
          (TimeNow.lt.e_EndTime(1))) then
          call umbrella_winds
          vx_pd(1:nxmax,1:nymax,ibase:itop) = vx_pd(1:nxmax,1:nymax,ibase:itop) + &
                                           uvx_pd(1:nxmax,1:nymax,ibase:itop)
          vy_pd(1:nxmax,1:nymax,ibase:itop) = vy_pd(1:nxmax,1:nymax,ibase:itop) + &
                                           uvy_pd(1:nxmax,1:nymax,ibase:itop)
      endif

!     CONVERT FROM M/S TO KM/HR.  
      vx_pd = vx_pd*MPS_2_KMPHR
      vy_pd = vy_pd*MPS_2_KMPHR
      vz_pd = vz_pd*MPS_2_KMPHR

      !If vxmax or vymax > 2000 km/hr, send a warning and see which values of vx and vy
      !are so high.
      if (abs(maxval(vx_pd(1:nxmax,1:nymax,1:nzmax))).gt.5.0e3_ip) then
        write(global_info,1041) nxmax, nymax, nzmax, SimStartHour, &
                      MR_MetStep_Hour_since_baseyear(MR_iMetStep_Now), &
                      MR_MetStep_Interval(MR_iMetStep_Now)
1041    format(4x,'Warning: max(vx) in Mesointerpolater is > 5.0e3.',&
                  '  Writing out locations of all ',/, &
                  'values greater than 5e3',/,4x,' nxmax=',i4,&
                  ' nymax=',i4,' nzmax=',i3,/, &
                  'SimStartHour=',f12.2,' MetStep_Hour_since_baseyear(MR_iMetStep_Now)',&
                  f12.2,' ForecastInterval=',f5.2)
        do iz=1,nzmax
          do ix=1,nxmax
            do iy=1,nymax
              if (abs(vx_pd(ix,iy,iz)).gt.5.0e3_ip) then
                write(global_info,1042) ix, iy, iz, vx_pd(ix,iy,iz), & 
                              vx_meso_last_step_sp(ix,iy,iz), vx_meso_next_step_sp(ix,iy,iz),&
                              HoursIntoInterval
1042            format(4x,'ix=',i4,' iy=',i4,' iz=',i4,' vx=',e12.4, &
                          ' vx_regrid(1)=',e12.4,' vx_regrid(2)=',e12.4,e12.4)
                !stop 1
              endif
            enddo
          enddo
        enddo
        write(global_info,*) 'Continue (y/n)?'
        read(5,'(a1)') answer
        if (answer.eq.'n') stop 1
      endif

      if (abs(maxval(vy_pd(1:nxmax,1:nymax,1:nzmax))).gt.5.0e3_ip) then
        write(global_info,1043) nxmax, nymax, nzmax, SimStartHour, &
                      MR_MetStep_Hour_since_baseyear(MR_iMetStep_Now), &
                      MR_MetStep_Interval(MR_iMetStep_Now)
1043    format(4x,'Warning: max(vy) in Mesointerpolater is > 5e3.',&
                  '  Writing out locations of all ',/, &
                  'values greater than 5e3',/,4x,' nxmax=',i4,&
                  ' nymax=',i4,' nzmax=',i3,/, &
                  'SimStartHour=',f12.2,' MetStep_Hour_since_baseyear(MR_iMetStep_Now)',&
                  f12.2,' ForecastInterval=',f5.2)

        do iz=1,nzmax
          do ix=1,nxmax
            do iy=1,nymax
              if (vy_pd(ix,iy,iz).gt.5.0e3_ip) then
                write(global_info,1044) ix, iy, iz, vy_pd(ix,iy,iz),  & 
                              vy_meso_last_step_sp(ix,iy,iz), vy_meso_next_step_sp(ix,iy,iz)
1044            format(4x,'ix=',i4,' iy=',i4,' iz=',i4,' vy=',e12.4,/, &
                          ' vy_regrid(1)=',e12.4,' vy_regrid(2)=',e12.4)
              endif
            enddo
          enddo
        enddo
        write(global_info,*) 'Continue (y/n)?'
        read(5,'(a1)') answer
        if (answer.eq.'n') stop 1
      endif

      return

      end subroutine MesoInterpolater
