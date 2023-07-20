      program Ash3d

      use precis_param

      use io_units

      use global_param,  only : &
         useCalcFallVel,useDiffusion,useHorzAdvect,useVertAdvect,&
         HR_2_S,useTemperature,DT_MIN,EPS_TINY,EPS_SMALL,&
         nmods,OPTMOD_names,StopConditions,CheckConditions      

      use mesh,          only : &
         ivent,jvent,nxmax,nymax,nzmax,nsmax,ts0,ts1

      use solution,      only : &
         concen_pd,DepositGranularity,StopValue,aloft_percent_remaining, &
         SourceCumulativeVol,dep_vol,aloft_vol,outflow_vol,tot_vol

      use Output_Vars,   only : &
         DepositAreaCovered,DepositThickness,LoadVal,CloudLoadArea,&
         Calculated_Cloud_Load,Calculated_AshThickness,Calc_vprofile, &
           Allocate_Output_UserVars, &
           Allocate_NTime,   &
           Allocate_Profile, &
           Gen_Output_Vars,&
           FirstAsh

      use io_data,       only : &
         Called_Gen_Output_Vars,isFinal_TS,LoadConcen,log_step,&
         Output_at_logsteps,Output_at_WriteTimes,Output_every_TS,&
         NextWriteTime,iTimeNext,nvprofiles,nWriteTimes,&
         Write_PT_Data,Write_PR_Data

      use time_data,     only : &
         time,dt,Simtime_in_hours,t0,t1,ntmax

      use Source,        only : &
         SourceNodeFlux,e_EndTime_final,e_Volume,MassFluxRate_now,&
         SourceType, &
           Allocate_Source_time,&
           MassFluxCalculator,&
           TephraSourceNodes,&
           SourceVolInc

      use Source_Umbrella, only : &
         ibase,itop,SourceNodeFlux_Umbrella, &
           Allocate_Source_Umbrella,&
           TephraSourceNodes_Umbrella,&
           SourceVolInc_Umbrella,&
           AvgCon_Umbrella

      use Tephra,        only : &
         n_gs_max,n_gs_aloft,&
           Allocate_Tephra,&
           Allocate_Tephra_Met,&
           Prune_GS
   
      use Atmosphere,    only : &
           Allocate_Atmosphere_Met

      use AdvectionHorz, only : &
           AdvectHorz

      use AdvectionVert_DCU, only : &
           advect_z

      use Diffusion,     only : &
           DiffuseHorz,&
           DiffuseVert

      use Airports,      only : &
         Airport_thickness_TS,Airport_thickness,nairports,&
           ReadAirports

#ifdef USENETCDF
      use Ash3d_Netcdf,  only : &
           NC_RestartFile_LoadConcen
#endif

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert 'use' statements here
!
!------------------------------------------------------------------------------

      implicit none

!      integer               :: iostatus
      integer               :: itime
      integer               :: i,j,k
      integer               :: iz,isize
      real(kind=ip)         :: Interval_Frac
      logical               :: Load_MesoSteps
      logical               :: StopTimeLoop   = .false.
      logical               :: first_time     = .true.
!      character(len=130)    :: tmp_str
      real(kind=ip)         :: MassConsErr

      INTERFACE
        subroutine Set_OS_Env
        end subroutine Set_OS_Env
        subroutine Read_Control_File
        end subroutine Read_Control_File
        subroutine input_data_ResetParams
        end subroutine input_data_ResetParams
        subroutine alloc_arrays
        end subroutine alloc_arrays
        subroutine calc_mesh_params
        end subroutine calc_mesh_params
        subroutine MesoInterpolater(TimeNow,Load_MesoSteps,Interval_Frac,first_time)
          integer,parameter  :: dp         = 8 ! Double precision
          real(kind=dp),intent(in)    :: TimeNow
          real(kind=dp),intent(out)   :: Interval_Frac
          logical      ,intent(inout) :: Load_MesoSteps
          logical      ,intent(in)    :: first_time
        end subroutine MesoInterpolater
        subroutine output_results
        end subroutine output_results
        subroutine Set_BC(bc_code)
          integer,intent(in) :: bc_code ! 1 for advection, 2 for diffusion
        end subroutine Set_BC
        subroutine vprofilewriter(itime)
          integer, intent(in) :: itime
        end subroutine vprofilewriter
        subroutine TimeStepTotals(itime)
          integer, intent(in) :: itime
        end subroutine TimeStepTotals
        subroutine dealloc_arrays
        end subroutine dealloc_arrays
      END INTERFACE

!      ! Before we do anything, start a log file
!      open(unit=global_log,file='Ash3d.lst',status='unknown')

      call Set_OS_Env

      aloft_percent_remaining = 1.0_ip
      SourceCumulativeVol     = 0.0_ip

      call cpu_time(t0) !time is a scaler real

        ! input data for ash transport
      call Read_Control_File

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to custom input blocks here
!
!  Loop through all the optional modules compiled and test against the list
!  from the input file (e.g. OPTMOD=TOPO), then call the special input reader
!  for that block
!  Do a sanity check on optional module requested in the input file v.s. those
!  compiled in this executable and for consistency among modules:
!  e.g. SRC_RESUSP will require the VARDIFF and LC be set
!
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Now looping through optional modules found in input file"
      endif;enddo
      do j=1,nmods
        do io=1,2;if(VB(io).le.verbosity_essential)then
          write(outlog(io),*)"Testing for ",OPTMOD_names(j),j
        endif;enddo
        if(OPTMOD_names(j).eq.'RESETPARAMS')then
          do io=1,2;if(VB(io).le.verbosity_essential)then
            write(outlog(io),*)"  Reading input block for RESETPARAMS"
          endif;enddo
          call input_data_ResetParams
        endif
      enddo
!
!------------------------------------------------------------------------------

        ! Read airports/POI and allocate/initilize arrays
        ! We only need to do this if an output variable demands it since this is
        ! a burden every time step
      if(Output_every_TS) &
        call ReadAirports

      call alloc_arrays
        ! Set up grids for solution and Met data
      call calc_mesh_params

      if (((SourceType.eq.'umbrella').or.(SourceType.eq.'umbrella_air')))then
        call Allocate_Source_Umbrella(nxmax,nymax,nzmax)
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Initialize concen and any special source terms here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(LoadConcen)then
        ! We are initializing the concentration and time from an output file
        ! Currently, Ash3d assumes the concentration file is compatible with
        ! the computational grid and grainsize distribution
#ifdef USENETCDF
        call NC_RestartFile_LoadConcen
#else
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: Loading concentration files requires previous netcdf"
          write(errlog(io),*)"       output.  This Ash3d executable was not compiled with"
          write(errlog(io),*)"       netcdf support.  Please recompile Ash3d with"
          write(errlog(io),*)"       USENETCDF=T, or select another source."
        endif;enddo
        stop 1
#endif
      else
        ! Initialize arrays if we haven't already loaded the concentration from
        ! a previous run
        concen_pd = 0.0_ip
        DepositGranularity = 0.0_ip
      endif
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert special source terms here
!
!------------------------------------------------------------------------------

      if(useTemperature)then
        call Allocate_Atmosphere_Met
      endif
      ! Now we can allocate the variables that live on the met grid
      if(useCalcFallVel)then
        call Allocate_Tephra_Met
      endif

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to optional variable allocation subroutines here
!
!------------------------------------------------------------------------------
      ! Allocate all the output variables
      call Allocate_Output_UserVars(nxmax,nymax,nzmax,nsmax)

      ! Now that we have the Met grids initialized, get the state variables
      ! interpoated on the start time
      time           = 0.0_ip
      Load_MesoSteps = .true.
      Interval_Frac  = 0.0_ip
      first_time     = .true.
      call MesoInterpolater(time , Load_MesoSteps , Interval_Frac, first_time)

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to special MesoInterpolaters subroutines here
!
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to prep user-specified output
!
!------------------------------------------------------------------------------

        ! Call output_results before time loop to create output files
      call output_results

      ntmax = max(1,3*int(Simtime_in_hours/dt))
      call Allocate_NTime(ntmax)
      if (Write_PR_Data)then
        call Allocate_Profile(nzmax,ntmax,nvprofiles)
      endif

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),7)
      endif;enddo

      call Allocate_Source_time

      !call MassFluxCalculator          !find current mass flux & plume height

      ! Write out starting volume, max time steps, and headers for the table that follows
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),1) tot_vol,ntmax
      endif;enddo

      ! ************************************************************************
      ! ****** begin time simulation *******************************************
      ! ************************************************************************
      itime = 0
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Starting time loop."
      endif;enddo

      do while (StopTimeLoop.eqv..false.)
        ! Note: stop conditions are evaluated at the end of the time loop

        ! Copy previous t=n+1 to t=n
        concen_pd(:,:,:,:,ts0) = concen_pd(:,:,:,:,ts1)
          ! re-initialize slice (1)
        concen_pd(:,:,:,1:nsmax,ts1) = 0.0_ip

        Called_Gen_Output_Vars  = .false.
        Calculated_Cloud_Load   = .false.
        Calculated_AshThickness = .false.

        itime = itime + 1
        if(itime.gt.ntmax)then
          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"WARNING: The number of time steps attempted exceeds 3x that anticipated."
            write(outlog(io),*)"         Check that the winds are stable"
            write(outlog(io),*)"        Simtime_in_hours = ",Simtime_in_hours
            write(outlog(io),*)"                   ntmax = ",ntmax
            write(outlog(io),*)"            current step = ",itime
          endif;enddo
        endif

          ! find the wind field at the current time
        first_time     = .false.
        call MesoInterpolater(time , Load_MesoSteps , Interval_Frac, first_time)

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to special MesoInterpolaters subroutines here
!
!------------------------------------------------------------------------------

        call MassFluxCalculator         ! call subroutine that determines mass flux & plume height

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to specialized MassFluxRate calculations here
!
!------------------------------------------------------------------------------

        ! Add source term
        if(MassFluxRate_now.gt.0.0_ip) then
          ! Check if the source type is one of the standard types with a 1-node column
          if ((SourceType.eq.'point')  .or. &
              (SourceType.eq.'line')   .or. &
              (SourceType.eq.'profile').or. &
              (SourceType.eq.'suzuki'))then

            ! Calculating the flux into the vent column
            call TephraSourceNodes

            ! Most standard source types (point, line, profile, suzuki) are
            ! integrated as follows.
            concen_pd(ivent,jvent,1:nzmax+1,1:n_gs_max,ts0) =  &
                concen_pd(ivent,jvent,1:nzmax+1,1:n_gs_max,ts0)    &
                  + dt*SourceNodeFlux(1:nzmax+1,1:n_gs_max)

            ! Keep track of the accumulated source inserted for mass conservation error-checking
            SourceCumulativeVol = SourceCumulativeVol + SourceVolInc(dt)

          elseif (SourceType.eq.'umbrella'.or. &
                 (SourceType.eq.'umbrella_air')) then
            ! Umbrella clouds have a special source insertion with a 3x3 column

            ! Umbrella sources still need the Suzuki distribution of mass above the vent
            ! Calculating the flux into the vent column
            call TephraSourceNodes

            ! Now modify the above result to distribute the total SourceNodeFlux
            ! into the umbrella form
            call TephraSourceNodes_Umbrella

            ! Before source insertion, we smooth over the concentration over the
            ! 3x3 patch within the umbrella zone using the weighted averaging stencil
            do iz=ibase,itop
              do isize=1,n_gs_max
                concen_pd(ivent-1:ivent+1,jvent-1:jvent+1,iz,isize,ts0) = &
                  AvgCon_Umbrella(concen_pd(ivent-1:ivent+1,jvent-1:jvent+1,iz,isize,ts0),iz)
              enddo
            enddo
            ! Here the integration is the same as for standard sources (point, line, profile, suzuki),
            ! except we use a 3x3 column around ivent,jvent instead of just a single-node column.
            concen_pd(ivent-1:ivent+1,jvent-1:jvent+1,1:nzmax+1,1:n_gs_max,ts0) =  &
                concen_pd(ivent-1:ivent+1,jvent-1:jvent+1,1:nzmax+1,1:n_gs_max,ts0)    &
                  + dt*SourceNodeFlux_Umbrella(1:3,1:3,1:nzmax+1,1:n_gs_max)

            ! Keep track of the accumulated source inserted for mass conservation error-checking
            SourceCumulativeVol = SourceCumulativeVol + SourceVolInc_Umbrella(dt)

          else
            ! This is not a standard source.
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"WARNING: source type is non-standard"
            endif;enddo
            stop 1
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to optional sources here
!         These subroutines need to calculate the mass of a species inserted
!         into specific cells and update concen accordingly
!
!          elseif () then
!
!------------------------------------------------------------------------------
          endif
        endif !MassFluxRate_now.gt.0.0_ip

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Set Boundary Conditions
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to optional boundary conditions here
!
!------------------------------------------------------------------------------
        call Set_BC(1)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Advection / Diffusion / Deposition
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to optional advection/diffusion routines here
!
!------------------------------------------------------------------------------

        if(useHorzAdvect) call AdvectHorz(itime)

        if(useVertAdvect) call advect_z

        if(useDiffusion)then
          call Set_BC(2)
          call DiffuseVert
          call DiffuseHorz(itime)
        endif

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to optional deposition routines here
!
!------------------------------------------------------------------------------

        ! Advance time
        time = time + dt

        ! If there is any time-series output requiring evaluation at each
        ! time-step, then extract output variables from concen here
        if(Output_every_TS)then
          call Gen_Output_Vars

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls output routines (every timestep) here
!
!------------------------------------------------------------------------------

            ! SEE WHETHER THE ASH HAS HIT ANY AIRPORTS
          call FirstAsh

            ! Track ash on vertical profiles
          if (Write_PR_Data)then
            call Calc_vprofile(itime)
            call vprofilewriter(itime)     !write out vertical profiles
          endif
        endif

        ! GO TO OUTPUT RESULTS IF WE'RE AT THE NEXT OUTPUT STAGE
        ! Note that dt was set in Adjust_DT so that it is no larger than
        ! DT_MIN, but may be adjusted down so as to land on the next
        ! output time.  time has already been integrated forward so
        ! NextWriteTime-time should be near zero for output steps.
        if(Output_at_WriteTimes.and.(NextWriteTime-time.lt.DT_MIN))then
            ! Generate output variables if we haven't already
          if(.not.Called_Gen_Output_Vars)then
            call Gen_Output_Vars
          endif
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls output routines (every output-step) here
!
!------------------------------------------------------------------------------
          call output_results
          !if ((WriteAirportFile_ASCII.or.WriteAirportFile_KML).and. &
          if (Write_PT_Data.and. &
              (iTimeNext.lt.nWriteTimes)) then
            do j=iTimeNext,nWriteTimes
              Airport_Thickness_TS(1:nairports,j) = Airport_Thickness(1:nairports)
            enddo
          endif
        endif

        if(Output_at_logsteps)then  
          !WRITE SUMMARY INFORMATION ON MASS CONSERVATION EVERY log_step TIME STEPS
          if(mod(itime,log_step).eq.0) then
            if(.not.Called_Gen_Output_Vars)then
              call Gen_Output_Vars
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls output routines (every log-step) here
!
!------------------------------------------------------------------------------
            endif
            call TimeStepTotals(itime)
          endif

            ! Only consider reducing GS bins at logsteps
          if(time.gt.e_EndTime_final)then
            ! Doesn't make sense to flag bins as flushed out while the eruption is on-going
            if(.not.Called_Gen_Output_Vars)then
              call Gen_Output_Vars
            endif
            call Prune_GS
          endif
        else
            ! If we are not monitoring deposits through logsteps, then set
            ! tot_vol to 0 and rely on the simulation ending through input
            ! duration values
          tot_vol = 0.0_ip
        endif

        if(tot_vol.gt.EPS_SMALL)then
          aloft_percent_remaining = aloft_vol/tot_vol
        else
          aloft_percent_remaining = 1.0_ip
        endif

        ! Check stop conditions
        !  If any of these is true, then the time loop will stop
           ! Stops if there is less than 1% of ash aloft in the domain
        StopConditions(1) = (aloft_percent_remaining.lt.(1.0_ip-StopValue))
           ! Normal stop condition if simulation exceeds alloted time
        StopConditions(2) = (time.ge.Simtime_in_hours)
           ! Normal stop condition when nothing is left to advect
        StopConditions(3) = (n_gs_aloft.eq.0)
        if(SourceCumulativeVol.gt.EPS_TINY)then
          MassConsErr = abs(SourceCumulativeVol-tot_vol)/SourceCumulativeVol
        endif
           ! Error stop condition if the concen and outflow do not match the source
        StopConditions(4) = (MassConsErr.gt.1.0e-3_ip)
        !StopConditions(4) = .false.  ! We override this condition until the umbrella source conserves mass
           ! Error stop condition if any volume measure is negative
        StopConditions(5) = (dep_vol.lt.-1.0_ip*EPS_SMALL).or.&
                            (aloft_vol.lt.-1.0_ip*EPS_SMALL).or.&
                            (outflow_vol.lt.-1.0_ip*EPS_SMALL).or.&
                            (SourceCumulativeVol.lt.-1.0_ip*EPS_SMALL)

        if((CheckConditions(1).eqv..true.).and.&
           (StopConditions(1).eqv..true.))then
          StopTimeLoop = .true.
        elseif((CheckConditions(2).eqv..true.).and.&
               (StopConditions(2).eqv..true.))then
          StopTimeLoop = .true.
        elseif((CheckConditions(3).eqv..true.).and.&
               (StopConditions(3).eqv..true.))then
          StopTimeLoop = .true.
        elseif((CheckConditions(4).eqv..true.).and.&
               (StopConditions(4).eqv..true.))then
          StopTimeLoop = .true.
        elseif((CheckConditions(5).eqv..true.).and.&
               (StopConditions(5).eqv..true.))then
          StopTimeLoop = .true.
        else
          StopTimeLoop = .false.
        endif
      enddo  !loop over itime
              !  ((dep_percent_accumulated.le.StopValue).and. &
              !    (time.lt.Simtime_in_hours)        .and. &
              !    (n_gs_aloft.gt.0))

      ! Reset ntmax to the actual number of time steps
      ntmax = itime

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Time integration completed for the following reason:"
      endif;enddo
      if((CheckConditions(1).eqv..true.).and.&
         (StopConditions(1).eqv..true.))then
        ! Normal stop condition set by user tracking the deposit
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Percent accumulated/exited exceeds ",StopValue
        endif;enddo
      endif
      if((CheckConditions(2).eqv..true.).and.&
         (StopConditions(2).eqv..true.))then
        ! Normal stop condition if simulation exceeds alloted time
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"time.ge.Simtime_in_hours"
          write(outlog(io),*)"              Time = ",real(time,kind=4)
          write(outlog(io),*)"  Simtime_in_hours = ",real(Simtime_in_hours,kind=4)
        endif;enddo
      endif
      if((CheckConditions(3).eqv..true.).and.&
         (StopConditions(3).eqv..true.))then
        ! Normal stop condition when nothing is left to advect
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"No ash species remain aloft."
        endif;enddo
      endif
      if((CheckConditions(4).eqv..true.).and.&
         (StopConditions(4).eqv..true.))then
        ! Error stop condition if the concen and outflow do not match the source
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"Cummulative source volume does not match aloft + outflow"
          write(errlog(io),*)" tot_vol = ",tot_vol
          write(errlog(io),*)" SourceCumulativeVol = ",SourceCumulativeVol
          write(errlog(io),*)" Abs. Error = ",&
                               abs((tot_vol-SourceCumulativeVol)/SourceCumulativeVol)
          write(errlog(io),*)" e_Volume = ",e_Volume
        endif;enddo
        stop 1
      endif
      if((CheckConditions(5).eqv..true.).and.&
         (StopConditions(5).eqv..true.))then
        ! Error stop condition if any volume measure is negative
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"One of the volume measures is negative."
          write(errlog(io),*)"        dep_vol = ",dep_vol
          write(errlog(io),*)"        aloft_vol = ",aloft_vol
          write(errlog(io),*)"        outflow_vol = ",outflow_vol
          write(errlog(io),*)"        SourceCumulativeVol = ",SourceCumulativeVol
        endif;enddo
        stop 1
      endif

      ! ************************************************************************
      ! ****** end time simulation *********************************************
      ! ************************************************************************

      isFinal_TS = .true.
      Called_Gen_Output_Vars  = .false.
      Calculated_Cloud_Load   = .false.
      Calculated_AshThickness = .false.

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),12)   !put footnotes below output table
        write(outlog(io),*)'time=',real(time,kind=4),',dt=',real(dt,kind=4)
        write(outlog(io),*)"Mass Conservation Error = ",MassConsErr
      endif;enddo

        ! Make sure we have the latest output variables and go to write routines
      call Gen_Output_Vars
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls output routines here
!
!------------------------------------------------------------------------------

      call output_results
      !WRITE RESULTS TO LOG AND STANDARD OUTPUT
      !TotalTime_sp = etime(elapsed_sp)
      call cpu_time(t1) !time is a scalar real
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),3) t1-t0, time*HR_2_S
        write(outlog(io),4) time*HR_2_S/(t1-t0)
      endif;enddo
      call TimeStepTotals(itime)
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),5) dep_vol
        write(outlog(io),6) tot_vol
        write(outlog(io),9) maxval(DepositThickness), DepositAreaCovered
        write(outlog(io),34)       !write out area of cloud at different thresholds
        do k=1,5
          write(outlog(io),35) LoadVal(k), CloudLoadArea(k)
        enddo
        write(outlog(io),33)    !write "normal completion"
      endif;enddo

      ! Format statements
1     format(/,5x,'Starting volume (km3 DRE)    = ',f11.4,       &
             /,5x,'maximum number of time steps = ',i8,          &
             //,21x,'Time',19x,                                  &
                '|--------------------Volume (km3 DRE)-------------------|', &
                3x,'Cloud Area',                                 &
              /,7x,'step',8x,'(hrs)',2x,'yyyymmddhh:mm',         &
                5x,'Source',6x,'Deposit',7x,'Aloft',5x,'Outflow',&
                7x,'Total',10x,'km2')

3     format(/,5x,'Execution time           = ',f15.4,' seconds',&
             /,5x,'Simulation time          = ',f15.4,' seconds')       
4     format(  5x,'Execution time/CPU time  = ',f15.4)       
5     format(  5x,'Ending deposit volume    = ',f15.4,' km3 DRE')       
6     format(  5x,'Ending total volume      = ',f15.4,' km3 DRE')       
7     format(  5x,'Building time array of plume height & eruption rate')
9     format(/,5x,'Maximum deposit thickness (mm)   = ',f10.4, &
             /,5x,'Area covered by >0.01 mm (km2)   = ',f10.1,/)
12    format(4x,'*=files written out')

33    format(/,5x,'Normal completion')
34    format(/,'  Ash load   cloud area',/, &
               '      T/km2         km2')
35    format(2f11.1)

      ! clean up memory
      call dealloc_arrays

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls deallocation routines here
!
!------------------------------------------------------------------------------

      close(global_log)       !close log file 

      end program Ash3d
