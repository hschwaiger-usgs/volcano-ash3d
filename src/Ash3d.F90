      program Ash3d

      use precis_param

      use io_units

      use global_param,  only : &
         useCalcFallVel,useDiffusion,useHorzAdvect,useVertAdvect,VERB,&
         HR_2_S,useTemperature,DT_MIN,KM3_2_M3,EPS_TINY,EPS_SMALL,&
         nmods,OPTMOD_names,StopConditions,CheckConditions

      use mesh,          only : &
         ivent,jvent,nxmax,nymax,nzmax,nsmax,ts0,ts1,kappa_pd

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
         ibase,itop,SourceNodeFlux,e_EndTime_final,e_Volume,MassFluxRate_now,&
         SourceType, &
           Allocate_Source_time,&
           MassFluxCalculator,&
           TephraSourceNodes

      use Tephra,        only : &
         n_gs_max,n_gs_aloft,MagmaDensity,&
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

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert 'use' statements here
!
!------------------------------------------------------------------------------

      implicit none

!      integer               :: iostatus
      integer               :: itime
      integer               :: j,k
      integer               :: ii,jj,iz,isize
      real(kind=ip)         :: avgcon        ! avg concen of cells in umbrella
      real(kind=ip)         :: Interval_Frac
      logical               :: Load_MesoSteps
      logical               :: StopTimeLoop   = .false.
      logical               :: first_time     = .true.
!      character(len=130)    :: tmp_str
      real(kind=ip)         :: MassConsErr

      INTERFACE
#ifdef USENETCDF
        !subroutine NC_RestartFile_LoadConcen
        !end subroutine NC_RestartFile_LoadConcen
#endif
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

      ! Before we do anything, start a log file
      open(unit=global_log,file='Ash3d.lst',status='unknown')

      call Set_OS_Env

      aloft_percent_remaining = 1.0_ip
      SourceCumulativeVol     = 0.0_ip

      call cpu_time(t0) !time is a scaler real

        ! input data for ash transport
      call Read_Control_File

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to custom input blocks here
      if(VERB.ge.1)write(global_info,*)&
        "Now looping through optional modules found in input file"
      do j=1,nmods
        if(VERB.ge.1)write(global_info,*)"Testing for ",OPTMOD_names(j),j
        if(OPTMOD_names(j).eq.'RESETPARAMS')then
          if(VERB.ge.1)write(global_info,*)"  Reading input block for RESETPARAMS"
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Initialize concen and any special source terms here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(LoadConcen)then
        ! We are initializing the concentration and time from an output file
        ! Currently, Ash3d assumes the concentration file is compatible with
        ! the computational grid and grainsize distribution
#ifdef USENETCDF
        !call NC_RestartFile_LoadConcen
#else
        write(global_error,*)"ERROR: Loading concentration files requires previous netcdf"
        write(global_error,*)"       output.  This Ash3d executable was not compiled with"
        write(global_error,*)"       netcdf support.  Please recompile Ash3d with"
        write(global_error,*)"       USENETCDF=T, or select another source."
        stop 1
#endif
      else
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

      ! write "Building time array of plume height & eruption rate"
      if(VERB.ge.1)write(global_info,7)
      if(VERB.ge.1)write(global_log ,7)

      call Allocate_Source_time

      !call MassFluxCalculator          !find current mass flux & plume height

      ! Write out starting volume, max time steps, and headers for the table that follows
      if(VERB.ge.1)write(global_info,1) tot_vol,ntmax
      if(VERB.ge.1)write(global_log ,1) tot_vol,ntmax

      ! ************************************************************************
      ! ****** begin time simulation *******************************************
      ! ************************************************************************
      itime = 0
      if(VERB.ge.1)write(global_info,*)"Starting time loop."

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
          if(VERB.ge.1)then
            write(global_info,*)"WARNING: The number of time steps attempted exceeds 3x that anticipated."
            write(global_info,*)"         Check that the winds are stable"
            write(global_info,*)"        Simtime_in_hours = ",Simtime_in_hours
            write(global_info,*)"                   ntmax = ",ntmax
            write(global_info,*)"            current step = ",itime
          endif
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
        ! erupt ash into column
        if(MassFluxRate_now.gt.0.0_ip) then
          ! Check if the source type is one of the standard types
          if ((SourceType.eq.'point')  .or. &
              (SourceType.eq.'line')   .or. &
              (SourceType.eq.'profile').or. &
              (SourceType.eq.'suzuki') .or. &
              (SourceType.eq.'umbrella') .or. &
              (SourceType.eq.'umbrella_air'))then
            ! Calculating the flux at the source nodes
            call TephraSourceNodes

            ! Now integrate the ash concentration with the SourceNodeFlux
            if (SourceType.eq.'umbrella'.or. &
               (SourceType.eq.'umbrella_air')) then
              ! Umbrella clouds have a special integration
              !  Below the umbrella cloud, add ash to vent nodes as above
              concen_pd(ivent,jvent,1:ibase-1,1:n_gs_max,ts0) =          & ! 
                       concen_pd(ivent,jvent,1:ibase-1,1:n_gs_max,ts0) + & ! kg/km3
                       dt                                              * & ! hr
                       SourceNodeFlux(1:ibase-1,1:n_gs_max)                ! kg/km3 hr
              do isize=1,n_gs_max
                do k=1,ibase-1
                  SourceCumulativeVol = SourceCumulativeVol + & ! final units is km3
                    dt                              * & ! hr
                    SourceNodeFlux(k,isize)         * & ! kg/km3 hr
                    kappa_pd(ivent,jvent,k)         / & ! km3
                    MagmaDensity                    / & ! kg/m3
                    KM3_2_M3                            ! m3/km3
                enddo
              enddo
              do iz=ibase,itop
                !Within the cloud: first, average the concentration that curently
                !exists in the 9 cells surrounding the vent
                do isize=1,n_gs_max
                  avgcon=sum(concen_pd(ivent-1:ivent+1,jvent-1:jvent+1,iz,isize,ts0))/9.0_ip
                  concen_pd(ivent-1:ivent+1,jvent-1:jvent+1,iz,isize,ts0)=avgcon
                enddo
              enddo
              !Then, add tephra to the 9 nodes surrounding the vent
              ! TephraSourceNodes has a special line to reduce SourceNodeFlux by a factor 9
              ! because it is applied 9 times here.  We need to be careful about mixing mass
              ! and concentration since cell volume differ in lat, but this should be minor
              do ii=ivent-1,ivent+1
                do jj=jvent-1,jvent+1
                  do iz=ibase,itop
                    concen_pd(ii,jj,iz,1:n_gs_max,ts0) =                &
                              concen_pd(ii,jj,iz,1:n_gs_max,ts0)        &
                                 + dt*SourceNodeFlux(iz,1:n_gs_max)
                    do isize=1,n_gs_max
                      SourceCumulativeVol = SourceCumulativeVol + & ! final units is km3
                        dt                              * & ! hr
                        SourceNodeFlux(iz,isize)         * & ! kg/km3 hr
                        kappa_pd(ivent,jvent,iz)         / & ! km3
                        MagmaDensity                    / & ! kg/m3
                        KM3_2_M3                            ! m3/km3
                    enddo
                  enddo
                enddo
              enddo
              !! this part is just for book-keeping and error checking
              !do isize=1,n_gs_max
              !  do k=1,nzmax+1
              !    SourceCumulativeVol = SourceCumulativeVol + & ! final units is km3
              !      dt                              * & ! hr
              !      SourceNodeFlux(k,isize)         * & ! kg/km3 hr
              !      kappa_pd(ivent,jvent,k)         / & ! km3
              !      MagmaDensity                    / & ! kg/m3
              !      KM3_2_M3                            ! m3/km3
              !  enddo
              !enddo
            else ! (SourceType.eq.'umbrella' or 'umbrella_air')
              ! All other standard source types (point,line,profile, suzuki) are
              ! integrated as follows.
              concen_pd(ivent,jvent,1:nzmax+1,1:n_gs_max,ts0) =  &
              concen_pd(ivent,jvent,1:nzmax+1,1:n_gs_max,ts0)    &
                + dt*SourceNodeFlux(1:nzmax+1,1:n_gs_max)
              ! this part is just for book-keeping and error checking
              do isize=1,n_gs_max
                do k=1,nzmax+1
                  SourceCumulativeVol = SourceCumulativeVol + & ! final units is km3
                    dt                              * & ! hr
                    SourceNodeFlux(k,isize)         * & ! kg/km3 hr
                    kappa_pd(ivent,jvent,k)         / & ! km3
                    MagmaDensity                    / & ! kg/m3
                    KM3_2_M3                            ! m3/km3
                enddo
              enddo
            endif
          else
            ! This is not a standard source.
            if(VERB.ge.1)write(global_info,*)"WARNING: source type is non-standard"
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
        StopConditions(4) = .false.
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

      if(VERB.ge.1)write(global_info,*)"Time integration completed for the following reason:"
      if((CheckConditions(1).eqv..true.).and.&
         (StopConditions(1).eqv..true.))then
        ! Normal stop condition set by user tracking the deposit
        if(VERB.ge.1)write(global_info,*)"Percent accumulated/exited exceeds ",StopValue
      endif
      if((CheckConditions(2).eqv..true.).and.&
         (StopConditions(2).eqv..true.))then
        ! Normal stop condition if simulation exceeds alloted time
        if(VERB.ge.1)then
          write(global_info,*)"time.ge.Simtime_in_hours"
          write(global_info,*)"              Time = ",real(time,kind=4)
          write(global_info,*)"  Simtime_in_hours = ",real(Simtime_in_hours,kind=4)
          write(global_log,*)"time.ge.Simtime_in_hours"
          write(global_log,*)"              Time = ",real(time,kind=4)
          write(global_log,*)"  Simtime_in_hours = ",real(Simtime_in_hours,kind=4)
        endif
      endif
      if((CheckConditions(3).eqv..true.).and.&
         (StopConditions(3).eqv..true.))then
        ! Normal stop condition when nothing is left to advect
        if(VERB.ge.1)then
          write(global_info,*)"No ash species remain aloft."
          write(global_log,*)"No ash species remain aloft."
        endif
      endif
      if((CheckConditions(4).eqv..true.).and.&
         (StopConditions(4).eqv..true.))then
        ! Error stop condition if the concen and outflow do not match the source
        write(global_error,*)"Cummulative source volume does not match aloft + outflow"
        write(global_error,*)" tot_vol = ",tot_vol
        write(global_error,*)" SourceCumulativeVol = ",SourceCumulativeVol
        write(global_error,*)" Abs. Error = ",&
                             abs((tot_vol-SourceCumulativeVol)/SourceCumulativeVol)
        write(global_error,*)" e_Volume = ",e_Volume
        stop 1
      endif
      if((CheckConditions(5).eqv..true.).and.&
         (StopConditions(5).eqv..true.))then
        ! Error stop condition if any volume measure is negative
        write(global_error,*)"One of the volume measures is negative."
        write(global_error,*)"        dep_vol = ",dep_vol
        write(global_error,*)"        aloft_vol = ",aloft_vol
        write(global_error,*)"        outflow_vol = ",outflow_vol
        write(global_error,*)"        SourceCumulativeVol = ",SourceCumulativeVol
        stop 1
      endif

      ! ************************************************************************
      ! ****** end time simulation *********************************************
      ! ************************************************************************

      isFinal_TS = .true.
      Called_Gen_Output_Vars  = .false.
      Calculated_Cloud_Load   = .false.
      Calculated_AshThickness = .false.

      if(VERB.ge.1)then
        write(global_info,12)   !put footnotes below output table
        write(global_log,12)   !put footnotes below output table
        write(global_log ,12)
        write(global_info,*)'time=',real(time,kind=4),',dt=',real(dt,kind=4)
        write(global_log,*)'time=',real(time,kind=4),',dt=',real(dt,kind=4)
        write(global_info,*)"Mass Conservation Error = ",MassConsErr
        write(global_log,*)"Mass Conservation Error = ",MassConsErr
      endif

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
      if(VERB.ge.1)then
        write(global_info,3) t1-t0, time*HR_2_S
        write(global_log ,3) t1-t0, time*HR_2_S
        write(global_info,4) time*HR_2_S/(t1-t0)
        write(global_log ,4) time*HR_2_S/(t1-t0)
      endif
      call TimeStepTotals(itime)
      if(VERB.ge.1)then
        write(global_info,5) dep_vol
        write(global_log ,5) dep_vol
        write(global_info,6) tot_vol
        write(global_log ,6) tot_vol
        write(global_info,9) maxval(DepositThickness), DepositAreaCovered
        write(global_log ,9) maxval(DepositThickness), DepositAreaCovered

        write(global_info,34)       !write out area of cloud at different thresholds
        write(global_log ,34)
        do k=1,5
          write(global_info,35) LoadVal(k), CloudLoadArea(k)
          write(global_log ,35) LoadVal(k), CloudLoadArea(k)
        enddo
        
        write(global_info,33)    !write "normal completion"
        write(global_log ,33)
      endif

      close(9)       !close log file 

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

      end program Ash3d
