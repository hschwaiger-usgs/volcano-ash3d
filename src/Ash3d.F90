      program Ash3d

      use precis_param

      use io_units

      use global_param,  only : &
         EPS_TINY,useCalcFallVel,useDiffusion,useHorzAdvect,useVertAdvect,VERB,&
         HR_2_S,useTemperature

      use mesh,          only : &
         ivent,jvent,nxmax,nymax,nzmax,nsmax,ts0,ts1

      use solution,      only : &
         concen_pd,dep_vol,tot_vol,DepositGranularity,StopValue,&
         vx_pd,vy_pd,vz_pd

      use Output_Vars,   only : &
         AreaCovered,DepositThickness,LoadVal,CloudLoadArea,&
           Allocate_Output_Vars, &
           Allocate_Output_UserVars, &
           Gen_Output_Vars,&
           FirstAsh

      use io_data,       only : &
         Called_Gen_Output_Vars,isFinal_TS,LoadConcen,log_step, &
         Output_at_logsteps,Output_at_WriteTimes,Output_every_TS,&
         NextWriteTime,nTimeNext,nvprofiles,nWriteTimes,&
         WriteAirportFile_ASCII,WriteAirportFile_KML

      use time_data,     only : &
         time,dt,Simtime_in_hours,t0,t1

      use Source,        only : &
         ibase,itop,SourceNodeFlux,e_EndTime_final,e_Volume,MassFluxRate_now,&
         SourceType, &
           Allocate_Source_time,&
           MassFluxCalculator,&
           TephraSourceNodes

      use Tephra,        only : &
         n_gs_max,ns_aloft,&
           Allocate_Tephra,&
           Allocate_Tephra_Met,&
           Collapse_GS

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

      integer               :: itime
      integer               :: j,k
      integer               :: ii,jj,iz,isize
      real(kind=ip)         :: avgcon        ! avg concen of cells in umbrella
      real(kind=ip)         :: percent_accumulated
      real(kind=ip)         :: Interval_Frac
      logical               :: Load_MesoSteps
      integer               :: ntmax
      logical, dimension(5) :: StopConditions = .false.
      logical               :: StopTimeLoop   = .false.
      logical               :: first_time     = .true.

      percent_accumulated = 0.0_ip

      call cpu_time(t0) !time is a scaler real

        ! input data for ash transport
      call Read_Control_File
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to custom input blocks here
!
!------------------------------------------------------------------------------

        ! Read airports/POI and allocate/initilize arrays
      if (WriteAirportFile_ASCII.or.WriteAirportFile_KML)  call ReadAirports

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
        call NC_RestartFile_LoadConcen
#else
        write(global_info,*)"ERROR: Loading concentration files requires previous netcdf"
        write(global_info,*)"       output.  This Ash3d executable was not compiled with"
        write(global_info,*)"       netcdf support.  Please recompile Ash3d with"
        write(global_info,*)"       USENETCDF=T, or select another source."
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

      call Adjust_DT

!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to prep user-specified output
!
!------------------------------------------------------------------------------

        ! Call output_results before time loop to create output files
      call output_results

      ntmax = max(1,int(Simtime_in_hours/dt))

      write(global_info,7)            !write "Building time array of plume height & eruption rate"
      write(global_log ,7)

      call Allocate_Source_time

      !call MassFluxCalculator          !find current mass flux & plume height

      ! Write out starting volume, max time steps, and headers for the table that follows
      write(global_info,1) tot_vol,ntmax
      write(global_log ,1) tot_vol,ntmax

      ! ************************************************************************
      ! ****** begin time simulation *******************************************
      ! ************************************************************************
      itime = 0
      write(global_info,*)"Starting time loop."
      do while (StopTimeLoop.eqv..false.)
        ! Note: stop conditions are evaluated at the end of the time loop

        ! Copy previous t=n+1 to t=n
        concen_pd(:,:,:,:,ts0) = concen_pd(:,:,:,:,ts1)
          ! re-initialize slice (1)
        concen_pd(:,:,:,1:nsmax,ts1) = 0.0_ip

        Called_Gen_Output_Vars = .false.

        itime = itime + 1
        if(itime.gt.3*ntmax)then
          write(global_info,*)"WARNING: The number of time steps attempted exceeds 3x that anticipated."
          write(global_info,*)"         Check that the winds are stable"
          write(global_info,*)"        Simtime_in_hours = ",Simtime_in_hours
          write(global_info,*)"                   ntmax = ",ntmax
          write(global_info,*)"            current step = ",itime
        endif

          ! find the wind field at the current time
        first_time     = .false.
        call MesoInterpolater(time , Load_MesoSteps , Interval_Frac, first_time)
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to special MesoInterpolaters subroutines here
!
!------------------------------------------------------------------------------

        call Adjust_DT

        if(VERB.gt.1)write(global_info,*)"Ash3d: Calling MassFluxCalculator"
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
              (SourceType.eq.'umbrella'))then
            ! Calculating the flux at the source nodes
            call TephraSourceNodes

            ! Now integrate the ash concentration with the SourceNodeFlux
            if (SourceType.eq.'umbrella') then
              ! Umbrella clouds have a special integration
              !  Below the umbrella cloud, add ash to vent nodes as above
              concen_pd(ivent,jvent,1:ibase-1,1:n_gs_max,ts0) =                        &
                                       concen_pd(ivent,jvent,1:ibase-1,1:n_gs_max,ts0) &
                                       + dt*SourceNodeFlux(1:ibase-1,1:n_gs_max)
              do iz=ibase,itop
                !Within the cloud: first, average the concentration that curently
                !exists in the 9 cells surrounding the vent
                do isize=1,n_gs_max
                  avgcon=sum(concen_pd(ivent-1:ivent+1,jvent-1:jvent+1,iz,isize,ts0))/9.0_ip
                  concen_pd(ivent-1:ivent+1,jvent-1:jvent+1,iz,isize,ts0)=avgcon
                enddo
              enddo
              !Then, add tephra to the 9 nodes surrounding the vent
              do ii=ivent-1,ivent+1
                do jj=jvent-1,jvent+1
                  do iz=ibase,itop
                    concen_pd(ii,jj,iz,1:nsmax,ts0) =                &
                              concen_pd(ii,jj,iz,1:n_gs_max,ts0)        &
                                 + dt*SourceNodeFlux(iz,1:nsmax)
                  enddo
                enddo
              enddo
            else
              ! All other standard source types (point,line,profile, suzuki) are
              ! integrated as follows.
              concen_pd(ivent,jvent,1:nzmax+1,1:n_gs_max,ts0) =  &
              concen_pd(ivent,jvent,1:nzmax+1,1:n_gs_max,ts0)    &
                + dt*SourceNodeFlux(1:nzmax+1,1:n_gs_max)
            endif
          else
            ! This is not a standard source.
            write(global_info,*)"WARNING: source type is non-standard"
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
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Set Boundary Conditions
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------
!       OPTIONAL MODULES
!         Insert calls to optional boundary conditions here
!
!------------------------------------------------------------------------------
        call Set_BC

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
          if (nvprofiles.gt.0) call vprofilewriter     !write out vertical profiles
        endif

        !GO TO OUTPUT RESULTS IF WE'RE AT THE NEXT OUTPUT STAGE
        if(Output_at_WriteTimes.and.(NextWriteTime-time.lt.1.0e-4_ip))then
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
          if ((WriteAirportFile_ASCII.or.WriteAirportFile_KML).and. &
              (nTimeNext.lt.nWriteTimes)) then
            do j=nTimeNext,nWriteTimes
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
            call Collapse_GS
          endif
        else
            ! If we are not monitoring deposits through logsteps, then set
            ! tot-vol to 0 and rely on the simulation ending through input
            ! duration values
          tot_vol = 0.0_ip
        endif

        if(tot_vol.gt.EPS_TINY)then
          percent_accumulated = dep_vol/tot_vol
        else
          percent_accumulated = 0.0_ip
        endif

        ! Check stop conditions
        !  If any of these is true, then the time loop will stop
        StopConditions(1) = (percent_accumulated.gt.StopValue)
        StopConditions(2) = (time.ge.Simtime_in_hours)
        StopConditions(3) = (ns_aloft.eq.0)
        StopConditions(4) = .false.  !tot_vol.le.(1.05_ip*sum(e_Volume))
        StopConditions(5) = .false.

        if(StopConditions(1).eqv..true.)then
          StopTimeLoop = .true.
        elseif(StopConditions(2).eqv..true.)then
          StopTimeLoop = .true.
        elseif(StopConditions(3).eqv..true.)then
          StopTimeLoop = .true.
        elseif(StopConditions(4).eqv..true.)then
          StopTimeLoop = .true.
        elseif(StopConditions(5).eqv..true.)then
          StopTimeLoop = .true.
        else
          StopTimeLoop = .false.
        endif

      enddo  !loop over itime
              !  ((percent_accumulated.le.StopValue).and. &
              !    (time.lt.Simtime_in_hours)        .and. &
              !    (ns_aloft.gt.0))

      write(global_info,*)"Time integration completed for the following reason:"
      if(.not.percent_accumulated.le.StopValue)then
        write(global_info,*)"Percent accumulated exceeds ",StopValue
      endif
      if(.not.time.lt.Simtime_in_hours)then
        write(global_info,*)"time.le.Simtime_in_hours"
        write(global_info,*)"              Time = ",time
        write(global_info,*)"  Simtime_in_hours = ",Simtime_in_hours
      endif
      if(ns_aloft.eq.0)then
        write(global_info,*)"ns_aloft = 0"
      endif
      if(tot_vol.gt.(1.05_ip*sum(e_Volume)))then
        write(global_info,*) "tot_vol>1.05*e_volume"
        write(global_info,*) " tot_vol = ",tot_vol
        write(global_info,*) "e_Volume = ",e_Volume
      endif
      ! ************************************************************************
      ! ****** end time simulation *********************************************
      ! ************************************************************************

      isFinal_TS = .true.
      write(global_info,12)   !put footnotes below output table
      write(global_log ,12)
      write(global_info,*) 'time=',time,', dt=',dt

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
      !write(global_info,*) elapsed_sp(2), time*3600.0_ip
      call cpu_time(t1) !time is a scalar real 
      write(global_info,3) t1-t0, time*HR_2_S
      write(global_log ,3) t1-t0, time*HR_2_S
      write(global_info,4) time*HR_2_S/(t1-t0)
      write(global_log ,4) time*HR_2_S/(t1-t0)
      call TimeStepTotals
      write(global_info,5) dep_vol
      write(global_log ,5) dep_vol
      write(global_info,6) tot_vol
      write(global_log ,6) tot_vol
      write(global_info,9) maxval(DepositThickness), AreaCovered
      write(global_log ,9) maxval(DepositThickness), AreaCovered

      write(global_info,34)       !write out area of cloud at different thresholds
      write(global_log ,34)
      do k=1,5
         write(global_info,35) LoadVal(k), CloudLoadArea(k)
         write(global_log ,35) LoadVal(k), CloudLoadArea(k)
      enddo
        
      write(global_info,33)    !write "normal completion"
      write(global_log ,33)
      
      close(9)       !close log file 

      ! Format statements
1     format(/,5x,'Starting volume (km3 DRE)    = ',f11.4,       &
             /,5x,'maximum number of time steps = ',i8,          &
             //,21x,'Time',19x,                                  &
                '|--------------Volume (km3 DRE)-------------|', &
                3x,'Cloud Area',                                 &
              /,7x,'step',8x,'(hrs)',2x,'yyyymmddhh:mm',         &
                5x,'Deposit',7x,'Aloft',5x,'Outflow',7x,          &
                   'Total',10x,'km2')
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
