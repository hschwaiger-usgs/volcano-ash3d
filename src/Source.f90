!##############################################################################
!
!  Source module
!
!  This module sets up the basic structure for all sources in Ash3d and the
!  standard sources in the column above the vent (line,Suzuki,profile,point).
!  Specialized sources, such the umbrella cloud or surface area sources,
!  can be included as separate modules, but generally rely on variables
!  provided by this module.
!
!      subroutine Allocate_Source_eruption
!      subroutine Allocate_Source_grid
!      subroutine Deallocate_Source
!      subroutine CheckEruptivePulses
!      subroutine TephraSourceNodes
!      function SourceVolInc
!
!  Note: source type is given on line 8 of block 1 of the control file
!        For a Suzuki source, we would have just the Suzuki constant:
!0.0       4.0                     # diffusion coefficient (m2/s), Suzuki constant
!        For a line, we would have:
!0.0       line                    # diffusion coefficient (m2/s), Suzuki constant
!        For a point, we would have:
!0.0       point                   # diffusion coefficient (m2/s), Suzuki constant
!        For an umbrella or umbrella air, we would have:
!0.0       umbrella[_air]            # diffusion coefficient (m2/s), Suzuki constant
!        And for a profile, we would have:
!0.0       profile                 # diffusion coefficient (m2/s), Suzuki constant
!
!  All standard sources have the following format:
!# ERUPTION LINES (number = neruptions)
!# In the following line, each line represents one eruptive pulse.
!# Parameters are (1-4) start time (yyyy mm dd h.hh (UT)); (5) duration (hrs);
!#                  (6) plume height (km);                 (7) eruped volume (km3)
!******************* BLOCK 2 ***************************************************
!2010 04 14   0.00   1.0     18.0  0.16
!
! The profile option has the additional two values for the dz and nz
!
!# Parameters are (1-4) start time (yyyy mm dd h.hh (UT)); (5) duration (hrs);
!#                  (6) plume height (km);                 (7) eruped volume (km3)
!#                  (8) dz of profile;                     (9) number of z segments
!******************* BLOCK 2 ***************************************************
!2010 04 14   0.00   1.0     18.0  0.16 1.0 18
!0.019 0.023 0.027 0.031 0.037 0.042 0.049 0.055 0.063 0.070 0.077 0.083 0.088 0.090 0.087 0.078 0.058 0.023
!
!##############################################################################

      module Source

      use precis_param

      use io_units

      implicit none

        ! Set everything to private by default
      private

        ! Publicly available subroutines/functions
      public Allocate_Source_eruption,  &
             Allocate_Source_grid,      &
             Deallocate_Source,         &
             Calc_Normalized_SourceCol, &
             EruptivePulse_MassFluxRate,&
             CheckEruptivePulses,       &
             TephraSourceNodes,         &
             SourceVolInc

        ! Publicly available variables

      real(kind=ip),public :: x_volcano, y_volcano      ! x & y points of volcano (km)
      real(kind=ip),public :: lat_volcano, lon_volcano  ! position of volcano in lat/lon
      real(kind=ip),public :: z_volcano                 ! vent elevation (km)

      integer,          public :: neruptions            ! number of eruptions or eruptive pulses
      character(len=12),public :: SourceType            ! may be 'point', 'line', or 'Suzuki' 
      real(kind=ip),    public :: Suzuki_A
      logical,          public :: IsCustom_SourceType = .false.

#ifdef USEPOINTERS
      real(kind=ip), dimension(:,:)  ,pointer,public :: NormSourceColumn    => null()
      real(kind=ip), dimension(:  )  ,pointer,public :: TephraFluxRate      => null()
      real(kind=ip), dimension(:,:)  ,pointer,public :: SourceNodeFlux      => null()
      real(kind=ip), dimension(:,:,:),pointer        :: SourceNodeFlux_Area => null()
#else
      real(kind=ip), dimension(:,:)  ,allocatable,public :: NormSourceColumn
      real(kind=ip), dimension(:)    ,allocatable,public :: TephraFluxRate
      real(kind=ip), dimension(:,:)  ,allocatable,public :: SourceNodeFlux
      real(kind=ip), dimension(:,:,:),allocatable        :: SourceNodeFlux_Area
#endif

        !The following arrays are used by MassFluxCalculator
      logical,       public :: Source_in_dt        ! true if any eruption contributes in this dt

        !The following arrays are of length neruptions
      real(kind=ip), dimension(:)        ,allocatable,public :: e_PlumeHeight
      real(kind=ip), dimension(:)        ,allocatable,public :: e_Volume
      real(kind=dp), dimension(:)        ,allocatable,public :: e_Duration   ! Time needs to be dp
      real(kind=dp), dimension(:)        ,allocatable,public :: e_StartTime  ! 
      real(kind=dp), dimension(:)        ,allocatable,public :: e_EndTime    ! 
      real(kind=ip), dimension(:)        ,allocatable,public :: e_prof_dz
      integer      , dimension(:)        ,allocatable,public :: e_prof_nzpoints
      real(kind=ip), dimension(:,:)      ,allocatable,public :: e_prof_Volume
      real(kind=ip), dimension(:,:)      ,allocatable,public :: e_prof_MassFluxRate
      real(kind=ip), dimension(:)        ,allocatable,public :: MassFluxRate

      real(kind=dp), dimension(:)        ,allocatable        :: dt_pulse_frac

      real(kind=dp), public :: e_EndTime_final

      real(kind=ip), public :: ESP_height        = 0.0_ip
      real(kind=dp), public :: ESP_duration      = 0.0_ip
      real(kind=ip), public :: ESP_MassFluxRate  = 0.0_ip
      real(kind=ip), public :: ESP_Vol           = 0.0_ip
      real(kind=ip), public :: ESP_massfracfine  = 0.0_ip

      integer, parameter :: MAXCustSrc = 10    ! The maximum number of custom
                                               ! source types that we will check for
      character(len=30),dimension(MAXCustSrc) :: SourceType_Custom = ""

      integer :: ieruption          ! eruption at the start of the time step
      integer :: jeruption          ! eruption at the end of the time step

      contains
      !------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Allocate_Source_eruption
!
!  Called from: Read_Control_File
!  Arguments:
!    none
!
!  This subroutine allocated the main variables needed to specify the source.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Allocate_Source_eruption

      allocate (e_StartTime(neruptions));            e_StartTime   = 0.0_ip
      allocate (e_Duration(neruptions));             e_Duration    = 0.0_ip
      allocate (e_PlumeHeight(neruptions));          e_PlumeHeight = 0.0_ip
      allocate (e_Volume(neruptions));               e_Volume      = 0.0_ip
      allocate (MassFluxRate(neruptions));           MassFluxRate  = 0.0_ip
      allocate (e_EndTime(neruptions));              e_EndTime     = 0.0_ip
      allocate (dt_pulse_frac(neruptions));          dt_pulse_frac = 0.0_dp

      if(SourceType.eq.'profile')then
        allocate (e_prof_dz(neruptions));             e_prof_dz       = 0.0_ip
        allocate (e_prof_nzpoints(neruptions));       e_prof_nzpoints = 0
          ! for profiles, assume 50 points
        allocate (e_prof_Volume(neruptions,50));      e_prof_Volume   = 0.0_ip
        allocate (e_prof_MassFluxRate(neruptions,50));e_prof_MassFluxRate = 0.0_ip
      endif

      ieruption = 1 ! Initialize eruption for the start of this dt to the starting eruption
      jeruption = 1 ! Initialize eruption for the end of this dt to the starting eruption

      end subroutine Allocate_Source_eruption

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Allocate_Source_grid
!
!  Called from: alloc_arrays
!  Arguments:
!    none
!
!  This subroutine source variables that are a function of z
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Allocate_Source_grid

      use mesh,          only : &
         nzmax,nsmax

      allocate(NormSourceColumn(neruptions,1:nzmax));    NormSourceColumn = 0.0_ip
      allocate(SourceNodeFlux(0:nzmax+1,1:nsmax));       SourceNodeFlux   = 0.0_ip
      allocate(TephraFluxRate(nzmax));                   TephraFluxRate   = 0.0_ip

      end subroutine Allocate_Source_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Deallocate_Source
!
!  Called from: dealloc_arrays
!  Arguments:
!    none
!
!  This subroutine deallocates variables allocated in Allocate_Source_eruption
!  and Allocate_Source_grid.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Deallocate_Source

      if(allocated(e_StartTime))      deallocate(e_StartTime)
      if(allocated(e_Duration))       deallocate(e_Duration)
      if(allocated(e_PlumeHeight))    deallocate(e_PlumeHeight)
      if(allocated(e_Volume))         deallocate(e_Volume)
      !if(allocated(MassFlux))         deallocate(MassFlux)
      if(allocated(MassFluxRate))     deallocate(MassFluxRate)
      if(allocated(e_EndTime))        deallocate(e_EndTime)

      ! SourceType.eq.'profile'
      if(allocated(e_prof_dz))           deallocate(e_prof_dz)
      if(allocated(e_prof_nzpoints))     deallocate(e_prof_nzpoints)
      if(allocated(e_prof_Volume))       deallocate(e_prof_Volume)
      if(allocated(e_prof_MassFluxRate)) deallocate(e_prof_MassFluxRate)

      if(allocated(NormSourceColumn))    deallocate(NormSourceColumn)
      if(allocated(SourceNodeFlux))      deallocate(SourceNodeFlux)
      if(allocated(TephraFluxRate))      deallocate(TephraFluxRate)

      end subroutine Deallocate_Source

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Calc_Normalized_SourceCol
!
!  Called from: Read_Control_File
!  Arguments:
!    none
!
!  This subroutine calculates the normalize representation of the source
!  column descretized to the z-grid of this Ash3d run.  This could be the
!  Suzuki distribution, a line or point source, or a profile specified on
!  a z-increment different than the dz of the simulation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Calc_Normalized_SourceCol

      use global_param,    only : &
         EPS_SMALL

      use mesh,          only : &
         nzmax,dz_vec_pd,z_lb_pd,z_cc_pd

      integer :: i
      integer :: k
      real(kind=ip) :: Suzuki_k     ! k factor in the Suzuki equation (see Hurst's Ashfall manual)
      real(kind=ip) :: z_cell_bot
      real(kind=ip) :: z_cell_top
      real(kind=ip) :: zbot_prof
      real(kind=ip) :: ztop_prof
      real(kind=ip) :: zground
      real(kind=ip) :: PlumeHeight_above_ground
      real(kind=ip) :: frac
      real(kind=ip) :: tot
      integer       :: kground
      integer       :: kk
      integer       :: kPlumeTop

      NormSourceColumn     = 0.0_ip

      do i=1,neruptions
        ! Get the cell containing the bottom of the source
        !  as well as the first cell above the plume top
        ! if(useTopo) kground = topo_indx(ivent,jvent)
        kground   = 1
        kPlumeTop = 0
        do k=1,nzmax+1
          if (z_volcano.ge.z_lb_pd(k  ).and. &
              z_volcano.lt.z_lb_pd(k+1))then
            kground = k
          endif
          if(e_PlumeHeight(i).ge.z_lb_pd(k  ).and. &
             e_PlumeHeight(i).lt.z_lb_pd(k+1))then
            kPlumeTop = k
          endif
        enddo
        ! Set ground level to be the lower boundary of cell kground
        zground = z_lb_pd(kground)

        if ((SourceType.eq.'suzuki')      .or. &
            (SourceType.eq.'umbrella')    .or. &
            (SourceType.eq.'umbrella_air')) then
          Suzuki_k = Suzuki_A/((e_PlumeHeight(i)-zground)* &
                    ((1.0_ip/Suzuki_A)-((Suzuki_A+1.0_ip)/Suzuki_A)* &
                    exp(-Suzuki_A)))
        endif

        ! Find the fraction of the erupted mass rate (kg/hr) that lies within
        ! each height interval dz
        do k=kground,kPlumeTop
          ! height at the top of this cell (or top of plume)
          z_cell_top  = min(z_cc_pd(k)+0.5_ip*dz_vec_pd(k),e_PlumeHeight(i))
          ! height at the bottom
          z_cell_bot = z_cc_pd(k)-0.5_ip*dz_vec_pd(k)

          ! First get the TephraFluxRate in kg/hr (total mass of tephra inserted per hour)
          PlumeHeight_above_ground = e_PlumeHeight(i)-zground
          if ((SourceType.eq.'suzuki')      .or. &
              (SourceType.eq.'umbrella')    .or. &
              (SourceType.eq.'umbrella_air')) then
            ! For Suzuki plumes and umbrella clouds
            ! It uses an equation obtained by integrating the 
            ! Suzuki equation given in Hurst.
            NormSourceColumn(i,k) = (Suzuki_k*                                  &
                               PlumeHeight_above_ground/Suzuki_A) *             &
                              ((1.0_ip+(1.0_ip/Suzuki_A)-((z_cell_top -zground)/&
                                PlumeHeight_above_ground)) *                    &
                                 exp(Suzuki_A *((z_cell_top -zground)/          &
                                  PlumeHeight_above_ground-1.0_ip)) -           &
                               (1.0_ip+(1.0_ip/Suzuki_A)-((z_cell_bot-zground)/ &
                                PlumeHeight_above_ground)) *                    &
                                 exp(Suzuki_A *((z_cell_bot-zground)/           &
                                  PlumeHeight_above_ground-1.0_ip)))
          elseif (SourceType.eq.'line') then
            ! For line sources, the fractional contribution of the cell
            ! at z is just the height of the cell over the length of the line
            NormSourceColumn(i,k) = (z_cell_top-z_cell_bot) / PlumeHeight_above_ground
          elseif (SourceType.eq.'point') then
            !for point sources, put all the contribution into the cell that contains
            ! the point
            if ((z_cell_top.ge.e_PlumeHeight(i)).and.    &
                (z_cell_bot.lt.e_PlumeHeight(i))) then
              NormSourceColumn(i,k) = 1.0_ip
            else
              NormSourceColumn(i,k) = 0.0_ip
            endif
          elseif (SourceType.eq.'profile') then
            ! loop over the points describing the eruption profile

            ! Recall that the column we are considering is within kground -> kPlumeTop
            ! The particular cell is z_cell_bot -> z_cell_top
            ! We need to loop over the z-steps of the profile find the weighted average
            ! of the z-steps that contribute to the cell in question.

            ! loop over the points describing the eruption profile
            do kk=1,e_prof_nzpoints(i)
              ! Find the bot/top altitude of this step of the eruption profile
              zbot_prof = e_prof_dz(i)*(kk-1)
              ztop_prof = e_prof_dz(i)*(kk)
              if(ztop_prof.lt.zground)then
                ! If line element is fully below zground
                frac = 0.0_ip
                if(e_prof_Volume(i,kk).gt.EPS_SMALL)then
                  do io=1,2;if(VB(io).le.verbosity_info)then
                    write(outlog(io),*)&
                      "WARNING: Eruption profile element is below vent elevation."
                    write(outlog(io),*)&
                      "         Profile line element will be disregarded and the profile renormalize."
                  endif;enddo
                endif
                cycle
              elseif((zbot_prof.le.z_cell_bot).and.&
                     (ztop_prof.ge.z_cell_bot).and.&
                     (ztop_prof.le.z_cell_top))then
                ! Add any fractional bit that straddles the bottom part of the z-cell
                frac = (ztop_prof-z_cell_bot)/dz_vec_pd(k)
              elseif((zbot_prof.ge.z_cell_bot).and.&
                     (ztop_prof.lt.z_cell_top))then
                ! Add any step fully within the z-cell
                frac = (ztop_prof-zbot_prof)/dz_vec_pd(k)
              elseif((zbot_prof.ge.z_cell_bot).and.&
                     (zbot_prof.le.z_cell_top).and.&
                     (ztop_prof.ge.z_cell_top))then
                ! Add any fractional bit that straddles the top part of the z-cell
                frac = (z_cell_top-zbot_prof)/dz_vec_pd(k)
              elseif(zbot_prof.lt.z_cell_bot.and.ztop_prof.gt.z_cell_top)then
                ! If the eruption profile step fully encompasses the z-cell, no frac needed
                frac = 1.0_ip
              else
                ! This profile step does not contribute to this z-cell
                frac = 0.0_ip
              endif
                NormSourceColumn(i,k) = NormSourceColumn(i,k) + &
                  e_prof_Volume(i,kk)*frac
            enddo
          else
            ! Source is none of suzuki,umbrella,umbrella_air,line,point,profile
            ! This is probably a non-tephra source or some custom source entered
            ! elsewhere. Set flux to zero for now.
            NormSourceColumn(i,k) = 0.0_ip
          endif
          ! Done with this k-cell; continue upwards to the top of the plume
        enddo ! nzmax+1

        tot = sum(NormSourceColumn(i,:))
        NormSourceColumn(i,:) = NormSourceColumn(i,:)/tot

      enddo ! neruptions

      end subroutine Calc_Normalized_SourceCol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  EruptivePulse_MassFluxRate
!
!  Called from: Ash3d.F90
!  Arguments:
!    none
!
!  This subroutine is called after the control file has been read and the
!  optional 'input_data_ResetParams' is called (MagmaDensity can be reset from
!  the control file).  For each eruptive pulse, the mass flux in kg/hr is
!  calculated from the specified volume (DRE) and duration.  The end time
!  of each eruptive pulse is also calculated here (in hours since BaseYear).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine EruptivePulse_MassFluxRate

      use global_param,  only : &
         KM3_2_M3

      use Tephra,        only : &
         MagmaDensity

      integer :: i

      ! Calculate mass flux and end times of each eruptive pulse
      do i=1,neruptions
             !mass flux in kg/hr
        if(SourceType.eq.'suzuki'      .or. &
           SourceType.eq.'point'       .or. &
           SourceType.eq.'line'        .or. &
           SourceType.eq.'umbrella'    .or. &
           SourceType.eq.'umbrella_air')then
          MassFluxRate(i)  = MagmaDensity * & ! kg/m3
                             e_Volume(i)  * & ! km3
                             KM3_2_M3     / & ! m3/km3
                             real(e_Duration(i),kind=ip)    ! hours  => kg/hr
          e_EndTime(i) = e_StartTime(i) + e_Duration(i)

          do io=1,2;if(VB(io).le.verbosity_info)then
            write(outlog(io),1023) MagmaDensity, e_Duration(i), &
                                   MassFluxRate(i), e_Volume(i)
          endif;enddo
1023      format('   Magma density (kg/m3) = ',f6.1,', Pulse Duration (hrs) = ',f6.3,/, &
                 '   Mass flux (kg/hr) = ',e12.4,', Pulse volume (km3 DRE)=',f8.4)
        elseif(SourceType.eq.'profile')then
          e_prof_MassFluxRate(i,1:e_prof_nzpoints(i)) =              &
                         MagmaDensity                          * & ! kg/m3
                         e_prof_Volume(i,1:e_prof_nzpoints(i)) * & ! km3
                         KM3_2_M3                              / & ! m3/km3
                         real(e_Duration(i),kind=ip)                             ! hours = kg/hr
          MassFluxRate(i) = sum(e_prof_MassFluxRate(i,1:e_prof_nzpoints(i)))
          e_EndTime(i) = e_StartTime(i) + e_Duration(i)
        else
          ! Custom source, initializing MassFluxRate and end time
          MassFluxRate(i)  = 0.0_ip
          e_EndTime(i) = 0.0_ip
        endif

      enddo

      e_EndTime_final = maxval(e_EndTime)  ! this marks the end of all eruptions

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),1024) e_EndTime_final-e_StartTime(1),&
                               sum(e_Volume(1:neruptions)), &
                               sum(e_Volume(1:neruptions))*MagmaDensity*KM3_2_M3*1.0e-9
      endif;enddo

1024  format('  Total Duration (hrs) = ',f6.3,/, &
             '  Total volume (km3 DRE) = ',f8.4,/,&
             '  Total mass (Tg) = ',f8.4)

      end subroutine EruptivePulse_MassFluxRate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  CheckEruptivePulses
!
!  Called from: Ash3d.F90
!  Arguments:
!    none
!
!  This subroutine is called in every step of the time integration.  Its job
!  is to test each of the eruptive pulse events to determine which are active
!  in the current dt.  In some cases, the dt of the time integration might overshoot
!  the end of an eruptive pulse, in which case, the mass flux is adjusted accordingly.
!  It might also involve multiple eruptive pulses.  If a pulse is active, the fraction
!  of the current time step (dt) that the pulse is active is logged to dt_pulse_frac()
!  and the logical flag, Source_in_dt, is set.  This flag is used in Ash3d.F90 to
!  actually insert the source term.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CheckEruptivePulses

      use time_data,     only : &
         time, dt

      real(kind=dp)    :: tstart, tend      ! start and end times of this time step
      logical          :: Pulse_contributes
      integer          :: i

      tstart = time
      tend   = time+dt

      Source_in_dt     = .false.
      dt_pulse_frac(:) = 0.0_dp

      if((SourceType.eq.'point')       .or. & ! profile is a branch below
         (SourceType.eq.'line')        .or. &
         (SourceType.eq.'suzuki')      .or. &
         (SourceType.eq.'umbrella')    .or. &
         (SourceType.eq.'umbrella_air').or. &
         (SourceType.eq.'profile'))then

        do i=ieruption,neruptions
          Pulse_contributes = .false.
          if((tstart.ge.e_StartTime(i)).and. & ! beginning of time step at or after pulse start
             (tstart.lt.e_EndTime(i)))then     ! beginning of time step is before same pulse ends
            ! This catches all pulses that touch the start of dt
            Pulse_contributes = .true.
          elseif((tend.gt.e_StartTime(i)).and. & ! end of time step at or after pulse start
                 (tend.le.e_EndTime(i)))then     ! end of time step is before same pulse ends
            ! This catches all pulses that touch the end of dt
            Pulse_contributes = .true.
            jeruption = i
            ! if we have found the eruptive pulse that touches the end of dt, then exit the
            ! do loop
            exit
          elseif(tstart.lt.e_StartTime(i).and.  &
                (tend.gt.e_EndTime(i)))then
            ! This catches all pulses that are wholely within dt
            ! If we are here, either the dt selected is huge (e.g. no winds) or the
            ! eruption duration is tiny.
            ! Issue a warning, but continue.
            do io=1,2;if(VB(io).le.verbosity_info)then
              write(outlog(io),1)i
            endif;enddo
            Pulse_contributes = .true.
          endif
          if(Pulse_contributes)then
            ! If any pulse contributes, update the global flag
            Source_in_dt = .true.
            ! calculate the sliver of time this pulse is active within dt
            dt_pulse_frac(i) = (min(tend,e_EndTime(i)) - &
                                max(tstart,e_StartTime(i)))/dt
          endif
        enddo
      else
        ! For all non-standard sources, assign the height given on the source
        ! line of the input file and assign a zero mass flux rate.
        Source_in_dt = .false.
        dt_pulse_frac(:) = 0.0_dp
        return
      endif

!     Format statements
 1     format(4x,'Warning.  Eruption ',i3,' is shorter than time steps dt.')

      end subroutine CheckEruptivePulses

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  TephraSourceNodes
!
!  Called from: Ash3d.F90
!  Arguments:
!    none
!
!  This subroutine is called for every step of the time integration in which
!  Source_in_dt = .true.. Its job is to calculate SourceNodeFlux(0:nzmax+1,1:nsmax)
!  which specifies the rate of tephra loading into each cell above the vent
!  in mass/vol/time (kg/km^3/hr) for each grainsize.  This is calculated for
!  all standard sources (point,line,suzuki,profile and umbrella/umbrella_air)
!  by using the normalized profile with the MassFluxRate for this dt (or the
!  weighted average is multiple eruptions are active in this dt.) to get
!  the total tephra flux rate as a function of height (TephraFluxRate(:)) in kg/hr.
!  This is then used with the GSD and cell volumes to get the concentration
!  inserted per hour for each cell and grainsize.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine TephraSourceNodes

      use global_param,    only : &
         EPS_SMALL

      use Tephra,        only : &
         Tephra_bin_mass,n_gs_max

      use mesh,          only : &
         nzmax,ivent,jvent,kappa_pd
 
      integer :: i,k
      real(kind=ip) :: z_cell_bot
      real(kind=ip) :: z_cell_top
      real(kind=ip) :: SumSourceNodeFlux      ! checking terms
      real(kind=ip) :: MassFluxRate_now

      SourceNodeFlux       = 0.0_ip           ! initialize SourceNodeFlux
      SumSourceNodeFlux    = 0.0_ip

      ! Loop over all the eruptive pulses active in this time step
      TephraFluxRate(:) = 0.0_ip
      MassFluxRate_now  = 0.0_ip
      do i=ieruption,jeruption
        TephraFluxRate(1:nzmax) = TephraFluxRate(1:nzmax)        + &
                                  MassFluxRate(i)                * &
                                  real(dt_pulse_frac(i),kind=ip) * &
                                  NormSourceColumn(ieruption,1:nzmax)
        MassFluxRate_now = MassFluxRate_now + MassFluxRate(i) * &
                           real(dt_pulse_frac(i),kind=ip)
      enddo

        ! Now that we have the TephraFluxRate as a function of k, convert it to mass
        ! over the grainsmax bins stored in SourceNodeFlux (kg/km3/hr)
        ! SumSourceNodeFlux is used for check the sum of all source nodes
        ! against the total MassFluxRate (i.e. should equal 1.0)
      do k=1,nzmax
        SourceNodeFlux(k,1:n_gs_max) =      & ! final units are kg/km3/hr
              Tephra_bin_mass(1:n_gs_max) * & ! fraction of total in bin
              TephraFluxRate(k)           / & ! kg/hr
              kappa_pd(ivent,jvent,k)         ! km3
        SumSourceNodeFlux = &
              SumSourceNodeFlux +        &         ! dimensionless
              sum(SourceNodeFlux(k,1:n_gs_max) * & ! kg/km3 hr
              kappa_pd(ivent,jvent,k)) / &         ! km3
              MassFluxRate_now
      enddo
      ! Make sure the sum of the fluxes in all the cells equals the total flux
      if (abs(SumSourceNodeFlux-1.0_ip).gt.EPS_SMALL) then
         do io=1,2;if(VB(io).le.verbosity_error)then
           write(errlog(io) ,2) SumSourceNodeFlux-1.0_ip
           write(errlog(io),*)"SourceType          = ",SourceType
           write(errlog(io),*)"z_cell_bot          = ",z_cell_bot
           write(errlog(io),*)"z_cell_top          = ",z_cell_top
           write(errlog(io),*)"MassFluxRate_now    = ",MassFluxRate_now
           write(errlog(io),*)"n_gs_max            = ",n_gs_max
           write(errlog(io),*)"SourceNodeFlux(1:nz)=",real(SourceNodeFlux(:,1),kind=4)
         endif;enddo
        stop 1
      endif

      ! Now that we have prepared the contributions, update the eruption index
      ! for the start of the next step
      ieruption = jeruption

!     Format statements
2     format(4x,'Source Node Flux does not agree with calculations.',/, &
              4x,'(Sum(SourceNodeFlux)/MassFluxRate)-1=',e12.5,/, &
              4x,'Program stopped')

      end subroutine TephraSourceNodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  SourceVolInc(dt)
!
!  Called from: Ash3d.F90
!  Arguments:
!    dt = time step in hours
!
!  This function calculates the total tephra volume inserted in this time step.
!  It is used only for mass-conservation error-checking.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function SourceVolInc(dt)

      use global_param,  only : &
         KM3_2_M3

      use Tephra,        only : &
         n_gs_max,MagmaDensity

      use mesh,          only : &
         nzmax,kappa_pd,ivent,jvent

      real(kind=ip) :: SourceVolInc
      real(kind=dp) :: dt

      real(kind=ip) :: tmp
      integer :: k,isize

      tmp = 0.0_ip

      do k=1,nzmax+1
        do isize=1,n_gs_max
          tmp = tmp                             + & ! final units is km3
                real(dt,kind=ip)                * & ! hr
                SourceNodeFlux(k,isize)         * & ! kg/km3 hr
                kappa_pd(ivent,jvent,k)         / & ! km3
                MagmaDensity                    / & ! kg/m3
                KM3_2_M3                            ! m3/km3
        enddo
      enddo

      SourceVolInc = tmp

      return

      end function SourceVolInc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module Source
!##############################################################################
