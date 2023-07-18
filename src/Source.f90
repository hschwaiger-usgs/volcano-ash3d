!      subroutine Allocate_Source_eruption
!      subroutine Allocate_Source_grid
!      subroutine Allocate_Source_time
!      subroutine Deallocate_Source
!      subroutine TephraSourceNodes
!      subroutine MassFluxCalculator
!      function SourceVolInc

      module Source

      use precis_param

      use io_units

      implicit none

        ! Set everything to private by default
      private

        ! Publicly available subroutines/functions
      public Allocate_Source_eruption, Allocate_Source_grid,Allocate_Source_time,&
             Deallocate_Source,MassFluxCalculator,TephraSourceNodes,SourceVolInc

        ! Publicly available variables


      integer, parameter :: MAXCustSrc = 10  ! The maximum number of custom
                                               !source types that we will check for

      real(kind=ip),public :: x_volcano, y_volcano      ! x & y points of volcano (m)
      real(kind=ip),public :: lat_volcano, lon_volcano  !position of volcano in lat/lon
      real(kind=ip),public :: z_volcano                 ! vent elevation

      integer,public :: neruptions                  ! number of eruptions or eruptive pulses
      character(len=12),public :: SourceType          !may be 'point', 'line', or 'Suzuki' 
      real(kind=ip),public :: Suzuki_A
      logical,public       :: IsCustom_SourceType = .false.
      character(len=30),dimension(MAXCustSrc) :: SourceType_Custom = ""

#ifdef USEPOINTERS
      real(kind=ip), dimension(:,:)  ,pointer,public :: SourceColumn
      real(kind=ip), dimension(:  )  ,pointer,public :: TephraFluxRate  => null()
      real(kind=ip), dimension(:,:)  ,pointer,public :: SourceNodeFlux      => null()
      real(kind=ip), dimension(:,:,:),pointer :: SourceNodeFlux_Area => null()
#else
      real(kind=ip), dimension(:)    ,allocatable,public     :: TephraFluxRate
      real(kind=ip), dimension(:,:)  ,allocatable,public     :: SourceNodeFlux
      real(kind=ip), dimension(:,:,:),allocatable     :: SourceNodeFlux_Area
#endif

        !The following arrays are used by MassFluxCalculator
      real(kind=ip),public :: MassFluxRate_now
      real(kind=ip),public :: Height_now
      integer :: ieruption !eruption we're currently on

        !The following arrays are of length neruptions
      real(kind=ip), dimension(:)        ,allocatable,public :: e_PlumeHeight
      real(kind=ip), dimension(:)        ,allocatable,public :: e_Volume
      real(kind=ip), dimension(:)        ,allocatable,public :: e_Duration
      real(kind=ip), dimension(:)        ,allocatable,public :: e_StartTime
      real(kind=ip), dimension(:)        ,allocatable,public :: e_EndTime
      real(kind=ip), dimension(:)        ,allocatable,public :: MassFlux
      real(kind=ip), dimension(:)        ,allocatable,public :: e_prof_dz
      integer      , dimension(:)        ,allocatable,public :: e_prof_zpoints
      real(kind=ip), dimension(:,:)      ,allocatable,public :: e_prof_Volume
      real(kind=ip), dimension(:,:)      ,allocatable,public :: e_prof_MassFlux
      real(kind=ip),public :: e_EndTime_final

      real(kind=ip),public :: ESP_height        = 0.0_ip
      real(kind=ip),public :: ESP_duration      = 0.0_ip
      real(kind=ip),public :: ESP_MassFluxRate  = 0.0_ip
      real(kind=ip),public :: ESP_Vol           = 0.0_ip
      real(kind=ip),public :: ESP_massfracfine  = 0.0_ip

      contains

!******************************************************************************

      subroutine Allocate_Source_eruption

      allocate (e_StartTime(neruptions));            e_StartTime   = 0.0_ip
      allocate (e_Duration(neruptions));             e_Duration    = 0.0_ip
      allocate (e_PlumeHeight(neruptions));          e_PlumeHeight = 0.0_ip
      allocate (e_Volume(neruptions));               e_Volume      = 0.0_ip
      allocate (MassFlux(neruptions));               MassFlux      = 0.0_ip
      allocate (e_EndTime(neruptions));              e_EndTime     = 0.0_ip

      if(SourceType.eq.'profile')then
        allocate (e_prof_dz(neruptions));             e_prof_dz      = 0.0_ip
        allocate (e_prof_zpoints(neruptions));        e_prof_zpoints = 0
          ! for profiles, assume 50 points
        allocate (e_prof_Volume(neruptions,50));      e_prof_Volume   = 0.0_ip
        allocate (e_prof_MassFlux(neruptions,50));    e_prof_MassFlux = 0.0_ip
      endif

      end subroutine Allocate_Source_eruption

!******************************************************************************

      subroutine Allocate_Source_grid(nx,ny,nz)

      use mesh,          only : &
         nsmax

      integer,intent(in) :: nx,ny,nz

      allocate(SourceNodeFlux(0:nz+1,1:nsmax));      SourceNodeFlux = 0.0_ip
      allocate(TephraFluxRate(nz));                  TephraFluxRate = 0.0_ip

      end subroutine Allocate_Source_grid

!******************************************************************************

      subroutine Allocate_Source_time

      MassFluxRate_now = 0.0_ip
      Height_now = 0.0_ip

      ieruption = 1 ! Initialize eruption to the starting eruption

      end subroutine Allocate_Source_time


!******************************************************************************

      subroutine Deallocate_Source

      deallocate (e_StartTime)
      deallocate (e_Duration)
      deallocate (e_PlumeHeight)
      deallocate (e_Volume)
      deallocate (MassFlux)
      deallocate (e_EndTime)

      if(SourceType.eq.'profile')then
        deallocate (e_prof_dz)
        deallocate (e_prof_zpoints)
        deallocate (e_prof_Volume)
        deallocate (e_prof_MassFlux)
      endif

      deallocate(SourceNodeFlux)

      end subroutine Deallocate_Source
!******************************************************************************

      subroutine TephraSourceNodes

!     subroutine that finds the source nodes and assigns an ash concentration to them

      use global_param,    only : &
         EPS_SMALL

      use Tephra,        only : &
         Tephra_bin_mass,n_gs_max

      use mesh,          only : &
         IsLatLon,nzmax,dx,dy,dz_vec_pd,ivent,jvent,z_lb_pd,z_cc_pd,kappa_pd
 
      integer :: k
      real(kind=ip) :: Suzuki_k     ! k factor in the Suzuki equation (see Hurst's Ashfall manual)
      real(kind=ip) :: z_cell_bot, z_cell_top
      real(kind=ip) :: SumSourceNodeFlux      ! checking terms
      real(kind=ip) :: zground, PlumeHeight_above_ground
      integer       :: kground
      real(kind=ip) :: ez
      integer       :: kk
      integer       :: kPlumeTop

      SourceNodeFlux       = 0.0_ip           !initialize SourceNodeFlux
      SumSourceNodeFlux    = 0.0_ip

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
        if(Height_now.ge.z_lb_pd(k  ).and. &
           Height_now.lt.z_lb_pd(k+1))then
          kPlumeTop = k
        endif
      enddo
      zground = z_cc_pd(kground) - 0.5_ip*dz_vec_pd(kground)

      if ((SourceType.eq.'suzuki')      .or. &
          (SourceType.eq.'umbrella')    .or. &
          (SourceType.eq.'umbrella_air')) then
        Suzuki_k = Suzuki_A/((Height_now-zground)* &
                  ((1.0_ip/Suzuki_A)-((Suzuki_A+1.0_ip)/Suzuki_A)* &
                  exp(-Suzuki_A)))
      endif

      ! FIND THE FRACTION OF THE ERUPTED MASS RATE (kg/hr) THAT LIES
      ! WITHIN EACH HEIGHT INTERVAL DZ
      TephraFluxRate(1:nzmax) = 0.0_ip ! kg/hr of total tephra as a function of z
      do k=kground,kPlumeTop
        ! height at the top of this cell (or top of plume)
        z_cell_top  = min(z_cc_pd(k)+0.5_ip*dz_vec_pd(k),Height_now)
        ! height at the bottom
        z_cell_bot = z_cc_pd(k)-0.5_ip*dz_vec_pd(k)

        ! First get the TephraFluxRate in kg/hr (total mass of tephra inserted)
        !TephraFluxRate = 0.0_ip
        PlumeHeight_above_ground = Height_now-zground
        if ((SourceType.eq.'suzuki')      .or. &
            (SourceType.eq.'umbrella')    .or. &
            (SourceType.eq.'umbrella_air')) then
          !For Suzuki plumes and umbrella clouds
          ! HANS: double-check units and cite equation from User's Guide
          ! The equation below calculates the tephra volume flux 
          ! (m3 DRE/s) in the height interval.  
          ! It uses an equation obtained by integrating the 
          ! Suzuki equation given in Hurst.
          TephraFluxRate(k) = (MassFluxRate_now*Suzuki_k*                     &
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
          !for line sources, the fractional Tephra Flux Rate into the cell
          ! at z is just the height of the cell over the length of the line
          TephraFluxRate(k) = MassFluxRate_now * &
                              (z_cell_top-z_cell_bot) / PlumeHeight_above_ground
        elseif (SourceType.eq.'point') then
          !for point sources, put all the MassFluxRate into the cell that contains
          ! the point
          if ((z_cell_top.ge.Height_now).and.    &
              (z_cell_bot.lt.Height_now)) then
            TephraFluxRate(k) = MassFluxRate_now
          else
            TephraFluxRate(k) = 0.0_ip
          endif
        elseif (SourceType.eq.'profile') then
          ! loop over the points describing the eruption profile
          do kk=1,e_prof_zpoints(ieruption)
            ! Find the top node in the eruption profile
            ez = e_PlumeHeight(ieruption) - &
                    e_prof_dz(ieruption)*e_prof_zpoints(ieruption) + &
                    (kk-1)*e_prof_dz(ieruption)
            if ((z_cell_top.ge.ez).and. & 
                (z_cell_bot.lt.ez)) then
              ! This assumes that the timestep is fully within the
              ! eruption
              TephraFluxRate(k) = TephraFluxRate(k) + e_prof_MassFlux(ieruption,kk)
            else
              TephraFluxRate(k) = 0.0_ip
            endif
          enddo
        else
          ! Source is none of suzuki,umbrella,umbrella_air,line,point,profile
          ! This is probably a non-tephra source or some custom source entered
          ! elsewhere. Set flux to zero for now.
          TephraFluxRate(k) = 0.0_ip
        endif
        ! Done with this k-cell; continue upwards to the top of the plume
      enddo

        ! Now that we have the TephraFluxRate as a function of k, convert it to mass
        ! over the grainsmax bins stored in SourceNodeFlux (kg/km3/hr)
        ! SumSourceNodeFlux is used for check the sum of all source nodes
        ! against the total MassFlux (i.e. should equal 1.0)
      do k=kground,kPlumeTop
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
           write(errlog(io),*)"SourceType = ",SourceType
           write(errlog(io),*)"Height_now =",Height_now
           write(errlog(io),*)"z_cell_top,z_cell_bot = ",z_cell_top,z_cell_bot
           write(errlog(io),*)"MassFluxRate_now = ",MassFluxRate_now
           write(errlog(io),*)"n_gs_max = ",n_gs_max
           write(errlog(io),*)"SourceNodeFlux(1:nz)=",SourceNodeFlux(:,1)
         endif;enddo
        stop 1
      endif

!     Format statements
2     format(4x,'Source Node Flux does not agree with calculations.',/, &
              4x,'(Sum(SourceNodeFlux)/MassFluxRate)-1=',e12.5,/, &
              4x,'Program stopped')

      end subroutine TephraSourceNodes

!******************************************************************************

      subroutine MassFluxCalculator

!     function that calculates the mass flux as a function of time

      use time_data,     only : &
         time, dt

      real(kind=ip)    :: tstart, tend    !start and end times of this time step

      tstart = time
      tend   = time+dt

      MassFluxRate_now = 0.0_ip  ! kg/hr
      Height_now       = 0.0_ip  ! km

      !if (itime.eq.0) return     !exit the subroutine if itime=0

      if((SourceType.eq.'point')       .or. & ! profile is a branch below
         (SourceType.eq.'line')        .or. &
         (SourceType.eq.'suzuki')      .or. &
         (SourceType.eq.'umbrella')    .or. &
         (SourceType.eq.'umbrella_air'))then

        !COMPARE THE TIME WITH THE START & END TIMES OF THE DIFFERENT ERUPTIONS
        !exit the subroutine if we're past the last eruption
        if (tstart.gt.e_EndTime(neruptions)) return
  
        if (ieruption.gt.neruptions) then
          ! If the time slice starts after the last eruption ends
          ! do nothing (MassFluxRate & Height are already 0)
          return
        elseif (tstart>=e_StartTime(ieruption)) then
          ! If the time slice starts after the eruption starts . . .
          if (tend.le.e_EndTime(ieruption))then
            ! . . . and ends before the eruption ends, MassFluxRate=avg.
            MassFluxRate_now = MassFlux(ieruption)
            Height_now       = e_PlumeHeight(ieruption)
          else
            ! . . . and ends after the eruption ends, do some more checking.
            if (ieruption.lt.neruptions) then
              !There are multiple eruptions. See if we need to construct
              !time averages
              if (tend.gt.e_EndTime(ieruption+1)) then
                ! If it ends AFTER the END of the NEXT eruption, stop the
                ! program.
                ! This means that the chosen dt is too large.
                do io=1,2;if(VB(io).le.verbosity_error)then
                  write(errlog(io),1)  ieruption+1
                endif;enddo
                stop 1
              elseif (tend.gt.e_StartTime(ieruption+1)) then
                ! If it ends AFTER the START of the NEXT eruption,
                ! construct weighted average of the two eruptions
                MassFluxRate_now    = MassFlux(ieruption)*(e_EndTime(ieruption)-tstart)  /dt + &
                                      MassFlux(ieruption+1)*(tend-e_StartTime(ieruption+1))/dt
                Height_now          = (e_PlumeHeight(ieruption) + e_PlumeHeight(ieruption+1))/2.0_ip
                ieruption = ieruption+1
              else
                ! If it ends BEFORE the START of the NEXT eruption,
                ! distribute the fraction of the time-step that is
                ! erupting over the whole dt
                MassFluxRate_now = MassFlux(ieruption)*(e_EndTime(ieruption)-tstart)/dt
                Height_now       = e_PlumeHeight(ieruption)
                ieruption = ieruption+1
              endif
            else
              !If this is either the only or the last eruption, distribute
              !the fraction of the time-step that is erupting over the
              !whole dt
  
              MassFluxRate_now = MassFlux(ieruption)*(e_EndTime(ieruption)-tstart)/dt
              Height_now       = e_PlumeHeight(ieruption)
              ieruption = ieruption+1
            endif
          endif  ! end of check on 
        else
          ! The time slice starts before the eruption starts . . .
          if (tend>e_StartTime(ieruption)) then
            if (tend.le.e_EndTime(ieruption)) then
              ! . . . and ends during the eruption, interpolate
  
              MassFluxRate_now = MassFlux(ieruption)*(tend-e_StartTime(ieruption))/dt
              Height_now       = e_PlumeHeight(ieruption)
            else
              !If it ends after the erupt ends, do some more checking
              if (ieruption.lt.neruptions) then
                !If it ends before the next eruption starts, interpolate
                if (tend.lt.e_StartTime(ieruption+1)) then
                  MassFluxRate_now = MassFlux(ieruption)*(e_Duration(ieruption)/dt)
                  Height_now       = e_PlumeHeight(ieruption)
                  ieruption=ieruption+1
                  !If it ends during the next eruption do another interpolation
                elseif (tend.le.e_EndTime(ieruption+1)) then
                  MassFluxRate_now = MassFlux(ieruption)*((e_EndTime(ieruption)-tstart)/dt) + &
                                     MassFlux(ieruption+1)*((tend-e_StartTime(ieruption))/dt)
                  Height_now       = (e_PlumeHeight(ieruption)+e_PlumeHeight(ieruption+1))/2.0_ip

                  ieruption=ieruption+1
                  !If it ends after the next eruption ends, stop the program
                else
                  do io=1,2;if(VB(io).le.verbosity_error)then
                    write(errlog(io),1)  ieruption+1
                  endif;enddo
                  stop 1
                endif
                !If if ends after the last eruption ends, interpolate and exit
              else
                Height_now       = e_PlumeHeight(ieruption)
                MassFluxRate_now = MassFlux(ieruption)*(e_Duration(ieruption)/dt)
                ieruption = ieruption + 1
              endif
            endif
          else                                 ! . . . and ends before the eruption starts, move on
            continue
          endif
        endif

        return

      elseif(SourceType.eq.'profile')then
        Height_now       = e_PlumeHeight(ieruption)
        MassFluxRate_now = sum(e_prof_MassFlux(ieruption,:))
        return

      else
        ! For all non-standard sources, assign the height given on the source
        ! line of the input file and assign a zero mass flux rate.
        Height_now       = e_PlumeHeight(ieruption)
        MassFluxRate_now = 0.0_ip
        return
      endif

!     Format statements
 1     format(4x,'Warning.  Eruption ',i3,' is shorter than time steps dt.  Program stopped.')

      end subroutine MassFluxCalculator

!******************************************************************************

      function SourceVolInc(dt)

      use global_param,  only : &
         KM3_2_M3

      use Tephra,        only : &
         n_gs_max,MagmaDensity

      use mesh,          only : &
         nzmax,kappa_pd,ivent,jvent

      real(kind=ip) :: SourceVolInc
      real(kind=ip) :: dt

      real(kind=ip) :: tmp
      integer :: k,isize

      tmp = 0.0_ip

      do k=1,nzmax+1
        do isize=1,n_gs_max
          tmp = tmp                             + & ! final units is km3
                dt                              * & ! hr
                SourceNodeFlux(k,isize)         * & ! kg/km3 hr
                kappa_pd(ivent,jvent,k)         / & ! km3
                MagmaDensity                    / & ! kg/m3
                KM3_2_M3                            ! m3/km3
        enddo
      enddo

      SourceVolInc = tmp

      return

      end function SourceVolInc

!******************************************************************************

      end module Source
