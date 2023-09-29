!##############################################################################
!
!  Source_Umbrella module
!
!      subroutine Allocate_Source_Umbrella
!      subroutine Deallocate_Source_Umbrella
!      subroutine umbrella_winds
!      subroutine TephraSourceNodes_Umbrella
!      function SourceVolInc_Umbrella
!      function AvgCon_Umbrella
!
!##############################################################################

      module Source_Umbrella

      use precis_param

      use io_units

      implicit none

        ! Set everything to private by default
      private

        ! Publicly available subroutines/functions
      public Allocate_Source_Umbrella,  &
             Deallocate_Source_Umbrella,&
             umbrella_winds,            &
             TephraSourceNodes_Umbrella,&
             SourceVolInc_Umbrella,     &
             AvgCon_Umbrella

      !components of the wind field used for umbrella clouds
#ifdef USEPOINTERS
      real(kind=ip),dimension(:,:,:,:),pointer,public :: SourceNodeFlux_Umbrella =>null()
      real(kind=ip),dimension(:,:,:)  ,pointer,public :: uvx_pd =>null() ! u (E) component of wind
      real(kind=ip),dimension(:,:,:)  ,pointer,public :: uvy_pd =>null() ! v (N) component of wind
#else
      real(kind=ip),dimension(:,:,:,:),allocatable,public :: SourceNodeFlux_Umbrella
      real(kind=ip),dimension(:,:,:)  ,allocatable,public :: uvx_pd ! u (E) component of wind
      real(kind=ip),dimension(:,:,:)  ,allocatable,public :: uvy_pd ! v (N) component of wind
#endif

      !The following are used by Mesointerpolator for umbrella clouds
      integer,public :: ibase     !z index of lowest node in the umbrella cloud
      integer,public :: itop     !z index of highest node in the umbrella cloud

      ! These are only public since we want to write these to the ouput file
      integer      ,public :: VelMod_umb         = 1       ! Velocity model to use (1=default)
      real(kind=ip),public :: k_entrainment_umb  = 0.1_ip  ! entrainment coefficient
      real(kind=ip),public :: lambda_umb         = 0.2_ip  ! umbrella cloud shape factor
      real(kind=ip),public :: N_BV_umb           = 0.02_ip ! Brunt-Vaisala frequency, 1/s
      real(kind=ip),public :: SuzK_umb           = 12.0_ip ! Suzuki parameter used for umb clouds

       !width & height of source nodes in km
      real(kind=ip) :: SourceNodeWidth_km
      real(kind=ip) :: SourceNodeHeight_km

      real(kind=ip),dimension(:,:,:)  ,allocatable :: AvgStenc_Umbrella
      real(kind=ip),dimension(:)      ,allocatable :: ScaleFac_Umbrella

      contains
      !------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Allocate_Source_Umbrella(nx,ny,nz)
!
!  Called from: Ash3d.F90
!  Arguments:
!    nx         = x length of computational grid (needed for supplemental wind)
!    ny         = y length of computational grid (needed for supplemental wind)
!    nz         = z length of computational grid (needed for modified source column)
!
!  This subroutine allocates the arrays needed for the umbrella source and
!  pre-calculates some geometry and normalizing terms.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Allocate_Source_Umbrella(nx,ny,nz)

      use mesh,          only : &
         z_vec_init,kappa_pd,ivent,jvent,de_km,dn_km

      use Source,        only : &
         e_PlumeHeight

      use Tephra,        only : &
         n_gs_max

      integer,intent(in) :: nx
      integer,intent(in) :: ny
      integer,intent(in) :: nz

      integer :: i,j,k,ii,jj
      real(kind=ip) :: tot_vol

      allocate(SourceNodeFlux_Umbrella(3,3,nz+1,n_gs_max)); SourceNodeFlux_Umbrella=0.0_ip
      allocate(ScaleFac_Umbrella(nz));                  ScaleFac_Umbrella=0.0_ip
      allocate(AvgStenc_Umbrella(3,3,nz));              AvgStenc_Umbrella=0.0_ip

      SourceNodeWidth_km  = (3.0_ip/2.0_ip)*de_km
      SourceNodeHeight_km = (3.0_ip/2.0_ip)*dn_km

      ! Get the weights for the averaging stencil for each z-level
      ! The averaging stencil is used to smooth the concentration in the 3x3
      ! zone around the vent to facilitate smooth radial spreading (otherwise
      ! it will be lumpy).
      ! The scale factor is used for converting mass inserted as a point above the
      ! vent to mass inserted in the same 3x3 patch.  These factors only need to
      ! be calculated once at the beginning of the run.
      do k = 1,nz
        tot_vol = 0.0_ip
        do i=1,3
          ii=ivent-2+i
          do j=1,3
            jj=jvent-2+j
            AvgStenc_Umbrella(i,j,k)=kappa_pd(ii,jj,k)
            tot_vol = tot_vol + kappa_pd(ii,jj,k)
          enddo
        enddo
        ScaleFac_Umbrella(k)     = kappa_pd(ivent,jvent,k)/tot_vol
        AvgStenc_Umbrella(:,:,k) = AvgStenc_Umbrella(:,:,k)/tot_vol
      enddo

      itop  = 0
      ibase = 0
      do k = 1,nz
        if(z_vec_init(k)  .ge.e_PlumeHeight(1).and.&
           z_vec_init(k-1).lt.e_PlumeHeight(1))then
          itop = k
        endif
        ! The Suzuki distribution for umbrella clouds is hardwired to k=12
        ! The 0.75 for the plume bottom just ensures that the significant part
        ! of the Suzuki mass loading profile will be within the ibase to itop nodes
        if(z_vec_init(k)  .ge.0.75_ip*e_PlumeHeight(1).and. &
           z_vec_init(k-1).lt.0.75_ip*e_PlumeHeight(1))then
          ibase = k
        endif
      enddo
      allocate(uvx_pd(-1:nx+2,-1:ny+2,ibase:itop));     uvx_pd = 0.0_ip
      allocate(uvy_pd(-1:nx+2,-1:ny+2,ibase:itop));     uvy_pd = 0.0_ip

      end subroutine Allocate_Source_Umbrella

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Deallocate_Source_Umbrella
!
!  Called from: dealloc_arrays
!  Arguments:
!    none
!
!  This subroutine deallocates the arrays allocated in Allocate_Source_Umbrella
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Deallocate_Source_Umbrella

        if(allocated(AvgStenc_Umbrella))       deallocate(AvgStenc_Umbrella)
        if(allocated(ScaleFac_Umbrella))       deallocate(ScaleFac_Umbrella)
        if(allocated(SourceNodeFlux_Umbrella)) deallocate(SourceNodeFlux_Umbrella)
        if(allocated(uvx_pd))                  deallocate(uvx_pd)
        if(allocated(uvx_pd))                  deallocate(uvy_pd)

      end subroutine Deallocate_Source_Umbrella

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  umbrella_winds(first_time)
!
!  Called from: MesoInterpolater
!  Arguments:
!    first_time = logical
!
!  This subroutine calculates the radial wind associated with the umbrella
!  spreading.  Equations implemented are described in:
!       Mastin, L.G. and Van Eaton, A.R., Comparing Simulations of Umbrella-Cloud Growth
!         and Ash Transport with Observations from Pinatubo, Kelud, and Calbuco Volcanoes
!         Atmosphere, 2020, 11, 1038, doi:10.3390/atmos11101038

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine umbrella_winds(first_time)

      use global_param,  only : &
         DEG2KMLAT,DEG2KMLON,DEG2RAD,KM_2_M,PI,HR_2_S,MPS_2_KMPHR,EPS_SMALL

      use mesh,          only : &
         nxmax,nymax,dn_km,de_km,IsLatLon,ivent,jvent,&
         lat_cc_pd,lon_cc_pd,de_km,dn_km

      use time_data,     only : &
         time,Simtime_in_hours,dt

      use Source,        only : &
         lat_volcano,lon_volcano,MassFluxRate,e_EndTime, SourceType

      use Tephra,        only : &
         n_gs_max

      logical, intent(in)  :: first_time

      real(kind=ip) :: avg_lat          ! avg latitude between vent & point
      real(kind=ip) :: C_Costa          ! C constant used in Costa et al., 2013
      real(kind=ip) :: cloud_radius     ! cloud radius, km
      real(kind=ip) :: cloudrad_raw     ! cloud radius, km, uncorrected
      real(kind=ip) :: edge_speed       ! expansion rate of cloud edge, m/s
      real(kind=dp) :: etime_s          ! time since eruption start, seconds
      real(kind=ip) :: ew_km,ns_km      ! distances between vent & point
      real(kind=ip) :: qnow             ! volume flow rate into umbrella cloud, m3/s
      real(kind=ip) :: radnow           ! radial distance from cloud center, km
      real(kind=ip) :: thetanow         ! angle of point CW from east
      real(kind=ip) :: windspeedhere    ! windspeed at this node
      real(kind=ip) :: rexp             ! exponent in radial vel term
      integer       :: ii,jj,iz         ! counters
      integer       :: ew_nodes,ns_nodes! radius of clouds in nodes
      integer       :: west_node,east_node
      integer       :: south_node,north_node
      real(kind=ip) :: MassFluxRateMKS_now ! current mass flux rate, kg/s

      uvx_pd(-1:nxmax+2,-1:nymax+2,ibase:itop) = 0.0_ip               !set umbrella winds to zero
      uvy_pd(-1:nxmax+2,-1:nymax+2,ibase:itop) = 0.0_ip    
 
      if(.not.IsLatLon)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*) 'Error: umbrella_winds is not yet set up to handle'
          write(errlog(io),*) 'projected coordinates.'
        endif;enddo
        stop 1
      endif

      ! Convert MassFluxRate from kg/hr to a local variable in kg/s
      MassFluxRateMKS_now = MassFluxRate(1)/HR_2_S        !mass flux rate, kg/s

      !If there is only one size class, assume it's an airborne run and
      !multiply the mass flux by 20
      if ((SourceType.eq.'umbrella_air').and.(n_gs_max.eq.1)) then
        MassFluxRateMKS_now = 20.0_ip*MassFluxRateMKS_now
      end if

      !set value of C based on latitude
      if (abs(lat_volcano).lt.23.0_ip) then
          ! m3 kg^(-3/4) s^(-7/8) for tropical eruptions
        C_Costa = 0.43e3_ip
      else
          ! m3 kg^(-3/4) s^(-7/8) for non-tropical eruptions
        C_Costa = 0.87e3_ip 
      endif

      ! Here is Eq. 2 of Mastin and Van Eaton, 2020 (m3/s)
      qnow  = C_Costa*sqrt(k_entrainment_umb)*MassFluxRateMKS_now**(3.0_ip/4.0_ip) / &
              N_BV_umb**(5.0_ip/4.0_ip)

      if(time.lt.EPS_SMALL) then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*) 
          write(outlog(io),*) 'in Umbrella_winds'
          write(outlog(io),*) '  massfluxnow (kg/s) = ',real(MassFluxRateMKS_now,kind=sp)
          write(outlog(io),*) '             C_Costa = ',real(C_Costa,kind=sp)
          write(outlog(io),*) '                N_BV = ',real(N_BV_umb,kind=sp)
          write(outlog(io),*) '       k_entrainment = ',real(k_entrainment_umb,kind=sp)
          write(outlog(io),*) '            Q (m3/s) = ',real(qnow,kind=sp)
          write(outlog(io),*)
          !If we're doing an airborne run and using only 1 grain size,
          !multiply the MER by 20 to make sure we're getting the right umbrella
          !growth rate
          if((SourceType.eq.'umbrella_air').and.(n_gs_max.eq.1))then
            write(outlog(io),*) 'n_gs_max=1, so we are assuming an airborne run'
            write(outlog(io),*) 'massflux has been multiplied by 20'
          endif
        endif;enddo
        if(VelMod_umb.ne.1.and.VelMod_umb.ne.2)then
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)" ERROR: Umb.Vel.Model unknown.  Should be 1 or 2."
            write(errlog(io),*)" VelMod_umb  = ",VelMod_umb
          endif;enddo
          stop 1
        endif
      endif ! time.eq.0.0_ip

      !convert from  hours to seconds
!      if (itime.gt.0) then
!        if (time.gt.0.0_ip) then
        if(.not.first_time)then
          ! For the first time step, time=0 -> cloudrad=0 and uR=Inf
          ! so just use dt to get values at the end of step 1
          etime_s      = max(time,dt)*real(HR_2_S,kind=ip)
        else
          etime_s      = min(Simtime_in_hours,e_EndTime(1))*real(HR_2_S,kind=ip)
          !return                      !return to Mesointerpolator if time=0
        endif
!      else
         ! If this is the first call to mesointerpolator before the beginning 
         ! of the simulation, the call is made simply to find the first value of 
         ! dt so that ntmax can be assigned using ntmax=int(SimTime_in_hours/dt).  
         ! If the initial value of dt underestimates the average time step used in 
         ! the simulation, it will under-allocate the array size. The radial wind 
         ! speeds in adjacent nodes increase with time during the eruption,
         ! meaning that dt should decrease.  Thus we want to estimate dt using radial 
         ! wind speeds at the end of the eruption, as below.  These wind speeds are 
         ! not actually used in the calculation of advecting ash.
!        etime_s      = min(3600.0_ip*Simtime_in_hours,3600.0_ip*e_EndTime(1))
!      endif

      ! Here is Eq. 1 of Mastin and Van Eaton, 2020
      !cloud radius, km
      cloudrad_raw = (3.0_ip*lambda_umb*N_BV_umb*qnow/(2.0_ip*PI))**(1.0_ip/3.0_ip) * &
                     real(etime_s,kind=ip)**(2.0_ip/3.0_ip) / KM_2_M
      !Make sure cloud radius extends beyond the source nodes
      !cloud_radius = max(cloudrad_raw,max(SourceNodeWidth_km,SourceNodeHeight_km))
      cloud_radius = cloudrad_raw            !for debugging

      ! This is the time derivitive of Eq. 1 of Mastin and Van Eaton, 2020
      !cloud expansion rate, m/s
      edge_speed   = (2.0_ip/3.0_ip)*(3.0_ip*lambda_umb*N_BV_umb*qnow/(2.0_ip*PI))**(1.0_ip/3.0_ip) * &
                    real(etime_s,kind=ip)**(-1.0_ip/3.0_ip)  

      if (cloud_radius.le.max(SourceNodeWidth_km,SourceNodeHeight_km)) then
        return
      else
        if (IsLatLon) then
          !calculate cloud size in nodes in x and y
          ew_nodes = int(cloud_radius/de_km)+1
          ns_nodes = int(cloud_radius/dn_km)+1
          west_node = max(1,ivent-ew_nodes)
          east_node = min(nxmax,ivent+ew_nodes)
          south_node = max(1,jvent-ns_nodes)
          north_node = min(nymax,jvent+ns_nodes)

          ! Find distance to volcano from  each node center
          do ii=west_node,east_node
            do jj=south_node,north_node
              !skip the source nodes
              if ((ii.eq.ivent).and.(jj.eq.jvent)) cycle
              !Calculate radial distance to the vent (km)
              avg_lat=(lat_cc_pd(jj)+lat_volcano)/2.
              !These formulas ensure that winds are added if the cell center is
              !within the umbrella radius
              ns_km  =(lat_cc_pd(jj)-lat_volcano)*DEG2KMLAT
              ew_km  =(lon_cc_pd(ii)-lon_volcano)* &
                        cos(avg_lat*DEG2RAD)*DEG2KMLON
              !These formulas ensure that winds are added if the all of the cell is
              !within the umbrella radius
              !ns_km  =(abs(lat_cc_pd(jj)-lat_volcano)+dn/2.0_ip)*DEG2KMLAT
              !ew_km  =(abs(lon_cc_pd(ii)-lon_volcano)+de/2.0_ip)* &
              !          cos(avg_lat*DEG2RAD)*DEG2KMLON
              radnow = sqrt(ns_km**2.0_ip+ew_km**2.0_ip)  !distance, km
              !make sure we're within the umbrella cloud
              if (radnow.lt.cloud_radius) then
                if(VelMod_umb.eq.1)then
                  ! This is Eq. 3 except the second term here has an extra r/R
                  ! Note that the radial windspeed function is just a non-dimensional radial
                  ! function scaling the of leading edge speed where the edge speed is from the
                  ! time derivitive of the edge position function.
                  rexp = 3.0_ip  ! This is from legacy code
                  !rexp = 2.0_ip  ! This is from Eq 3 of paper
                  windspeedhere = edge_speed  *                          &  ! m/s
                                  MPS_2_KMPHR *                          &  ! km/hr
                                  (3.0_ip/4.0_ip)*(cloud_radius/radnow)* &  ! Start of scaling term
                                  (1.0_ip                              + &  ! Cons. of Vol term
                                   (1.0_ip/3.0_ip)*(radnow/cloud_radius)**rexp) ! outward spreading,time flattening
                elseif(VelMod_umb.eq.2)then
                  ! This is Eq. 4
                  windspeedhere = edge_speed  *                          &  ! m/s
                                  MPS_2_KMPHR *                          &  ! km/hr
                                  (cloud_radius/radnow)**0.5_ip             ! 
                endif

                thetanow = atan2(ns_km,ew_km)     !angle CW from E
                do iz=ibase,itop
                  ! These velocities are in km/hr
                  uvx_pd(ii,jj,iz)=windspeedhere*cos(thetanow)
                  uvy_pd(ii,jj,iz)=windspeedhere*sin(thetanow)
                enddo
              endif
            enddo
          enddo
        endif
      endif

      return

      end subroutine umbrella_winds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  TephraSourceNodes_Umbrella
!
!  Called from: Ash3d.F90
!  Arguments:
!    none
!
!  This subroutine takes SourceNodeFlux calculated for the vent column in
!  TephraSourceNodes and spreads it over the 3x3 source patch for the umbrella.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine TephraSourceNodes_Umbrella

      use mesh,          only : &
         nzmax

      use Source,        only : &
         SourceNodeFlux

      use Tephra,        only : &
         n_gs_max

      integer :: i,j,k

      SourceNodeFlux_Umbrella(1:3,1:3,1:nzmax,1:n_gs_max) = 0.0_ip
      do k=1,nzmax+1
        if(k.lt.ibase)then
          ! Below the cloud, use the normal Suzuki profile
          SourceNodeFlux_Umbrella(2,2,k,1:n_gs_max)=SourceNodeFlux(k,1:n_gs_max)
        elseif(k.le.itop)then
          ! Above the umbrella base, but below its top, spread the source over the
          ! 3x3 patch
          do i=1,3
            do j=1,3
              SourceNodeFlux_Umbrella(i,j,k,1:n_gs_max) = &
                SourceNodeFlux(k,1:n_gs_max)*ScaleFac_Umbrella(k)
            enddo
          enddo
        else
          ! Above the umbrella top, zero out source
          SourceNodeFlux_Umbrella(2,2,k,1:n_gs_max) = 0.0_ip
        endif
      enddo

      return

      end subroutine TephraSourceNodes_Umbrella

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  SourceVolInc_Umbrella(dt)
!
!  Called from: Ash3d.F90
!  Arguments:
!    dt = time step in hours
!
!  This function calculates the total tephra volume inserted in this time step.
!  It is used only for mass-conservation error-checking.  It does the same
!  work as SourceVolInc but also adds the umbrella bit over the 3x3 zone.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function SourceVolInc_Umbrella(dt)

      use global_param,  only : &
         KM3_2_M3

      use mesh,          only : &
         kappa_pd,ivent,jvent

      use Tephra,        only : &
         n_gs_max,MagmaDensity

      use Source,        only : &
         SourceNodeFlux

      real(kind=ip) :: SourceVolInc_Umbrella
      real(kind=dp) :: dt

      real(kind=ip) :: tmp
      integer :: i,j,ii,jj,k,isize

      tmp = 0.0_ip

      do isize=1,n_gs_max
        do k=1,ibase-1
          tmp = tmp                             + & ! final units is km3
                real(dt,kind=ip)                * & ! hr
                SourceNodeFlux(k,isize)         * & ! kg/km3 hr
                kappa_pd(ivent,jvent,k)         / & ! km3
                MagmaDensity                    / & ! kg/m3
                KM3_2_M3                            ! m3/km3
        enddo
      enddo

      do i=1,3
        ii=ivent-2+i
        do j=1,3
          jj=jvent-2+j
          do k=ibase,itop
            do isize=1,n_gs_max
              tmp= tmp                                           + & ! final units is km3
                real(dt,kind=ip)                                 * & ! hr      
                SourceNodeFlux(k,isize)*AvgStenc_Umbrella(i,j,k) * & ! kg/km3 hr
                kappa_pd(ivent,jvent,k)                          / & ! km3
                MagmaDensity                                     / & ! kg/m3
                KM3_2_M3                                             ! m3/km3
            enddo
          enddo
        enddo
      enddo

      SourceVolInc_Umbrella = tmp

      return

      end function SourceVolInc_Umbrella

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  AvgCon_Umbrella(conpatch,klevel)
!
!  Called from: Ash3d.F90
!  Arguments:
!    conpatch = 3x3 patch over which we are averaging concentrations
!    klevel   = z-index of patch (used to get volumes)
!
!  This function is used for averaging the concentration in the 3x3 patch
!  surrounding the vent within the umbrella cloud.  This is applied before
!  source insertion and helps smooth the solution.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function AvgCon_Umbrella(conpatch,klevel)

      real(kind=ip) :: AvgCon_Umbrella
      real(kind=ip),dimension(3,3) :: conpatch
      integer :: klevel

      integer :: i,j

      AvgCon_Umbrella = 0.0_ip
      do i=1,3
        do j=1,3
          AvgCon_Umbrella = AvgCon_Umbrella + conpatch(i,j)*AvgStenc_Umbrella(i,j,klevel)
        enddo
      enddo

      return

      end function AvgCon_Umbrella

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module Source_Umbrella
!##############################################################################

