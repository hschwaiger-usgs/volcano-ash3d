      subroutine umbrella_winds(first_time)

      !subroutine that calculates radial winds from the center of an umbrella cloud

      use precis_param

      use io_units

      use global_param,  only : &
         DEG2KMLAT,DEG2KMLON,DEG2RAD,KM_2_M,PI,HR_2_S,MPS_2_KMPHR,EPS_SMALL

      use time_data,     only : &
         time,Simtime_in_hours,dt

      use mesh,          only : &
         nxmax,nymax,dn_km,de_km,IsLatLon,ivent,jvent,&
         lat_cc_pd,lon_cc_pd

      use Source,        only : &
         uvx_pd,uvy_pd,ibase,itop,lat_volcano,lon_volcano,&
         MassFlux,SourceNodeWidth_km, SourceNodeHeight_km,&
         e_EndTime, SourceType

      use Tephra,        only : &
         n_gs_max

      implicit none

      logical, intent(in)  :: first_time

      real(kind=ip):: avg_lat          !avg latitude between vent & point
      real(kind=ip):: C_Costa          !C constant used in Costa et al., 2013
      real(kind=ip):: cloud_radius     !cloud radius, km
      real(kind=ip):: cloudrad_raw     !cloud radius, km, uncorrected
      real(kind=ip):: edge_speed       !expansion rate of cloud edge, m/s
      real(kind=ip):: etime_s          !time since eruption start, seconds
      real(kind=ip):: ew_km,ns_km      !distances between vent & point
      real(kind=ip):: k_entrainment    !entrainment coefficient
      real(kind=ip):: lambda           !umbrella cloud shape factor
      !real(kind=ip):: latnow, lonnow   !present latitude, longitude
      real(kind=ip):: massfluxnow      !current mass flux, kg/s
      real(kind=ip):: N_BV             !Brunt-Vaisala frequency, 1/s
      real(kind=ip):: qnow             !volume flow rate into umbrella cloud, m3/s
      real(kind=ip):: radnow           !radial distance from cloud center, km
      real(kind=ip):: thetanow         !angle of point CW from east
      real(kind=ip) :: windspeedhere    !windspeed at this node
      !real(kind=ip):: xyspacing        !average spacing between nodes, in km
      integer      :: ii,jj,iz         !counters
      integer      :: ew_nodes,ns_nodes!radius of clouds in nodes
      integer      :: west_node,east_node
      integer      :: south_node,north_node
      !character    :: answer*1

      !Set standard values
      lambda                        = 0.2_ip
      N_BV                          = 0.02_ip
      k_entrainment                 = 0.1_ip
      uvx_pd(-1:nxmax+2,-1:nymax+2,ibase:itop) = 0.0_ip               !set umbrella winds to zero
      uvy_pd(-1:nxmax+2,-1:nymax+2,ibase:itop) = 0.0_ip    
 
      if (.not.IsLatLon) then
        write(global_error,*) 'Error: umbrella_winds is not yet set up to handle'
        write(global_error,*) 'projected coordinates.'
        stop 1
      endif

      !call MassFluxCalculator
      ! Convert MassFlux from kg/hr to a local variable in kg/s
      massfluxnow = MassFlux(1)/HR_2_S        !mass flux rate, kg/s

      !If there is only one size class, assume it's an airborne run and
      !multiply the mass flux by 20
      if ((SourceType.eq.'umbrella_air').and.(n_gs_max.eq.1)) then
             massfluxnow = 20.0_ip*massfluxnow
      end if

      !set value of C based on latitude
      if (abs(lat_volcano).lt.23.0_ip) then
          !m3 kg^(-3/4) s^(-7/8) for tropical eruptions
        C_Costa = 0.43e3_ip
      else
          !m3 kg^(-3/4) s^(-7/8) for non-tropical eruptions
        C_Costa = 0.87e3_ip 
      endif

      ! Here is Eq. 2 of Mastin and Van Eaton, 2020 (m3/s)
      qnow  = C_Costa*sqrt(k_entrainment)*massfluxnow**(3.0_ip/4.0_ip) / &
              N_BV**(5.0_ip/4.0_ip)

      if (time.lt.EPS_SMALL)  then
          write(global_info,*) 
          write(global_info,*) 'in Umbrella_winds'
          write(global_info,*) '  massfluxnow (kg/s) = ',real(massfluxnow,kind=sp)
          write(global_info,*) '             C_Costa = ',real(C_Costa,kind=sp)
          write(global_info,*) '                N_BV = ',real(N_BV,kind=sp)
          write(global_info,*) '       k_entrainment = ',real(k_entrainment,kind=sp)
          write(global_info,*) '            Q (m3/s) = ',real(qnow,kind=sp)
          write(global_info,*)
          write(global_log,*) 
          write(global_log,*)  'in Umbrella_winds'
          write(global_log,*)  '  massfluxnow (kg/s) = ',real(massfluxnow,kind=sp)
          write(global_log,*)  '             C_Costa = ',real(C_Costa,kind=sp)
          write(global_log,*)  '                N_BV = ',real(N_BV,kind=sp)
          write(global_log,*)  '       k_entrainment = ',real(k_entrainment,kind=sp)
          write(global_log,*)  '            Q (m3/s) = ',real(qnow,kind=sp)
          write(global_log,*)
          !If we're doing an airborne run and using only 1 grain size,
          !multiply the MER by 20 to make sure we're getting the right umbrella
          !growth rate
          if ((SourceType.eq.'umbrella_air').and.(n_gs_max.eq.1)) then
              write(global_info,*) 'n_gs_max=1, so we are assuming an airborne run'
              write(global_info,*) 'massflux has been multiplied by 20'
              write(global_log,*) 'n_gs_max=1, so we are assuming an airborne run'
              write(global_log,*) 'massflux has been multiplied by 20'
          end if
       end if ! time.eq.0.0_ip

      !convert from  hours to seconds
!      if (itime.gt.0) then
!        if (time.gt.0.0_ip) then
        if(.not.first_time)then
          ! For the first time step, time=0 -> cloudrad=0 and uR=Inf
          ! so just use dt to get values at the end of step 1
          etime_s      = max(time,dt)*HR_2_S
        else
          etime_s      = min(Simtime_in_hours,e_EndTime(1))*HR_2_S
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
      cloudrad_raw = (3.0_ip*lambda*N_BV*qnow/(2.0_ip*PI))**(1.0_ip/3.0_ip) * &
                     etime_s**(2.0_ip/3.0_ip) / KM_2_M
      !Make sure cloud radius extends beyond the source nodes
      !cloud_radius = max(cloudrad_raw,max(SourceNodeWidth_km,SourceNodeHeight_km))
      cloud_radius = cloudrad_raw            !for debugging

      ! This is the time derivitive of Eq. 1 of Mastin and Van Eaton, 2020
      !cloud expansion rate, m/s
      edge_speed   = (2.0_ip/3.0_ip)*(3.0_ip*lambda*N_BV*qnow/(2.0_ip*PI))**(1.0_ip/3.0_ip) * &
                    etime_s**(-1.0_ip/3.0_ip)  

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
                ! This is Eq. 3 except the second term here has an extra r/R
                ! Note that the radial windspeed function is just a non-dimensional radial
                ! function scaling the of leading edge speed where the edge speed is from the
                ! time derivitive of the edge position function.
                windspeedhere = edge_speed  *                          &  ! m/s
                                MPS_2_KMPHR *                          &  ! km/hr
                                (3.0_ip/4.0_ip)*(cloud_radius/radnow)* &  ! Start of scaling term
                                (1.0_ip                              + &  ! Cons. of Vol term
                                 (1.0_ip/3.0_ip)*(radnow/cloud_radius)**3.0_ip) ! outward spreading,time flattening

!                windspeedhere = (3.0_ip/4.0_ip)*edge_speed* &          !m/s
!                    (cloud_radius/radnow)* &
!                    (1.0_ip+(1.0_ip/3.0_ip)*(radnow**3.0_ip/cloud_radius**3.0_ip))
                thetanow = atan2(ns_km,ew_km)     !angle CW from E
                do iz=ibase,itop
                  ! These velocities are in km/hr
                  uvx_pd(ii,jj,iz)=windspeedhere*cos(thetanow)
                  uvy_pd(ii,jj,iz)=windspeedhere*sin(thetanow)
                enddo
              endif
              !write(global_info,13) ii,jj,gridlat(jj),gridlon(ii), &
              !            radnow/cloud_radius, &
              !            uvx(ii,jj,ibase),uvy(ii,jj,ibase)
!13            !format(2i4,3f8.3,2f7.1)
              !if (radnow.lt.cloud_radius) then
              !   write(global_info,*) 'Continue?'
              !   read(5,'(a1)') answer
              !   if (answer.eq.'n') stop 1
              !endif
!              if (jj.eq.jvent)write(*,*)"UMBR",ii,ew_km,thetanow,uvx_pd(ii,jj,ibase)
            enddo
          enddo
        endif
      endif

      return

      end subroutine umbrella_winds
