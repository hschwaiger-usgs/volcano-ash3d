      subroutine calc_mesh_params

      ! set up grids for solution and wind arrays

      use precis_param

      use io_units

      use global_param,  only : &
         DEG2KMLON,DEG2RAD,DEG2KMLAT,RAD_EARTH,PI

      use mesh,          only : &
         nxmax,nymax,nzmax,x_cc_pd,y_cc_pd,dx,dy,&
         lon_cc_pd,lat_cc_pd,de,dn,de_km,dn_km, &
         z_lb_pd,z_vec_init,z_cc_pd,dz_vec_pd, &
         sigma_nx_pd,sigma_ny_pd,sigma_nz_pd,kappa_pd,&
         xLL,yLL,latLL,lonLL,rdphi_pd,rdlambda_pd, &
         A3d_iprojflag,A3d_k0_scale,A3d_phi0,A3d_lam0,A3d_phi1,&
         A3d_phi2,A3d_radius_earth,IsLatLon,IsPeriodic

      use time_data,     only : &
         SimStartHour,Simtime_in_hours

      use Source,        only : &
         lat_volcano,SourceNodeHeight_km,SourceType,SourceNodeWidth_km

      use MetReader,     only : &
           MR_Set_CompProjection, &
           MR_Initialize_Met_Grids, &
           MR_Set_Met_Times

      implicit none

      integer :: i,j,k

      real(kind=ip) :: r_1,r_2,rr_2,drr,drrr
      real(kind=ip) :: phi_1,phi_2,del_phi,del_cphi
      real(kind=ip) :: del_lam
      real(kind=ip) :: phi

      write(global_production,*)"--------------------------------------------------"
      write(global_production,*)"---------- CALC_MESH_PARAMS ----------------------"
      write(global_production,*)"--------------------------------------------------"

!     SET UP CONCENTRATION GRID.
!     THE CONCENTRATION IS CALCULATED AT CELL CENTERS
      write(global_info,30)
      write(global_log,30)

      do k=1,nzmax
        ! Note: z_vec_init is the array (0:nz_init) for the initial input z profile.
        !       This could be much larger than PlumeHeight since it could describe
        !       a bunch of linear segments up to 50 km or something.  We will use
        !       the cell-centered value and the lower boundary value for all subsequent
        !       calculations, limited to just the height needed (PlumeHeight*ZPADDING)
        dz_vec_pd(k) = z_vec_init(k) - z_vec_init(k-1)
        z_cc_pd(k)   = z_vec_init(k) - 0.5_ip*dz_vec_pd(k)
        z_lb_pd(k)   = z_cc_pd(k) - 0.5_ip*dz_vec_pd(k)  ! This is essentially just z_vec_init(k-1)
      enddo
      z_cc_pd( 0)      = z_cc_pd(1)     - dz_vec_pd(1)
      z_lb_pd( 0)      = z_cc_pd(0)     - dz_vec_pd(1)
      z_cc_pd(-1)      = z_cc_pd(0)     - dz_vec_pd(1)
      z_lb_pd(-1)      = z_cc_pd(-1)    - dz_vec_pd(1)
      z_cc_pd(nzmax+1) = z_cc_pd(nzmax) + dz_vec_pd(nzmax)
      z_cc_pd(nzmax+2) = z_cc_pd(nzmax) + dz_vec_pd(nzmax)*2.0_ip
      z_lb_pd(nzmax+1) = z_lb_pd(nzmax) + dz_vec_pd(nzmax)
      z_lb_pd(nzmax+2) = z_lb_pd(nzmax) + dz_vec_pd(nzmax)*2.0_ip
      dz_vec_pd(-1)      = dz_vec_pd(1)
      dz_vec_pd( 0)      = dz_vec_pd(1)
      dz_vec_pd(nzmax+1) = dz_vec_pd(nzmax)
      dz_vec_pd(nzmax+2) = dz_vec_pd(nzmax)

      if (IsLatLon) then
        !find width and height of a node (km) at the volcano's location
        de_km = de*cos(lat_volcano*DEG2RAD)*DEG2KMLON
        dn_km = dn*DEG2KMLAT
        if (SourceType.eq.'umbrella') then
           !calculate radius of umbrella cloud source nodes
           SourceNodeWidth_km  = (3./2.)*de_km
           SourceNodeHeight_km = (3./2.)*dn_km
           write(global_info,142) SourceNodeWidth_km, SourceNodeHeight_km
142        format(/,'Calculating width of source for umbrella.',/, &
                  'SourceNodeWidth_km=',f5.1,/, &
                  'SourceNodeHeight_km=',f5.1)
        endif
        do i=-1,nxmax+2
          lon_cc_pd(i) = lonLL + de*real(i,kind=ip) - de*0.5_ip
        enddo
        do j=-1,nymax+2
          lat_cc_pd(j) = latLL + dn*real(j,kind=ip) - dn*0.5_ip
        enddo

        !*********************************************************************************
        do k=-1,nzmax+2
          rdphi_pd(k) = dn*DEG2RAD * (RAD_EARTH + z_cc_pd(k))
          r_1  = RAD_EARTH+z_cc_pd(k)               ! r at bottom of cell
          r_2  = RAD_EARTH+z_cc_pd(k)+dz_vec_pd(k)  ! r at top of cell

          rr_2 = r_2*r_2
          drr  = (    r_2*r_2 -     r_1*r_1)  ! difference of squares
          drrr = (r_2*r_2*r_2 - r_1*r_1*r_1)  ! difference of cubes
          do j=-1,nymax+2
            !rdlambda(j,k) = (RAD_EARTH+k*dz-0.5_ip*dz)*cos(DEG2RAD*(gridlat(j)-0.5_ip*dn))*de*DEG2RAD
            rdlambda_pd(j,k) = (RAD_EARTH+z_cc_pd(k))*cos(DEG2RAD*(lat_cc_pd(j)-0.5_ip*dn))*de*DEG2RAD
            phi_1 = DEG2RAD*(lat_cc_pd(j)-0.5_ip*dn)
            phi_2 = DEG2RAD*(lat_cc_pd(j)+0.5_ip*dn)
            del_phi  = dn*DEG2RAD
            del_cphi = cos(0.5_ip*PI-phi_2)-cos(0.5_ip*PI-phi_1)
            ! No need to loop over lambda since each slice is identical
            del_lam = de*DEG2RAD
              ! Area of face at i-1/2,j,k
            sigma_nx_pd(-1:nxmax+2,j,k) = 0.5_ip*del_phi*drr
              ! Area of face at i,j-1/2,k
            sigma_ny_pd(-1:nxmax+2,j,k) = 0.5_ip*sin(0.5_ip*PI-phi_2)*del_lam*drr
              ! Area of face at i,j,k-1/2
            sigma_nz_pd(-1:nxmax+2,j,k) = rr_2*del_lam*del_cphi
              ! Volume of cell
            kappa_pd(-1:nxmax+2,j,k)=del_lam*drrr*del_cphi/3.0_ip
          enddo
        enddo
       !*********************************************************************************
        ! For lon/lat cases, we will need the dx,dy that is most restrictive on
        ! dt calculations
        phi = maxval(abs(lat_cc_pd(1:nymax)))
        dx = RAD_EARTH*cos(DEG2RAD*phi)*de*DEG2RAD
        dy = RAD_EARTH*dn*PI/180.0_ip
      else
          ! Set up cell-centered coordinates
        do i=-1,nxmax+2
          x_cc_pd(i) = xLL + dx*real(i,kind=ip) - dx*0.5_ip
        enddo
        !xmin = minval(x_cc_pd)
        !xmax = maxval(x_cc_pd)
        do j=-1,nymax+2
          y_cc_pd(j) = yLL + dy*real(j,kind=ip) - dy*0.5_ip
        enddo
        !ymin = minval(y_cc_pd)
        !ymax = maxval(y_cc_pd)
        do k=-1,nzmax+2
          kappa_pd(-1:nxmax+2,-1:nymax+2,k)=dx*dy*dz_vec_pd(k)
        enddo

      endif

      write(global_info,*)"    Cell-centered computational grid extends from:"
      if (IsLatLon) then
        write(global_info,*)"     in lon: ",real(lon_cc_pd(1),kind=sp),real(lon_cc_pd(nxmax),kind=sp)
        write(global_info,*)"     in lat: ",real(lat_cc_pd(1),kind=sp),real(lat_cc_pd(nymax),kind=sp)
        write(global_info,*)"     in   z: ",real(   z_cc_pd(1),kind=sp),real(   z_cc_pd(nzmax),kind=sp)
      else
        write(global_info,*)"     in x: ",real(x_cc_pd(1),kind=sp),real(x_cc_pd(nxmax),kind=sp)
        write(global_info,*)"     in y: ",real(y_cc_pd(1),kind=sp),real(y_cc_pd(nymax),kind=sp)
        write(global_info,*)"     in z: ",real(z_cc_pd(1),kind=sp),real(z_cc_pd(nzmax),kind=sp)
      endif
      ! Evaluate wind files for time/space consistency with model requests
        ! Initialize the grids needed for met data
      write(global_info,*)"Calling MR_Set_CompProjection"
      call MR_Set_CompProjection(IsLatLon,A3d_iprojflag,A3d_lam0, &
                                 A3d_phi0,A3d_phi1,A3d_phi2,       &
                                 A3d_k0_scale,A3d_radius_earth)

      write(global_info,*)"Calling MR_Initialize_Met_Grids"
      if(IsLatLon)then
        call MR_Initialize_Met_Grids(nxmax,nymax,nzmax,          &
                                real(lon_cc_pd(1:nxmax),kind=sp), &
                                real(lat_cc_pd(1:nymax),kind=sp), &
                                real(z_cc_pd(1:nzmax) ,kind=sp), &
                                IsPeriodic)
      else
        call MR_Initialize_Met_Grids(nxmax,nymax,nzmax,          &
                                real(x_cc_pd(1:nxmax)    ,kind=sp), &
                                real(y_cc_pd(1:nymax)    ,kind=sp), &
                                real(z_cc_pd(1:nzmax)    ,kind=sp), &
                                IsPeriodic)
      endif
      write(global_info,*)"Finished initializing Met Grids"
      write(global_info,*)"Now determining which NWP files and steps needed."
      call MR_Set_Met_Times(SimStartHour, Simtime_in_hours)

30    format(/,4x,'Calculating the locations of each cell-centered node in the grid.')
!31    format(/,4x,&
!             'interpolating to find winds at regular grid elevations')
!32    format(/,4x,'Finding source nodes')
  
      end subroutine calc_mesh_params

!##############################################################################

      subroutine get_minmax_lonlat(lonmin,lonmax,latmin,latmax)

      use precis_param

      use io_units

      use global_param,      only : &
         VERB

      use mesh,              only : &
         nxmax,nymax,x_cc_pd,y_cc_pd,xy2ll_xlon,xy2ll_ylat,&
         A3d_iprojflag,A3d_k0_scale,A3d_phi0,A3d_lam0,A3d_phi1,&
         A3d_phi2,A3d_radius_earth

      use projection,        only : &
           PJ_proj_inv

      implicit none

      real(kind=ip),intent(out) :: lonmin
      real(kind=ip),intent(out) :: lonmax
      real(kind=ip),intent(out) :: latmin
      real(kind=ip),intent(out) :: latmax

      integer :: i,j
      real(kind=ip)      :: olam_ip,ophi_ip ! using internal precision
      real(kind=dp)       :: xout,yout

      if(VERB.gt.1)write(global_info,*)"Inside get_minmax_lonlat"
      if(VERB.gt.1)write(global_info,*)"Allocating of size: ",nxmax+2,nymax+2
      if(.not.allocated(xy2ll_ylat))allocate(xy2ll_ylat(0:nxmax+1,0:nymax+1))
      if(.not.allocated(xy2ll_xlon))allocate(xy2ll_xlon(0:nxmax+1,0:nymax+1))

      ! This block calculates the lon/lat for each computational grid
      ! point
      ! Note:  All we need here is just the min/max for lat/lon so that
      ! we can generate our own, regular lat/lon grid filled with
      ! interpolated values.
      do i=0,nxmax+1
        do j=0,nymax+1
          xout = x_cc_pd(i)
          yout = y_cc_pd(j)
          call PJ_proj_inv(xout,yout, &
                         A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
                         A3d_k0_scale,A3d_radius_earth, &
                         olam_ip,ophi_ip)
          xy2ll_ylat(i,j)=ophi_ip
          if(olam_ip.lt.0.0_ip)olam_ip=olam_ip+360.0_ip
          xy2ll_xlon(i,j)=olam_ip
        enddo
      enddo

      ! Get the extremal extents in lat/lon space
      latmin = minval(minval(xy2ll_ylat,1),1)
      latmax = maxval(maxval(xy2ll_ylat,1),1)
      lonmin = minval(minval(xy2ll_xlon,1),1)
      lonmax = maxval(maxval(xy2ll_xlon,1),1)

      end subroutine get_minmax_lonlat

!##############################################################################

