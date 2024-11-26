!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  calc_mesh_params
!
!  Called from: Ash3d.F90 just after arrays are allocated
!  Arguments:
!    none
!
!  This subroutine defines the mesh used by Ash3d and fills all the variables from
!  the mesh data module.  Cartesian grids are treated separately with the
!  horizontal cell-centered coordinate variables x_cc_pd and y_cc_pd.  Lon/Lat grids
!  use lon_cc_pd and lat_cc_pd.  Additionally, both coordinate systems use the same
!  geometry arrays for volume (kappa_pd) and surface area (sigma_[nx,ny,nz]_pd).
!  The vertical coordinate is calculated from z_vec_init (built in Read_Control_File),
!  populating the variables z_cc_pd (cell-centered vert. grid), dz_vec_pd (thickness
!  of each cell), and z_lb_pd (height of lower boundary).  Next, the projection
!  information for the computational grid is sent to MetReader and the corresponding
!  grids are built with a call to call MR_Initialize_Met_Grids.  Lastly, a call
!  to MR_Set_Met_Times informs MetReader of the start and simulation time.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine calc_mesh_params

      ! set up grids for solution and wind arrays using x,y (or lon,lat) and z in km

      use precis_param

      use io_units

      use global_param,  only : &
         DEG2KMLON,DEG2RAD,DEG2KMLAT,RAD_EARTH,PI

      use mesh,          only : &
         nxmax,nymax,nzmax,x_cc_pd,y_cc_pd,dx,dy,&
         lon_cc_pd,lat_cc_pd,de,dn,de_km,dn_km, &
         z_lb_pd,z_cc_pd,dz_vec_pd,z_vec_init, &
         sigma_nx_pd,sigma_ny_pd,sigma_nz_pd,kappa_pd,&
         xLL,yLL,latLL,lonLL, &
         A3d_iprojflag,A3d_k0_scale,A3d_phi0,A3d_lam0,A3d_phi1,&
         A3d_phi2,A3d_Re,IsLatLon,IsPeriodic

      use time_data,     only : &
         SimStartHour,Simtime_in_hours

      use Source,        only : &
         lat_volcano

      use MetReader,     only : &
           MR_Set_CompProjection, &
           MR_Initialize_Met_Grids, &
           MR_Set_Met_Times, &
           MR_Set_SigmaAlt_Scaling

      implicit none

      integer :: i,j,k

      real(kind=ip) :: r_1,r_2,rr_2,rr_1,drr,drrr
      real(kind=ip) :: phi_1,phi_2
      real(kind=ip) :: theta_1,theta_2,del_theta,del_costheta
      real(kind=ip) :: del_lam
      real(kind=ip) :: phi_bot,phi_top,phi
      real(kind=sp),allocatable,dimension(:) :: dumx_sp,dumy_sp,dumz_sp

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"--------------------------------------------------"
        write(outlog(io),*)"---------- CALC_MESH_PARAMS ----------------------"
        write(outlog(io),*)"--------------------------------------------------"
      endif;enddo

     ! Set up concentration grid.
     ! The concentration is calculated at cell centers
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),30)
      endif;enddo

      do k=1,nzmax
        ! Note: z_vec_init is the array (0:nz_init) for the initial input z profile.
        !       This could be much larger than e_PlumeHeight since it could describe
        !       a bunch of linear segments up to 50 km or something.  We will use
        !       the cell-centered value and the lower boundary value for all subsequent
        !       calculations, limited to just the height needed (e_PlumeHeight*ZPADDING)
        dz_vec_pd(k) = z_vec_init(k) - z_vec_init(k-1)
        z_cc_pd(k)   = z_vec_init(k) - 0.5_ip*dz_vec_pd(k)
        z_lb_pd(k)   = z_cc_pd(k)    - 0.5_ip*dz_vec_pd(k)  ! This is essentially just z_vec_init(k-1)
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
        ! Find width and height of a node (km) at the volcano's location
        de_km = de*cos(lat_volcano*DEG2RAD)*DEG2KMLON
        dn_km = dn*DEG2KMLAT
        do i=-1,nxmax+2
          lon_cc_pd(i) = lonLL + de*real(i,kind=ip) - de*0.5_ip
        enddo
        do j=-1,nymax+2
          lat_cc_pd(j) = latLL + dn*real(j,kind=ip) - dn*0.5_ip
        enddo

        !*********************************************************************************
        do k=-1,nzmax+2
          r_1  = RAD_EARTH+z_cc_pd(k)-0.5_ip*dz_vec_pd(k)  ! r at bottom of cell
          r_2  = RAD_EARTH+z_cc_pd(k)+0.5_ip*dz_vec_pd(k)  ! r at top of cell

          rr_1 =      r_1*r_1
          rr_2 =      r_2*r_2
          drr  = (    r_2*r_2 -     r_1*r_1)  ! difference of squares
          drrr = (r_2*r_2*r_2 - r_1*r_1*r_1)  ! difference of cubes
          do j=-1,nymax+2
            phi_1 = DEG2RAD*(lat_cc_pd(j)-0.5_ip*dn) ! phi at lower y of cell
            phi_2 = DEG2RAD*(lat_cc_pd(j)+0.5_ip*dn) ! phi at upper y of cell
            theta_1 = 0.5_ip*PI-phi_1
            theta_2 = 0.5_ip*PI-phi_2
            del_theta  = dn*DEG2RAD
            del_costheta = cos(theta_2)-cos(theta_1)
            ! No need to loop over lambda (longitude) since each slice is identical
            del_lam = de*DEG2RAD

              ! Area of face at i-1/2,j,k
            sigma_nx_pd(-1:nxmax+2,j,k) = 0.5_ip*del_theta*drr
              ! Area of face at i,j-1/2,k
            sigma_ny_pd(-1:nxmax+2,j,k) = 0.5_ip*sin(theta_1)*del_lam*drr
              ! Area of face at i,j,k-1/2
            sigma_nz_pd(-1:nxmax+2,j,k) = rr_1*del_lam*del_costheta
              ! Volume of cell
            kappa_pd(-1:nxmax+2,j,k)=del_lam*drrr*del_costheta/3.0_ip

          enddo
        enddo
        !*********************************************************************************
        ! For lon/lat cases, we will need the dx,dy that is most restrictive on
        ! dt calculations
        phi_bot = DEG2RAD*(lat_cc_pd(      1) - dn*0.5_ip)
        phi_top = DEG2RAD*(lat_cc_pd(nymax+1) + dn*0.5_ip)
        phi = max(abs(phi_bot),abs(phi_top))
        dx = RAD_EARTH*de*DEG2RAD*cos(phi)
        dy = RAD_EARTH*dn*DEG2RAD
      else ! This is the .not.IsLatLon case
          ! Set up cell-centered coordinates
        do i=-1,nxmax+2
          x_cc_pd(i) = xLL + dx*real(i,kind=ip) - dx*0.5_ip
        enddo
        do j=-1,nymax+2
          y_cc_pd(j) = yLL + dy*real(j,kind=ip) - dy*0.5_ip
        enddo
          ! Area of face at i,j,k-1/2
        sigma_nz_pd(:,:,:) = dy*dx
        do k=-1,nzmax+2
            ! Area of face at i-1/2,j,k
          sigma_nx_pd(:,:,k) = dz_vec_pd(k)*dy
            ! Area of face at i,j-1/2,k
          sigma_ny_pd(:,:,k) = dz_vec_pd(k)*dx
            ! Volume of cell
          kappa_pd(-1:nxmax+2,-1:nymax+2,k)=dx*dy*dz_vec_pd(k)
        enddo

      endif ! IsLatLon

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"    Cell-centered computational grid extends from:"
        if (IsLatLon) then
          write(outlog(io),'(a13,2f10.3)')"     in lon: ",lon_cc_pd(1),lon_cc_pd(nxmax)
          write(outlog(io),'(a13,2f10.3)')"     in lat: ",lat_cc_pd(1),lat_cc_pd(nymax)
          write(outlog(io),'(a13,2f10.3)')"     in   z: ",   z_cc_pd(1), z_cc_pd(nzmax)
        else
          write(outlog(io),'(a11,2f10.3)')"     in x: ",x_cc_pd(1),x_cc_pd(nxmax)
          write(outlog(io),'(a11,2f10.3)')"     in y: ",y_cc_pd(1),y_cc_pd(nymax)
          write(outlog(io),'(a11,2f10.3)')"     in z: ",z_cc_pd(1),z_cc_pd(nzmax)
        endif
      endif;enddo

      ! Now setting up projection information
      ! Evaluate wind files for time/space consistency with model requests
        ! Initialize the grids needed for met data
      call MR_Set_CompProjection(IsLatLon,A3d_iprojflag,A3d_lam0, &
                                 A3d_phi0,A3d_phi1,A3d_phi2,       &
                                 A3d_k0_scale,A3d_Re)
      allocate(dumx_sp(nxmax))
      allocate(dumy_sp(nymax))
      allocate(dumz_sp(nzmax))
      if(IsLatLon)then
        dumx_sp(1:nxmax) = real(lon_cc_pd(1:nxmax),kind=sp)
        dumy_sp(1:nymax) = real(lat_cc_pd(1:nymax),kind=sp)
        dumz_sp(1:nzmax) = real(  z_cc_pd(1:nzmax),kind=sp)
      else
        dumx_sp(1:nxmax) = real(  x_cc_pd(1:nxmax),kind=sp)
        dumy_sp(1:nymax) = real(  y_cc_pd(1:nymax),kind=sp)
        dumz_sp(1:nzmax) = real(  z_cc_pd(1:nzmax),kind=sp)
      endif
      call MR_Initialize_Met_Grids(nxmax,nymax,nzmax,             &
                              dumx_sp,dumy_sp,dumz_sp,            &
                              IsPeriodic)

      deallocate(dumx_sp,dumy_sp,dumz_sp)

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Finished initializing Met Grids"
      endif;enddo

      ! Informing MetReader of the start and simulation time
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Now determining which NWP files and steps needed."
      endif;enddo
      call MR_Set_Met_Times(SimStartHour, Simtime_in_hours)

30    format(/,4x,'Calculating the locations of each cell-centered node in the grid.')
  
      end subroutine calc_mesh_params

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  calc_s_mesh
!
!  Called from: Ash3d.F90 just after topography is loaded. If topography is not
!  used (ZScaling_ID=0; default), then s = z. If ZScaling_ID=1 (shifted topo) or
!  ZScaling_ID=2 (scaled topo), then s is a bit different.
!
!  Arguments:
!    none
!
!  Assigns:
!    s_cc_pd,s_lb_pd,ds_vec_pd and Met variables through MR_Set_SigmaAlt_Scaling
!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine calc_s_mesh

      ! set up grids for solution and wind arrays

      use precis_param

      use global_param,  only : &
         CFL

      use io_units

      use mesh,          only : &
         nxmax,nymax,nzmax,j_cc_pd,sigma_nx_pd,sigma_ny_pd,kappa_pd,&
         z_lb_pd,z_cc_pd,dz_vec_pd, &
         s_lb_pd,s_cc_pd,ds_vec_pd,Ztop,ZScaling_ID, &
         IsLatLon,Zsurf,ivent,jvent

      use MetReader,     only : &
           MR_Set_CompProjection, &
           MR_Initialize_Met_Grids, &
           MR_Set_Met_Times, &
           MR_Set_SigmaAlt_Scaling

      implicit none

      real(kind=sp),allocatable,dimension(:) :: dums_sp
      integer       :: i,j
      real(kind=ip) :: j_max

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"--------------------------------------------------"
        write(outlog(io),*)"---------- CALC_S_MESH ---------------------------"
        write(outlog(io),*)"--------------------------------------------------"
      endif;enddo

      ! s_cc is the replacement z-coordinate: either s=z or s=(z-zsurf)/(ztop-zsurf)
      ! Here we set up the s-values with zsurf=0 and might apply topography later
      if(ZScaling_ID.eq.0)then
        ! No topography modification; s=z in km
        s_cc_pd(:) = z_cc_pd(:)
        s_lb_pd(:) = z_lb_pd(:)
      elseif(ZScaling_ID.eq.1)then
        ! Coordinate system will be shifted by topography (surf -> surf + Ztop in km)
        ! Since we want the s array to start at the surface (not just the s-value of
        ! the z-array), we use the same values as z with the understanding that it is
        ! shifted
        s_cc_pd(:) = z_cc_pd(:)
        s_lb_pd(:) = z_lb_pd(:)
      elseif(ZScaling_ID.eq.2)then
        ! Coordinate system will be scaled by topography (0->Ztop in km)
        ! Same logic as above; surface to top should be 0->Ztop so use Zsurf of 0
        ! for the purpose of indexing (otherwise we just have the s-values of
        ! the Cartesian z-array)
        s_cc_pd(:) = z_cc_pd(:)
        s_lb_pd(:) = z_lb_pd(:)
      endif

      if(ZScaling_ID.eq.0.or.ZScaling_ID.eq.1)then
        ds_vec_pd(:) = dz_vec_pd(:)
      else
        ds_vec_pd(:) = (dz_vec_pd(:))*Ztop/(Ztop-Zsurf(ivent,jvent))
      endif

      allocate(dums_sp(nzmax))
      if(IsLatLon)then
        dums_sp(1:nzmax) = real(  s_cc_pd(1:nzmax),kind=sp)
      else
        dums_sp(1:nzmax) = real(  s_cc_pd(1:nzmax),kind=sp)
      endif
      call MR_Set_SigmaAlt_Scaling(nzmax,dums_sp)

      deallocate(dums_sp)

      ! Now modify all the cell face area and volume measures if we are using
      ! the scaled coordinates. Note: these variables are in km2 and km3
      ! since the jacobian for all cases is unitless
      if(ZScaling_ID.eq.2)then
        ! The Jacobian measures the change in volumen of the cell which is
        ! accommodated by a change in thickness. We have to scale the CFL to
        ! be stable with these squeezed cells
        j_max = maxval(j_cc_pd(:,:))
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Cells are thinned by up to ",real(j_max,kind=4)
          write(outlog(io),*)"Resetting CFL from ",real(CFL,kind=4)
          write(outlog(io),*)"  to ",real(CFL * j_max,kind=4)
        endif;enddo
        CFL = CFL * j_max
        do i=-1,nxmax+2
          do j=-1,nymax+2
            ! Note that sigma_nz_pd is unaffected
              ! Area of face at i-1/2,j,k
            sigma_nx_pd(i,j,:) = sigma_nx_pd(i,j,:)*j_cc_pd(i,j)
              ! Area of face at i,j-1/2,k
            sigma_ny_pd(i,j,:) = sigma_ny_pd(i,j,:)*j_cc_pd(i,j)
              ! Volume of cell
            kappa_pd(i,j,:)    = kappa_pd(i,j,:)*j_cc_pd(i,j)
          enddo
        enddo
      endif
      end subroutine calc_s_mesh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  get_minmax_lonlat
!
!  Called from: Not called from standard Ash3d, but used in some optional modules
!  Arguments:
!    lonmin,lonmax : minimum and maximum longitude of the computational grid
!    latmin,latmax : minimum and maximum latitude of the computational grid
!
!  This subroutine has no input arguments, but returns min/max values for lon/lat
!  for a projected grid.  Projected grids will have curved boundaries in lon/lat
!  coordinates so the lon/lat minima/maxima might be along these boundaries.  This
!  subroutine loops through all the cells, inverse-projects the coordinates, then
!  returns the min/max values.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_minmax_lonlat(lonmin,lonmax,latmin,latmax)

      use precis_param

      use io_units

      use mesh,              only : &
         nxmax,nymax,x_cc_pd,y_cc_pd,xy2ll_xlon,xy2ll_ylat,&
         A3d_iprojflag,A3d_k0_scale,A3d_phi0,A3d_lam0,A3d_phi1,&
         A3d_phi2,A3d_Re

      use projection,        only : &
           PJ_proj_inv

      implicit none

      real(kind=8),intent(out) :: lonmin
      real(kind=8),intent(out) :: lonmax
      real(kind=8),intent(out) :: latmin
      real(kind=8),intent(out) :: latmax

      integer       :: i,j
      real(kind=8)  :: olam,ophi ! using precision needed by libprojection
      real(kind=8)  :: xin,yin

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"Inside get_minmax_lonlat"
        write(outlog(io),*)"Allocating of size: ",nxmax+2,nymax+2
      endif;enddo
#ifdef USEPOINTERS
      if(.not.associated(xy2ll_ylat))then
#else
      if(.not.allocated(xy2ll_ylat))then
#endif
        allocate(xy2ll_ylat(0:nxmax+1,0:nymax+1))
      endif
#ifdef USEPOINTERS
      if(.not.associated(xy2ll_xlon))then
#else
      if(.not.allocated(xy2ll_xlon))then
#endif
        allocate(xy2ll_xlon(-1:nxmax+2,-1:nymax+2))
      endif
      ! This block calculates the lon/lat for each computational grid point
      ! Note:  All we need here is just the min/max for lat/lon so that we
      !        can generate our own, regular lat/lon grid filled with
      !        interpolated values.
      latmax =  -90.0_8
      latmin =   90.0_8
      lonmin =  360.0_8
      lonmax = -360.0_8
      do i=-1,nxmax+2
        do j=-1,nymax+2
          xin = real(x_cc_pd(i),kind=dp)  ! Projection routines use kind=8
          yin = real(y_cc_pd(j),kind=dp)
          call PJ_proj_inv(xin,yin, &
                         A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
                         A3d_k0_scale,A3d_Re, &
                         olam,ophi)
          if(olam.lt.lonmin)lonmin=olam
          if(olam.gt.lonmax)lonmax=olam
          if(ophi.lt.latmin)latmin=ophi
          if(ophi.gt.latmax)latmax=ophi
          xy2ll_ylat(i,j) = real(ophi,kind=ip)
          xy2ll_xlon(i,j) = real(olam,kind=ip)
          if(xy2ll_xlon(i,j).lt.0.0_ip) xy2ll_xlon(i,j) = xy2ll_xlon(i,j) + 360.0_ip
        enddo
      enddo

      ! Get the extremal extents in lat/lon space
      !latmin = minval(minval(xy2ll_ylat,1),1)
      !latmax = maxval(maxval(xy2ll_ylat,1),1)
      !lonmin = minval(minval(xy2ll_xlon,1),1)
      !lonmax = maxval(maxval(xy2ll_xlon,1),1)

      end subroutine get_minmax_lonlat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  get_minmax_index
!
!  Called from: AdvectHorz
!  Arguments:
!    none
!
!  This subroutine is called every time the horizontal advection routine is
!  called if the executable was compiled with FAST_SUBGRID set.  This subroutine
!  loops through all the cells in the domain and tracks the maximal and minimal
!  extents of the airborne ash cloud, setting imin,imax,jmin,jmax,kmin,kmax.
!  These values define the range of the do-loops for this time step, greatly
!  reducing the computation time for small clouds in large domains.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_minmax_index

      use precis_param

      use io_units

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,ts0,ivent,jvent,kvent,IsPeriodic

      use solution,      only : &
         concen_pd,imin,imax,jmin,jmax,kmin,kmax

      use Output_Vars,  only : &
         CLOUDCON_GRID_THRESH
 
      implicit none

      integer :: i,j,k
      real(kind=ip) :: tmp_flt

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine get_minmax_index"
      endif;enddo

      ! Find extent of ash cloud making sure to include the vent
      !  First in x
      imin = 1
      imax = nxmax
      do i=2,ivent-1
        tmp_flt=maxval(concen_pd(i,1:nymax,1:nzmax,1:nsmax,ts0))
        if(tmp_flt.lt.CLOUDCON_GRID_THRESH)then
          imin = i
        else
          exit
        endif
      enddo
      do i=nxmax,imin+2,-1
        tmp_flt=maxval(concen_pd(i,1:nymax,1:nzmax,1:nsmax,ts0))
        if(tmp_flt.lt.CLOUDCON_GRID_THRESH)then
          imax = i
        else
          exit
        endif
      enddo
      if(IsPeriodic)then
        ! if this is a periodic case and the ash touches a boundary, set the min and max
        ! to the full domain
        if(imin.eq.1.or.imax.eq.nxmax)then
          imin = 1
          imax = nxmax
        endif
      endif
      !  Now in y
      jmin = 1
      jmax = nymax
      do j=2,jvent-1
        tmp_flt=maxval(concen_pd(1:nxmax,j,1:nzmax,1:nsmax,ts0))
        if(tmp_flt.lt.CLOUDCON_GRID_THRESH)then
          jmin = j
        else
          exit
        endif
      enddo
      do j=nymax,jmin+2,-1
        tmp_flt=maxval(concen_pd(1:nxmax,j,1:nzmax,1:nsmax,ts0))
        if(tmp_flt.lt.CLOUDCON_GRID_THRESH)then
          jmax = j
        else
          exit
        endif
      enddo
      !  Now in z
      kmin = 1
      kmax = nzmax
      do k=2,kvent-1
        tmp_flt=maxval(concen_pd(1:nxmax,1:nymax,k,1:nsmax,ts0))
        if(tmp_flt.lt.CLOUDCON_GRID_THRESH)then
          kmin = k
        else
          exit
        endif
      enddo
      do k=nzmax,kmin+2,-1
        tmp_flt=maxval(concen_pd(1:nxmax,1:nymax,k,1:nsmax,ts0))
        if(tmp_flt.lt.CLOUDCON_GRID_THRESH)then
          kmax = k
        else
          exit
        endif
      enddo
      if(kmax.lt.kmin)then
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"ERROR: kmax<kmin in get_minmax_index"
        endif;enddo
        stop 1
      endif
      ! Advection routines expect at least three cells. Adjust i,j,k as needed
      if(imax-imin.lt.2)then
        imax = min(imax+1,nxmax)
        imin = imax-2
      endif
      if(jmax-jmin.lt.2)then
        jmax = min(jmax+1,nymax)
        jmin = jmax-2
      endif
      if(kmax-kmin.lt.2)then
        kmax = min(kmax+1,nzmax)
        kmin = kmax-2
      endif

      end subroutine get_minmax_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

