      module Tephra

      use precis_param

      use io_units

      use global_param,  only : &
         useCalcFallVel,useVariableGSbins,PI,MPS_2_KMPHR,GRAV

      implicit none

      real(kind=ip) :: MagmaDensity   = 2500.0_ip  !density of magma, in kg/m3
      real(kind=ip) :: DepositDensity = 1000.0_ip  !deposit density, in kg/m3
      real(kind=ip) :: LAM_GS_THRESH  =  250.0_ip  ! Invokes Cslip once effect is 1%
                                     !=  125.0_ip  ! Invokes Cslip once effect is 2%
                                     !=   50.0_ip  ! Invokes Cslip once effect is 5%
      real(kind=ip) :: AIRBORNE_THRESH = 1.0e-3_ip ! Mass threshold for flagging bin as empty (1 gram)

      integer :: n_gs_max                    ! # size classes of particles 
      integer :: ns_aloft                    ! max gs bin still aloft

      real(kind=ip), dimension(:)  ,allocatable  :: Tephra_v_s         ! Settling vel
      real(kind=ip), dimension(:)  ,allocatable  :: Tephra_gsdiam      ! Grain-size diameter (read in mmm, stored in m)
      real(kind=ip), dimension(:)  ,allocatable  :: Tephra_bin_mass    ! mass    (kg)
      real(kind=ip), dimension(:)  ,allocatable  :: Tephra_rho_m       ! density (kg/m3)
      real(kind=ip), dimension(:)  ,allocatable  :: Tephra_gsF         ! Grain-size shape (b+c)/2a
      real(kind=ip), dimension(:,:),allocatable  :: Tephra_gsF_fac     ! Precalculated shape factors
                                                                !  i=1 WH Stokes fac
                                                                !  i=2 WH Newton fac
                                                                !  i=3 Gans Stokes fac
                                                                !  i=4 Gans Newton fac
                                                                !  i=5 slip adjustment to diameter
              ! Fall velocity model
              !   0 = No fall (just tracer)
              !   1 = Wilson and Huang
              !   2 = Wilson and Huang + Cunningham slip
              !   3 = Wilson and Huang + Mod by Pfeiffer Et al.
              !   4 = Ganser
              !   5 = Stokes flow for spherical particles + slip
      integer                                    :: FV_ID

#ifdef USEPOINTERS
      real(kind=sp),dimension(:,:,:,:),pointer :: vf_meso_last_step_MetP_sp => null()
      real(kind=sp),dimension(:,:,:,:),pointer :: vf_meso_next_step_MetP_sp => null()
#else
      real(kind=sp),dimension(:,:,:,:),allocatable :: vf_meso_last_step_MetP_sp
      real(kind=sp),dimension(:,:,:,:),allocatable :: vf_meso_next_step_MetP_sp
#endif

      real(kind=ip) :: phi_mean, phi_stddev

      contains

!******************************************************************************

      subroutine Allocate_Tephra

      implicit none

      allocate(Tephra_v_s(n_gs_max))
      allocate(Tephra_gsdiam(n_gs_max))
      allocate(Tephra_bin_mass(n_gs_max))
      allocate(Tephra_rho_m(n_gs_max))
      allocate(Tephra_gsF(n_gs_max))
      allocate(Tephra_gsF_fac(n_gs_max,5))

      end subroutine Allocate_Tephra

!******************************************************************************

      subroutine Allocate_Tephra_Met

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet

      implicit none

      write(global_production,*)"--------------------------------------------------"
      write(global_production,*)"---------- ALLOCATE_TEPHRA_MET -------------------"
      write(global_production,*)"--------------------------------------------------"

      allocate(vf_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet,n_gs_max))
      allocate(vf_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet,n_gs_max))

      end subroutine Allocate_Tephra_Met

!******************************************************************************

      subroutine Deallocate_Tephra

      implicit none

      if(allocated(Tephra_v_s))        deallocate(Tephra_v_s)
      if(allocated(Tephra_gsdiam))     deallocate(Tephra_gsdiam)
      if(allocated(Tephra_bin_mass))   deallocate(Tephra_bin_mass)
      if(allocated(Tephra_rho_m))      deallocate(Tephra_rho_m)
      if(allocated(Tephra_gsF))        deallocate(Tephra_gsF)
      if(allocated(Tephra_gsF_fac))    deallocate(Tephra_gsF_fac)

      end subroutine Deallocate_Tephra

!******************************************************************************

      subroutine Deallocate_Tephra_Met

      implicit none

#ifdef USEPOINTERS
      if(associated(vf_meso_last_step_MetP_sp)) deallocate(vf_meso_last_step_MetP_sp)
      if(associated(vf_meso_next_step_MetP_sp)) deallocate(vf_meso_next_step_MetP_sp)
#else
      if(allocated(vf_meso_last_step_MetP_sp)) deallocate(vf_meso_last_step_MetP_sp)
      if(allocated(vf_meso_next_step_MetP_sp)) deallocate(vf_meso_next_step_MetP_sp)
#endif

      end subroutine Deallocate_Tephra_Met

!******************************************************************************

      subroutine Set_Vf_Meso(Load_MesoSteps,Interval_Frac)

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax

      use solution,      only : &
         vf_pd,IsAloft

      use wind_grid,     only : &
         vf_meso_last_step_sp,vf_meso_next_step_sp

      use atmosphere,    only : &
         AirDens_meso_last_step_MetP_sp,AirDens_meso_next_step_MetP_sp,&
         AirVisc_meso_last_step_MetP_sp,AirVisc_meso_next_step_MetP_sp,&
         AirLamb_meso_last_step_MetP_sp,AirLamb_meso_next_step_MetP_sp

      use MetReader,     only : &
         nx_submet,ny_submet,np_fullmet,MR_dum3d_metP,MR_dum3d_compH,MR_iMetStep_Now,&
           MR_Regrid_MetP_to_CompGrid

      implicit none

      logical      ,intent(in) :: Load_MesoSteps
      real(kind=ip),intent(in) :: Interval_Frac

      real(kind=sp) :: vf_sp(n_gs_max)
      integer :: i,j,k,is,l
      real(kind=ip) :: dens,visc,lambda,Kna
      real(kind=ip) :: v_grav_set

      if(Load_MesoSteps)then
        do i=1,nx_submet
         do j=1,ny_submet
            do k=1,np_fullmet
               do is = 1,2  ! loop over the two braketing meso steps
                 if(is.eq.1)then
                   dens   = AirDens_meso_last_step_MetP_sp(i,j,k)
                   visc   = AirVisc_meso_last_step_MetP_sp(i,j,k)
                   lambda = AirLamb_meso_last_step_MetP_sp(i,j,k)
                 else
                   dens   = AirDens_meso_next_step_MetP_sp(i,j,k)
                   visc   = AirVisc_meso_next_step_MetP_sp(i,j,k)
                   lambda = AirLamb_meso_next_step_MetP_sp(i,j,k)
                 endif
                 select case (FV_ID)
                   ! Get settling velocity in m/s (will be converted to
                   ! km/hr and negated later)
                 case (0)  ! No fall
                   vf_sp(:) = 0.0_ip
                 case (1)  ! Wilson and Huang
                   do l=1,nsmax
                     if(.not.IsAloft(l)) cycle
                     v_grav_set = vset_WH(dens,Tephra_rho_m(l),visc, &
                                   Tephra_gsdiam(l),Tephra_gsF_fac(l,1),Tephra_gsF_fac(l,2))
                     vf_sp(l) = real(v_grav_set,kind=sp)
                   enddo
                 case (2)  ! Wilson and Huang + Cunningham slip
                   do l=1,nsmax
                     if(.not.IsAloft(l)) cycle
                     Kna = 2.0_ip*lambda/(Tephra_gsdiam(l)*Tephra_gsF_fac(l,5))
                     ! The non-continuum effect are > 1% when
                     ! gs<250*lam_col_windp(k)
                     if(Tephra_gsdiam(l).gt.LAM_GS_THRESH*lambda)then
                       v_grav_set = vset_WH(dens,Tephra_rho_m(l),visc, &
                                     Tephra_gsdiam(l),Tephra_gsF_fac(l,1),Tephra_gsF_fac(l,2))
                     else
                       v_grav_set = vset_WH_slip(dens,Tephra_rho_m(l),visc, &
                                     Tephra_gsdiam(l),Tephra_gsF_fac(l,1),Tephra_gsF_fac(l,2), &
                                     Kna)
                     endif
                     vf_sp(l) = real(v_grav_set,kind=sp)
                   enddo
                 case (3)  ! Wilson and Huang + Mod by Pfeiffer Et al.
                   do l=1,nsmax
                     if(.not.IsAloft(l)) cycle
                     v_grav_set = vset_WH_PCM(dens,Tephra_rho_m(l),visc, &
                                   Tephra_gsdiam(l),Tephra_gsF_fac(l,1),Tephra_gsF_fac(l,2))
                     vf_sp(l) = real(v_grav_set,kind=sp)
                   enddo
                 case (4)  ! Ganser
                   do l=1,nsmax
                     if(.not.IsAloft(l)) cycle
                     v_grav_set = vset_Gans(dens,Tephra_rho_m(l),&
                                            visc,Tephra_gsdiam(l), &
                                            Tephra_gsF_fac(l,3),Tephra_gsF_fac(l,4))
                     vf_sp(l) = real(v_grav_set,kind=sp)
                   enddo
                 case (5)  ! Stokes flow for spherical particles + slip
                   do l=1,nsmax
                     if(.not.IsAloft(l)) cycle
                     Kna = 2.0_ip*lambda/(Tephra_gsdiam(l))
                     v_grav_set = vset_Stokesslip(Tephra_rho_m(l),visc,Tephra_gsdiam(l),Kna)
                     vf_sp(l) = real(v_grav_set,kind=sp)
                   enddo
                 case default ! Wilson and Huang
                   do l=1,nsmax
                     if(.not.IsAloft(l)) cycle
                     v_grav_set = vset_WH(dens,Tephra_rho_m(l),visc, &
                                   Tephra_gsdiam(l),Tephra_gsF_fac(l,1),Tephra_gsF_fac(l,2))
                     vf_sp(l) = real(v_grav_set,kind=sp)
                   enddo
                 end select
                 if(is.eq.1)then
                   vf_meso_last_step_MetP_sp(i,j,k,:) = vf_sp
                 else
                   vf_meso_next_step_MetP_sp(i,j,k,:) = vf_sp
                 endif
              enddo !is
            enddo !k
          enddo !j
        enddo !i

        ! Now need to interpolate vs_meso_last_step_MetP_sp onto vf_meso_1_sp
        do l=1,nsmax
          if(.not.IsAloft(l)) cycle
          MR_dum3d_metP(:,:,:) = vf_meso_last_step_MetP_sp(:,:,:,l)
          call MR_Regrid_MetP_to_CompGrid(MR_iMetStep_Now)
          vf_meso_last_step_sp(:,:,:,l) = MR_dum3d_compH(:,:,:)*real(MPS_2_KMPHR,kind=sp)*(-1.0_sp)

          MR_dum3d_metP(:,:,:) = vf_meso_next_step_MetP_sp(:,:,:,l)
          call MR_Regrid_MetP_to_CompGrid(MR_iMetStep_Now+1)
          vf_meso_next_step_sp(:,:,:,l) = MR_dum3d_compH(:,:,:)*real(MPS_2_KMPHR,kind=sp)*(-1.0_sp)
        enddo
      endif

      vf_pd = 0.0_ip
      ! Now interpolate onto current time
      vf_pd(1:nxmax,1:nymax,1:nzmax,:) = real( vf_meso_last_step_sp(:,:,:,:),kind=ip) + &
                                         real((vf_meso_next_step_sp(:,:,:,:) - &
                                               vf_meso_last_step_sp(:,:,:,:)),kind=ip) * &
                                              Interval_Frac

      ! We only need actual fall velocity values in ghosts cells on the top and bottom
      vf_pd(-1:nxmax+2,-1:nymax+2,      -1,:) = vf_pd(-1:nxmax+2,-1:nymax+2,    1,:)
      vf_pd(-1:nxmax+2,-1:nymax+2,       0,:) = vf_pd(-1:nxmax+2,-1:nymax+2,    1,:)
      vf_pd(-1:nxmax+2,-1:nymax+2, nzmax+1,:) = vf_pd(-1:nxmax+2,-1:nymax+2,nzmax,:)
      vf_pd(-1:nxmax+2,-1:nymax+2, nzmax+2,:) = vf_pd(-1:nxmax+2,-1:nymax+2,nzmax,:)


      end subroutine Set_Vf_Meso

!******************************************************************************

      subroutine Calculate_Tephra_Shape!(ngs)

      implicit none

      !integer :: ngs

      integer :: i,j
      real(kind=ip) :: ellipse_area,ellipse_vol,equiv_rad
      real(kind=ip) :: p_exp,sphere_area,phi_sphere

      real(kind=ip) :: Dahneke_LD(13), Dahneke_RL(13),onF
      integer       :: Dahni

!     Calculating terms for Cunningham slip coefficients for non-spherical
!     particles
      Dahneke_LD(1)  = 1.0_ip
      Dahneke_LD(2)  = 2.0_ip
      Dahneke_LD(3)  = 3.0_ip
      Dahneke_LD(4)  = 4.0_ip
      Dahneke_LD(5)  = 6.0_ip
      Dahneke_LD(6)  = 8.0_ip
      Dahneke_LD(7)  = 10.0_ip
      Dahneke_LD(8)  = 15.0_ip
      Dahneke_LD(9)  = 20.0_ip
      Dahneke_LD(10) = 30.0_ip
      Dahneke_LD(11) = 50.0_ip
      Dahneke_LD(12) = 75.0_ip
      Dahneke_LD(13) = 100.0_ip
      Dahneke_RL(1)  = 1.0_ip
      Dahneke_RL(2)  = 1.23_ip
      Dahneke_RL(3)  = 1.39_ip
      Dahneke_RL(4)  = 1.51_ip
      Dahneke_RL(5)  = 1.72_ip
      Dahneke_RL(6)  = 1.87_ip
      Dahneke_RL(7)  = 2.0_ip
      Dahneke_RL(8)  = 2.25_ip
      Dahneke_RL(9)  = 2.43_ip
      Dahneke_RL(10) = 2.69_ip
      Dahneke_RL(11) = 3.02_ip
      Dahneke_RL(12) = 3.28_ip
      Dahneke_RL(13) = 3.47_ip

      do i=1,n_gs_max
        ! Precalculate the WilsonHuang factors based on particle shapes
        Tephra_gsF_fac(i,1) = Tephra_gsF(i)**(-0.828_ip)   ! WH Stokes factor
        Tephra_gsF_fac(i,2) = sqrt(1.07_ip-Tephra_gsF(i))  ! WH Newton factor
        ! and precalculate the K1 and K2 if we use the Ganser model
        ! First get sphericity.
        ! Note: sphericity is the ratio of the surface area of a sphere with
        !       equivalent volume to the actual surface area of the particle
        ! following http://en.wikipedia.org/wiki/Ellipsoid
        ellipse_vol  = 4.0_ip*PI*Tephra_gsF(i)*Tephra_gsF(i)/3.0_ip
        equiv_rad    = (ellipse_vol*3.0_ip/(4.0_ip*PI))**(1.0_ip/3.0_ip)
        sphere_area  = 4.0_ip*PI*equiv_rad*equiv_rad
        p_exp = 1.6075_ip
        ellipse_area = 4.0_ip*PI*(((Tephra_gsF(i)*Tephra_gsF(i))**p_exp + &
                            2.0_ip*(Tephra_gsF(i))**p_exp)/3.0_ip)**(1.0_ip/p_exp)
        phi_sphere = sphere_area/ellipse_area
        Tephra_gsF_fac(i,3) = 3.0_ip/(1.0_ip+2.0_ip*phi_sphere**(-0.5_ip))
        Tephra_gsF_fac(i,4) = 1.84148_ip*(-log10(phi_sphere))**0.5743_ip
        Tephra_gsF_fac(i,4) = 10.0_ip**Tephra_gsF_fac(i,4);
          ! Now calculate the aspherical correction to the Cunningham
          ! slip correction using Table 1 from Dahneke, Aerosol Sci, v4p163
          ! 1973.
        onF = 1.0_ip/Tephra_gsF(i)
        Dahni = 0
        do j = 1,12
          if(onF.ge.Dahneke_LD(j).and.onF.lt.Dahneke_LD(j+1)) Dahni = j
        enddo
        if(Dahni.eq.0)then
          write(global_info,*)"ERROR: could not determine slip correction factor."
          Tephra_gsF_fac(i,5) = 1.0_ip
        else
          Tephra_gsF_fac(i,5) = Dahneke_RL(Dahni) + (onF-Dahneke_LD(Dahni))* &
                         (Dahneke_RL(Dahni+1)-Dahneke_RL(Dahni))/     &
                          (Dahneke_LD(Dahni+1)-Dahneke_LD(Dahni))
        endif
      enddo

      end subroutine Calculate_Tephra_Shape

!******************************************************************************

      subroutine Sort_Tephra_Size

      implicit none

      integer :: i,j,l
      real(kind=ip) :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
      real(kind=ip) :: temp_a(5)
      real(kind=ip) :: viscnow, densnow, vfnow
      real(kind=ip), dimension(:), allocatable      :: vf_now

      densnow = 1.2_ip           !air density at STP (approximate)
      viscnow = 1.0e-05_ip       !air viscosity at STP (approximate)

      !Calculate fall velocity at 1 atmosphere and assume the relative values
      !are the same at higher elevation
      allocate(vf_now(n_gs_max))
      write(global_info,*) 'GSD before sorting:'
      do l=1,n_gs_max
         if (useCalcFallVel) then
            vfnow = vset_WH(densnow,Tephra_rho_m(l),viscnow, &
                                Tephra_gsdiam(l),Tephra_gsF_fac(l,1),Tephra_gsF_fac(l,2))
            vf_now(l) = vfnow
            write(global_info,*) 'l = ',l,', gsdiam = ',real(Tephra_gsdiam(l),kind=sp),&
                       ', rho_m = ',real(Tephra_rho_m(l),kind=sp),&
                       ', vf_now(l) = ',real(vf_now(l),kind=sp)
         else
            vf_now(l) = Tephra_v_s(l)
         endif
      enddo

        ! Using insertion sort on p321 of Numerical Recipes
      do j=2,n_gs_max
        tmp1   = vf_now(j)
        tmp2   = Tephra_gsdiam(j)
        tmp3   = Tephra_bin_mass(j)
        tmp4   = Tephra_rho_m(j)
        tmp5   = Tephra_v_s(j)
        tmp6   = Tephra_gsF(j)
        temp_a = Tephra_gsF_fac(j,:)
        do i=j-1,1,-1
            ! sort on grain-size
          if (vf_now(i).le.tmp1) goto 101
          vf_now(i+1)           = vf_now(i)
          Tephra_gsdiam(i+1)    = Tephra_gsdiam(i)
          Tephra_bin_mass(i+1)  = Tephra_bin_mass(i)
          Tephra_rho_m(i+1)     = Tephra_rho_m(i)
          Tephra_v_s(i+1)       = Tephra_v_s(i)
          Tephra_gsF(i+1)       = Tephra_gsF(i)
          Tephra_gsF_fac(i+1,:) = Tephra_gsF_fac(i,:)
        enddo
        i=0
 101    vf_now(i+1)          = tmp1
        Tephra_gsdiam(i+1)   = tmp2
        Tephra_bin_mass(i+1) = tmp3
        Tephra_rho_m(i+1)    = tmp4
        Tephra_v_s(i+1)      = tmp5
        Tephra_gsF(i+1)      = tmp6
        Tephra_gsF_fac(i+1,:)= temp_a
      enddo

      if (useVariableGSbins) call partition_gsbins(phi_mean,phi_stddev)

      end subroutine Sort_Tephra_Size

!******************************************************************************

      subroutine partition_gsbins(mu,sigma)

      implicit none

      integer :: i

      real(kind=ip) :: mu,sigma

      real(kind=ip),dimension(n_gs_max)   :: phi
      real(kind=ip),dimension(n_gs_max-1) :: phi_boundaries
      real(kind=ip),dimension(n_gs_max)   :: suppl_frac
      integer       :: mid_bin
      real(kind=ip) :: mid_bin_neg_half
      real(kind=ip) :: mid_bin_pos_half
      real(kind=ip) :: erf_at_a,erf_at_b
      real(kind=ip) :: fac1,frac_to_distrib

      suppl_frac = 0.0_ip
      mid_bin = 1
      fac1 = 0.5_ip

      !Note: This subroutine assumes that the grainsmax bins have been
      !sorted from smallest to largest

      ! We have n_gs_max bins with real grainsmax info.  The last bin read
      ! in (init_n_gs_max = n_gs_max+1) is
      ! the remaining fraction of the total volume that will be
      ! redistributed among the first n_gs_max bins with a normal
      ! distribution in phi.

      if(sigma.le.0.0_ip)then
        write(global_info,*)"ERROR: StdDev <= 0.0 for supplimental GS"
        write(global_info,*)"       distribution"
        stop 1
      endif

      ! Get the phi of the other bins
      do i = 1,n_gs_max
        phi(i) = -log(Tephra_gsdiam(i))/log(2.0_ip)
      enddo

      ! boundary between bins is the average phi
      do i = 1,n_gs_max-1
        phi_boundaries(i) = (phi(i)+phi(i+1))*0.5_ip
      enddo

      if (n_gs_max.eq.0)then
        write(global_info,*)"ERROR: Must have n_gs_max.ge.1"
        stop 1
      else if (n_gs_max.eq.1)then
        suppl_frac(1)=1.0_ip
      else
        ! Convert boundaries between grainsmaxs given to locations on a
        ! standard normal Gaussian, given mu and sigma
        do i = 1,n_gs_max-1
          phi_boundaries(i) = (phi_boundaries(i) - mu)/sigma
        enddo
        if(phi_boundaries(1).lt.0.0_ip)then
          write(global_info,*)"ERROR: Mean is not within GS-distribution."
            write(global_info,*)&
             "        Bin    :       gsdiam        :",&
             "        phi    :       right bin boundary (phi)"
          do i = 1,n_gs_max
            write(global_info,*)i,real(Tephra_gsdiam(i),kind=sp),real(phi(i),kind=sp),&
                        real(phi_boundaries(i),kind=sp)
          enddo
          write(global_info,*)"Mean = ",mu
          stop 1
        endif
        if(n_gs_max.eq.2)then
          ! We're dividing the surplus into two bins
           erf_at_a = 0.5_ip - fac1*erf(abs(phi_boundaries(1)))
           if(phi_boundaries(1).lt.0.0_ip)then
             suppl_frac(1) = erf_at_a
             suppl_frac(2) = 1.0_ip - erf_at_a
           else
             suppl_frac(2) = erf_at_a
             suppl_frac(1) = 1.0_ip - erf_at_a
           endif
        else if(n_gs_max.ge.3)then
          do i = 2,n_gs_max-1
            ! Find the bin that contains the mean of the distribution
            if(phi_boundaries(i).le.0.0_ip.and.phi_boundaries(i-1).ge.0.0_ip)then
              mid_bin = i
            endif
          enddo
        endif

        ! Integrating the positive half (lower index) of phi
        !   integrate from zero to the first boundary
        mid_bin_pos_half = fac1*erf(abs(phi_boundaries(mid_bin-1)))
        suppl_frac(mid_bin) = mid_bin_pos_half
        ! Now integrate the remaining bins on the positive side of mu
        do i = mid_bin-1,2,-1
          erf_at_a = fac1*erf(abs(phi_boundaries(i)))
          erf_at_b = fac1*erf(abs(phi_boundaries(i-1)))
          suppl_frac(i) = erf_at_b - erf_at_a
        enddo
        ! The last positive bin is from phi_boundaries(1) to Inf
        suppl_frac(1) = 0.5_ip - fac1*erf(abs(phi_boundaries(1)))

        !Now integrate the negative side
        mid_bin_neg_half = fac1*erf(abs(phi_boundaries(mid_bin)))
        suppl_frac(mid_bin) = mid_bin_neg_half
        ! Now integrate the remaining bins on the negative side of mu
        do i = mid_bin+1,n_gs_max-1
          erf_at_a = fac1*erf(abs(phi_boundaries(i-1)))
          erf_at_b = fac1*erf(abs(phi_boundaries(i)))
          suppl_frac(i) = erf_at_b - erf_at_a
        enddo
          ! The last negative bin is from phi_boundaries(1) to -Inf
        suppl_frac(n_gs_max) = 0.5_ip - fac1*erf(abs(phi_boundaries(n_gs_max-1)))
          ! Reconstruct middle bin from the two one-sided integrations
        suppl_frac(mid_bin) = mid_bin_pos_half + mid_bin_neg_half
      endif

      frac_to_distrib = 1.0_ip - sum(Tephra_bin_mass(1:n_gs_max))
      do i = 1,n_gs_max
        Tephra_bin_mass(i) = Tephra_bin_mass(i) + suppl_frac(i)*frac_to_distrib
      enddo

      end subroutine partition_gsbins

!******************************************************************************

      subroutine Collapse_GS

      use solution,      only : &
           mass_aloft,IsAloft

      implicit none

      integer :: n

      ! Loop through the grain sizes starting with the largest (fastest falling)
      ! and check if any have completely flushed out.  If so, then
      ! modify the max index of the grain-size loop (ns_aloft)
      ! The threshold for collapsing GS array is if there is less than a gram in
      ! that size bin aloft
      !do n = 1,ns_aloft
      do n = 1,n_gs_max
        if(IsAloft(n).and. &                     ! if bin is currently flagged as aloft
           mass_aloft(n).lt.AIRBORNE_THRESH)then ! but the mass is less than the thresh
          IsAloft(n) = .false.
          ns_aloft = ns_aloft-1
          write(global_info,*)"Grainsmax bin ",ns_aloft+1," has fully deposited."
          write(global_log ,*)"Grainsmax bin ",ns_aloft+1," has fully deposited."
        endif
      enddo

      end subroutine Collapse_GS

!******************************************************************************

      function vset_WH(rho_air,rho_m,eta,diam,Ffac1,Ffac2)

      implicit none

      real(kind=ip) :: vset_WH      ! Settling velocity in m/s
      real(kind=ip) :: rho_air   ! density of air in kg/m3
      real(kind=ip) :: rho_m     ! density of the particle in km/m3
      real(kind=ip) :: eta       ! dynamic viscosity of air in (kg/(m s))
      real(kind=ip) :: diam      ! diameter of the particle in m
      real(kind=ip) :: Ffac1     ! ! = F**(-0.828)
      real(kind=ip) :: Ffac2     ! = sqrt(1.07-F)

      real(kind=ip) :: vset2

      !Wilson and Huang (1979, EPSL v. 44, ppp. 311-324) give the drag coefficient of tephra 
      !as a function of its shape by defining tephra fragments as ellipsoidal with axes a, b, and c.
      !they also define F=(b+c)/2a, and give the drag coefficient as:
      !Cd = (24/Re)*F**(-0.828) + 2*sqrt(1.07-F)
      !They have a table in which a, b, and c were measured for pumice, glass, and feldspar fragments.
      !The average values of F for these three fragment types ranges from 0.39 to 0.42.  So I use 0.4
      !here.

      !F = 0.4
      !Ffac1 = 2.13547419198160     ! = F**(-0.828)  where F=0.4
      !Ffac1 = 1.34072838100802     ! = F**(-0.32)  where F=0.4; Suzuki's mod
      !Ffac2 = 0.818535277187245    ! = sqrt(1.07-F) where F=0.4

      !F = 0.6666
      !Ffac1 = 1.39896599613964      ! = F**(-0.828)  where F=0.6666
      !Ffac1 = 1.13854602830969      ! = F**(-0.32)  where F=0.6666; Suzuki's mod
      !Ffac2 = 0.635137780328017     ! = sqrt(1.07-F) where F=0.6666

      ! The following roots are solved by the maxima script:
      !  load(f90)$
      !  eq_1 : Vs*Vs = (4*diam*rho_m*g/(3*Cd*rho_air))$
      !  eq_2 : Cd = (24/Re)*Ffac1 + 2*Ffac2$
      !  eq_3 : Re = (rho_air*Vs*diam)/eta$
      !  eq_4:subst(eq_3,eq_2)$
      !  eq_5:subst(%,eq_1)$
      !  VsSoln : solve(%,Vs);
      !  f90(VsSoln[2]);    /* We want the second (positive) root */

      !vset1 =-(sqrt(6.)*sqrt(diam**3.*Ffac2*g*rho_air*rho_m + &
      !         54.*eta**2.*Ffac1**2.)+18.*eta*Ffac1)/(diam*Ffac2*rho_air)/3.0
      !vset2 = (sqrt(6.)*sqrt(diam**3.*Ffac2*g*rho_air*rho_m + &
      !         54.*eta**2.*Ffac1**2.)-18.*eta*Ffac1)/(diam*Ffac2*rho_air)/3.0

      ! Here's expanding out exponents and precalculating roots
      vset2 = (2.44948974278318_ip*sqrt(diam*diam*diam*Ffac2*GRAV*rho_air*rho_m + &
               54.0_ip*eta*eta*Ffac1*Ffac1)-18.0_ip*eta*Ffac1)/(diam*Ffac2*rho_air)/3.0_ip
      vset_WH = vset2

      end function vset_WH

!******************************************************************************

      function vset_WH_slip(rho_air,rho_m,eta,diam,Ffac1,Ffac2,Kna)
      ! Modification to Wilson and Huang's with Cunningham slip

      implicit none

      real(kind=ip) :: vset_WH_slip      ! Settling velocity in m/s
      real(kind=ip) :: rho_air   ! density of air in kg/m3
      real(kind=ip) :: rho_m     ! density of the particle in km/m3
      real(kind=ip) :: eta       ! dynamic viscosity of air in (kg/(m s))
      real(kind=ip) :: diam      ! diameter of the particle in m
      real(kind=ip) :: Ffac1     ! ! = F**(-0.828)
      real(kind=ip) :: Ffac2     ! = sqrt(1.07-F)
      real(kind=ip) :: Kna       ! adjusted Knudsen number 

      real(kind=ip) :: vnew, vold, Re                      ! old and new settling velocity
      real(kind=ip) :: Cd                                  ! drag coefficient
      real(kind=ip) :: Cslip                               ! Slip correction

      vold = 1.0_ip                               !assume an initial settling velocity of 1 m/s

      Re = rho_air*vold*diam/eta                       ! Eq. 4  of Wilson79

      ! Eq. 9.34 of Seinfeld and Pandis
      Cslip = 1.0_ip+ Kna * (1.257_ip + 0.4_ip*exp(-1.1_ip/Kna))

      ! Get initial Cd from Wilson79, Eq 12
      Cd = (24.0_ip/Re)*Ffac1 + Ffac2

      vnew = sqrt((4.0_ip*rho_m*diam*GRAV)/(3.0_ip*rho_air*Cd))   ! Eq. 15 of WoodsBursik91

      do while ((abs(vold-vnew)/vnew).gt.0.001_ip)
        vold = vnew
        Re = rho_air*vold*diam/eta
        Cd = (24.0_ip/Re)*Ffac1 + Ffac2
        Cd = Cd/Cslip
        vnew = sqrt((4.0_ip*rho_m*diam*GRAV)/(3.0_ip*rho_air*Cd))
      enddo

      vset_WH_slip = vnew

      return

      end function vset_WH_slip

!******************************************************************************

      function vset_WH_PCM(rho_air,rho_m,eta,diam,Ffac1,Ffac2)
      ! Modification to Wilson and Huang's model suggested by:
      ! T. Pfeiffer and A. Costa and G. Macedonio, JVGR, v140n4p273 2005
      ! DOI:10.1016/j.jvolgeores.2004.09.001

      implicit none

      real(kind=ip) :: vset_WH_PCM      ! Settling velocity in m/s
      real(kind=ip) :: rho_air   ! density of air in kg/m3
      real(kind=ip) :: rho_m     ! density of the particle in km/m3
      real(kind=ip) :: eta       ! dynamic viscosity of air in (kg/(m s))
      real(kind=ip) :: diam      ! diameter of the particle in m
      real(kind=ip) :: Ffac1     ! ! = F**(-0.828)
      real(kind=ip) :: Ffac2     ! = sqrt(1.07-F)

      real(kind=ip) :: vnew, vold, Re                      ! old and new settling velocity
      real(kind=ip) :: Cd                                  ! drag coefficient

      real(kind=ip) :: Cd100

      !Cd100 = (24./100.0)*F**(-0.828) + 2.*sqrt(1.07-F)   ! Cd at Re=100
      Cd100 = (24.0_ip/100.0_ip)*Ffac1 + Ffac2

      vold = 1.0_ip                               !assume an initial settling velocity of 1 m/s

      Re = rho_air*vold*diam/eta                       ! Eq. 4  of Wilson79

      ! Get initial Cd from Wilson79, Eq 12
      Cd = (24.0_ip/Re)*Ffac1 + Ffac2
      vnew = sqrt((4.0_ip*rho_m*diam*GRAV)/(3.0_ip*rho_air*Cd))   ! Eq. 15 of WoodsBursik91

      do while ((abs(vold-vnew)/vnew).gt.0.001_ip)
        vold = vnew
        Re = rho_air*vold*diam/eta
        if (Re.lt.100.0_ip)then
          Cd = (24.0_ip/Re)*Ffac1 + Ffac2
        elseif (Re.gt.1000.0_ip)then
          !  
          Cd = 1.0_ip
        else
          !  !interpolate between values at Re=1000 and 100
          Cd = 1.0_ip-(1.0_ip-Cd100)/900.0_ip *(1000.0_ip-Re)
        endif
        vnew = sqrt((4.0_ip*rho_m*diam*GRAV)/(3.0_ip*rho_air*Cd))
      enddo

      vset_WH_PCM = vnew

      return

      end function vset_WH_PCM

!******************************************************************************

      function vset_Gans(rho_air,rho_m,eta,diam,K1,K2)
      ! Fall velocity as calculated from
      ! Ganser, Powder Tech., v77,2p143, 1993
      ! DOI:10.1016/0032-5910(93)80051-B

      implicit none

      real(kind=ip) :: vset_Gans      ! Settling velocity in m/s
      real(kind=ip) :: rho_air   ! density of air in kg/m3
      real(kind=ip) :: rho_m     ! density of the particle in km/m3
      real(kind=ip) :: eta       ! dynamic viscosity of air in (kg/(m s))
      real(kind=ip) :: diam      ! diameter of the particle in m
      real(kind=ip) :: K1,K2     ! Stokes shape factor and Newtons shape factor

      real(kind=ip) :: vnew, vold, Re                      ! old and new settling velocity
      real(kind=ip) :: Cd                                  ! drag coefficient

      real(kind=ip) :: Cd1,Cd2

      vold = 1.0_ip                               !assume an initial settling velocity of 1 m/s

      Re = rho_air*vold*diam/eta                       ! Eq. 4  of Wilson79

      Cd1 = 1.0_ip + 0.1118_ip * (Re*K1*K2)**0.6567_ip
      Cd1 = Cd1 * 24.0_ip /(K1*Re)
      Cd2 = 0.4305_ip*K2/(1.0_ip+3305.0_ip/(Re*K1*K2))
      Cd = Cd1+Cd2;
      vnew = sqrt((4.0_ip*rho_m*diam*GRAV)/(3.0_ip*rho_air*Cd))   ! Eq. 15 of WoodsBursik91

      do while ((abs(vold-vnew)/vnew).gt.0.001_ip)
        vold = vnew
        Re = rho_air*vold*diam/eta
        Cd1 = 1.0_ip+0.1118_ip*(Re*K1*K2)**0.6567_ip
        Cd1 = Cd1 * 24.0_ip /(K1*Re)
        Cd2 = 0.4305_ip*K2/(1.0_ip+3305.0_ip/(Re*K1*K2))
        Cd = Cd1+Cd2;
        vnew = sqrt((4.0_ip*rho_m*diam*GRAV)/(3.0_ip*rho_air*Cd))
      enddo

      vset_Gans = vnew

      return

      end function vset_Gans

!******************************************************************************

      function vset_Stokesslip(rho_m,eta,diam,Kna)
      ! Fall velocity as calculated from
      ! Stokes flow plus Cunningham slip

      implicit none

      real(kind=ip) :: vset_Stokesslip      ! Settling velocity in m/s
      real(kind=ip) :: rho_m     ! density of the particle in km/m3
      real(kind=ip) :: eta       ! dynamic viscosity of air in (kg/(m s))
      real(kind=ip) :: diam      ! diameter of the particle in m
      real(kind=ip) :: Kna       ! adjusted Knudsen number

      real(kind=ip) :: Cslip                                  ! drag coefficient

      ! Eq. 9.34 of Seinfeld and Pandis
      Cslip = 1.0_ip+ Kna * (1.257_ip + 0.4_ip*exp(-1.1_ip/Kna))

      vset_Stokesslip = diam*diam*rho_m*GRAV*Cslip/(18.0_ip*eta)

      return

      end function vset_Stokesslip

!******************************************************************************

      function vset_Stokes_Cloud(rho_m,eta,diam)
      ! Fall velocity as calculated from Stokes flow
      ! for a sub-10 um particle acting as a cloud condensation nucleus
      ! to produce a 10 um ash/water particle

      implicit none

      real(kind=ip) vset_Stokes_Cloud  ! Settling velocity in m/s
      real(kind=ip) rho_m              ! density of the particle in kg/m3
      real(kind=ip) eta                ! dynamic viscosity of air in (kg/(m s))
      real(kind=ip) diam               ! diameter of the particle in m

      real(kind=ip) rho_cloud ! density of the cloud particle in kg/m3

      rho_cloud = (diam/1.0e-5_ip)**3.0_ip *(rho_m-1000.0_ip) + 1000.0_ip

      vset_Stokes_Cloud = 1.0_ip-10.0_ip*rho_cloud*GRAV/(18.0_ip*eta)

      return

      end function vset_Stokes_Cloud

!******************************************************************************

      end module Tephra

