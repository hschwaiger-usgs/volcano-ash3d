      module Output_Vars

      use precis_param

      use io_units

      use global_param,  only : &
         PI, M_2_MM, KM2_2_M2, KM3_2_M3, EPS_SMALL

      real(kind=ip)                :: DEPO_THRESH      = 1.0e-1_ip  !threshold deposit thickness (mm)
      real(kind=ip)                :: DEPRATE_THRESH   = 1.0e-2_ip  !threshold deposition rate (mm/hr)
      real(kind=ip)                :: CLOUDLOAD_THRESH = 2.0e-1_ip  !threshold cloud load (t/km2)
      real(kind=ip)                :: THICKNESS_THRESH = 1.0e-2_ip  !threshold thickness for start of deposition (mm)
      real(kind=ip)                :: TOTCON_THRESH    = 1.0e-3_ip  !threshold concentration (kg/km3)
      real(kind=ip)                :: DBZ_THRESH       =-2.0e+1_ip  !threshold dbZ

        ! These are the initialized values
      real(kind=ip)                :: DepositThickness_FillValue   =  0.0_ip
      real(kind=ip)                :: MaxConcentration_FillValue   =  0.0_ip
      real(kind=ip)                :: DepArrivalTime_FillValue     = -9999.0_ip
      real(kind=ip)                :: CloudArrivalTime_FillValue   = -9999.0_ip
      real(kind=ip)                :: CloudLoad_FillValue          = -9999.0_ip
      real(kind=ip)                :: MaxHeight_FillValue          = -9999.0_ip
      real(kind=ip)                :: MinHeight_FillValue          = -9999.0_ip
      real(kind=ip)                :: dbZCol_FillValue             = -100.0_ip

        ! Set this parameter if you want to include the standard derived variables
        ! in the output file (this should nearly always be true, downstream
        ! products rely on these variables)
      logical, parameter :: USE_OUTPROD_VARS  = .true.

        ! Set this parameter if you want to include velocities in the output file
      logical, parameter :: USE_WIND_VARS  = .false.

        ! Set this parameter to false if you do not want raw concentration values
        ! exported (only derived products and deposits)
      logical, parameter :: USE_RESTART_VARS  = .false.

        ! Set this to true if you want the extra output variables defined in the
        ! optional modules
      logical, parameter :: USE_OPTMOD_VARS  = .false.

        ! Set this parameter if you want to export additional variables
        ! to the netcdf file
      !logical, parameter :: USE_ADDITIONAL_VARS  = .false.
      !logical, parameter :: log_2d1_tmp          = .false.
      !logical, parameter :: log_3d1_tmp          = .false.
      !logical, parameter :: log_4d1_tmp          = .false.


      real(kind=ip) :: CloudArea                               ! area of ash cloud at a given time
      real(kind=ip) :: LoadVal(5), CloudLoadArea(5)            ! CloudLoadArea(i)=area of ash cloud exceeding LoadVal(i)
      real(kind=ip) :: AreaCovered                             ! area covered by ash deposit

        ! 1-D variables (in z)
      !real(kind=ip), dimension(:),allocatable   :: TotalConcentration ! concentration from all grain sizes in a vertical column

#ifdef USEPOINTERS
        ! 2-D variables (in x,y)
      real(kind=ip), dimension(:,:),pointer :: DepositThickness => null() ! accumulated ash thickness on ground in mm (x,y)
      real(kind=ip), dimension(:,:),pointer :: MaxConcentration => null() ! max concentration in the cloud at any i,j node
      real(kind=ip), dimension(:,:),pointer :: DepArrivalTime   => null() 
      real(kind=ip), dimension(:,:),pointer :: CloudArrivalTime => null() 
      real(kind=ip), dimension(:,:),pointer :: CloudLoad        => null() ! Ash load in cloud, tonnes/km2
      real(kind=ip), dimension(:,:),pointer :: CloudLoadLast    => null() ! Ash load at last time step, tonnes/km2
      real(kind=ip), dimension(:,:),pointer :: MaxHeight        => null() ! maximum cloud height
      real(kind=ip), dimension(:,:),pointer :: MinHeight        => null() ! cloud bottom height
      real(kind=ip), dimension(:,:),pointer :: dbZCol           => null() ! max reflectivity in a vertical column

        ! 3-D variables
        !   (in x,y,z)
      real(kind=ip), dimension(:,:,:),pointer :: dbZ => null()               ! radar reflectivty at time t (dbZ)
#else
        ! 2-D variables (in x,y)
      real(kind=ip), dimension(:,:),allocatable :: DepositThickness   ! accumulated ash thickness on ground in mm (x,y)
      real(kind=ip), dimension(:,:),allocatable :: MaxConcentration   ! max concentration in the cloud at any i,j node
      real(kind=ip), dimension(:,:),allocatable :: DepArrivalTime
      real(kind=ip), dimension(:,:),allocatable :: CloudArrivalTime
      real(kind=ip), dimension(:,:),allocatable :: CloudLoad          ! Ash load in cloud, tonnes/km2
      real(kind=ip), dimension(:,:),allocatable :: CloudLoadLast      ! Ash load at last time step, tonnes/km2
      real(kind=ip), dimension(:,:),allocatable :: MaxHeight          ! maximum cloud height
      real(kind=ip), dimension(:,:),allocatable :: MinHeight          ! cloud bottom height
      real(kind=ip), dimension(:,:),allocatable :: dbZCol             ! max reflectivity in a vertical column

        ! 3-D variables
        !   (in x,y,z)
      real(kind=ip), dimension(:,:,:),allocatable :: dbZ                ! radar reflectivty at time t (dbZ)
#endif

        ! User-defined 2-D static variables (in x,y)
      character(len=30), dimension(:),    allocatable :: var_User2d_static_XY_name
      character(len=30), dimension(:),    allocatable :: var_User2d_static_XY_unit
      character(len=30), dimension(:),    allocatable :: var_User2d_static_XY_lname
      real(kind=op),     dimension(:),    allocatable :: var_User2d_static_XY_MissVal
      real(kind=op),     dimension(:),    allocatable :: var_User2d_static_XY_FillVal
      real(kind=ip), dimension(:,:,:),    allocatable :: var_User2d_static_XY
        ! User-defined 2-D variables (in x,y)
      character(len=30), dimension(:),    allocatable :: var_User2d_XY_name
      character(len=30), dimension(:),    allocatable :: var_User2d_XY_unit
      character(len=30), dimension(:),    allocatable :: var_User2d_XY_lname
      real(kind=op),     dimension(:),    allocatable :: var_User2d_XY_MissVal
      real(kind=op),     dimension(:),    allocatable :: var_User2d_XY_FillVal
      real(kind=ip), dimension(:,:,:),    allocatable :: var_User2d_XY
        ! User-defined 3-D variables (in x,y,gs)
      character(len=30), dimension(:),    allocatable :: var_User3d_XYGs_name
      character(len=30), dimension(:),    allocatable :: var_User3d_XYGs_unit
      character(len=30), dimension(:),    allocatable :: var_User3d_XYGs_lname
      real(kind=op),     dimension(:),    allocatable :: var_User3d_XYGs_MissVal
      real(kind=op),     dimension(:),    allocatable :: var_User3d_XYGs_FillVal
      real(kind=ip), dimension(:,:,:,:),  allocatable :: var_User3d_XYGs
        ! User-defined 3-D variables (in x,y,z)
      character(len=30), dimension(:),    allocatable :: var_User3d_XYZ_name
      character(len=30), dimension(:),    allocatable :: var_User3d_XYZ_unit
      character(len=30), dimension(:),    allocatable :: var_User3d_XYZ_lname
      real(kind=op),     dimension(:),    allocatable :: var_User3d_XYZ_MissVal
      real(kind=op),     dimension(:),    allocatable :: var_User3d_XYZ_FillVal
      real(kind=ip), dimension(:,:,:,:),  allocatable :: var_User3d_XYZ
        ! User-defined 4-D variables (in x,y,z,gs)
      character(len=30), dimension(:),allocatable     :: var_User4d_XYZGs_name
      character(len=30), dimension(:),allocatable     :: var_User4d_XYZGs_unit
      character(len=30), dimension(:),    allocatable :: var_User4d_XYZGs_lname
      real(kind=op),     dimension(:),    allocatable :: var_User4d_XYZGs_MissVal
      real(kind=op),     dimension(:),    allocatable :: var_User4d_XYZGs_FillVal
      real(kind=ip), dimension(:,:,:,:,:),allocatable :: var_User4d_XYZGs

      contains

!******************************************************************************

!******************************************************************************

      subroutine Allocate_Output_Vars(nx,ny,nz)

      implicit none

      integer :: nx,ny,nz

      allocate(DepositThickness(nx,ny))
      DepositThickness = DepositThickness_FillValue
      allocate(MaxConcentration(nx,ny))
      MaxConcentration = MaxConcentration_FillValue
      allocate(DepArrivalTime(nx,ny))
      DepArrivalTime = DepArrivalTime_FillValue
      allocate(CloudArrivalTime(nx,ny))                    ! time of arrival of ash cloud
      CloudArrivalTime = CloudArrivalTime_FillValue
      allocate(CloudLoad(nx,ny))
      CloudLoad    = CloudLoad_FillValue
      allocate(CloudLoadLast(nx,ny))
      CloudLoadLast = CloudLoad_FillValue
      allocate(MaxHeight(nx,ny))                           ! maximum height (top) of ash cloud
      MaxHeight = MaxHeight_FillValue
      allocate(MinHeight(nx,ny))                           ! minimum height (bottom) of ash cloud
      MinHeight = MinHeight_FillValue
      allocate(dbZCol(nx,ny))                              ! reflectivity in a column of nodes
      dbZCol = dbZCol_FillValue

      allocate(dbZ(nx,ny,nz))                              ! radar reflectivity (dbZ)
      dbZ = dbZCol_FillValue

      end subroutine Allocate_Output_Vars

!******************************************************************************

!******************************************************************************

      subroutine Allocate_Output_UserVars(nx,ny,nz,ns)

      use io_data,       only : &
         nvar_User2d_static_XY,nvar_User2d_XY,&
         nvar_User3d_XYGs,nvar_User3d_XYZ,    &
         nvar_User4d_XYZGs

      implicit none

      integer :: nx,ny,nz,ns

        ! User-defined 2-D static variables (in x,y)
      write(global_info,*)"Allocating var_User2d_static_XY: ",nx,ny,nvar_User2d_static_XY
      allocate(var_User2d_static_XY_name(nvar_User2d_static_XY))
      allocate(var_User2d_static_XY_unit(nvar_User2d_static_XY))
      allocate(var_User2d_static_XY_lname(nvar_User2d_static_XY))
      allocate(var_User2d_static_XY_MissVal(nvar_User2d_static_XY))
      allocate(var_User2d_static_XY_FillVal(nvar_User2d_static_XY))
      allocate(var_User2d_static_XY(nx,ny,nvar_User2d_static_XY))
        ! User-defined 2-D variables (in x,y)
      write(global_info,*)"Allocating var_User2d_XY: ",nx,ny,nvar_User2d_XY
      allocate(var_User2d_XY_name(nvar_User2d_XY))
      allocate(var_User2d_XY_unit(nvar_User2d_XY))
      allocate(var_User2d_XY_lname(nvar_User2d_XY))
      allocate(var_User2d_XY_MissVal(nvar_User2d_XY))
      allocate(var_User2d_XY_FillVal(nvar_User2d_XY))
      allocate(var_User2d_XY(nx,ny,nvar_User2d_XY))
        ! User-defined 3-D variables (in x,y,gs)
      write(global_info,*)"Allocating var_User3d_XYGs: ",nx,ny,ns,nvar_User3d_XYGs
      allocate(var_User3d_XYGs_name(nvar_User3d_XYGs))
      allocate(var_User3d_XYGs_unit(nvar_User3d_XYGs))
      allocate(var_User3d_XYGs_lname(nvar_User3d_XYGs))
      allocate(var_User3d_XYGs_MissVal(nvar_User3d_XYGs))
      allocate(var_User3d_XYGs_FillVal(nvar_User3d_XYGs))
      allocate(var_User3d_XYGs(nx,ny,ns,nvar_User3d_XYGs))
        ! User-defined 3-D variables (in x,y,z)
      write(global_info,*)"Allocating var_User3d_XYZ: ",nx,ny,nz,nvar_User3d_XYZ
      allocate(var_User3d_XYZ_name(nvar_User3d_XYZ))
      allocate(var_User3d_XYZ_unit(nvar_User3d_XYZ))
      allocate(var_User3d_XYZ_lname(nvar_User3d_XYZ))
      allocate(var_User3d_XYZ_MissVal(nvar_User3d_XYZ))
      allocate(var_User3d_XYZ_FillVal(nvar_User3d_XYZ))
      allocate(var_User3d_XYZ(nx,ny,nz,nvar_User3d_XYZ))
        ! User-defined 4-D variables (in x,y,z,gs)
      write(global_info,*)"Allocating var_User4d_XYZGs: ",nx,ny,ns,nvar_User4d_XYZGs
      allocate(var_User4d_XYZGs_name(nvar_User4d_XYZGs))
      allocate(var_User4d_XYZGs_unit(nvar_User4d_XYZGs))
      allocate(var_User4d_XYZGs_lname(nvar_User4d_XYZGs))
      allocate(var_User4d_XYZGs_MissVal(nvar_User4d_XYZGs))
      allocate(var_User4d_XYZGs_FillVal(nvar_User4d_XYZGs))
      allocate(var_User4d_XYZGs(nx,ny,nz,ns,nvar_User4d_XYZGs))

      end subroutine Allocate_Output_UserVars

!******************************************************************************

      subroutine Deallocate_Output_Vars()

      implicit none

#ifndef USEPOINTERS
      deallocate(DepositThickness)
      deallocate(MaxConcentration)
      deallocate(DepArrivalTime)
      deallocate(CloudArrivalTime)
      deallocate(CloudLoad)
      deallocate(CloudLoadLast)
      deallocate(MaxHeight)
      deallocate(MinHeight)
      deallocate(dbZCol)
      deallocate(dbZ)
#endif

      end subroutine Deallocate_Output_Vars


!******************************************************************************

      subroutine Deallocate_Output_UserVars()

      implicit none

      deallocate(var_User2d_static_XY_name)
      deallocate(var_User2d_static_XY)
      deallocate(var_User2d_XY_name)
      deallocate(var_User2d_XY)
      deallocate(var_User3d_XYGs_name)
      deallocate(var_User3d_XYGs)
      deallocate(var_User3d_XYZ_name)
      deallocate(var_User3d_XYZ)
      deallocate(var_User4d_XYZGs_name)
      deallocate(var_User4d_XYZGs)

      end subroutine Deallocate_Output_UserVars

!******************************************************************************

      subroutine AshThicknessCalculator

      use Tephra,        only : &
         DepositDensity,n_gs_max

      use mesh,          only : &
         IsLatLon,nxmax,nymax,dx,dy,dz_vec_pd,kappa_pd

      use solution,      only : &
         DepositGranularity

      implicit none
      
      integer :: i,j
      real(kind=ip) :: dvol

      !calculate deposit thickness in mm, and area covered
      DepositThickness(:,:)     = DepositThickness_FillValue
      AreaCovered               = 0.0_ip
      dvol = dx*dy*dz_vec_pd(1)

      do i=1,nxmax
        do j=1,nymax
          DepositThickness(i,j) = sum(DepositGranularity(i,j,1:n_gs_max)) * &  ! in kg/km^3
                                  dz_vec_pd(1)                            / &  ! convert to kg/km^2
                                  KM2_2_M2                                / &  ! from kg/km^2 to kg/m^2
                                  DepositDensity                          * &  ! from kg/m^2 to m
                                  M_2_MM                                       ! from m to mm
          if(IsLatLon)then
            dvol = kappa_pd(i,j,0)
          endif
          if (DepositThickness(i,j).gt.DEPO_THRESH)then
            AreaCovered = AreaCovered + dvol/dz_vec_pd(1)
          endif
        enddo
      enddo

      end subroutine AshThicknessCalculator
      
!******************************************************************************

      subroutine dbZCalculator
      
      use Tephra,        only : &
         n_gs_max,Tephra_gsdiam,Tephra_rho_m

      use mesh,          only : &
         nxmax,nymax,nzmax,ts1

      use solution,      only : &
         concen_pd

      implicit none

      integer :: i,j,k,l
      real(kind=ip) :: NumDens          !number densities (#/m3) of particles
      real(kind=ip) :: zcol             !z value of cell

      call AshLoadCalculator         !calculate cloud load, T/km2

      dbZCol(:,:) = dbZCol_FillValue
      dbZ(:,:,:)  = dbZCol_FillValue

      !calculate particle collision rate between two particle sizes
      !note: this requires that only two particle sizes be used as input

      do i=1,nxmax
        do j=1,nymax
          !if (CloudLoad(i,j).lt.CLOUDLOAD_THRESH) cycle
          do k=1,nzmax
            zcol = 0.0_ip
            do l=1,n_gs_max
              !convert concentration (kg/km3) to number density (#/m3)
              NumDens = concen_pd(i,j,k,l,ts1)/(Tephra_rho_m(l)*PI* &
                            Tephra_gsdiam(l)**3.0_ip/6.0_ip)/KM3_2_M3     !particles/m3
              zcol    = zcol + NumDens*(1000.0_ip*Tephra_gsdiam(l))**6.0_ip
            enddo
            if(zcol.gt.EPS_SMALL)then
              dbZ(i,j,k) = 10.0_ip*log10(zcol)
            else
              dbZ(i,j,k) = dbZCol_FillValue
            endif
          enddo
          dbZCol(i,j) = maxval(dbZ(i,j,1:nzmax))
        enddo
      enddo

!      if (WriteKMLreflectivity) then
!        write(306,2) time, maxval(dbZ)
!2       format(f8.3,e14.4)
!      endif
    
      end subroutine dbZCalculator      

!******************************************************************************

      subroutine ConcentrationCalculator(nz)
      
      use mesh,          only : &
         IsLatLon,nxmax,nymax,nzmax,dx,dy,dz_vec_pd,kappa_pd,z_cc_pd,ts1

      use solution,      only : &
         concen_pd

      use Tephra,        only : &
         n_gs_max

      implicit none

      integer :: nz
      integer :: i,j,k,kk
      real(kind=ip),dimension(nz)  :: TotalConcentration
      real(kind=ip) :: dvol
 
      !calculate cloud concentration and cloud height

      CloudArea             = 0.0_ip
      MaxConcentration(1:nxmax,1:nymax) = MaxConcentration_FillValue
      MaxHeight(1:nxmax,1:nymax)        = MaxHeight_FillValue
      MinHeight(1:nxmax,1:nymax)        = MinHeight_FillValue

      do i=1,nxmax
        do j=1,nymax
          TotalConcentration = 0.0_ip
          do k=1,nzmax
            TotalConcentration(k) = sum(concen_pd(i,j,k,1:n_gs_max,ts1))
            if (TotalConcentration(k)>TOTCON_THRESH) then 
              if(IsLatLon)then
                !dvol = kappa(i,j,1)
                dvol = kappa_pd(i,j,k)
              else
                dvol = dx*dy*dz_vec_pd(k)
              endif
              !set height only if load>CLOUDLOAD_THRESH
              if (CloudLoad(i,j).gt.CLOUDLOAD_THRESH) &
                !MaxHeight(i,j)=z_cc_pd(k)+0.5_ip*dz
                MaxHeight(i,j)=z_cc_pd(k)+0.5_ip*dz_vec_pd(k)
              CloudArea = CloudArea + dvol/dz_vec_pd(k)
            endif
          enddo

          ! Now get cloud bottom
          if(MaxHeight(i,j).gt.MaxHeight_FillValue)then
            !kk = floor(MaxHeight(i,j)/dz)+1
            kk = nzmax
            do k=kk,1,-1
              if(TotalConcentration(k)>TOTCON_THRESH)then
                MinHeight(i,j)=z_cc_pd(k)-0.5_ip*dz_vec_pd(k)
              endif
            enddo
            !if cloud goes all the way to the ground, set min to 0.0
            if(MinHeight(i,j).lt.0.0_ip.and.MaxHeight(i,j).gt.0.0_ip)MinHeight(i,j)=0.0_ip
            ! Double-check that min doesn't exceed max
            MinHeight(i,j)=min(MaxHeight(i,j),MinHeight(i,j))

          endif

          MaxConcentration(i,j) = maxval(TotalConcentration)          !concentration in kg/km3
        enddo
      enddo
            
      end subroutine ConcentrationCalculator      

!******************************************************************************

      subroutine AshLoadCalculator
      
      use Tephra,        only : &
         n_gs_max

      use mesh,          only : &
         IsLatlon,nxmax,nymax,nzmax,dx,dy,dz_vec_pd,kappa_pd,ts1

      use solution,      only : &
         concen_pd

      implicit none

      integer :: i,j,k

      real(kind=ip) :: CellArea
 
      !calculate cloud concentration and cloud height
      LoadVal(1)    = 0.24_ip
      LoadVal(2)    = 1.0_ip
      LoadVal(3)    = 2.0_ip
      LoadVal(4)    = 4.0_ip
      LoadVal(5)    = 6.0_ip
      CloudLoadArea = 0.0_ip

      CloudLoad = 0.0_ip

      CellArea = dx*dy
      do i=1,nxmax
        do j=1,nymax
          if (IsLatlon) CellArea=kappa_pd(i,j,1)/dz_vec_pd(1)
          do k=1,nzmax
            ! Increment the cloud load for this cell
            CloudLoad(i,j) = CloudLoad(i,j) + &
                              sum(concen_pd(i,j,k,1:n_gs_max,ts1)) * & ! in kg/km^3
                                 dz_vec_pd(k)                      / & ! convert to kg/km^2
                               1.0e3_ip                            ! tonnes/km2
          enddo
          do k=1,5
            if (CloudLoad(i,j).gt.LoadVal(k))     CloudLoadArea(k) = CloudLoadArea(k) + CellArea
          enddo
        enddo
      enddo
            
      end subroutine AshLoadCalculator      

!******************************************************************************

      subroutine Gen_Output_Vars

      use mesh,          only : &
         nxmax, nymax, nzmax

      use time_data,     only : &
         time

      use io_data,       only : &
         Called_Gen_Output_Vars

      implicit none

      integer :: i,j

      real(kind=ip) :: thickness
      real(kind=ip) :: CloudLoadHere

      call AshThicknessCalculator

      call AshLoadCalculator

      call ConcentrationCalculator(nzmax)

      ! Mark the arrival time of any new deposit
      do i=1,nxmax
        do j=1,nymax
          thickness = DepositThickness(i,j)
          if(thickness.gt.DEPO_THRESH.and.DepArrivalTime(i,j).lt.0.0_ip)then
            DepArrivalTime(i,j)=time
          endif
          CloudLoadHere = CloudLoad(i,j)
          !0.2 T/km2 is roughly the detection limit of Pavolonis's SEVIRI
          !satellite retrievals
          if((CloudLoadHere.gt.CLOUDLOAD_THRESH).and.(CloudArrivalTime(i,j).lt.0.0_ip))then
            CloudArrivalTime(i,j)=time
          endif
        enddo
      enddo

      Called_Gen_Output_Vars = .true.

      end subroutine Gen_Output_Vars

!******************************************************************************

      subroutine Calc_AshVol_Aloft(vol)

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,kappa_pd,ts1

      use solution,      only : &
         concen_pd,mass_aloft

      use Tephra,        only : &
         MagmaDensity,n_gs_max

      implicit none

      real(kind=ip),intent(out) :: vol ! Total volume of ash still airborne
      integer :: isize

      vol = 0.0_ip

      ! First calculate mass of all species
      do isize=1,nsmax
        mass_aloft(isize) = sum(concen_pd(1:nxmax,1:nymax,1:nzmax,isize,ts1) *   & ! in kg/km^3
                                 kappa_pd(1:nxmax,1:nymax,1:nzmax))                ! convert to kg
      enddo

      ! Now loop over the tephra bins and calculate volume
      do isize=1,n_gs_max
          ! Increment total ash in air
        vol = vol + mass_aloft(isize)               /   & ! in kg
                    MagmaDensity                    /   & ! convert to m3
                    KM3_2_M3                              ! convert to km3
      enddo

      end subroutine Calc_AshVol_Aloft

!******************************************************************************

      subroutine Calc_AshVol_Deposit(vol)

      use mesh,          only : &
         nxmax,nymax,kappa_pd,ts1

      use solution,      only : &
         DepositGranularity

      use Tephra,        only : &
         MagmaDensity,n_gs_max

      implicit none

      real(kind=ip),intent(out) :: vol ! Total volume of ash in deposit
      integer :: isize

      vol = 0.0_ip

      do isize=1,n_gs_max
        vol = vol + sum(DepositGranularity(1:nxmax,1:nymax,isize) *   & ! in kg/km^3
                         kappa_pd(1:nxmax,1:nymax,0))          /   & ! convert to kg
                    MagmaDensity                      /   & ! convert to m3
                    KM3_2_M3                                ! convert to km3
      enddo
      end subroutine Calc_AshVol_Deposit

!******************************************************************************

      subroutine Calc_AshVol_Outflow(vol)

      use mesh,          only : &
         nxmax,nymax,nzmax,IsLatLon,kappa_pd,ts1

      use solution,      only : &
         outflow_xz1_pd,outflow_xz2_pd,outflow_yz1_pd,outflow_yz2_pd,outflow_xy2_pd

      use Tephra,        only : &
         MagmaDensity,n_gs_max

      implicit none

      real(kind=ip),intent(out) :: vol ! Total volume of ash flowing out of gird

      integer :: isize

      vol = 0.0_ip

      do isize=1,n_gs_max
        if(IsLatLon)then
          vol = vol + (                                &
                        sum(outflow_yz1_pd(        1:nymax,1:nzmax,isize)*   &
                                  kappa_pd(      1,1:nymax,1:nzmax)) +       &
                        sum(outflow_yz2_pd(        1:nymax,1:nzmax,isize)*   &
                                  kappa_pd(  nxmax,1:nymax,1:nzmax)) +       &
                        sum(outflow_xz1_pd(1:nxmax,        1:nzmax,isize)   *   &
                                  kappa_pd(1:nxmax,      1,1:nzmax)) +          &
                        sum(outflow_xz2_pd(1:nxmax,        1:nzmax,isize)  *   &
                                  kappa_pd(1:nxmax,  nymax,1:nzmax)) +         &
                        sum(outflow_xy2_pd(1:nxmax,1:nymax        ,isize)   *   &
                                  kappa_pd(1:nxmax,1:nymax,  nzmax)) )     /   & ! convert to kg
                        MagmaDensity                            /   & ! convert to m3
                        KM3_2_M3                                      ! convert to km3
        else
          vol = vol +     (                              &
                        sum(outflow_xz1_pd(1:nxmax,1:nzmax,isize))  +   &
                        sum(outflow_xz2_pd(1:nxmax,1:nzmax,isize))  +   &
                        sum(outflow_yz1_pd(1:nymax,1:nzmax,isize))  +   &
                        sum(outflow_yz2_pd(1:nymax,1:nzmax,isize))  +   &
                        sum(outflow_xy2_pd(1:nxmax,1:nymax,isize)) )/   & ! convert to kg
                        MagmaDensity                       /   & ! convert to m3
                        KM3_2_M3                                 ! convert to km3

        endif
      enddo

      end subroutine Calc_AshVol_Outflow


!******************************************************************************

      subroutine FirstAsh
      
!     Subroutine that determines whether the ash has yet hit any airports

      use mesh,          only : &
         IsLatLon
 
      use Airports,      only : &
         Airport_CloudHere,Airport_thickness,Airport_depRate,Airport_AshDuration,&
         Airport_thicknessLast,Airport_depRateLast,Airport_CloudHereLast,&
         Airport_CloudDuration,Airport_CloudArrived,Airport_CloudArrivalTime,&
         Airport_AshArrived,Airport_AshArrivalTime,nairports,Airport_i,Airport_j,&
           bilinear_thickness

      use time_data,     only : &
         time, dt
      
      integer :: i

      ! Requires that AshThicknessCalculator and AshLoadCalculator were called

        !record thickness & cloud concentration in last time step      
      Airport_CloudHereLast  = Airport_CloudHere
      Airport_ThicknessLast  = Airport_thickness
      Airport_depRateLast    = Airport_depRate
      do i=1,nairports
        Airport_CloudHere(i)  = CloudLoad(Airport_i(i),Airport_j(i))
        Airport_thickness(i)  = bilinear_thickness(i,DepositThickness)
        Airport_depRate(i)    = (Airport_thickness(i)-Airport_thicknessLast(i))/dt       !dep. rate, mm/hr

        !For airports where ash has arrived . . .
        !mark fall duration if ash has stopped falling
        if (Airport_AshArrived(i).eqv..true.) then
          if ((Airport_depRate(i).lt.DEPRATE_THRESH).and.(Airport_depRateLast(i).ge.DEPRATE_THRESH)) then
            Airport_AshDuration(i) = time-Airport_AshArrivalTime(i)
          endif
        endif

        !mark cloud duration is cloud has passed
        if (Airport_CloudArrived(i).eqv..true.) then
          if ((Airport_CloudHere(i).le.CLOUDLOAD_THRESH).and.(Airport_CloudHereLast(i).gt.CLOUDLOAD_THRESH)) then
            Airport_CloudDuration(i) = time-Airport_CloudArrivalTime(i)
          endif
        endif

        !For airports where ash has not yet arrived . . .
        !if ash load>0.1 T/km2 , call it "arrived"
        if ((Airport_CloudArrived(i).eqv..false.).and.(Airport_CloudHere(i).gt.CLOUDLOAD_THRESH)) then
          Airport_CloudArrived(i) = .true.
          Airport_CloudArrivalTime(i) = time
        endif
        !if ash thickness>0.1 mm, call it "arrived"
        if ((Airport_AshArrived(i).eqv..false.).and.(Airport_thickness(i).gt.THICKNESS_THRESH)) then
          Airport_AshArrived(i) = .true.
          Airport_AshArrivalTime(i) = time
        endif
      enddo

      return

      end subroutine FirstAsh

      end module Output_Vars

