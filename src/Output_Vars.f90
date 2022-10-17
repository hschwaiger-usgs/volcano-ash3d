      module Output_Vars

      use precis_param

      use io_units

      use global_param,  only : &
         PI, M_2_MM, KM2_2_M2, KG_2_MG, KM3_2_M3, EPS_SMALL

      use time_data,     only : &
         time_native

        ! PRODUCT ID: KML var ID:: variable                     : units  : dims        : Output kml file name
        !----------------------------------------------------------------------------------------------------------
        ! iprod =  1:ivar =     :: full concentration array     : kg/km3 : (x,y,z,t,g) :
        ! iprod =  2:ivar =     :: deposit granularity          : kg/m2  : (x,y,t,g)   :
        ! iprod =  3:ivar =   7 :: deposit                      : mm     : (x,y,t)     : deposit_thickness_mm.kml
        ! iprod =  4:ivar =   8 :: deposit (NWS)                : inches : (x,y,t)     : deposit_thickness_inches.kml
        ! iprod =  5:ivar =  *7 :: deposit final                : mm     : (x,y)       :
        ! iprod =  6:ivar =  *8 :: deposit final (NWS)          : inches : (x,y)       :
        ! iprod =  7:ivar =   9 :: ashfall arrival time         : hours  : (x,y)       : ashfall_arrivaltimes_hours.kml
        ! iprod =  8:ivar =     :: ashfall at airports          : mm     : (pt,tn)     : ash_arrivaltimes_airports.kml

        ! iprod =  9:ivar =   1 :: ash-cloud concentration      : mg/m3  : (x,y,t)     : CloudConcentration.kml
        ! iprod = 10:ivar =   2 :: ash-cloud height (top)       : km     : (x,y,t)     : CloudHeight.kml
        ! iprod = 11:ivar =   3 :: ash-cloud height (bot)       : km     : (x,y,t)     : CloudBottom.kml
        ! iprod = 12:ivar =   4 :: ash-cloud load               : T/km2  : (x,y,t)     : CloudLoad.kml
        ! iprod = 13:ivar =   6 :: ash-cloud radar reflectivity : dBz    : (x,y,t)     : reflectivity.kml
        ! iprod = 14:ivar =   5 :: ash-cloud arrival time       : hours  : (x,y)       : cloud_arrivaltimes_hours.kml

        ! iprod = 15:ivar =  10 :: topography                   : km     : (x,y)       : Topography.kml
        ! iprod = 16:ivar =     :: profiles concentration?      :

      real(kind=ip)                :: DEPO_THRESH           = 1.0e-2_ip  ! threshold deposit thickness (mm)
      real(kind=ip)                :: DEPRATE_THRESH        = 1.0e-2_ip  ! threshold deposition rate (mm/hr)
      real(kind=ip)                :: CLOUDCON_THRESH       = 1.0e-3_ip  ! threshold cloud concentration (kg/km3) for output
      real(kind=ip)                :: CLOUDCON_GRID_THRESH  = 1.0e-7_ip  ! threshold cloud concentration (kg/km3) for subgrid
      real(kind=ip)                :: CLOUDLOAD_THRESH      = 2.0e-1_ip  ! threshold cloud load (t/km2)
                                       ! 0.2 T/km2 is roughly the detection
                                       ! limit of Pavolonis's SEVIRI satellite retrievals

      real(kind=ip)                :: THICKNESS_THRESH = 1.0e-2_ip  !threshold thickness for start of deposition (mm)
      real(kind=ip)                :: DBZ_THRESH       =-2.0e+1_ip  !threshold dbZ

      ! These are the initialized values
      real(kind=op)                :: DepositThickness_FillValue   =     0.0_op
      real(kind=op)                :: MaxConcentration_FillValue   = -9999.0_op
      real(kind=op)                :: DepArrivalTime_FillValue     = -9999.0_op
      real(kind=op)                :: CloudArrivalTime_FillValue   = -9999.0_op
      real(kind=op)                :: CloudLoad_FillValue          = -9999.0_op
      real(kind=op)                :: MaxHeight_FillValue          = -9999.0_op
      real(kind=op)                :: MinHeight_FillValue          = -9999.0_op
      real(kind=op)                :: dbZCol_FillValue             =  -100.0_op
      real(kind=op)                :: pr_ash_FillValue             =     0.0_op

      logical            :: Calculated_Cloud_Load
      logical            :: Calculated_AshThickness
        ! Set this parameter if you want to include velocities in the output file
      logical, parameter :: USE_WIND_VARS  = .false.
      !logical, parameter :: USE_WIND_VARS  = .true.

        ! Set this to true if you want the extra output variables defined in the
        ! optional modules
      logical, parameter :: USE_OPTMOD_VARS  = .true.

        ! Set this parameter if you want to export additional variables
        ! to the netcdf file
      !logical, parameter :: USE_ADDITIONAL_VARS  = .false.
      !logical, parameter :: log_2d1_tmp          = .false.
      !logical, parameter :: log_3d1_tmp          = .false.
      !logical, parameter :: log_4d1_tmp          = .false.

        ! This variable should be set to true if you want to include
        ! the standard derived variables
        ! in the output file (this should nearly always be true, downstream
        ! products rely on these variables)
      logical :: USE_OUTPROD_VARS  = .true.

        ! Set this variable to false if you do not want raw concentration
        ! values exported (only derived products and deposits)
      logical :: USE_RESTART_VARS  = .true.

      real(kind=ip) :: CloudArea                ! area of ash cloud at a given time
      real(kind=ip) :: LoadVal(5)               ! 5 threshold values for area calculations
      real(kind=ip) :: CloudLoadArea(5)         ! Corresponding areas where ash cloud exceeds LoadVal(i)
      real(kind=ip) :: DepositAreaCovered       ! area covered by ash deposit

#ifdef USEPOINTERS
        ! 2-D variables (in x,y)
      logical, dimension(:,:),pointer :: Mask_Cloud
      logical, dimension(:,:),pointer :: Mask_Deposit
      real(kind=ip), dimension(:,:),pointer :: DepositThickness => null() ! accumulated ash thickness on ground (mm)
      real(kind=ip), dimension(:,:),pointer :: MaxConcentration => null() ! max concentration in the cloud at any i,j node (mg/m3)
      real(kind=ip), dimension(:,:),pointer :: DepArrivalTime   => null() ! (hours)
      real(kind=ip), dimension(:,:),pointer :: CloudArrivalTime => null() ! (hours)
      real(kind=ip), dimension(:,:),pointer :: CloudLoad        => null() ! Ash load in cloud, (tonnes/km2)
      real(kind=ip), dimension(:,:),pointer :: CloudLoadLast    => null() ! Ash load at last time step, (tonnes/km2)
      real(kind=ip), dimension(:,:),pointer :: MaxHeight        => null() ! maximum cloud height (km)
      real(kind=ip), dimension(:,:),pointer :: MinHeight        => null() ! cloud bottom height (km)
      real(kind=ip), dimension(:,:),pointer :: dbZCol           => null() ! max reflectivity in a vertical column (dB)

      real(kind=ip), dimension(:,:,:),pointer :: pr_ash         => null() ! concentration profile

        ! 3-D variables
        !   (in x,y,z)
      real(kind=ip), dimension(:,:,:),pointer :: dbZ => null()               ! radar reflectivty at time t (dbZ)
#else
        ! 2-D variables (in x,y)
      logical, dimension(:,:),allocatable :: Mask_Cloud
      logical, dimension(:,:),allocatable :: Mask_Deposit
      real(kind=ip), dimension(:,:),allocatable :: DepositThickness   ! accumulated ash thickness on ground in mm (x,y)
      real(kind=ip), dimension(:,:),allocatable :: MaxConcentration   ! max concentration in the cloud at any i,j node (mg/m3)
      real(kind=ip), dimension(:,:),allocatable :: DepArrivalTime     ! (hours)
      real(kind=ip), dimension(:,:),allocatable :: CloudArrivalTime   ! (hours)
      real(kind=ip), dimension(:,:),allocatable :: CloudLoad          ! Ash load in cloud, tonnes/km2
      real(kind=ip), dimension(:,:),allocatable :: CloudLoadLast      ! Ash load at last time step, tonnes/km2
      real(kind=ip), dimension(:,:),allocatable :: MaxHeight          ! maximum cloud height
      real(kind=ip), dimension(:,:),allocatable :: MinHeight          ! cloud bottom height
      real(kind=ip), dimension(:,:),allocatable :: dbZCol             ! max reflectivity in a vertical column

      real(kind=ip), dimension(:,:,:),allocatable :: pr_ash           ! concentration profile

        ! These arrays are only used when reading an output file of unknown size
      real(kind=ip), dimension(:,:),allocatable :: R_XY
      integer       :: R_nx,R_ny
      real(kind=ip) :: R_xll,R_yll
      real(kind=ip) :: R_dx,R_dy
      real(kind=ip) :: R_Fill

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
      real(kind=op), dimension(:,:,:),    allocatable :: var_User2d_static_XY
        ! User-defined 2-D variables (in x,y)
      character(len=30), dimension(:),    allocatable :: var_User2d_XY_name
      character(len=30), dimension(:),    allocatable :: var_User2d_XY_unit
      character(len=30), dimension(:),    allocatable :: var_User2d_XY_lname
      real(kind=op),     dimension(:),    allocatable :: var_User2d_XY_MissVal
      real(kind=op),     dimension(:),    allocatable :: var_User2d_XY_FillVal
      real(kind=op), dimension(:,:,:),    allocatable :: var_User2d_XY
        ! User-defined 3-D variables (in x,y,gs)
      character(len=30), dimension(:),    allocatable :: var_User3d_XYGs_name
      character(len=30), dimension(:),    allocatable :: var_User3d_XYGs_unit
      character(len=30), dimension(:),    allocatable :: var_User3d_XYGs_lname
      real(kind=op),     dimension(:),    allocatable :: var_User3d_XYGs_MissVal
      real(kind=op),     dimension(:),    allocatable :: var_User3d_XYGs_FillVal
      real(kind=op), dimension(:,:,:,:),  allocatable :: var_User3d_XYGs
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
      real(kind=op), dimension(:,:,:,:,:),allocatable :: var_User4d_XYZGs

      contains

!******************************************************************************

!******************************************************************************

      subroutine Allocate_Output_Vars(nx,ny,nz)

      implicit none

      integer :: nx,ny,nz

      allocate(Mask_Cloud(nx,ny))
      Mask_Cloud = .false.
      allocate(Mask_Deposit(nx,ny))
      Mask_Deposit = .false.
      allocate(DepositThickness(nx,ny))
      DepositThickness = 0.0_ip
      allocate(MaxConcentration(nx,ny))
      MaxConcentration = 0.0_ip
      allocate(DepArrivalTime(nx,ny))
      DepArrivalTime = DepArrivalTime_FillValue
      allocate(CloudArrivalTime(nx,ny))                    ! time of arrival of ash cloud
      CloudArrivalTime = CloudArrivalTime_FillValue
      allocate(CloudLoad(nx,ny))
      CloudLoad    = 0.0_ip
      allocate(CloudLoadLast(nx,ny))
      CloudLoadLast = 0.0_ip
      allocate(MaxHeight(nx,ny))                           ! maximum height (top) of ash cloud
      MaxHeight = 0.0_ip
      allocate(MinHeight(nx,ny))                           ! minimum height (bottom) of ash cloud
      MinHeight = 0.0_ip
      allocate(dbZCol(nx,ny))                              ! reflectivity in a column of nodes
      dbZCol = 0.0_ip

      allocate(dbZ(nx,ny,nz))                              ! radar reflectivity (dbZ)
      dbZ = 0.0_ip

      end subroutine Allocate_Output_Vars

!******************************************************************************

!******************************************************************************

      subroutine Allocate_Profile(nz,nt,nv)

      implicit none

      integer :: nz,nt,nv

      if(nv.gt.0)then
        allocate(time_native(nt))
        allocate(pr_ash(nz,nt,nv))                       ! vertical ash profile
        pr_ash = 0.0_op
      endif

      end subroutine Allocate_Profile

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
      if(allocated(Mask_Cloud))       deallocate(Mask_Cloud)
      if(allocated(Mask_Deposit))     deallocate(Mask_Deposit)
      if(allocated(DepositThickness)) deallocate(DepositThickness)
      if(allocated(MaxConcentration)) deallocate(MaxConcentration)
      if(allocated(DepArrivalTime))   deallocate(DepArrivalTime)
      if(allocated(CloudArrivalTime)) deallocate(CloudArrivalTime)
      if(allocated(CloudLoad))        deallocate(CloudLoad)
      if(allocated(CloudLoadLast))    deallocate(CloudLoadLast)
      if(allocated(MaxHeight))        deallocate(MaxHeight)
      if(allocated(MinHeight))        deallocate(MinHeight)
      if(allocated(dbZCol))           deallocate(dbZCol)
      if(allocated(dbZ))              deallocate(dbZ)

      if(allocated(pr_ash))           deallocate(pr_ash)
#endif

      end subroutine Deallocate_Output_Vars


!******************************************************************************

      subroutine Deallocate_Output_UserVars()

      implicit none

      if(allocated(var_User2d_static_XY_name))    deallocate(var_User2d_static_XY_name)
      if(allocated(var_User2d_static_XY_unit))    deallocate(var_User2d_static_XY_unit)
      if(allocated(var_User2d_static_XY_lname))   deallocate(var_User2d_static_XY_lname)
      if(allocated(var_User2d_static_XY_MissVal)) deallocate(var_User2d_static_XY_MissVal)
      if(allocated(var_User2d_static_XY_FillVal)) deallocate(var_User2d_static_XY_FillVal)
      if(allocated(var_User2d_static_XY))         deallocate(var_User2d_static_XY)
      if(allocated(var_User2d_XY_name))           deallocate(var_User2d_XY_name)
      if(allocated(var_User2d_XY_unit))           deallocate(var_User2d_XY_unit)
      if(allocated(var_User2d_XY_lname))          deallocate(var_User2d_XY_lname)
      if(allocated(var_User2d_XY_MissVal))        deallocate(var_User2d_XY_MissVal)
      if(allocated(var_User2d_XY_FillVal))        deallocate(var_User2d_XY_FillVal)
      if(allocated(var_User2d_XY))                deallocate(var_User2d_XY)
      if(allocated(var_User3d_XYGs_name))         deallocate(var_User3d_XYGs_name)
      if(allocated(var_User3d_XYGs_unit))         deallocate(var_User3d_XYGs_unit)
      if(allocated(var_User3d_XYGs_lname))        deallocate(var_User3d_XYGs_lname)
      if(allocated(var_User3d_XYGs_MissVal))      deallocate(var_User3d_XYGs_MissVal)
      if(allocated(var_User3d_XYGs_FillVal))      deallocate(var_User3d_XYGs_FillVal)
      if(allocated(var_User3d_XYGs))              deallocate(var_User3d_XYGs)
      if(allocated(var_User3d_XYZ_name))          deallocate(var_User3d_XYZ_name)
      if(allocated(var_User3d_XYZ_unit))          deallocate(var_User3d_XYZ_unit)
      if(allocated(var_User3d_XYZ_lname))         deallocate(var_User3d_XYZ_lname)
      if(allocated(var_User3d_XYZ_MissVal))       deallocate(var_User3d_XYZ_MissVal)
      if(allocated(var_User3d_XYZ_FillVal))       deallocate(var_User3d_XYZ_FillVal)
      if(allocated(var_User3d_XYZ))               deallocate(var_User3d_XYZ)
      if(allocated(var_User4d_XYZGs_name))        deallocate(var_User4d_XYZGs_name)
      if(allocated(var_User4d_XYZGs_unit))        deallocate(var_User4d_XYZGs_unit)
      if(allocated(var_User4d_XYZGs_lname))       deallocate(var_User4d_XYZGs_lname)
      if(allocated(var_User4d_XYZGs_MissVal))     deallocate(var_User4d_XYZGs_MissVal)
      if(allocated(var_User4d_XYZGs_FillVal))     deallocate(var_User4d_XYZGs_FillVal)
      if(allocated(var_User4d_XYZGs))             deallocate(var_User4d_XYZGs)

      end subroutine Deallocate_Output_UserVars

!******************************************************************************

      subroutine AshThicknessCalculator

      use Tephra,        only : &
         DepositDensity,n_gs_max

      use mesh,          only : &
         nxmax,nymax,dz_vec_pd,sigma_nz_pd

      use solution,      only : &
         DepositGranularity

      implicit none
      
      integer :: i,j

      !calculate deposit thickness in mm, and area covered
      Mask_Deposit(:,:)         = .false.
      DepositThickness(:,:)     = 0.0_ip
      DepositAreaCovered        = 0.0_ip

      if(n_gs_max.gt.0)then
        do i=1,nxmax
          do j=1,nymax
            DepositThickness(i,j) = sum(DepositGranularity(i,j,1:n_gs_max)) * &  ! in kg/km^3
                                    dz_vec_pd(1)                            / &  ! convert to kg/km^2
                                    KM2_2_M2                                / &  ! from kg/km^2 to kg/m^2
                                    DepositDensity                          * &  ! from kg/m^2 to m
                                    M_2_MM                                       ! from m to mm
            if (DepositThickness(i,j).gt.DEPO_THRESH)then
              Mask_Deposit(i,j) = .true.
              DepositAreaCovered = DepositAreaCovered + sigma_nz_pd(i,j,1)
            endif
          enddo
        enddo
      endif

      Calculated_AshThickness = .true.

      end subroutine AshThicknessCalculator
      
!******************************************************************************

      subroutine dbZCalculator
      
      use Tephra,        only : &
         n_gs_max,Tephra_gsdiam,Tephra_rho_m

      use mesh,          only : &
         nzmax,ts1

      use solution,      only : &
         concen_pd,imin,imax,jmin,jmax,kmin,kmax

      implicit none

      integer :: i,j,k,l
      real(kind=ip) :: NumDens          !number densities (#/m3) of particles
      real(kind=ip) :: zcol             !z value of cell

      dbZCol(:,:) = 0.0_op
      dbZ(:,:,:)  = 0.0_op

      !calculate particle collision rate between two particle sizes
      !note: this requires that only two particle sizes be used as input

      if(n_gs_max.gt.0)then
        do i=imin,imax
          do j=jmin,jmax
            !if (CloudLoad(i,j).lt.CLOUDLOAD_THRESH) cycle
            do k=kmin,kmax
              zcol = 0.0_ip
              do l=1,n_gs_max
                !convert concentration (kg/km3) to number density (#/m3)
                NumDens = concen_pd(i,j,k,l,ts1) / &
                            (Tephra_rho_m(l)*PI*Tephra_gsdiam(l)**3.0_ip/6.0_ip) / &
                            KM3_2_M3                                  !particles/m3
                zcol    = zcol + NumDens*(1000.0_ip*Tephra_gsdiam(l))**6.0_ip
              enddo
              if(zcol.gt.EPS_SMALL)then
                dbZ(i,j,k) = 10.0_ip*log10(zcol)
              else
                dbZ(i,j,k) = 0.0_op
              endif
            enddo
            dbZCol(i,j) = maxval(dbZ(i,j,1:nzmax))
          enddo
        enddo
      endif

      end subroutine dbZCalculator      

!******************************************************************************

      subroutine ConcentrationCalculator

      use mesh,          only : &
         nxmax,nymax,nzmax,dz_vec_pd,z_cc_pd,ts1,sigma_nz_pd

      use solution,      only : &
         concen_pd,imin,imax,jmin,jmax

      use Tephra,        only : &
         n_gs_max

      implicit none

      integer :: i,j,k
      real(kind=ip) :: CellArea

      ! Both these concentration variables are in the 'natural' units of kg/km3
      real(kind=ip),dimension(nzmax) :: TotalConcentration ! concentration from all grain sizes as a vertical column
      real(kind=ip)                  :: MaxTotalConcentration
 
      ! calculate cloud concentration,
      ! cloud height, and cloud bottom

      CloudArea                         = 0.0_ip
      CloudLoad(1:nxmax,1:nymax)        = 0.0_ip
      MaxConcentration(1:nxmax,1:nymax) = 0.0_ip
      MaxHeight(1:nxmax,1:nymax)        = 0.0_ip
      MinHeight(1:nxmax,1:nymax)        = 100.0_ip

      Mask_Cloud(:,:) = .false.

      if(n_gs_max.gt.0)then
        do i=imin,imax
          do j=jmin,jmax
            CellArea = sigma_nz_pd(i,j,1)
            TotalConcentration = 0.0_ip
            do k=1,nzmax
               ! Increment the cloud load for this column
              CloudLoad(i,j) = CloudLoad(i,j) + &
                                sum(concen_pd(i,j,k,1:n_gs_max,ts1)) * & ! in kg/km^3
                                   dz_vec_pd(k)                      / & ! convert to kg/km^2
                                1.0e3_ip                                 ! tonnes/km2
               ! Pull out column total concentration
              TotalConcentration(k) = sum(concen_pd(i,j,k,1:n_gs_max,ts1))

               ! Note max height
              if(TotalConcentration(k).ge.CLOUDCON_THRESH)then
                ! this is overwritten if the k+1 is also identified as an ash cloud
                MaxHeight(i,j)=z_cc_pd(k)+0.5_ip*dz_vec_pd(k)
              endif
            enddo

             ! Now go from the top down and find the bottom
            do k=nzmax,1,-1
             if(TotalConcentration(k).le.CLOUDCON_THRESH)then
                MinHeight(i,j)=z_cc_pd(k)-0.5_ip*dz_vec_pd(k)
              endif
            enddo

            MaxTotalConcentration = maxval(TotalConcentration(:))

            ! Generate the cloud mask base on integrated cloud load
            if(CloudLoad(i,j).ge.CLOUDLOAD_THRESH)then
              Mask_Cloud(i,j) = .true.
            else
              Mask_Cloud(i,j) = .false.
            endif

            ! Finally, modify the output variables based on the mask
            if (Mask_Cloud(i,j))then
              ! This is a cloud, accumulate area and error-check heights
              if(MinHeight(i,j).lt.0.0_ip.and. &
                 MaxHeight(i,j).gt.0.0_ip) &
                   MinHeight(i,j) = 0.0_ip
              ! Double-check that min doesn't exceed max
              MinHeight(i,j)=min(MaxHeight(i,j),MinHeight(i,j))

              CloudArea = CloudArea + sigma_nz_pd(i,j,1)
                ! Set the output variable for max concentration in mg/m3
              MaxConcentration(i,j) = MaxTotalConcentration * KG_2_MG/KM3_2_M3
            else
              ! Not an ash cloud column so write 0.0
              ! This will be set to FillValue in netcdf when writing out
              CloudLoad(i,j)        = 0.0_ip
              MaxConcentration(i,j) = 0.0_ip
              MaxHeight(i,j)        = 0.0_ip
              MinHeight(i,j)        = 0.0_ip
            endif

          enddo
        enddo
      endif

      Calculated_Cloud_Load = .true.

      end subroutine ConcentrationCalculator      

!******************************************************************************

      subroutine CloudAreaCalculator

      use Tephra,        only : &
         n_gs_max

      use mesh,          only : &
         sigma_nz_pd,ts1

      use solution,      only : &
         imin,imax,jmin,jmax

      implicit none

      integer :: i,j,k

      real(kind=ip) :: CellArea
 
      !calculate cloud concentration
      LoadVal(1)         = 0.24_ip
      LoadVal(2)         = 1.0_ip
      LoadVal(3)         = 2.0_ip
      LoadVal(4)         = 4.0_ip
      LoadVal(5)         = 6.0_ip
      CloudLoadArea(1:5) = 0.0_ip

      if(n_gs_max.gt.0)then
        do i=imin,imax
          do j=jmin,jmax
            CellArea = sigma_nz_pd(i,j,1)
            do k=1,5
              if (CloudLoad(i,j).gt.LoadVal(k)) &
                CloudLoadArea(k) = CloudLoadArea(k) + CellArea
            enddo
          enddo
        enddo
      endif

      end subroutine CloudAreaCalculator      

!******************************************************************************

      subroutine Gen_Output_Vars

      use mesh,          only : &
         nxmax, nymax

      use time_data,     only : &
         time

      use io_data,       only : &
         Called_Gen_Output_Vars

      implicit none

      integer :: i,j

      !real(kind=ip) :: thickness
      !real(kind=ip) :: CloudLoadHere

      Mask_Cloud   = .false.
      Mask_Deposit = .false.

      if(.not.Calculated_AshThickness)then
        ! Need to calculate deposit thickness since this is a proxy for
        ! where the deposit is located and where FillValues should used
        call AshThicknessCalculator
      endif

      if(.not.Calculated_Cloud_Load)then
        ! Need to calculate cloud load since this is a proxy for
        ! where the cloud is located and where FillValues should used
        call ConcentrationCalculator
      endif

      call CloudAreaCalculator

      ! Mark the arrival time of any new deposit
      do i=1,nxmax
        do j=1,nymax
          !thickness = DepositThickness(i,j)
          ! Note: this check below requires that the fill value be < 0
          if(Mask_Deposit(i,j).and.DepArrivalTime(i,j).lt.0.0_ip)then
            DepArrivalTime(i,j)=time
          endif
          !CloudLoadHere = CloudLoad(i,j)
          if((Mask_Cloud(i,j)).and.(CloudArrivalTime(i,j).lt.0.0_ip))then
            CloudArrivalTime(i,j)=time
          endif
        enddo
      enddo

      Called_Gen_Output_Vars = .true.

      end subroutine Gen_Output_Vars

!******************************************************************************

      subroutine Calc_AshVol_Aloft(vol)

      use mesh,          only : &
         nxmax,nymax,nzmax,nsmax,kappa_pd,ts0

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
        mass_aloft(isize) = sum(concen_pd(1:nxmax,1:nymax,1:nzmax,isize,ts0) *   & ! in kg/km^3
                                 kappa_pd(1:nxmax,1:nymax,1:nzmax))                ! convert to kg
      enddo

      ! Now loop over just the tephra bins (first n_gs_max bins) and
      ! calculate volume
      if(n_gs_max.gt.0)then
        do isize=1,n_gs_max
            ! Increment total ash in air
          vol = vol + mass_aloft(isize)               /   & ! in kg
                      MagmaDensity                    /   & ! convert to m3
                      KM3_2_M3                              ! convert to km3

        enddo
      endif

      end subroutine Calc_AshVol_Aloft

!******************************************************************************

      subroutine Calc_vprofile(itime)

      use global_param,  only : &
         EPS_THRESH

      use io_data,       only : &
         nvprofiles,i_vprofile,j_vprofile

      use mesh,          only : &
         nzmax,ts1

      use solution,      only : &
         concen_pd

      use time_data,     only : &
         time

      use Tephra,         only : &
         n_gs_max

      use time_data,      only : &
         ntmax

      implicit none

      integer, intent(in) :: itime

      integer :: i,k
      real(kind=ip) :: totalash

      if(itime.gt.ntmax)then
        write(global_error)"ERROR: itime is greater than ntmax"
        write(global_error)"       cannot write to profile"
        return
      else
        time_native(itime) = time
      endif


      do i=1,nvprofiles
        ! Get the total ash aloft in the coloumn at this point in kg/km3
        totalash = sum(concen_pd(i_vprofile(i),j_vprofile(i),1:nzmax,1:n_gs_max,ts1))
        ! don't write if there's no ash
        !if(totalash.lt.EPS_THRESH) cycle
        if(totalash.lt.CLOUDCON_THRESH) cycle
        do k=1,nzmax
          pr_ash(k,itime,i) = sum(concen_pd(i_vprofile(i),j_vprofile(i),k,1:n_gs_max,ts1))&
                              /1000.0_ip     !convert from kg/km3 to mg/m3
        enddo
      enddo

      return

      end subroutine Calc_vprofile

!******************************************************************************

      subroutine Calc_AshVol_Deposit(vol)

      use mesh,          only : &
         nxmax,nymax,kappa_pd

      use solution,      only : &
         DepositGranularity

      use Tephra,        only : &
         MagmaDensity,n_gs_max

      implicit none

      real(kind=ip),intent(out) :: vol ! Total volume of ash in deposit
      integer :: isize

      vol = 0.0_ip

      if(n_gs_max.gt.0)then
        do isize=1,n_gs_max
          vol = vol + sum(DepositGranularity(1:nxmax,1:nymax,isize) *   & ! in kg/km^3
                           kappa_pd(1:nxmax,1:nymax,0))          /   & ! convert to kg
                      MagmaDensity                      /   & ! convert to m3
                      KM3_2_M3                                ! convert to km3
        enddo
      endif

      end subroutine Calc_AshVol_Deposit

!******************************************************************************

      subroutine Calc_AshVol_Outflow(vol)

      use mesh,          only : &
         nxmax,nymax,nzmax,kappa_pd

      use solution,      only : &
         outflow_xz1_pd,outflow_xz2_pd,outflow_yz1_pd,outflow_yz2_pd,outflow_xy2_pd

      use Tephra,        only : &
         MagmaDensity,n_gs_max

      implicit none

      real(kind=ip),intent(out) :: vol ! Total volume of ash flowing out of gird

      integer :: isize

      vol = 0.0_ip

      if(n_gs_max.gt.0)then
        do isize=1,n_gs_max
          vol = vol + (                                &
                        sum(outflow_yz1_pd(        1:nymax,1:nzmax,isize)*   &
                                  kappa_pd(      0,1:nymax,1:nzmax)) +       &
                        sum(outflow_yz2_pd(        1:nymax,1:nzmax,isize)*   &
                                  kappa_pd(nxmax+1,1:nymax,1:nzmax)) +       &
                        sum(outflow_xz1_pd(1:nxmax,        1:nzmax,isize)   *   &
                                  kappa_pd(1:nxmax,        0,1:nzmax)) +          &
                        sum(outflow_xz2_pd(1:nxmax,        1:nzmax,isize)  *   &
                                  kappa_pd(1:nxmax,nymax+1,1:nzmax)) +         &
                        sum(outflow_xy2_pd(1:nxmax,1:nymax          ,isize)   *   &
                                  kappa_pd(1:nxmax,1:nymax,  nzmax+1)) )     /   & ! convert to kg
                        MagmaDensity                            /   & ! convert to m3
                        KM3_2_M3                                      ! convert to km3

        enddo
      endif

      end subroutine Calc_AshVol_Outflow


!******************************************************************************

      subroutine FirstAsh

!     Subroutine that determines whether the ash has yet hit any airports

      use Airports,      only : &
         Airport_CloudHere,Airport_thickness,Airport_depRate,Airport_AshDuration,&
         Airport_thicknessLast,Airport_depRateLast,Airport_CloudHereLast,&
         Airport_CloudDuration,Airport_CloudArrived,Airport_CloudArrivalTime,&
         Airport_AshArrived,Airport_AshArrivalTime,nairports,Airport_i,Airport_j,&
           bilinear_thickness

      use time_data,     only : &
         time, dt

      integer :: i

      ! Requires that AshThicknessCalculator and ConcentrationCalculator were called
      if(.not.Calculated_Cloud_Load)then
        ! Need to calculate cloud load since this is a proxy for
        ! where the cloud is located and where FillValues should used
        call ConcentrationCalculator
      endif
      if(.not.Calculated_AshThickness)then
        ! Need to calculate deposit thickness since this is a proxy for
        ! where the deposit is located and where FillValues should used
        call AshThicknessCalculator
      endif

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
          if ((Airport_depRate(i).lt.DEPRATE_THRESH).and.&
              (Airport_depRateLast(i).ge.DEPRATE_THRESH)) then
            Airport_AshDuration(i) = time-Airport_AshArrivalTime(i)
          endif
        endif

        !mark cloud duration if cloud has passed
        if (Airport_CloudArrived(i).eqv..true.) then
          if ((Airport_CloudHere(i).le.CLOUDLOAD_THRESH).and.&
              (Airport_CloudHereLast(i).gt.CLOUDLOAD_THRESH)) then
            Airport_CloudDuration(i) = time-Airport_CloudArrivalTime(i)
          endif
        endif

        !For airports where ash has not yet arrived . . .
        !if ash load>CLOUDLOAD_THRESH T/km2 , call it "arrived"
        if ((Airport_CloudArrived(i).eqv..false.).and.&
            (Airport_CloudHere(i).gt.CLOUDLOAD_THRESH)) then
          Airport_CloudArrived(i) = .true.
          Airport_CloudArrivalTime(i) = time
        endif
        !if ash thickness>THICKNESS_THRESH mm, call it "arrived"
        if ((Airport_AshArrived(i).eqv..false.).and.&
            (Airport_thickness(i).gt.THICKNESS_THRESH)) then
          Airport_AshArrived(i) = .true.
          Airport_AshArrivalTime(i) = time
          ! Some cases with high eruptive volume might have a deposit that
          ! arrived before the 'cloud' is triggered.  Force the cloud to be
          ! marked
          if (Airport_CloudArrived(i).eqv..false.) then
            Airport_CloudArrived(i) = .true.
            Airport_CloudArrivalTime(i) = time
          endif
        endif
      enddo

      return

      end subroutine FirstAsh

      end module Output_Vars
