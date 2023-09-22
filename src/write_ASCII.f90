!##############################################################################
!
!      Ash3d_ASCII_IO module
!
!  This module manages all output to ASCII files
!
!      subroutine vprofileopener
!      subroutine vprofilewriter(itime)
!      subroutine vprofilecloser
!      subroutine write_2D_ASCII(nx,ny,OutVar,VarMask,Fill_Value,filename_root)
!      subroutine read_2D_ASCII(filename)
!      subroutine write_3D_ASCII(cio)
!      subroutine Write_PointData_Airports_ASCII
!
!##############################################################################

      module Ash3d_ASCII_IO

      use precis_param

      use io_units

      implicit none

        ! Set everything to private by default
      private

        ! Publicly available subroutines/functions
      public vprofileopener,     &
             vprofilewriter,     &
             vprofilecloser,     &
             write_2D_ASCII,     &
             read_2D_ASCII,      &
             write_3D_ASCII,     &
             Write_PointData_Airports_ASCII

        ! Publicly available variables

        ! These arrays are only used when reading an output file of unknown size
      real(kind=ip), dimension(:,:),allocatable,public :: R_XY
      integer      ,public :: R_nx
      integer      ,public :: R_ny
      real(kind=ip),public :: R_xll
      real(kind=ip),public :: R_yll
      real(kind=ip),public :: R_dx
      real(kind=ip),public :: R_dy
      real(kind=ip),public :: R_Fill

      contains
      !------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  vprofileopener
!
!  Called from: output_results and Ash3d_PostProc.F90
!  Arguments:
!    none
!
!  This subroutine opens files for vertical profiles
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine vprofileopener

      use io_data,       only : &
         MAXPROFILES,nvprofiles,Site_vprofile,x_vprofile, y_vprofile

      use mesh,          only : &
         nzmax,z_cc_pd

      integer :: i,j

      integer  ::  ionumber            !number of output file
      character(len=14)  :: cio

      do i=1,nvprofiles
        ionumber = fid_vprofbase + i-1
        if (i.lt.10) then
          write(cio,1) i
1         format('vprofile0',i1,'.txt')
        elseif (i.lt.MAXPROFILES) then
          write(cio,2) i
2         format('vprofile',i2,'.txt')
        else
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"ERROR: Too many vertical profiles."
            write(errlog(io),*)"nvprofiles must be < ",MAXPROFILES
          endif;enddo
          stop 1
        endif
        open(unit=ionumber,file=cio,status='replace',action='write')
        write(ionumber,*)'Vertical profile data for location: ',trim(adjustl(Site_vprofile(i)))
        write(ionumber,3)x_vprofile(i), y_vprofile(i)
3       format( &
       'x:',f10.3,/, &
       'y:',f10.3,/, &
       '                          Output is ash concentration in mg/m3',/, &
       '                          elevation (km) ---->',/, &
       'date-time         hrs')
        write(ionumber,4) (z_cc_pd(j), j=1,nzmax)
4       format(50f16.3)
      enddo

      return

      end subroutine vprofileopener

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  vprofilewriter
!
!  Called from: Ash3d.F90 and Ash3d_PostProc.F90
!  Arguments:
!    itime = time step index
!
!  This subroutine writes data on vertical profiles
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine vprofilewriter(itime)

      use global_param,  only : &
         KG_2_MG,KM3_2_M3

      use io_data,       only : &
         nvprofiles

      use mesh,          only : &
         nzmax

      use time_data,     only : &
         SimStartHour,time,BaseYear,useLeap,OutputOffset

      use Output_Vars,    only : &
         CLOUDCON_THRESH,pr_ash

      integer, intent(in) :: itime

      integer :: i,k
      integer   :: ionumber
      character(len=13)  :: cio

      INTERFACE
        character (len=13) function HS_yyyymmddhh_since(HoursSince,byear,useLeaps)
          real(kind=8)               ::  HoursSince
          integer                    ::  byear
          logical                    ::  useLeaps
        end function HS_yyyymmddhh_since
      END INTERFACE

      do i=1,nvprofiles
        ! don't write if there's no ash
        if(maxval(pr_ash(1:nzmax,itime,i)).lt.CLOUDCON_THRESH*KG_2_MG/KM3_2_M3) cycle
        ionumber = fid_vprofbase + i-1
        cio = HS_yyyymmddhh_since(SimStartHour+time+OutputOffset,&
                                  BaseYear,useLeap)
        write(ionumber,1) cio, time, (pr_ash(k,itime,i), k=1,nzmax) ! write tot. concen in mg/m3

1       format(a13,',',f10.3,',',50(e15.3,','))
      enddo

      return

      end subroutine vprofilewriter


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  vprofilecloser
!
!  Called from: output_results and Ash3d_PostProc.F90
!  Arguments:
!    none
!
!  This subroutine closes vertical profile files
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine vprofilecloser

      use io_data,       only : &
         nvprofiles

      integer  ::  i, ionumber

      do i=1,nvprofiles
        ionumber = fid_vprofbase + i-1
        close(ionumber)
      enddo

      return

      end subroutine vprofilecloser

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  write_2D_ASCII
!
!  Called from: output_results and Ash3d_PostProc.F90
!  Arguments:
!    nx            = x length of output array OutVar
!    ny            = y length of output array OutVar
!    OutVar        = 2-d array to be written to ASCII file
!    VarMask       = 2-d logical that toggles data v.s. Fill_Value
!    Fill_Value    = number used for No-data
!    filename_root = root name of file (20 characters)
!
!  Subroutine that writes out 2-D arrays in ESRI ASCII raster format
!  This format can be post-processed with gmt converting to grid files with
!  gmt grdconvert out.dat=ef out.grd
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine write_2D_ASCII(nx,ny,OutVar,VarMask,Fill_Value,filename_root)

      use global_param,  only  : &
         KM_2_M

      use mesh,          only : &
         dx,dy,de,dn,IsLatLon,latLL,lonLL,xLL,yLL

      use io_data,       only : &
         isFinal_TS,iout3d,WriteTimes

      integer          ,intent(in) :: nx,ny
      real(kind=ip)    ,intent(in) :: OutVar(nx,ny)
      logical          ,intent(in) :: VarMask(nx,ny)
      character(len=6) ,intent(in) :: Fill_Value
      character(len=20),intent(in) :: filename_root

      real(kind=op)  :: OVar(nx,ny)
      real(kind=op)  :: FValue
      integer :: i,j
      character (len=9)  :: cio
      character(len=50)  :: filename_out

      read(Fill_Value,*)FValue

      !fid = 30

      if(isFinal_TS)then
        cio='____final'
      else
        if (WriteTimes(iout3d).lt.10.0_ip) then
          write(cio,1) WriteTimes(iout3d)
1         format('00',f4.2,'hrs')
        elseif (WriteTimes(iout3d).lt.100.0_ip) then
          write(cio,2) WriteTimes(iout3d)
2         format('0',f5.2,'hrs')
        else
          write(cio,3) WriteTimes(iout3d)
3         format(f6.2,'hrs')
        endif
      endif
      if(index(filename_root,'ArrivalTime').gt.0)then
          ! For the special cases of DepositArrivalTime.dat and
          ! CloudArrivalTime.dat
        write(filename_out,*)trim(adjustl(filename_root)),'.dat'
      else
        write(filename_out,*)trim(adjustl(filename_root)),cio,'.dat'
      endif

      open(unit=fid_ascii2dout,file=trim(adjustl(filename_out)), status='replace',action='write',err=2500)

      write(fid_ascii2dout,3000) nx        ! write header values
      write(fid_ascii2dout,3001) ny
      if (IsLatLon) then
        write(fid_ascii2dout,3002) lonLL    
        write(fid_ascii2dout,3003) latLL    
        write(fid_ascii2dout,3004) de,dn
      else
        write(fid_ascii2dout,3002) xLL*KM_2_M    ! convert xLL from km to meters so ArcMap can read it
        write(fid_ascii2dout,3003) yLL*KM_2_M    ! same with yLL
        write(fid_ascii2dout,3004) dx*KM_2_M,dy*KM_2_M    ! and with dx and dy
      endif
      write(fid_ascii2dout,3005)Fill_Value

      ! Apply threshold mask to output variable
      do i=1,nx
        do j=1,ny
          if(VarMask(i,j))then
            OVar(i,j) = real(OutVar(i,j),kind=op)
          else
            OVar(i,j) = FValue
          endif
        enddo
      enddo

      !Write out arrays of maximum concentration and maximum height
      do j=ny,1,-1
        write(fid_ascii2dout,3006) (OVar(i,j), i=1,nx)
        write(fid_ascii2dout,*)                                         !make a blank line between rows
      enddo
      
      close(fid_ascii2dout)

!     format statements
3000  format('NCOLS ',i5)      
3001  format('NROWS ',i5)
3002  format('XLLCORNER ',f15.3)
3003  format('YLLCORNER ',f15.3)
3004  format('CELLSIZE ',2f15.3)
3005  format('NODATA_VALUE ',a6)
3006  format(10f15.3)
      return

!     Error traps
2500  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*) 'Error opening output file ASCII_output_file.txt.  Program stopped'
      endif;enddo
      stop 1
      
      end subroutine write_2D_ASCII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  read_2D_ASCII
!
!  Called from: Ash3d_ASCII_check.f90
!  Arguments:
!    filename = root name of file (20 characters)
!
!  Subroutine that reads in 2-D arrays in ESRI ASCII raster format and
!  populates R_nx,R_ny,R_XY,R_xll,R_yll,R_dx,R_dy,R_Fill from Output_Vars
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine read_2D_ASCII(filename)

      character(len=80),intent(in) :: filename

      !integer :: fid
      integer :: i,j
      integer :: iost

      !fid = 40

      open(unit=fid_ascii2din,file=trim(adjustl(filename)), status='old',action='read',err=2500)

      read(fid_ascii2din,3000,iostat=iost,err=2600) R_nx        ! read header values
      if(iost.gt.0)then
        ! We might have an empty file
        ! Issue warning and return
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*) 'Error reading file ',trim(adjustl(filename))
          write(errlog(io),*) 'Check for zero-length file.'
        endif;enddo
        return
      endif
      read(fid_ascii2din,3001,err=2600) R_ny
      allocate(R_XY(R_nx,R_ny))
      read(fid_ascii2din,3002,err=2600) R_xll
      read(fid_ascii2din,3003,err=2600) R_yll
      read(fid_ascii2din,3004,err=2600) R_dx,R_dy
      read(fid_ascii2din,3005,err=2600) R_Fill

      do j=R_ny,1,-1
        read(fid_ascii2din,3006,err=2600) (R_XY(i,j), i=1,R_nx)
        read(fid_ascii2din,*)
      enddo

      close(fid_ascii2din)

!     format statements
3000  format(6x,i5)
3001  format(6x,i5)
3002  format(10x,f15.3)
3003  format(10x,f15.3)
3004  format(10x,2f15.3)
3005  format(13x,a6)
3006  format(10f15.3)
      return

!     Error traps
2500  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*) 'Error opening ASCII file. Program stopped'
      endif;enddo
      stop 1

2600  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*) 'Error reading from ASCII file.'
      endif;enddo
      stop 1

      end subroutine read_2D_ASCII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  write_3D_ASCII
!
!  Called from: output_results
!  Arguments:
!    cio = time string to be inserted into filename; either '________final'
!          or yyyymmddhh.h
!
!  Subroutine that writes out 3-D arrays in ESRI ASCII raster format
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine write_3D_ASCII(cio)

      use mesh,          only : &
         nxmax,nymax,nzmax,lon_cc_pd,lat_cc_pd,IsLatLon,&
         x_cc_pd,y_cc_pd,z_cc_pd,ts1

      use solution,      only : &
         concen_pd

      character(len=13) ,intent(in) :: cio

      integer :: i,j,k
      real(kind=ip)     :: rhom
      character(len=32) :: DepOutfileName
      logical,save :: first_time = .true.

      ! Output data in ASCII format

      DepOutfileName='3d_tephra_fall_'//cio//'.dat'
      open(unit=fid_ascii3dout,file=DepOutfileName,status='replace',action='write')
      write(fid_ascii3dout,*)&
      'VARIABLES = "X","Y","Z","AshConc"'
      if(first_time)then
        write(fid_ascii3dout,*) 'ZONE I = ',nxmax,' J = ',nymax,' K = ',nzmax
        first_time = .false.
      else
        write(fid_ascii3dout,*) 'ZONE '
      endif

      do k=1,nzmax
        do j=1,nymax
          do i=1,nxmax
            rhom = sum(concen_pd(i,j,k,:,ts1)) !kg/km3
            if (IsLatLon) then
              write(fid_ascii3dout,'(3(4x,f20.3),g20.5)') &
                lon_cc_pd(i), lat_cc_pd(j), z_cc_pd(k), rhom
            else
              write(fid_ascii3dout,'(3(4x,f20.3),g20.5)') &
                x_cc_pd(i), y_cc_pd(j), z_cc_pd(k), rhom
            endif
          enddo
        enddo
      enddo
      close(fid_ascii3dout)

      end subroutine write_3D_ASCII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Write_PointData_Airports_ASCII
!
!  Called from: output_results and Ash3d_PostProc.f90
!  Arguments:
!    none
!
!  Subroutine that writes the arrival time at airports to an ASCII file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Write_PointData_Airports_ASCII

      use global_param,  only : &
         M_2_MM,UseCalcFallVel

      use io_data,       only : &
         infile,WriteGSD,WriteAirportFile_ASCII,VolcanoName

      use mesh,          only : &
         nsmax

       use solution,      only : &
         DepositGranularity

      use time_data,     only : &
         time,SimStartHour,OutputOffset,BaseYear,useLeap,RunStartMinute,&
         RunStartYear,RunStartMonth,RunStartDay,RunStartHr

      use Airports,      only : &
         Airport_AshArrivalTime,Airport_CloudArrivalTime, &
         Airport_thickness,Airport_AshDuration, &
         Airport_CloudDuration,nairports,Airport_CloudArrived,&
         Airport_AshArrived,Airport_depRate,Airport_i,Airport_j,&
         Airport_Longitude,Airport_Latitude,Airport_Name,&
         bilinear_thickness

      use Output_Vars,   only : &
         DepositThickness,DEPRATE_THRESH, &
         CloudLoad,CLOUDLOAD_THRESH

      use Tephra,        only : &
         n_gs_max,Tephra_gsdiam,Tephra_rho_m,Tephra_v_s

      use Source,        only : &
         neruptions,e_Duration,e_Volume,e_PlumeHeight

      integer             :: i
      integer             :: nWrittenOut
      character (len=20)  :: yyyymmddhh_ash, yyyymmddhh_cloud
      character (len=1)   :: cloud_morethan, deposit_morethan      !equals ">" if cloud is still overhead or ash is still falling
      character (len=13)  :: nwsthickness                !nwp ashfall  terms (trace, minor, substantial, heavy, severe)
      integer             :: isize
      real(kind=ip)       :: longitude_now

      INTERFACE
        character (len=20) function HS_xmltime(HoursSince,byear,useLeaps)
          real(kind=8)              :: HoursSince
          integer                   :: byear
          logical                   :: useLeaps
        end function HS_xmltime
      END INTERFACE

      ! Write values out to ASCII file
      if (WriteAirportFile_ASCII) then

        ! Write out source parameters for airports file
        open(unit=fid_asharrive,file='ash_arrivaltimes_airports.txt',status='replace',action='write',err=2000)
        write(fid_asharrive,98)  infile, RunStartYear,RunStartMonth,RunStartDay,RunStartHr, &
                              RunStartMinute, VolcanoName  !write infile, simulation time
        do i=1,neruptions  ! write source parameters
          write(fid_asharrive,99) i, HS_xmltime(SimStartHour+OutputOffset,&
                                           BaseYear,useLeap), &
                           e_Duration(i), e_PlumeHeight(i), e_PlumeHeight(i)*3280.8_ip, e_Volume(i)
        enddo
        write(fid_asharrive,995)
        if (WriteGSD) then
          ! if we're writing out grain sizes.
          if (UseCalcFallVel) then
            ! if fall velocity is calculated
            write(fid_asharrive,100) (Tephra_gsdiam(isize)*M_2_MM, isize=1,n_gs_max)
            write(fid_asharrive,101)
            write(fid_asharrive,102) (Tephra_rho_m(isize), isize=1,n_gs_max)
            write(fid_asharrive,103)
          else
            ! if fall velocity is specified
            write(fid_asharrive,110) (Tephra_v_s(isize), isize=1,n_gs_max)
          endif
        else
          ! if not writing out grain sizes
          write(fid_asharrive,1)
        endif
        nWrittenOut = 0
        do i=1,nairports                      ! write out the airports that are hit.
          if (Airport_AshArrived(i).or.Airport_CloudArrived(i)) then
            ! rank ash thickness in NWS rank
            Airport_thickness(i) = bilinear_thickness(i,DepositThickness) ! interpolate to find thickness
            if (Airport_thickness(i).le.0.7935_ip) then         ! <1/32" thickness
              nwsthickness="trace or less"
            elseif (Airport_thickness(i).le.6.35_ip) then       ! <=1/4"
              nwsthickness="minor"
            elseif (Airport_thickness(i).le.25.4_ip) then       ! <=1"
              nwsthickness="substantial"
            elseif (Airport_thickness(i).le.101.6_ip) then      ! <=4"
              nwsthickness="heavy"
            else
              nwsthickness="severe"
            endif

            ! get yyyymmddhh of arrival
            yyyymmddhh_cloud = HS_xmltime(Airport_CloudArrivalTime(i)+SimStartHour+OutputOffset,&
                                          BaseYear,useLeap)
            if (Airport_AshArrived(i)) then
              yyyymmddhh_ash = HS_xmltime(Airport_AshArrivaltime(i)+SimStartHour+OutputOffset,&
                                          BaseYear,useLeap)
            else
              yyyymmddhh_ash = '0000-00-00T00:00:00Z'
            endif

            ! See whether cloud is still overhead, or whether ash is still
            ! falling
            if((Airport_AshArrived(i)).and.(Airport_depRate(i).gt.DEPRATE_THRESH)) then
              Airport_AshDuration(i) = time-Airport_AshArrivalTime(i)
              deposit_morethan = '>'
            else
              deposit_morethan = ' '
            endif
            if (CloudLoad(Airport_i(i),Airport_j(i)).gt.CLOUDLOAD_THRESH)then
              Airport_CloudDuration(i) = time-Airport_CloudArrivalTime(i)
              cloud_morethan = '>'
            else
              cloud_morethan = ' '
            endif

            if (Airport_Longitude(i).gt.180.0_ip) then
              longitude_now = Airport_Longitude(i)-360.0_ip
            else
              longitude_now = Airport_Longitude(i)
            endif
            if (WriteGSD) then
              write(fid_asharrive,20) Airport_Name(i), Airport_Latitude(i), longitude_now, &
                       yyyymmddhh_cloud, Airport_CloudArrivalTime(i), cloud_morethan, Airport_CloudDuration(i), &
                       yyyymmddhh_ash, Airport_AshArrivalTime(i), deposit_morethan, Airport_AshDuration(i), &
                       Airport_thickness(i), nwsthickness, &
                       ((DepositGranularity(Airport_i(i),Airport_j(i),isize)/ &
                         sum(DepositGranularity(Airport_i(i),Airport_j(i),:))),isize=1,nsmax)  ! mass fraction of size i
            else
              write(fid_asharrive,2) Airport_Name(i), Airport_Latitude(i), longitude_now, &
                       yyyymmddhh_cloud, Airport_CloudArrivalTime(i), cloud_morethan, Airport_CloudDuration(i), &
                       yyyymmddhh_ash, Airport_AshArrivalTime(i), deposit_morethan, Airport_AshDuration(i), &
                       Airport_thickness(i), nwsthickness
            endif
            nWrittenOut = nWrittenOut + 1
          endif
        enddo
        if (nWrittenOut.eq.0) write(fid_asharrive,3)    ! if no airports are hit, say it in the file.
        write(fid_asharrive,120)                        ! write footnotes & caveats
        close(fid_asharrive)
      endif

      return

!     Error traps
2000  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'Error opening ash_arrivaltimes_airports.txt.  Program stopped.'
      endif;enddo
      stop 1

!     Format statements
1     format('---------------------------------------------------------------------------------------------------', &
             '----------------------------------------------------------------',/, &
             '                      LOCATION                          |                  CLOUD                  |', &
             '                               DEPOSIT                         |',/, &
             '                                                        |                                         |', &
             '                                                               |',/, &
             '                                                        |                                         |', &
             '                                                               |',/, &
             '(Airport code & ) Place name         Latitude Longitude |        Cloud Arrival Time      Duration |', &
             '      Deposit Arrival Time        Duration      Thickness      |',/, &
             '                                                        |                        hrs after        |', &
             '                          hrs after                            |',/, &
             '                                                        |      date/time UTC    start      hrs    |', &
             '       date/time UTC      start     hrs      mm   NWS rank     |')
98    format('ARRIVAL TIME OF ASH IN AREA MODELED BY ASH3D',/, &  !if WriteGSD=.true. and CalcFallVel=.true.
             'simulation using input file: ',A130,/,                        &
             'Model run date: ',i4,'.',i2.2,'.',i2,', time UTC: ',i4,':',i2.2,//, &
             '---------------------------------------------------------------------------------',/, &
             '             SOURCE PARAMETERS FOR SIMULATION OF: ',a30,'|',/, &
             'Pulse       start time                Duration       Plume height      volume   |',/, &
             '               UTC                      hrs          km     feet         km3    |')
99    format(i3,4x,a20,f16.2,f13.1,f9.0,f13.5,'  |')
995   format('---------------------------------------------------------------------------------',/)
100   format('---------------------------------------------------------------------------------------------------', &
             '---------------------------------------------------------------',/, &
             '                      LOCATION                          |                  CLOUD                  |', &
             '                               DEPOSIT                         |  GRAIN SIZE, MM',/, &
             '                                                        |                                         |', &
             '                                                               |',30f12.4)                                    !grain size, mm
101   format('                                                        |                                         |', &
             '                                                               |  DENSITY, KG/M3')
102   format('(Airport code & ) Place name         Latitude Longitude |        Cloud Arrival Time      Duration |', &
             '      Deposit Arrival Time        Duration      Thickness      |',30f12.1)                                    !density (kg/m3)
103   format('                                                        |                      hrs after          |', &
             '                                                        |',/, &
             '                                                        |      date/time UTC    start      hrs    |', &
             '       date/time UTC       hrs      hrs      mm   NWS rank     |  MASS FRACTION AT THE GIVEN GRAIN SIZE')

110   format('---------------------------------------------------------------------------------------------------', &
             '----------------------------------------------------------------',/, &
             '                      LOCATION                          |                  CLOUD                  |', &
             '                               DEPOSIT                         |', &
             '    MASS PER UNIT AREA AT THE GIVEN FALL SPEED (M/S)',/, &
             '                                                        |                                         |', &
             '                                                               |',/, &
             '                                                        |                                         |', &
             '                                                               |',/, &
             '(Airport code & ) Place name         Latitude Longitude |        Cloud Arrival Time      Duration |', &
             '      Deposit Arrival Time        Duration      Thickness      |',/, &
             '                                                        |                                         |', &
             '                                                               |',/, &
             '                                                        |      date/time UTC     hrs       hrs    |', &
             '       date/time UTC       hrs      hrs      mm   NWS rank     |',30f12.4)
120   format(//, &
             '---------------------------------------------------------------------------------------------------', &
             '----------------------------------------------------------------',//, &
             'NOTES ON ITEMS IN THIS TABLE:',/, &
             'LOCATION: If the location is an airport, the first three letters are the ICAO airport code',/, &
             'CLOUD DATA: Cloud arrival time is given in hours after the eruption start and the date and time', &
                          ' in UTC,in the format yyyy-mm-ddThh:mm:ssZ).  Duration is the',/, &
             '  number of hours during which the cloud is overhead (or, if there is a break in the cloud, the ', &
                          'time from the clouds first arrival until it last passes).',/, &
             '  A character ">" in front of the duration indicates that the cloud was still overhead at the end ', &
                          'of the simulation. The vertically integrated cloud mass must',/, &
             '  exceed a threshold value of 0.2 tonnes per square kilometer to be considered in these calculations.',/, &
             'DEPOSIT DATA:  The deposit arrival time is given in hours since eruption start or in the date and time ', &
                          'UTC as formatted in the cloud output (The cloud and ',/, &
             '   deposit may arrive at different times). "Deposit arrival time" is the time of arrival of the deposit ', &
                          'at a thickness exceeding 0.01 mm (0.0004 inches).',/, &
             '   An arrival time of -9999.00 hrs indicates the deposit did not arrive at this location.  Deposit ', &
                          'duration is the time period (hrs) over which the deposit',/, &
             '   was falling at a rate exceeding 0.01 mm/hr.  A ">" character before this number indicates that the ', &
                          'deposit was still falling at the end of the simulation.',/, &
             '   The thickness of the deposit is given in millimeters (left column) and as ranked according to the ', &
                          'following system devised by the U.S. National Weather Service:',/, &
             '     NWS Rank         Thickness',/, &
             '                        up to',/, &
             '                     mm      in.',/, &
             '     trace           0.8     1/32"',/, &
             '     minor           6.3     1/4"',/, &
             '     substantial    25.4     1"',/, &
             '     heavy          100      4"',/, &
             '     severe        >100     >4"',//, &
             'NOTE: This table is the estimate at time of issuance: changing conditions at the volcano may require ', &
                         'updating the forecast.')
2     format(a35,2f10.4,' |',2x,a20,f7.2,3x,a1,f5.2,'   |',2x,a20,f9.2,3x,a1,f5.2,f8.2,2x,a13,'|')
20    format(a35,2f10.4,' |',2x,a20,f7.2,3x,a1,f5.2,'   |',2x,a20,f9.2,3x,a1,f5.2,f8.2,2x,a13,'|',30e12.4)
3     format(/,'No airports affected by ash')

      end subroutine Write_PointData_Airports_ASCII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module Ash3d_ASCII_IO

!##############################################################################
