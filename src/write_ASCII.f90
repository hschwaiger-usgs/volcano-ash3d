!      module Ash3d_ASCII_IO
!
!      use precis_param
!
!      use io_units
!
!      implicit none
!
!        ! These arrays are only used when reading an output file of unknown size
!      real(kind=ip), dimension(:,:),allocatable :: R_XY
!      integer       :: R_nx,R_ny
!      real(kind=ip) :: R_xll,R_yll
!      real(kind=ip) :: R_dx,R_dy
!      real(kind=ip) :: R_Fill
!
!      contains

!******************************************************************************

      subroutine vprofileopener

!     subroutine that opens files for vertical profiles

      use io_units

      use io_data,       only : &
         nvprofiles,x_vprofile, y_vprofile

      use mesh,          only : &
         nzmax,z_cc_pd

      implicit none

      integer :: i,j

      integer  ::  ionumber            !number of output file
      character(len=14)  :: cio

      do i=1,nvprofiles
        ionumber = 200+i
        if (i.lt.10) then
          write(cio,1) i
1         format('vprofile0',i1,'.txt')
        elseif (i.lt.100) then
          write(cio,2) i
2         format('vprofile',i2,'.txt')
        else
          write(global_info,*)"Too many vertical profiles."
          write(global_info,*)"nvprofiles must be < 100."
          stop 1
        endif
        open(unit=ionumber,file=cio)
        write(ionumber,3) x_vprofile(i), y_vprofile(i)
3       format('Vertical profile data for location',/, &
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

!******************************************************************************

      subroutine vprofilewriter(itime)

!     subroutine that writes data on vertical profiles

      use precis_param

      use global_param,  only : &
         KG_2_MG,KM3_2_M3

      use io_data,       only : &
         nvprofiles

      use mesh,          only : &
         nzmax,ts1

      use time_data,     only : &
         SimStartHour,time,BaseYear,useLeap,OutputOffset

      use Output_Vars,    only : &
         CLOUDCON_THRESH,pr_ash

      implicit none

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
        ionumber = 200+i
        cio = HS_yyyymmddhh_since(SimStartHour+time+OutputOffset,&
                                  BaseYear,useLeap)
        write(ionumber,1) cio, time, (pr_ash(k,itime,i), k=1,nzmax)

1       format(a13,',',f10.3,',',50(e15.3,','))
      enddo

      return

      end subroutine vprofilewriter

!******************************************************************************

      subroutine vprofilecloser

!     subroutine that closes vertical profile files

      use io_data,       only : &
         nvprofiles

      implicit none
      integer  ::  i, ionumber

      do i=1,nvprofiles
        ionumber = 200+i
        close(ionumber)
      enddo

      return

      end subroutine vprofilecloser

!******************************************************************************

      subroutine write_2D_ASCII(nx,ny,OutVar,VarMask,Fill_Value,filename_root)

!     Subroutine that writes out 2-D arrays in ESRI ASCII raster format

      use precis_param

      use io_units

      use global_param,  only  : &
         KM_2_M

      use mesh,          only : &
         dx,dy,de,dn,IsLatLon,latLL,lonLL,xLL,yLL

      use io_data,       only : &
         isFinal_TS,iout3d,WriteTimes

      implicit none

      integer          ,intent(in) :: nx,ny
      real(kind=ip)    ,intent(in) :: OutVar(nx,ny)
      logical          ,intent(in) :: VarMask(nx,ny)
      character(len=6) ,intent(in) :: Fill_Value
      character(len=20),intent(in) :: filename_root

      real(kind=op)  :: OVar(nx,ny)
      real(kind=op)  :: FValue
      integer :: fid
      integer :: i,j
      character (len=9)  :: cio
      character(len=50)  :: filename_out

      read(Fill_Value,*)FValue

      fid = 30

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
      if(INDEX(filename_root,'ArrivalTime').gt.0)then
          ! For the special cases of DepositArrivalTime.dat and
          ! CloudArrivalTime.dat
        write(filename_out,*)trim(adjustl(filename_root)),'.dat'
      else
        write(filename_out,*)trim(adjustl(filename_root)),cio,'.dat'
      endif

      open(unit=fid,file=trim(adjustl(filename_out)), status='unknown',err=2500)

      write(fid,3000) nx        ! write header values
      write(fid,3001) ny
      if (IsLatLon) then
        write(fid,3002) lonLL    
        write(fid,3003) latLL    
        write(fid,3004) de,dn
      else
        write(fid,3002) xLL*KM_2_M    ! convert xLL from km to meters so ArcMap can read it
        write(fid,3003) yLL*KM_2_M    ! same with yLL
        write(fid,3004) dx*KM_2_M,dy*KM_2_M    ! and with dx and dy
      endif
      write(fid,3005)Fill_Value

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
        write(fid,3006) (OVar(i,j), i=1,nx)
        write(fid,*)                                         !make a blank line between rows
      enddo
      
      close(fid)

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
2500  write(global_info,*) 'Error opening output file ASCII_output_file.txt.  Program stopped'            
      write(global_log ,*) 'Error opening output file ASCII_output_file.txt.  Program stopped'
      stop 1
      
      end subroutine write_2D_ASCII

!******************************************************************************

      subroutine read_2D_ASCII(filename)

!     Subroutine that reads in 2-D arrays in ESRI ASCII raster format

      use precis_param

      use io_units

      use Output_Vars

      implicit none

      character(len=50),intent(in) :: filename

      integer :: fid
      integer :: i,j

      fid = 40

      open(unit=fid,file=trim(adjustl(filename)), status='old',err=2500)

      read(fid,3000) R_nx        ! read header values
      read(fid,3001) R_ny
      allocate(R_XY(R_nx,R_ny))
      read(fid,3002) R_xll
      read(fid,3003) R_yll
      read(fid,3004) R_dx,R_dy
      read(fid,3005) R_Fill

      do j=R_ny,1,-1
        read(fid,3006) (R_XY(i,j), i=1,R_nx)
        read(fid,*)
      enddo

      close(fid)

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
2500  write(global_info,*) 'Error opening ASCII file. Program stopped'
      write(global_log ,*) 'Error opening ASCII file. Program stopped'
      stop 1

      end subroutine read_2D_ASCII

!******************************************************************************
 
      subroutine write_3D_ASCII(cio)

!     Subroutine that writes out 3-D arrays in ESRI ASCII raster format

      use precis_param

      use mesh,          only : &
         nxmax,nymax,nzmax,lon_cc_pd,lat_cc_pd,IsLatLon,&
         x_cc_pd,y_cc_pd,z_cc_pd,ts1

      use solution,      only : &
         concen_pd

      implicit none

      character(len=13) ,intent(in) :: cio

      integer :: i,j,k
      real(kind=ip)     :: rhom
      character(len=32) :: DepOutfileName
      logical,save :: first_time = .true.

      ! Output data in ASCII format

      DepOutfileName='3d_tephra_fall_'//cio//'.dat'
      open(unit=100,file=DepOutfileName,status='replace')
      write(100,*)&
      'VARIABLES = "X","Y","Z","AshConc"'
      if(first_time)then
        write(100,*) 'ZONE I = ',nxmax,' J = ',nymax,' K = ',nzmax
        first_time = .false.
      else
        write(100,*) 'ZONE '
      endif

      do k=1,nzmax
        do j=1,nymax
          do i=1,nxmax
            rhom = sum(concen_pd(i,j,k,:,ts1)) !kg/km3
            if (IsLatLon) then
              write(100,'(3(4x,f20.3),g20.5)') &
                lon_cc_pd(i), lat_cc_pd(j), z_cc_pd(k), rhom
            else
              write(100,'(3(4x,f20.3),g20.5)') &
                x_cc_pd(i), y_cc_pd(j), z_cc_pd(k), rhom
            endif
          enddo
        enddo
      enddo
      close(100)

      end subroutine write_3D_ASCII

!******************************************************************************

      subroutine Write_PointData_Airports_ASCII

!     Subroutine that writes the arrival time at airports to an ASCII file

      use precis_param

      use io_units

      use global_param,  only : &
         M_2_MM,UseCalcFallVel

      use Airports,      only : &
         Airport_AshArrivalTime,Airport_CloudArrivalTime, &
         Airport_thickness,Airport_AshDuration, &
         Airport_CloudDuration,nairports,Airport_CloudArrived,&
         Airport_AshArrived,Airport_depRate,Airport_i,Airport_j,&
         Airport_Longitude,Airport_Latitude,Airport_Name,&
         bilinear_thickness

      use solution,      only : &
         DepositGranularity

      use Output_Vars,   only : &
         DepositThickness,DEPRATE_THRESH, &
         CloudLoad,CLOUDLOAD_THRESH

      use io_data,       only : &
         infile,WriteGSD,WriteAirportFile_ASCII,VolcanoName

      use time_data,     only : &
         time,SimStartHour,OutputOffset,BaseYear,useLeap,RunStartMinute,&
         RunStartYear,RunStartMonth,RunStartDay,RunStartHr

      use Tephra,        only : &
         n_gs_max,Tephra_gsdiam,Tephra_rho_m,Tephra_v_s

      use mesh,          only : &
         nsmax

      use Source,        only : &
         neruptions,e_Duration,e_Volume,e_PlumeHeight

      implicit none

      integer             :: i
      integer             :: nWrittenOut
      character (len=20)  :: yyyymmddhh_ash, yyyymmddhh_cloud
      character (len=1)   :: cloud_morethan, deposit_morethan      !equals ">" if cloud is still overhead or ash is still falling
      character (len=13)  :: nwsthickness                !nwp ashfall  terms (trace, minor, substantial, heavy, severe)
      integer             :: isize
      real(kind=ip)       :: longitude_now

      integer             :: out_unit

      INTERFACE
        character (len=20) function HS_xmltime(HoursSince,byear,useLeaps)
          real(kind=8)              :: HoursSince
          integer                   :: byear
          logical                   :: useLeaps
        end function HS_xmltime
      END INTERFACE

      out_unit = 18

      !WRITE VALUES OUT TO ASCII FILE
      if (WriteAirportFile_ASCII) then
!        open(unit=28,file='dep_thick.txt',err=2000)
!
!        !WRITE OUT SOURCE PARAMETERS FOR AIRPORTS FILE
        open(unit=out_unit,file='ash_arrivaltimes_airports.txt',err=2000)
        write(out_unit,98)  infile, RunStartYear,RunStartMonth,RunStartDay,RunStartHr, &
                              RunStartMinute, VolcanoName  !write infile, simulation time
        do i=1,neruptions  !write source parameters
          write(out_unit,99) i, HS_xmltime(SimStartHour+OutputOffset,&
                                           BaseYear,useLeap), &
                           e_Duration(i), e_PlumeHeight(i), e_PlumeHeight(i)*3280.8_ip, e_Volume(i)
        enddo
        write(out_unit,995)
        if (WriteGSD) then                 !If we're writing out grain sizes.
          if (UseCalcFallVel) then              !if fall velocity is calculated
            write(out_unit,100) (Tephra_gsdiam(isize)*M_2_MM, isize=1,n_gs_max)
            write(out_unit,101)
            write(out_unit,102) (Tephra_rho_m(isize), isize=1,n_gs_max)
            write(out_unit,103)
          else                                !if fall velocity is specified
            write(out_unit,110) (Tephra_v_s(isize), isize=1,n_gs_max)
          endif
        else                             !If not
          write(out_unit,1)
        endif
        nWrittenOut = 0
        do i=1,nairports                      !write out the airports that are hit.
          if (Airport_AshArrived(i).or.Airport_CloudArrived(i)) then

            !rank ash thickness in NWS rank
            Airport_thickness(i) = bilinear_thickness(i,DepositThickness) !interpolate to find thickness
            if (Airport_thickness(i).le.0.7935_ip) then         !<1/32" thickness
              nwsthickness="trace or less"
            elseif (Airport_thickness(i).le.6.35_ip) then      !<=1/4"
              nwsthickness="minor"
            elseif (Airport_thickness(i).le.25.4_ip) then      !<=1"
              nwsthickness="substantial"
            elseif (Airport_thickness(i).le.101.6_ip) then     !<=4"
              nwsthickness="heavy"
            else
              nwsthickness="severe"
            endif

            !get yyyymmddhh of arrival
            yyyymmddhh_cloud = HS_xmltime(Airport_CloudArrivalTime(i)+SimStartHour+OutputOffset,&
                                          BaseYear,useLeap)
            if (Airport_AshArrived(i)) then
              yyyymmddhh_ash = HS_xmltime(Airport_AshArrivaltime(i)+SimStartHour+OutputOffset,&
                                          BaseYear,useLeap)
            else
              yyyymmddhh_ash = '0000-00-00T00:00:00Z'
            endif

            !See whether cloud is still overhead, or whether ash is still
            !falling
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
              write(out_unit,20) Airport_Name(i), Airport_Latitude(i), longitude_now, &
                       yyyymmddhh_cloud, Airport_CloudArrivalTime(i), cloud_morethan, Airport_CloudDuration(i), &
                       yyyymmddhh_ash, Airport_AshArrivalTime(i), deposit_morethan, Airport_AshDuration(i), &
                       Airport_thickness(i), nwsthickness, &
                       ((DepositGranularity(Airport_i(i),Airport_j(i),isize)/ &
                         sum(DepositGranularity(Airport_i(i),Airport_j(i),:))),isize=1,nsmax)  !mass fraction of size i
                      !(DepositGranularity(Airport_i(i),Airport_j(i),isize)*dz_vec_pd(0)/1.0e6_ip,isize=1,nsmax) !mpua, kg/m2 of size i
            else
              write(out_unit,2) Airport_Name(i), Airport_Latitude(i), longitude_now, &
                       yyyymmddhh_cloud, Airport_CloudArrivalTime(i), cloud_morethan, Airport_CloudDuration(i), &
                       yyyymmddhh_ash, Airport_AshArrivalTime(i), deposit_morethan, Airport_AshDuration(i), &
                       Airport_thickness(i), nwsthickness
            endif
            nWrittenOut = nWrittenOut + 1
          endif
        enddo
        if (nWrittenOut.eq.0) write(out_unit,3)    !If no airports are hit, say it in the file.
        write(out_unit,120)                        !write footnotes & caveats
        close(out_unit)
      endif

      return

!     Error traps
2000   write(global_info,*)  'Error opening ash_arrivaltimes_airports.txt.  Program stopped.'
       write(global_info,*)  'Error opening ash_arrivaltimes_airports.txt.  Program stopped.'
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

!100   format('---------------------------------------------------------------------------------------------------', &
!             '---------------------------------------------------------------',/, &
!             '                      LOCATION                          |                  CLOUD                  |', &
!             '                               DEPOSIT                         |', &
!             '    MASS PER UNIT AREA AT THE GIVEN GRAIN SIZE (MM)',/, &
!             '                                                        |                                         |', &
!             '                                                               |',/, &
!             '                                                        |                                         |', &
!             '                                                               |',/, &
!             '(Airport code & ) Place name         Latitude Longitude |        Cloud Arrival Time      Duration |', &
!             '      Deposit Arrival Time        Duration      Thickness      |',/, &
!             '                                                        |                      hrs after          |', &
!             '                                                               |',/, &
!             '                                                        |      date/time UTC    start      hrs    |', &
!             '       date/time UTC       hrs      hrs      mm   NWS rank     |',30f12.4)
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

!      end module Ash3d_ASCII_IO

