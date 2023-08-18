!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  output_results
!
!  Called from: Ash3d.F90 at three points; first before the time loop, once within
!               it, and once after.  This is also called once from the post-
!               procession tool
!  Arguments:
!    none
!
!  This subroutine is primarily used to write out the requested output data on
!  the output timesteps specified in the control file.  Profile data, or other
!  data on a different output interval are exported elsewhere. There are three
!  blocks in this subroutine: first a block that prepares the output files which
!  is only executed the first time this subroutine is called (prior to the time
!  loop in Ash3d.F90), next the block that writes out the data for the output
!  step, and finally a block that closes the files at the end of the simulation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine output_results

      use precis_param

      use io_units

      use global_param,  only : &
         EPS_SMALL,MM_2_IN

      use io_data,       only : &
         iout3d,ioutputFormat,WriteTimes,nWriteTimes,isFinal_TS,&
         NextWriteTime,iTimeNext,OutputStep_Marker,&
         Write3dFiles,WriteAirportFile_ASCII,&
         WriteCloudConcentration_ASCII,WriteAirportFile_KML,&
         WriteCloudConcentration_KML,WriteCloudHeight_ASCII,&
         WriteCloudHeight_KML,WriteCloudLoad_ASCII,WriteCloudLoad_KML,&
         WriteCloudTime_ASCII,WriteCloudTime_KML,&
         WriteDepositFinal_ASCII,WriteDepositTS_KML,WriteDepositFinal_KML,&
         WriteDepositTime_ASCII,WriteDepositTime_KML,WriteDepositTS_ASCII,&
         WriteReflectivity_KML,Write_PR_Data

      use mesh,          only : &
         nxmax,nymax,nzmax

      use Output_Vars,   only : &
         DepositThickness,MaxConcentration,MaxHeight,&
         Mask_Cloud,Mask_Deposit,ashcon_tot, &
         CloudLoad,DepArrivalTime,CloudArrivalTime,dbZCol,&
         AshTotalCalculator

      use time_data,     only : &
         BaseYear,useLeap,dt,time,SimStartHour,Simtime_in_hours,&
         xmlTimeSpanStart,xmlTimeSpanEnd

      use Output_KML,    only : &
           OpenFile_KML,&
           Close_KML,&
           Write_2D_KML,&
           Write_PointData_Airports_KML,&
           Set_OutVar_Specs

#ifdef USENETCDF
      use Ash3d_Netcdf
#endif

      implicit none

      character(len=13) :: cio
      real(kind=ip)     :: timestart
      real(kind=ip)     :: timeend
      logical,save      :: first_time = .true.

      INTERFACE
        subroutine vprofileopener
        end subroutine vprofileopener
        subroutine vprofilecloser
        end subroutine vprofilecloser
        subroutine Write_PointData_Airports_ASCII
        end subroutine Write_PointData_Airports_ASCII
        subroutine write_2D_ASCII(nx,ny,OutVar,VarMask,Fill_Value,filename_root)
          integer,parameter  :: ip         = 8 ! Internal precision
          integer          ,intent(in) :: nx,ny
          real(kind=ip)    ,intent(in) :: OutVar(nx,ny)
          logical          ,intent(in) :: VarMask(nx,ny)
          character(len=6) ,intent(in) :: Fill_Value
          character(len=20),intent(in) :: filename_root
        end subroutine write_2D_ASCII
        subroutine write_3D_ASCII(cio)
          character(len=13) ,intent(in) :: cio
        end subroutine write_3D_ASCII
        subroutine write_3D_Binary(cio,nx,ny,nz,ashcon_tot)
          integer,parameter  :: op         = 4 ! Output precision
          character(len=13) ,intent(in) :: cio
          integer           ,intent(in) :: nx,ny,nz
          real(kind=op)     ,intent(in) :: ashcon_tot(nx,ny,nz)
        end subroutine write_3D_Binary
        subroutine write_2D_Binary(nx,ny,OutVar,VarMask,Fill_Value,filename_root)
          integer,parameter  :: ip         = 8 ! Internal precision
          integer          ,intent(in) :: nx,ny
          real(kind=ip)    ,intent(in) :: OutVar(nx,ny)
          logical          ,intent(in) :: VarMask(nx,ny)
          character(len=6) ,intent(in) :: Fill_Value
          character(len=20),intent(in) :: filename_root
        end subroutine write_2D_Binary
        character (len=13) function HS_yyyymmddhh_since(HoursSince,byear,useLeaps)
          real(kind=8)               ::  HoursSince
          integer                    ::  byear
          logical                    ::  useLeaps
        end function HS_yyyymmddhh_since
        character (len=20) function HS_xmltime(HoursSince,byear,useLeaps)
          real(kind=8)              :: HoursSince
          integer                   :: byear
          logical                   :: useLeaps
        end function HS_xmltime
      END INTERFACE

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Open/Create output files
      !   This is only executed the first time output_results is called
      !   before the main time loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(first_time)then

        ! Note: 2d and 3d ASCII/binary files are written as individual files
        !       as are Airport/POI files and cloud/deposit arrival times.
        !       Profiles, KML time-series, and netcdf files need to be created
        !       and left opened throughout the run.
        iout3d = 0

        ! Open vertical profiles files
        if (Write_PR_Data) call vprofileopener

        ! Open KML files
        call Set_OutVar_Specs          ! Initialize variables local to the Output_KML module
        if (WriteCloudConcentration_KML)  call OpenFile_KML(1) ! Cloud Concentration
        if (WriteCloudHeight_KML)         call OpenFile_KML(2) ! Cloud Top Height
        if (WriteCloudLoad_KML)           call OpenFile_KML(4) ! Cloud Load
        if (WriteReflectivity_KML)        call OpenFile_KML(6) ! Reflectivity
        if (WriteDepositTS_KML.or.WriteDepositFinal_KML)  then
          call OpenFile_KML(7) ! Deposit
          call OpenFile_KML(8) ! Deposit (NWS)
        endif

        if(Write3dFiles)then
          if(ioutputFormat.eq.3)then
#ifdef USENETCDF
            ! Create the netcdf file and define dimensions/variables
            call NC_create_netcdf_file
#else
            do io=1,2;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"ERROR: Ash3d was not compiled with netcdf libraries, but netcdf"
              write(errlog(io),*)"       output is requested.  Please recompile Ash3d with"
              write(errlog(io),*)"       USENETCDF=T, or select another output format."
            endif;enddo
            stop 1
#endif
          endif
        endif

        ! This is an not output step, so toggle off the marker for the output log
        OutputStep_Marker = ' '

        first_time = .false.

        if(NextWriteTime.gt.EPS_SMALL)then
          ! If the first output timestep is essentially t=0, then continue in
          !  this subroutine, else return to Ash3d.F90
          return
        endif

      endif  ! first_time
      ! Finished opening all the netcdf and kml output files.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !construct text string for timespan written to KML files
      if (iTimeNext.gt.0) then
         if (iTimeNext.eq.1) then
             !if this is the first time, set timestart equal to the eruption time
             timestart = SimStartHour
           else
             !otherwise, set it equal to the midpoint between now and the last write time
             timestart = SimStartHour + (WriteTimes(iTimeNext-1)+WriteTimes(iTimeNext))/2.0_ip
         end if
         if (iTimeNext.lt.nWriteTimes) then
             !Set timeend to the midpoint between now and the next write time
             timeend = SimStartHour + (WriteTimes(iTimeNext)+WriteTimes(iTimeNext+1))/2.0_ip
           else
             !If this is the last write time, set timeend to the end of the simulation
             timeend = SimStartHour + Simtime_in_hours
         endif
         xmlTimeSpanStart = HS_xmltime(timestart,BaseYear,useLeap)
         xmlTimeSpanEnd   = HS_xmltime(timeend,BaseYear,useLeap)
      endif

      !increment iTimeNext
      iout3d = iout3d + 1    ! increment the counter for the output step
      if (iTimeNext.lt.nWriteTimes) then   !adjust next write time
         iTimeNext = iTimeNext + 1
         NextWriteTime = WriteTimes(iTimeNext)
        else
          ! Otherwise, set the next time to very large number
         NextWriteTime = 1.0_ip/EPS_SMALL
      endif

      !get data string of current date & time
      if (isFinal_TS) then
        cio='________final'
      else
        cio = HS_yyyymmddhh_since(SimStartHour+time,BaseYear,useLeap)
      endif
      ! Write time-series data at this output step
      if (.not.isFinal_TS) then
          ! First ascii files
        if (WriteDepositTS_ASCII)            &
          call write_2D_ASCII(nxmax,nymax,&
                              DepositThickness(1:nxmax,1:nymax), &
                              Mask_Deposit(1:nxmax,1:nymax),&
                              ' 0.000','DepositFile_        ')
        if (WriteCloudConcentration_ASCII)  &
          call write_2D_ASCII(nxmax,nymax,&
                              MaxConcentration(1:nxmax,1:nymax), &
                              Mask_Cloud(1:nxmax,1:nymax),&
                              ' 0.000','CloudConcentration_ ')
        if (WriteCloudHeight_ASCII)         &
          call write_2D_ASCII(nxmax,nymax,&
                              MaxHeight(1:nxmax,1:nymax), &
                              Mask_Cloud(1:nxmax,1:nymax),&
                              ' 0.000','CloudHeight_        ')
        if (WriteCloudLoad_ASCII)           &
          call write_2D_ASCII(nxmax,nymax,&
                              CloudLoad(1:nxmax,1:nymax), &
                              Mask_Cloud(1:nxmax,1:nymax),&
                              ' 0.000','CloudLoad_          ')

          ! Now KML files
        if (WriteCloudConcentration_KML)   call Write_2D_KML(1,MaxConcentration,1,1) ! Cloud Concentration
        if (WriteCloudHeight_KML)          call Write_2D_KML(2,MaxHeight, 1,1)       ! Cloud Top Height
        if (WriteCloudLoad_KML)            call Write_2D_KML(4,CloudLoad, 1,1)       ! Cloud Load
        if (WriteReflectivity_KML)         call Write_2D_KML(6,dbZCol, 1,1)          ! Reflectivity
        if (WriteDepositTS_KML.or.WriteDepositFinal_KML)  then
          call Write_2D_KML(7,DepositThickness,0,1)         ! Deposit (mm)
          call Write_2D_KML(8,DepositThickness*MM_2_IN,0,1) ! Deposit (NWS thresholds in inches)
        endif

      endif !.not.isFinal_TS

      !************************************************************************
      !  WRITE OUT 3D CONCENTRATION FILES      
      if (Write3dFiles) then
        if(ioutputFormat.eq.1)then
          call write_3D_ASCII(cio)
        elseif(ioutputFormat.eq.2)then
          allocate(ashcon_tot(nxmax,nymax,nzmax))
          ashcon_tot = 0.0_op
          call AshTotalCalculator
          call write_3D_Binary(cio,nxmax,nymax,nzmax,ashcon_tot)
          deallocate(ashcon_tot)
          call write_2D_Binary(nxmax,nymax,&
                              DepositThickness(1:nxmax,1:nymax), &
                              Mask_Deposit(1:nxmax,1:nymax),&
                              ' 0.000','DepositFile_        ')
        elseif(ioutputFormat.eq.3)then
#ifdef USENETCDF
          call NC_append_to_netcdf
#endif
        endif
      endif
      !************************************************************************

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Write final data and close output files
      !   This is only executed the last time output_results is called
      !   after the main time loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (isFinal_TS) then

        ! close output files for vertical profiles
        if (Write_PR_Data) call vprofilecloser

        !if files of deposit arrival time are to be written out
        if (WriteDepositTime_KML) then
          call OpenFile_KML(9) ! Deposit Arrival Time
          call Write_2D_KML(9,DepArrivalTime,0,0) ! Deposit
          call Close_KML(9,0)
        endif
        if (WriteDepositTime_ASCII)     &
          call write_2D_ASCII(nxmax,nymax,&
                              DepArrivalTime(1:nxmax,1:nymax), &
                              Mask_Deposit(1:nxmax,1:nymax),&
                              '-1.000','DepositArrivalTime  ')
        if (WriteCloudTime_KML) then
          call OpenFile_KML(5) ! Cloud Arrival Time
          call Write_2D_KML(5,CloudArrivalTime,0,0) ! Deposit
          call Close_KML(5,0)
        endif
        if (WriteCloudTime_ASCII)       &
          call write_2D_ASCII(nxmax,nymax,&
                              CloudArrivalTime(1:nxmax,1:nymax), &
                              Mask_Cloud(1:nxmax,1:nymax),&
                              '-1.000','CloudArrivalTime    ')
        ! Write Final deposit file
        if ((WriteDepositFinal_ASCII).and.                       &
            (((nWriteTimes.gt.0).and.(abs(time-dt-WriteTimes(nWriteTimes)).gt.1.0e-4_ip)).or. &
            (nWriteTimes.eq.0))) then
          call write_2D_ASCII(nxmax,nymax,&
                              DepositThickness(1:nxmax,1:nymax), &
                              Mask_Deposit(1:nxmax,1:nymax),&
                              ' 0.000','DepositFile_        ')
        endif
        if (WriteDepositFinal_KML) then
          call Write_2D_KML(7,DepositThickness,0,0) ! Deposit
          call Write_2D_KML(8,DepositThickness*MM_2_IN,0,0) ! Deposit (NWS)
        endif

        ! Close KML files
        if (WriteCloudConcentration_KML) call Close_KML(1,1)
        if (WriteCloudHeight_KML)        call Close_KML(2,1)
        if (WriteCloudLoad_KML)          call Close_KML(4,1)
        if (WriteReflectivity_KML)       call Close_KML(6,1)
        if (WriteDepositTS_KML) then
            call Close_KML(7,1)
            call Close_KML(8,1)
          else if (WriteDepositFinal_KML) then
            call Close_KML(7,0)
            call Close_KML(8,0)
        endif

        ! Airport / Point-of-interest files
        if (WriteAirportFile_KML)   call Write_PointData_Airports_KML
        if (WriteAirportFile_ASCII) call Write_PointData_Airports_ASCII

      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! This is an output step, so toggle on the marker for the output log
      OutputStep_Marker = '*'

      return

      end subroutine output_results

!******************************************************************************
