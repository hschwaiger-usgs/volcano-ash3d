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

      use time_data,     only : &
         BaseYear,useLeap,dt,time,SimStartHour,Simtime_in_hours,&
         xmlTimeSpanStart,xmlTimeSpanEnd

      use Output_Vars,   only : &
         DepositThickness,MaxConcentration,MaxHeight,&
         Mask_Cloud,Mask_Deposit,ashcon_tot, &
         CloudLoad,DepArrivalTime,CloudArrivalTime,dbZCol,&
         AshTotalCalculator

      use Ash3d_ASCII_IO,  only : &
           write_2D_ASCII, &
           write_3D_ASCII, &
           vprofileopener, &
           vprofilecloser, &
           Write_PointData_Airports_ASCII

      use Ash3d_Binary_IO, only : &
           write_2D_Binary, &
           write_3D_Binary

      use Ash3d_KML_IO,  only : &
           OpenFile_KML,&
           Close_KML,&
           Write_2D_KML,&
           Write_PointData_Airports_KML,&
           Set_OutVar_Specs

#ifdef USENETCDF
      use Ash3d_Netcdf_IO
#endif

      implicit none

      character(len=13) :: cio
      real(kind=8)      :: timestart
      real(kind=8)      :: timeend
      logical           :: Mask(nxmax,nymax)
      integer           :: i,j
      logical,save      :: first_time = .true.

      INTERFACE
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

      ! Construct text string for timespan written to KML files
      if (iTimeNext.gt.0) then
         if (iTimeNext.eq.1) then
             ! If this is the first time, set timestart equal to the eruption time
             timestart = SimStartHour
           else
             ! otherwise, set it equal to the midpoint between now and the last write time
             timestart = SimStartHour + (WriteTimes(iTimeNext-1)+WriteTimes(iTimeNext))/2.0_ip
         end if
         if (iTimeNext.lt.nWriteTimes) then
             ! Set timeend to the midpoint between now and the next write time
             timeend = SimStartHour + (WriteTimes(iTimeNext)+WriteTimes(iTimeNext+1))/2.0_ip
           else
             ! If this is the last write time, set timeend to the end of the simulation
             timeend = SimStartHour + Simtime_in_hours
         endif
         xmlTimeSpanStart = HS_xmltime(timestart,BaseYear,useLeap)
         xmlTimeSpanEnd   = HS_xmltime(timeend,BaseYear,useLeap)
      endif

      ! increment counter for the output step (iTimeNext)
      iout3d = iout3d + 1
      if (iTimeNext.lt.nWriteTimes) then   ! adjust next write time
         iTimeNext = iTimeNext + 1
         NextWriteTime = WriteTimes(iTimeNext)
        else
          ! Otherwise, set the next time to very large number
         NextWriteTime = 1.0_ip/EPS_SMALL
      endif

      ! get data string of current date and time
      if (isFinal_TS) then
        cio='________final'
      else
        cio = HS_yyyymmddhh_since(SimStartHour+time,BaseYear,useLeap)
      endif
      ! Write time-series data at this output step
      if (.not.isFinal_TS) then
          ! First ascii files
        if (WriteDepositTS_ASCII)then
          Mask(1:nxmax,1:nymax) = .true.
          call write_2D_ASCII(nxmax,nymax,                       &
                              DepositThickness(1:nxmax,1:nymax), &
                              Mask(1:nxmax,1:nymax),&
                              '-9999.','DepositFile_        ')
        elseif (WriteCloudConcentration_ASCII)then
          Mask(1:nxmax,1:nymax) = .true.
          call write_2D_ASCII(nxmax,nymax,                       &
                              MaxConcentration(1:nxmax,1:nymax), &
                              Mask(1:nxmax,1:nymax),             &
                              '-9999.','CloudConcentration_ ')
        elseif (WriteCloudHeight_ASCII)then
          Mask(1:nxmax,1:nymax) = Mask_Cloud(1:nxmax,1:nymax)
          call write_2D_ASCII(nxmax,nymax,                &
                              MaxHeight(1:nxmax,1:nymax), &
                              Mask(1:nxmax,1:nymax),      &
                              '-9999.','CloudHeight_        ')
        elseif (WriteCloudLoad_ASCII)then
          Mask(1:nxmax,1:nymax) = .true.
          call write_2D_ASCII(nxmax,nymax,                &
                              CloudLoad(1:nxmax,1:nymax), &
                              Mask(1:nxmax,1:nymax),      &
                              '-9999.','CloudLoad_          ')
        endif

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
          allocate(ashcon_tot(nxmax,nymax,nzmax))
          ashcon_tot = 0.0_op
          call AshTotalCalculator
          call write_3D_ASCII(cio,nxmax,nymax,nzmax,ashcon_tot)
          deallocate(ashcon_tot)
        elseif(ioutputFormat.eq.2)then
          allocate(ashcon_tot(nxmax,nymax,nzmax))
          ashcon_tot = 0.0_op
          call AshTotalCalculator
          call write_3D_Binary(cio,nxmax,nymax,nzmax,ashcon_tot)
          deallocate(ashcon_tot)
          Mask(1:nxmax,1:nymax) = .true.
          call write_2D_Binary(nxmax,nymax,                      &
                              DepositThickness(1:nxmax,1:nymax), &
                              Mask(1:nxmax,1:nymax),             &
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
          call Write_2D_KML(9,real(DepArrivalTime,kind=ip),0,0) ! Deposit
          call Close_KML(9,0)
        endif
        if (WriteDepositTime_ASCII)then
          Mask(1:nxmax,1:nymax) = Mask_Deposit(1:nxmax,1:nymax)
          call write_2D_ASCII(nxmax,nymax,                                   &
                              real(DepArrivalTime(1:nxmax,1:nymax),kind=ip), &
                              Mask(1:nxmax,1:nymax),                         &
                              '-9999.','DepositArrivalTime  ')
        endif
        if (WriteCloudTime_KML) then
          call OpenFile_KML(5) ! Cloud Arrival Time
          call Write_2D_KML(5,real(CloudArrivalTime,kind=ip),0,0) ! Deposit
          call Close_KML(5,0)
        endif
        if (WriteCloudTime_ASCII)then
          ! cloud mask based on cloud load does not work in this case the cloud load mask
          ! is a function of time
          Mask(1:nxmax,1:nymax) = .true.
          do i=1,nxmax
            do j=1,nymax
              if(CloudArrivalTime(i,j).lt.0.0_ip)Mask(i,j) = .false.
            enddo
          enddo
          call write_2D_ASCII(nxmax,nymax,&
                              real(CloudArrivalTime(1:nxmax,1:nymax),kind=ip), &
                              Mask(1:nxmax,1:nymax),&
                              '-9999.','CloudArrivalTime    ')
        endif
        ! Write Final deposit file
        if ((WriteDepositFinal_ASCII).and.                       &
            (((nWriteTimes.gt.0).and.&
             (abs(time-dt-WriteTimes(nWriteTimes)).gt.1.0e-4_ip)).or. &
             (nWriteTimes.eq.0))) then
          Mask(1:nxmax,1:nymax) = .true.
          call write_2D_ASCII(nxmax,nymax,&
                              DepositThickness(1:nxmax,1:nymax), &
                              Mask(1:nxmax,1:nymax),&
                              '-9999.','DepositFile_        ')
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
