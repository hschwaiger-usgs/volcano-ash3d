      subroutine output_results

      use precis_param

      use io_units

      use global_param,  only : &
         EPS_SMALL

      use io_data,       only : &
         io,ioutputFormat,WriteTimes,nWriteTimes,isFinal_TS,&
         NextWriteTime,nTimeNext,nvprofiles,OutputStep_Marker,&
         Write3dFiles,WriteAirportFile_ASCII,&
         WriteCloudConcentration_ASCII,WriteAirportFile_KML,&
         WriteCloudConcentration_KML,WriteCloudHeight_ASCII,&
         WriteCloudHeight_KML,WriteCloudLoad_ASCII,WriteCloudLoad_KML,&
         WriteCloudTime_ASCII,WriteCloudTime_KML,&
         WriteDepositFinal_ASCII,WriteDepositTS_KML,WriteDepositFinal_KML,&
         WriteDepositTime_ASCII,WriteDepositTime_KML,WriteDepositTS_ASCII,&
         Writereflectivity_KML

      use mesh,          only : &
         nxmax,nymax,nsmax

      use Output_Vars,   only : &
         DepositThickness,MaxConcentration,MaxHeight,MinHeight,&
         CloudLoad,DepArrivalTime,CloudArrivalTime,dbZCol

      use time_data,     only : &
         BaseYear,useLeap,dt,time,SimStartHour,Simtime_in_hours,&
         xmlTimeSpanStart,xmlTimeSpanEnd

      use Output_KML,    only : &
           OpenFile_KML,&
           Close_KML,&
           Write_2D_KML,&
           Write_PointData_Airports_KML,&
           Set_OutVar_Specs

      implicit none

      character(len=13) :: cio
      character(len=13) :: HS_yyyymmddhh_since
      character(len=20) :: HS_xmltime
      real(kind=ip)     :: timestart
      real(kind=ip)     :: timeend
      logical,save      :: first_time = .true.

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

        ! Open vertical profiles files
        if (nvprofiles.gt.0) call vprofileopener

        ! Open KML files
        call Set_OutVar_Specs          ! Initialize variables local to the Output_KLM module
        if (WriteCloudConcentration_KML)  call OpenFile_KML(1) ! Cloud Concentration
        if (WriteCloudHeight_KML)         call OpenFile_KML(2) ! Cloud Top Height
        if (WriteCloudHeight_KML)         call OpenFile_KML(3) ! Cloud Bottom Height
        if (WriteCloudLoad_KML)           call OpenFile_KML(4) ! Cloud Load
        if (Writereflectivity_KML)        call OpenFile_KML(6) ! Reflectivity
        if (WriteDepositTS_KML.or.WriteDepositFinal_KML)  then
             call OpenFile_KML(7) ! Deposit
             call OpenFile_KML(8) ! Deposit (NWS)
        endif

        if(Write3dFiles)then
          if(ioutputFormat.eq.3)then
#ifdef USENETCDF
            ! Create the netcdf file and define dimensions/variables
            call create_netcdf_file
#else
            write(global_info,*)"ERROR: Ash3d was not compiled with netcdf libraries, but netcdf"
            write(global_info,*)"       output is requested.  Please recompile Ash3d with"
            write(global_info,*)"       USENETCDF=T, or select another output format."
            stop 1
#endif
          endif
        endif

        OutputStep_Marker = ' '

        first_time = .false.

        !if(time+dt.lt.NextWriteTime)then
        if(NextWriteTime.gt.EPS_SMALL)then
          ! If the first output timestep is essentially t=0, then continue in
          !  this subroutine, else return to Ash3d.F90
          return
        endif

      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      io = io+1
      if (nTimeNext.lt.nWriteTimes) then   !adjust next write time
         nTimeNext = nTimeNext + 1
         NextWriteTime = WriteTimes(nTimeNext)
        else
         NextWriteTime = 1.0e10_ip
      endif
      !construct text string for timespan written to KML files
      if (nTimeNext.gt.0) then
         !timenow = SimStartHour+time
         !xmlTimeNow = HS_xmltime(timenow,BaseYear,useLeap)
         timestart = SimStartHour + (WriteTimes(nTimeNext-1)+WriteTimes(nTimeNext))/2.0_ip
         if (nTimeNext.lt.nWriteTimes) then
             timeend = SimStartHour + (WriteTimes(nTimeNext+1)+WriteTimes(nTimeNext))/2.0_ip
           else
             timeend = SimStartHour + Simtime_in_hours
         endif
         xmlTimeSpanStart = HS_xmltime(timestart,BaseYear,useLeap)
         xmlTimeSpanEnd   = HS_xmltime(timeend,BaseYear,useLeap)
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
          call write_2D_ASCII(DepositThickness(1:nxmax,1:nymax), &
                              ' 0.000','DepositFile_        ')
        if (WriteCloudConcentration_ASCII)  &
          call write_2D_ASCII(MaxConcentration(1:nxmax,1:nymax)/1.0e3_ip, &
                              ' 0.000','CloudConcentration_ ')
        if (WriteCloudHeight_ASCII)         &
          call write_2D_ASCII(MaxHeight(1:nxmax,1:nymax), &
                              '-9999.','CloudHeight_        ')
        if (WriteCloudHeight_ASCII)         &
          call write_2D_ASCII(MinHeight(1:nxmax,1:nymax), &
                              '-9999.','CloudHeightBot_     ')
        if (WriteCloudLoad_ASCII)           &
          call write_2D_ASCII(CloudLoad(1:nxmax,1:nymax), &
                              ' 0.000','CloudLoad_          ')

          ! Now KML files
        if (WriteCloudConcentration_KML)   call Write_2D_KML(1,MaxConcentration/1.0e3_ip,1,1) ! Cloud Concentration
        if (WriteCloudHeight_KML)          call Write_2D_KML(2,MaxHeight, 1,1) ! Cloud Top Height
        if (WriteCloudHeight_KML)          call Write_2D_KML(3,MinHeight,-1,1) ! Cloud Bottom Height
        if (WriteCloudLoad_KML)            call Write_2D_KML(4,CloudLoad, 1,1) ! Cloud Load
        if (Writereflectivity_KML)         call Write_2D_KML(6,dbZCol, 1,1) ! Reflectivity
        !if (WriteDepositTS_KML.or.WriteDepositFinal_KML)  then
        if (WriteDepositTS_KML)  then
             call Write_2D_KML(7,DepositThickness,0,1) ! Deposit
             call Write_2D_KML(8,DepositThickness,0,1) ! Deposit (NWS)
        endif

      endif

      !************************************************************************
      !  WRITE OUT 3D CONCENTRATION FILES      
      if (Write3dFiles) then
        if(ioutputFormat.eq.1)then
          call write_3D_ASCII(cio)
        elseif(ioutputFormat.eq.2)then
          call write_3D_Binary(cio,nsmax)
        elseif(ioutputFormat.eq.3)then
#ifdef USENETCDF
          call append_to_netcdf
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
        if (nvprofiles.gt.0) call vprofilecloser

        !if files of deposit arrival time are to be written out
        if (WriteDepositTime_KML) then
          call OpenFile_KML(9) ! Deposit Arrival Time
          call Write_2D_KML(9,DepArrivalTime,0,0) ! Deposit
          call Close_KML(9,0)
        endif
        if (WriteDepositTime_ASCII)     &
          call write_2D_ASCII(DepArrivalTime(1:nxmax,1:nymax), &
                                '-1.000','DepositArrivalTime  ')
        if (WriteCloudTime_KML) then
          call OpenFile_KML(5) ! Cloud Arrival Time
          call Write_2D_KML(5,CloudArrivalTime,0,0) ! Deposit
          call Close_KML(5,0)
        endif
        if (WriteCloudTime_ASCII)       &
          call write_2D_ASCII(CloudArrivalTime(1:nxmax,1:nymax), &
                                '-1.000','CloudArrivalTime    ')
        ! Write Final deposit file
        if ((WriteDepositFinal_ASCII).and.                       &
            (((nWriteTimes.gt.0).and.(abs(time-dt-WriteTimes(nWriteTimes)).gt.1.0e-4_ip)).or. &
            (nWriteTimes.eq.0))) then
                 call write_2D_ASCII(DepositThickness(1:nxmax,1:nymax), &
                                ' 0.000','DepositFile_        ')
        endif
        if (WriteDepositFinal_KML) then
          call Write_2D_KML(7,DepositThickness,0,0) ! Deposit
          call Write_2D_KML(8,DepositThickness,0,0) ! Deposit (NWS)
        endif

        ! Close KML files
        if (WriteCloudConcentration_KML) call Close_KML(1,1)
        if (WriteCloudHeight_KML)        call Close_KML(2,1)
        if (WriteCloudHeight_KML)        call Close_KML(3,1)
        if (WriteCloudLoad_KML)          call Close_KML(4,1)
        if (Writereflectivity_KML)       call Close_KML(6,1)
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

      OutputStep_Marker = '*'

      return

      end subroutine output_results

!******************************************************************************
