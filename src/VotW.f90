!##############################################################################
!
!  VotW_ESP module
!
!      subroutine get_ESP
!      subroutine VotW_v12
!
!##############################################################################

      module VotW_ESP

      use precis_param

      use io_units

      use io_data,       only : &
         Ash3dHome

      implicit none

        ! Set everything to private by default
      private

        ! Publicly available subroutines/functions
      public get_ESP

        ! Publicly available variables
      !-- None --


      character(len=130)  :: VotWMasterFile  !Only needed if USEEXTDATA=T

      contains
      !------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  get_ESP(volc_code)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_ESP(volc_code)

      use io_data,       only : &
         VolcanoName
      use Source,        only : &
         ESP_duration,ESP_height,ESP_MassFluxRate,ESP_massfracfine,ESP_Vol

      integer,parameter  :: MAXVOLCS = 1535
      
      character(len=8),intent(in) :: volc_code

      real(kind=ip)      :: volcLat(MAXVOLCS)
      real(kind=ip)      :: volcLon(MAXVOLCS)
      integer            :: volcElev(MAXVOLCS)
      character(len=30)  :: volcLoc(MAXVOLCS)
      character(len=2 )  :: volcESP_Code(MAXVOLCS)
      character(len=8 )  :: volcID(MAXVOLCS)
      character(len=42)  :: volcName(MAXVOLCS)

      integer            :: volcESP(MAXVOLCS)

      integer            :: i,Volcano_ID

      !"http://www.volcano.si.edu/world/volcano.cfm?vnum="
      ! http://dx.doi.org/10.5479/si.GVP.VOTW4-2013

      ! Populate the Smithsonian lists
      call VotW_v12(MAXVOLCS,volcLat,volcLon,volcElev,volcLoc,volcESP_Code,volcID,volcName)

      do i = 1, MAXVOLCS
        volcName(i) = adjustl(volcName(i))
        volcLoc(i)  = adjustl(volcLoc(i))

        ! Assign ESP's from Mastin et al., JVGR 2009
        volcESP(i) = 0
        if(volcESP_Code(i).eq."M0")then
          ! M0 -- Mafic Standard (Cerro Negro, Nicaragua 4/13/1992)
          volcESP(i)  = 1
        elseif(volcESP_Code(i).eq."M1")then
          ! M1 -- Mafic Small (Etna, Italy 7/19-24/2001)
          volcESP(i) = 2
        elseif(volcESP_Code(i).eq."M2")then
          ! M2 -- Mafic Medium (Cerro Negro, Nicaragua 4/13/1992)
          volcESP(i) = 3
        elseif(volcESP_Code(i).eq."M3")then
          ! M3 -- Mafic Large (Fuego, Guatemala 10/14/1974)
          volcESP(i) = 4
        elseif(volcESP_Code(i).eq."S0")then
          ! S0 -- Silicic Standard (Spurr, USA 8/18/1992)
          volcESP(i) = 5
        elseif(volcESP_Code(i).eq."S1")then
          ! S1 -- Silicic Small (Ruapehu, New Zealand 6/17/1996)
          volcESP(i) = 6
        elseif(volcESP_Code(i).eq."S2")then
          ! S2 -- Silicic Medium (Spurr, USA 8/18/1992)
          volcESP(i) = 7
        elseif(volcESP_Code(i).eq."S3")then
          ! S3 -- Silicic Large (Mt. St. Helens, USA 5/18/1980)
          volcESP(i) = 8
        elseif(volcESP_Code(i).eq."S8")then
          ! S8 -- Silicic co-ignimbrite cloud (Mt. St. Helens, USA 5/18/1980 pre-9am)
          volcESP(i) = 9
        elseif(volcESP_Code(i).eq."S9")then
          ! S9 -- Silicic brief (Soufriere Hills, Montserrat)
          volcESP(i) = 10
        elseif(volcESP_Code(i).eq."U0")then
          ! U0 -- Submarine (No example)
          volcESP(i) = 11
        endif

      enddo

      ! Read all the eruption source parameters from the Smithsonian data
      ! base.  Now check if 'volc_code' is on the list and return the
      ! ESP values
      Volcano_ID = 0
      do i = 1,MAXVOLCS
        if(volc_code.eq.volcID(i))then
          Volcano_ID = i
          cycle
        endif
      enddo
      if (Volcano_ID.gt.0)then
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Found volcano in database"
          write(outlog(io),*)"  Volcano ID    : ",volcID(Volcano_ID)
          write(outlog(io),*)"  Name          : ",volcName(Volcano_ID)
          write(outlog(io),*)"  Location      : ",volcLoc(Volcano_ID)
          write(outlog(io),*)"  Latitude      : ",volcLat(Volcano_ID)
          write(outlog(io),*)"  Longitude     : ",volcLon(Volcano_ID)
          write(outlog(io),*)"  Elevation     : ",volcElev(Volcano_ID)
          write(outlog(io),*)"  Eruption type : ",volcESP_Code(Volcano_ID)
          write(outlog(io),*)"  "
        endif;enddo

        !Change the input name to the EPS database name
        write(VolcanoName,*)trim(volcName(Volcano_ID)),&
                            " (",volcID(Volcano_ID),")"

        if(volcESP(Volcano_ID).eq.1)then
          ! M0 -- Mafic Standard (Cerro Negro, Nicaragua 4/13/1992)
          ESP_height       = 7.0_ip
          ESP_duration     = 60.0_ip
          ESP_MassFluxRate = 1.0e5_ip
          ESP_Vol          = 1.0e-2_ip
          ESP_massfracfine = 5.0e-2_ip
        elseif(volcESP(Volcano_ID).eq.2)then
          ! M1 -- Mafic Small (Etna, Italy 7/19-24/2001)
          ESP_height       = 2.0_ip
          ESP_duration     = 100.0_ip
          ESP_MassFluxRate = 5.0e3_ip
          ESP_Vol          = 1.0e-3_ip
          ESP_massfracfine = 2.0e-2_ip
        elseif(volcESP(Volcano_ID).eq.3)then
          ! M2 -- Mafic Medium (Cerro Negro, Nicaragua 4/13/1992)
          ESP_height       = 7.0_ip
          ESP_duration     = 60.0_ip
          ESP_MassFluxRate = 1.0e5_ip
          ESP_Vol          = 1.0e-2_ip
          ESP_massfracfine = 5.0e-2_ip
        elseif(volcESP(Volcano_ID).eq.4)then
          ! M3 -- Mafic Large (Fuego, Guatemala 10/14/1974)
          ESP_height       = 10.0_ip
          ESP_duration     = 5.0_ip
          ESP_MassFluxRate = 1.0e6_ip
          ESP_Vol          = 1.7e-1_ip
          ESP_massfracfine = 1.0e-1_ip
        elseif(volcESP(Volcano_ID).eq.5)then
          ! S0 -- Silicic Standard (Spurr, USA 8/18/1992)
          ESP_height       = 11.0_ip
          ESP_duration     = 3.0_ip
          ESP_MassFluxRate = 4.0e6_ip
          ESP_Vol          = 1.5e-2_ip
          ESP_massfracfine = 4.0e-1_ip
        elseif(volcESP(Volcano_ID).eq.6)then
          ! S1 -- Silicic Small (Ruapehu, New Zealand 6/17/1996)
          volcESP(i) = 6
          ESP_height       = 5.0_ip
          ESP_duration     = 12.0_ip
          ESP_MassFluxRate = 2.0e5_ip
          ESP_Vol          = 3.0e-3_ip
          ESP_massfracfine = 1.0e-1_ip
        elseif(volcESP(Volcano_ID).eq.7)then
          ! S2 -- Silicic Medium (Spurr, USA 8/18/1992)
          ESP_height       = 11.0_ip
          ESP_duration     = 3.0_ip
          ESP_MassFluxRate = 4.0e6_ip
          ESP_Vol          = 1.5e-2_ip
          ESP_massfracfine = 4.0e-1_ip
        elseif(volcESP(Volcano_ID).eq.8)then
          ! S3 -- Silicic Large (Mt. St. Helens, USA 5/18/1980)
          ESP_height       = 15.0_ip
          ESP_duration     = 8.0_ip
          ESP_MassFluxRate = 1.0e7_ip
          ESP_Vol          = 1.5e-1_ip
          ESP_massfracfine = 5.0e-1_ip
        elseif(volcESP(Volcano_ID).eq.9)then
          ! S8 -- Silicic co-ignimbrite cloud (Mt. St. Helens, USA
          ! 5/18/1980 pre-9am)
          ESP_height       = 25.0_ip
          ESP_duration     = 0.5_ip
          ESP_MassFluxRate = 1.0e8_ip
          ESP_Vol          = 5.0e-2_ip
          ESP_massfracfine = 5.0e-1_ip
        elseif(volcESP(Volcano_ID).eq.10)then
          ! S9 -- Silicic brief (Soufriere Hills, Montserrat)
          ESP_height       = 10.0_ip
          ESP_duration     = 0.01_ip
          ESP_MassFluxRate = 3.0e6_ip
          ESP_Vol          = 3.0e-4_ip
          ESP_massfracfine = 6.0e-1_ip
        elseif(volcESP(Volcano_ID).eq.11)then
          ! U0 -- Submarine (No example)
          ESP_height       = 0.0_ip
          ESP_duration     = 0.0_ip
          ESP_MassFluxRate = 0.0_ip
          ESP_Vol          = 0.0_ip
          ESP_massfracfine = 0.0_ip
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"Not set up for submarine ESP"
          endif;enddo
          stop 1
        else
          do io=1,2;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"Could not read the eruption style."
          endif;enddo
          stop 1
        endif
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"The following ESPs will be used IF variables are not"
          write(outlog(io),*)"assigned in the input file:"
          write(outlog(io),*)"  Plume Height (km)     : ",ESP_height
          write(outlog(io),*)"  Duration (h)          : ",ESP_duration
          write(outlog(io),*)"  Mass Flux Rate (kg/s) : ",ESP_MassFluxRate
          write(outlog(io),*)"  Volume (km^3)         : ",ESP_Vol
          write(outlog(io),*)"  Mass frac of fines    : ",ESP_massfracfine
        endif;enddo
      else ! if (Volcano_ID.gt.0)
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Did not find volcano in database"
        endif;enddo
      endif

      return

      
      !FORMAT STATEMENTS
!1     format(/,5x,'Reading from ESP file')      
!2     format(49x,a35)
!5     format(5x,'Error.  cannot find input file ',a130,/,5x,&
!                'Program stopped')
!50    format(a8,a42,a31,f21.3,a18,f19.3,a18,a15,a18)

      end subroutine get_ESP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  VotW_v12(MAXVOLCS,volcLat,volcLon,volcElev,volcLoc,volcESP_Code,volcID,volcName)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef USEEXTDATA
      subroutine VotW_v12(MAXVOLCS,volcLat,volcLon,volcElev,volcLoc,volcESP_Code,volcID,volcName)

      integer           ,intent(in)    :: MAXVOLCS
      real(kind=ip)     ,intent(inout) :: volcLat(MAXVOLCS)
      real(kind=ip)     ,intent(inout) :: volcLon(MAXVOLCS)
      integer           ,intent(inout) :: volcElev(MAXVOLCS)
      character(len=30) ,intent(inout) :: volcLoc(MAXVOLCS)
      character(len=2 ) ,intent(inout) :: volcESP_Code(MAXVOLCS)
      character(len=8 ) ,intent(inout) :: volcID(MAXVOLCS)
      character(len=42) ,intent(inout) :: volcName(MAXVOLCS)
      logical              :: IsThere
      integer              :: Iostatus    = 1
      character(len=20)    ::  temp1
      character(len=20)    ::  temp2
      character(len=20)    ::  temp3
      character(len=30)    ::  volcElev_c(MAXVOLCS)
      character            ::  volcNS(MAXVOLCS)
      character            ::  volcWE(MAXVOLCS)

      character(len=195) :: linebuffer195
      integer         :: nvolcs
      !integer         :: i,Volcano_ID
      character       :: testkey

      VotWMasterFile = trim(Ash3dHome) // &
                          '/share/VotW_ESP_v12_csv.txt'
      ! Test for existance of the VotW file
      inquire( file=trim(adjustl(VotWMasterFile)), exist=IsThere )
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"     ",trim(adjustl(VotWMasterFile)),IsThere
      endif;enddo
      if(.not.IsThere)then
        do io=1,2;if(VB(io).le.verbosity_error)then          
          write(errlog(io),*)"ERROR: Could not find VotW file."
          write(errlog(io),*)"       Please copy file to this location:"
          write(errlog(io),*)VotWMasterFile
        endif;enddo
        stop 1 
      endif

      !"http://www.volcano.si.edu/world/volcano.cfm?vnum="
      open(unit=20,file=VotWMasterFile,status='old',err=3000)
      read(20,'(a195)')linebuffer195
      nvolcs = 0
      read(20,'(a195)',IOSTAT=Iostatus) linebuffer195
      do while (Iostatus.ge.0)
        nvolcs = nvolcs + 1
        read(linebuffer195,50)volcID(nvolcs),     &
                              volcName(nvolcs),   &
                              volcLoc(nvolcs),    &
                              volcLat(nvolcs),    &
                              temp1,              &
                              volcLon(nvolcs),    &
                              temp2,              &
                              volcElev_c(nvolcs), &
                              temp3
        volcName(nvolcs) = adjustl(volcName(nvolcs))
        volcLoc(nvolcs)  = adjustl(volcLoc(nvolcs))
        volcElev_c(nvolcs) = adjustl(volcElev_c(nvolcs))
        temp1              = adjustl(temp1)
        volcNS(nvolcs)     = temp1(1:1)
        temp2              = adjustl(temp2)
        volcWE(nvolcs)     = temp2(1:1)
        temp3              = adjustl(temp3)
        volcESP_Code(nvolcs) = temp3(1:1)

        read(volcElev_c(nvolcs),*)testkey
        if(testkey.eq.'U'.or.testkey.eq.'v')then
          volcElev(nvolcs) = 0
        else
          read(volcElev_c(nvolcs),*)volcElev(nvolcs)
        endif
        if(volcNS(nvolcs).eq.'S') volcLat(nvolcs) = -volcLat(nvolcs)
        if(volcWE(nvolcs).eq.'W') volcLon(nvolcs) = -volcLon(nvolcs) + 360.0_ip
        do io=1,2;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)&
             volcID(nvolcs),volcName(nvolcs),volcLoc(nvolcs),volcLat(nvolcs),volcLon(nvolcs)
        endif;enddo
        read(20,'(a195)',IOSTAT=Iostatus) linebuffer195
      enddo

      close(20)

      return

!     ERROR TRAPS
3000  do io=1,2;if(VB(io).le.verbosity_error)then          
        write(errlog(io),5) "VotW_ESP_v12.esp"
      endif;enddo
      stop 1

      !FORMAT STATEMENTS
5     format(5x,'Error.  cannot find input file ',a130,/,5x,&
                'Program stopped')
50    format(a8,a42,a31,f21.3,a18,f19.3,a18,a15,a18)


      end subroutine VotW_v12
#else
      subroutine VotW_v12(MAXVOLCS,volcLat,volcLon,volcElev,volcLoc,volcESP_Code,volcID,volcName)

      integer           ,intent(in)    :: MAXVOLCS
      real(kind=ip)     ,intent(inout) :: volcLat(MAXVOLCS)
      real(kind=ip)     ,intent(inout) :: volcLon(MAXVOLCS)
      integer           ,intent(inout) :: volcElev(MAXVOLCS)
      character(len=30) ,intent(inout) :: volcLoc(MAXVOLCS)
      character(len=2 ) ,intent(inout) :: volcESP_Code(MAXVOLCS)
      character(len=8 ) ,intent(inout) :: volcID(MAXVOLCS)
      character(len=42) ,intent(inout) :: volcName(MAXVOLCS)

      integer :: i

      i = 0;
      i=i+1; volcLat(i)= 50.170; volcLon(i)=  6.850; volcElev(i)=  600;  volcLoc(i)="Germany                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0100-01-"; volcName(i)="West Eifel Volc Field          "; 
      i=i+1; volcLat(i)= 45.775; volcLon(i)=  2.970; volcElev(i)= 1464;  volcLoc(i)="France                        "; 
             volcESP_Code(i)="M0"; volcID(i)="0100-02-"; volcName(i)="Chaîne des Puys                "; 
      i=i+1; volcLat(i)= 42.170; volcLon(i)=  2.530; volcElev(i)=  893;  volcLoc(i)="Spain                         "; 
             volcESP_Code(i)="M0"; volcID(i)="0100-03-"; volcName(i)="Olot Volc Field                "; 
      i=i+1; volcLat(i)= 38.870; volcLon(i)=355.980; volcElev(i)= 1117;  volcLoc(i)="Spain                         "; 
             volcESP_Code(i)="M0"; volcID(i)="0100-04-"; volcName(i)="Calatrava Volc Field           "; 
      i=i+1; volcLat(i)= 43.250; volcLon(i)= 10.870; volcElev(i)=  500;  volcLoc(i)="Italy                         "; 
             volcESP_Code(i)="S0"; volcID(i)="0101-001"; volcName(i)="Larderello                     "; 
      i=i+1; volcLat(i)= 42.600; volcLon(i)= 11.930; volcElev(i)=  800;  volcLoc(i)="Italy                         "; 
             volcESP_Code(i)="S0"; volcID(i)="0101-003"; volcName(i)="Vulsini                        "; 
      i=i+1; volcLat(i)= 41.730; volcLon(i)= 12.700; volcElev(i)=  949;  volcLoc(i)="Italy                         "; 
             volcESP_Code(i)="S0"; volcID(i)="0101-004"; volcName(i)="Alban Hills                    "; 
      i=i+1; volcLat(i)= 40.827; volcLon(i)= 14.139; volcElev(i)=  458;  volcLoc(i)="Italy                         "; 
             volcESP_Code(i)="S0"; volcID(i)="0101-01="; volcName(i)="Campi Flegrei                  "; 
      i=i+1; volcLat(i)= 40.821; volcLon(i)= 14.426; volcElev(i)= 1281;  volcLoc(i)="Italy                         "; 
             volcESP_Code(i)="S2"; volcID(i)="0101-02="; volcName(i)="Vesuvius                       "; 
      i=i+1; volcLat(i)= 40.730; volcLon(i)= 13.897; volcElev(i)=  789;  volcLoc(i)="Italy                         "; 
             volcESP_Code(i)="S0"; volcID(i)="0101-03="; volcName(i)="Ischia                         "; 
      i=i+1; volcLat(i)= 38.630; volcLon(i)= 15.070; volcElev(i)=  421;  volcLoc(i)="Italy                         "; 
             volcESP_Code(i)="S0"; volcID(i)="0101-041"; volcName(i)="Panarea                        "; 
      i=i+1; volcLat(i)= 38.480; volcLon(i)= 14.950; volcElev(i)=  602;  volcLoc(i)="Italy                         "; 
             volcESP_Code(i)="S0"; volcID(i)="0101-042"; volcName(i)="Lipari                         "; 
      i=i+1; volcLat(i)= 38.789; volcLon(i)= 15.213; volcElev(i)=  924;  volcLoc(i)="Italy                         "; 
             volcESP_Code(i)="M1"; volcID(i)="0101-04="; volcName(i)="Stromboli                      "; 
      i=i+1; volcLat(i)= 38.404; volcLon(i)= 14.962; volcElev(i)=  500;  volcLoc(i)="Italy                         "; 
             volcESP_Code(i)="S2"; volcID(i)="0101-05="; volcName(i)="Vulcano                        "; 
      i=i+1; volcLat(i)= 37.734; volcLon(i)= 15.004; volcElev(i)= 3330;  volcLoc(i)="Italy                         "; 
             volcESP_Code(i)="M1"; volcID(i)="0101-06="; volcName(i)="Etna                           "; 
      i=i+1; volcLat(i)= 36.770; volcLon(i)= 12.020; volcElev(i)=  836;  volcLoc(i)="Italy                         "; 
             volcESP_Code(i)="S0"; volcID(i)="0101-071"; volcName(i)="Pantelleria                    "; 
      i=i+1; volcLat(i)= 37.100; volcLon(i)= 12.700; volcElev(i)=   -8;  volcLoc(i)="Italy                         "; 
             volcESP_Code(i)="M0"; volcID(i)="0101-07="; volcName(i)="Campi Flegrei Mar Sicilia      "; 
      i=i+1; volcLat(i)= 37.615; volcLon(i)= 23.336; volcElev(i)=  760;  volcLoc(i)="Greece                        "; 
             volcESP_Code(i)="S0"; volcID(i)="0102-02="; volcName(i)="Methana                        "; 
      i=i+1; volcLat(i)= 36.699; volcLon(i)= 24.439; volcElev(i)=  751;  volcLoc(i)="Greece                        "; 
             volcESP_Code(i)="S0"; volcID(i)="0102-03="; volcName(i)="Mílos                          "; 
      i=i+1; volcLat(i)= 36.404; volcLon(i)= 25.396; volcElev(i)=  367;  volcLoc(i)="Greece                        "; 
             volcESP_Code(i)="S1"; volcID(i)="0102-04="; volcName(i)="Santorini                      "; 
      i=i+1; volcLat(i)= 36.671; volcLon(i)= 27.140; volcElev(i)=  180;  volcLoc(i)="Greece                        "; 
             volcESP_Code(i)="S0"; volcID(i)="0102-051"; volcName(i)="Yali                           "; 
      i=i+1; volcLat(i)= 36.586; volcLon(i)= 27.160; volcElev(i)=  698;  volcLoc(i)="Greece                        "; 
             volcESP_Code(i)="S0"; volcID(i)="0102-05="; volcName(i)="Nisyros                        "; 
      i=i+1; volcLat(i)= 38.580; volcLon(i)= 28.520; volcElev(i)=  750;  volcLoc(i)="Turkey                        "; 
             volcESP_Code(i)="M0"; volcID(i)="0103-00-"; volcName(i)="Kula                           "; 
      i=i+1; volcLat(i)= 37.670; volcLon(i)= 33.650; volcElev(i)= 1302;  volcLoc(i)="Turkey                        "; 
             volcESP_Code(i)="M0"; volcID(i)="0103-001"; volcName(i)="Karapinar Field                "; 
      i=i+1; volcLat(i)= 38.130; volcLon(i)= 34.170; volcElev(i)= 3253;  volcLoc(i)="Turkey                        "; 
             volcESP_Code(i)="S0"; volcID(i)="0103-002"; volcName(i)="Hasan Dagi                     "; 
      i=i+1; volcLat(i)= 38.250; volcLon(i)= 34.570; volcElev(i)= 2143;  volcLoc(i)="Turkey                        "; 
             volcESP_Code(i)="S0"; volcID(i)="0103-003"; volcName(i)="Göllü Dag                      "; 
      i=i+1; volcLat(i)= 38.570; volcLon(i)= 34.520; volcElev(i)= 1689;  volcLoc(i)="Turkey                        "; 
             volcESP_Code(i)="S0"; volcID(i)="0103-004"; volcName(i)="Acigöl-Nevsehir                "; 
      i=i+1; volcLat(i)= 37.670; volcLon(i)= 39.830; volcElev(i)= 1957;  volcLoc(i)="Turkey                        "; 
             volcESP_Code(i)="M0"; volcID(i)="0103-011"; volcName(i)="Karaca Dag                     "; 
      i=i+1; volcLat(i)= 38.520; volcLon(i)= 35.480; volcElev(i)= 3916;  volcLoc(i)="Turkey                        "; 
             volcESP_Code(i)="S0"; volcID(i)="0103-01="; volcName(i)="Erciyes Dagi                   "; 
      i=i+1; volcLat(i)= 38.920; volcLon(i)= 42.820; volcElev(i)= 4158;  volcLoc(i)="Turkey                        "; 
             volcESP_Code(i)="S0"; volcID(i)="0103-021"; volcName(i)="Süphan Dagi                    "; 
      i=i+1; volcLat(i)= 39.170; volcLon(i)= 43.330; volcElev(i)=    0;  volcLoc(i)="Turkey                        "; 
             volcESP_Code(i)="S0"; volcID(i)="0103-022"; volcName(i)="Girekol                        "; 
      i=i+1; volcLat(i)= 38.650; volcLon(i)= 42.230; volcElev(i)= 2948;  volcLoc(i)="Turkey                        "; 
             volcESP_Code(i)="S0"; volcID(i)="0103-02="; volcName(i)="Nemrut Dagi                    "; 
      i=i+1; volcLat(i)= 39.370; volcLon(i)= 43.870; volcElev(i)= 3584;  volcLoc(i)="Turkey                        "; 
             volcESP_Code(i)="M0"; volcID(i)="0103-03="; volcName(i)="Tendürek Dagi                  "; 
      i=i+1; volcLat(i)= 39.700; volcLon(i)= 44.300; volcElev(i)= 5165;  volcLoc(i)="Turkey                        "; 
             volcESP_Code(i)="S0"; volcID(i)="0103-04-"; volcName(i)="Ararat                         "; 
      i=i+1; volcLat(i)= 40.750; volcLon(i)= 42.900; volcElev(i)= 3000;  volcLoc(i)="Turkey                        "; 
             volcESP_Code(i)="S0"; volcID(i)="0103-05-"; volcName(i)="Kars Plateau                   "; 
      i=i+1; volcLat(i)= 43.330; volcLon(i)= 42.450; volcElev(i)= 5633;  volcLoc(i)="Russia-SW                     "; 
             volcESP_Code(i)="S0"; volcID(i)="0104-01-"; volcName(i)="Elbrus                         "; 
      i=i+1; volcLat(i)= 42.700; volcLon(i)= 44.500; volcElev(i)= 5050;  volcLoc(i)="Georgia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0104-02-"; volcName(i)="Kasbek                         "; 
      i=i+1; volcLat(i)= 42.550; volcLon(i)= 44.000; volcElev(i)= 3650;  volcLoc(i)="Georgia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0104-03-"; volcName(i)="Kabargin Oth Group             "; 
      i=i+1; volcLat(i)= 42.450; volcLon(i)= 44.250; volcElev(i)= 3750;  volcLoc(i)="Georgia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0104-04-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 41.550; volcLon(i)= 43.600; volcElev(i)= 3400;  volcLoc(i)="Georgia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0104-05-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 40.530; volcLon(i)= 44.200; volcElev(i)= 4095;  volcLoc(i)="Armenia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0104-06-"; volcName(i)="Aragats                        "; 
      i=i+1; volcLat(i)= 40.275; volcLon(i)= 44.750; volcElev(i)= 3597;  volcLoc(i)="Armenia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0104-07-"; volcName(i)="Ghegam Ridge                   "; 
      i=i+1; volcLat(i)= 39.700; volcLon(i)= 45.542; volcElev(i)= 3329;  volcLoc(i)="Armenia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0104-08-"; volcName(i)="Dar-Alages                     "; 
      i=i+1; volcLat(i)= 40.020; volcLon(i)= 45.780; volcElev(i)= 2800;  volcLoc(i)="Armenia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0104-09-"; volcName(i)="Porak                          "; 
      i=i+1; volcLat(i)= 39.730; volcLon(i)= 46.020; volcElev(i)= 3000;  volcLoc(i)="Armenia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0104-10-"; volcName(i)="Tskhouk-Karckar                "; 
      i=i+1; volcLat(i)= 15.550; volcLon(i)= 41.830; volcElev(i)=  244;  volcLoc(i)="Red Sea                       "; 
             volcESP_Code(i)="M1"; volcID(i)="0201-01="; volcName(i)="Tair, Jebel at                 "; 
      i=i+1; volcLat(i)= 14.020; volcLon(i)= 42.750; volcElev(i)=  624;  volcLoc(i)="Red Sea                       "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-021"; volcName(i)="Zukur                          "; 
      i=i+1; volcLat(i)= 13.720; volcLon(i)= 42.730; volcElev(i)=  422;  volcLoc(i)="Red Sea                       "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-022"; volcName(i)="Hanish                         "; 
      i=i+1; volcLat(i)= 15.050; volcLon(i)= 42.180; volcElev(i)=  191;  volcLoc(i)="Red Sea                       "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-02="; volcName(i)="Zubair, Jebel                  "; 
      i=i+1; volcLat(i)= 15.042; volcLon(i)= 39.820; volcElev(i)=  713;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-03="; volcName(i)="Jalua                          "; 
      i=i+1; volcLat(i)= 14.242; volcLon(i)= 40.300; volcElev(i)=  -48;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-041"; volcName(i)="Dallol                         "; 
      i=i+1; volcLat(i)= 14.880; volcLon(i)= 39.920; volcElev(i)=  904;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-04="; volcName(i)="Alid                           "; 
      i=i+1; volcLat(i)= 13.975; volcLon(i)= 40.408; volcElev(i)=  287;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-05="; volcName(i)="Gada Ale                       "; 
      i=i+1; volcLat(i)= 13.825; volcLon(i)= 40.508; volcElev(i)=  429;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-06="; volcName(i)="Alu                            "; 
      i=i+1; volcLat(i)= 13.725; volcLon(i)= 40.600; volcElev(i)=  668;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-071"; volcName(i)="Borale Ale                     "; 
      i=i+1; volcLat(i)= 13.792; volcLon(i)= 40.550; volcElev(i)=  613;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-07="; volcName(i)="Dalaffilla                     "; 
      i=i+1; volcLat(i)= 13.600; volcLon(i)= 40.670; volcElev(i)=  613;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M1"; volcID(i)="0201-08="; volcName(i)="Erta Ale                       "; 
      i=i+1; volcLat(i)= 13.500; volcLon(i)= 40.720; volcElev(i)=  521;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-091"; volcName(i)="Hayli Gubbi                    "; 
      i=i+1; volcLat(i)= 13.520; volcLon(i)= 40.630; volcElev(i)= 1031;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-09="; volcName(i)="Ale Bagu                       "; 
      i=i+1; volcLat(i)= 13.370; volcLon(i)= 41.700; volcElev(i)= 2218;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-101"; volcName(i)="Nabro                          "; 
      i=i+1; volcLat(i)= 13.270; volcLon(i)= 41.650; volcElev(i)= 1875;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-102"; volcName(i)="Mallahle                       "; 
      i=i+1; volcLat(i)= 13.180; volcLon(i)= 41.725; volcElev(i)= 1611;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-103"; volcName(i)="Sork Ale                       "; 
      i=i+1; volcLat(i)= 13.070; volcLon(i)= 41.600; volcElev(i)= 1200;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-104"; volcName(i)="Asavyo                         "; 
      i=i+1; volcLat(i)= 13.100; volcLon(i)= 41.150; volcElev(i)=  523;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-105"; volcName(i)="Mat Ala                        "; 
      i=i+1; volcLat(i)= 13.280; volcLon(i)= 41.070; volcElev(i)=  700;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-106"; volcName(i)="Tat Ali                        "; 
      i=i+1; volcLat(i)= 13.300; volcLon(i)= 40.980; volcElev(i)=  812;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-107"; volcName(i)="Borawli                        "; 
      i=i+1; volcLat(i)= 13.580; volcLon(i)= 41.808; volcElev(i)= 1625;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-10="; volcName(i)="Dubbi                          "; 
      i=i+1; volcLat(i)= 13.020; volcLon(i)= 40.200; volcElev(i)= 1815;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-111"; volcName(i)="Ma Alalta                      "; 
      i=i+1; volcLat(i)= 12.880; volcLon(i)= 40.570; volcElev(i)= 1501;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-112"; volcName(i)="Alayta                         "; 
      i=i+1; volcLat(i)= 12.600; volcLon(i)= 40.480; volcElev(i)= 1442;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-113"; volcName(i)="Dabbahu                        "; 
      i=i+1; volcLat(i)= 12.380; volcLon(i)= 40.070; volcElev(i)= 1302;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-114"; volcName(i)="Dabbayra                       "; 
      i=i+1; volcLat(i)= 12.170; volcLon(i)= 40.820; volcElev(i)=  600;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-115"; volcName(i)="Manda Hararo                   "; 
      i=i+1; volcLat(i)= 11.730; volcLon(i)= 40.250; volcElev(i)=  930;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-116"; volcName(i)="Groppo                         "; 
      i=i+1; volcLat(i)= 13.080; volcLon(i)= 40.850; volcElev(i)= 1295;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-11="; volcName(i)="Afderà                         "; 
      i=i+1; volcLat(i)= 11.630; volcLon(i)= 41.450; volcElev(i)=  875;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-121"; volcName(i)="Borawli                        "; 
      i=i+1; volcLat(i)= 12.380; volcLon(i)= 42.200; volcElev(i)=  600;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-122"; volcName(i)="Manda-Inakir                   "; 
      i=i+1; volcLat(i)= 12.470; volcLon(i)= 42.400; volcElev(i)= 2028;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-123"; volcName(i)="Mousa Alli                     "; 
      i=i+1; volcLat(i)= 12.550; volcLon(i)= 42.530; volcElev(i)=  600;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-124"; volcName(i)="Gufa                           "; 
      i=i+1; volcLat(i)= 12.950; volcLon(i)= 42.430; volcElev(i)=  987;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-125"; volcName(i)="Assab Volc Field               "; 
      i=i+1; volcLat(i)= 11.580; volcLon(i)= 42.470; volcElev(i)=  298;  volcLoc(i)="Djibouti                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-126"; volcName(i)="Ardoukôba                      "; 
      i=i+1; volcLat(i)= 11.880; volcLon(i)= 41.208; volcElev(i)=  625;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-12="; volcName(i)="Kurub                          "; 
      i=i+1; volcLat(i)= 11.280; volcLon(i)= 41.630; volcElev(i)= 1068;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-141"; volcName(i)="Dama Ali                       "; 
      i=i+1; volcLat(i)= 10.580; volcLon(i)= 41.042; volcElev(i)= 1383;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-151"; volcName(i)="Yangudi                        "; 
      i=i+1; volcLat(i)= 11.080; volcLon(i)= 41.270; volcElev(i)= 1459;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-15="; volcName(i)="Gabillema                      "; 
      i=i+1; volcLat(i)= 10.082; volcLon(i)= 40.702; volcElev(i)= 2145;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-16="; volcName(i)="Ayelu                          "; 
      i=i+1; volcLat(i)=  9.780; volcLon(i)= 40.330; volcElev(i)=  900;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-171"; volcName(i)="Hertali                        "; 
      i=i+1; volcLat(i)=  9.570; volcLon(i)= 40.280; volcElev(i)=  878;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-172"; volcName(i)="Liado Hayk                     "; 
      i=i+1; volcLat(i)= 10.070; volcLon(i)= 40.840; volcElev(i)= 1733;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-17="; volcName(i)="Adwa                           "; 
      i=i+1; volcLat(i)=  9.350; volcLon(i)= 40.130; volcElev(i)= 1151;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-18="; volcName(i)="Dofen                          "; 
      i=i+1; volcLat(i)=  8.950; volcLon(i)= 39.750; volcElev(i)= 1100;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-191"; volcName(i)="Beru                           "; 
      i=i+1; volcLat(i)=  8.975; volcLon(i)= 39.930; volcElev(i)= 2007;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-19="; volcName(i)="Fentale                        "; 
      i=i+1; volcLat(i)=  8.800; volcLon(i)= 39.692; volcElev(i)= 1619;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-20-"; volcName(i)="Kone                           "; 
      i=i+1; volcLat(i)=  8.700; volcLon(i)= 39.630; volcElev(i)= 1300;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-201"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=  8.558; volcLon(i)= 39.475; volcElev(i)= 2447;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-21-"; volcName(i)="Boset-Bericha                  "; 
      i=i+1; volcLat(i)=  8.780; volcLon(i)= 38.980; volcElev(i)= 1850;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-22-"; volcName(i)="Bishoftu Volc Field            "; 
      i=i+1; volcLat(i)=  8.620; volcLon(i)= 38.950; volcElev(i)= 1800;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-221"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=  8.430; volcLon(i)= 39.350; volcElev(i)= 1765;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-222"; volcName(i)="Sodore                         "; 
      i=i+1; volcLat(i)=  8.350; volcLon(i)= 39.180; volcElev(i)= 1984;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-23-"; volcName(i)="Gedamsa Caldera                "; 
      i=i+1; volcLat(i)=  8.270; volcLon(i)= 39.030; volcElev(i)= 2285;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-24-"; volcName(i)="Bora-Bericcio                  "; 
      i=i+1; volcLat(i)=  8.158; volcLon(i)= 39.130; volcElev(i)= 2349;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-25-"; volcName(i)="Tullu Moje                     "; 
      i=i+1; volcLat(i)=  8.070; volcLon(i)= 39.070; volcElev(i)= 1800;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-251"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=  7.950; volcLon(i)= 38.930; volcElev(i)= 1889;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-252"; volcName(i)="East Zway                      "; 
      i=i+1; volcLat(i)=  8.050; volcLon(i)= 38.350; volcElev(i)= 2281;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-26-"; volcName(i)="Butajiri-Silti Field           "; 
      i=i+1; volcLat(i)=  7.770; volcLon(i)= 38.780; volcElev(i)= 2335;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-27-"; volcName(i)="Alutu                          "; 
      i=i+1; volcLat(i)=  7.470; volcLon(i)= 38.580; volcElev(i)= 2075;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-28-"; volcName(i)="O'a Caldera                    "; 
      i=i+1; volcLat(i)=  7.180; volcLon(i)= 38.430; volcElev(i)= 2320;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-29-"; volcName(i)="Corbetti Caldera               "; 
      i=i+1; volcLat(i)=  7.070; volcLon(i)= 38.100; volcElev(i)= 1700;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-291"; volcName(i)="Bilate River Field             "; 
      i=i+1; volcLat(i)=  7.420; volcLon(i)= 35.430; volcElev(i)= 2728;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-292"; volcName(i)="Tepi                           "; 
      i=i+1; volcLat(i)=  6.780; volcLon(i)= 37.830; volcElev(i)= 1800;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-293"; volcName(i)="Hobicha Caldera                "; 
      i=i+1; volcLat(i)=  6.650; volcLon(i)= 38.120; volcElev(i)= 1650;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0201-30-"; volcName(i)="Chiracha                       "; 
      i=i+1; volcLat(i)=  5.930; volcLon(i)= 37.570; volcElev(i)= 1650;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-31-"; volcName(i)="Tosa Sucha                     "; 
      i=i+1; volcLat(i)=  5.650; volcLon(i)= 37.670; volcElev(i)= 1200;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-311"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=  5.100; volcLon(i)= 35.880; volcElev(i)=  912;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-32-"; volcName(i)="Korath Range                   "; 
      i=i+1; volcLat(i)=  4.080; volcLon(i)= 37.420; volcElev(i)= 1067;  volcLoc(i)="Ethiopia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0201-33-"; volcName(i)="Mega Basalt Field              "; 
      i=i+1; volcLat(i)=  4.070; volcLon(i)= 36.050; volcElev(i)=  520;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0202-001"; volcName(i)="North Island                   "; 
      i=i+1; volcLat(i)=  3.500; volcLon(i)= 36.042; volcElev(i)=  550;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0202-01="; volcName(i)="Central Island                 "; 
      i=i+1; volcLat(i)=  2.320; volcLon(i)= 37.970; volcElev(i)= 1707;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0202-021"; volcName(i)="Marsabit                       "; 
      i=i+1; volcLat(i)=  2.630; volcLon(i)= 36.600; volcElev(i)=  800;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0202-02="; volcName(i)="South Island                   "; 
      i=i+1; volcLat(i)=  2.320; volcLon(i)= 36.570; volcElev(i)= 1032;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="S1"; volcID(i)="0202-03="; volcName(i)="Barrier, The                   "; 
      i=i+1; volcLat(i)=  1.980; volcLon(i)= 36.430; volcElev(i)=  817;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0202-04-"; volcName(i)="Namarunu                       "; 
      i=i+1; volcLat(i)=  1.570; volcLon(i)= 37.900; volcElev(i)=  699;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0202-05-"; volcName(i)="Segererua Plateau              "; 
      i=i+1; volcLat(i)=  1.500; volcLon(i)= 36.330; volcElev(i)= 1328;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0202-051"; volcName(i)="Emuruangogolak                 "; 
      i=i+1; volcLat(i)=  1.150; volcLon(i)= 36.230; volcElev(i)= 1528;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0202-052"; volcName(i)="Silali                         "; 
      i=i+1; volcLat(i)=  0.920; volcLon(i)= 36.180; volcElev(i)= 1697;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0202-053"; volcName(i)="Paka                           "; 
      i=i+1; volcLat(i)=  0.770; volcLon(i)= 36.120; volcElev(i)= 1446;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0202-054"; volcName(i)="Korosi                         "; 
      i=i+1; volcLat(i)=  0.620; volcLon(i)= 36.075; volcElev(i)= 1130;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0202-055"; volcName(i)="Ol Kokwe                       "; 
      i=i+1; volcLat(i)=  0.230; volcLon(i)= 37.870; volcElev(i)=  750;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0202-056"; volcName(i)="Nyambeni Hills                 "; 
      i=i+1; volcLat(i)= -0.200; volcLon(i)= 36.070; volcElev(i)= 2278;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0202-06="; volcName(i)="Menengai                       "; 
      i=i+1; volcLat(i)= -0.520; volcLon(i)= 36.270; volcElev(i)= 2126;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0202-071"; volcName(i)="Elmenteita Badlands            "; 
      i=i+1; volcLat(i)= -0.380; volcLon(i)= 34.500; volcElev(i)= 1751;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0202-07="; volcName(i)="Homa Mountain                  "; 
      i=i+1; volcLat(i)= -0.650; volcLon(i)= 36.220; volcElev(i)= 2856;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0202-08="; volcName(i)="Eburru, Ol Doinyo              "; 
      i=i+1; volcLat(i)= -0.904; volcLon(i)= 36.292; volcElev(i)= 2434;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0202-09="; volcName(i)="Olkaria                        "; 
      i=i+1; volcLat(i)= -0.914; volcLon(i)= 36.446; volcElev(i)= 2776;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0202-10="; volcName(i)="Longonot                       "; 
      i=i+1; volcLat(i)= -1.175; volcLon(i)= 36.350; volcElev(i)= 2356;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0202-11="; volcName(i)="Suswa                          "; 
      i=i+1; volcLat(i)= -2.764; volcLon(i)= 35.914; volcElev(i)= 2962;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="S2"; volcID(i)="0202-12="; volcName(i)="Lengai, Ol Doinyo              "; 
      i=i+1; volcLat(i)= -2.680; volcLon(i)= 37.880; volcElev(i)= 2188;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0202-13="; volcName(i)="Chyulu Hills                   "; 
      i=i+1; volcLat(i)= -3.070; volcLon(i)= 37.350; volcElev(i)= 5895;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0202-15="; volcName(i)="Kilimanjaro                    "; 
      i=i+1; volcLat(i)= -4.870; volcLon(i)= 31.920; volcElev(i)=    0;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0202-161"; volcName(i)="Igwisi Hills                   "; 
      i=i+1; volcLat(i)= -8.630; volcLon(i)= 33.570; volcElev(i)=    0;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0202-162"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= -8.750; volcLon(i)= 33.800; volcElev(i)= 2179;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0202-163"; volcName(i)="SW Usangu Basin                "; 
      i=i+1; volcLat(i)= -8.970; volcLon(i)= 33.570; volcElev(i)= 2622;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0202-164"; volcName(i)="Ngozi                          "; 
      i=i+1; volcLat(i)= -8.930; volcLon(i)= 33.400; volcElev(i)= 1568;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0202-165"; volcName(i)="Izumbwe-Mpoli                  "; 
      i=i+1; volcLat(i)= -9.130; volcLon(i)= 33.670; volcElev(i)= 2961;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0202-166"; volcName(i)="Rungwe                         "; 
      i=i+1; volcLat(i)= -3.250; volcLon(i)= 36.750; volcElev(i)= 4565;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0202-16="; volcName(i)="Meru                           "; 
      i=i+1; volcLat(i)= -9.230; volcLon(i)= 33.780; volcElev(i)= 2175;  volcLoc(i)="Africa-E                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0202-17="; volcName(i)="Kieyo                          "; 
      i=i+1; volcLat(i)=  0.700; volcLon(i)= 30.250; volcElev(i)= 1615;  volcLoc(i)="Africa-C                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0203-001"; volcName(i)="Fort Portal                    "; 
      i=i+1; volcLat(i)=  0.450; volcLon(i)= 30.250; volcElev(i)= 1430;  volcLoc(i)="Africa-C                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0203-002"; volcName(i)="Kyatwa                         "; 
      i=i+1; volcLat(i)= -0.080; volcLon(i)= 29.920; volcElev(i)= 1067;  volcLoc(i)="Africa-C                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0203-003"; volcName(i)="Katwe-Kikorongo                "; 
      i=i+1; volcLat(i)= -0.200; volcLon(i)= 30.080; volcElev(i)= 1554;  volcLoc(i)="Africa-C                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0203-004"; volcName(i)="Bunyaruguru                    "; 
      i=i+1; volcLat(i)= -0.471; volcLon(i)= 30.191; volcElev(i)= 1707;  volcLoc(i)="Africa-C                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0203-005"; volcName(i)="Katunga                        "; 
      i=i+1; volcLat(i)= -0.930; volcLon(i)= 29.330; volcElev(i)=  950;  volcLoc(i)="Africa-C                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0203-01="; volcName(i)="May-ya-moto                    "; 
      i=i+1; volcLat(i)= -1.408; volcLon(i)= 29.200; volcElev(i)= 3058;  volcLoc(i)="Africa-C                      "; 
             volcESP_Code(i)="M1"; volcID(i)="0203-02="; volcName(i)="Nyamuragira                    "; 
      i=i+1; volcLat(i)= -1.520; volcLon(i)= 29.250; volcElev(i)= 3470;  volcLoc(i)="Africa-C                      "; 
             volcESP_Code(i)="M1"; volcID(i)="0203-03="; volcName(i)="Nyiragongo                     "; 
      i=i+1; volcLat(i)= -1.500; volcLon(i)= 29.450; volcElev(i)= 4507;  volcLoc(i)="Africa-C                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0203-04-"; volcName(i)="Karisimbi                      "; 
      i=i+1; volcLat(i)= -1.470; volcLon(i)= 29.492; volcElev(i)= 3711;  volcLoc(i)="Africa-C                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0203-05-"; volcName(i)="Visoke                         "; 
      i=i+1; volcLat(i)= -1.380; volcLon(i)= 29.670; volcElev(i)= 4127;  volcLoc(i)="Africa-C                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0203-06-"; volcName(i)="Muhavura                       "; 
      i=i+1; volcLat(i)= -1.230; volcLon(i)= 29.720; volcElev(i)= 2440;  volcLoc(i)="Africa-C                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0203-07-"; volcName(i)="Bufumbira                      "; 
      i=i+1; volcLat(i)= -2.320; volcLon(i)= 28.750; volcElev(i)= 1460;  volcLoc(i)="Africa-C                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0203-08-"; volcName(i)="Tshibinda                      "; 
      i=i+1; volcLat(i)=  0.200; volcLon(i)=  6.580; volcElev(i)= 2024;  volcLoc(i)="Africa-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0204-001"; volcName(i)="Sao Tome                       "; 
      i=i+1; volcLat(i)=  3.350; volcLon(i)=  8.520; volcElev(i)= 2260;  volcLoc(i)="Africa-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0204-002"; volcName(i)="San Carlos                     "; 
      i=i+1; volcLat(i)=  3.350; volcLon(i)=  8.630; volcElev(i)= 2009;  volcLoc(i)="Africa-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0204-003"; volcName(i)="San Joaquin                    "; 
      i=i+1; volcLat(i)=  3.580; volcLon(i)=  8.750; volcElev(i)= 3007;  volcLoc(i)="Africa-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0204-004"; volcName(i)="Santa Isabel                   "; 
      i=i+1; volcLat(i)=  4.750; volcLon(i)=  9.670; volcElev(i)=  500;  volcLoc(i)="Africa-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0204-011"; volcName(i)="Tombel Graben                  "; 
      i=i+1; volcLat(i)=  4.203; volcLon(i)=  9.170; volcElev(i)= 4095;  volcLoc(i)="Africa-W                      "; 
             volcESP_Code(i)="M1"; volcID(i)="0204-01="; volcName(i)="Cameroon                       "; 
      i=i+1; volcLat(i)=  5.030; volcLon(i)=  9.830; volcElev(i)= 2411;  volcLoc(i)="Africa-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0204-02-"; volcName(i)="Manengouba                     "; 
      i=i+1; volcLat(i)=  6.250; volcLon(i)= 10.500; volcElev(i)= 3011;  volcLoc(i)="Africa-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0204-03-"; volcName(i)="Oku Volc Field                 "; 
      i=i+1; volcLat(i)=  7.250; volcLon(i)= 13.670; volcElev(i)=    0;  volcLoc(i)="Africa-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0204-04-"; volcName(i)="Ngaoundere Plateau             "; 
      i=i+1; volcLat(i)= 10.750; volcLon(i)= 12.000; volcElev(i)=    0;  volcLoc(i)="Africa-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0204-05-"; volcName(i)="Biu Plateau                    "; 
      i=i+1; volcLat(i)= 17.680; volcLon(i)=  8.500; volcElev(i)= 1780;  volcLoc(i)="Africa-N                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0205-001"; volcName(i)="Todra Volc Field               "; 
      i=i+1; volcLat(i)= 19.830; volcLon(i)=  2.830; volcElev(i)=    0;  volcLoc(i)="Africa-N                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0205-002"; volcName(i)="Tin Zaouatene Volc Field       "; 
      i=i+1; volcLat(i)= 23.000; volcLon(i)= 10.830; volcElev(i)=    0;  volcLoc(i)="Africa-N                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0205-003"; volcName(i)="In Ezzane Volc Field           "; 
      i=i+1; volcLat(i)= 22.670; volcLon(i)=  5.000; volcElev(i)= 1467;  volcLoc(i)="Africa-N                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0205-004"; volcName(i)="Tahalra Volc Field             "; 
      i=i+1; volcLat(i)= 23.330; volcLon(i)=  5.830; volcElev(i)= 2918;  volcLoc(i)="Africa-N                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0205-005"; volcName(i)="Atakor Volc Field              "; 
      i=i+1; volcLat(i)= 23.920; volcLon(i)=  5.830; volcElev(i)= 1672;  volcLoc(i)="Africa-N                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0205-006"; volcName(i)="Manzaz Volc Field              "; 
      i=i+1; volcLat(i)= 27.250; volcLon(i)= 17.500; volcElev(i)= 1200;  volcLoc(i)="Africa-N                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0205-007"; volcName(i)="Haruj                          "; 
      i=i+1; volcLat(i)= 25.050; volcLon(i)= 17.550; volcElev(i)=  547;  volcLoc(i)="Africa-N                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0205-008"; volcName(i)="Wau-en-Namus                   "; 
      i=i+1; volcLat(i)= 21.330; volcLon(i)= 16.330; volcElev(i)= 2000;  volcLoc(i)="Africa-N                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0205-009"; volcName(i)="Tôh, Tarso                     "; 
      i=i+1; volcLat(i)= 21.030; volcLon(i)= 16.450; volcElev(i)= 3265;  volcLoc(i)="Africa-N                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0205-01="; volcName(i)="Toussidé, Tarso                "; 
      i=i+1; volcLat(i)= 19.800; volcLon(i)= 18.530; volcElev(i)= 3415;  volcLoc(i)="Africa-N                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0205-021"; volcName(i)="Koussi, Emi                    "; 
      i=i+1; volcLat(i)= 20.920; volcLon(i)= 17.280; volcElev(i)= 3100;  volcLoc(i)="Africa-N                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0205-02="; volcName(i)="Voon, Tarso                    "; 
      i=i+1; volcLat(i)= 12.950; volcLon(i)= 24.270; volcElev(i)= 3042;  volcLoc(i)="Africa-N                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0205-03-"; volcName(i)="Marra, Jebel                   "; 
      i=i+1; volcLat(i)= 14.570; volcLon(i)= 25.850; volcElev(i)=    0;  volcLoc(i)="Africa-N                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0205-04-"; volcName(i)="Kutum Volc Field               "; 
      i=i+1; volcLat(i)= 15.320; volcLon(i)= 26.470; volcElev(i)= 2000;  volcLoc(i)="Africa-N                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0205-05-"; volcName(i)="Meidob Volc Field              "; 
      i=i+1; volcLat(i)= 18.330; volcLon(i)= 32.750; volcElev(i)=  670;  volcLoc(i)="Africa-N                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0205-06-"; volcName(i)="Bayuda Volc Field              "; 
      i=i+1; volcLat(i)= 18.170; volcLon(i)= 33.830; volcElev(i)=    0;  volcLoc(i)="Africa-N                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0205-07-"; volcName(i)="Umm Arafieb, Jebel             "; 
      i=i+1; volcLat(i)= 36.530; volcLon(i)= 40.850; volcElev(i)=  534;  volcLoc(i)="Syria                         "; 
             volcESP_Code(i)="M0"; volcID(i)="0300-01-"; volcName(i)="Sharat Kovakab                 "; 
      i=i+1; volcLat(i)= 36.670; volcLon(i)= 37.000; volcElev(i)=    0;  volcLoc(i)="Syria                         "; 
             volcESP_Code(i)="M0"; volcID(i)="0300-02-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 33.100; volcLon(i)= 35.970; volcElev(i)= 1197;  volcLoc(i)="Syria                         "; 
             volcESP_Code(i)="M0"; volcID(i)="0300-03-"; volcName(i)="Golan Heights                  "; 
      i=i+1; volcLat(i)= 33.000; volcLon(i)= 36.430; volcElev(i)=  945;  volcLoc(i)="Syria                         "; 
             volcESP_Code(i)="M0"; volcID(i)="0300-04-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 33.250; volcLon(i)= 37.070; volcElev(i)=  979;  volcLoc(i)="Syria                         "; 
             volcESP_Code(i)="M0"; volcID(i)="0300-05-"; volcName(i)="Es Safa                        "; 
      i=i+1; volcLat(i)= 32.658; volcLon(i)= 36.425; volcElev(i)= 1803;  volcLoc(i)="Syria                         "; 
             volcESP_Code(i)="M0"; volcID(i)="0300-06-"; volcName(i)="Druze, Jabal ad                "; 
      i=i+1; volcLat(i)= 31.080; volcLon(i)= 38.420; volcElev(i)= 1100;  volcLoc(i)="Arabia-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0301-001"; volcName(i)="Harrah, Al                     "; 
      i=i+1; volcLat(i)= 27.800; volcLon(i)= 36.170; volcElev(i)= 1950;  volcLoc(i)="Arabia-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0301-01="; volcName(i)="Rahah, Harrat ar               "; 
      i=i+1; volcLat(i)= 27.080; volcLon(i)= 37.250; volcElev(i)= 1920;  volcLoc(i)="Arabia-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0301-02="; volcName(i)="'Uwayrid, Harrat               "; 
      i=i+1; volcLat(i)= 25.170; volcLon(i)= 37.750; volcElev(i)= 1370;  volcLoc(i)="Arabia-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0301-04-"; volcName(i)="Lunayyir, Harrat               "; 
      i=i+1; volcLat(i)= 26.580; volcLon(i)= 40.200; volcElev(i)= 1625;  volcLoc(i)="Arabia-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0301-05="; volcName(i)="Ithnayn, Harrat                "; 
      i=i+1; volcLat(i)= 25.000; volcLon(i)= 39.920; volcElev(i)= 2093;  volcLoc(i)="Arabia-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0301-06="; volcName(i)="Khaybar, Harrat                "; 
      i=i+1; volcLat(i)= 22.800; volcLon(i)= 41.380; volcElev(i)= 1475;  volcLoc(i)="Arabia-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0301-071"; volcName(i)="Kishb, Harrat                  "; 
      i=i+1; volcLat(i)= 18.370; volcLon(i)= 41.630; volcElev(i)=  381;  volcLoc(i)="Arabia-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0301-072"; volcName(i)="Birk, Harrat al                "; 
      i=i+1; volcLat(i)= 23.080; volcLon(i)= 39.780; volcElev(i)= 1744;  volcLoc(i)="Arabia-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0301-07="; volcName(i)="Rahat, Harrat                  "; 
      i=i+1; volcLat(i)= 17.050; volcLon(i)= 42.830; volcElev(i)=  305;  volcLoc(i)="Arabia-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0301-08-"; volcName(i)="Yar, Jabal                     "; 
      i=i+1; volcLat(i)= 15.630; volcLon(i)= 44.080; volcElev(i)= 3100;  volcLoc(i)="Arabia-S                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0301-09-"; volcName(i)="Arhab, Harra of                "; 
      i=i+1; volcLat(i)= 15.245; volcLon(i)= 44.236; volcElev(i)= 2506;  volcLoc(i)="Arabia-S                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0301-10-"; volcName(i)="Marha, Jabal el-               "; 
      i=i+1; volcLat(i)= 15.430; volcLon(i)= 44.780; volcElev(i)= 1550;  volcLoc(i)="Arabia-S                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0301-11-"; volcName(i)="Haylan, Jabal                  "; 
      i=i+1; volcLat(i)= 14.570; volcLon(i)= 44.670; volcElev(i)= 3500;  volcLoc(i)="Arabia-S                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0301-12-"; volcName(i)="Dhamar, Harras of              "; 
      i=i+1; volcLat(i)= 12.250; volcLon(i)= 45.000; volcElev(i)=    0;  volcLoc(i)="Arabia-S                      "; 
             volcESP_Code(i)="U0"; volcID(i)="0301-15-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 13.580; volcLon(i)= 46.120; volcElev(i)= 1737;  volcLoc(i)="Arabia-S                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0301-16-"; volcName(i)="Sawâd, Harra es-               "; 
      i=i+1; volcLat(i)= 14.050; volcLon(i)= 48.330; volcElev(i)=  233;  volcLoc(i)="Arabia-S                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0301-17-"; volcName(i)="Bal Haf, Harra of              "; 
      i=i+1; volcLat(i)= 15.550; volcLon(i)= 50.630; volcElev(i)=    0;  volcLoc(i)="Arabia-S                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0301-18-"; volcName(i)="Bir Borhut                     "; 
      i=i+1; volcLat(i)= 39.330; volcLon(i)= 45.170; volcElev(i)=    0;  volcLoc(i)="Iran                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0302-00-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 37.750; volcLon(i)= 46.430; volcElev(i)= 3707;  volcLoc(i)="Iran                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0302-001"; volcName(i)="Sahand                         "; 
      i=i+1; volcLat(i)= 38.250; volcLon(i)= 47.920; volcElev(i)= 4811;  volcLoc(i)="Iran                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0302-002"; volcName(i)="Sabalan                        "; 
      i=i+1; volcLat(i)= 35.951; volcLon(i)= 52.109; volcElev(i)= 5670;  volcLoc(i)="Iran                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0302-01-"; volcName(i)="Damavand                       "; 
      i=i+1; volcLat(i)= 29.400; volcLon(i)= 57.570; volcElev(i)=    0;  volcLoc(i)="Iran                          "; 
             volcESP_Code(i)="M0"; volcID(i)="0302-02-"; volcName(i)="Qal'eh Hasan Ali               "; 
      i=i+1; volcLat(i)= 28.070; volcLon(i)= 60.000; volcElev(i)= 3490;  volcLoc(i)="Iran                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0302-03-"; volcName(i)="Bazman                         "; 
      i=i+1; volcLat(i)= 28.170; volcLon(i)= 60.670; volcElev(i)=    0;  volcLoc(i)="Iran                          "; 
             volcESP_Code(i)="M0"; volcID(i)="0302-04-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 28.600; volcLon(i)= 61.130; volcElev(i)= 3940;  volcLoc(i)="Iran                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0302-05-"; volcName(i)="Taftan                         "; 
      i=i+1; volcLat(i)= 33.950; volcLon(i)= 67.920; volcElev(i)= 3800;  volcLoc(i)="Afghanistan                   "; 
             volcESP_Code(i)="S0"; volcID(i)="0302-06-"; volcName(i)="Dacht-i-Navar Group            "; 
      i=i+1; volcLat(i)= 34.250; volcLon(i)= 67.970; volcElev(i)= 3190;  volcLoc(i)="Afghanistan                   "; 
             volcESP_Code(i)="S0"; volcID(i)="0302-07-"; volcName(i)="Vakak Group                    "; 
      i=i+1; volcLat(i)=-11.470; volcLon(i)= 43.330; volcElev(i)= 1087;  volcLoc(i)="Indian O.-W                   "; 
             volcESP_Code(i)="M0"; volcID(i)="0303-001"; volcName(i)="Grille, La                     "; 
      i=i+1; volcLat(i)=-12.600; volcLon(i)= 49.150; volcElev(i)= 1475;  volcLoc(i)="Madagascar                    "; 
             volcESP_Code(i)="M0"; volcID(i)="0303-011"; volcName(i)="Ambre-Bobaomby                 "; 
      i=i+1; volcLat(i)=-13.320; volcLon(i)= 48.480; volcElev(i)=  214;  volcLoc(i)="Madagascar                    "; 
             volcESP_Code(i)="S0"; volcID(i)="0303-012"; volcName(i)="Nosy-Be                        "; 
      i=i+1; volcLat(i)=-14.300; volcLon(i)= 48.670; volcElev(i)= 2878;  volcLoc(i)="Madagascar                    "; 
             volcESP_Code(i)="M0"; volcID(i)="0303-013"; volcName(i)="Ankaizina Field                "; 
      i=i+1; volcLat(i)=-19.000; volcLon(i)= 46.770; volcElev(i)= 1800;  volcLoc(i)="Madagascar                    "; 
             volcESP_Code(i)="S0"; volcID(i)="0303-014"; volcName(i)="Itasy Volc Field               "; 
      i=i+1; volcLat(i)=-19.400; volcLon(i)= 47.200; volcElev(i)= 2644;  volcLoc(i)="Madagascar                    "; 
             volcESP_Code(i)="S0"; volcID(i)="0303-015"; volcName(i)="Ankaratra Field                "; 
      i=i+1; volcLat(i)=-11.750; volcLon(i)= 43.380; volcElev(i)= 2361;  volcLoc(i)="Indian O.-W                   "; 
             volcESP_Code(i)="M1"; volcID(i)="0303-01="; volcName(i)="Karthala                       "; 
      i=i+1; volcLat(i)=-21.231; volcLon(i)= 55.713; volcElev(i)= 2632;  volcLoc(i)="Indian O.-W                   "; 
             volcESP_Code(i)="M1"; volcID(i)="0303-02="; volcName(i)="Fournaise, Piton de la         "; 
      i=i+1; volcLat(i)=-37.721; volcLon(i)= 77.825; volcElev(i)= -650;  volcLoc(i)="Indian O.-S                   "; 
             volcESP_Code(i)="U0"; volcID(i)="0304-00-"; volcName(i)="Boomerang Seamount             "; 
      i=i+1; volcLat(i)=-37.830; volcLon(i)= 77.520; volcElev(i)=  881;  volcLoc(i)="Indian O.-S                   "; 
             volcESP_Code(i)="M0"; volcID(i)="0304-001"; volcName(i)="Amsterdam Island               "; 
      i=i+1; volcLat(i)=-38.720; volcLon(i)= 77.530; volcElev(i)=  268;  volcLoc(i)="Indian O.-S                   "; 
             volcESP_Code(i)="M0"; volcID(i)="0304-002"; volcName(i)="St. Paul                       "; 
      i=i+1; volcLat(i)=-53.030; volcLon(i)= 72.600; volcElev(i)=  230;  volcLoc(i)="Indian O.-S                   "; 
             volcESP_Code(i)="S1"; volcID(i)="0304-011"; volcName(i)="McDonald Islands               "; 
      i=i+1; volcLat(i)=-53.106; volcLon(i)= 73.513; volcElev(i)= 2745;  volcLoc(i)="Indian O.-S                   "; 
             volcESP_Code(i)="M1"; volcID(i)="0304-01="; volcName(i)="Heard                          "; 
      i=i+1; volcLat(i)=-49.580; volcLon(i)= 69.500; volcElev(i)= 1840;  volcLoc(i)="Indian O.-S                   "; 
             volcESP_Code(i)="S0"; volcID(i)="0304-02="; volcName(i)="Kerguelen Islands              "; 
      i=i+1; volcLat(i)=-46.430; volcLon(i)= 52.200; volcElev(i)= 1090;  volcLoc(i)="Indian O.-S                   "; 
             volcESP_Code(i)="M0"; volcID(i)="0304-03-"; volcName(i)="Est, Ile de l'                 "; 
      i=i+1; volcLat(i)=-46.420; volcLon(i)= 51.750; volcElev(i)=  934;  volcLoc(i)="Indian O.-S                   "; 
             volcESP_Code(i)="M0"; volcID(i)="0304-04-"; volcName(i)="Possession, Ile de la          "; 
      i=i+1; volcLat(i)=-46.100; volcLon(i)= 50.230; volcElev(i)=  775;  volcLoc(i)="Indian O.-S                   "; 
             volcESP_Code(i)="M0"; volcID(i)="0304-05-"; volcName(i)="Cochons, Ile aux               "; 
      i=i+1; volcLat(i)=-46.630; volcLon(i)= 37.950; volcElev(i)=  672;  volcLoc(i)="Indian O.                     "; 
             volcESP_Code(i)="M0"; volcID(i)="0304-06-"; volcName(i)="Prince Edward Island           "; 
      i=i+1; volcLat(i)=-46.900; volcLon(i)= 37.750; volcElev(i)= 1230;  volcLoc(i)="Indian O.-S                   "; 
             volcESP_Code(i)="M1"; volcID(i)="0304-07-"; volcName(i)="Marion Island                  "; 
      i=i+1; volcLat(i)= 11.750; volcLon(i)= 80.750; volcElev(i)=    0;  volcLoc(i)="Indian O.-E                   "; 
             volcESP_Code(i)="U0"; volcID(i)="0305-01="; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=-35.750; volcLon(i)=174.270; volcElev(i)=  397;  volcLoc(i)="New Zealand                   "; 
             volcESP_Code(i)="M0"; volcID(i)="0401-011"; volcName(i)="Whangarei                      "; 
      i=i+1; volcLat(i)=-35.300; volcLon(i)=173.900; volcElev(i)=  388;  volcLoc(i)="New Zealand                   "; 
             volcESP_Code(i)="M0"; volcID(i)="0401-01="; volcName(i)="Kaikohe-Bay of Islands         "; 
      i=i+1; volcLat(i)=-37.280; volcLon(i)=176.250; volcElev(i)=  355;  volcLoc(i)="New Zealand                   "; 
             volcESP_Code(i)="S0"; volcID(i)="0401-021"; volcName(i)="Mayor Island                   "; 
      i=i+1; volcLat(i)=-36.900; volcLon(i)=174.870; volcElev(i)=  260;  volcLoc(i)="New Zealand                   "; 
             volcESP_Code(i)="M0"; volcID(i)="0401-02="; volcName(i)="Auckland Field                 "; 
      i=i+1; volcLat(i)=-39.300; volcLon(i)=174.070; volcElev(i)= 2518;  volcLoc(i)="New Zealand                   "; 
             volcESP_Code(i)="S0"; volcID(i)="0401-03="; volcName(i)="Taranaki                       "; 
      i=i+1; volcLat(i)=-37.520; volcLon(i)=177.180; volcElev(i)=  321;  volcLoc(i)="New Zealand                   "; 
             volcESP_Code(i)="S1"; volcID(i)="0401-04="; volcName(i)="White Island                   "; 
      i=i+1; volcLat(i)=-38.120; volcLon(i)=176.500; volcElev(i)= 1111;  volcLoc(i)="New Zealand                   "; 
             volcESP_Code(i)="S0"; volcID(i)="0401-05="; volcName(i)="Okataina                       "; 
      i=i+1; volcLat(i)=-38.420; volcLon(i)=176.330; volcElev(i)=  592;  volcLoc(i)="New Zealand                   "; 
             volcESP_Code(i)="S0"; volcID(i)="0401-06-"; volcName(i)="Reporoa                        "; 
      i=i+1; volcLat(i)=-38.420; volcLon(i)=176.080; volcElev(i)= 1156;  volcLoc(i)="New Zealand                   "; 
             volcESP_Code(i)="S0"; volcID(i)="0401-061"; volcName(i)="Maroa                          "; 
      i=i+1; volcLat(i)=-38.820; volcLon(i)=176.000; volcElev(i)=  760;  volcLoc(i)="New Zealand                   "; 
             volcESP_Code(i)="S0"; volcID(i)="0401-07="; volcName(i)="Taupo                          "; 
      i=i+1; volcLat(i)=-39.130; volcLon(i)=175.642; volcElev(i)= 1978;  volcLoc(i)="New Zealand                   "; 
             volcESP_Code(i)="S1"; volcID(i)="0401-08="; volcName(i)="Tongariro                      "; 
      i=i+1; volcLat(i)=-36.446; volcLon(i)=177.839; volcElev(i)= -860;  volcLoc(i)="New Zealand                   "; 
             volcESP_Code(i)="U0"; volcID(i)="0401-101"; volcName(i)="Clark                          "; 
      i=i+1; volcLat(i)=-36.321; volcLon(i)=178.028; volcElev(i)=  600;  volcLoc(i)="New Zealand                   "; 
             volcESP_Code(i)="U0"; volcID(i)="0401-102"; volcName(i)="Tangaroa                       "; 
      i=i+1; volcLat(i)=-39.280; volcLon(i)=175.570; volcElev(i)= 2797;  volcLoc(i)="New Zealand                   "; 
             volcESP_Code(i)="S1"; volcID(i)="0401-10="; volcName(i)="Ruapehu                        "; 
      i=i+1; volcLat(i)=-36.142; volcLon(i)=178.196; volcElev(i)=  400;  volcLoc(i)="New Zealand                   "; 
             volcESP_Code(i)="U0"; volcID(i)="0401-11-"; volcName(i)="Rumble V                       "; 
      i=i+1; volcLat(i)=-36.130; volcLon(i)=178.050; volcElev(i)=  500;  volcLoc(i)="New Zealand                   "; 
             volcESP_Code(i)="U0"; volcID(i)="0401-12-"; volcName(i)="Rumble IV                      "; 
      i=i+1; volcLat(i)=-35.745; volcLon(i)=178.478; volcElev(i)= -220;  volcLoc(i)="New Zealand                   "; 
             volcESP_Code(i)="U0"; volcID(i)="0401-13-"; volcName(i)="Rumble III                     "; 
      i=i+1; volcLat(i)=-35.353; volcLon(i)=178.527; volcElev(i)= 1200;  volcLoc(i)="New Zealand                   "; 
             volcESP_Code(i)="U0"; volcID(i)="0401-131"; volcName(i)="Rumble II West                 "; 
      i=i+1; volcLat(i)=-35.004; volcLon(i)=178.973; volcElev(i)=  980;  volcLoc(i)="New Zealand                   "; 
             volcESP_Code(i)="U0"; volcID(i)="0401-14-"; volcName(i)="Healy                          "; 
      i=i+1; volcLat(i)=-34.875; volcLon(i)=179.075; volcElev(i)=-1350;  volcLoc(i)="New Zealand                   "; 
             volcESP_Code(i)="U0"; volcID(i)="0401-15-"; volcName(i)="Brothers                       "; 
      i=i+1; volcLat(i)=-31.850; volcLon(i)=180.820; volcElev(i)= -900;  volcLoc(i)="Kermadec Is                   "; 
             volcESP_Code(i)="U0"; volcID(i)="0402-001"; volcName(i)="Volcano W                      "; 
      i=i+1; volcLat(i)=-30.542; volcLon(i)=181.439; volcElev(i)=  137;  volcLoc(i)="Kermadec Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="0402-01="; volcName(i)="Curtis Island                  "; 
      i=i+1; volcLat(i)=-30.200; volcLon(i)=181.530; volcElev(i)=  238;  volcLoc(i)="Kermadec Is                   "; 
             volcESP_Code(i)="M0"; volcID(i)="0402-021"; volcName(i)="Macauley Island                "; 
      i=i+1; volcLat(i)=-30.036; volcLon(i)=181.288; volcElev(i)=  -65;  volcLoc(i)="Kermadec Is                   "; 
             volcESP_Code(i)="U0"; volcID(i)="0402-022"; volcName(i)="Giggenbach                     "; 
      i=i+1; volcLat(i)=-29.270; volcLon(i)=182.080; volcElev(i)=  516;  volcLoc(i)="Kermadec Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="0402-03="; volcName(i)="Raoul Island                   "; 
      i=i+1; volcLat(i)=-25.887; volcLon(i)=182.812; volcElev(i)= -132;  volcLoc(i)="Kermadec Is                   "; 
             volcESP_Code(i)="U0"; volcID(i)="0402-05-"; volcName(i)="Monowai Seamount               "; 
      i=i+1; volcLat(i)=-24.800; volcLon(i)=182.980; volcElev(i)= -385;  volcLoc(i)="Tonga-SW Pacific              "; 
             volcESP_Code(i)="U0"; volcID(i)="0403-001"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=-21.150; volcLon(i)=184.250; volcElev(i)=  -65;  volcLoc(i)="Tonga-SW Pacific              "; 
             volcESP_Code(i)="U0"; volcID(i)="0403-011"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=-21.380; volcLon(i)=184.350; volcElev(i)= -500;  volcLoc(i)="Tonga-SW Pacific              "; 
             volcESP_Code(i)="U0"; volcID(i)="0403-01="; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=-20.850; volcLon(i)=184.470; volcElev(i)=  -13;  volcLoc(i)="Tonga-SW Pacific              "; 
             volcESP_Code(i)="S0"; volcID(i)="0403-03="; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=-20.570; volcLon(i)=184.620; volcElev(i)=  149;  volcLoc(i)="Tonga-SW Pacific              "; 
             volcESP_Code(i)="S0"; volcID(i)="0403-04="; volcName(i)="Hunga Tonga-Hunga Ha'apai      "; 
      i=i+1; volcLat(i)=-20.320; volcLon(i)=184.580; volcElev(i)=  -17;  volcLoc(i)="Tonga-SW Pacific              "; 
             volcESP_Code(i)="S0"; volcID(i)="0403-05="; volcName(i)="Falcon Island                  "; 
      i=i+1; volcLat(i)=-19.670; volcLon(i)=184.970; volcElev(i)= 1030;  volcLoc(i)="Tonga-SW Pacific              "; 
             volcESP_Code(i)="S0"; volcID(i)="0403-061"; volcName(i)="Kao                            "; 
      i=i+1; volcLat(i)=-19.750; volcLon(i)=184.930; volcElev(i)=  515;  volcLoc(i)="Tonga-SW Pacific              "; 
             volcESP_Code(i)="S0"; volcID(i)="0403-06="; volcName(i)="Tofua                          "; 
      i=i+1; volcLat(i)=-19.180; volcLon(i)=185.130; volcElev(i)=   43;  volcLoc(i)="Tonga-SW Pacific              "; 
             volcESP_Code(i)="S0"; volcID(i)="0403-07="; volcName(i)="Metis Shoal                    "; 
      i=i+1; volcLat(i)=-18.992; volcLon(i)=185.225; volcElev(i)=   -2;  volcLoc(i)="Tonga-SW Pacific              "; 
             volcESP_Code(i)="S0"; volcID(i)="0403-08="; volcName(i)="Home Reef                      "; 
      i=i+1; volcLat(i)=-18.325; volcLon(i)=185.635; volcElev(i)=  -40;  volcLoc(i)="Tonga-SW Pacific              "; 
             volcESP_Code(i)="S0"; volcID(i)="0403-091"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=-18.806; volcLon(i)=185.350; volcElev(i)=  540;  volcLoc(i)="Tonga-SW Pacific              "; 
             volcESP_Code(i)="S0"; volcID(i)="0403-09="; volcName(i)="Late                           "; 
      i=i+1; volcLat(i)=-15.850; volcLon(i)=186.280; volcElev(i)=  560;  volcLoc(i)="Tonga-SW Pacific              "; 
             volcESP_Code(i)="S0"; volcID(i)="0403-101"; volcName(i)="Tafahi                         "; 
      i=i+1; volcLat(i)=-15.620; volcLon(i)=186.330; volcElev(i)=  -33;  volcLoc(i)="Tonga-SW Pacific              "; 
             volcESP_Code(i)="S0"; volcID(i)="0403-102"; volcName(i)="Curacoa                        "; 
      i=i+1; volcLat(i)=-18.020; volcLon(i)=185.675; volcElev(i)=  180;  volcLoc(i)="Tonga-SW Pacific              "; 
             volcESP_Code(i)="S0"; volcID(i)="0403-10="; volcName(i)="Fonualei                       "; 
      i=i+1; volcLat(i)=-15.600; volcLon(i)=184.370; volcElev(i)=  260;  volcLoc(i)="Tonga-SW Pacific              "; 
             volcESP_Code(i)="M0"; volcID(i)="0403-11="; volcName(i)="Niuafo'ou                      "; 
      i=i+1; volcLat(i)=-14.215; volcLon(i)=190.942; volcElev(i)= -592;  volcLoc(i)="Samoa-SW Pacific              "; 
             volcESP_Code(i)="U0"; volcID(i)="0404-00-"; volcName(i)="Vailulu'u                      "; 
      i=i+1; volcLat(i)=-14.230; volcLon(i)=190.546; volcElev(i)=  931;  volcLoc(i)="Samoa-SW Pacific              "; 
             volcESP_Code(i)="M0"; volcID(i)="0404-001"; volcName(i)="Ta'u                           "; 
      i=i+1; volcLat(i)=-14.175; volcLon(i)=190.382; volcElev(i)=  639;  volcLoc(i)="Samoa-SW Pacific              "; 
             volcESP_Code(i)="M0"; volcID(i)="0404-01="; volcName(i)="Ofu-Olosega                    "; 
      i=i+1; volcLat(i)=-14.295; volcLon(i)=189.300; volcElev(i)=  653;  volcLoc(i)="Samoa-SW Pacific              "; 
             volcESP_Code(i)="M0"; volcID(i)="0404-02-"; volcName(i)="Tutuila                        "; 
      i=i+1; volcLat(i)=-13.935; volcLon(i)=188.280; volcElev(i)= 1100;  volcLoc(i)="Samoa-SW Pacific              "; 
             volcESP_Code(i)="M0"; volcID(i)="0404-03-"; volcName(i)="Upolu                          "; 
      i=i+1; volcLat(i)=-13.612; volcLon(i)=187.475; volcElev(i)= 1858;  volcLoc(i)="Samoa-SW Pacific              "; 
             volcESP_Code(i)="M0"; volcID(i)="0404-04="; volcName(i)="Savai'i                        "; 
      i=i+1; volcLat(i)=-13.300; volcLon(i)=183.830; volcElev(i)=  143;  volcLoc(i)="SW Pacific                    "; 
             volcESP_Code(i)="M0"; volcID(i)="0404-05-"; volcName(i)="Wallis Islands                 "; 
      i=i+1; volcLat(i)=-16.820; volcLon(i)=180.030; volcElev(i)= 1241;  volcLoc(i)="Fiji Is-SW Pacific            "; 
             volcESP_Code(i)="M0"; volcID(i)="0405-01-"; volcName(i)="Taveuni                        "; 
      i=i+1; volcLat(i)=-17.320; volcLon(i)=179.400; volcElev(i)=  522;  volcLoc(i)="Fiji Is-SW Pacific            "; 
             volcESP_Code(i)="M0"; volcID(i)="0405-02-"; volcName(i)="Koro                           "; 
      i=i+1; volcLat(i)=-19.120; volcLon(i)=177.980; volcElev(i)=  805;  volcLoc(i)="Fiji Is-SW Pacific            "; 
             volcESP_Code(i)="S0"; volcID(i)="0405-03-"; volcName(i)="Nabukelevu                     "; 
      i=i+1; volcLat(i)= -2.380; volcLon(i)=147.350; volcElev(i)=  270;  volcLoc(i)="Admiralty Is-SW Pacific       "; 
             volcESP_Code(i)="S0"; volcID(i)="0500-01="; volcName(i)="St. Andrew Strait              "; 
      i=i+1; volcLat(i)= -2.570; volcLon(i)=147.280; volcElev(i)=  254;  volcLoc(i)="Admiralty Is-SW Pacific       "; 
             volcESP_Code(i)="M0"; volcID(i)="0500-02-"; volcName(i)="Baluan                         "; 
      i=i+1; volcLat(i)= -3.030; volcLon(i)=147.780; volcElev(i)=-1300;  volcLoc(i)="Admiralty Is-SW Pacific       "; 
             volcESP_Code(i)="U0"; volcID(i)="0500-03-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= -3.507; volcLon(i)=144.605; volcElev(i)=  402;  volcLoc(i)="New Guinea-NE of              "; 
             volcESP_Code(i)="S0"; volcID(i)="0501-001"; volcName(i)="Blup Blup                      "; 
      i=i+1; volcLat(i)= -3.630; volcLon(i)=144.631; volcElev(i)=  365;  volcLoc(i)="New Guinea-NE of              "; 
             volcESP_Code(i)="S0"; volcID(i)="0501-002"; volcName(i)="Kadovar                        "; 
      i=i+1; volcLat(i)= -3.994; volcLon(i)=144.963; volcElev(i)=  240;  volcLoc(i)="New Guinea-NE of              "; 
             volcESP_Code(i)="M0"; volcID(i)="0501-011"; volcName(i)="Boisa                          "; 
      i=i+1; volcLat(i)= -3.613; volcLon(i)=144.818; volcElev(i)=  685;  volcLoc(i)="New Guinea-NE of              "; 
             volcESP_Code(i)="S2"; volcID(i)="0501-01="; volcName(i)="Bam                            "; 
      i=i+1; volcLat(i)= -4.080; volcLon(i)=145.037; volcElev(i)= 1807;  volcLoc(i)="New Guinea-NE of              "; 
             volcESP_Code(i)="M1"; volcID(i)="0501-02="; volcName(i)="Manam                          "; 
      i=i+1; volcLat(i)= -4.649; volcLon(i)=145.964; volcElev(i)= 1839;  volcLoc(i)="New Guinea-NE of              "; 
             volcESP_Code(i)="S1"; volcID(i)="0501-03="; volcName(i)="Karkar                         "; 
      i=i+1; volcLat(i)= -4.900; volcLon(i)=146.750; volcElev(i)=    0;  volcLoc(i)="New Guinea-NE of              "; 
             volcESP_Code(i)="U0"; volcID(i)="0501-041"; volcName(i)="Yomba                          "; 
      i=i+1; volcLat(i)= -4.311; volcLon(i)=146.256; volcElev(i)=-2000;  volcLoc(i)="New Guinea-NE of              "; 
             volcESP_Code(i)="U0"; volcID(i)="0501-04="; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= -5.358; volcLon(i)=147.120; volcElev(i)= 1280;  volcLoc(i)="New Guinea-NE of              "; 
             volcESP_Code(i)="S1"; volcID(i)="0501-05="; volcName(i)="Long Island                    "; 
      i=i+1; volcLat(i)= -5.589; volcLon(i)=147.875; volcElev(i)= 1548;  volcLoc(i)="New Guinea-NE of              "; 
             volcESP_Code(i)="M0"; volcID(i)="0501-06="; volcName(i)="Umboi                          "; 
      i=i+1; volcLat(i)= -5.520; volcLon(i)=148.121; volcElev(i)=  140;  volcLoc(i)="New Guinea-NE of              "; 
             volcESP_Code(i)="M1"; volcID(i)="0501-07="; volcName(i)="Ritter Island                  "; 
      i=i+1; volcLat(i)= -5.414; volcLon(i)=148.094; volcElev(i)=  992;  volcLoc(i)="New Guinea-NE of              "; 
             volcESP_Code(i)="M0"; volcID(i)="0501-08="; volcName(i)="Sakar                          "; 
      i=i+1; volcLat(i)= -5.200; volcLon(i)=148.570; volcElev(i)=    0;  volcLoc(i)="New Britain-SW Pac            "; 
             volcESP_Code(i)="U0"; volcID(i)="0502-001"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= -5.525; volcLon(i)=148.420; volcElev(i)= 1330;  volcLoc(i)="New Britain-SW Pac            "; 
             volcESP_Code(i)="M1"; volcID(i)="0502-01="; volcName(i)="Langila                        "; 
      i=i+1; volcLat(i)= -4.630; volcLon(i)=149.350; volcElev(i)=  179;  volcLoc(i)="New Britain-SW Pac            "; 
             volcESP_Code(i)="M0"; volcID(i)="0502-021"; volcName(i)="Mundua                         "; 
      i=i+1; volcLat(i)= -4.692; volcLon(i)=149.500; volcElev(i)=  368;  volcLoc(i)="New Britain-SW Pac            "; 
             volcESP_Code(i)="S0"; volcID(i)="0502-03="; volcName(i)="Garove                         "; 
      i=i+1; volcLat(i)= -5.056; volcLon(i)=150.108; volcElev(i)=  400;  volcLoc(i)="New Britain-SW Pac            "; 
             volcESP_Code(i)="S0"; volcID(i)="0502-04="; volcName(i)="Dakataua                       "; 
      i=i+1; volcLat(i)= -5.150; volcLon(i)=150.030; volcElev(i)= 1155;  volcLoc(i)="New Britain-SW Pac            "; 
             volcESP_Code(i)="S0"; volcID(i)="0502-05="; volcName(i)="Bola                           "; 
      i=i+1; volcLat(i)= -5.300; volcLon(i)=150.070; volcElev(i)=  565;  volcLoc(i)="New Britain-SW Pac            "; 
             volcESP_Code(i)="S0"; volcID(i)="0502-06="; volcName(i)="Garua Harbour                  "; 
      i=i+1; volcLat(i)= -5.468; volcLon(i)=150.507; volcElev(i)=  805;  volcLoc(i)="New Britain-SW Pac            "; 
             volcESP_Code(i)="S0"; volcID(i)="0502-071"; volcName(i)="Lolo                           "; 
      i=i+1; volcLat(i)= -5.450; volcLon(i)=150.030; volcElev(i)=  564;  volcLoc(i)="New Britain-SW Pac            "; 
             volcESP_Code(i)="S0"; volcID(i)="0502-07="; volcName(i)="Garbuna Group                  "; 
      i=i+1; volcLat(i)= -5.580; volcLon(i)=150.520; volcElev(i)=  742;  volcLoc(i)="New Britain-SW Pac            "; 
             volcESP_Code(i)="S1"; volcID(i)="0502-08="; volcName(i)="Pago                           "; 
      i=i+1; volcLat(i)= -5.500; volcLon(i)=150.942; volcElev(i)=  610;  volcLoc(i)="New Britain-SW Pac            "; 
             volcESP_Code(i)="S0"; volcID(i)="0502-09="; volcName(i)="Sulu Range                     "; 
      i=i+1; volcLat(i)= -5.330; volcLon(i)=151.100; volcElev(i)= 1148;  volcLoc(i)="New Britain-SW Pac            "; 
             volcESP_Code(i)="S0"; volcID(i)="0502-10="; volcName(i)="Hargy                          "; 
      i=i+1; volcLat(i)= -5.200; volcLon(i)=151.230; volcElev(i)= 2248;  volcLoc(i)="New Britain-SW Pac            "; 
             volcESP_Code(i)="S0"; volcID(i)="0502-11="; volcName(i)="Bamus                          "; 
      i=i+1; volcLat(i)= -5.050; volcLon(i)=151.330; volcElev(i)= 2334;  volcLoc(i)="New Britain-SW Pac            "; 
             volcESP_Code(i)="S2"; volcID(i)="0502-12="; volcName(i)="Ulawun                         "; 
      i=i+1; volcLat(i)= -4.750; volcLon(i)=150.850; volcElev(i)=    0;  volcLoc(i)="New Britain-SW Pac            "; 
             volcESP_Code(i)="U0"; volcID(i)="0502-131"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= -4.920; volcLon(i)=151.158; volcElev(i)=  858;  volcLoc(i)="New Britain-SW Pac            "; 
             volcESP_Code(i)="M3"; volcID(i)="0502-13="; volcName(i)="Lolobau                        "; 
      i=i+1; volcLat(i)= -4.271; volcLon(i)=152.203; volcElev(i)=  688;  volcLoc(i)="New Britain-SW Pac            "; 
             volcESP_Code(i)="S1"; volcID(i)="0502-14="; volcName(i)="Rabaul                         "; 
      i=i+1; volcLat(i)= -4.120; volcLon(i)=152.200; volcElev(i)=  200;  volcLoc(i)="New Britain-SW Pac            "; 
             volcESP_Code(i)="S0"; volcID(i)="0502-15-"; volcName(i)="Tavui                          "; 
      i=i+1; volcLat(i)= -5.900; volcLon(i)=143.150; volcElev(i)= 3568;  volcLoc(i)="New Guinea                    "; 
             volcESP_Code(i)="S0"; volcID(i)="0503-00-"; volcName(i)="Doma Peaks                     "; 
      i=i+1; volcLat(i)= -6.580; volcLon(i)=145.080; volcElev(i)= 3233;  volcLoc(i)="New Guinea                    "; 
             volcESP_Code(i)="M0"; volcID(i)="0503-001"; volcName(i)="Crater Mountain                "; 
      i=i+1; volcLat(i)= -7.050; volcLon(i)=145.858; volcElev(i)= 3384;  volcLoc(i)="New Guinea                    "; 
             volcESP_Code(i)="S0"; volcID(i)="0503-002"; volcName(i)="Yelia                          "; 
      i=i+1; volcLat(i)= -7.330; volcLon(i)=146.708; volcElev(i)= 1500;  volcLoc(i)="New Guinea                    "; 
             volcESP_Code(i)="S0"; volcID(i)="0503-003"; volcName(i)="Koranga                        "; 
      i=i+1; volcLat(i)= -9.200; volcLon(i)=147.570; volcElev(i)=  850;  volcLoc(i)="New Guinea                    "; 
             volcESP_Code(i)="M0"; volcID(i)="0503-004"; volcName(i)="Madilogo                       "; 
      i=i+1; volcLat(i)= -9.000; volcLon(i)=148.370; volcElev(i)= 1915;  volcLoc(i)="New Guinea                    "; 
             volcESP_Code(i)="S0"; volcID(i)="0503-011"; volcName(i)="Hydrographers Range            "; 
      i=i+1; volcLat(i)= -8.950; volcLon(i)=148.150; volcElev(i)= 1680;  volcLoc(i)="New Guinea                    "; 
             volcESP_Code(i)="S2"; volcID(i)="0503-01="; volcName(i)="Lamington                      "; 
      i=i+1; volcLat(i)= -9.080; volcLon(i)=148.330; volcElev(i)= 1342;  volcLoc(i)="New Guinea                    "; 
             volcESP_Code(i)="M0"; volcID(i)="0503-021"; volcName(i)="Managlase Plateau              "; 
      i=i+1; volcLat(i)= -9.308; volcLon(i)=148.130; volcElev(i)=  808;  volcLoc(i)="New Guinea                    "; 
             volcESP_Code(i)="S0"; volcID(i)="0503-02="; volcName(i)="Musa River                     "; 
      i=i+1; volcLat(i)= -9.480; volcLon(i)=149.130; volcElev(i)=  370;  volcLoc(i)="New Guinea                    "; 
             volcESP_Code(i)="S0"; volcID(i)="0503-031"; volcName(i)="Sessagara                      "; 
      i=i+1; volcLat(i)= -9.200; volcLon(i)=149.070; volcElev(i)= 1925;  volcLoc(i)="New Guinea                    "; 
             volcESP_Code(i)="S1"; volcID(i)="0503-03="; volcName(i)="Victory                        "; 
      i=i+1; volcLat(i)= -9.480; volcLon(i)=150.350; volcElev(i)=  220;  volcLoc(i)="D'Entrecasteaux Is            "; 
             volcESP_Code(i)="S0"; volcID(i)="0503-041"; volcName(i)="Goodenough                     "; 
      i=i+1; volcLat(i)= -9.570; volcLon(i)=149.075; volcElev(i)=  640;  volcLoc(i)="New Guinea                    "; 
             volcESP_Code(i)="S2"; volcID(i)="0503-04="; volcName(i)="Waiowa                         "; 
      i=i+1; volcLat(i)= -9.520; volcLon(i)=150.530; volcElev(i)=  200;  volcLoc(i)="D'Entrecasteaux Is            "; 
             volcESP_Code(i)="S0"; volcID(i)="0503-05="; volcName(i)="Iamalele                       "; 
      i=i+1; volcLat(i)= -9.620; volcLon(i)=150.880; volcElev(i)=  500;  volcLoc(i)="D'Entrecasteaux Is            "; 
             volcESP_Code(i)="S0"; volcID(i)="0503-06="; volcName(i)="Dawson Strait Group            "; 
      i=i+1; volcLat(i)= -3.125; volcLon(i)=152.642; volcElev(i)=  700;  volcLoc(i)="New Ireland-SW Pacific        "; 
             volcESP_Code(i)="M0"; volcID(i)="0504-01="; volcName(i)="Lihir                          "; 
      i=i+1; volcLat(i)= -4.080; volcLon(i)=153.650; volcElev(i)=  450;  volcLoc(i)="New Ireland-SW Pacific        "; 
             volcESP_Code(i)="S0"; volcID(i)="0504-02="; volcName(i)="Ambitle                        "; 
      i=i+1; volcLat(i)= -5.830; volcLon(i)=154.930; volcElev(i)= 2200;  volcLoc(i)="Bougainville-SW Pacific       "; 
             volcESP_Code(i)="S0"; volcID(i)="0505-00-"; volcName(i)="Tore                           "; 
      i=i+1; volcLat(i)= -6.092; volcLon(i)=155.225; volcElev(i)= 1544;  volcLoc(i)="Bougainville-SW Pacific       "; 
             volcESP_Code(i)="S0"; volcID(i)="0505-011"; volcName(i)="Billy Mitchell                 "; 
      i=i+1; volcLat(i)= -5.920; volcLon(i)=154.980; volcElev(i)= 2715;  volcLoc(i)="Bougainville-SW Pacific       "; 
             volcESP_Code(i)="S0"; volcID(i)="0505-01="; volcName(i)="Balbi                          "; 
      i=i+1; volcLat(i)= -6.442; volcLon(i)=155.608; volcElev(i)= 2210;  volcLoc(i)="Bougainville-SW Pacific       "; 
             volcESP_Code(i)="S0"; volcID(i)="0505-021"; volcName(i)="Takuan Group                   "; 
      i=i+1; volcLat(i)= -6.140; volcLon(i)=155.195; volcElev(i)= 1750;  volcLoc(i)="Bougainville-SW Pacific       "; 
             volcESP_Code(i)="S1"; volcID(i)="0505-02="; volcName(i)="Bagana                         "; 
      i=i+1; volcLat(i)= -6.520; volcLon(i)=155.620; volcElev(i)= 1887;  volcLoc(i)="Bougainville-SW Pacific       "; 
             volcESP_Code(i)="S0"; volcID(i)="0505-03="; volcName(i)="Loloru                         "; 
      i=i+1; volcLat(i)= -8.750; volcLon(i)=157.030; volcElev(i)= -700;  volcLoc(i)="Solomon Is-SW Pacific         "; 
             volcESP_Code(i)="U0"; volcID(i)="0505-052"; volcName(i)="Kana Keoki                     "; 
      i=i+1; volcLat(i)= -8.830; volcLon(i)=157.170; volcElev(i)=    0;  volcLoc(i)="Solomon Is-SW Pacific         "; 
             volcESP_Code(i)="U0"; volcID(i)="0505-053"; volcName(i)="Coleman Seamount               "; 
      i=i+1; volcLat(i)= -8.292; volcLon(i)=156.520; volcElev(i)=  335;  volcLoc(i)="Solomon Is-SW Pacific         "; 
             volcESP_Code(i)="S0"; volcID(i)="0505-05="; volcName(i)="Simbo                          "; 
      i=i+1; volcLat(i)= -8.920; volcLon(i)=158.030; volcElev(i)= -240;  volcLoc(i)="Solomon Is-SW Pacific         "; 
             volcESP_Code(i)="U0"; volcID(i)="0505-061"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= -9.350; volcLon(i)=159.730; volcElev(i)= 1000;  volcLoc(i)="Solomon Is-SW Pacific         "; 
             volcESP_Code(i)="S0"; volcID(i)="0505-062"; volcName(i)="Gallego                        "; 
      i=i+1; volcLat(i)= -9.020; volcLon(i)=157.950; volcElev(i)=  -20;  volcLoc(i)="Solomon Is-SW Pacific         "; 
             volcESP_Code(i)="M1"; volcID(i)="0505-06="; volcName(i)="Kavachi                        "; 
      i=i+1; volcLat(i)= -9.130; volcLon(i)=159.820; volcElev(i)=  485;  volcLoc(i)="Solomon Is-SW Pacific         "; 
             volcESP_Code(i)="S2"; volcID(i)="0505-07="; volcName(i)="Savo                           "; 
      i=i+1; volcLat(i)=-10.380; volcLon(i)=165.800; volcElev(i)=  851;  volcLoc(i)="Santa Cruz Is-SW Pacific      "; 
             volcESP_Code(i)="S1"; volcID(i)="0506-01="; volcName(i)="Tinakula                       "; 
      i=i+1; volcLat(i)=-13.670; volcLon(i)=167.670; volcElev(i)=  411;  volcLoc(i)="Vanuatu-SW Pacific            "; 
             volcESP_Code(i)="M0"; volcID(i)="0507-001"; volcName(i)="Motlav                         "; 
      i=i+1; volcLat(i)=-13.800; volcLon(i)=167.470; volcElev(i)=  921;  volcLoc(i)="Vanuatu-SW Pacific            "; 
             volcESP_Code(i)="S0"; volcID(i)="0507-01="; volcName(i)="Suretamatai                    "; 
      i=i+1; volcLat(i)=-14.450; volcLon(i)=168.050; volcElev(i)= 1028;  volcLoc(i)="Vanuatu-SW Pacific            "; 
             volcESP_Code(i)="M0"; volcID(i)="0507-021"; volcName(i)="Mere Lava                      "; 
      i=i+1; volcLat(i)=-14.270; volcLon(i)=167.500; volcElev(i)=  797;  volcLoc(i)="Vanuatu-SW Pacific            "; 
             volcESP_Code(i)="S1"; volcID(i)="0507-02="; volcName(i)="Gaua                           "; 
      i=i+1; volcLat(i)=-15.400; volcLon(i)=167.830; volcElev(i)= 1496;  volcLoc(i)="Vanuatu-SW Pacific            "; 
             volcESP_Code(i)="M0"; volcID(i)="0507-03="; volcName(i)="Aoba                           "; 
      i=i+1; volcLat(i)=-16.250; volcLon(i)=168.120; volcElev(i)= 1334;  volcLoc(i)="Vanuatu-SW Pacific            "; 
             volcESP_Code(i)="M1"; volcID(i)="0507-04="; volcName(i)="Ambrym                         "; 
      i=i+1; volcLat(i)=-16.507; volcLon(i)=168.346; volcElev(i)= 1413;  volcLoc(i)="Vanuatu-SW Pacific            "; 
             volcESP_Code(i)="M2"; volcID(i)="0507-05="; volcName(i)="Lopevi                         "; 
      i=i+1; volcLat(i)=-16.680; volcLon(i)=168.370; volcElev(i)=  833;  volcLoc(i)="Vanuatu-SW Pacific            "; 
             volcESP_Code(i)="M0"; volcID(i)="0507-06="; volcName(i)="Epi                            "; 
      i=i+1; volcLat(i)=-16.829; volcLon(i)=168.536; volcElev(i)=   -2;  volcLoc(i)="Vanuatu-SW Pacific            "; 
             volcESP_Code(i)="M0"; volcID(i)="0507-07="; volcName(i)="Kuwae                          "; 
      i=i+1; volcLat(i)=-16.992; volcLon(i)=168.592; volcElev(i)=  216;  volcLoc(i)="Vanuatu-SW Pacific            "; 
             volcESP_Code(i)="M0"; volcID(i)="0507-08-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=-17.470; volcLon(i)=168.353; volcElev(i)=  594;  volcLoc(i)="Vanuatu-SW Pacific            "; 
             volcESP_Code(i)="M0"; volcID(i)="0507-081"; volcName(i)="North Vate                     "; 
      i=i+1; volcLat(i)=-18.750; volcLon(i)=169.230; volcElev(i)=  837;  volcLoc(i)="Vanuatu-SW Pacific            "; 
             volcESP_Code(i)="M0"; volcID(i)="0507-09="; volcName(i)="Traitor's Head                 "; 
      i=i+1; volcLat(i)=-19.530; volcLon(i)=169.442; volcElev(i)=  361;  volcLoc(i)="Vanuatu-SW Pacific            "; 
             volcESP_Code(i)="M1"; volcID(i)="0507-10="; volcName(i)="Yasur                          "; 
      i=i+1; volcLat(i)=-20.200; volcLon(i)=169.780; volcElev(i)=  852;  volcLoc(i)="Vanuatu-SW Pacific            "; 
             volcESP_Code(i)="M0"; volcID(i)="0507-11-"; volcName(i)="Aneityum                       "; 
      i=i+1; volcLat(i)=-20.980; volcLon(i)=170.280; volcElev(i)=  -80;  volcLoc(i)="SW Pacific                    "; 
             volcESP_Code(i)="U0"; volcID(i)="0508-001"; volcName(i)="Eastern Gemini Seamount        "; 
      i=i+1; volcLat(i)=-22.330; volcLon(i)=171.320; volcElev(i)=  177;  volcLoc(i)="SW Pacific                    "; 
             volcESP_Code(i)="S0"; volcID(i)="0508-01="; volcName(i)="Matthew Island                 "; 
      i=i+1; volcLat(i)=-22.400; volcLon(i)=172.050; volcElev(i)=  297;  volcLoc(i)="SW Pacific                    "; 
             volcESP_Code(i)="S0"; volcID(i)="0508-02="; volcName(i)="Hunter Island                  "; 
      i=i+1; volcLat(i)=-25.780; volcLon(i)=168.630; volcElev(i)=-2400;  volcLoc(i)="SW Pacific                    "; 
             volcESP_Code(i)="U0"; volcID(i)="0508-03-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=-37.770; volcLon(i)=142.500; volcElev(i)= 1011;  volcLoc(i)="Australia                     "; 
             volcESP_Code(i)="M0"; volcID(i)="0509-01-"; volcName(i)="Newer Volcanics Prov           "; 
      i=i+1; volcLat(i)= 13.430; volcLon(i)= 94.280; volcElev(i)=  710;  volcLoc(i)="Andaman Is-Indian O           "; 
             volcESP_Code(i)="S0"; volcID(i)="0600-001"; volcName(i)="Narcondum                      "; 
      i=i+1; volcLat(i)= 12.278; volcLon(i)= 93.858; volcElev(i)=  354;  volcLoc(i)="Andaman Is-Indian O           "; 
             volcESP_Code(i)="S2"; volcID(i)="0600-01="; volcName(i)="Barren Island                  "; 
      i=i+1; volcLat(i)=  5.448; volcLon(i)= 95.658; volcElev(i)= 1810;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0601-02="; volcName(i)="Seulawah Agam                  "; 
      i=i+1; volcLat(i)=  4.914; volcLon(i)= 96.329; volcElev(i)= 2801;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S1"; volcID(i)="0601-03="; volcName(i)="Peuet Sague                    "; 
      i=i+1; volcLat(i)=  4.769; volcLon(i)= 96.821; volcElev(i)= 2617;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S1"; volcID(i)="0601-05="; volcName(i)="Telong, Bur ni                 "; 
      i=i+1; volcLat(i)=  3.230; volcLon(i)= 98.520; volcElev(i)= 2212;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0601-07="; volcName(i)="Sibayak                        "; 
      i=i+1; volcLat(i)=  3.170; volcLon(i)= 98.392; volcElev(i)= 2460;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0601-08="; volcName(i)="Sinabung                       "; 
      i=i+1; volcLat(i)=  2.580; volcLon(i)= 98.830; volcElev(i)= 2157;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0601-09="; volcName(i)="Toba                           "; 
      i=i+1; volcLat(i)=  2.158; volcLon(i)= 98.930; volcElev(i)= 1505;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0601-101"; volcName(i)="Imun                           "; 
      i=i+1; volcLat(i)=  1.478; volcLon(i)= 99.209; volcElev(i)= 1862;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0601-111"; volcName(i)="Lubukraya                      "; 
      i=i+1; volcLat(i)=  1.556; volcLon(i)= 99.255; volcElev(i)= 1819;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0601-11="; volcName(i)="Sibualbuali                    "; 
      i=i+1; volcLat(i)=  0.686; volcLon(i)= 99.539; volcElev(i)= 2145;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S1"; volcID(i)="0601-12="; volcName(i)="Sorikmarapi                    "; 
      i=i+1; volcLat(i)=  0.080; volcLon(i)=100.200; volcElev(i)=    0;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0601-131"; volcName(i)="Sarik-Gajah                    "; 
      i=i+1; volcLat(i)=  0.079; volcLon(i)= 99.980; volcElev(i)= 2919;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0601-13="; volcName(i)="Talakmau                       "; 
      i=i+1; volcLat(i)= -0.381; volcLon(i)=100.473; volcElev(i)= 2891;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S1"; volcID(i)="0601-14="; volcName(i)="Marapi                         "; 
      i=i+1; volcLat(i)= -0.433; volcLon(i)=100.317; volcElev(i)= 2438;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S1"; volcID(i)="0601-15="; volcName(i)="Tandikat                       "; 
      i=i+1; volcLat(i)= -0.978; volcLon(i)=100.679; volcElev(i)= 2597;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S1"; volcID(i)="0601-16="; volcName(i)="Talang                         "; 
      i=i+1; volcLat(i)= -2.274; volcLon(i)=101.483; volcElev(i)= 2151;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0601-171"; volcName(i)="Kunyit                         "; 
      i=i+1; volcLat(i)= -2.330; volcLon(i)=101.600; volcElev(i)= 2021;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0601-172"; volcName(i)="Hutapanjang                    "; 
      i=i+1; volcLat(i)= -1.697; volcLon(i)=101.264; volcElev(i)= 3800;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S1"; volcID(i)="0601-17="; volcName(i)="Kerinci                        "; 
      i=i+1; volcLat(i)= -2.414; volcLon(i)=101.728; volcElev(i)= 2507;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S1"; volcID(i)="0601-18="; volcName(i)="Sumbing                        "; 
      i=i+1; volcLat(i)= -2.820; volcLon(i)=102.020; volcElev(i)=    0;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0601-191"; volcName(i)="Pendan                         "; 
      i=i+1; volcLat(i)= -2.820; volcLon(i)=102.180; volcElev(i)= 1958;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="M0"; volcID(i)="0601-20="; volcName(i)="Belirang-Beriti                "; 
      i=i+1; volcLat(i)= -3.380; volcLon(i)=102.370; volcElev(i)= 2467;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="M0"; volcID(i)="0601-21="; volcName(i)="Daun, Bukit                    "; 
      i=i+1; volcLat(i)= -3.520; volcLon(i)=102.620; volcElev(i)= 1952;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S1"; volcID(i)="0601-22="; volcName(i)="Kaba                           "; 
      i=i+1; volcLat(i)= -4.270; volcLon(i)=103.300; volcElev(i)= 2817;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0601-231"; volcName(i)="Patah                          "; 
      i=i+1; volcLat(i)= -4.030; volcLon(i)=103.130; volcElev(i)= 3173;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S1"; volcID(i)="0601-23="; volcName(i)="Dempo                          "; 
      i=i+1; volcLat(i)= -4.220; volcLon(i)=103.620; volcElev(i)= 2055;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0601-24="; volcName(i)="Lumut Balai, Bukit             "; 
      i=i+1; volcLat(i)= -4.830; volcLon(i)=103.920; volcElev(i)= 1881;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0601-251"; volcName(i)="Ranau                          "; 
      i=i+1; volcLat(i)= -4.430; volcLon(i)=103.670; volcElev(i)= 1899;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0601-25="; volcName(i)="Besar                          "; 
      i=i+1; volcLat(i)= -5.120; volcLon(i)=104.320; volcElev(i)= 1719;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0601-26="; volcName(i)="Sekincau Belirang              "; 
      i=i+1; volcLat(i)= -5.250; volcLon(i)=104.270; volcElev(i)= 1000;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0601-27="; volcName(i)="Suoh                           "; 
      i=i+1; volcLat(i)= -5.350; volcLon(i)=104.600; volcElev(i)= 1040;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0601-28="; volcName(i)="Hulubelu                       "; 
      i=i+1; volcLat(i)= -5.780; volcLon(i)=105.625; volcElev(i)= 1281;  volcLoc(i)="Sumatra                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0601-29="; volcName(i)="Rajabasa                       "; 
      i=i+1; volcLat(i)= -6.102; volcLon(i)=105.423; volcElev(i)=  813;  volcLoc(i)="Indonesia                     "; 
             volcESP_Code(i)="S1"; volcID(i)="0602-00="; volcName(i)="Krakatau                       "; 
      i=i+1; volcLat(i)= -6.342; volcLon(i)=105.975; volcElev(i)= 1346;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0603-01="; volcName(i)="Pulosari                       "; 
      i=i+1; volcLat(i)= -6.270; volcLon(i)=106.042; volcElev(i)= 1778;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0603-02="; volcName(i)="Karang                         "; 
      i=i+1; volcLat(i)= -6.750; volcLon(i)=106.700; volcElev(i)= 1699;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0603-04="; volcName(i)="Perbakti-Gagak                 "; 
      i=i+1; volcLat(i)= -6.720; volcLon(i)=106.730; volcElev(i)= 2211;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S1"; volcID(i)="0603-05="; volcName(i)="Salak                          "; 
      i=i+1; volcLat(i)= -6.780; volcLon(i)=106.980; volcElev(i)= 2958;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S2"; volcID(i)="0603-06="; volcName(i)="Gede                           "; 
      i=i+1; volcLat(i)= -7.160; volcLon(i)=107.400; volcElev(i)= 2434;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0603-07="; volcName(i)="Patuha                         "; 
      i=i+1; volcLat(i)= -7.130; volcLon(i)=107.650; volcElev(i)= 2343;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0603-081"; volcName(i)="Malabar                        "; 
      i=i+1; volcLat(i)= -7.208; volcLon(i)=107.630; volcElev(i)= 2182;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0603-08="; volcName(i)="Wayang-Windu                   "; 
      i=i+1; volcLat(i)= -6.770; volcLon(i)=107.600; volcElev(i)= 2084;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S1"; volcID(i)="0603-09="; volcName(i)="Tangkubanparahu                "; 
      i=i+1; volcLat(i)= -7.320; volcLon(i)=107.730; volcElev(i)= 2665;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S1"; volcID(i)="0603-10="; volcName(i)="Papandayan                     "; 
      i=i+1; volcLat(i)= -7.230; volcLon(i)=107.720; volcElev(i)= 2608;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0603-11="; volcName(i)="Kendang                        "; 
      i=i+1; volcLat(i)= -6.770; volcLon(i)=107.950; volcElev(i)= 1684;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0603-131"; volcName(i)="Tampomas                       "; 
      i=i+1; volcLat(i)= -7.143; volcLon(i)=107.840; volcElev(i)= 2249;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S1"; volcID(i)="0603-13="; volcName(i)="Guntur                         "; 
      i=i+1; volcLat(i)= -7.250; volcLon(i)=108.058; volcElev(i)= 2168;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="M2"; volcID(i)="0603-14="; volcName(i)="Galunggung                     "; 
      i=i+1; volcLat(i)= -7.208; volcLon(i)=108.070; volcElev(i)= 2201;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="M0"; volcID(i)="0603-15="; volcName(i)="Talagabodas                    "; 
      i=i+1; volcLat(i)= -7.120; volcLon(i)=108.080; volcElev(i)= 1155;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0603-16="; volcName(i)="Karaha, Kawah                  "; 
      i=i+1; volcLat(i)= -6.892; volcLon(i)=108.400; volcElev(i)= 3078;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S1"; volcID(i)="0603-17="; volcName(i)="Cereme                         "; 
      i=i+1; volcLat(i)= -7.242; volcLon(i)=109.208; volcElev(i)= 3428;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="M1"; volcID(i)="0603-18="; volcName(i)="Slamet                         "; 
      i=i+1; volcLat(i)= -7.200; volcLon(i)=109.920; volcElev(i)= 2565;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S1"; volcID(i)="0603-20="; volcName(i)="Dieng Volc Complex             "; 
      i=i+1; volcLat(i)= -7.300; volcLon(i)=109.992; volcElev(i)= 3136;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S1"; volcID(i)="0603-21="; volcName(i)="Sundoro                        "; 
      i=i+1; volcLat(i)= -7.384; volcLon(i)=110.070; volcElev(i)= 3371;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0603-22="; volcName(i)="Sumbing                        "; 
      i=i+1; volcLat(i)= -7.370; volcLon(i)=110.400; volcElev(i)= 1894;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0603-231"; volcName(i)="Telomoyo                       "; 
      i=i+1; volcLat(i)= -7.180; volcLon(i)=110.330; volcElev(i)= 2050;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0603-23="; volcName(i)="Ungaran                        "; 
      i=i+1; volcLat(i)= -7.450; volcLon(i)=110.430; volcElev(i)= 3145;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="M1"; volcID(i)="0603-24="; volcName(i)="Merbabu                        "; 
      i=i+1; volcLat(i)= -6.620; volcLon(i)=110.880; volcElev(i)= 1625;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="M0"; volcID(i)="0603-251"; volcName(i)="Muria                          "; 
      i=i+1; volcLat(i)= -7.542; volcLon(i)=110.442; volcElev(i)= 2968;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S9"; volcID(i)="0603-25="; volcName(i)="Merapi                         "; 
      i=i+1; volcLat(i)= -7.625; volcLon(i)=111.192; volcElev(i)= 3265;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S1"; volcID(i)="0603-26="; volcName(i)="Lawu                           "; 
      i=i+1; volcLat(i)= -7.808; volcLon(i)=111.758; volcElev(i)= 2563;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0603-27="; volcName(i)="Wilis                          "; 
      i=i+1; volcLat(i)= -7.920; volcLon(i)=112.450; volcElev(i)= 2651;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0603-281"; volcName(i)="Kawi-Butak                     "; 
      i=i+1; volcLat(i)= -7.930; volcLon(i)=112.308; volcElev(i)= 1731;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S2"; volcID(i)="0603-28="; volcName(i)="Kelut                          "; 
      i=i+1; volcLat(i)= -7.620; volcLon(i)=112.630; volcElev(i)= 1653;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0603-291"; volcName(i)="Penanggungan                   "; 
      i=i+1; volcLat(i)= -8.020; volcLon(i)=112.680; volcElev(i)=  680;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0603-292"; volcName(i)="Malang Plain                   "; 
      i=i+1; volcLat(i)= -7.725; volcLon(i)=112.580; volcElev(i)= 3339;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S1"; volcID(i)="0603-29="; volcName(i)="Arjuno-Welirang                "; 
      i=i+1; volcLat(i)= -8.108; volcLon(i)=112.920; volcElev(i)= 3676;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S1"; volcID(i)="0603-30="; volcName(i)="Semeru                         "; 
      i=i+1; volcLat(i)= -7.942; volcLon(i)=112.950; volcElev(i)= 2329;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S1"; volcID(i)="0603-31="; volcName(i)="Tengger Caldera                "; 
      i=i+1; volcLat(i)= -7.730; volcLon(i)=113.580; volcElev(i)=  539;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0603-321"; volcName(i)="Lurus                          "; 
      i=i+1; volcLat(i)= -7.979; volcLon(i)=113.342; volcElev(i)= 1651;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="M1"; volcID(i)="0603-32="; volcName(i)="Lamongan                       "; 
      i=i+1; volcLat(i)= -7.970; volcLon(i)=113.570; volcElev(i)= 3088;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0603-33="; volcName(i)="Iyang-Argapura                 "; 
      i=i+1; volcLat(i)= -8.125; volcLon(i)=114.042; volcElev(i)= 3332;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S1"; volcID(i)="0603-34="; volcName(i)="Raung                          "; 
      i=i+1; volcLat(i)= -7.850; volcLon(i)=114.370; volcElev(i)= 1247;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S0"; volcID(i)="0603-351"; volcName(i)="Baluran                        "; 
      i=i+1; volcLat(i)= -8.058; volcLon(i)=114.242; volcElev(i)= 2799;  volcLoc(i)="Java                          "; 
             volcESP_Code(i)="S1"; volcID(i)="0603-35="; volcName(i)="Ijen                           "; 
      i=i+1; volcLat(i)= -8.280; volcLon(i)=115.130; volcElev(i)= 2276;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="S0"; volcID(i)="0604-001"; volcName(i)="Bratan                         "; 
      i=i+1; volcLat(i)= -8.242; volcLon(i)=115.375; volcElev(i)= 1717;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="S1"; volcID(i)="0604-01="; volcName(i)="Batur                          "; 
      i=i+1; volcLat(i)= -8.342; volcLon(i)=115.508; volcElev(i)= 3142;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="S0"; volcID(i)="0604-02="; volcName(i)="Agung                          "; 
      i=i+1; volcLat(i)= -8.420; volcLon(i)=116.470; volcElev(i)= 3726;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="S1"; volcID(i)="0604-03="; volcName(i)="Rinjani                        "; 
      i=i+1; volcLat(i)= -8.250; volcLon(i)=118.000; volcElev(i)= 2850;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="S0"; volcID(i)="0604-04="; volcName(i)="Tambora                        "; 
      i=i+1; volcLat(i)= -8.200; volcLon(i)=119.070; volcElev(i)= 1949;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="M2"; volcID(i)="0604-05="; volcName(i)="Sangeang Api                   "; 
      i=i+1; volcLat(i)= -8.720; volcLon(i)=120.020; volcElev(i)=  903;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="S0"; volcID(i)="0604-06="; volcName(i)="Sano, Wai                      "; 
      i=i+1; volcLat(i)= -8.620; volcLon(i)=120.520; volcElev(i)= 2350;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="S0"; volcID(i)="0604-071"; volcName(i)="Ranakah                        "; 
      i=i+1; volcLat(i)= -8.680; volcLon(i)=120.480; volcElev(i)= 1675;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="S0"; volcID(i)="0604-07="; volcName(i)="Poco Leok                      "; 
      i=i+1; volcLat(i)= -8.875; volcLon(i)=120.950; volcElev(i)= 2245;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="S0"; volcID(i)="0604-08="; volcName(i)="Inierie                        "; 
      i=i+1; volcLat(i)= -8.730; volcLon(i)=120.980; volcElev(i)= 1559;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="S1"; volcID(i)="0604-09="; volcName(i)="Inielika                       "; 
      i=i+1; volcLat(i)= -8.820; volcLon(i)=121.180; volcElev(i)= 2124;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="S1"; volcID(i)="0604-10="; volcName(i)="Ebulobo                        "; 
      i=i+1; volcLat(i)= -8.897; volcLon(i)=121.645; volcElev(i)=  637;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="M1"; volcID(i)="0604-11="; volcName(i)="Iya                            "; 
      i=i+1; volcLat(i)= -8.792; volcLon(i)=121.770; volcElev(i)= 1500;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="S0"; volcID(i)="0604-12="; volcName(i)="Sukaria Caldera                "; 
      i=i+1; volcLat(i)= -8.720; volcLon(i)=121.780; volcElev(i)=  750;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="S0"; volcID(i)="0604-13="; volcName(i)="Ndete Napu                     "; 
      i=i+1; volcLat(i)= -8.770; volcLon(i)=121.820; volcElev(i)= 1639;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="S1"; volcID(i)="0604-14="; volcName(i)="Kelimutu                       "; 
      i=i+1; volcLat(i)= -8.320; volcLon(i)=121.708; volcElev(i)=  875;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="S2"; volcID(i)="0604-15="; volcName(i)="Paluweh                        "; 
      i=i+1; volcLat(i)= -8.670; volcLon(i)=122.450; volcElev(i)= 1703;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="S1"; volcID(i)="0604-16="; volcName(i)="Egon                           "; 
      i=i+1; volcLat(i)= -8.478; volcLon(i)=122.671; volcElev(i)= 1100;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="S0"; volcID(i)="0604-17="; volcName(i)="Ilimuda                        "; 
      i=i+1; volcLat(i)= -8.542; volcLon(i)=122.775; volcElev(i)= 1703;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="S1"; volcID(i)="0604-18="; volcName(i)="Lewotobi                       "; 
      i=i+1; volcLat(i)= -8.358; volcLon(i)=122.842; volcElev(i)= 1117;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="S1"; volcID(i)="0604-20="; volcName(i)="Leroboleng                     "; 
      i=i+1; volcLat(i)= -8.342; volcLon(i)=123.258; volcElev(i)= 1659;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="M1"; volcID(i)="0604-22="; volcName(i)="Iliboleng                      "; 
      i=i+1; volcLat(i)= -8.272; volcLon(i)=123.505; volcElev(i)= 1423;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="S1"; volcID(i)="0604-23="; volcName(i)="Lewotolo                       "; 
      i=i+1; volcLat(i)= -8.550; volcLon(i)=123.380; volcElev(i)= 1018;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="M0"; volcID(i)="0604-24="; volcName(i)="Ililabalekan                   "; 
      i=i+1; volcLat(i)= -8.530; volcLon(i)=123.570; volcElev(i)= 1018;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="M1"; volcID(i)="0604-25="; volcName(i)="Iliwerung                      "; 
      i=i+1; volcLat(i)= -7.792; volcLon(i)=123.579; volcElev(i)=  748;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="M1"; volcID(i)="0604-26="; volcName(i)="Tara, Batu                     "; 
      i=i+1; volcLat(i)= -8.508; volcLon(i)=124.130; volcElev(i)=  862;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="S1"; volcID(i)="0604-27="; volcName(i)="Sirung                         "; 
      i=i+1; volcLat(i)= -7.530; volcLon(i)=123.950; volcElev(i)=-3800;  volcLoc(i)="Lesser Sunda Is               "; 
             volcESP_Code(i)="U0"; volcID(i)="0604-28="; volcName(i)="Yersey                         "; 
      i=i+1; volcLat(i)= -6.620; volcLon(i)=124.220; volcElev(i)=-2850;  volcLoc(i)="Banda Sea                     "; 
             volcESP_Code(i)="U0"; volcID(i)="0605-01="; volcName(i)="Emperor of China               "; 
      i=i+1; volcLat(i)= -6.600; volcLon(i)=124.675; volcElev(i)=-2285;  volcLoc(i)="Banda Sea                     "; 
             volcESP_Code(i)="U0"; volcID(i)="0605-02="; volcName(i)="Nieuwerkerk                    "; 
      i=i+1; volcLat(i)= -6.642; volcLon(i)=126.650; volcElev(i)=  282;  volcLoc(i)="Banda Sea                     "; 
             volcESP_Code(i)="S2"; volcID(i)="0605-03="; volcName(i)="Gunungapi Wetar                "; 
      i=i+1; volcLat(i)= -7.125; volcLon(i)=128.675; volcElev(i)=  868;  volcLoc(i)="Banda Sea                     "; 
             volcESP_Code(i)="S1"; volcID(i)="0605-04="; volcName(i)="Wurlali                        "; 
      i=i+1; volcLat(i)= -6.920; volcLon(i)=129.125; volcElev(i)=  655;  volcLoc(i)="Banda Sea                     "; 
             volcESP_Code(i)="S2"; volcID(i)="0605-05="; volcName(i)="Teon                           "; 
      i=i+1; volcLat(i)= -6.730; volcLon(i)=129.500; volcElev(i)=  781;  volcLoc(i)="Banda Sea                     "; 
             volcESP_Code(i)="S1"; volcID(i)="0605-06="; volcName(i)="Nila                           "; 
      i=i+1; volcLat(i)= -6.300; volcLon(i)=130.000; volcElev(i)=  641;  volcLoc(i)="Banda Sea                     "; 
             volcESP_Code(i)="S1"; volcID(i)="0605-07="; volcName(i)="Serua                          "; 
      i=i+1; volcLat(i)= -5.530; volcLon(i)=130.292; volcElev(i)=  282;  volcLoc(i)="Banda Sea                     "; 
             volcESP_Code(i)="S0"; volcID(i)="0605-08="; volcName(i)="Manuk                          "; 
      i=i+1; volcLat(i)= -4.525; volcLon(i)=129.871; volcElev(i)=  640;  volcLoc(i)="Banda Sea                     "; 
             volcESP_Code(i)="S1"; volcID(i)="0605-09="; volcName(i)="Banda Api                      "; 
      i=i+1; volcLat(i)= -0.170; volcLon(i)=121.608; volcElev(i)=  507;  volcLoc(i)="Sulawesi-Indonesia            "; 
             volcESP_Code(i)="S2"; volcID(i)="0606-01="; volcName(i)="Colo                           "; 
      i=i+1; volcLat(i)=  0.750; volcLon(i)=124.420; volcElev(i)= 1795;  volcLoc(i)="Sulawesi-Indonesia            "; 
             volcESP_Code(i)="S0"; volcID(i)="0606-02="; volcName(i)="Ambang                         "; 
      i=i+1; volcLat(i)=  1.108; volcLon(i)=124.730; volcElev(i)= 1784;  volcLoc(i)="Sulawesi-Indonesia            "; 
             volcESP_Code(i)="S1"; volcID(i)="0606-03="; volcName(i)="Soputan                        "; 
      i=i+1; volcLat(i)=  1.130; volcLon(i)=124.758; volcElev(i)= 1549;  volcLoc(i)="Sulawesi-Indonesia            "; 
             volcESP_Code(i)="S0"; volcID(i)="0606-04="; volcName(i)="Sempu                          "; 
      i=i+1; volcLat(i)=  1.230; volcLon(i)=124.830; volcElev(i)= 1202;  volcLoc(i)="Sulawesi-Indonesia            "; 
             volcESP_Code(i)="S0"; volcID(i)="0606-07-"; volcName(i)="Tondano Caldera                "; 
      i=i+1; volcLat(i)=  1.358; volcLon(i)=124.792; volcElev(i)= 1580;  volcLoc(i)="Sulawesi-Indonesia            "; 
             volcESP_Code(i)="S1"; volcID(i)="0606-10="; volcName(i)="Lokon-Empung                   "; 
      i=i+1; volcLat(i)=  1.358; volcLon(i)=124.858; volcElev(i)= 1324;  volcLoc(i)="Sulawesi-Indonesia            "; 
             volcESP_Code(i)="S1"; volcID(i)="0606-11="; volcName(i)="Mahawu                         "; 
      i=i+1; volcLat(i)=  1.470; volcLon(i)=125.030; volcElev(i)= 1995;  volcLoc(i)="Sulawesi-Indonesia            "; 
             volcESP_Code(i)="S0"; volcID(i)="0606-12="; volcName(i)="Klabat                         "; 
      i=i+1; volcLat(i)=  1.520; volcLon(i)=125.200; volcElev(i)= 1149;  volcLoc(i)="Sulawesi-Indonesia            "; 
             volcESP_Code(i)="S1"; volcID(i)="0606-13="; volcName(i)="Tongkoko                       "; 
      i=i+1; volcLat(i)=  2.300; volcLon(i)=125.370; volcElev(i)=  725;  volcLoc(i)="Sangihe Is-Indonesia          "; 
             volcESP_Code(i)="S1"; volcID(i)="0607-01="; volcName(i)="Ruang                          "; 
      i=i+1; volcLat(i)=  2.780; volcLon(i)=125.400; volcElev(i)= 1784;  volcLoc(i)="Sangihe Is-Indonesia          "; 
             volcESP_Code(i)="S1"; volcID(i)="0607-02="; volcName(i)="Karangetang                    "; 
      i=i+1; volcLat(i)=  3.138; volcLon(i)=125.491; volcElev(i)=   -5;  volcLoc(i)="Sangihe Is-Indonesia          "; 
             volcESP_Code(i)="S1"; volcID(i)="0607-03="; volcName(i)="Banua Wuhu                     "; 
      i=i+1; volcLat(i)=  3.670; volcLon(i)=125.500; volcElev(i)= 1320;  volcLoc(i)="Sangihe Is-Indonesia          "; 
             volcESP_Code(i)="S1"; volcID(i)="0607-04="; volcName(i)="Awu                            "; 
      i=i+1; volcLat(i)=  3.970; volcLon(i)=124.170; volcElev(i)=-5000;  volcLoc(i)="Sangihe Is-Indonesia          "; 
             volcESP_Code(i)="U0"; volcID(i)="0607-05="; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=  1.830; volcLon(i)=127.830; volcElev(i)=  318;  volcLoc(i)="Halmahera-Indonesia           "; 
             volcESP_Code(i)="M0"; volcID(i)="0608-001"; volcName(i)="Tarakan                        "; 
      i=i+1; volcLat(i)=  1.680; volcLon(i)=127.880; volcElev(i)= 1335;  volcLoc(i)="Halmahera-Indonesia           "; 
             volcESP_Code(i)="S1"; volcID(i)="0608-01="; volcName(i)="Dukono                         "; 
      i=i+1; volcLat(i)=  1.630; volcLon(i)=127.670; volcElev(i)= 1035;  volcLoc(i)="Halmahera-Indonesia           "; 
             volcESP_Code(i)="S0"; volcID(i)="0608-02-"; volcName(i)="Tobaru                         "; 
      i=i+1; volcLat(i)=  1.488; volcLon(i)=127.630; volcElev(i)= 1325;  volcLoc(i)="Halmahera-Indonesia           "; 
             volcESP_Code(i)="S1"; volcID(i)="0608-03="; volcName(i)="Ibu                            "; 
      i=i+1; volcLat(i)=  1.380; volcLon(i)=127.530; volcElev(i)= 1635;  volcLoc(i)="Halmahera-Indonesia           "; 
             volcESP_Code(i)="S1"; volcID(i)="0608-04="; volcName(i)="Gamkonora                      "; 
      i=i+1; volcLat(i)=  1.080; volcLon(i)=127.420; volcElev(i)= 1130;  volcLoc(i)="Halmahera-Indonesia           "; 
             volcESP_Code(i)="S0"; volcID(i)="0608-051"; volcName(i)="Jailolo                        "; 
      i=i+1; volcLat(i)=  0.900; volcLon(i)=127.320; volcElev(i)=  630;  volcLoc(i)="Halmahera-Indonesia           "; 
             volcESP_Code(i)="M0"; volcID(i)="0608-052"; volcName(i)="Hiri                           "; 
      i=i+1; volcLat(i)=  1.250; volcLon(i)=127.470; volcElev(i)=  979;  volcLoc(i)="Halmahera-Indonesia           "; 
             volcESP_Code(i)="S0"; volcID(i)="0608-05="; volcName(i)="Todoko-Ranu                    "; 
      i=i+1; volcLat(i)=  0.658; volcLon(i)=127.400; volcElev(i)= 1730;  volcLoc(i)="Halmahera-Indonesia           "; 
             volcESP_Code(i)="S0"; volcID(i)="0608-061"; volcName(i)="Tidore                         "; 
      i=i+1; volcLat(i)=  0.570; volcLon(i)=127.400; volcElev(i)=  308;  volcLoc(i)="Halmahera-Indonesia           "; 
             volcESP_Code(i)="S0"; volcID(i)="0608-062"; volcName(i)="Mare                           "; 
      i=i+1; volcLat(i)=  0.450; volcLon(i)=127.400; volcElev(i)=  950;  volcLoc(i)="Halmahera-Indonesia           "; 
             volcESP_Code(i)="S0"; volcID(i)="0608-063"; volcName(i)="Moti                           "; 
      i=i+1; volcLat(i)=  0.800; volcLon(i)=127.330; volcElev(i)= 1715;  volcLoc(i)="Halmahera-Indonesia           "; 
             volcESP_Code(i)="S1"; volcID(i)="0608-06="; volcName(i)="Gamalama                       "; 
      i=i+1; volcLat(i)=  0.070; volcLon(i)=127.420; volcElev(i)=  422;  volcLoc(i)="Halmahera-Indonesia           "; 
             volcESP_Code(i)="S0"; volcID(i)="0608-071"; volcName(i)="Tigalalu                       "; 
      i=i+1; volcLat(i)= -0.530; volcLon(i)=127.480; volcElev(i)= 1030;  volcLoc(i)="Halmahera-Indonesia           "; 
             volcESP_Code(i)="S0"; volcID(i)="0608-072"; volcName(i)="Amasing                        "; 
      i=i+1; volcLat(i)= -0.770; volcLon(i)=127.720; volcElev(i)=  900;  volcLoc(i)="Halmahera-Indonesia           "; 
             volcESP_Code(i)="S0"; volcID(i)="0608-073"; volcName(i)="Bibinoi                        "; 
      i=i+1; volcLat(i)=  0.320; volcLon(i)=127.400; volcElev(i)= 1357;  volcLoc(i)="Halmahera-Indonesia           "; 
             volcESP_Code(i)="S2"; volcID(i)="0608-07="; volcName(i)="Makian                         "; 
      i=i+1; volcLat(i)=  4.400; volcLon(i)=117.880; volcElev(i)=  531;  volcLoc(i)="Borneo                        "; 
             volcESP_Code(i)="M0"; volcID(i)="0610-01-"; volcName(i)="Bombalai                       "; 
      i=i+1; volcLat(i)=  6.013; volcLon(i)=121.057; volcElev(i)=  811;  volcLoc(i)="Sulu Is-Philippines           "; 
             volcESP_Code(i)="M0"; volcID(i)="0700-01="; volcName(i)="Jolo                           "; 
      i=i+1; volcLat(i)=  6.113; volcLon(i)=124.892; volcElev(i)= 1824;  volcLoc(i)="Mindanao-Philippines          "; 
             volcESP_Code(i)="S0"; volcID(i)="0701-011"; volcName(i)="Parker                         "; 
      i=i+1; volcLat(i)=  5.400; volcLon(i)=125.375; volcElev(i)=  862;  volcLoc(i)="Mindanao-Philippines          "; 
             volcESP_Code(i)="S0"; volcID(i)="0701-01="; volcName(i)="Balut                          "; 
      i=i+1; volcLat(i)=  6.370; volcLon(i)=125.070; volcElev(i)= 2286;  volcLoc(i)="Mindanao-Philippines          "; 
             volcESP_Code(i)="S0"; volcID(i)="0701-02="; volcName(i)="Matutum                        "; 
      i=i+1; volcLat(i)=  7.382; volcLon(i)=126.047; volcElev(i)= 1080;  volcLoc(i)="Mindanao-Philippines          "; 
             volcESP_Code(i)="S0"; volcID(i)="0701-031"; volcName(i)="Leonard Range                  "; 
      i=i+1; volcLat(i)=  6.989; volcLon(i)=125.269; volcElev(i)= 2938;  volcLoc(i)="Mindanao-Philippines          "; 
             volcESP_Code(i)="S0"; volcID(i)="0701-03="; volcName(i)="Apo                            "; 
      i=i+1; volcLat(i)=  7.647; volcLon(i)=124.320; volcElev(i)= 1940;  volcLoc(i)="Mindanao-Philippines          "; 
             volcESP_Code(i)="S0"; volcID(i)="0701-04="; volcName(i)="Makaturing                     "; 
      i=i+1; volcLat(i)=  7.650; volcLon(i)=124.450; volcElev(i)= 2338;  volcLoc(i)="Mindanao-Philippines          "; 
             volcESP_Code(i)="S0"; volcID(i)="0701-05="; volcName(i)="Latukan                        "; 
      i=i+1; volcLat(i)=  7.950; volcLon(i)=124.800; volcElev(i)= 2824;  volcLoc(i)="Mindanao-Philippines          "; 
             volcESP_Code(i)="S0"; volcID(i)="0701-061"; volcName(i)="Kalatungan                     "; 
      i=i+1; volcLat(i)=  7.700; volcLon(i)=124.500; volcElev(i)= 2815;  volcLoc(i)="Mindanao-Philippines          "; 
             volcESP_Code(i)="S1"; volcID(i)="0701-06="; volcName(i)="Ragang                         "; 
      i=i+1; volcLat(i)=  8.220; volcLon(i)=123.630; volcElev(i)= 2404;  volcLoc(i)="Mindanao-Philippines          "; 
             volcESP_Code(i)="S0"; volcID(i)="0701-071"; volcName(i)="Malindang                      "; 
      i=i+1; volcLat(i)=  8.770; volcLon(i)=124.980; volcElev(i)= 2450;  volcLoc(i)="Mindanao-Philippines          "; 
             volcESP_Code(i)="M0"; volcID(i)="0701-072"; volcName(i)="Balatukan                      "; 
      i=i+1; volcLat(i)=  7.877; volcLon(i)=125.068; volcElev(i)=  646;  volcLoc(i)="Mindanao-Philippines          "; 
             volcESP_Code(i)="S0"; volcID(i)="0701-07="; volcName(i)="Musuan                         "; 
      i=i+1; volcLat(i)=  9.203; volcLon(i)=124.673; volcElev(i)= 1552;  volcLoc(i)="Mindanao-Philippines          "; 
             volcESP_Code(i)="S0"; volcID(i)="0701-08="; volcName(i)="Camiguin                       "; 
      i=i+1; volcLat(i)=  9.593; volcLon(i)=125.520; volcElev(i)=  524;  volcLoc(i)="Mindanao-Philippines          "; 
             volcESP_Code(i)="S0"; volcID(i)="0701-09-"; volcName(i)="Paco                           "; 
      i=i+1; volcLat(i)=  9.250; volcLon(i)=123.170; volcElev(i)= 1862;  volcLoc(i)="Philippines-C                 "; 
             volcESP_Code(i)="S0"; volcID(i)="0702-01="; volcName(i)="Cuernos de Negros              "; 
      i=i+1; volcLat(i)= 10.412; volcLon(i)=123.132; volcElev(i)= 2435;  volcLoc(i)="Philippines-C                 "; 
             volcESP_Code(i)="S1"; volcID(i)="0702-02="; volcName(i)="Kanlaon                        "; 
      i=i+1; volcLat(i)= 10.650; volcLon(i)=123.250; volcElev(i)= 1885;  volcLoc(i)="Philippines-C                 "; 
             volcESP_Code(i)="S0"; volcID(i)="0702-03="; volcName(i)="Mandalagan                     "; 
      i=i+1; volcLat(i)= 10.770; volcLon(i)=123.230; volcElev(i)= 1510;  volcLoc(i)="Philippines-C                 "; 
             volcESP_Code(i)="S0"; volcID(i)="0702-04="; volcName(i)="Silay                          "; 
      i=i+1; volcLat(i)= 10.287; volcLon(i)=125.221; volcElev(i)=  945;  volcLoc(i)="Philippines-C                 "; 
             volcESP_Code(i)="S0"; volcID(i)="0702-05="; volcName(i)="Cabalían                       "; 
      i=i+1; volcLat(i)= 10.896; volcLon(i)=125.870; volcElev(i)=  860;  volcLoc(i)="Philippines-C                 "; 
             volcESP_Code(i)="S0"; volcID(i)="0702-07="; volcName(i)="Mahagnao                       "; 
      i=i+1; volcLat(i)= 11.523; volcLon(i)=124.535; volcElev(i)= 1301;  volcLoc(i)="Philippines-C                 "; 
             volcESP_Code(i)="S0"; volcID(i)="0702-08="; volcName(i)="Biliran                        "; 
      i=i+1; volcLat(i)= 12.770; volcLon(i)=124.050; volcElev(i)= 1565;  volcLoc(i)="Luzon-Philippines             "; 
             volcESP_Code(i)="S1"; volcID(i)="0703-01="; volcName(i)="Bulusan                        "; 
      i=i+1; volcLat(i)= 13.050; volcLon(i)=123.958; volcElev(i)= 1102;  volcLoc(i)="Luzon-Philippines             "; 
             volcESP_Code(i)="S0"; volcID(i)="0703-02="; volcName(i)="Pocdol Mountains               "; 
      i=i+1; volcLat(i)= 13.320; volcLon(i)=123.600; volcElev(i)= 1328;  volcLoc(i)="Luzon-Philippines             "; 
             volcESP_Code(i)="S0"; volcID(i)="0703-031"; volcName(i)="Masaraga                       "; 
      i=i+1; volcLat(i)= 13.257; volcLon(i)=123.685; volcElev(i)= 2462;  volcLoc(i)="Luzon-Philippines             "; 
             volcESP_Code(i)="S2"; volcID(i)="0703-03="; volcName(i)="Mayon                          "; 
      i=i+1; volcLat(i)= 13.457; volcLon(i)=123.457; volcElev(i)= 1196;  volcLoc(i)="Luzon-Philippines             "; 
             volcESP_Code(i)="S0"; volcID(i)="0703-041"; volcName(i)="Iriga                          "; 
      i=i+1; volcLat(i)= 13.658; volcLon(i)=123.380; volcElev(i)= 1966;  volcLoc(i)="Luzon-Philippines             "; 
             volcESP_Code(i)="S0"; volcID(i)="0703-042"; volcName(i)="Isarog                         "; 
      i=i+1; volcLat(i)= 13.240; volcLon(i)=122.018; volcElev(i)= 1157;  volcLoc(i)="Luzon-Philippines             "; 
             volcESP_Code(i)="S0"; volcID(i)="0703-044"; volcName(i)="Malindig                       "; 
      i=i+1; volcLat(i)= 14.070; volcLon(i)=121.480; volcElev(i)= 2158;  volcLoc(i)="Luzon-Philippines             "; 
             volcESP_Code(i)="S0"; volcID(i)="0703-05="; volcName(i)="Banahaw                        "; 
      i=i+1; volcLat(i)= 14.120; volcLon(i)=121.300; volcElev(i)= 1090;  volcLoc(i)="Luzon-Philippines             "; 
             volcESP_Code(i)="S0"; volcID(i)="0703-06="; volcName(i)="San Pablo Volc Field           "; 
      i=i+1; volcLat(i)= 14.002; volcLon(i)=120.993; volcElev(i)=  311;  volcLoc(i)="Luzon-Philippines             "; 
             volcESP_Code(i)="S2"; volcID(i)="0703-07="; volcName(i)="Taal                           "; 
      i=i+1; volcLat(i)= 14.520; volcLon(i)=120.470; volcElev(i)= 1388;  volcLoc(i)="Luzon-Philippines             "; 
             volcESP_Code(i)="S0"; volcID(i)="0703-081"; volcName(i)="Mariveles                      "; 
      i=i+1; volcLat(i)= 14.720; volcLon(i)=120.400; volcElev(i)= 1253;  volcLoc(i)="Luzon-Philippines             "; 
             volcESP_Code(i)="S0"; volcID(i)="0703-082"; volcName(i)="Natib                          "; 
      i=i+1; volcLat(i)= 15.130; volcLon(i)=120.350; volcElev(i)= 1486;  volcLoc(i)="Luzon-Philippines             "; 
             volcESP_Code(i)="S3"; volcID(i)="0703-083"; volcName(i)="Pinatubo                       "; 
      i=i+1; volcLat(i)= 15.200; volcLon(i)=120.742; volcElev(i)= 1026;  volcLoc(i)="Luzon-Philippines             "; 
             volcESP_Code(i)="S0"; volcID(i)="0703-084"; volcName(i)="Arayat                         "; 
      i=i+1; volcLat(i)= 15.828; volcLon(i)=120.805; volcElev(i)=  376;  volcLoc(i)="Luzon-Philippines             "; 
             volcESP_Code(i)="S0"; volcID(i)="0703-085"; volcName(i)="Amorong                        "; 
      i=i+1; volcLat(i)= 16.330; volcLon(i)=120.550; volcElev(i)= 2260;  volcLoc(i)="Luzon-Philippines             "; 
             volcESP_Code(i)="S0"; volcID(i)="0703-086"; volcName(i)="Santo Tomas                    "; 
      i=i+1; volcLat(i)= 17.147; volcLon(i)=120.980; volcElev(i)= 1865;  volcLoc(i)="Luzon-Philippines             "; 
             volcESP_Code(i)="S0"; volcID(i)="0703-087"; volcName(i)="Patoc                          "; 
      i=i+1; volcLat(i)= 17.320; volcLon(i)=121.100; volcElev(i)= 2329;  volcLoc(i)="Luzon-Philippines             "; 
             volcESP_Code(i)="S0"; volcID(i)="0703-088"; volcName(i)="Ambalatungan Group             "; 
      i=i+1; volcLat(i)= 14.420; volcLon(i)=121.270; volcElev(i)=  743;  volcLoc(i)="Luzon-Philippines             "; 
             volcESP_Code(i)="S0"; volcID(i)="0703-08="; volcName(i)="Laguna Caldera                 "; 
      i=i+1; volcLat(i)= 18.222; volcLon(i)=122.123; volcElev(i)= 1133;  volcLoc(i)="Luzon-Philippines             "; 
             volcESP_Code(i)="S0"; volcID(i)="0703-09="; volcName(i)="Cagua                          "; 
      i=i+1; volcLat(i)= 18.830; volcLon(i)=121.860; volcElev(i)=  712;  volcLoc(i)="Luzon-N of                    "; 
             volcESP_Code(i)="S0"; volcID(i)="0704-01="; volcName(i)="Camiguin de Babuyanes          "; 
      i=i+1; volcLat(i)= 19.077; volcLon(i)=122.202; volcElev(i)=  228;  volcLoc(i)="Luzon-N of                    "; 
             volcESP_Code(i)="S0"; volcID(i)="0704-02="; volcName(i)="Didicas                        "; 
      i=i+1; volcLat(i)= 19.523; volcLon(i)=121.940; volcElev(i)= 1080;  volcLoc(i)="Luzon-N of                    "; 
             volcESP_Code(i)="S1"; volcID(i)="0704-03="; volcName(i)="Babuyan Claro                  "; 
      i=i+1; volcLat(i)= 20.330; volcLon(i)=121.750; volcElev(i)=  -24;  volcLoc(i)="Luzon-N of                    "; 
             volcESP_Code(i)="S0"; volcID(i)="0704-05="; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 20.469; volcLon(i)=122.010; volcElev(i)= 1009;  volcLoc(i)="Luzon-N of                    "; 
             volcESP_Code(i)="S0"; volcID(i)="0704-06-"; volcName(i)="Iraya                          "; 
      i=i+1; volcLat(i)= 19.700; volcLon(i)=110.100; volcElev(i)=    0;  volcLoc(i)="SE Asia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0705-001"; volcName(i)="Hainan Dao                     "; 
      i=i+1; volcLat(i)= 20.780; volcLon(i)=110.170; volcElev(i)=  259;  volcLoc(i)="SE Asia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0705-01-"; volcName(i)="Leizhou Bandao                 "; 
      i=i+1; volcLat(i)= 15.380; volcLon(i)=109.120; volcElev(i)=  181;  volcLoc(i)="SE Asia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0705-02-"; volcName(i)="Cù-Lao Ré Group                "; 
      i=i+1; volcLat(i)= 14.930; volcLon(i)=108.000; volcElev(i)=  800;  volcLoc(i)="SE Asia                       "; 
             volcESP_Code(i)="M0"; volcID(i)="0705-03-"; volcName(i)="Toroeng Prong                  "; 
      i=i+1; volcLat(i)= 11.600; volcLon(i)=108.200; volcElev(i)= 1000;  volcLoc(i)="SE Asia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0705-04-"; volcName(i)="Haut Dong Nai                  "; 
      i=i+1; volcLat(i)= 10.800; volcLon(i)=107.200; volcElev(i)=  392;  volcLoc(i)="SE Asia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0705-05-"; volcName(i)="Bas Dong Nai                   "; 
      i=i+1; volcLat(i)= 10.158; volcLon(i)=109.014; volcElev(i)=  -20;  volcLoc(i)="SE Asia                       "; 
             volcESP_Code(i)="M0"; volcID(i)="0705-06-"; volcName(i)="Cendres, Ile des               "; 
      i=i+1; volcLat(i)=  9.830; volcLon(i)=109.050; volcElev(i)=    0;  volcLoc(i)="SE Asia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0705-07-"; volcName(i)="Veteran                        "; 
      i=i+1; volcLat(i)= 20.920; volcLon(i)= 95.250; volcElev(i)= 1518;  volcLoc(i)="SE Asia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0705-08-"; volcName(i)="Popa                           "; 
      i=i+1; volcLat(i)= 22.280; volcLon(i)= 95.100; volcElev(i)=  385;  volcLoc(i)="SE Asia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0705-09-"; volcName(i)="Lower Chindwin                 "; 
      i=i+1; volcLat(i)= 22.700; volcLon(i)= 95.980; volcElev(i)=  507;  volcLoc(i)="SE Asia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0705-10-"; volcName(i)="Singu Plateau                  "; 
      i=i+1; volcLat(i)= 25.230; volcLon(i)= 98.500; volcElev(i)= 2865;  volcLoc(i)="SE Asia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="0705-11-"; volcName(i)="Tengchong                      "; 
      i=i+1; volcLat(i)= 19.170; volcLon(i)=132.250; volcElev(i)=  -10;  volcLoc(i)="Taiwan-E of                   "; 
             volcESP_Code(i)="S0"; volcID(i)="0801-011"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 20.930; volcLon(i)=134.750; volcElev(i)=-6000;  volcLoc(i)="Taiwan-E of                   "; 
             volcESP_Code(i)="U0"; volcID(i)="0801-01="; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 21.830; volcLon(i)=121.180; volcElev(i)= -115;  volcLoc(i)="Taiwan-E of                   "; 
             volcESP_Code(i)="U0"; volcID(i)="0801-02="; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 24.850; volcLon(i)=121.920; volcElev(i)=  401;  volcLoc(i)="Taiwan-E of                   "; 
             volcESP_Code(i)="S0"; volcID(i)="0801-031"; volcName(i)="Kueishantao                    "; 
      i=i+1; volcLat(i)= 24.000; volcLon(i)=121.830; volcElev(i)=    0;  volcLoc(i)="Taiwan-E of                   "; 
             volcESP_Code(i)="S0"; volcID(i)="0801-03="; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 25.400; volcLon(i)=122.200; volcElev(i)= -100;  volcLoc(i)="Taiwan-N of                   "; 
             volcESP_Code(i)="U0"; volcID(i)="0801-04="; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 26.180; volcLon(i)=122.458; volcElev(i)= -418;  volcLoc(i)="Taiwan-N of                   "; 
             volcESP_Code(i)="U0"; volcID(i)="0801-05="; volcName(i)="Zengyu                         "; 
      i=i+1; volcLat(i)= 24.558; volcLon(i)=124.000; volcElev(i)= -200;  volcLoc(i)="Ryukyu Is                     "; 
             volcESP_Code(i)="U0"; volcID(i)="0802-01="; volcName(i)="Iriomote-jima                  "; 
      i=i+1; volcLat(i)= 28.797; volcLon(i)=128.997; volcElev(i)=  495;  volcLoc(i)="Ryukyu Is                     "; 
             volcESP_Code(i)="S0"; volcID(i)="0802-021"; volcName(i)="Yokoate-jima                   "; 
      i=i+1; volcLat(i)= 29.461; volcLon(i)=129.597; volcElev(i)=  584;  volcLoc(i)="Ryukyu Is                     "; 
             volcESP_Code(i)="S0"; volcID(i)="0802-022"; volcName(i)="Akuseki-jima                   "; 
      i=i+1; volcLat(i)= 27.877; volcLon(i)=128.224; volcElev(i)=  212;  volcLoc(i)="Ryukyu Is                     "; 
             volcESP_Code(i)="S1"; volcID(i)="0802-02="; volcName(i)="Iwo-Tori-shima                 "; 
      i=i+1; volcLat(i)= 29.635; volcLon(i)=129.716; volcElev(i)=  799;  volcLoc(i)="Ryukyu Is                     "; 
             volcESP_Code(i)="S1"; volcID(i)="0802-03="; volcName(i)="Suwanose-jima                  "; 
      i=i+1; volcLat(i)= 29.879; volcLon(i)=129.625; volcElev(i)=  301;  volcLoc(i)="Ryukyu Is                     "; 
             volcESP_Code(i)="S0"; volcID(i)="0802-041"; volcName(i)="Kogaja-jima                    "; 
      i=i+1; volcLat(i)= 29.964; volcLon(i)=129.927; volcElev(i)=  628;  volcLoc(i)="Ryukyu Is                     "; 
             volcESP_Code(i)="S0"; volcID(i)="0802-043"; volcName(i)="Kuchino-shima                  "; 
      i=i+1; volcLat(i)= 29.856; volcLon(i)=129.859; volcElev(i)=  979;  volcLoc(i)="Ryukyu Is                     "; 
             volcESP_Code(i)="S0"; volcID(i)="0802-04="; volcName(i)="Nakano-shima                   "; 
      i=i+1; volcLat(i)= 30.440; volcLon(i)=130.219; volcElev(i)=  657;  volcLoc(i)="Ryukyu Is                     "; 
             volcESP_Code(i)="S2"; volcID(i)="0802-05="; volcName(i)="Kuchinoerabu-jima              "; 
      i=i+1; volcLat(i)= 30.789; volcLon(i)=130.308; volcElev(i)=  704;  volcLoc(i)="Ryukyu Is                     "; 
             volcESP_Code(i)="S1"; volcID(i)="0802-06="; volcName(i)="Kikai                          "; 
      i=i+1; volcLat(i)= 31.220; volcLon(i)=130.570; volcElev(i)=  922;  volcLoc(i)="Kyushu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0802-07="; volcName(i)="Ibusuki Volc Field             "; 
      i=i+1; volcLat(i)= 31.768; volcLon(i)=130.594; volcElev(i)=   15;  volcLoc(i)="Kyushu-Japan                  "; 
             volcESP_Code(i)="M0"; volcID(i)="0802-081"; volcName(i)="Sumiyoshi-ike                  "; 
      i=i+1; volcLat(i)= 31.585; volcLon(i)=130.657; volcElev(i)= 1117;  volcLoc(i)="Kyushu-Japan                  "; 
             volcESP_Code(i)="S1"; volcID(i)="0802-08="; volcName(i)="Sakura-jima                    "; 
      i=i+1; volcLat(i)= 32.653; volcLon(i)=128.851; volcElev(i)=  317;  volcLoc(i)="Kyushu-Japan                  "; 
             volcESP_Code(i)="M0"; volcID(i)="0802-091"; volcName(i)="Fukue-jima                     "; 
      i=i+1; volcLat(i)= 31.931; volcLon(i)=130.864; volcElev(i)= 1700;  volcLoc(i)="Kyushu-Japan                  "; 
             volcESP_Code(i)="S1"; volcID(i)="0802-09="; volcName(i)="Kirishima                      "; 
      i=i+1; volcLat(i)= 32.757; volcLon(i)=130.294; volcElev(i)= 1500;  volcLoc(i)="Kyushu-Japan                  "; 
             volcESP_Code(i)="S1"; volcID(i)="0802-10="; volcName(i)="Unzen                          "; 
      i=i+1; volcLat(i)= 32.881; volcLon(i)=131.106; volcElev(i)= 1592;  volcLoc(i)="Kyushu-Japan                  "; 
             volcESP_Code(i)="M1"; volcID(i)="0802-11="; volcName(i)="Aso                            "; 
      i=i+1; volcLat(i)= 33.083; volcLon(i)=131.251; volcElev(i)= 1791;  volcLoc(i)="Kyushu-Japan                  "; 
             volcESP_Code(i)="M1"; volcID(i)="0802-12="; volcName(i)="Kuju                           "; 
      i=i+1; volcLat(i)= 33.280; volcLon(i)=131.432; volcElev(i)= 1584;  volcLoc(i)="Kyushu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0802-13="; volcName(i)="Tsurumi                        "; 
      i=i+1; volcLat(i)= 34.500; volcLon(i)=131.600; volcElev(i)=  641;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="M0"; volcID(i)="0803-001"; volcName(i)="Abu                            "; 
      i=i+1; volcLat(i)= 35.130; volcLon(i)=132.620; volcElev(i)= 1126;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0803-002"; volcName(i)="Sanbe                          "; 
      i=i+1; volcLat(i)= 36.176; volcLon(i)=133.334; volcElev(i)=  151;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="M0"; volcID(i)="0803-003"; volcName(i)="Oki-Dogo                       "; 
      i=i+1; volcLat(i)= 34.900; volcLon(i)=139.098; volcElev(i)= 1406;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="M0"; volcID(i)="0803-01="; volcName(i)="Izu-Tobu                       "; 
      i=i+1; volcLat(i)= 35.230; volcLon(i)=139.024; volcElev(i)= 1438;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0803-02="; volcName(i)="Hakone                         "; 
      i=i+1; volcLat(i)= 36.100; volcLon(i)=138.300; volcElev(i)= 2530;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0803-031"; volcName(i)="Kita Yatsuga-take              "; 
      i=i+1; volcLat(i)= 35.358; volcLon(i)=138.731; volcElev(i)= 3776;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S2"; volcID(i)="0803-03="; volcName(i)="Fuji                           "; 
      i=i+1; volcLat(i)= 35.890; volcLon(i)=137.480; volcElev(i)= 3063;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0803-04="; volcName(i)="On-take                        "; 
      i=i+1; volcLat(i)= 36.152; volcLon(i)=136.774; volcElev(i)= 2702;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S2"; volcID(i)="0803-05="; volcName(i)="Haku-san                       "; 
      i=i+1; volcLat(i)= 36.103; volcLon(i)=137.557; volcElev(i)= 3026;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0803-06="; volcName(i)="Norikura                       "; 
      i=i+1; volcLat(i)= 36.408; volcLon(i)=137.594; volcElev(i)= 2924;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="M0"; volcID(i)="0803-071"; volcName(i)="Washiba-Kumonotaira            "; 
      i=i+1; volcLat(i)= 36.224; volcLon(i)=137.590; volcElev(i)= 2455;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0803-07="; volcName(i)="Yake-dake                      "; 
      i=i+1; volcLat(i)= 36.568; volcLon(i)=137.593; volcElev(i)= 2621;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0803-08="; volcName(i)="Tate-yama                      "; 
      i=i+1; volcLat(i)= 36.918; volcLon(i)=138.039; volcElev(i)= 2400;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S1"; volcID(i)="0803-09="; volcName(i)="Niigata-Yake-yama              "; 
      i=i+1; volcLat(i)= 36.888; volcLon(i)=138.120; volcElev(i)= 2446;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0803-10="; volcName(i)="Myoko                          "; 
      i=i+1; volcLat(i)= 36.403; volcLon(i)=138.526; volcElev(i)= 2568;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S1"; volcID(i)="0803-11="; volcName(i)="Asama                          "; 
      i=i+1; volcLat(i)= 36.688; volcLon(i)=138.519; volcElev(i)= 2041;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="M0"; volcID(i)="0803-121"; volcName(i)="Shiga                          "; 
      i=i+1; volcLat(i)= 36.474; volcLon(i)=138.881; volcElev(i)= 1449;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0803-122"; volcName(i)="Haruna                         "; 
      i=i+1; volcLat(i)= 36.620; volcLon(i)=138.535; volcElev(i)= 2171;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S1"; volcID(i)="0803-12="; volcName(i)="Kusatsu-Shirane                "; 
      i=i+1; volcLat(i)= 36.952; volcLon(i)=139.289; volcElev(i)= 2356;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0803-131"; volcName(i)="Hiuchi                         "; 
      i=i+1; volcLat(i)= 36.557; volcLon(i)=139.196; volcElev(i)= 1828;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0803-13="; volcName(i)="Akagi                          "; 
      i=i+1; volcLat(i)= 36.762; volcLon(i)=139.494; volcElev(i)= 2486;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0803-141"; volcName(i)="Nantai                         "; 
      i=i+1; volcLat(i)= 36.792; volcLon(i)=139.510; volcElev(i)= 2367;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0803-142"; volcName(i)="Omanago Group                  "; 
      i=i+1; volcLat(i)= 36.897; volcLon(i)=139.780; volcElev(i)= 1795;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0803-143"; volcName(i)="Takahara                       "; 
      i=i+1; volcLat(i)= 36.796; volcLon(i)=139.379; volcElev(i)= 2578;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S1"; volcID(i)="0803-14="; volcName(i)="Nikko-Shirane                  "; 
      i=i+1; volcLat(i)= 37.450; volcLon(i)=139.579; volcElev(i)= 1100;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="M0"; volcID(i)="0803-151"; volcName(i)="Numazawa                       "; 
      i=i+1; volcLat(i)= 37.122; volcLon(i)=139.966; volcElev(i)= 1915;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S1"; volcID(i)="0803-15="; volcName(i)="Nasu                           "; 
      i=i+1; volcLat(i)= 37.598; volcLon(i)=140.076; volcElev(i)= 1819;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S1"; volcID(i)="0803-16="; volcName(i)="Bandai                         "; 
      i=i+1; volcLat(i)= 37.644; volcLon(i)=140.286; volcElev(i)= 1718;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S1"; volcID(i)="0803-17="; volcName(i)="Adatara                        "; 
      i=i+1; volcLat(i)= 37.732; volcLon(i)=140.248; volcElev(i)= 2035;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S1"; volcID(i)="0803-18="; volcName(i)="Azuma                          "; 
      i=i+1; volcLat(i)= 38.606; volcLon(i)=140.178; volcElev(i)=  516;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0803-191"; volcName(i)="Hijiori                        "; 
      i=i+1; volcLat(i)= 38.141; volcLon(i)=140.443; volcElev(i)= 1841;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S1"; volcID(i)="0803-19="; volcName(i)="Zao                            "; 
      i=i+1; volcLat(i)= 38.733; volcLon(i)=140.732; volcElev(i)=  470;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0803-20="; volcName(i)="Narugo                         "; 
      i=i+1; volcLat(i)= 38.958; volcLon(i)=140.792; volcElev(i)= 1628;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S1"; volcID(i)="0803-21="; volcName(i)="Kurikoma                       "; 
      i=i+1; volcLat(i)= 39.096; volcLon(i)=140.052; volcElev(i)= 2233;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S1"; volcID(i)="0803-22="; volcName(i)="Chokai                         "; 
      i=i+1; volcLat(i)= 39.758; volcLon(i)=140.803; volcElev(i)= 1637;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S1"; volcID(i)="0803-23="; volcName(i)="Akita-Komaga-take              "; 
      i=i+1; volcLat(i)= 39.850; volcLon(i)=141.004; volcElev(i)= 2041;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0803-24="; volcName(i)="Iwate                          "; 
      i=i+1; volcLat(i)= 39.955; volcLon(i)=140.857; volcElev(i)= 1614;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0803-25="; volcName(i)="Hachimantai                    "; 
      i=i+1; volcLat(i)= 39.950; volcLon(i)=139.730; volcElev(i)=  291;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="M0"; volcID(i)="0803-262"; volcName(i)="Megata                         "; 
      i=i+1; volcLat(i)= 39.961; volcLon(i)=140.761; volcElev(i)= 1366;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0803-26="; volcName(i)="Akita-Yake-yama                "; 
      i=i+1; volcLat(i)= 40.470; volcLon(i)=140.920; volcElev(i)= 1159;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0803-271"; volcName(i)="Towada                         "; 
      i=i+1; volcLat(i)= 40.653; volcLon(i)=140.307; volcElev(i)= 1625;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S1"; volcID(i)="0803-27="; volcName(i)="Iwaki                          "; 
      i=i+1; volcLat(i)= 40.656; volcLon(i)=140.881; volcElev(i)= 1585;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0803-28="; volcName(i)="Hakkoda Group                  "; 
      i=i+1; volcLat(i)= 41.276; volcLon(i)=141.124; volcElev(i)=  879;  volcLoc(i)="Honshu-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0803-29="; volcName(i)="Osore-yama                     "; 
      i=i+1; volcLat(i)= 34.517; volcLon(i)=139.283; volcElev(i)=  508;  volcLoc(i)="Izu Is-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0804-011"; volcName(i)="To-shima                       "; 
      i=i+1; volcLat(i)= 34.721; volcLon(i)=139.398; volcElev(i)=  764;  volcLoc(i)="Izu Is-Japan                  "; 
             volcESP_Code(i)="S1"; volcID(i)="0804-01="; volcName(i)="Oshima                         "; 
      i=i+1; volcLat(i)= 34.393; volcLon(i)=139.273; volcElev(i)=  432;  volcLoc(i)="Izu Is-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0804-02="; volcName(i)="Nii-jima                       "; 
      i=i+1; volcLat(i)= 34.216; volcLon(i)=139.156; volcElev(i)=  572;  volcLoc(i)="Izu Is-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0804-03="; volcName(i)="Kozu-shima                     "; 
      i=i+1; volcLat(i)= 33.871; volcLon(i)=139.605; volcElev(i)=  851;  volcLoc(i)="Izu Is-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0804-041"; volcName(i)="Mikura-jima                    "; 
      i=i+1; volcLat(i)= 33.400; volcLon(i)=139.680; volcElev(i)= -107;  volcLoc(i)="Izu Is-Japan                  "; 
             volcESP_Code(i)="U0"; volcID(i)="0804-042"; volcName(i)="Kurose Hole                    "; 
      i=i+1; volcLat(i)= 34.079; volcLon(i)=139.529; volcElev(i)=  815;  volcLoc(i)="Izu Is-Japan                  "; 
             volcESP_Code(i)="M2"; volcID(i)="0804-04="; volcName(i)="Miyake-jima                    "; 
      i=i+1; volcLat(i)= 33.130; volcLon(i)=139.769; volcElev(i)=  854;  volcLoc(i)="Izu Is-Japan                  "; 
             volcESP_Code(i)="M0"; volcID(i)="0804-05="; volcName(i)="Hachijo-jima                   "; 
      i=i+1; volcLat(i)= 32.100; volcLon(i)=139.850; volcElev(i)=  360;  volcLoc(i)="Izu Is-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0804-061"; volcName(i)="Myojin Knoll                   "; 
      i=i+1; volcLat(i)= 32.454; volcLon(i)=139.762; volcElev(i)=  423;  volcLoc(i)="Izu Is-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0804-06="; volcName(i)="Aoga-shima                     "; 
      i=i+1; volcLat(i)= 31.880; volcLon(i)=139.920; volcElev(i)=   11;  volcLoc(i)="Izu Is-Japan                  "; 
             volcESP_Code(i)="S1"; volcID(i)="0804-07="; volcName(i)="Bayonnaise Rocks               "; 
      i=i+1; volcLat(i)= 31.436; volcLon(i)=140.054; volcElev(i)=  136;  volcLoc(i)="Izu Is-Japan                  "; 
             volcESP_Code(i)="S0"; volcID(i)="0804-08="; volcName(i)="Smith Rock                     "; 
      i=i+1; volcLat(i)= 29.789; volcLon(i)=140.345; volcElev(i)=   99;  volcLoc(i)="Izu Is-Japan                  "; 
             volcESP_Code(i)="M0"; volcID(i)="0804-091"; volcName(i)="Sofugan                        "; 
      i=i+1; volcLat(i)= 28.600; volcLon(i)=140.630; volcElev(i)=-1418;  volcLoc(i)="Izu Is-Japan                  "; 
             volcESP_Code(i)="U0"; volcID(i)="0804-093"; volcName(i)="Suiyo Seamount                 "; 
      i=i+1; volcLat(i)= 28.320; volcLon(i)=140.570; volcElev(i)= -920;  volcLoc(i)="Izu Is-Japan                  "; 
             volcESP_Code(i)="U0"; volcID(i)="0804-094"; volcName(i)="Mokuyo Seamount                "; 
      i=i+1; volcLat(i)= 27.680; volcLon(i)=140.800; volcElev(i)= -860;  volcLoc(i)="Izu Is-Japan                  "; 
             volcESP_Code(i)="U0"; volcID(i)="0804-095"; volcName(i)="Doyo Seamount                  "; 
      i=i+1; volcLat(i)= 27.274; volcLon(i)=140.882; volcElev(i)=   38;  volcLoc(i)="Volcano Is-Japan              "; 
             volcESP_Code(i)="S0"; volcID(i)="0804-096"; volcName(i)="Nishino-shima                  "; 
      i=i+1; volcLat(i)= 26.670; volcLon(i)=141.000; volcElev(i)= -162;  volcLoc(i)="Volcano Is-Japan              "; 
             volcESP_Code(i)="U0"; volcID(i)="0804-097"; volcName(i)="Kaikata Seamount               "; 
      i=i+1; volcLat(i)= 30.480; volcLon(i)=140.306; volcElev(i)=  394;  volcLoc(i)="Izu Is-Japan                  "; 
             volcESP_Code(i)="M0"; volcID(i)="0804-09="; volcName(i)="Tori-shima                     "; 
      i=i+1; volcLat(i)= 26.130; volcLon(i)=144.480; volcElev(i)=-3200;  volcLoc(i)="Volcano Is-Japan              "; 
             volcESP_Code(i)="U0"; volcID(i)="0804-101"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 26.122; volcLon(i)=141.102; volcElev(i)= -103;  volcLoc(i)="Volcano Is-Japan              "; 
             volcESP_Code(i)="U0"; volcID(i)="0804-10="; volcName(i)="Kaitoku Seamount               "; 
      i=i+1; volcLat(i)= 25.424; volcLon(i)=141.284; volcElev(i)=  792;  volcLoc(i)="Volcano Is-Japan              "; 
             volcESP_Code(i)="M0"; volcID(i)="0804-11="; volcName(i)="Kita-Iwo-jima                  "; 
      i=i+1; volcLat(i)= 24.414; volcLon(i)=141.419; volcElev(i)=  -73;  volcLoc(i)="Volcano Is-Japan              "; 
             volcESP_Code(i)="U0"; volcID(i)="0804-121"; volcName(i)="Kita-Fukutokutai               "; 
      i=i+1; volcLat(i)= 24.754; volcLon(i)=141.290; volcElev(i)=  161;  volcLoc(i)="Volcano Is-Japan              "; 
             volcESP_Code(i)="S1"; volcID(i)="0804-12="; volcName(i)="Ioto                           "; 
      i=i+1; volcLat(i)= 23.497; volcLon(i)=141.940; volcElev(i)=  -30;  volcLoc(i)="Volcano Is-Japan              "; 
             volcESP_Code(i)="M0"; volcID(i)="0804-131"; volcName(i)="Minami-Hiyoshi                 "; 
      i=i+1; volcLat(i)= 23.075; volcLon(i)=142.308; volcElev(i)= -391;  volcLoc(i)="Volcano Is-Japan              "; 
             volcESP_Code(i)="U0"; volcID(i)="0804-132"; volcName(i)="Nikko                          "; 
      i=i+1; volcLat(i)= 21.930; volcLon(i)=143.470; volcElev(i)= -217;  volcLoc(i)="Volcano Is-Japan              "; 
             volcESP_Code(i)="U0"; volcID(i)="0804-133"; volcName(i)="Fukujin                        "; 
      i=i+1; volcLat(i)= 21.765; volcLon(i)=143.710; volcElev(i)= -598;  volcLoc(i)="Volcano Is-Japan              "; 
             volcESP_Code(i)="U0"; volcID(i)="0804-134"; volcName(i)="Kasuga                         "; 
      i=i+1; volcLat(i)= 21.600; volcLon(i)=143.637; volcElev(i)= -274;  volcLoc(i)="Volcano Is-Japan              "; 
             volcESP_Code(i)="U0"; volcID(i)="0804-135"; volcName(i)="Minami Kasuga                  "; 
      i=i+1; volcLat(i)= 21.485; volcLon(i)=144.043; volcElev(i)=-1535;  volcLoc(i)="Volcano Is-Japan              "; 
             volcESP_Code(i)="U0"; volcID(i)="0804-136"; volcName(i)="NW Eifuku                      "; 
      i=i+1; volcLat(i)= 21.324; volcLon(i)=144.194; volcElev(i)= -323;  volcLoc(i)="Volcano Is-Japan              "; 
             volcESP_Code(i)="U0"; volcID(i)="0804-137"; volcName(i)="Daikoku                        "; 
      i=i+1; volcLat(i)= 21.000; volcLon(i)=142.900; volcElev(i)=    0;  volcLoc(i)="Mariana Is-C Pacific          "; 
             volcESP_Code(i)="U0"; volcID(i)="0804-138"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 20.300; volcLon(i)=143.200; volcElev(i)=    0;  volcLoc(i)="Mariana Is-C Pacific          "; 
             volcESP_Code(i)="U0"; volcID(i)="0804-139"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 24.280; volcLon(i)=141.485; volcElev(i)=  -14;  volcLoc(i)="Volcano Is-Japan              "; 
             volcESP_Code(i)="S1"; volcID(i)="0804-13="; volcName(i)="Fukutoku-Okanoba               "; 
      i=i+1; volcLat(i)= 20.420; volcLon(i)=145.030; volcElev(i)= -137;  volcLoc(i)="Mariana Is-C Pacific          "; 
             volcESP_Code(i)="U0"; volcID(i)="0804-141"; volcName(i)="Ahyi                           "; 
      i=i+1; volcLat(i)= 20.130; volcLon(i)=145.100; volcElev(i)=   -8;  volcLoc(i)="Mariana Is-C Pacific          "; 
             volcESP_Code(i)="S0"; volcID(i)="0804-142"; volcName(i)="Supply Reef                    "; 
      i=i+1; volcLat(i)= 20.020; volcLon(i)=145.220; volcElev(i)=  227;  volcLoc(i)="Mariana Is-C Pacific          "; 
             volcESP_Code(i)="M0"; volcID(i)="0804-143"; volcName(i)="Maug Islands                   "; 
      i=i+1; volcLat(i)= 20.538; volcLon(i)=144.896; volcElev(i)=  360;  volcLoc(i)="Mariana Is-C Pacific          "; 
             volcESP_Code(i)="S1"; volcID(i)="0804-14="; volcName(i)="Farallon de Pajaros            "; 
      i=i+1; volcLat(i)= 19.671; volcLon(i)=145.406; volcElev(i)=  857;  volcLoc(i)="Mariana Is-C Pacific          "; 
             volcESP_Code(i)="S0"; volcID(i)="0804-15="; volcName(i)="Asuncion                       "; 
      i=i+1; volcLat(i)= 18.770; volcLon(i)=145.670; volcElev(i)=  965;  volcLoc(i)="Mariana Is-C Pacific          "; 
             volcESP_Code(i)="M0"; volcID(i)="0804-16="; volcName(i)="Agrigan                        "; 
      i=i+1; volcLat(i)= 18.130; volcLon(i)=145.800; volcElev(i)=  570;  volcLoc(i)="Mariana Is-C Pacific          "; 
             volcESP_Code(i)="M2"; volcID(i)="0804-17="; volcName(i)="Pagan                          "; 
      i=i+1; volcLat(i)= 17.600; volcLon(i)=145.830; volcElev(i)=  744;  volcLoc(i)="Mariana Is-C Pacific          "; 
             volcESP_Code(i)="S0"; volcID(i)="0804-18="; volcName(i)="Alamagan                       "; 
      i=i+1; volcLat(i)= 16.880; volcLon(i)=145.850; volcElev(i)=    0;  volcLoc(i)="Mariana Is-C Pacific          "; 
             volcESP_Code(i)="S0"; volcID(i)="0804-191"; volcName(i)="Zealandia Bank                 "; 
      i=i+1; volcLat(i)= 16.708; volcLon(i)=145.780; volcElev(i)=  538;  volcLoc(i)="Mariana Is-C Pacific          "; 
             volcESP_Code(i)="S0"; volcID(i)="0804-192"; volcName(i)="Sarigan                        "; 
      i=i+1; volcLat(i)= 17.307; volcLon(i)=145.845; volcElev(i)=  287;  volcLoc(i)="Mariana Is-C Pacific          "; 
             volcESP_Code(i)="S0"; volcID(i)="0804-19="; volcName(i)="Guguan                         "; 
      i=i+1; volcLat(i)= 15.930; volcLon(i)=145.670; volcElev(i)= -127;  volcLoc(i)="Mariana Is-C Pacific          "; 
             volcESP_Code(i)="U0"; volcID(i)="0804-201"; volcName(i)="East Diamante                  "; 
      i=i+1; volcLat(i)= 15.620; volcLon(i)=145.570; volcElev(i)= -230;  volcLoc(i)="Mariana Is-C Pacific          "; 
             volcESP_Code(i)="U0"; volcID(i)="0804-202"; volcName(i)="Ruby                           "; 
      i=i+1; volcLat(i)= 16.350; volcLon(i)=145.670; volcElev(i)=  790;  volcLoc(i)="Mariana Is-C Pacific          "; 
             volcESP_Code(i)="S2"; volcID(i)="0804-20="; volcName(i)="Anatahan                       "; 
      i=i+1; volcLat(i)= 14.601; volcLon(i)=144.775; volcElev(i)= -517;  volcLoc(i)="Mariana Is-C Pacific          "; 
             volcESP_Code(i)="U0"; volcID(i)="0804-211"; volcName(i)="NW Rota-1                      "; 
      i=i+1; volcLat(i)= 15.000; volcLon(i)=145.250; volcElev(i)=  -43;  volcLoc(i)="Mariana Is-C Pacific          "; 
             volcESP_Code(i)="M0"; volcID(i)="0804-21="; volcName(i)="Esmeralda Bank                 "; 
      i=i+1; volcLat(i)= 13.400; volcLon(i)=143.920; volcElev(i)=    0;  volcLoc(i)="Mariana Is-C Pacific          "; 
             volcESP_Code(i)="M0"; volcID(i)="0804-22-"; volcName(i)="Forecast Seamount              "; 
      i=i+1; volcLat(i)= 13.250; volcLon(i)=144.020; volcElev(i)=-1230;  volcLoc(i)="Mariana Is-C Pacific          "; 
             volcESP_Code(i)="M0"; volcID(i)="0804-23-"; volcName(i)="Seamount X                     "; 
      i=i+1; volcLat(i)= 41.802; volcLon(i)=141.170; volcElev(i)=  618;  volcLoc(i)="Hokkaido-Japan                "; 
             volcESP_Code(i)="S0"; volcID(i)="0805-011"; volcName(i)="E-san                          "; 
      i=i+1; volcLat(i)= 41.507; volcLon(i)=139.371; volcElev(i)=  737;  volcLoc(i)="Hokkaido-Japan                "; 
             volcESP_Code(i)="M0"; volcID(i)="0805-01="; volcName(i)="Oshima-Oshima                  "; 
      i=i+1; volcLat(i)= 42.061; volcLon(i)=140.681; volcElev(i)= 1131;  volcLoc(i)="Hokkaido-Japan                "; 
             volcESP_Code(i)="S2"; volcID(i)="0805-02="; volcName(i)="Komaga-take                    "; 
      i=i+1; volcLat(i)= 42.880; volcLon(i)=140.630; volcElev(i)= 1154;  volcLoc(i)="Hokkaido-Japan                "; 
             volcESP_Code(i)="S0"; volcID(i)="0805-031"; volcName(i)="Niseko                         "; 
      i=i+1; volcLat(i)= 42.830; volcLon(i)=140.815; volcElev(i)= 1898;  volcLoc(i)="Hokkaido-Japan                "; 
             volcESP_Code(i)="S0"; volcID(i)="0805-032"; volcName(i)="Yotei                          "; 
      i=i+1; volcLat(i)= 42.489; volcLon(i)=141.163; volcElev(i)=  581;  volcLoc(i)="Hokkaido-Japan                "; 
             volcESP_Code(i)="S0"; volcID(i)="0805-034"; volcName(i)="Kuttara                        "; 
      i=i+1; volcLat(i)= 42.541; volcLon(i)=140.843; volcElev(i)=  737;  volcLoc(i)="Hokkaido-Japan                "; 
             volcESP_Code(i)="S0"; volcID(i)="0805-03="; volcName(i)="Usu                            "; 
      i=i+1; volcLat(i)= 45.180; volcLon(i)=141.250; volcElev(i)= 1721;  volcLoc(i)="Hokkaido-Japan                "; 
             volcESP_Code(i)="S0"; volcID(i)="0805-041"; volcName(i)="Rishiri                        "; 
      i=i+1; volcLat(i)= 42.688; volcLon(i)=141.380; volcElev(i)= 1320;  volcLoc(i)="Hokkaido-Japan                "; 
             volcESP_Code(i)="S0"; volcID(i)="0805-04="; volcName(i)="Shikotsu                       "; 
      i=i+1; volcLat(i)= 43.416; volcLon(i)=142.690; volcElev(i)= 2077;  volcLoc(i)="Hokkaido-Japan                "; 
             volcESP_Code(i)="S0"; volcID(i)="0805-05="; volcName(i)="Tokachi                        "; 
      i=i+1; volcLat(i)= 43.453; volcLon(i)=143.036; volcElev(i)= 2013;  volcLoc(i)="Hokkaido-Japan                "; 
             volcESP_Code(i)="S0"; volcID(i)="0805-061"; volcName(i)="Nipesotsu-Maruyama             "; 
      i=i+1; volcLat(i)= 43.312; volcLon(i)=143.096; volcElev(i)= 1401;  volcLoc(i)="Hokkaido-Japan                "; 
             volcESP_Code(i)="S0"; volcID(i)="0805-062"; volcName(i)="Shikaribetsu Group             "; 
      i=i+1; volcLat(i)= 43.661; volcLon(i)=142.858; volcElev(i)= 2290;  volcLoc(i)="Hokkaido-Japan                "; 
             volcESP_Code(i)="S0"; volcID(i)="0805-06="; volcName(i)="Daisetsu                       "; 
      i=i+1; volcLat(i)= 43.384; volcLon(i)=144.013; volcElev(i)= 1499;  volcLoc(i)="Hokkaido-Japan                "; 
             volcESP_Code(i)="S1"; volcID(i)="0805-07="; volcName(i)="Akan                           "; 
      i=i+1; volcLat(i)= 43.570; volcLon(i)=144.565; volcElev(i)=  855;  volcLoc(i)="Hokkaido-Japan                "; 
             volcESP_Code(i)="S0"; volcID(i)="0805-081"; volcName(i)="Mashu                          "; 
      i=i+1; volcLat(i)= 44.073; volcLon(i)=145.126; volcElev(i)= 1660;  volcLoc(i)="Hokkaido-Japan                "; 
             volcESP_Code(i)="S0"; volcID(i)="0805-082"; volcName(i)="Rausu                          "; 
      i=i+1; volcLat(i)= 43.608; volcLon(i)=144.443; volcElev(i)=  999;  volcLoc(i)="Hokkaido-Japan                "; 
             volcESP_Code(i)="S0"; volcID(i)="0805-08="; volcName(i)="Kutcharo                       "; 
      i=i+1; volcLat(i)= 44.131; volcLon(i)=145.165; volcElev(i)= 1563;  volcLoc(i)="Hokkaido-Japan                "; 
             volcESP_Code(i)="S0"; volcID(i)="0805-09="; volcName(i)="Shiretoko-Iwo-zan              "; 
      i=i+1; volcLat(i)= 43.841; volcLon(i)=145.509; volcElev(i)=  543;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-01="; volcName(i)="Golovnin                       "; 
      i=i+1; volcLat(i)= 44.420; volcLon(i)=146.135; volcElev(i)= 1189;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-021"; volcName(i)="Smirnov                        "; 
      i=i+1; volcLat(i)= 43.976; volcLon(i)=145.736; volcElev(i)=  888;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-02="; volcName(i)="Mendeleev                      "; 
      i=i+1; volcLat(i)= 44.351; volcLon(i)=146.256; volcElev(i)= 1819;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0900-03="; volcName(i)="Tiatia                         "; 
      i=i+1; volcLat(i)= 44.608; volcLon(i)=146.994; volcElev(i)=  528;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-041"; volcName(i)="Lvinaya Past                   "; 
      i=i+1; volcLat(i)= 44.459; volcLon(i)=146.936; volcElev(i)= 1221;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-04="; volcName(i)="Berutarube                     "; 
      i=i+1; volcLat(i)= 44.805; volcLon(i)=147.135; volcElev(i)= 1206;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-05="; volcName(i)="Atsonupuri                     "; 
      i=i+1; volcLat(i)= 44.833; volcLon(i)=147.342; volcElev(i)= 1634;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-06-"; volcName(i)="Bogatyr Ridge                  "; 
      i=i+1; volcLat(i)= 45.030; volcLon(i)=147.208; volcElev(i)= -930;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="U0"; volcID(i)="0900-061"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 45.026; volcLon(i)=147.922; volcElev(i)= 1211;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S1"; volcID(i)="0900-07="; volcName(i)="Grozny Group                   "; 
      i=i+1; volcLat(i)= 45.097; volcLon(i)=148.024; volcElev(i)= 1132;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-08="; volcName(i)="Baransky                       "; 
      i=i+1; volcLat(i)= 45.250; volcLon(i)=148.350; volcElev(i)=  442;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-091"; volcName(i)="Golets-Tornyi Group            "; 
      i=i+1; volcLat(i)= 45.338; volcLon(i)=147.925; volcElev(i)= 1587;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0900-09="; volcName(i)="Chirip                         "; 
      i=i+1; volcLat(i)= 45.387; volcLon(i)=148.843; volcElev(i)= 1125;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-10="; volcName(i)="Medvezhia                      "; 
      i=i+1; volcLat(i)= 45.500; volcLon(i)=148.850; volcElev(i)= 1205;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-11-"; volcName(i)="Demon                          "; 
      i=i+1; volcLat(i)= 45.770; volcLon(i)=149.680; volcElev(i)= 1426;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-111"; volcName(i)="Ivao Group                     "; 
      i=i+1; volcLat(i)= 45.880; volcLon(i)=149.830; volcElev(i)=  542;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-112"; volcName(i)="Rudakov                        "; 
      i=i+1; volcLat(i)= 45.930; volcLon(i)=149.920; volcElev(i)=  998;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-113"; volcName(i)="Tri Sestry                     "; 
      i=i+1; volcLat(i)= 46.042; volcLon(i)=150.050; volcElev(i)= 1328;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-12="; volcName(i)="Kolokol Group                  "; 
      i=i+1; volcLat(i)= 46.100; volcLon(i)=150.500; volcElev(i)= -100;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="U0"; volcID(i)="0900-13-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 46.525; volcLon(i)=150.875; volcElev(i)=  742;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-15="; volcName(i)="Chirpoi                        "; 
      i=i+1; volcLat(i)= 46.470; volcLon(i)=151.280; volcElev(i)= -502;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="U0"; volcID(i)="0900-16-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 46.820; volcLon(i)=151.780; volcElev(i)= 1540;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-161"; volcName(i)="Milne                          "; 
      i=i+1; volcLat(i)= 46.830; volcLon(i)=151.750; volcElev(i)=  891;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-17="; volcName(i)="Goriaschaia Sopka              "; 
      i=i+1; volcLat(i)= 46.925; volcLon(i)=151.950; volcElev(i)=  624;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-18="; volcName(i)="Zavaritzki Caldera             "; 
      i=i+1; volcLat(i)= 47.120; volcLon(i)=152.250; volcElev(i)=  678;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-191"; volcName(i)="Urataman                       "; 
      i=i+1; volcLat(i)= 47.020; volcLon(i)=152.120; volcElev(i)= 1360;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0900-19="; volcName(i)="Prevo Peak                     "; 
      i=i+1; volcLat(i)= 47.350; volcLon(i)=152.475; volcElev(i)= 1172;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-20="; volcName(i)="Ketoi                          "; 
      i=i+1; volcLat(i)= 47.600; volcLon(i)=152.920; volcElev(i)=   36;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-211"; volcName(i)="Srednii                        "; 
      i=i+1; volcLat(i)= 47.520; volcLon(i)=152.800; volcElev(i)=  401;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-21="; volcName(i)="Ushishur                       "; 
      i=i+1; volcLat(i)= 47.770; volcLon(i)=153.020; volcElev(i)=  956;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-22="; volcName(i)="Rasshua                        "; 
      i=i+1; volcLat(i)= 48.080; volcLon(i)=153.330; volcElev(i)= -150;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="U0"; volcID(i)="0900-23="; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 48.092; volcLon(i)=153.200; volcElev(i)= 1496;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S1"; volcID(i)="0900-24="; volcName(i)="Sarychev Peak                  "; 
      i=i+1; volcLat(i)= 48.292; volcLon(i)=153.250; volcElev(i)=  551;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="M0"; volcID(i)="0900-25="; volcName(i)="Raikoke                        "; 
      i=i+1; volcLat(i)= 48.980; volcLon(i)=153.480; volcElev(i)=  724;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-26="; volcName(i)="Chirinkotan                    "; 
      i=i+1; volcLat(i)= 48.958; volcLon(i)=153.930; volcElev(i)= 1170;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-27="; volcName(i)="Ekarma                         "; 
      i=i+1; volcLat(i)= 48.875; volcLon(i)=154.175; volcElev(i)=  934;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-29="; volcName(i)="Sinarka                        "; 
      i=i+1; volcLat(i)= 49.120; volcLon(i)=154.508; volcElev(i)= 1145;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-30="; volcName(i)="Kharimkotan                    "; 
      i=i+1; volcLat(i)= 49.350; volcLon(i)=154.700; volcElev(i)= 1325;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-31="; volcName(i)="Tao-Rusyr Caldera              "; 
      i=i+1; volcLat(i)= 49.570; volcLon(i)=154.808; volcElev(i)= 1018;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-32="; volcName(i)="Nemo Peak                      "; 
      i=i+1; volcLat(i)= 50.200; volcLon(i)=154.980; volcElev(i)=  761;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-331"; volcName(i)="Shirinki                       "; 
      i=i+1; volcLat(i)= 50.270; volcLon(i)=155.250; volcElev(i)= 1772;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-34="; volcName(i)="Fuss Peak                      "; 
      i=i+1; volcLat(i)= 50.250; volcLon(i)=155.430; volcElev(i)= 1681;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-351"; volcName(i)="Lomonosov Group                "; 
      i=i+1; volcLat(i)= 50.130; volcLon(i)=155.370; volcElev(i)= 1345;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-35="; volcName(i)="Karpinsky Group                "; 
      i=i+1; volcLat(i)= 50.325; volcLon(i)=155.458; volcElev(i)= 1816;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="M2"; volcID(i)="0900-36="; volcName(i)="Chikurachki                    "; 
      i=i+1; volcLat(i)= 50.550; volcLon(i)=155.970; volcElev(i)= 1183;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S0"; volcID(i)="0900-37-"; volcName(i)="Vernadskii Ridge               "; 
      i=i+1; volcLat(i)= 50.680; volcLon(i)=156.020; volcElev(i)= 1156;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="S1"; volcID(i)="0900-38="; volcName(i)="Ebeko                          "; 
      i=i+1; volcLat(i)= 50.858; volcLon(i)=155.550; volcElev(i)= 2339;  volcLoc(i)="Kuril Is                      "; 
             volcESP_Code(i)="M1"; volcID(i)="0900-39="; volcName(i)="Alaid                          "; 
      i=i+1; volcLat(i)= 51.100; volcLon(i)=156.720; volcElev(i)=  503;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-001"; volcName(i)="Mashkovtsev                    "; 
      i=i+1; volcLat(i)= 51.300; volcLon(i)=156.870; volcElev(i)= 2156;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-01="; volcName(i)="Kambalny                       "; 
      i=i+1; volcLat(i)= 51.570; volcLon(i)=156.600; volcElev(i)=  705;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-021"; volcName(i)="Yavinsky                       "; 
      i=i+1; volcLat(i)= 51.450; volcLon(i)=156.970; volcElev(i)= 1070;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-022"; volcName(i)="Diky Greben                    "; 
      i=i+1; volcLat(i)= 51.450; volcLon(i)=157.120; volcElev(i)=   81;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-023"; volcName(i)="Kurile Lake                    "; 
      i=i+1; volcLat(i)= 51.357; volcLon(i)=156.750; volcElev(i)= 1812;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-02="; volcName(i)="Koshelev                       "; 
      i=i+1; volcLat(i)= 51.490; volcLon(i)=157.200; volcElev(i)= 1578;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-03="; volcName(i)="Ilyinsky                       "; 
      i=i+1; volcLat(i)= 51.650; volcLon(i)=157.350; volcElev(i)=  900;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-041"; volcName(i)="Kell                           "; 
      i=i+1; volcLat(i)= 51.750; volcLon(i)=157.270; volcElev(i)=  892;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-042"; volcName(i)="Belenkaya                      "; 
      i=i+1; volcLat(i)= 51.570; volcLon(i)=157.323; volcElev(i)= 1953;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-04="; volcName(i)="Zheltovsky                     "; 
      i=i+1; volcLat(i)= 51.880; volcLon(i)=157.380; volcElev(i)=  562;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-051"; volcName(i)="Ozernoy                        "; 
      i=i+1; volcLat(i)= 52.020; volcLon(i)=157.530; volcElev(i)=  681;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-052"; volcName(i)="Olkoviy Volc Group             "; 
      i=i+1; volcLat(i)= 52.063; volcLon(i)=157.703; volcElev(i)= 2090;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-053"; volcName(i)="Khodutka                       "; 
      i=i+1; volcLat(i)= 52.113; volcLon(i)=157.849; volcElev(i)= 1322;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-054"; volcName(i)="Piratkovsky                    "; 
      i=i+1; volcLat(i)= 52.146; volcLon(i)=157.322; volcElev(i)=  719;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-055"; volcName(i)="Ostanets                       "; 
      i=i+1; volcLat(i)= 52.220; volcLon(i)=157.428; volcElev(i)=  791;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-056"; volcName(i)="Otdelniy                       "; 
      i=i+1; volcLat(i)= 52.263; volcLon(i)=157.787; volcElev(i)=  858;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-057"; volcName(i)="Golaya                         "; 
      i=i+1; volcLat(i)= 52.355; volcLon(i)=157.827; volcElev(i)= 1910;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-058"; volcName(i)="Asacha                         "; 
      i=i+1; volcLat(i)= 52.430; volcLon(i)=157.930; volcElev(i)= 1234;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-059"; volcName(i)="Visokiy                        "; 
      i=i+1; volcLat(i)= 51.800; volcLon(i)=157.530; volcElev(i)= 1079;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-05="; volcName(i)="Ksudach                        "; 
      i=i+1; volcLat(i)= 52.453; volcLon(i)=158.195; volcElev(i)= 2322;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M1"; volcID(i)="1000-06="; volcName(i)="Mutnovsky                      "; 
      i=i+1; volcLat(i)= 52.558; volcLon(i)=158.030; volcElev(i)= 1829;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-07="; volcName(i)="Gorely                         "; 
      i=i+1; volcLat(i)= 52.570; volcLon(i)=157.020; volcElev(i)=  610;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-081"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 52.630; volcLon(i)=157.580; volcElev(i)= 1021;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-082"; volcName(i)="Tolmachev Dol                  "; 
      i=i+1; volcLat(i)= 52.700; volcLon(i)=158.280; volcElev(i)= 2173;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-083"; volcName(i)="Vilyuchik                      "; 
      i=i+1; volcLat(i)= 52.823; volcLon(i)=158.270; volcElev(i)=  870;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-084"; volcName(i)="Barkhatnaya Sopka              "; 
      i=i+1; volcLat(i)= 52.920; volcLon(i)=158.520; volcElev(i)=  450;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-085"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 52.880; volcLon(i)=158.300; volcElev(i)=  700;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-086"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 52.900; volcLon(i)=157.780; volcElev(i)= 1200;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-087"; volcName(i)="Bolshe-Bannaya                 "; 
      i=i+1; volcLat(i)= 52.543; volcLon(i)=157.335; volcElev(i)= 2475;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-08="; volcName(i)="Opala                          "; 
      i=i+1; volcLat(i)= 53.320; volcLon(i)=158.688; volcElev(i)= 3456;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-09="; volcName(i)="Koryaksky                      "; 
      i=i+1; volcLat(i)= 53.255; volcLon(i)=158.830; volcElev(i)= 2741;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-10="; volcName(i)="Avachinsky                     "; 
      i=i+1; volcLat(i)= 53.637; volcLon(i)=158.922; volcElev(i)= 2285;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-11="; volcName(i)="Dzenzursky                     "; 
      i=i+1; volcLat(i)= 53.750; volcLon(i)=158.450; volcElev(i)=  520;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-121"; volcName(i)="Veer                           "; 
      i=i+1; volcLat(i)= 53.830; volcLon(i)=158.050; volcElev(i)= 1150;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-122"; volcName(i)="Kostakan                       "; 
      i=i+1; volcLat(i)= 53.905; volcLon(i)=158.070; volcElev(i)= 2278;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-123"; volcName(i)="Bakening                       "; 
      i=i+1; volcLat(i)= 53.905; volcLon(i)=158.385; volcElev(i)= 1567;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-124"; volcName(i)="Zavaritsky                     "; 
      i=i+1; volcLat(i)= 53.980; volcLon(i)=159.450; volcElev(i)= 1180;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-125"; volcName(i)="Akademia Nauk                  "; 
      i=i+1; volcLat(i)= 53.590; volcLon(i)=159.147; volcElev(i)= 2958;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-12="; volcName(i)="Zhupanovsky                    "; 
      i=i+1; volcLat(i)= 54.050; volcLon(i)=159.450; volcElev(i)= 1536;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S1"; volcID(i)="1000-13="; volcName(i)="Karymsky                       "; 
      i=i+1; volcLat(i)= 54.130; volcLon(i)=159.670; volcElev(i)= 1560;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-14="; volcName(i)="Maly Semiachik                 "; 
      i=i+1; volcLat(i)= 54.320; volcLon(i)=160.020; volcElev(i)= 1720;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-15="; volcName(i)="Bolshoi Semiachik              "; 
      i=i+1; volcLat(i)= 54.530; volcLon(i)=159.800; volcElev(i)= 2353;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-16-"; volcName(i)="Taunshits                      "; 
      i=i+1; volcLat(i)= 54.500; volcLon(i)=159.970; volcElev(i)= 1617;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-17="; volcName(i)="Uzon                           "; 
      i=i+1; volcLat(i)= 54.487; volcLon(i)=160.253; volcElev(i)= 1552;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-18="; volcName(i)="Kikhpinych                     "; 
      i=i+1; volcLat(i)= 54.593; volcLon(i)=160.273; volcElev(i)= 1856;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-19="; volcName(i)="Krasheninnikov                 "; 
      i=i+1; volcLat(i)= 54.920; volcLon(i)=160.630; volcElev(i)= 2020;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-201"; volcName(i)="Schmidt                        "; 
      i=i+1; volcLat(i)= 54.753; volcLon(i)=160.527; volcElev(i)= 3528;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-20="; volcName(i)="Kronotsky                      "; 
      i=i+1; volcLat(i)= 54.973; volcLon(i)=160.702; volcElev(i)= 2576;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-21="; volcName(i)="Gamchen                        "; 
      i=i+1; volcLat(i)= 55.070; volcLon(i)=160.770; volcElev(i)= 2161;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-221"; volcName(i)="Vysoky                         "; 
      i=i+1; volcLat(i)= 55.032; volcLon(i)=160.720; volcElev(i)= 2070;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-22="; volcName(i)="Komarov                        "; 
      i=i+1; volcLat(i)= 55.920; volcLon(i)=161.750; volcElev(i)=    0;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-232"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 55.130; volcLon(i)=160.320; volcElev(i)= 2376;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-23="; volcName(i)="Kizimen                        "; 
      i=i+1; volcLat(i)= 55.755; volcLon(i)=160.527; volcElev(i)= 2923;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-241"; volcName(i)="Udina                          "; 
      i=i+1; volcLat(i)= 55.862; volcLon(i)=160.603; volcElev(i)= 3081;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-242"; volcName(i)="Zimina                         "; 
      i=i+1; volcLat(i)= 55.830; volcLon(i)=160.330; volcElev(i)= 3682;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M1"; volcID(i)="1000-24="; volcName(i)="Tolbachik                      "; 
      i=i+1; volcLat(i)= 56.020; volcLon(i)=160.593; volcElev(i)= 4585;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-251"; volcName(i)="Kamen                          "; 
      i=i+1; volcLat(i)= 55.978; volcLon(i)=160.587; volcElev(i)= 2882;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S2"; volcID(i)="1000-25="; volcName(i)="Bezymianny                     "; 
      i=i+1; volcLat(i)= 56.070; volcLon(i)=160.470; volcElev(i)= 3943;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-261"; volcName(i)="Ushkovsky                      "; 
      i=i+1; volcLat(i)= 56.057; volcLon(i)=160.638; volcElev(i)= 4835;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M1"; volcID(i)="1000-26="; volcName(i)="Kliuchevskoi                   "; 
      i=i+1; volcLat(i)= 55.420; volcLon(i)=167.330; volcElev(i)= -300;  volcLoc(i)="Kamchatka-E of                "; 
             volcESP_Code(i)="U0"; volcID(i)="1000-271"; volcName(i)="Piip                           "; 
      i=i+1; volcLat(i)= 54.750; volcLon(i)=157.380; volcElev(i)= 2000;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-272"; volcName(i)="Khangar                        "; 
      i=i+1; volcLat(i)= 55.550; volcLon(i)=157.470; volcElev(i)= 1868;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-273"; volcName(i)="Cherpuk Group                  "; 
      i=i+1; volcLat(i)= 56.653; volcLon(i)=161.360; volcElev(i)= 3283;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S1"; volcID(i)="1000-27="; volcName(i)="Shiveluch                      "; 
      i=i+1; volcLat(i)= 55.680; volcLon(i)=157.730; volcElev(i)= 3621;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-28="; volcName(i)="Ichinsky                       "; 
      i=i+1; volcLat(i)= 55.820; volcLon(i)=157.980; volcElev(i)= 1802;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-29-"; volcName(i)="Maly Payalpan                  "; 
      i=i+1; volcLat(i)= 55.880; volcLon(i)=157.780; volcElev(i)= 1906;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-30-"; volcName(i)="Bolshoi Payalpan               "; 
      i=i+1; volcLat(i)= 55.200; volcLon(i)=158.470; volcElev(i)= 1236;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-31-"; volcName(i)="Plosky                         "; 
      i=i+1; volcLat(i)= 55.430; volcLon(i)=158.650; volcElev(i)= 1956;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-32-"; volcName(i)="Akhtang                        "; 
      i=i+1; volcLat(i)= 55.580; volcLon(i)=158.380; volcElev(i)= 2016;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-33-"; volcName(i)="Kozyrevsky                     "; 
      i=i+1; volcLat(i)= 55.650; volcLon(i)=158.800; volcElev(i)= 1442;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-34-"; volcName(i)="Romanovka                      "; 
      i=i+1; volcLat(i)= 56.080; volcLon(i)=158.380; volcElev(i)= 1692;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-35-"; volcName(i)="Uksichan                       "; 
      i=i+1; volcLat(i)= 56.470; volcLon(i)=157.800; volcElev(i)= 1401;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-36-"; volcName(i)="Bolshoi-Kekuknaysky            "; 
      i=i+1; volcLat(i)= 56.370; volcLon(i)=158.370; volcElev(i)=  915;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-37-"; volcName(i)="Kulkev                         "; 
      i=i+1; volcLat(i)= 56.330; volcLon(i)=158.670; volcElev(i)= 1170;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-38-"; volcName(i)="Geodesistoy                    "; 
      i=i+1; volcLat(i)= 56.320; volcLon(i)=158.830; volcElev(i)= 1828;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-39-"; volcName(i)="Anaun                          "; 
      i=i+1; volcLat(i)= 56.370; volcLon(i)=159.030; volcElev(i)= 1554;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-40-"; volcName(i)="Krainy                         "; 
      i=i+1; volcLat(i)= 56.400; volcLon(i)=158.850; volcElev(i)= 1377;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-41-"; volcName(i)="Kekurny                        "; 
      i=i+1; volcLat(i)= 56.570; volcLon(i)=158.520; volcElev(i)= 1046;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-42-"; volcName(i)="Eggella                        "; 
      i=i+1; volcLat(i)= 56.820; volcLon(i)=158.950; volcElev(i)= 1185;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-43-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 56.520; volcLon(i)=159.530; volcElev(i)= 1400;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-44-"; volcName(i)="Verkhovoy                      "; 
      i=i+1; volcLat(i)= 56.700; volcLon(i)=159.650; volcElev(i)= 2598;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-45-"; volcName(i)="Alney-Chashakondzha            "; 
      i=i+1; volcLat(i)= 56.820; volcLon(i)=159.670; volcElev(i)= 1778;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-46-"; volcName(i)="Cherny                         "; 
      i=i+1; volcLat(i)= 56.850; volcLon(i)=159.800; volcElev(i)= 1427;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-47-"; volcName(i)="Pogranychny                    "; 
      i=i+1; volcLat(i)= 56.880; volcLon(i)=159.950; volcElev(i)= 1349;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-48-"; volcName(i)="Zaozerny                       "; 
      i=i+1; volcLat(i)= 56.970; volcLon(i)=159.780; volcElev(i)= 1244;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-49-"; volcName(i)="Bliznets                       "; 
      i=i+1; volcLat(i)= 57.100; volcLon(i)=159.930; volcElev(i)= 1527;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-50-"; volcName(i)="Kebeney                        "; 
      i=i+1; volcLat(i)= 57.130; volcLon(i)=160.400; volcElev(i)=  965;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-51-"; volcName(i)="Fedotych                       "; 
      i=i+1; volcLat(i)= 57.150; volcLon(i)=161.080; volcElev(i)=  379;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-511"; volcName(i)="Shisheika                      "; 
      i=i+1; volcLat(i)= 57.200; volcLon(i)=159.830; volcElev(i)=  765;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-512"; volcName(i)="Terpuk                         "; 
      i=i+1; volcLat(i)= 57.270; volcLon(i)=160.080; volcElev(i)= 1241;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-52-"; volcName(i)="Sedankinsky                    "; 
      i=i+1; volcLat(i)= 57.300; volcLon(i)=159.830; volcElev(i)= 1333;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-53-"; volcName(i)="Leutongey                      "; 
      i=i+1; volcLat(i)= 57.320; volcLon(i)=159.970; volcElev(i)= 1533;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-54-"; volcName(i)="Tuzovsky                       "; 
      i=i+1; volcLat(i)= 57.330; volcLon(i)=160.200; volcElev(i)= 2125;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-55-"; volcName(i)="Gorny Institute                "; 
      i=i+1; volcLat(i)= 57.350; volcLon(i)=160.970; volcElev(i)=  583;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-551"; volcName(i)="Kinenin                        "; 
      i=i+1; volcLat(i)= 57.350; volcLon(i)=161.370; volcElev(i)=  265;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-552"; volcName(i)="Bliznetsy                      "; 
      i=i+1; volcLat(i)= 57.400; volcLon(i)=160.100; volcElev(i)= 1559;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-56-"; volcName(i)="Titila                         "; 
      i=i+1; volcLat(i)= 57.470; volcLon(i)=160.250; volcElev(i)= 1641;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-57-"; volcName(i)="Mezhdusopochny                 "; 
      i=i+1; volcLat(i)= 57.450; volcLon(i)=160.370; volcElev(i)= 2525;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-58-"; volcName(i)="Shishel                        "; 
      i=i+1; volcLat(i)= 57.550; volcLon(i)=160.530; volcElev(i)= 1381;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-59-"; volcName(i)="Elovsky                        "; 
      i=i+1; volcLat(i)= 57.700; volcLon(i)=160.400; volcElev(i)= 1853;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-60-"; volcName(i)="Alngey                         "; 
      i=i+1; volcLat(i)= 57.700; volcLon(i)=160.580; volcElev(i)= 1643;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-61-"; volcName(i)="Uka                            "; 
      i=i+1; volcLat(i)= 57.800; volcLon(i)=160.670; volcElev(i)= 1582;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-62-"; volcName(i)="Kaileney                       "; 
      i=i+1; volcLat(i)= 57.830; volcLon(i)=160.250; volcElev(i)= 1255;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-63-"; volcName(i)="Plosky                         "; 
      i=i+1; volcLat(i)= 57.880; volcLon(i)=160.530; volcElev(i)= 2080;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-64-"; volcName(i)="Bely                           "; 
      i=i+1; volcLat(i)= 57.970; volcLon(i)=160.650; volcElev(i)= 1764;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-65-"; volcName(i)="Nylgimelkin                    "; 
      i=i+1; volcLat(i)= 58.020; volcLon(i)=160.800; volcElev(i)= 2169;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-66-"; volcName(i)="Snezhniy                       "; 
      i=i+1; volcLat(i)= 58.080; volcLon(i)=160.770; volcElev(i)= 2300;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-67-"; volcName(i)="Iktunup                        "; 
      i=i+1; volcLat(i)= 58.130; volcLon(i)=160.820; volcElev(i)= 2171;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1000-671"; volcName(i)="Spokoiny                       "; 
      i=i+1; volcLat(i)= 58.180; volcLon(i)=160.820; volcElev(i)= 2552;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-68-"; volcName(i)="Ostry                          "; 
      i=i+1; volcLat(i)= 58.200; volcLon(i)=160.970; volcElev(i)= 2169;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-69-"; volcName(i)="Snegovoy                       "; 
      i=i+1; volcLat(i)= 58.280; volcLon(i)=160.870; volcElev(i)= 1936;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-70-"; volcName(i)="Severny                        "; 
      i=i+1; volcLat(i)= 58.400; volcLon(i)=161.080; volcElev(i)= 1340;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-71-"; volcName(i)="Iettunup                       "; 
      i=i+1; volcLat(i)= 58.370; volcLon(i)=160.620; volcElev(i)= 1225;  volcLoc(i)="Kamchatka                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1000-72-"; volcName(i)="Voyampolsky                    "; 
      i=i+1; volcLat(i)= 47.000; volcLon(i)=137.500; volcElev(i)=    0;  volcLoc(i)="Russia-SE                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1002-01-"; volcName(i)="Sikhote-Alin                   "; 
      i=i+1; volcLat(i)= 56.280; volcLon(i)=117.770; volcElev(i)= 2180;  volcLoc(i)="Russia-SE                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1002-03-"; volcName(i)="Udokan Plateau                 "; 
      i=i+1; volcLat(i)= 53.700; volcLon(i)=113.300; volcElev(i)= 1250;  volcLoc(i)="Russia-SE                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1002-04-"; volcName(i)="Vitim Plateau                  "; 
      i=i+1; volcLat(i)= 51.500; volcLon(i)=102.500; volcElev(i)= 1200;  volcLoc(i)="Russia-SE                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1002-05-"; volcName(i)="Tunkin Depression              "; 
      i=i+1; volcLat(i)= 52.700; volcLon(i)= 98.980; volcElev(i)= 2077;  volcLoc(i)="Russia-SE                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1002-06-"; volcName(i)="Oka Plateau                    "; 
      i=i+1; volcLat(i)= 52.520; volcLon(i)= 98.600; volcElev(i)= 2765;  volcLoc(i)="Russia-SE                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1002-07-"; volcName(i)="Azas Plateau                   "; 
      i=i+1; volcLat(i)= 48.170; volcLon(i)= 99.700; volcElev(i)= 2400;  volcLoc(i)="Mongolia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1003-01-"; volcName(i)="Taryatu-Chulutu                "; 
      i=i+1; volcLat(i)= 48.670; volcLon(i)=102.750; volcElev(i)= 1886;  volcLoc(i)="Mongolia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1003-02-"; volcName(i)="Khanuy Gol                     "; 
      i=i+1; volcLat(i)= 47.120; volcLon(i)=109.080; volcElev(i)= 1162;  volcLoc(i)="Mongolia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1003-03-"; volcName(i)="Bus-Obo                        "; 
      i=i+1; volcLat(i)= 45.330; volcLon(i)=114.000; volcElev(i)= 1778;  volcLoc(i)="Mongolia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1003-04-"; volcName(i)="Dariganga Volc Field           "; 
      i=i+1; volcLat(i)= 45.280; volcLon(i)=106.700; volcElev(i)= 1120;  volcLoc(i)="Mongolia                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1003-05-"; volcName(i)="Middle Gobi                    "; 
      i=i+1; volcLat(i)= 42.900; volcLon(i)= 89.250; volcElev(i)=    0;  volcLoc(i)="China-W                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1004-01-"; volcName(i)="Turfan                         "; 
      i=i+1; volcLat(i)= 42.500; volcLon(i)= 82.500; volcElev(i)=    0;  volcLoc(i)="China-W                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1004-02-"; volcName(i)="Tianshan Volc Group            "; 
      i=i+1; volcLat(i)= 35.520; volcLon(i)= 80.200; volcElev(i)= 5808;  volcLoc(i)="China-W                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1004-03-"; volcName(i)="Kunlun Volc Group              "; 
      i=i+1; volcLat(i)= 35.850; volcLon(i)= 91.700; volcElev(i)= 5400;  volcLoc(i)="China-W                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1004-04-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 41.470; volcLon(i)=113.000; volcElev(i)= 1700;  volcLoc(i)="China-E                       "; 
             volcESP_Code(i)="M0"; volcID(i)="1005-01-"; volcName(i)="Honggeertu                     "; 
      i=i+1; volcLat(i)= 47.500; volcLon(i)=120.700; volcElev(i)=    0;  volcLoc(i)="China-E                       "; 
             volcESP_Code(i)="M0"; volcID(i)="1005-011"; volcName(i)="Arshan                         "; 
      i=i+1; volcLat(i)= 49.370; volcLon(i)=125.920; volcElev(i)=  670;  volcLoc(i)="China-E                       "; 
             volcESP_Code(i)="M0"; volcID(i)="1005-02-"; volcName(i)="Keluo Group                    "; 
      i=i+1; volcLat(i)= 48.720; volcLon(i)=126.120; volcElev(i)=  597;  volcLoc(i)="China-E                       "; 
             volcESP_Code(i)="M0"; volcID(i)="1005-03-"; volcName(i)="Wudalianchi                    "; 
      i=i+1; volcLat(i)= 44.080; volcLon(i)=128.830; volcElev(i)= 1000;  volcLoc(i)="China-E                       "; 
             volcESP_Code(i)="M0"; volcID(i)="1005-04-"; volcName(i)="Jingbo                         "; 
      i=i+1; volcLat(i)= 42.330; volcLon(i)=126.500; volcElev(i)= 1000;  volcLoc(i)="China-E                       "; 
             volcESP_Code(i)="M0"; volcID(i)="1005-05-"; volcName(i)="Longgang Group                 "; 
      i=i+1; volcLat(i)= 41.980; volcLon(i)=128.080; volcElev(i)= 2744;  volcLoc(i)="China-E                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1005-06-"; volcName(i)="Changbaishan                   "; 
      i=i+1; volcLat(i)= 41.330; volcLon(i)=128.000; volcElev(i)=    0;  volcLoc(i)="Korea                         "; 
             volcESP_Code(i)="S0"; volcID(i)="1006-01-"; volcName(i)="Xianjindao                     "; 
      i=i+1; volcLat(i)= 38.330; volcLon(i)=127.330; volcElev(i)=  452;  volcLoc(i)="Korea                         "; 
             volcESP_Code(i)="M0"; volcID(i)="1006-02-"; volcName(i)="Ch'uga-ryong                   "; 
      i=i+1; volcLat(i)= 37.500; volcLon(i)=130.870; volcElev(i)=  984;  volcLoc(i)="Korea                         "; 
             volcESP_Code(i)="S0"; volcID(i)="1006-03-"; volcName(i)="Ulreung                        "; 
      i=i+1; volcLat(i)= 33.370; volcLon(i)=126.530; volcElev(i)= 1950;  volcLoc(i)="Korea                         "; 
             volcESP_Code(i)="M0"; volcID(i)="1006-04-"; volcName(i)="Halla                          "; 
      i=i+1; volcLat(i)= 52.350; volcLon(i)=175.911; volcElev(i)=  656;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-01-"; volcName(i)="Buldir                         "; 
      i=i+1; volcLat(i)= 52.103; volcLon(i)=177.602; volcElev(i)= 1220;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="M0"; volcID(i)="1101-02-"; volcName(i)="Kiska                          "; 
      i=i+1; volcLat(i)= 52.015; volcLon(i)=178.136; volcElev(i)= 1160;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-03-"; volcName(i)="Segula                         "; 
      i=i+1; volcLat(i)= 51.970; volcLon(i)=178.330; volcElev(i)=  328;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="M0"; volcID(i)="1101-04-"; volcName(i)="Davidof                        "; 
      i=i+1; volcLat(i)= 51.950; volcLon(i)=178.543; volcElev(i)= 1174;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-05-"; volcName(i)="Little Sitkin                  "; 
      i=i+1; volcLat(i)= 51.930; volcLon(i)=179.580; volcElev(i)= 1221;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-06-"; volcName(i)="Semisopochnoi                  "; 
      i=i+1; volcLat(i)= 51.790; volcLon(i)=181.206; volcElev(i)= 1573;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-07-"; volcName(i)="Gareloi                        "; 
      i=i+1; volcLat(i)= 51.885; volcLon(i)=181.854; volcElev(i)= 1806;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-08-"; volcName(i)="Tanaga                         "; 
      i=i+1; volcLat(i)= 51.873; volcLon(i)=181.994; volcElev(i)= 1449;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-09-"; volcName(i)="Takawangha                     "; 
      i=i+1; volcLat(i)= 51.910; volcLon(i)=182.562; volcElev(i)=  738;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-10-"; volcName(i)="Bobrof                         "; 
      i=i+1; volcLat(i)= 51.923; volcLon(i)=182.832; volcElev(i)= 1307;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="M0"; volcID(i)="1101-11-"; volcName(i)="Kanaga                         "; 
      i=i+1; volcLat(i)= 51.944; volcLon(i)=183.253; volcElev(i)= 1196;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-111"; volcName(i)="Moffett                        "; 
      i=i+1; volcLat(i)= 52.076; volcLon(i)=183.870; volcElev(i)= 1740;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S1"; volcID(i)="1101-12-"; volcName(i)="Great Sitkin                   "; 
      i=i+1; volcLat(i)= 52.177; volcLon(i)=184.492; volcElev(i)=  314;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S2"; volcID(i)="1101-13-"; volcName(i)="Kasatochi                      "; 
      i=i+1; volcLat(i)= 52.220; volcLon(i)=184.870; volcElev(i)=  273;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-14-"; volcName(i)="Koniuji                        "; 
      i=i+1; volcLat(i)= 52.050; volcLon(i)=185.050; volcElev(i)=  560;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-15-"; volcName(i)="Sergief                        "; 
      i=i+1; volcLat(i)= 52.332; volcLon(i)=185.863; volcElev(i)= 1451;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S1"; volcID(i)="1101-16-"; volcName(i)="Atka                           "; 
      i=i+1; volcLat(i)= 52.381; volcLon(i)=185.846; volcElev(i)= 1533;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S1"; volcID(i)="1101-161"; volcName(i)="Korovin                        "; 
      i=i+1; volcLat(i)= 52.315; volcLon(i)=187.490; volcElev(i)= 1054;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-18-"; volcName(i)="Seguam                         "; 
      i=i+1; volcLat(i)= 52.500; volcLon(i)=188.748; volcElev(i)= 1066;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S1"; volcID(i)="1101-19-"; volcName(i)="Amukta                         "; 
      i=i+1; volcLat(i)= 52.577; volcLon(i)=188.870; volcElev(i)= 1142;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-20-"; volcName(i)="Chagulak                       "; 
      i=i+1; volcLat(i)= 52.643; volcLon(i)=189.371; volcElev(i)=  550;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S2"; volcID(i)="1101-21-"; volcName(i)="Yunaska                        "; 
      i=i+1; volcLat(i)= 52.742; volcLon(i)=189.889; volcElev(i)= 1280;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-22-"; volcName(i)="Herbert                        "; 
      i=i+1; volcLat(i)= 52.894; volcLon(i)=189.946; volcElev(i)= 1620;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S1"; volcID(i)="1101-23-"; volcName(i)="Carlisle                       "; 
      i=i+1; volcLat(i)= 52.825; volcLon(i)=190.056; volcElev(i)= 1730;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S2"; volcID(i)="1101-24-"; volcName(i)="Cleveland                      "; 
      i=i+1; volcLat(i)= 52.830; volcLon(i)=190.230; volcElev(i)= 1170;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-241"; volcName(i)="Tana                           "; 
      i=i+1; volcLat(i)= 53.065; volcLon(i)=190.230; volcElev(i)=  888;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-25-"; volcName(i)="Uliaga                         "; 
      i=i+1; volcLat(i)= 52.974; volcLon(i)=190.280; volcElev(i)=  893;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-26-"; volcName(i)="Kagamil                        "; 
      i=i+1; volcLat(i)= 53.130; volcLon(i)=191.307; volcElev(i)= 2149;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-27-"; volcName(i)="Vsevidof                       "; 
      i=i+1; volcLat(i)= 53.157; volcLon(i)=191.461; volcElev(i)= 1984;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-28-"; volcName(i)="Recheschnoi                    "; 
      i=i+1; volcLat(i)= 53.430; volcLon(i)=191.870; volcElev(i)= 1073;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S2"; volcID(i)="1101-29-"; volcName(i)="Okmok                          "; 
      i=i+1; volcLat(i)= 53.930; volcLon(i)=191.970; volcElev(i)=  150;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S2"; volcID(i)="1101-30-"; volcName(i)="Bogoslof                       "; 
      i=i+1; volcLat(i)= 53.891; volcLon(i)=193.077; volcElev(i)= 1800;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S1"; volcID(i)="1101-31-"; volcName(i)="Makushin                       "; 
      i=i+1; volcLat(i)= 54.134; volcLon(i)=194.014; volcElev(i)= 1303;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S1"; volcID(i)="1101-32-"; volcName(i)="Akutan                         "; 
      i=i+1; volcLat(i)= 54.518; volcLon(i)=195.350; volcElev(i)= 1654;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S2"; volcID(i)="1101-34-"; volcName(i)="Westdahl                       "; 
      i=i+1; volcLat(i)= 54.650; volcLon(i)=195.570; volcElev(i)= 1112;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-35-"; volcName(i)="Fisher                         "; 
      i=i+1; volcLat(i)= 54.756; volcLon(i)=196.030; volcElev(i)= 2857;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S2"; volcID(i)="1101-36-"; volcName(i)="Shishaldin                     "; 
      i=i+1; volcLat(i)= 54.765; volcLon(i)=196.277; volcElev(i)= 2446;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-37-"; volcName(i)="Isanotski                      "; 
      i=i+1; volcLat(i)= 54.800; volcLon(i)=196.411; volcElev(i)= 1871;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-38-"; volcName(i)="Roundtop                       "; 
      i=i+1; volcLat(i)= 55.424; volcLon(i)=196.851; volcElev(i)=  488;  volcLoc(i)="Aleutian Is                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1101-39-"; volcName(i)="Amak                           "; 
      i=i+1; volcLat(i)= 55.082; volcLon(i)=197.186; volcElev(i)= 2012;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-01-"; volcName(i)="Frosty                         "; 
      i=i+1; volcLat(i)= 55.168; volcLon(i)=197.728; volcElev(i)= 1506;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-011"; volcName(i)="Dutton                         "; 
      i=i+1; volcLat(i)= 55.341; volcLon(i)=197.921; volcElev(i)= 1436;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-02-"; volcName(i)="Emmons Lake                    "; 
      i=i+1; volcLat(i)= 55.420; volcLon(i)=198.113; volcElev(i)= 2519;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S2"; volcID(i)="1102-03-"; volcName(i)="Pavlof                         "; 
      i=i+1; volcLat(i)= 55.453; volcLon(i)=198.157; volcElev(i)= 2142;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-04-"; volcName(i)="Pavlof Sister                  "; 
      i=i+1; volcLat(i)= 55.641; volcLon(i)=198.786; volcElev(i)= 1354;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-05-"; volcName(i)="Dana                           "; 
      i=i+1; volcLat(i)= 55.913; volcLon(i)=199.959; volcElev(i)= 1323;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-051"; volcName(i)="Stepovak Bay 2                 "; 
      i=i+1; volcLat(i)= 55.929; volcLon(i)=199.998; volcElev(i)= 1555;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-052"; volcName(i)="Stepovak Bay 3                 "; 
      i=i+1; volcLat(i)= 55.954; volcLon(i)=200.046; volcElev(i)= 1557;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-053"; volcName(i)="Stepovak Bay 4                 "; 
      i=i+1; volcLat(i)= 56.011; volcLon(i)=200.203; volcElev(i)= 1895;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-06-"; volcName(i)="Kupreanof                      "; 
      i=i+1; volcLat(i)= 56.170; volcLon(i)=200.620; volcElev(i)= 2507;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S1"; volcID(i)="1102-07-"; volcName(i)="Veniaminof                     "; 
      i=i+1; volcLat(i)= 56.552; volcLon(i)=201.215; volcElev(i)= 1032;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-08-"; volcName(i)="Black Peak                     "; 
      i=i+1; volcLat(i)= 56.880; volcLon(i)=201.830; volcElev(i)= 1341;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S1"; volcID(i)="1102-09-"; volcName(i)="Aniakchak                      "; 
      i=i+1; volcLat(i)= 57.019; volcLon(i)=202.815; volcElev(i)= 1345;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-10-"; volcName(i)="Yantarni                       "; 
      i=i+1; volcLat(i)= 57.135; volcLon(i)=203.010; volcElev(i)= 2221;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-11-"; volcName(i)="Chiginagak                     "; 
      i=i+1; volcLat(i)= 57.203; volcLon(i)=203.255; volcElev(i)= 1677;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-12-"; volcName(i)="Kialagvik                      "; 
      i=i+1; volcLat(i)= 57.751; volcLon(i)=203.632; volcElev(i)= 1474;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-13-"; volcName(i)="Ugashik-Peulik                 "; 
      i=i+1; volcLat(i)= 57.832; volcLon(i)=203.490; volcElev(i)=   91;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S1"; volcID(i)="1102-131"; volcName(i)="Ukinrek Maars                  "; 
      i=i+1; volcLat(i)= 57.870; volcLon(i)=204.580; volcElev(i)=  300;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-132"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 58.172; volcLon(i)=204.639; volcElev(i)= 1863;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-14-"; volcName(i)="Martin                         "; 
      i=i+1; volcLat(i)= 58.195; volcLon(i)=204.747; volcElev(i)= 2165;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-15-"; volcName(i)="Mageik                         "; 
      i=i+1; volcLat(i)= 58.236; volcLon(i)=204.900; volcElev(i)= 1864;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S2"; volcID(i)="1102-16-"; volcName(i)="Trident                        "; 
      i=i+1; volcLat(i)= 58.280; volcLon(i)=205.037; volcElev(i)= 2047;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-17-"; volcName(i)="Katmai                         "; 
      i=i+1; volcLat(i)= 58.270; volcLon(i)=204.843; volcElev(i)=  841;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-18-"; volcName(i)="Novarupta                      "; 
      i=i+1; volcLat(i)= 58.354; volcLon(i)=204.908; volcElev(i)= 2317;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-19-"; volcName(i)="Griggs                         "; 
      i=i+1; volcLat(i)= 58.336; volcLon(i)=205.318; volcElev(i)= 2162;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-20-"; volcName(i)="Snowy Mountain                 "; 
      i=i+1; volcLat(i)= 58.418; volcLon(i)=205.551; volcElev(i)= 2287;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-21-"; volcName(i)="Denison                        "; 
      i=i+1; volcLat(i)= 58.430; volcLon(i)=205.600; volcElev(i)= 2272;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-22-"; volcName(i)="Steller                        "; 
      i=i+1; volcLat(i)= 58.453; volcLon(i)=205.645; volcElev(i)= 2043;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-23-"; volcName(i)="Kukak                          "; 
      i=i+1; volcLat(i)= 58.608; volcLon(i)=205.972; volcElev(i)=  901;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-25-"; volcName(i)="Kaguyak                        "; 
      i=i+1; volcLat(i)= 58.770; volcLon(i)=206.328; volcElev(i)= 2105;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S1"; volcID(i)="1102-26-"; volcName(i)="Fourpeaked                     "; 
      i=i+1; volcLat(i)= 58.855; volcLon(i)=206.458; volcElev(i)= 2140;  volcLoc(i)="Alaska Peninsula              "; 
             volcESP_Code(i)="S0"; volcID(i)="1102-27-"; volcName(i)="Douglas                        "; 
      i=i+1; volcLat(i)= 59.363; volcLon(i)=206.570; volcElev(i)= 1252;  volcLoc(i)="Alaska-SW                     "; 
             volcESP_Code(i)="S2"; volcID(i)="1103-01-"; volcName(i)="Augustine                      "; 
      i=i+1; volcLat(i)= 60.032; volcLon(i)=206.910; volcElev(i)= 3053;  volcLoc(i)="Alaska-SW                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1103-02-"; volcName(i)="Iliamna                        "; 
      i=i+1; volcLat(i)= 60.485; volcLon(i)=207.258; volcElev(i)= 3108;  volcLoc(i)="Alaska-SW                     "; 
             volcESP_Code(i)="S2"; volcID(i)="1103-03-"; volcName(i)="Redoubt                        "; 
      i=i+1; volcLat(i)= 61.299; volcLon(i)=207.749; volcElev(i)= 3374;  volcLoc(i)="Alaska-SW                     "; 
             volcESP_Code(i)="S2"; volcID(i)="1103-04-"; volcName(i)="Spurr                          "; 
      i=i+1; volcLat(i)= 61.640; volcLon(i)=207.589; volcElev(i)= 3034;  volcLoc(i)="Alaska-SW                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1103-05-"; volcName(i)="Hayes                          "; 
      i=i+1; volcLat(i)= 57.180; volcLon(i)=189.700; volcElev(i)=  203;  volcLoc(i)="Alaska-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1104-01-"; volcName(i)="St. Paul Island                "; 
      i=i+1; volcLat(i)= 60.020; volcLon(i)=193.670; volcElev(i)=  511;  volcLoc(i)="Alaska-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1104-02-"; volcName(i)="Nunivak Island                 "; 
      i=i+1; volcLat(i)= 61.430; volcLon(i)=195.530; volcElev(i)=  190;  volcLoc(i)="Alaska-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1104-03-"; volcName(i)="Ingakslugwat Hills             "; 
      i=i+1; volcLat(i)= 63.450; volcLon(i)=197.880; volcElev(i)=  715;  volcLoc(i)="Alaska-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1104-04-"; volcName(i)="St. Michael                    "; 
      i=i+1; volcLat(i)= 63.600; volcLon(i)=189.570; volcElev(i)=  673;  volcLoc(i)="Alaska-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1104-05-"; volcName(i)="Kookooligit Mountains          "; 
      i=i+1; volcLat(i)= 65.600; volcLon(i)=196.080; volcElev(i)=  610;  volcLoc(i)="Alaska-W                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1104-06-"; volcName(i)="Imuruk Lake                    "; 
      i=i+1; volcLat(i)= 64.070; volcLon(i)=211.580; volcElev(i)=  830;  volcLoc(i)="Alaska-E                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1105-001"; volcName(i)="Buzzard Creek                  "; 
      i=i+1; volcLat(i)= 62.220; volcLon(i)=215.870; volcElev(i)= 4949;  volcLoc(i)="Alaska-E                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1105-01-"; volcName(i)="Sanford                        "; 
      i=i+1; volcLat(i)= 62.000; volcLon(i)=215.980; volcElev(i)= 4317;  volcLoc(i)="Alaska-E                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1105-02-"; volcName(i)="Wrangell                       "; 
      i=i+1; volcLat(i)= 62.130; volcLon(i)=216.920; volcElev(i)= 2755;  volcLoc(i)="Alaska-E                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1105-021"; volcName(i)="Gordon                         "; 
      i=i+1; volcLat(i)= 61.380; volcLon(i)=218.250; volcElev(i)= 5005;  volcLoc(i)="Alaska-E                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1105-03-"; volcName(i)="Churchill                      "; 
      i=i+1; volcLat(i)= 57.050; volcLon(i)=224.250; volcElev(i)=  970;  volcLoc(i)="Alaska-E                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1105-04-"; volcName(i)="Edgecumbe                      "; 
      i=i+1; volcLat(i)= 56.500; volcLon(i)=226.900; volcElev(i)=   15;  volcLoc(i)="Alaska-E                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1105-05-"; volcName(i)="Duncan Canal                   "; 
      i=i+1; volcLat(i)= 55.250; volcLon(i)=226.700; volcElev(i)=   50;  volcLoc(i)="Alaska-E                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1105-06-"; volcName(i)="Tlevak Strait-Suemez Is.       "; 
      i=i+1; volcLat(i)= 55.320; volcLon(i)=228.950; volcElev(i)=  500;  volcLoc(i)="Alaska-E                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1105-07-"; volcName(i)="Behm Canal-Rudyerd Bay         "; 
      i=i+1; volcLat(i)= 62.930; volcLon(i)=222.620; volcElev(i)= 1239;  volcLoc(i)="Canada                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1200-01-"; volcName(i)="Fort Selkirk                   "; 
      i=i+1; volcLat(i)= 60.420; volcLon(i)=224.580; volcElev(i)= 2217;  volcLoc(i)="Canada                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1200-02-"; volcName(i)="Alligator Lake                 "; 
      i=i+1; volcLat(i)= 59.680; volcLon(i)=226.680; volcElev(i)= 1880;  volcLoc(i)="Canada                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1200-03-"; volcName(i)="Atlin Volc Field               "; 
      i=i+1; volcLat(i)= 59.370; volcLon(i)=229.420; volcElev(i)= 2123;  volcLoc(i)="Canada                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1200-031"; volcName(i)="Tuya Volc Field                "; 
      i=i+1; volcLat(i)= 58.600; volcLon(i)=228.030; volcElev(i)= 2012;  volcLoc(i)="Canada                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1200-04-"; volcName(i)="Heart Peaks                    "; 
      i=i+1; volcLat(i)= 58.420; volcLon(i)=228.650; volcElev(i)= 2190;  volcLoc(i)="Canada                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1200-05-"; volcName(i)="Level Mountain                 "; 
      i=i+1; volcLat(i)= 57.720; volcLon(i)=229.370; volcElev(i)= 2786;  volcLoc(i)="Canada                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1200-06-"; volcName(i)="Edziza                         "; 
      i=i+1; volcLat(i)= 57.430; volcLon(i)=229.320; volcElev(i)= 2430;  volcLoc(i)="Canada                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1200-07-"; volcName(i)="Spectrum Range                 "; 
      i=i+1; volcLat(i)= 56.780; volcLon(i)=228.720; volcElev(i)= 1850;  volcLoc(i)="Canada                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1200-08-"; volcName(i)="Hoodoo Mountain                "; 
      i=i+1; volcLat(i)= 56.580; volcLon(i)=229.450; volcElev(i)= 1880;  volcLoc(i)="Canada                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1200-09-"; volcName(i)="Iskut-Unuk River Cones         "; 
      i=i+1; volcLat(i)= 55.120; volcLon(i)=231.100; volcElev(i)=  609;  volcLoc(i)="Canada                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1200-10-"; volcName(i)="Tseax River Cone               "; 
      i=i+1; volcLat(i)= 54.700; volcLon(i)=229.770; volcElev(i)=  335;  volcLoc(i)="Canada                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1200-11-"; volcName(i)="Crow Lagoon                    "; 
      i=i+1; volcLat(i)= 52.500; volcLon(i)=231.270; volcElev(i)=  335;  volcLoc(i)="Canada                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1200-12-"; volcName(i)="Milbanke Sound Group           "; 
      i=i+1; volcLat(i)= 52.470; volcLon(i)=235.300; volcElev(i)= 1921;  volcLoc(i)="Canada                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1200-13-"; volcName(i)="Satah Mountain                 "; 
      i=i+1; volcLat(i)= 52.900; volcLon(i)=236.270; volcElev(i)= 1230;  volcLoc(i)="Canada                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1200-14-"; volcName(i)="Nazko                          "; 
      i=i+1; volcLat(i)= 52.330; volcLon(i)=239.430; volcElev(i)= 2015;  volcLoc(i)="Canada                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1200-15-"; volcName(i)="Wells Gray-Clearwater          "; 
      i=i+1; volcLat(i)= 51.430; volcLon(i)=233.700; volcElev(i)= 3160;  volcLoc(i)="Canada                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1200-16-"; volcName(i)="Silverthrone                   "; 
      i=i+1; volcLat(i)= 50.800; volcLon(i)=236.600; volcElev(i)= 2500;  volcLoc(i)="Canada                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1200-17-"; volcName(i)="Bridge River Cones             "; 
      i=i+1; volcLat(i)= 50.630; volcLon(i)=236.500; volcElev(i)= 2680;  volcLoc(i)="Canada                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1200-18-"; volcName(i)="Meager                         "; 
      i=i+1; volcLat(i)= 49.920; volcLon(i)=236.970; volcElev(i)= 2316;  volcLoc(i)="Canada                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1200-19-"; volcName(i)="Garibaldi Lake                 "; 
      i=i+1; volcLat(i)= 49.850; volcLon(i)=237.000; volcElev(i)= 2678;  volcLoc(i)="Canada                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1200-20-"; volcName(i)="Garibaldi                      "; 
      i=i+1; volcLat(i)= 48.777; volcLon(i)=238.187; volcElev(i)= 3285;  volcLoc(i)="US-Washington                 "; 
             volcESP_Code(i)="S1"; volcID(i)="1201-01="; volcName(i)="Baker                          "; 
      i=i+1; volcLat(i)= 48.112; volcLon(i)=238.887; volcElev(i)= 3213;  volcLoc(i)="US-Washington                 "; 
             volcESP_Code(i)="S0"; volcID(i)="1201-02-"; volcName(i)="Glacier Peak                   "; 
      i=i+1; volcLat(i)= 46.853; volcLon(i)=238.240; volcElev(i)= 4392;  volcLoc(i)="US-Washington                 "; 
             volcESP_Code(i)="S1"; volcID(i)="1201-03-"; volcName(i)="Rainier                        "; 
      i=i+1; volcLat(i)= 46.206; volcLon(i)=238.510; volcElev(i)= 3742;  volcLoc(i)="US-Washington                 "; 
             volcESP_Code(i)="S0"; volcID(i)="1201-04-"; volcName(i)="Adams                          "; 
      i=i+1; volcLat(i)= 46.200; volcLon(i)=237.820; volcElev(i)= 2549;  volcLoc(i)="US-Washington                 "; 
             volcESP_Code(i)="S2"; volcID(i)="1201-05-"; volcName(i)="St. Helens                     "; 
      i=i+1; volcLat(i)= 45.880; volcLon(i)=237.920; volcElev(i)= 1329;  volcLoc(i)="US-Washington                 "; 
             volcESP_Code(i)="M0"; volcID(i)="1201-06-"; volcName(i)="West Crater                    "; 
      i=i+1; volcLat(i)= 45.930; volcLon(i)=238.180; volcElev(i)= 1806;  volcLoc(i)="US-Washington                 "; 
             volcESP_Code(i)="M0"; volcID(i)="1201-07-"; volcName(i)="Indian Heaven                  "; 
      i=i+1; volcLat(i)= 45.374; volcLon(i)=238.305; volcElev(i)= 3426;  volcLoc(i)="US-Oregon                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1202-01-"; volcName(i)="Hood                           "; 
      i=i+1; volcLat(i)= 44.674; volcLon(i)=238.200; volcElev(i)= 3199;  volcLoc(i)="US-Oregon                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1202-02-"; volcName(i)="Jefferson                      "; 
      i=i+1; volcLat(i)= 44.411; volcLon(i)=238.226; volcElev(i)= 1230;  volcLoc(i)="US-Oregon                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1202-03-"; volcName(i)="Blue Lake Crater               "; 
      i=i+1; volcLat(i)= 44.380; volcLon(i)=238.070; volcElev(i)= 1664;  volcLoc(i)="US-Oregon                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1202-04-"; volcName(i)="Sand Mountain Field            "; 
      i=i+1; volcLat(i)= 44.285; volcLon(i)=238.159; volcElev(i)= 2095;  volcLoc(i)="US-Oregon                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1202-06-"; volcName(i)="Belknap                        "; 
      i=i+1; volcLat(i)= 44.170; volcLon(i)=238.230; volcElev(i)= 3074;  volcLoc(i)="US-Oregon                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1202-07-"; volcName(i)="North Sister Field             "; 
      i=i+1; volcLat(i)= 44.103; volcLon(i)=238.232; volcElev(i)= 3157;  volcLoc(i)="US-Oregon                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1202-08-"; volcName(i)="South Sister                   "; 
      i=i+1; volcLat(i)= 43.979; volcLon(i)=238.312; volcElev(i)= 2763;  volcLoc(i)="US-Oregon                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1202-09-"; volcName(i)="Bachelor                       "; 
      i=i+1; volcLat(i)= 43.570; volcLon(i)=238.180; volcElev(i)= 2163;  volcLoc(i)="US-Oregon                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1202-10-"; volcName(i)="Davis Lake                     "; 
      i=i+1; volcLat(i)= 43.722; volcLon(i)=238.771; volcElev(i)= 2434;  volcLoc(i)="US-Oregon                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1202-11-"; volcName(i)="Newberry                       "; 
      i=i+1; volcLat(i)= 43.512; volcLon(i)=239.139; volcElev(i)= 1698;  volcLoc(i)="US-Oregon                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1202-12-"; volcName(i)="Devils Garden                  "; 
      i=i+1; volcLat(i)= 43.472; volcLon(i)=239.246; volcElev(i)= 1711;  volcLoc(i)="US-Oregon                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1202-13-"; volcName(i)="Squaw Ridge Lava Field         "; 
      i=i+1; volcLat(i)= 43.361; volcLon(i)=239.331; volcElev(i)= 1501;  volcLoc(i)="US-Oregon                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1202-14-"; volcName(i)="Four Craters Lava Field        "; 
      i=i+1; volcLat(i)= 43.241; volcLon(i)=237.892; volcElev(i)= 1956;  volcLoc(i)="US-Oregon                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1202-15-"; volcName(i)="Cinnamon Butte                 "; 
      i=i+1; volcLat(i)= 42.930; volcLon(i)=237.880; volcElev(i)= 2487;  volcLoc(i)="US-Oregon                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1202-16-"; volcName(i)="Crater Lake                    "; 
      i=i+1; volcLat(i)= 43.100; volcLon(i)=241.250; volcElev(i)= 1435;  volcLoc(i)="US-Oregon                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1202-17-"; volcName(i)="Diamond Craters                "; 
      i=i+1; volcLat(i)= 43.147; volcLon(i)=242.540; volcElev(i)= 1473;  volcLoc(i)="US-Oregon                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1202-19-"; volcName(i)="Jordan Craters                 "; 
      i=i+1; volcLat(i)= 41.409; volcLon(i)=237.807; volcElev(i)= 4317;  volcLoc(i)="US-California                 "; 
             volcESP_Code(i)="S2"; volcID(i)="1203-01-"; volcName(i)="Shasta                         "; 
      i=i+1; volcLat(i)= 41.611; volcLon(i)=238.446; volcElev(i)= 2412;  volcLoc(i)="US-California                 "; 
             volcESP_Code(i)="S0"; volcID(i)="1203-02-"; volcName(i)="Medicine Lake                  "; 
      i=i+1; volcLat(i)= 41.178; volcLon(i)=238.557; volcElev(i)= 1174;  volcLoc(i)="US-California                 "; 
             volcESP_Code(i)="M0"; volcID(i)="1203-03-"; volcName(i)="Brushy Butte                   "; 
      i=i+1; volcLat(i)= 40.777; volcLon(i)=238.409; volcElev(i)= 1631;  volcLoc(i)="US-California                 "; 
             volcESP_Code(i)="M0"; volcID(i)="1203-04-"; volcName(i)="Twin Buttes                    "; 
      i=i+1; volcLat(i)= 40.731; volcLon(i)=238.159; volcElev(i)= 1535;  volcLoc(i)="US-California                 "; 
             volcESP_Code(i)="M0"; volcID(i)="1203-05-"; volcName(i)="Silver Lake                    "; 
      i=i+1; volcLat(i)= 40.680; volcLon(i)=238.450; volcElev(i)= 2191;  volcLoc(i)="US-California                 "; 
             volcESP_Code(i)="M0"; volcID(i)="1203-06-"; volcName(i)="Tumble Buttes                  "; 
      i=i+1; volcLat(i)= 40.492; volcLon(i)=238.492; volcElev(i)= 3187;  volcLoc(i)="US-California                 "; 
             volcESP_Code(i)="S2"; volcID(i)="1203-08-"; volcName(i)="Lassen Volc Center             "; 
      i=i+1; volcLat(i)= 40.630; volcLon(i)=239.170; volcElev(i)= 1652;  volcLoc(i)="US-California                 "; 
             volcESP_Code(i)="M0"; volcID(i)="1203-09-"; volcName(i)="Eagle Lake Field               "; 
      i=i+1; volcLat(i)= 38.970; volcLon(i)=237.230; volcElev(i)= 1439;  volcLoc(i)="US-California                 "; 
             volcESP_Code(i)="S0"; volcID(i)="1203-10-"; volcName(i)="Clear Lake                     "; 
      i=i+1; volcLat(i)= 38.000; volcLon(i)=240.970; volcElev(i)= 2121;  volcLoc(i)="US-California                 "; 
             volcESP_Code(i)="S0"; volcID(i)="1203-11-"; volcName(i)="Mono Lake Volc Field           "; 
      i=i+1; volcLat(i)= 37.880; volcLon(i)=241.000; volcElev(i)= 2796;  volcLoc(i)="US-California                 "; 
             volcESP_Code(i)="S0"; volcID(i)="1203-12-"; volcName(i)="Mono Craters                   "; 
      i=i+1; volcLat(i)= 37.692; volcLon(i)=240.980; volcElev(i)= 2629;  volcLoc(i)="US-California                 "; 
             volcESP_Code(i)="S0"; volcID(i)="1203-13-"; volcName(i)="Inyo Craters                   "; 
      i=i+1; volcLat(i)= 37.631; volcLon(i)=240.968; volcElev(i)= 3369;  volcLoc(i)="US-California                 "; 
             volcESP_Code(i)="S0"; volcID(i)="1203-15-"; volcName(i)="Mammoth Mountain               "; 
      i=i+1; volcLat(i)= 37.020; volcLon(i)=242.550; volcElev(i)=  752;  volcLoc(i)="US-California                 "; 
             volcESP_Code(i)="M0"; volcID(i)="1203-16-"; volcName(i)="Ubehebe Craters                "; 
      i=i+1; volcLat(i)= 36.358; volcLon(i)=241.680; volcElev(i)= 2886;  volcLoc(i)="US-California                 "; 
             volcESP_Code(i)="M0"; volcID(i)="1203-17-"; volcName(i)="Golden Trout Creek             "; 
      i=i+1; volcLat(i)= 36.030; volcLon(i)=242.180; volcElev(i)= 2400;  volcLoc(i)="US-California                 "; 
             volcESP_Code(i)="S0"; volcID(i)="1203-18-"; volcName(i)="Coso Volc Field                "; 
      i=i+1; volcLat(i)= 34.750; volcLon(i)=243.375; volcElev(i)= 1495;  volcLoc(i)="US-California                 "; 
             volcESP_Code(i)="M0"; volcID(i)="1203-19-"; volcName(i)="Lavic Lake                     "; 
      i=i+1; volcLat(i)= 43.180; volcLon(i)=245.650; volcElev(i)= 1478;  volcLoc(i)="US-Idaho                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1204-01-"; volcName(i)="Shoshone Lava Field            "; 
      i=i+1; volcLat(i)= 43.420; volcLon(i)=246.500; volcElev(i)= 2005;  volcLoc(i)="US-Idaho                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1204-02-"; volcName(i)="Craters of the Moon            "; 
      i=i+1; volcLat(i)= 42.880; volcLon(i)=246.780; volcElev(i)= 1604;  volcLoc(i)="US-Idaho                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1204-03-"; volcName(i)="Wapi Lava Field                "; 
      i=i+1; volcLat(i)= 43.500; volcLon(i)=247.550; volcElev(i)= 1631;  volcLoc(i)="US-Idaho                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1204-04-"; volcName(i)="Hell's Half Acre               "; 
      i=i+1; volcLat(i)= 44.430; volcLon(i)=249.330; volcElev(i)= 2805;  volcLoc(i)="US-Wyoming                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1205-01-"; volcName(i)="Yellowstone                    "; 
      i=i+1; volcLat(i)= 39.530; volcLon(i)=241.130; volcElev(i)= 1251;  volcLoc(i)="US-Nevada                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1206-01-"; volcName(i)="Soda Lakes                     "; 
      i=i+1; volcLat(i)= 37.257; volcLon(i)=246.375; volcElev(i)= 1465;  volcLoc(i)="US-Utah                       "; 
             volcESP_Code(i)="M0"; volcID(i)="1207-01-"; volcName(i)="Santa Clara                    "; 
      i=i+1; volcLat(i)= 37.328; volcLon(i)=247.592; volcElev(i)= 2135;  volcLoc(i)="US-Utah                       "; 
             volcESP_Code(i)="M0"; volcID(i)="1207-03-"; volcName(i)="Bald Knoll                     "; 
      i=i+1; volcLat(i)= 37.580; volcLon(i)=247.330; volcElev(i)= 2840;  volcLoc(i)="US-Utah                       "; 
             volcESP_Code(i)="M0"; volcID(i)="1207-04-"; volcName(i)="Markagunt Plateau              "; 
      i=i+1; volcLat(i)= 38.970; volcLon(i)=247.500; volcElev(i)= 1800;  volcLoc(i)="US-Utah                       "; 
             volcESP_Code(i)="M0"; volcID(i)="1207-05-"; volcName(i)="Black Rock Desert              "; 
      i=i+1; volcLat(i)= 39.661; volcLon(i)=252.965; volcElev(i)= 2230;  volcLoc(i)="US-Colorado                   "; 
             volcESP_Code(i)="M0"; volcID(i)="1208-01-"; volcName(i)="Dotsero                        "; 
      i=i+1; volcLat(i)= 36.380; volcLon(i)=246.870; volcElev(i)= 1555;  volcLoc(i)="US-Arizona                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1209-01-"; volcName(i)="Uinkaret Field                 "; 
      i=i+1; volcLat(i)= 35.370; volcLon(i)=248.500; volcElev(i)= 2447;  volcLoc(i)="US-Arizona                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1209-02-"; volcName(i)="Sunset Crater                  "; 
      i=i+1; volcLat(i)= 33.780; volcLon(i)=254.070; volcElev(i)= 1731;  volcLoc(i)="US-New Mexico                 "; 
             volcESP_Code(i)="M0"; volcID(i)="1210-01-"; volcName(i)="Carrizozo                      "; 
      i=i+1; volcLat(i)= 34.800; volcLon(i)=252.000; volcElev(i)= 2550;  volcLoc(i)="US-New Mexico                 "; 
             volcESP_Code(i)="M0"; volcID(i)="1210-02-"; volcName(i)="Zuni-Bandera                   "; 
      i=i+1; volcLat(i)= 47.950; volcLon(i)=230.900; volcElev(i)=-2050;  volcLoc(i)="Pacific-NE                    "; 
             volcESP_Code(i)="U0"; volcID(i)="1301-01-"; volcName(i)="Endeavour Ridge                "; 
      i=i+1; volcLat(i)= 46.880; volcLon(i)=230.670; volcElev(i)=-2100;  volcLoc(i)="Pacific-NE                    "; 
             volcESP_Code(i)="U0"; volcID(i)="1301-011"; volcName(i)="Cobb Segment                   "; 
      i=i+1; volcLat(i)= 46.520; volcLon(i)=230.420; volcElev(i)=-2400;  volcLoc(i)="Pacific-NE                    "; 
             volcESP_Code(i)="U0"; volcID(i)="1301-02-"; volcName(i)="CoAxial Segment                "; 
      i=i+1; volcLat(i)= 45.950; volcLon(i)=230.000; volcElev(i)=-1410;  volcLoc(i)="Pacific-NE                    "; 
             volcESP_Code(i)="U0"; volcID(i)="1301-021"; volcName(i)="Axial Seamount                 "; 
      i=i+1; volcLat(i)= 44.830; volcLon(i)=229.700; volcElev(i)=-2140;  volcLoc(i)="Pacific-NE                    "; 
             volcESP_Code(i)="U0"; volcID(i)="1301-03-"; volcName(i)="Cleft Segment                  "; 
      i=i+1; volcLat(i)= 42.670; volcLon(i)=233.220; volcElev(i)=-3000;  volcLoc(i)="Pacific-NE                    "; 
             volcESP_Code(i)="U0"; volcID(i)="1301-031"; volcName(i)="North Gorda Ridge              "; 
      i=i+1; volcLat(i)= 40.980; volcLon(i)=232.500; volcElev(i)=-1700;  volcLoc(i)="Pacific-NE                    "; 
             volcESP_Code(i)="U0"; volcID(i)="1301-04-"; volcName(i)="Escanaba Segment               "; 
      i=i+1; volcLat(i)= 31.750; volcLon(i)=235.750; volcElev(i)=-2533;  volcLoc(i)="Pacific-NE                    "; 
             volcESP_Code(i)="U0"; volcID(i)="1301-05-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 18.920; volcLon(i)=204.730; volcElev(i)= -975;  volcLoc(i)="Hawaiian Is                   "; 
             volcESP_Code(i)="U0"; volcID(i)="1302-00-"; volcName(i)="Loihi                          "; 
      i=i+1; volcLat(i)= 19.421; volcLon(i)=204.713; volcElev(i)= 1222;  volcLoc(i)="Hawaiian Is                   "; 
             volcESP_Code(i)="M1"; volcID(i)="1302-01-"; volcName(i)="Kilauea                        "; 
      i=i+1; volcLat(i)= 19.475; volcLon(i)=204.392; volcElev(i)= 4170;  volcLoc(i)="Hawaiian Is                   "; 
             volcESP_Code(i)="M1"; volcID(i)="1302-02="; volcName(i)="Mauna Loa                      "; 
      i=i+1; volcLat(i)= 19.820; volcLon(i)=204.530; volcElev(i)= 4205;  volcLoc(i)="Hawaiian Is                   "; 
             volcESP_Code(i)="M0"; volcID(i)="1302-03-"; volcName(i)="Mauna Kea                      "; 
      i=i+1; volcLat(i)= 19.692; volcLon(i)=204.130; volcElev(i)= 2523;  volcLoc(i)="Hawaiian Is                   "; 
             volcESP_Code(i)="M0"; volcID(i)="1302-04-"; volcName(i)="Hualalai                       "; 
      i=i+1; volcLat(i)= 20.708; volcLon(i)=203.750; volcElev(i)= 3055;  volcLoc(i)="Hawaiian Is                   "; 
             volcESP_Code(i)="M0"; volcID(i)="1302-06-"; volcName(i)="Haleakala                      "; 
      i=i+1; volcLat(i)= 21.750; volcLon(i)=201.250; volcElev(i)=-3000;  volcLoc(i)="Hawaiian Is                   "; 
             volcESP_Code(i)="U0"; volcID(i)="1302-08-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 23.580; volcLon(i)=196.170; volcElev(i)=-4000;  volcLoc(i)="Hawaiian Is                   "; 
             volcESP_Code(i)="U0"; volcID(i)="1302-09-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=-17.570; volcLon(i)=211.150; volcElev(i)=-1400;  volcLoc(i)="Society Is-C Pacific          "; 
             volcESP_Code(i)="U0"; volcID(i)="1303-01-"; volcName(i)="Teahitia                       "; 
      i=i+1; volcLat(i)=-17.642; volcLon(i)=211.400; volcElev(i)=-2100;  volcLoc(i)="Society Is-C Pacific          "; 
             volcESP_Code(i)="U0"; volcID(i)="1303-02-"; volcName(i)="Rocard                         "; 
      i=i+1; volcLat(i)=-18.320; volcLon(i)=211.330; volcElev(i)= -160;  volcLoc(i)="Society Is-C Pacific          "; 
             volcESP_Code(i)="U0"; volcID(i)="1303-03-"; volcName(i)="Moua Pihaa                     "; 
      i=i+1; volcLat(i)=-17.870; volcLon(i)=211.930; volcElev(i)=  435;  volcLoc(i)="Society Is-C Pacific          "; 
             volcESP_Code(i)="U0"; volcID(i)="1303-04-"; volcName(i)="Mehetia                        "; 
      i=i+1; volcLat(i)=-25.370; volcLon(i)=230.730; volcElev(i)=  -39;  volcLoc(i)="Pacific-C                     "; 
             volcESP_Code(i)="U0"; volcID(i)="1303-05-"; volcName(i)="Adams Seamount                 "; 
      i=i+1; volcLat(i)=-28.980; volcLon(i)=219.750; volcElev(i)=  -39;  volcLoc(i)="Austral Is-C Pacific          "; 
             volcESP_Code(i)="S0"; volcID(i)="1303-06-"; volcName(i)="Macdonald                      "; 
      i=i+1; volcLat(i)= 16.550; volcLon(i)=254.680; volcElev(i)=-2700;  volcLoc(i)="Pacific-E                     "; 
             volcESP_Code(i)="U0"; volcID(i)="1304-02-"; volcName(i)="Northern EPR-Segment RO2       "; 
      i=i+1; volcLat(i)= 15.830; volcLon(i)=254.570; volcElev(i)=-2300;  volcLoc(i)="Pacific-E                     "; 
             volcESP_Code(i)="U0"; volcID(i)="1304-021"; volcName(i)="Northern EPR-Segment RO3       "; 
      i=i+1; volcLat(i)= 10.730; volcLon(i)=256.420; volcElev(i)=    0;  volcLoc(i)="Pacific-E                     "; 
             volcESP_Code(i)="U0"; volcID(i)="1304-04-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=  9.830; volcLon(i)=255.700; volcElev(i)=-2500;  volcLoc(i)="Pacific-E                     "; 
             volcESP_Code(i)="U0"; volcID(i)="1304-05-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=  0.792; volcLon(i)=273.850; volcElev(i)=-2430;  volcLoc(i)="Pacific-E                     "; 
             volcESP_Code(i)="U0"; volcID(i)="1304-07-"; volcName(i)="Galápagos Rift                 "; 
      i=i+1; volcLat(i)= -8.270; volcLon(i)=252.050; volcElev(i)=-2800;  volcLoc(i)="Pacific-E                     "; 
             volcESP_Code(i)="U0"; volcID(i)="1304-10-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=-17.436; volcLon(i)=246.794; volcElev(i)=-2566;  volcLoc(i)="Pacific-E                     "; 
             volcESP_Code(i)="U0"; volcID(i)="1304-12-"; volcName(i)="Southern EPR-Segment K         "; 
      i=i+1; volcLat(i)=-18.175; volcLon(i)=246.650; volcElev(i)=-2650;  volcLoc(i)="Pacific-E                     "; 
             volcESP_Code(i)="U0"; volcID(i)="1304-13-"; volcName(i)="Southern EPR-Segment J         "; 
      i=i+1; volcLat(i)=-18.530; volcLon(i)=246.580; volcElev(i)=-2600;  volcLoc(i)="Pacific-E                     "; 
             volcESP_Code(i)="U0"; volcID(i)="1304-14-"; volcName(i)="Southern EPR-Segment I         "; 
      i=i+1; volcLat(i)=-49.680; volcLon(i)=178.770; volcElev(i)=  402;  volcLoc(i)="Pacific-S                     "; 
             volcESP_Code(i)="U0"; volcID(i)="1305-01-"; volcName(i)="Antipodes Island               "; 
      i=i+1; volcLat(i)=-53.900; volcLon(i)=219.700; volcElev(i)=-1000;  volcLoc(i)="Pacific-S                     "; 
             volcESP_Code(i)="U0"; volcID(i)="1305-02-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=-55.970; volcLon(i)=216.830; volcElev(i)=    0;  volcLoc(i)="Pacific-S                     "; 
             volcESP_Code(i)="U0"; volcID(i)="1305-03-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 32.418; volcLon(i)=244.695; volcElev(i)=  223;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1401-00-"; volcName(i)="Prieto, Cerro                  "; 
      i=i+1; volcLat(i)= 31.772; volcLon(i)=246.502; volcElev(i)= 1200;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1401-001"; volcName(i)="Pinacate                       "; 
      i=i+1; volcLat(i)= 30.468; volcLon(i)=244.004; volcElev(i)=  260;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1401-002"; volcName(i)="San Quintín Volc Field         "; 
      i=i+1; volcLat(i)= 29.970; volcLon(i)=245.600; volcElev(i)=  180;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1401-003"; volcName(i)="San Luis, Isla                 "; 
      i=i+1; volcLat(i)= 29.330; volcLon(i)=245.500; volcElev(i)=  960;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1401-004"; volcName(i)="Jaraguay Volc Field            "; 
      i=i+1; volcLat(i)= 29.080; volcLon(i)=246.487; volcElev(i)=  440;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1401-005"; volcName(i)="Coronado                       "; 
      i=i+1; volcLat(i)= 29.070; volcLon(i)=241.720; volcElev(i)= 1100;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1401-006"; volcName(i)="Guadalupe                      "; 
      i=i+1; volcLat(i)= 28.500; volcLon(i)=246.250; volcElev(i)= 1360;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1401-007"; volcName(i)="San Borja Volc Field           "; 
      i=i+1; volcLat(i)= 28.000; volcLon(i)=245.000; volcElev(i)=    0;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="U0"; volcID(i)="1401-008"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 27.430; volcLon(i)=248.120; volcElev(i)=  210;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1401-011"; volcName(i)="Tortuga, Isla                  "; 
      i=i+1; volcLat(i)= 26.000; volcLon(i)=248.080; volcElev(i)=  780;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1401-012"; volcName(i)="Comondú-La Purísima            "; 
      i=i+1; volcLat(i)= 27.470; volcLon(i)=247.409; volcElev(i)= 1940;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1401-01="; volcName(i)="Tres Vírgenes                  "; 
      i=i+1; volcLat(i)= 18.780; volcLon(i)=249.050; volcElev(i)= 1050;  volcLoc(i)="Mexico-Is                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1401-021"; volcName(i)="Socorro                        "; 
      i=i+1; volcLat(i)= 24.150; volcLon(i)=255.550; volcElev(i)= 2075;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1401-022"; volcName(i)="Durango Volc Field             "; 
      i=i+1; volcLat(i)= 21.450; volcLon(i)=255.270; volcElev(i)= 2340;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1401-023"; volcName(i)="Sangangüey                     "; 
      i=i+1; volcLat(i)= 19.300; volcLon(i)=249.180; volcElev(i)=  332;  volcLoc(i)="Mexico-Is                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1401-02="; volcName(i)="Bárcena                        "; 
      i=i+1; volcLat(i)= 20.620; volcLon(i)=255.170; volcElev(i)= 2560;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1401-031"; volcName(i)="Mascota Volc Field             "; 
      i=i+1; volcLat(i)= 21.125; volcLon(i)=255.492; volcElev(i)= 2280;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="S2"; volcID(i)="1401-03="; volcName(i)="Ceboruco                       "; 
      i=i+1; volcLat(i)= 19.514; volcLon(i)=256.380; volcElev(i)= 3850;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1401-04="; volcName(i)="Colima                         "; 
      i=i+1; volcLat(i)= 19.400; volcLon(i)=259.750; volcElev(i)= 3500;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1401-061"; volcName(i)="Zitácuaro-Valle de Bravo       "; 
      i=i+1; volcLat(i)= 19.730; volcLon(i)=260.242; volcElev(i)= 3900;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1401-062"; volcName(i)="Jocotitlán                     "; 
      i=i+1; volcLat(i)= 19.850; volcLon(i)=258.250; volcElev(i)= 3860;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="M2"; volcID(i)="1401-06="; volcName(i)="Michoacán-Guanajuato           "; 
      i=i+1; volcLat(i)= 19.108; volcLon(i)=260.242; volcElev(i)= 4680;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="S2"; volcID(i)="1401-07-"; volcName(i)="Toluca, Nevado de              "; 
      i=i+1; volcLat(i)= 19.308; volcLon(i)=261.300; volcElev(i)= 3600;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1401-081"; volcName(i)="Papayo                         "; 
      i=i+1; volcLat(i)= 19.179; volcLon(i)=261.358; volcElev(i)= 5230;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1401-082"; volcName(i)="Iztaccíhuatl                   "; 
      i=i+1; volcLat(i)= 19.080; volcLon(i)=260.870; volcElev(i)= 3930;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1401-08="; volcName(i)="Chichinautzin                  "; 
      i=i+1; volcLat(i)= 19.231; volcLon(i)=261.968; volcElev(i)= 4461;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1401-091"; volcName(i)="Malinche, La                   "; 
      i=i+1; volcLat(i)= 19.270; volcLon(i)=262.530; volcElev(i)= 3485;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1401-092"; volcName(i)="Serdán-Oriental                "; 
      i=i+1; volcLat(i)= 19.680; volcLon(i)=262.550; volcElev(i)= 3150;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1401-093"; volcName(i)="Humeros, Los                   "; 
      i=i+1; volcLat(i)= 19.809; volcLon(i)=263.474; volcElev(i)=  800;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1401-094"; volcName(i)="Atlixcos, Los                  "; 
      i=i+1; volcLat(i)= 19.670; volcLon(i)=263.250; volcElev(i)= 2000;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1401-095"; volcName(i)="Naolinco Volc Field            "; 
      i=i+1; volcLat(i)= 19.492; volcLon(i)=262.850; volcElev(i)= 4282;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1401-096"; volcName(i)="Cofre de Perote                "; 
      i=i+1; volcLat(i)= 19.330; volcLon(i)=262.750; volcElev(i)= 3500;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1401-097"; volcName(i)="Gloria, La                     "; 
      i=i+1; volcLat(i)= 19.150; volcLon(i)=262.730; volcElev(i)= 3940;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1401-098"; volcName(i)="Cumbres,  Las                  "; 
      i=i+1; volcLat(i)= 19.023; volcLon(i)=261.378; volcElev(i)= 5426;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="S1"; volcID(i)="1401-09="; volcName(i)="Popocatépetl                   "; 
      i=i+1; volcLat(i)= 19.030; volcLon(i)=262.732; volcElev(i)= 5675;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1401-10="; volcName(i)="Orizaba, Pico de               "; 
      i=i+1; volcLat(i)= 18.570; volcLon(i)=264.800; volcElev(i)= 1650;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1401-11="; volcName(i)="San Martín                     "; 
      i=i+1; volcLat(i)= 17.360; volcLon(i)=266.772; volcElev(i)= 1150;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1401-12="; volcName(i)="Chichón, El                    "; 
      i=i+1; volcLat(i)= 15.130; volcLon(i)=267.888; volcElev(i)= 4060;  volcLoc(i)="Mexico                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1401-13="; volcName(i)="Tacaná                         "; 
      i=i+1; volcLat(i)= 15.034; volcLon(i)=268.097; volcElev(i)= 4220;  volcLoc(i)="Guatemala                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1402-02="; volcName(i)="Tajumulco                      "; 
      i=i+1; volcLat(i)= 14.756; volcLon(i)=268.448; volcElev(i)= 3772;  volcLoc(i)="Guatemala                     "; 
             volcESP_Code(i)="S9"; volcID(i)="1402-03="; volcName(i)="Santa María                    "; 
      i=i+1; volcLat(i)= 14.820; volcLon(i)=268.520; volcElev(i)= 3197;  volcLoc(i)="Guatemala                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1402-04="; volcName(i)="Almolonga                      "; 
      i=i+1; volcLat(i)= 14.583; volcLon(i)=268.814; volcElev(i)= 3535;  volcLoc(i)="Guatemala                     "; 
             volcESP_Code(i)="S1"; volcID(i)="1402-06="; volcName(i)="Atitlán                        "; 
      i=i+1; volcLat(i)= 14.612; volcLon(i)=268.811; volcElev(i)= 3158;  volcLoc(i)="Guatemala                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1402-07="; volcName(i)="Tolimán                        "; 
      i=i+1; volcLat(i)= 14.501; volcLon(i)=269.124; volcElev(i)= 3976;  volcLoc(i)="Guatemala                     "; 
             volcESP_Code(i)="S1"; volcID(i)="1402-08="; volcName(i)="Acatenango                     "; 
      i=i+1; volcLat(i)= 14.473; volcLon(i)=269.120; volcElev(i)= 3763;  volcLoc(i)="Guatemala                     "; 
             volcESP_Code(i)="M2"; volcID(i)="1402-09="; volcName(i)="Fuego                          "; 
      i=i+1; volcLat(i)= 14.465; volcLon(i)=269.257; volcElev(i)= 3760;  volcLoc(i)="Guatemala                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1402-10="; volcName(i)="Agua                           "; 
      i=i+1; volcLat(i)= 14.330; volcLon(i)=269.600; volcElev(i)= 1454;  volcLoc(i)="Guatemala                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1402-111"; volcName(i)="Cuilapa-Barbarena              "; 
      i=i+1; volcLat(i)= 14.381; volcLon(i)=269.399; volcElev(i)= 2552;  volcLoc(i)="Guatemala                     "; 
             volcESP_Code(i)="M1"; volcID(i)="1402-11="; volcName(i)="Pacaya                         "; 
      i=i+1; volcLat(i)= 14.336; volcLon(i)=269.731; volcElev(i)= 1815;  volcLoc(i)="Guatemala                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1402-121"; volcName(i)="Jumaytepeque                   "; 
      i=i+1; volcLat(i)= 14.156; volcLon(i)=269.593; volcElev(i)= 1845;  volcLoc(i)="Guatemala                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1402-12="; volcName(i)="Tecuamburro                    "; 
      i=i+1; volcLat(i)= 14.030; volcLon(i)=269.900; volcElev(i)= 1662;  volcLoc(i)="Guatemala                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1402-13-"; volcName(i)="Moyuta                         "; 
      i=i+1; volcLat(i)= 14.308; volcLon(i)=270.008; volcElev(i)= 1600;  volcLoc(i)="Guatemala                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1402-14-"; volcName(i)="Flores                         "; 
      i=i+1; volcLat(i)= 14.430; volcLon(i)=270.100; volcElev(i)= 1716;  volcLoc(i)="Guatemala                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1402-141"; volcName(i)="Tahual                         "; 
      i=i+1; volcLat(i)= 14.330; volcLon(i)=270.130; volcElev(i)= 1192;  volcLoc(i)="Guatemala                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1402-15-"; volcName(i)="Santiago, Cerro                "; 
      i=i+1; volcLat(i)= 14.400; volcLon(i)=270.220; volcElev(i)= 2042;  volcLoc(i)="Guatemala                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1402-16-"; volcName(i)="Suchitán                       "; 
      i=i+1; volcLat(i)= 14.120; volcLon(i)=270.270; volcElev(i)= 1775;  volcLoc(i)="Guatemala                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1402-17-"; volcName(i)="Chingo                         "; 
      i=i+1; volcLat(i)= 14.420; volcLon(i)=270.320; volcElev(i)= 1292;  volcLoc(i)="Guatemala                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1402-18-"; volcName(i)="Ixtepeque                      "; 
      i=i+1; volcLat(i)= 14.550; volcLon(i)=270.370; volcElev(i)= 1650;  volcLoc(i)="Guatemala                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1402-19-"; volcName(i)="Ipala                          "; 
      i=i+1; volcLat(i)= 14.830; volcLon(i)=270.450; volcElev(i)= 1192;  volcLoc(i)="Guatemala                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1402-20-"; volcName(i)="Chiquimula Volc Field          "; 
      i=i+1; volcLat(i)= 14.570; volcLon(i)=270.550; volcElev(i)= 1200;  volcLoc(i)="Guatemala                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1402-21-"; volcName(i)="Quezaltepeque                  "; 
      i=i+1; volcLat(i)= 14.270; volcLon(i)=270.520; volcElev(i)=  781;  volcLoc(i)="El Salvador                   "; 
             volcESP_Code(i)="M0"; volcID(i)="1403-001"; volcName(i)="San Diego                      "; 
      i=i+1; volcLat(i)= 14.050; volcLon(i)=270.350; volcElev(i)=  957;  volcLoc(i)="El Salvador                   "; 
             volcESP_Code(i)="M0"; volcID(i)="1403-002"; volcName(i)="Singüil, Cerro                 "; 
      i=i+1; volcLat(i)= 13.891; volcLon(i)=270.214; volcElev(i)= 2036;  volcLoc(i)="El Salvador                   "; 
             volcESP_Code(i)="M0"; volcID(i)="1403-01="; volcName(i)="Apaneca Range                  "; 
      i=i+1; volcLat(i)= 13.853; volcLon(i)=270.370; volcElev(i)= 2381;  volcLoc(i)="El Salvador                   "; 
             volcESP_Code(i)="M2"; volcID(i)="1403-02="; volcName(i)="Santa Ana                      "; 
      i=i+1; volcLat(i)= 13.813; volcLon(i)=270.367; volcElev(i)= 1950;  volcLoc(i)="El Salvador                   "; 
             volcESP_Code(i)="M1"; volcID(i)="1403-03="; volcName(i)="Izalco                         "; 
      i=i+1; volcLat(i)= 13.870; volcLon(i)=270.450; volcElev(i)=  746;  volcLoc(i)="El Salvador                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1403-041"; volcName(i)="Coatepeque Caldera             "; 
      i=i+1; volcLat(i)= 14.020; volcLon(i)=270.750; volcElev(i)=  665;  volcLoc(i)="El Salvador                   "; 
             volcESP_Code(i)="M0"; volcID(i)="1403-051"; volcName(i)="Cinotepeque, Cerro             "; 
      i=i+1; volcLat(i)= 13.900; volcLon(i)=270.880; volcElev(i)= 1438;  volcLoc(i)="El Salvador                   "; 
             volcESP_Code(i)="M0"; volcID(i)="1403-052"; volcName(i)="Guazapa                        "; 
      i=i+1; volcLat(i)= 13.734; volcLon(i)=270.706; volcElev(i)= 1893;  volcLoc(i)="El Salvador                   "; 
             volcESP_Code(i)="M2"; volcID(i)="1403-05="; volcName(i)="San Salvador                   "; 
      i=i+1; volcLat(i)= 13.672; volcLon(i)=270.947; volcElev(i)=  450;  volcLoc(i)="El Salvador                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1403-06="; volcName(i)="Ilopango                       "; 
      i=i+1; volcLat(i)= 13.720; volcLon(i)=271.230; volcElev(i)=  700;  volcLoc(i)="El Salvador                   "; 
             volcESP_Code(i)="M0"; volcID(i)="1403-071"; volcName(i)="Apastepeque Field              "; 
      i=i+1; volcLat(i)= 13.435; volcLon(i)=271.468; volcElev(i)= 1172;  volcLoc(i)="El Salvador                   "; 
             volcESP_Code(i)="M0"; volcID(i)="1403-072"; volcName(i)="Taburete                       "; 
      i=i+1; volcLat(i)= 13.595; volcLon(i)=271.163; volcElev(i)= 2182;  volcLoc(i)="El Salvador                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1403-07="; volcName(i)="San Vicente                    "; 
      i=i+1; volcLat(i)= 13.419; volcLon(i)=271.529; volcElev(i)= 1449;  volcLoc(i)="El Salvador                   "; 
             volcESP_Code(i)="M0"; volcID(i)="1403-081"; volcName(i)="Usulután                       "; 
      i=i+1; volcLat(i)= 13.470; volcLon(i)=271.570; volcElev(i)= 1640;  volcLoc(i)="El Salvador                   "; 
             volcESP_Code(i)="M0"; volcID(i)="1403-082"; volcName(i)="Tigre, El                      "; 
      i=i+1; volcLat(i)= 13.494; volcLon(i)=271.498; volcElev(i)= 1593;  volcLoc(i)="El Salvador                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1403-08="; volcName(i)="Tecapa                         "; 
      i=i+1; volcLat(i)= 13.478; volcLon(i)=271.670; volcElev(i)= 1300;  volcLoc(i)="El Salvador                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1403-09="; volcName(i)="Chinameca                      "; 
      i=i+1; volcLat(i)= 13.428; volcLon(i)=271.895; volcElev(i)=  181;  volcLoc(i)="El Salvador                   "; 
             volcESP_Code(i)="M0"; volcID(i)="1403-101"; volcName(i)="Aramuaca, Laguna               "; 
      i=i+1; volcLat(i)= 13.434; volcLon(i)=271.731; volcElev(i)= 2130;  volcLoc(i)="El Salvador                   "; 
             volcESP_Code(i)="S1"; volcID(i)="1403-10="; volcName(i)="San Miguel                     "; 
      i=i+1; volcLat(i)= 13.275; volcLon(i)=272.155; volcElev(i)= 1225;  volcLoc(i)="El Salvador                   "; 
             volcESP_Code(i)="S0"; volcID(i)="1403-11="; volcName(i)="Conchagua                      "; 
      i=i+1; volcLat(i)= 13.229; volcLon(i)=272.233; volcElev(i)=  505;  volcLoc(i)="El Salvador                   "; 
             volcESP_Code(i)="M0"; volcID(i)="1403-12="; volcName(i)="Conchagüita                    "; 
      i=i+1; volcLat(i)= 13.272; volcLon(i)=272.359; volcElev(i)=  783;  volcLoc(i)="Honduras                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1403-13-"; volcName(i)="Tigre, Isla el                 "; 
      i=i+1; volcLat(i)= 13.330; volcLon(i)=272.370; volcElev(i)=  640;  volcLoc(i)="Honduras                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1403-14-"; volcName(i)="Zacate Grande, Isla            "; 
      i=i+1; volcLat(i)= 14.980; volcLon(i)=272.020; volcElev(i)= 1090;  volcLoc(i)="Honduras                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1403-15-"; volcName(i)="Yojoa, Lake                    "; 
      i=i+1; volcLat(i)= 16.100; volcLon(i)=273.100; volcElev(i)=   74;  volcLoc(i)="Honduras                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1403-16-"; volcName(i)="Utila Island                   "; 
      i=i+1; volcLat(i)= 12.980; volcLon(i)=272.430; volcElev(i)=  872;  volcLoc(i)="Nicaragua                     "; 
             volcESP_Code(i)="S2"; volcID(i)="1404-01="; volcName(i)="Cosigüina                      "; 
      i=i+1; volcLat(i)= 12.702; volcLon(i)=272.996; volcElev(i)= 1745;  volcLoc(i)="Nicaragua                     "; 
             volcESP_Code(i)="S1"; volcID(i)="1404-02="; volcName(i)="San Cristóbal                  "; 
      i=i+1; volcLat(i)= 12.602; volcLon(i)=273.155; volcElev(i)= 1061;  volcLoc(i)="Nicaragua                     "; 
             volcESP_Code(i)="S1"; volcID(i)="1404-04="; volcName(i)="Telica                         "; 
      i=i+1; volcLat(i)= 12.550; volcLon(i)=273.250; volcElev(i)=  832;  volcLoc(i)="Nicaragua                     "; 
             volcESP_Code(i)="S1"; volcID(i)="1404-06-"; volcName(i)="Rota                           "; 
      i=i+1; volcLat(i)= 12.506; volcLon(i)=273.298; volcElev(i)=  728;  volcLoc(i)="Nicaragua                     "; 
             volcESP_Code(i)="M2"; volcID(i)="1404-07="; volcName(i)="Negro, Cerro                   "; 
      i=i+1; volcLat(i)= 12.495; volcLon(i)=273.312; volcElev(i)= 1088;  volcLoc(i)="Nicaragua                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1404-08="; volcName(i)="Pilas, Las                     "; 
      i=i+1; volcLat(i)= 12.242; volcLon(i)=273.658; volcElev(i)=  518;  volcLoc(i)="Nicaragua                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1404-091"; volcName(i)="Apoyeque                       "; 
      i=i+1; volcLat(i)= 12.120; volcLon(i)=273.680; volcElev(i)=  360;  volcLoc(i)="Nicaragua                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1404-092"; volcName(i)="Nejapa-Miraflores              "; 
      i=i+1; volcLat(i)= 12.422; volcLon(i)=273.460; volcElev(i)= 1297;  volcLoc(i)="Nicaragua                     "; 
             volcESP_Code(i)="M2"; volcID(i)="1404-09="; volcName(i)="Momotombo                      "; 
      i=i+1; volcLat(i)= 11.920; volcLon(i)=274.020; volcElev(i)=  300;  volcLoc(i)="Nicaragua                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1404-101"; volcName(i)="Granada                        "; 
      i=i+1; volcLat(i)= 11.984; volcLon(i)=273.839; volcElev(i)=  635;  volcLoc(i)="Nicaragua                     "; 
             volcESP_Code(i)="M1"; volcID(i)="1404-10="; volcName(i)="Masaya                         "; 
      i=i+1; volcLat(i)= 11.730; volcLon(i)=274.180; volcElev(i)=  629;  volcLoc(i)="Nicaragua                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1404-111"; volcName(i)="Zapatera                       "; 
      i=i+1; volcLat(i)= 11.826; volcLon(i)=274.032; volcElev(i)= 1344;  volcLoc(i)="Nicaragua                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1404-11="; volcName(i)="Mombacho                       "; 
      i=i+1; volcLat(i)= 11.538; volcLon(i)=274.378; volcElev(i)= 1700;  volcLoc(i)="Nicaragua                     "; 
             volcESP_Code(i)="S1"; volcID(i)="1404-12="; volcName(i)="Concepción                     "; 
      i=i+1; volcLat(i)= 11.446; volcLon(i)=274.485; volcElev(i)= 1394;  volcLoc(i)="Nicaragua                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1404-13-"; volcName(i)="Maderas                        "; 
      i=i+1; volcLat(i)= 13.170; volcLon(i)=273.600; volcElev(i)=  899;  volcLoc(i)="Nicaragua                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1404-131"; volcName(i)="Estelí                         "; 
      i=i+1; volcLat(i)= 12.530; volcLon(i)=273.858; volcElev(i)=  603;  volcLoc(i)="Nicaragua                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1404-132"; volcName(i)="Ciguatepe, Cerro el            "; 
      i=i+1; volcLat(i)= 12.300; volcLon(i)=274.270; volcElev(i)=  926;  volcLoc(i)="Nicaragua                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1404-133"; volcName(i)="Lajas, Las                     "; 
      i=i+1; volcLat(i)= 12.530; volcLon(i)=276.130; volcElev(i)=  201;  volcLoc(i)="Nicaragua                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1404-14-"; volcName(i)="Azul, Volcán                   "; 
      i=i+1; volcLat(i)= 10.980; volcLon(i)=274.527; volcElev(i)= 1659;  volcLoc(i)="Costa Rica                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1405-01="; volcName(i)="Orosí                          "; 
      i=i+1; volcLat(i)= 10.830; volcLon(i)=274.676; volcElev(i)= 1916;  volcLoc(i)="Costa Rica                    "; 
             volcESP_Code(i)="S2"; volcID(i)="1405-02="; volcName(i)="Rincón de la Vieja             "; 
      i=i+1; volcLat(i)= 10.673; volcLon(i)=274.985; volcElev(i)= 1916;  volcLoc(i)="Costa Rica                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1405-031"; volcName(i)="Tenorio                        "; 
      i=i+1; volcLat(i)= 10.463; volcLon(i)=275.297; volcElev(i)= 1670;  volcLoc(i)="Costa Rica                    "; 
             volcESP_Code(i)="S1"; volcID(i)="1405-033"; volcName(i)="Arenal                         "; 
      i=i+1; volcLat(i)= 10.300; volcLon(i)=275.634; volcElev(i)= 2267;  volcLoc(i)="Costa Rica                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1405-034"; volcName(i)="Platanar                       "; 
      i=i+1; volcLat(i)= 10.748; volcLon(i)=274.847; volcElev(i)= 2028;  volcLoc(i)="Costa Rica                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1405-03="; volcName(i)="Miravalles                     "; 
      i=i+1; volcLat(i)= 10.200; volcLon(i)=275.767; volcElev(i)= 2708;  volcLoc(i)="Costa Rica                    "; 
             volcESP_Code(i)="S1"; volcID(i)="1405-04="; volcName(i)="Poás                           "; 
      i=i+1; volcLat(i)= 10.135; volcLon(i)=275.900; volcElev(i)= 2906;  volcLoc(i)="Costa Rica                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1405-05="; volcName(i)="Barva                          "; 
      i=i+1; volcLat(i)=  9.979; volcLon(i)=276.148; volcElev(i)= 3432;  volcLoc(i)="Costa Rica                    "; 
             volcESP_Code(i)="S1"; volcID(i)="1405-06="; volcName(i)="Irazú                          "; 
      i=i+1; volcLat(i)= 10.025; volcLon(i)=276.233; volcElev(i)= 3340;  volcLoc(i)="Costa Rica                    "; 
             volcESP_Code(i)="S1"; volcID(i)="1405-07="; volcName(i)="Turrialba                      "; 
      i=i+1; volcLat(i)=  8.808; volcLon(i)=277.457; volcElev(i)= 3474;  volcLoc(i)="Panama                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1406-01-"; volcName(i)="Barú                           "; 
      i=i+1; volcLat(i)=  8.470; volcLon(i)=279.180; volcElev(i)= 1297;  volcLoc(i)="Panama                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1406-02-"; volcName(i)="Yeguada, La                    "; 
      i=i+1; volcLat(i)=  8.580; volcLon(i)=279.830; volcElev(i)= 1185;  volcLoc(i)="Panama                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1406-03-"; volcName(i)="Valle, El                      "; 
      i=i+1; volcLat(i)=  5.206; volcLon(i)=284.636; volcElev(i)= 3858;  volcLoc(i)="Colombia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1501-011"; volcName(i)="Romeral                        "; 
      i=i+1; volcLat(i)=  5.092; volcLon(i)=284.700; volcElev(i)= 4000;  volcLoc(i)="Colombia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1501-012"; volcName(i)="Bravo, Cerro                   "; 
      i=i+1; volcLat(i)=  4.820; volcLon(i)=284.630; volcElev(i)= 4950;  volcLoc(i)="Colombia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1501-021"; volcName(i)="Santa Isabel                   "; 
      i=i+1; volcLat(i)=  4.895; volcLon(i)=284.678; volcElev(i)= 5321;  volcLoc(i)="Colombia                      "; 
             volcESP_Code(i)="S1"; volcID(i)="1501-02="; volcName(i)="Ruiz, Nevado del               "; 
      i=i+1; volcLat(i)=  4.670; volcLon(i)=284.670; volcElev(i)= 5200;  volcLoc(i)="Colombia                      "; 
             volcESP_Code(i)="S2"; volcID(i)="1501-03="; volcName(i)="Tolima, Nevado del             "; 
      i=i+1; volcLat(i)=  4.480; volcLon(i)=284.608; volcElev(i)= 2650;  volcLoc(i)="Colombia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1501-04="; volcName(i)="Machín                         "; 
      i=i+1; volcLat(i)=  2.930; volcLon(i)=283.970; volcElev(i)= 5364;  volcLoc(i)="Colombia                      "; 
             volcESP_Code(i)="S1"; volcID(i)="1501-05="; volcName(i)="Huila, Nevado del              "; 
      i=i+1; volcLat(i)=  2.108; volcLon(i)=283.408; volcElev(i)= 4400;  volcLoc(i)="Colombia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1501-061"; volcName(i)="Sotará                         "; 
      i=i+1; volcLat(i)=  1.570; volcLon(i)=283.220; volcElev(i)= 4054;  volcLoc(i)="Colombia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1501-062"; volcName(i)="Petacas                        "; 
      i=i+1; volcLat(i)=  2.320; volcLon(i)=283.600; volcElev(i)= 4650;  volcLoc(i)="Colombia                      "; 
             volcESP_Code(i)="S1"; volcID(i)="1501-06="; volcName(i)="Puracé                         "; 
      i=i+1; volcLat(i)=  1.470; volcLon(i)=283.080; volcElev(i)= 4150;  volcLoc(i)="Colombia                      "; 
             volcESP_Code(i)="S2"; volcID(i)="1501-07="; volcName(i)="Doña Juana                     "; 
      i=i+1; volcLat(i)=  1.220; volcLon(i)=282.630; volcElev(i)= 4276;  volcLoc(i)="Colombia                      "; 
             volcESP_Code(i)="S1"; volcID(i)="1501-08="; volcName(i)="Galeras                        "; 
      i=i+1; volcLat(i)=  1.080; volcLon(i)=282.320; volcElev(i)= 4070;  volcLoc(i)="Colombia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1501-09="; volcName(i)="Azufral                        "; 
      i=i+1; volcLat(i)=  0.950; volcLon(i)=282.130; volcElev(i)= 4764;  volcLoc(i)="Colombia                      "; 
             volcESP_Code(i)="S1"; volcID(i)="1501-10="; volcName(i)="Cumbal                         "; 
      i=i+1; volcLat(i)=  0.828; volcLon(i)=282.036; volcElev(i)= 4445;  volcLoc(i)="Colombia                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1501-11="; volcName(i)="Negro de Mayasquer, Cerro      "; 
      i=i+1; volcLat(i)=  0.552; volcLon(i)=282.420; volcElev(i)= 3955;  volcLoc(i)="Ecuador                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1502-001"; volcName(i)="Soche                          "; 
      i=i+1; volcLat(i)=  0.468; volcLon(i)=281.713; volcElev(i)= 4106;  volcLoc(i)="Ecuador                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1502-002"; volcName(i)="Chachimbiro                    "; 
      i=i+1; volcLat(i)=  0.308; volcLon(i)=281.636; volcElev(i)= 3246;  volcLoc(i)="Ecuador                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1502-003"; volcName(i)="Cuicocha                       "; 
      i=i+1; volcLat(i)=  0.258; volcLon(i)=281.817; volcElev(i)= 4609;  volcLoc(i)="Ecuador                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1502-004"; volcName(i)="Imbabura                       "; 
      i=i+1; volcLat(i)=  0.130; volcLon(i)=281.730; volcElev(i)= 4263;  volcLoc(i)="Ecuador                       "; 
             volcESP_Code(i)="S2"; volcID(i)="1502-005"; volcName(i)="Mojanda                        "; 
      i=i+1; volcLat(i)=  0.029; volcLon(i)=282.014; volcElev(i)= 5790;  volcLoc(i)="Ecuador                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1502-006"; volcName(i)="Cayambe                        "; 
      i=i+1; volcLat(i)=  0.038; volcLon(i)=281.537; volcElev(i)= 3356;  volcLoc(i)="Ecuador                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1502-011"; volcName(i)="Pululagua                      "; 
      i=i+1; volcLat(i)= -0.077; volcLon(i)=282.344; volcElev(i)= 3562;  volcLoc(i)="Ecuador                       "; 
             volcESP_Code(i)="S2"; volcID(i)="1502-01="; volcName(i)="Reventador                     "; 
      i=i+1; volcLat(i)= -0.353; volcLon(i)=281.383; volcElev(i)= 4463;  volcLoc(i)="Ecuador                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1502-021"; volcName(i)="Atacazo                        "; 
      i=i+1; volcLat(i)= -0.375; volcLon(i)=281.750; volcElev(i)= 4643;  volcLoc(i)="Ecuador                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1502-022"; volcName(i)="Chacana                        "; 
      i=i+1; volcLat(i)= -0.171; volcLon(i)=281.402; volcElev(i)= 4784;  volcLoc(i)="Ecuador                       "; 
             volcESP_Code(i)="S1"; volcID(i)="1502-02="; volcName(i)="Guagua Pichincha               "; 
      i=i+1; volcLat(i)= -0.481; volcLon(i)=281.859; volcElev(i)= 5753;  volcLoc(i)="Ecuador                       "; 
             volcESP_Code(i)="S2"; volcID(i)="1502-03="; volcName(i)="Antisana                       "; 
      i=i+1; volcLat(i)= -0.659; volcLon(i)=281.286; volcElev(i)= 5248;  volcLoc(i)="Ecuador                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1502-041"; volcName(i)="Illiniza                       "; 
      i=i+1; volcLat(i)= -0.538; volcLon(i)=282.374; volcElev(i)= 3990;  volcLoc(i)="Ecuador                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1502-04="; volcName(i)="Sumaco                         "; 
      i=i+1; volcLat(i)= -0.677; volcLon(i)=281.564; volcElev(i)= 5911;  volcLoc(i)="Ecuador                       "; 
             volcESP_Code(i)="S2"; volcID(i)="1502-05="; volcName(i)="Cotopaxi                       "; 
      i=i+1; volcLat(i)= -0.850; volcLon(i)=281.100; volcElev(i)= 3914;  volcLoc(i)="Ecuador                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1502-06="; volcName(i)="Quilotoa                       "; 
      i=i+1; volcLat(i)= -1.464; volcLon(i)=281.185; volcElev(i)= 6310;  volcLoc(i)="Ecuador                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1502-071"; volcName(i)="Chimborazo                     "; 
      i=i+1; volcLat(i)= -1.780; volcLon(i)=281.387; volcElev(i)= 3336;  volcLoc(i)="Ecuador                       "; 
             volcESP_Code(i)="M0"; volcID(i)="1502-081"; volcName(i)="Licto                          "; 
      i=i+1; volcLat(i)= -1.467; volcLon(i)=281.558; volcElev(i)= 5023;  volcLoc(i)="Ecuador                       "; 
             volcESP_Code(i)="S1"; volcID(i)="1502-08="; volcName(i)="Tungurahua                     "; 
      i=i+1; volcLat(i)= -2.002; volcLon(i)=281.659; volcElev(i)= 5230;  volcLoc(i)="Ecuador                       "; 
             volcESP_Code(i)="S9"; volcID(i)="1502-09="; volcName(i)="Sangay                         "; 
      i=i+1; volcLat(i)= -0.020; volcLon(i)=268.454; volcElev(i)=  790;  volcLoc(i)="Galapagos                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1503-011"; volcName(i)="Ecuador                        "; 
      i=i+1; volcLat(i)= -0.370; volcLon(i)=268.450; volcElev(i)= 1476;  volcLoc(i)="Galapagos                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1503-01="; volcName(i)="Fernandina                     "; 
      i=i+1; volcLat(i)=  0.020; volcLon(i)=268.650; volcElev(i)= 1710;  volcLoc(i)="Galapagos                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1503-02="; volcName(i)="Wolf                           "; 
      i=i+1; volcLat(i)= -0.180; volcLon(i)=268.720; volcElev(i)= 1330;  volcLoc(i)="Galapagos                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1503-03="; volcName(i)="Darwin                         "; 
      i=i+1; volcLat(i)= -0.430; volcLon(i)=268.880; volcElev(i)= 1130;  volcLoc(i)="Galapagos                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1503-04="; volcName(i)="Alcedo                         "; 
      i=i+1; volcLat(i)= -0.830; volcLon(i)=268.830; volcElev(i)= 1124;  volcLoc(i)="Galapagos                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1503-05="; volcName(i)="Negra, Sierra                  "; 
      i=i+1; volcLat(i)= -0.920; volcLon(i)=268.592; volcElev(i)= 1640;  volcLoc(i)="Galapagos                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1503-06="; volcName(i)="Azul, Cerro                    "; 
      i=i+1; volcLat(i)=  0.580; volcLon(i)=269.250; volcElev(i)=  780;  volcLoc(i)="Galapagos                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1503-07="; volcName(i)="Pinta                          "; 
      i=i+1; volcLat(i)=  0.320; volcLon(i)=270.042; volcElev(i)=   64;  volcLoc(i)="Galapagos                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1503-081"; volcName(i)="Genovesa                       "; 
      i=i+1; volcLat(i)=  0.330; volcLon(i)=269.530; volcElev(i)=  343;  volcLoc(i)="Galapagos                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1503-08="; volcName(i)="Marchena                       "; 
      i=i+1; volcLat(i)= -0.620; volcLon(i)=269.670; volcElev(i)=  864;  volcLoc(i)="Galapagos                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1503-091"; volcName(i)="Santa Cruz                     "; 
      i=i+1; volcLat(i)= -0.220; volcLon(i)=269.230; volcElev(i)=  920;  volcLoc(i)="Galapagos                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1503-09="; volcName(i)="Santiago                       "; 
      i=i+1; volcLat(i)= -0.880; volcLon(i)=270.500; volcElev(i)=  759;  volcLoc(i)="Galapagos                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1503-12-"; volcName(i)="San Cristóbal                  "; 
      i=i+1; volcLat(i)=-14.200; volcLon(i)=288.670; volcElev(i)= 3923;  volcLoc(i)="Peru                          "; 
             volcESP_Code(i)="S0"; volcID(i)="1504-00-"; volcName(i)="Quimsachata                    "; 
      i=i+1; volcLat(i)=-15.070; volcLon(i)=286.820; volcElev(i)= 4980;  volcLoc(i)="Peru                          "; 
             volcESP_Code(i)="M0"; volcID(i)="1504-001"; volcName(i)="Auquihuato, Cerro              "; 
      i=i+1; volcLat(i)=-15.330; volcLon(i)=286.550; volcElev(i)= 5522;  volcLoc(i)="Peru                          "; 
             volcESP_Code(i)="S0"; volcID(i)="1504-002"; volcName(i)="Sara Sara                      "; 
      i=i+1; volcLat(i)=-15.520; volcLon(i)=287.350; volcElev(i)= 6377;  volcLoc(i)="Peru                          "; 
             volcESP_Code(i)="S0"; volcID(i)="1504-003"; volcName(i)="Coropuna                       "; 
      i=i+1; volcLat(i)=-15.420; volcLon(i)=287.670; volcElev(i)= 4713;  volcLoc(i)="Peru                          "; 
             volcESP_Code(i)="M0"; volcID(i)="1504-004"; volcName(i)="Andahua-Orcopampa              "; 
      i=i+1; volcLat(i)=-15.830; volcLon(i)=287.870; volcElev(i)= 4550;  volcLoc(i)="Peru                          "; 
             volcESP_Code(i)="S0"; volcID(i)="1504-005"; volcName(i)="Huambo                         "; 
      i=i+1; volcLat(i)=-15.780; volcLon(i)=288.150; volcElev(i)= 5967;  volcLoc(i)="Peru                          "; 
             volcESP_Code(i)="S0"; volcID(i)="1504-006"; volcName(i)="Sabancaya                      "; 
      i=i+1; volcLat(i)=-16.191; volcLon(i)=288.470; volcElev(i)= 6057;  volcLoc(i)="Peru                          "; 
             volcESP_Code(i)="S0"; volcID(i)="1504-007"; volcName(i)="Chachani, Nevado               "; 
      i=i+1; volcLat(i)=-16.261; volcLon(i)=288.270; volcElev(i)= 2520;  volcLoc(i)="Peru                          "; 
             volcESP_Code(i)="M0"; volcID(i)="1504-008"; volcName(i)="Nicholson, Cerro               "; 
      i=i+1; volcLat(i)=-16.294; volcLon(i)=288.591; volcElev(i)= 5822;  volcLoc(i)="Peru                          "; 
             volcESP_Code(i)="S0"; volcID(i)="1504-01="; volcName(i)="Misti, El                      "; 
      i=i+1; volcLat(i)=-16.355; volcLon(i)=289.097; volcElev(i)= 5672;  volcLoc(i)="Peru                          "; 
             volcESP_Code(i)="S1"; volcID(i)="1504-02="; volcName(i)="Ubinas                         "; 
      i=i+1; volcLat(i)=-16.755; volcLon(i)=289.405; volcElev(i)= 5408;  volcLoc(i)="Peru                          "; 
             volcESP_Code(i)="S0"; volcID(i)="1504-031"; volcName(i)="Ticsani                        "; 
      i=i+1; volcLat(i)=-16.608; volcLon(i)=289.150; volcElev(i)= 4850;  volcLoc(i)="Peru                          "; 
             volcESP_Code(i)="S3"; volcID(i)="1504-03="; volcName(i)="Huaynaputina                   "; 
      i=i+1; volcLat(i)=-17.025; volcLon(i)=289.642; volcElev(i)= 5815;  volcLoc(i)="Peru                          "; 
             volcESP_Code(i)="S0"; volcID(i)="1504-04="; volcName(i)="Tutupaca                       "; 
      i=i+1; volcLat(i)=-17.180; volcLon(i)=289.800; volcElev(i)= 5550;  volcLoc(i)="Peru                          "; 
             volcESP_Code(i)="S0"; volcID(i)="1504-05-"; volcName(i)="Yucamane                       "; 
      i=i+1; volcLat(i)=-17.470; volcLon(i)=290.187; volcElev(i)= 5650;  volcLoc(i)="Peru                          "; 
             volcESP_Code(i)="S0"; volcID(i)="1504-06-"; volcName(i)="Casiri, Nevados                "; 
      i=i+1; volcLat(i)=-18.100; volcLon(i)=290.500; volcElev(i)= 5860;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-011"; volcName(i)="Taapaca                        "; 
      i=i+1; volcLat(i)=-18.170; volcLon(i)=290.850; volcElev(i)= 6348;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-012"; volcName(i)="Parinacota                     "; 
      i=i+1; volcLat(i)=-17.720; volcLon(i)=290.230; volcElev(i)= 5980;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-01="; volcName(i)="Tacora                         "; 
      i=i+1; volcLat(i)=-18.620; volcLon(i)=291.250; volcElev(i)= 4215;  volcLoc(i)="Bolivia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-021"; volcName(i)="Tambo Quemado                  "; 
      i=i+1; volcLat(i)=-18.420; volcLon(i)=290.908; volcElev(i)= 6071;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S1"; volcID(i)="1505-02="; volcName(i)="Guallatiri                     "; 
      i=i+1; volcLat(i)=-19.130; volcLon(i)=291.470; volcElev(i)= 5430;  volcLoc(i)="Bolivia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-032"; volcName(i)="Tata Sabaya                    "; 
      i=i+1; volcLat(i)=-19.450; volcLon(i)=292.580; volcElev(i)= 3650;  volcLoc(i)="Bolivia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-035"; volcName(i)="Jayu Khota, Laguna             "; 
      i=i+1; volcLat(i)=-19.780; volcLon(i)=293.520; volcElev(i)= 5438;  volcLoc(i)="Bolivia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-036"; volcName(i)="Nuevo Mundo                    "; 
      i=i+1; volcLat(i)=-19.150; volcLon(i)=291.170; volcElev(i)= 5550;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S1"; volcID(i)="1505-03="; volcName(i)="Isluga                         "; 
      i=i+1; volcLat(i)=-20.850; volcLon(i)=291.800; volcElev(i)= 5543;  volcLoc(i)="Bolivia                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-042"; volcName(i)="Pampa Luxsar                   "; 
      i=i+1; volcLat(i)=-20.730; volcLon(i)=291.450; volcElev(i)= 5163;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-04="; volcName(i)="Irruputuncu                    "; 
      i=i+1; volcLat(i)=-20.930; volcLon(i)=291.520; volcElev(i)= 5407;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-05="; volcName(i)="Olca-Paruma                    "; 
      i=i+1; volcLat(i)=-21.787; volcLon(i)=291.763; volcElev(i)= 5846;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-061"; volcName(i)="Azufre, Cerro del              "; 
      i=i+1; volcLat(i)=-21.300; volcLon(i)=291.820; volcElev(i)= 5868;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-06="; volcName(i)="Ollagüe                        "; 
      i=i+1; volcLat(i)=-21.880; volcLon(i)=291.600; volcElev(i)= 6145;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S1"; volcID(i)="1505-07="; volcName(i)="San Pedro                      "; 
      i=i+1; volcLat(i)=-22.720; volcLon(i)=292.108; volcElev(i)= 5971;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-091"; volcName(i)="Sairecabur                     "; 
      i=i+1; volcLat(i)=-22.830; volcLon(i)=292.120; volcElev(i)= 5916;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-092"; volcName(i)="Licancabur                     "; 
      i=i+1; volcLat(i)=-22.895; volcLon(i)=292.434; volcElev(i)= 5598;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-093"; volcName(i)="Guayaques                      "; 
      i=i+1; volcLat(i)=-23.000; volcLon(i)=292.250; volcElev(i)= 5703;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-094"; volcName(i)="Purico Complex                 "; 
      i=i+1; volcLat(i)=-23.236; volcLon(i)=292.355; volcElev(i)= 5631;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-095"; volcName(i)="Colachi                        "; 
      i=i+1; volcLat(i)=-23.300; volcLon(i)=292.380; volcElev(i)= 6046;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-096"; volcName(i)="Acamarachi                     "; 
      i=i+1; volcLat(i)=-23.520; volcLon(i)=292.330; volcElev(i)= 4555;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-097"; volcName(i)="Overo, Cerro                   "; 
      i=i+1; volcLat(i)=-23.580; volcLon(i)=292.300; volcElev(i)= 5778;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-098"; volcName(i)="Chiliques                      "; 
      i=i+1; volcLat(i)=-22.550; volcLon(i)=292.150; volcElev(i)= 5890;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-09="; volcName(i)="Putana                         "; 
      i=i+1; volcLat(i)=-23.743; volcLon(i)=292.466; volcElev(i)= 5852;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-101"; volcName(i)="Cordón de Puntas Negras        "; 
      i=i+1; volcLat(i)=-23.820; volcLon(i)=292.230; volcElev(i)= 5910;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-102"; volcName(i)="Miñiques                       "; 
      i=i+1; volcLat(i)=-23.830; volcLon(i)=292.050; volcElev(i)= 3550;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="M0"; volcID(i)="1505-103"; volcName(i)="Tujle, Cerro                   "; 
      i=i+1; volcLat(i)=-23.950; volcLon(i)=292.270; volcElev(i)= 4450;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-104"; volcName(i)="Caichinque                     "; 
      i=i+1; volcLat(i)=-23.970; volcLon(i)=291.870; volcElev(i)= 3116;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-105"; volcName(i)="Tilocalar                      "; 
      i=i+1; volcLat(i)=-24.180; volcLon(i)=291.750; volcElev(i)= 3500;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="M0"; volcID(i)="1505-106"; volcName(i)="Negrillar, El                  "; 
      i=i+1; volcLat(i)=-24.188; volcLon(i)=291.946; volcElev(i)= 6233;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-107"; volcName(i)="Pular                          "; 
      i=i+1; volcLat(i)=-24.280; volcLon(i)=291.400; volcElev(i)= 4109;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="M0"; volcID(i)="1505-108"; volcName(i)="Negrillar, La                  "; 
      i=i+1; volcLat(i)=-24.400; volcLon(i)=291.750; volcElev(i)= 6051;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-109"; volcName(i)="Socompa                        "; 
      i=i+1; volcLat(i)=-23.370; volcLon(i)=292.270; volcElev(i)= 5592;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S2"; volcID(i)="1505-10="; volcName(i)="Láscar                         "; 
      i=i+1; volcLat(i)=-25.080; volcLon(i)=291.630; volcElev(i)= 5447;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-112"; volcName(i)="Escorial, Cerro                "; 
      i=i+1; volcLat(i)=-24.720; volcLon(i)=291.470; volcElev(i)= 6739;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-11="; volcName(i)="Llullaillaco                   "; 
      i=i+1; volcLat(i)=-25.330; volcLon(i)=291.480; volcElev(i)= 5463;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-121"; volcName(i)="Cordón del Azufre              "; 
      i=i+1; volcLat(i)=-25.420; volcLon(i)=291.420; volcElev(i)= 5401;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-122"; volcName(i)="Bayo, Cerro                    "; 
      i=i+1; volcLat(i)=-26.480; volcLon(i)=291.420; volcElev(i)= 6127;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-123"; volcName(i)="Nevada, Sierra                 "; 
      i=i+1; volcLat(i)=-26.800; volcLon(i)=291.630; volcElev(i)= 5890;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-124"; volcName(i)="Falso Azufre                   "; 
      i=i+1; volcLat(i)=-27.042; volcLon(i)=291.720; volcElev(i)= 6621;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-125"; volcName(i)="Incahuasi, Nevado de           "; 
      i=i+1; volcLat(i)=-25.170; volcLon(i)=291.500; volcElev(i)= 5697;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-12="; volcName(i)="Lastarria                      "; 
      i=i+1; volcLat(i)=-27.108; volcLon(i)=291.280; volcElev(i)= 6190;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-131"; volcName(i)="Solo, El                       "; 
      i=i+1; volcLat(i)=-27.120; volcLon(i)=291.450; volcElev(i)= 6887;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-13="; volcName(i)="Ojos del Salado, Nevados       "; 
      i=i+1; volcLat(i)=-27.300; volcLon(i)=290.870; volcElev(i)= 6052;  volcLoc(i)="Chile-N                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-14-"; volcName(i)="Copiapó                        "; 
      i=i+1; volcLat(i)=-24.050; volcLon(i)=293.520; volcElev(i)= 5500;  volcLoc(i)="Argentina                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-15-"; volcName(i)="Tuzgle, Cerro                  "; 
      i=i+1; volcLat(i)=-24.250; volcLon(i)=292.230; volcElev(i)= 6082;  volcLoc(i)="Argentina                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-16-"; volcName(i)="Aracar                         "; 
      i=i+1; volcLat(i)=-25.100; volcLon(i)=291.730; volcElev(i)=    0;  volcLoc(i)="Argentina                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1505-161"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=-26.080; volcLon(i)=292.500; volcElev(i)= 4000;  volcLoc(i)="Argentina                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1505-18-"; volcName(i)="Antofagasta de la Sierra       "; 
      i=i+1; volcLat(i)=-26.620; volcLon(i)=291.650; volcElev(i)= 6532;  volcLoc(i)="Argentina                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-19-"; volcName(i)="Cóndor, Cerro el               "; 
      i=i+1; volcLat(i)=-26.620; volcLon(i)=291.850; volcElev(i)= 5740;  volcLoc(i)="Argentina                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-20-"; volcName(i)="Peinado                        "; 
      i=i+1; volcLat(i)=-26.770; volcLon(i)=292.280; volcElev(i)= 4400;  volcLoc(i)="Argentina                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-21-"; volcName(i)="Robledo                        "; 
      i=i+1; volcLat(i)=-27.200; volcLon(i)=291.450; volcElev(i)= 6660;  volcLoc(i)="Argentina                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1505-22-"; volcName(i)="Tipas                          "; 
      i=i+1; volcLat(i)=-27.150; volcLon(i)=250.620; volcElev(i)=  511;  volcLoc(i)="Chile-Is                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1506-011"; volcName(i)="Easter Island                  "; 
      i=i+1; volcLat(i)=-26.280; volcLon(i)=279.880; volcElev(i)=  193;  volcLoc(i)="Chile-Is                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1506-01="; volcName(i)="San Félix                      "; 
      i=i+1; volcLat(i)=-33.658; volcLon(i)=281.150; volcElev(i)=  922;  volcLoc(i)="Chile-Is                      "; 
             volcESP_Code(i)="M0"; volcID(i)="1506-02="; volcName(i)="Robinson Crusoe                "; 
      i=i+1; volcLat(i)=-33.620; volcLon(i)=283.170; volcElev(i)= -642;  volcLoc(i)="Chile-Is                      "; 
             volcESP_Code(i)="U0"; volcID(i)="1506-04="; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=-33.400; volcLon(i)=290.200; volcElev(i)= 6000;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S1"; volcID(i)="1507-01="; volcName(i)="Tupungatito                    "; 
      i=i+1; volcLat(i)=-34.161; volcLon(i)=290.167; volcElev(i)= 5264;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S1"; volcID(i)="1507-021"; volcName(i)="Maipo                          "; 
      i=i+1; volcLat(i)=-34.608; volcLon(i)=289.705; volcElev(i)= 4860;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-022"; volcName(i)="Palomo                         "; 
      i=i+1; volcLat(i)=-34.650; volcLon(i)=289.950; volcElev(i)= 5189;  volcLoc(i)="Argentina                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-023"; volcName(i)="Atuel, Caldera del             "; 
      i=i+1; volcLat(i)=-34.930; volcLon(i)=290.000; volcElev(i)= 4999;  volcLoc(i)="Argentina                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-024"; volcName(i)="Risco Plateado                 "; 
      i=i+1; volcLat(i)=-33.782; volcLon(i)=290.103; volcElev(i)= 5856;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S1"; volcID(i)="1507-02="; volcName(i)="San José                       "; 
      i=i+1; volcLat(i)=-34.814; volcLon(i)=289.648; volcElev(i)= 4280;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-03="; volcName(i)="Tinguiririca                   "; 
      i=i+1; volcLat(i)=-35.558; volcLon(i)=289.504; volcElev(i)= 3508;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-042"; volcName(i)="Calabozos                      "; 
      i=i+1; volcLat(i)=-35.240; volcLon(i)=289.430; volcElev(i)= 4107;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S1"; volcID(i)="1507-04="; volcName(i)="Planchón-Peteroa               "; 
      i=i+1; volcLat(i)=-35.580; volcLon(i)=289.250; volcElev(i)= 3953;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-05="; volcName(i)="Descabezado Grande             "; 
      i=i+1; volcLat(i)=-36.020; volcLon(i)=289.420; volcElev(i)= 3092;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-061"; volcName(i)="Maule, Laguna del              "; 
      i=i+1; volcLat(i)=-35.989; volcLon(i)=289.151; volcElev(i)= 3621;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-062"; volcName(i)="San Pedro-Pellado              "; 
      i=i+1; volcLat(i)=-36.193; volcLon(i)=288.839; volcElev(i)= 3242;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-063"; volcName(i)="Longaví, Nevado de             "; 
      i=i+1; volcLat(i)=-36.286; volcLon(i)=288.991; volcElev(i)= 2268;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-064"; volcName(i)="Blancas, Lomas                 "; 
      i=i+1; volcLat(i)=-36.450; volcLon(i)=289.080; volcElev(i)= 1890;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="M0"; volcID(i)="1507-065"; volcName(i)="Resago                         "; 
      i=i+1; volcLat(i)=-36.420; volcLon(i)=290.800; volcElev(i)= 3680;  volcLoc(i)="Argentina                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1507-066"; volcName(i)="Payún Matru                    "; 
      i=i+1; volcLat(i)=-36.580; volcLon(i)=289.580; volcElev(i)= 4709;  volcLoc(i)="Argentina                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-067"; volcName(i)="Domuyo                         "; 
      i=i+1; volcLat(i)=-35.653; volcLon(i)=289.239; volcElev(i)= 3788;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S2"; volcID(i)="1507-06="; volcName(i)="Azul, Cerro                    "; 
      i=i+1; volcLat(i)=-36.770; volcLon(i)=290.180; volcElev(i)= 1435;  volcLoc(i)="Argentina                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-071"; volcName(i)="Cochiquito Volc Group          "; 
      i=i+1; volcLat(i)=-37.142; volcLon(i)=289.970; volcElev(i)= 3978;  volcLoc(i)="Argentina                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-072"; volcName(i)="Tromen                         "; 
      i=i+1; volcLat(i)=-37.570; volcLon(i)=290.380; volcElev(i)=  970;  volcLoc(i)="Argentina                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1507-073"; volcName(i)="Puesto Cortaderas              "; 
      i=i+1; volcLat(i)=-36.863; volcLon(i)=288.623; volcElev(i)= 3212;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S1"; volcID(i)="1507-07="; volcName(i)="Chillán, Nevados de            "; 
      i=i+1; volcLat(i)=-37.730; volcLon(i)=289.100; volcElev(i)= 2500;  volcLoc(i)="Argentina                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-081"; volcName(i)="Trocon                         "; 
      i=i+1; volcLat(i)=-37.406; volcLon(i)=288.651; volcElev(i)= 2979;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S1"; volcID(i)="1507-08="; volcName(i)="Antuco                         "; 
      i=i+1; volcLat(i)=-37.920; volcLon(i)=288.550; volcElev(i)= 3164;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S1"; volcID(i)="1507-091"; volcName(i)="Callaqui                       "; 
      i=i+1; volcLat(i)=-38.270; volcLon(i)=288.900; volcElev(i)= 2143;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="M0"; volcID(i)="1507-092"; volcName(i)="Mariñaqui, Laguna              "; 
      i=i+1; volcLat(i)=-38.310; volcLon(i)=288.355; volcElev(i)= 2806;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-093"; volcName(i)="Tolguaca                       "; 
      i=i+1; volcLat(i)=-37.850; volcLon(i)=288.830; volcElev(i)= 2997;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S1"; volcID(i)="1507-09="; volcName(i)="Copahue                        "; 
      i=i+1; volcLat(i)=-38.377; volcLon(i)=288.420; volcElev(i)= 2865;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S2"; volcID(i)="1507-10="; volcName(i)="Lonquimay                      "; 
      i=i+1; volcLat(i)=-38.970; volcLon(i)=288.480; volcElev(i)= 2282;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-111"; volcName(i)="Sollipulli                     "; 
      i=i+1; volcLat(i)=-39.250; volcLon(i)=288.300; volcElev(i)= 1496;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="M0"; volcID(i)="1507-112"; volcName(i)="Caburgua-Huelemolle            "; 
      i=i+1; volcLat(i)=-38.692; volcLon(i)=288.271; volcElev(i)= 3125;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S1"; volcID(i)="1507-11="; volcName(i)="Llaima                         "; 
      i=i+1; volcLat(i)=-39.500; volcLon(i)=288.300; volcElev(i)= 2360;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-121"; volcName(i)="Quetrupillan                   "; 
      i=i+1; volcLat(i)=-39.633; volcLon(i)=288.500; volcElev(i)= 3747;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-122"; volcName(i)="Lanín                          "; 
      i=i+1; volcLat(i)=-39.880; volcLon(i)=288.420; volcElev(i)= 2139;  volcLoc(i)="Argentina                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-123"; volcName(i)="Huanquihue Group               "; 
      i=i+1; volcLat(i)=-39.420; volcLon(i)=288.070; volcElev(i)= 2847;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S1"; volcID(i)="1507-12="; volcName(i)="Villarrica                     "; 
      i=i+1; volcLat(i)=-39.927; volcLon(i)=287.973; volcElev(i)= 2422;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-13="; volcName(i)="Mocho-Choshuenco               "; 
      i=i+1; volcLat(i)=-40.350; volcLon(i)=287.930; volcElev(i)= 1114;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="M2"; volcID(i)="1507-14="; volcName(i)="Carrán-Los Venados             "; 
      i=i+1; volcLat(i)=-40.770; volcLon(i)=288.050; volcElev(i)= 2024;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-152"; volcName(i)="Pantoja, Cerro                 "; 
      i=i+1; volcLat(i)=-40.771; volcLon(i)=287.847; volcElev(i)= 1990;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-153"; volcName(i)="Antillanca Group               "; 
      i=i+1; volcLat(i)=-40.590; volcLon(i)=287.883; volcElev(i)= 2236;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-15="; volcName(i)="Puyehue-Cordón Caulle          "; 
      i=i+1; volcLat(i)=-40.969; volcLon(i)=287.736; volcElev(i)= 2493;  volcLoc(i)="Chile-C                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1507-16-"; volcName(i)="Puntiagudo-Cordón Cenizos      "; 
      i=i+1; volcLat(i)=-41.157; volcLon(i)=288.115; volcElev(i)= 3491;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1508-011"; volcName(i)="Tronador                       "; 
      i=i+1; volcLat(i)=-41.250; volcLon(i)=287.730; volcElev(i)=  506;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="M0"; volcID(i)="1508-012"; volcName(i)="Cayutué-La Viguería            "; 
      i=i+1; volcLat(i)=-41.100; volcLon(i)=287.507; volcElev(i)= 2652;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S1"; volcID(i)="1508-01="; volcName(i)="Osorno                         "; 
      i=i+1; volcLat(i)=-41.400; volcLon(i)=288.000; volcElev(i)= 1862;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1508-021"; volcName(i)="Cuernos del Diablo             "; 
      i=i+1; volcLat(i)=-41.755; volcLon(i)=287.604; volcElev(i)= 2187;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1508-022"; volcName(i)="Yate                           "; 
      i=i+1; volcLat(i)=-41.874; volcLon(i)=287.569; volcElev(i)= 1572;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1508-023"; volcName(i)="Hornopirén                     "; 
      i=i+1; volcLat(i)=-41.880; volcLon(i)=287.420; volcElev(i)= 1210;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1508-024"; volcName(i)="Apagado                        "; 
      i=i+1; volcLat(i)=-42.020; volcLon(i)=289.820; volcElev(i)= 1359;  volcLoc(i)="Chile-S/Argentina             "; 
             volcESP_Code(i)="M0"; volcID(i)="1508-025"; volcName(i)="Crater Basalt Volc Field       "; 
      i=i+1; volcLat(i)=-41.326; volcLon(i)=287.386; volcElev(i)= 2003;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S2"; volcID(i)="1508-02="; volcName(i)="Calbuco                        "; 
      i=i+1; volcLat(i)=-42.377; volcLon(i)=287.422; volcElev(i)= 1318;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S1"; volcID(i)="1508-03="; volcName(i)="Huequi                         "; 
      i=i+1; volcLat(i)=-42.833; volcLon(i)=287.354; volcElev(i)= 1122;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S2"; volcID(i)="1508-041"; volcName(i)="Chaitén                        "; 
      i=i+1; volcLat(i)=-42.793; volcLon(i)=287.561; volcElev(i)= 2404;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S1"; volcID(i)="1508-04="; volcName(i)="Minchinmávida                  "; 
      i=i+1; volcLat(i)=-43.500; volcLon(i)=287.200; volcElev(i)= 2042;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1508-050"; volcName(i)="Yanteles                       "; 
      i=i+1; volcLat(i)=-43.780; volcLon(i)=287.530; volcElev(i)=    0;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="M0"; volcID(i)="1508-051"; volcName(i)="Palena Volc Group              "; 
      i=i+1; volcLat(i)=-44.080; volcLon(i)=287.120; volcElev(i)= 2400;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1508-052"; volcName(i)="Melimoyu                       "; 
      i=i+1; volcLat(i)=-44.300; volcLon(i)=287.470; volcElev(i)=  524;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="M0"; volcID(i)="1508-053"; volcName(i)="Puyuhuapi                      "; 
      i=i+1; volcLat(i)=-44.700; volcLon(i)=286.920; volcElev(i)= 1660;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1508-054"; volcName(i)="Mentolat                       "; 
      i=i+1; volcLat(i)=-45.059; volcLon(i)=287.016; volcElev(i)= 2090;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1508-055"; volcName(i)="Cay                            "; 
      i=i+1; volcLat(i)=-45.100; volcLon(i)=286.830; volcElev(i)= 2960;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1508-056"; volcName(i)="Maca                           "; 
      i=i+1; volcLat(i)=-45.900; volcLon(i)=287.030; volcElev(i)= 1905;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S3"; volcID(i)="1508-057"; volcName(i)="Hudson, Cerro                  "; 
      i=i+1; volcLat(i)=-46.170; volcLon(i)=287.330; volcElev(i)=    0;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="M0"; volcID(i)="1508-058"; volcName(i)="Río Murta                      "; 
      i=i+1; volcLat(i)=-47.200; volcLon(i)=286.520; volcElev(i)= 3437;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1508-059"; volcName(i)="Arenales                       "; 
      i=i+1; volcLat(i)=-43.180; volcLon(i)=287.200; volcElev(i)= 2300;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1508-05="; volcName(i)="Corcovado                      "; 
      i=i+1; volcLat(i)=-49.358; volcLon(i)=286.720; volcElev(i)= 1500;  volcLoc(i)="Argentina                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1508-061"; volcName(i)="Viedma                         "; 
      i=i+1; volcLat(i)=-50.330; volcLon(i)=286.250; volcElev(i)= 2546;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1508-062"; volcName(i)="Aguilera                       "; 
      i=i+1; volcLat(i)=-50.964; volcLon(i)=286.415; volcElev(i)= 1000;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1508-063"; volcName(i)="Reclus                         "; 
      i=i+1; volcLat(i)=-49.020; volcLon(i)=286.450; volcElev(i)= 3607;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S1"; volcID(i)="1508-06="; volcName(i)="Lautaro                        "; 
      i=i+1; volcLat(i)=-52.330; volcLon(i)=286.600; volcElev(i)= 1758;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1508-07="; volcName(i)="Burney, Monte                  "; 
      i=i+1; volcLat(i)=-52.000; volcLon(i)=290.000; volcElev(i)=  282;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="M0"; volcID(i)="1508-08-"; volcName(i)="Palei-Aike Volc Field          "; 
      i=i+1; volcLat(i)=-54.950; volcLon(i)=289.750; volcElev(i)=  150;  volcLoc(i)="Chile-S                       "; 
             volcESP_Code(i)="S0"; volcID(i)="1508-09-"; volcName(i)="Fueguino                       "; 
      i=i+1; volcLat(i)= 17.630; volcLon(i)=296.770; volcElev(i)=  887;  volcLoc(i)="W Indies                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1600-01="; volcName(i)="Saba                           "; 
      i=i+1; volcLat(i)= 17.478; volcLon(i)=297.040; volcElev(i)=  601;  volcLoc(i)="W Indies                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1600-02="; volcName(i)="Quill, The                     "; 
      i=i+1; volcLat(i)= 17.370; volcLon(i)=297.200; volcElev(i)= 1156;  volcLoc(i)="W Indies                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1600-03="; volcName(i)="Liamuiga                       "; 
      i=i+1; volcLat(i)= 17.150; volcLon(i)=297.420; volcElev(i)=  985;  volcLoc(i)="W Indies                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1600-04="; volcName(i)="Nevis Peak                     "; 
      i=i+1; volcLat(i)= 16.720; volcLon(i)=297.820; volcElev(i)=  915;  volcLoc(i)="W Indies                      "; 
             volcESP_Code(i)="S9"; volcID(i)="1600-05="; volcName(i)="Soufrière Hills                "; 
      i=i+1; volcLat(i)= 16.050; volcLon(i)=298.330; volcElev(i)= 1467;  volcLoc(i)="W Indies                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1600-06="; volcName(i)="Soufrière Guadeloupe           "; 
      i=i+1; volcLat(i)= 15.612; volcLon(i)=298.570; volcElev(i)=  861;  volcLoc(i)="W Indies                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1600-08="; volcName(i)="Diables, Morne aux             "; 
      i=i+1; volcLat(i)= 15.503; volcLon(i)=298.603; volcElev(i)= 1430;  volcLoc(i)="W Indies                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1600-09="; volcName(i)="Diablotins, Morne              "; 
      i=i+1; volcLat(i)= 15.307; volcLon(i)=298.695; volcElev(i)= 1224;  volcLoc(i)="W Indies                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1600-101"; volcName(i)="Watt, Morne                    "; 
      i=i+1; volcLat(i)= 15.370; volcLon(i)=298.670; volcElev(i)= 1387;  volcLoc(i)="W Indies                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1600-10="; volcName(i)="Trois Pitons, Morne            "; 
      i=i+1; volcLat(i)= 15.255; volcLon(i)=298.659; volcElev(i)=  940;  volcLoc(i)="W Indies                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1600-11="; volcName(i)="Plat Pays, Morne               "; 
      i=i+1; volcLat(i)= 14.820; volcLon(i)=298.830; volcElev(i)= 1397;  volcLoc(i)="W Indies                      "; 
             volcESP_Code(i)="S2"; volcID(i)="1600-12="; volcName(i)="Pelée                          "; 
      i=i+1; volcLat(i)= 13.830; volcLon(i)=298.950; volcElev(i)=  777;  volcLoc(i)="W Indies                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1600-14="; volcName(i)="Qualibou                       "; 
      i=i+1; volcLat(i)= 13.330; volcLon(i)=298.820; volcElev(i)= 1220;  volcLoc(i)="W Indies                      "; 
             volcESP_Code(i)="S2"; volcID(i)="1600-15="; volcName(i)="Soufrière St. Vincent          "; 
      i=i+1; volcLat(i)= 12.300; volcLon(i)=298.360; volcElev(i)= -185;  volcLoc(i)="W Indies                      "; 
             volcESP_Code(i)="U0"; volcID(i)="1600-16="; volcName(i)="Kick 'em Jenny                 "; 
      i=i+1; volcLat(i)= 12.150; volcLon(i)=298.330; volcElev(i)=  840;  volcLoc(i)="W Indies                      "; 
             volcESP_Code(i)="S0"; volcID(i)="1600-17="; volcName(i)="St. Catherine                  "; 
      i=i+1; volcLat(i)= 64.800; volcLon(i)=336.220; volcElev(i)= 1448;  volcLoc(i)="Iceland-W                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1700-01="; volcName(i)="Snaefellsjökull                "; 
      i=i+1; volcLat(i)= 64.870; volcLon(i)=336.750; volcElev(i)=  647;  volcLoc(i)="Iceland-W                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1700-02="; volcName(i)="Helgrindur                     "; 
      i=i+1; volcLat(i)= 64.870; volcLon(i)=337.770; volcElev(i)= 1063;  volcLoc(i)="Iceland-W                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1700-03="; volcName(i)="Ljósufjöll                     "; 
      i=i+1; volcLat(i)= 63.880; volcLon(i)=337.500; volcElev(i)=  230;  volcLoc(i)="Iceland-SW                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1701-02="; volcName(i)="Reykjanes                      "; 
      i=i+1; volcLat(i)= 63.930; volcLon(i)=337.900; volcElev(i)=  379;  volcLoc(i)="Iceland-SW                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1701-03="; volcName(i)="Krísuvík                       "; 
      i=i+1; volcLat(i)= 63.920; volcLon(i)=338.170; volcElev(i)=  621;  volcLoc(i)="Iceland-SW                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1701-04="; volcName(i)="Brennisteinsfjöll              "; 
      i=i+1; volcLat(i)= 64.073; volcLon(i)=338.798; volcElev(i)=  540;  volcLoc(i)="Iceland-S                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1701-051"; volcName(i)="Hrómundartindur                "; 
      i=i+1; volcLat(i)= 64.080; volcLon(i)=338.680; volcElev(i)=  803;  volcLoc(i)="Iceland-SW                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1701-05="; volcName(i)="Hengill                        "; 
      i=i+1; volcLat(i)= 64.030; volcLon(i)=339.130; volcElev(i)=  214;  volcLoc(i)="Iceland-SW                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1701-06="; volcName(i)="Grímsnes                       "; 
      i=i+1; volcLat(i)= 64.600; volcLon(i)=339.420; volcElev(i)= 1400;  volcLoc(i)="Iceland-SW                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1701-07="; volcName(i)="Prestahnukur                   "; 
      i=i+1; volcLat(i)= 64.750; volcLon(i)=340.020; volcElev(i)= 1360;  volcLoc(i)="Iceland-SW                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1701-08="; volcName(i)="Hveravellir                    "; 
      i=i+1; volcLat(i)= 64.780; volcLon(i)=341.080; volcElev(i)= 1782;  volcLoc(i)="Iceland-SW                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1701-09="; volcName(i)="Hofsjökull                     "; 
      i=i+1; volcLat(i)= 63.430; volcLon(i)=339.720; volcElev(i)=  279;  volcLoc(i)="Iceland-S                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1702-01="; volcName(i)="Vestmannaeyjar                 "; 
      i=i+1; volcLat(i)= 63.630; volcLon(i)=340.380; volcElev(i)= 1666;  volcLoc(i)="Iceland-S                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1702-02="; volcName(i)="Eyjafjöll                      "; 
      i=i+1; volcLat(i)= 63.630; volcLon(i)=340.950; volcElev(i)= 1512;  volcLoc(i)="Iceland-S                     "; 
             volcESP_Code(i)="M3"; volcID(i)="1702-03="; volcName(i)="Katla                          "; 
      i=i+1; volcLat(i)= 63.780; volcLon(i)=340.430; volcElev(i)= 1463;  volcLoc(i)="Iceland-S                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1702-04="; volcName(i)="Tindfjallajökull               "; 
      i=i+1; volcLat(i)= 63.920; volcLon(i)=340.830; volcElev(i)= 1259;  volcLoc(i)="Iceland-S                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1702-05="; volcName(i)="Torfajökull                    "; 
      i=i+1; volcLat(i)= 63.980; volcLon(i)=340.300; volcElev(i)= 1491;  volcLoc(i)="Iceland-S                     "; 
             volcESP_Code(i)="S2"; volcID(i)="1702-07="; volcName(i)="Hekla                          "; 
      i=i+1; volcLat(i)= 64.420; volcLon(i)=342.670; volcElev(i)= 1725;  volcLoc(i)="Iceland-NE                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1703-01="; volcName(i)="Grímsvötn                      "; 
      i=i+1; volcLat(i)= 64.630; volcLon(i)=342.470; volcElev(i)= 2009;  volcLoc(i)="Iceland-NE                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1703-03="; volcName(i)="Bárdarbunga                    "; 
      i=i+1; volcLat(i)= 64.730; volcLon(i)=342.080; volcElev(i)= 1535;  volcLoc(i)="Iceland-NE                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1703-04="; volcName(i)="Tungnafellsjökull              "; 
      i=i+1; volcLat(i)= 64.650; volcLon(i)=343.280; volcElev(i)= 1929;  volcLoc(i)="Iceland-NE                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1703-05="; volcName(i)="Kverkfjöll                     "; 
      i=i+1; volcLat(i)= 65.030; volcLon(i)=343.250; volcElev(i)= 1516;  volcLoc(i)="Iceland-NE                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1703-06="; volcName(i)="Askja                          "; 
      i=i+1; volcLat(i)= 65.430; volcLon(i)=343.350; volcElev(i)=  939;  volcLoc(i)="Iceland-NE                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1703-07="; volcName(i)="Fremrinamur                    "; 
      i=i+1; volcLat(i)= 65.730; volcLon(i)=343.220; volcElev(i)=  818;  volcLoc(i)="Iceland-NE                    "; 
             volcESP_Code(i)="M1"; volcID(i)="1703-08="; volcName(i)="Krafla                         "; 
      i=i+1; volcLat(i)= 65.880; volcLon(i)=343.170; volcElev(i)=  564;  volcLoc(i)="Iceland-NE                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1703-09="; volcName(i)="Theistareykjarbunga            "; 
      i=i+1; volcLat(i)= 66.300; volcLon(i)=342.900; volcElev(i)=    0;  volcLoc(i)="Iceland-N of                  "; 
             volcESP_Code(i)="U0"; volcID(i)="1703-10="; volcName(i)="Tjörnes Fracture Zone          "; 
      i=i+1; volcLat(i)= 64.000; volcLon(i)=343.350; volcElev(i)= 2119;  volcLoc(i)="Iceland-SE                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1704-01="; volcName(i)="Öraefajökull                   "; 
      i=i+1; volcLat(i)= 64.270; volcLon(i)=343.350; volcElev(i)= 1760;  volcLoc(i)="Iceland-SE                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1704-02="; volcName(i)="Esjufjöll                      "; 
      i=i+1; volcLat(i)= 66.670; volcLon(i)=341.500; volcElev(i)=    5;  volcLoc(i)="Iceland-N of                  "; 
             volcESP_Code(i)="U0"; volcID(i)="1705-01="; volcName(i)="Kolbeinsey Ridge               "; 
      i=i+1; volcLat(i)= 71.080; volcLon(i)=351.830; volcElev(i)= 2277;  volcLoc(i)="Atlantic-N-Jan Mayen          "; 
             volcESP_Code(i)="M0"; volcID(i)="1706-01="; volcName(i)="Jan Mayen                      "; 
      i=i+1; volcLat(i)= 88.270; volcLon(i)=294.400; volcElev(i)=-1500;  volcLoc(i)="Arctic Ocean                  "; 
             volcESP_Code(i)="U0"; volcID(i)="1707-01-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 85.580; volcLon(i)= 85.000; volcElev(i)=-3800;  volcLoc(i)="Arctic Ocean                  "; 
             volcESP_Code(i)="U0"; volcID(i)="1707-02-"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 49.000; volcLon(i)=325.500; volcElev(i)=-1650;  volcLoc(i)="Atlantic-N                    "; 
             volcESP_Code(i)="U0"; volcID(i)="1801-02="; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 39.950; volcLon(i)=334.170; volcElev(i)=-2835;  volcLoc(i)="Atlantic-N                    "; 
             volcESP_Code(i)="U0"; volcID(i)="1801-03="; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 38.750; volcLon(i)=321.920; volcElev(i)=-4200;  volcLoc(i)="Atlantic-N                    "; 
             volcESP_Code(i)="U0"; volcID(i)="1801-04="; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= 39.462; volcLon(i)=328.784; volcElev(i)=  914;  volcLoc(i)="Azores                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1802-001"; volcName(i)="Flores                         "; 
      i=i+1; volcLat(i)= 39.699; volcLon(i)=328.889; volcElev(i)=  718;  volcLoc(i)="Azores                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1802-002"; volcName(i)="Corvo                          "; 
      i=i+1; volcLat(i)= 38.600; volcLon(i)=331.270; volcElev(i)= 1043;  volcLoc(i)="Azores                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1802-01="; volcName(i)="Fayal                          "; 
      i=i+1; volcLat(i)= 38.470; volcLon(i)=331.600; volcElev(i)= 2351;  volcLoc(i)="Azores                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1802-02="; volcName(i)="Pico                           "; 
      i=i+1; volcLat(i)= 38.650; volcLon(i)=331.920; volcElev(i)= 1053;  volcLoc(i)="Azores                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1802-03="; volcName(i)="San Jorge                      "; 
      i=i+1; volcLat(i)= 39.020; volcLon(i)=332.030; volcElev(i)=  402;  volcLoc(i)="Azores                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1802-04="; volcName(i)="Graciosa                       "; 
      i=i+1; volcLat(i)= 38.730; volcLon(i)=332.680; volcElev(i)= 1023;  volcLoc(i)="Azores                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1802-05="; volcName(i)="Terceira                       "; 
      i=i+1; volcLat(i)= 38.230; volcLon(i)=333.370; volcElev(i)=  -13;  volcLoc(i)="Azores                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1802-07="; volcName(i)="Don Joao de Castro Bank        "; 
      i=i+1; volcLat(i)= 37.780; volcLon(i)=334.330; volcElev(i)=  350;  volcLoc(i)="Azores                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1802-081"; volcName(i)="Picos Volc System              "; 
      i=i+1; volcLat(i)= 37.870; volcLon(i)=334.220; volcElev(i)=  856;  volcLoc(i)="Azores                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1802-08="; volcName(i)="Sete Cidades                   "; 
      i=i+1; volcLat(i)= 37.770; volcLon(i)=334.530; volcElev(i)=  947;  volcLoc(i)="Azores                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1802-09="; volcName(i)="Agua de Pau                    "; 
      i=i+1; volcLat(i)= 37.770; volcLon(i)=334.680; volcElev(i)=  805;  volcLoc(i)="Azores                        "; 
             volcESP_Code(i)="S0"; volcID(i)="1802-10="; volcName(i)="Furnas                         "; 
      i=i+1; volcLat(i)= 37.600; volcLon(i)=334.120; volcElev(i)= -197;  volcLoc(i)="Azores                        "; 
             volcESP_Code(i)="U0"; volcID(i)="1802-11="; volcName(i)="Monaco Bank                    "; 
      i=i+1; volcLat(i)= 32.730; volcLon(i)=343.030; volcElev(i)= 1862;  volcLoc(i)="Azores                        "; 
             volcESP_Code(i)="M0"; volcID(i)="1802-12-"; volcName(i)="Madeira                        "; 
      i=i+1; volcLat(i)= 28.570; volcLon(i)=342.170; volcElev(i)= 2426;  volcLoc(i)="Canary Is                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1803-01-"; volcName(i)="La Palma                       "; 
      i=i+1; volcLat(i)= 27.730; volcLon(i)=341.970; volcElev(i)= 1500;  volcLoc(i)="Canary Is                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1803-02-"; volcName(i)="Hierro                         "; 
      i=i+1; volcLat(i)= 28.271; volcLon(i)=343.359; volcElev(i)= 3715;  volcLoc(i)="Canary Is                     "; 
             volcESP_Code(i)="S0"; volcID(i)="1803-03-"; volcName(i)="Tenerife                       "; 
      i=i+1; volcLat(i)= 28.000; volcLon(i)=344.420; volcElev(i)= 1950;  volcLoc(i)="Canary Is                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1803-04-"; volcName(i)="Gran Canaria                   "; 
      i=i+1; volcLat(i)= 28.358; volcLon(i)=345.980; volcElev(i)=  529;  volcLoc(i)="Canary Is                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1803-05-"; volcName(i)="Fuerteventura                  "; 
      i=i+1; volcLat(i)= 29.030; volcLon(i)=346.370; volcElev(i)=  670;  volcLoc(i)="Canary Is                     "; 
             volcESP_Code(i)="M0"; volcID(i)="1803-06-"; volcName(i)="Lanzarote                      "; 
      i=i+1; volcLat(i)= 14.950; volcLon(i)=335.650; volcElev(i)= 2829;  volcLoc(i)="Cape Verde Is                 "; 
             volcESP_Code(i)="S0"; volcID(i)="1804-01="; volcName(i)="Fogo                           "; 
      i=i+1; volcLat(i)= 14.850; volcLon(i)=335.280; volcElev(i)=  900;  volcLoc(i)="Cape Verde Is                 "; 
             volcESP_Code(i)="S0"; volcID(i)="1804-02-"; volcName(i)="Brava                          "; 
      i=i+1; volcLat(i)= 16.850; volcLon(i)=335.030; volcElev(i)=  725;  volcLoc(i)="Cape Verde Is                 "; 
             volcESP_Code(i)="S0"; volcID(i)="1804-03-"; volcName(i)="Sao Vicente                    "; 
      i=i+1; volcLat(i)=  7.000; volcLon(i)=338.170; volcElev(i)=-1415;  volcLoc(i)="Atlantic-C                    "; 
             volcESP_Code(i)="U0"; volcID(i)="1805-01="; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=  4.200; volcLon(i)=338.550; volcElev(i)=-2900;  volcLoc(i)="Atlantic-C                    "; 
             volcESP_Code(i)="U0"; volcID(i)="1805-02="; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= -0.720; volcLon(i)=339.470; volcElev(i)=-1528;  volcLoc(i)="Atlantic-C                    "; 
             volcESP_Code(i)="U0"; volcID(i)="1805-03="; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= -3.500; volcLon(i)=335.500; volcElev(i)=-5300;  volcLoc(i)="Atlantic-C                    "; 
             volcESP_Code(i)="U0"; volcID(i)="1805-04="; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)= -7.950; volcLon(i)=345.630; volcElev(i)=  858;  volcLoc(i)="Atlantic-C                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1805-05-"; volcName(i)="Ascensión                      "; 
      i=i+1; volcLat(i)=-20.514; volcLon(i)=330.669; volcElev(i)=  600;  volcLoc(i)="Atlantic-C                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1805-051"; volcName(i)="Trindade                       "; 
      i=i+1; volcLat(i)=-37.420; volcLon(i)=347.520; volcElev(i)=  365;  volcLoc(i)="Atlantic-S                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1806-011"; volcName(i)="Nightingale Island             "; 
      i=i+1; volcLat(i)=-37.092; volcLon(i)=347.720; volcElev(i)= 2060;  volcLoc(i)="Atlantic-S                    "; 
             volcESP_Code(i)="M1"; volcID(i)="1806-01="; volcName(i)="Tristan da Cunha               "; 
      i=i+1; volcLat(i)=-54.420; volcLon(i)=  3.350; volcElev(i)=  780;  volcLoc(i)="Atlantic-S                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1806-02-"; volcName(i)="Bouvet                         "; 
      i=i+1; volcLat(i)=-53.930; volcLon(i)=  5.500; volcElev(i)=    0;  volcLoc(i)="Atlantic-S                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1806-03-"; volcName(i)="Thompson Island                "; 
      i=i+1; volcLat(i)=-66.420; volcLon(i)=162.470; volcElev(i)= 1340;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1900-011"; volcName(i)="Young Island                   "; 
      i=i+1; volcLat(i)=-67.400; volcLon(i)=164.830; volcElev(i)= 1167;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1900-012"; volcName(i)="Sturge Island                  "; 
      i=i+1; volcLat(i)=-72.670; volcLon(i)=165.500; volcElev(i)= 3040;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1900-013"; volcName(i)="Pleiades, The                  "; 
      i=i+1; volcLat(i)=-73.450; volcLon(i)=164.580; volcElev(i)= 2987;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1900-014"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=-74.350; volcLon(i)=164.700; volcElev(i)= 2732;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1900-015"; volcName(i)="Melbourne                      "; 
      i=i+1; volcLat(i)=-76.830; volcLon(i)=163.000; volcElev(i)= -500;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="U0"; volcID(i)="1900-016"; volcName(i)="Unnamed                        "; 
      i=i+1; volcLat(i)=-66.780; volcLon(i)=163.250; volcElev(i)= 1239;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1900-01="; volcName(i)="Buckle Island                  "; 
      i=i+1; volcLat(i)=-78.250; volcLon(i)=163.330; volcElev(i)= 3000;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1900-021"; volcName(i)="Royal Society Range            "; 
      i=i+1; volcLat(i)=-76.050; volcLon(i)=224.000; volcElev(i)= 3478;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1900-022"; volcName(i)="Berlin                         "; 
      i=i+1; volcLat(i)=-75.800; volcLon(i)=227.670; volcElev(i)= 2978;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1900-023"; volcName(i)="Andrus                         "; 
      i=i+1; volcLat(i)=-77.170; volcLon(i)=233.120; volcElev(i)= 3292;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1900-024"; volcName(i)="Waesche                        "; 
      i=i+1; volcLat(i)=-73.430; volcLon(i)=233.330; volcElev(i)= 3110;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1900-025"; volcName(i)="Siple                          "; 
      i=i+1; volcLat(i)=-75.800; volcLon(i)=244.170; volcElev(i)= 3595;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1900-026"; volcName(i)="Toney Mountain                 "; 
      i=i+1; volcLat(i)=-76.280; volcLon(i)=247.920; volcElev(i)= 3460;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1900-027"; volcName(i)="Takahe                         "; 
      i=i+1; volcLat(i)=-74.330; volcLon(i)=260.580; volcElev(i)=  749;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1900-028"; volcName(i)="Hudson Mountains               "; 
      i=i+1; volcLat(i)=-68.850; volcLon(i)=269.420; volcElev(i)= 1640;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1900-029"; volcName(i)="Peter I Island                 "; 
      i=i+1; volcLat(i)=-77.530; volcLon(i)=167.170; volcElev(i)= 3794;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="S1"; volcID(i)="1900-02="; volcName(i)="Erebus                         "; 
      i=i+1; volcLat(i)=-62.100; volcLon(i)=302.070; volcElev(i)=  180;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1900-031"; volcName(i)="Penguin Island                 "; 
      i=i+1; volcLat(i)=-62.970; volcLon(i)=299.350; volcElev(i)=  576;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="S2"; volcID(i)="1900-03="; volcName(i)="Deception Island               "; 
      i=i+1; volcLat(i)=-63.580; volcLon(i)=304.230; volcElev(i)=  353;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1900-041"; volcName(i)="Paulet                         "; 
      i=i+1; volcLat(i)=-62.050; volcLon(i)=303.250; volcElev(i)=  240;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1900-04="; volcName(i)="Bridgeman Island               "; 
      i=i+1; volcLat(i)=-65.030; volcLon(i)=299.950; volcElev(i)=  368;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="M0"; volcID(i)="1900-05="; volcName(i)="Seal Nunataks Group            "; 
      i=i+1; volcLat(i)=-59.450; volcLon(i)=332.630; volcElev(i)= 1075;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1900-07="; volcName(i)="Thule Islands                  "; 
      i=i+1; volcLat(i)=-58.420; volcLon(i)=333.670; volcElev(i)= 1370;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1900-081"; volcName(i)="Montagu Island                 "; 
      i=i+1; volcLat(i)=-59.030; volcLon(i)=333.420; volcElev(i)= 1100;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="S1"; volcID(i)="1900-08="; volcName(i)="Bristol Island                 "; 
      i=i+1; volcLat(i)=-57.780; volcLon(i)=333.550; volcElev(i)=  990;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="S1"; volcID(i)="1900-09="; volcName(i)="Michael                        "; 
      i=i+1; volcLat(i)=-57.080; volcLon(i)=333.330; volcElev(i)=  550;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1900-10="; volcName(i)="Candlemas Island               "; 
      i=i+1; volcLat(i)=-56.700; volcLon(i)=332.850; volcElev(i)= 1005;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1900-11="; volcName(i)="Hodson                         "; 
      i=i+1; volcLat(i)=-56.670; volcLon(i)=331.870; volcElev(i)=  190;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1900-12="; volcName(i)="Leskov Island                  "; 
      i=i+1; volcLat(i)=-56.300; volcLon(i)=332.430; volcElev(i)=  551;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1900-13="; volcName(i)="Zavodovski                     "; 
      i=i+1; volcLat(i)=-55.920; volcLon(i)=331.920; volcElev(i)=  -27;  volcLoc(i)="Antarctica                    "; 
             volcESP_Code(i)="S0"; volcID(i)="1900-14-"; volcName(i)="Protector Shoal                "; 

      end subroutine VotW_v12
#endif

      end module VotW_ESP
!##############################################################################

