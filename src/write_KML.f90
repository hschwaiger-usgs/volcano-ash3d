!##############################################################################
!
!  Ash3d_KML_IO module
!
!  This module manages all output to kml files
!
!      subroutine Set_OutVar_Specs
!      subroutine OpenFile_KML
!      subroutine Write_2D_KML
!      subroutine Write_PointData_Airports_KML
!      subroutine Close_KML
!      subroutine PlotModelBoundary
!
!##############################################################################

      module Ash3d_KML_IO

      use precis_param

      use io_units

      implicit none

        ! Set everything to private by default
      private

        ! Publicly available subroutines/functions
      public OpenFile_KML,              &
             Close_KML,                 &
             Write_2D_KML,              &
             Write_PointData_Airports_KML,&
             Set_OutVar_Specs

        ! Publicly available variables

      integer, parameter :: nvars      = 10  ! Number of output variables with style profiles
      integer, parameter :: max_nclrmp = 11   ! Max number of colormap points

      ! Most variables are written to a file with the generic writers:
      !  Set_OutVar_Specs
      !  OpenFile_KML(ivar)
      !  Write_2D_KML(ivar,data,groundflag,TS_flag)
      !  Close_KML(ivar,TS_flag)
      !ivar = 1 :: cloud concentration (mg/m3)
      !ivar = 2 :: cloud height (top)  (km)
      !ivar = 3 :: cloud height (bot)  (km)
      !ivar = 4 :: cloud load          (tonnes/km2)
      !ivar = 5 :: cloud arrival time  (hours)
      !ivar = 6 :: cloud reflectivity  (dBz)
      !ivar = 7 :: deposit             (mm)
      !ivar = 8 :: deposit (NWS)       (inches)
      !ivar = 9 :: deposit time        (hours)
      !ivar =10 :: topography          (km)

      ! Note file ash_arrivaltimes_airports.kml is written by the custom
      ! subroutine Write_PointData_Airports_KML

      character(len=30),dimension(nvars           ),public :: KMZ_filename
      character(len=30),dimension(nvars           ),public :: KML_filename
      character(len=5) ,dimension(nvars           ) :: KML_units
      integer          ,dimension(nvars           ) :: KML_fid
      integer          ,dimension(nvars           ) :: KML_n_clrmp
      real(kind=ip)    ,dimension(nvars,max_nclrmp) :: KML_color_map
      character(len=9) ,dimension(nvars,max_nclrmp) :: KML_Styles
      character(len=6) ,dimension(nvars,max_nclrmp) :: KML_Colors
      character(len=30),dimension(nvars           ) :: KML_description
      character(len=30),dimension(nvars           ) :: KML_legend
      character(len=3) ,dimension(nvars           ) :: KML_overlayX
      character(len=3) ,dimension(nvars           ) :: KML_overlayY
      character(len=3) ,dimension(nvars           ) :: KML_screenX
      character(len=3) ,dimension(nvars           ) :: KML_screenY
      character(len=3) ,dimension(nvars           ) :: KML_sizeX
      character(len=3) ,dimension(nvars           ) :: KML_sizeY
      character(len=13),dimension(nvars           ) :: KML_AltMode

      contains
      !------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Set_OutVar_Specs
!
!  Called from: output_results and Ash3d_PostProc.F90
!  Arguments:
!    none
!
!  This subroutine essentially initializes variables local to the Output_KML module
!  Many of these variables are slight modifications of varibles in Output_Vars
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Set_OutVar_Specs

      integer :: ivar

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Set_OutVar_Specs"
      endif;enddo

      ivar = 1 ! cloud concentration
      KMZ_filename(ivar)      = 'CloudConcentration.kmz       '
      KML_filename(ivar)      = 'CloudConcentration.kml       '
      KML_units(ivar)         = 'mg/m3'
      KML_fid(ivar)           = fid_kmlbase + ivar
      KML_n_clrmp(ivar)       = 8
      KML_color_map(ivar,:) = (/ 0.1_ip,   0.3_ip,   1.0_ip, 2.0_ip, 10.0_ip, &
                             30.0_ip, 100.0_ip, 300.0_ip, 0.0_ip, 0.0_ip, 0.0_ip/)
      KML_Styles(ivar,1) = 'pale_blue';  KML_Colors(ivar,1) = 'ffe5e5'
      KML_Styles(ivar,2) = 'lite_blue';  KML_Colors(ivar,2) = 'ffcccc'
      KML_Styles(ivar,3) = 'med__blue';  KML_Colors(ivar,3) = 'ffb2b2'
      KML_Styles(ivar,4) = 'pale_pink';  KML_Colors(ivar,4) = 'ff99ff'
      KML_Styles(ivar,5) = 'lite_pink';  KML_Colors(ivar,5) = 'ff7fff'
      KML_Styles(ivar,6) = 'med__pink';  KML_Colors(ivar,6) = 'ff66ff'
      KML_Styles(ivar,7) = 'strg_pink';  KML_Colors(ivar,7) = 'ff4cff'
      KML_Styles(ivar,8) = 'dark_pink';  KML_Colors(ivar,8) = 'ff33ff'
      KML_Styles(ivar,9) = '         ';  KML_Colors(ivar,9) = '      '
      KML_Styles(ivar,10) = '         ';  KML_Colors(ivar,10) = '      '
      KML_Styles(ivar,11) = '         ';  KML_Colors(ivar,11) = '      '
      KML_description(ivar)   = 'cloud concentration shading  '
      KML_legend(ivar)        = 'concentration_legend.png     '
      KML_overlayX(ivar)      = '0.0'
      KML_overlayY(ivar)      = '0.0'
      KML_screenX(ivar)       = '0.0'
      KML_screenY(ivar)       = '0.1'
      KML_sizeX(ivar)         = ' 77'
      KML_sizeY(ivar)         = '305'
      KML_AltMode(ivar)       = 'absolute'

      ivar = 2 ! cloud height (top)
      KMZ_filename(ivar)      = 'CloudHeight.kmz              '
      KML_filename(ivar)      = 'CloudHeight.kml              '
      KML_units(ivar)         = ' km  '
      KML_fid(ivar)           = fid_kmlbase + ivar
      KML_n_clrmp(ivar)       = 9
      KML_color_map(ivar,:) = (/0.24_ip,  3.0_ip,  6.0_ip, 10.0_ip, 13.0_ip, &
                             16.0_ip, 20.0_ip, 25.0_ip, 30.0_ip, 0.0_ip, 0.0_ip/)
      KML_Styles(ivar,1) = '00.-03.km';  KML_Colors(ivar,1) = '800080'
      KML_Styles(ivar,2) = '03.-06.km';  KML_Colors(ivar,2) = 'ff0000'
      KML_Styles(ivar,3) = '06.-10.km';  KML_Colors(ivar,3) = 'ff8000'
      KML_Styles(ivar,4) = '10.-13.km';  KML_Colors(ivar,4) = 'ffff00'
      KML_Styles(ivar,5) = '13.-16.km';  KML_Colors(ivar,5) = '80ff80'
      KML_Styles(ivar,6) = '16.-20.km';  KML_Colors(ivar,6) = '00ffff'
      KML_Styles(ivar,7) = '20.-25.km';  KML_Colors(ivar,7) = '0080ff'
      KML_Styles(ivar,8) = '25.-30.km';  KML_Colors(ivar,8) = '0000ff'
      KML_Styles(ivar,9) = '>>>30.0km';  KML_Colors(ivar,9) = '000080'
      KML_Styles(ivar,10)= '         ';  KML_Colors(ivar,10)= '      '
      KML_Styles(ivar,11)= '         ';  KML_Colors(ivar,11)= '      '
      KML_description(ivar)   = 'cloud height shading         '
      KML_legend(ivar)        = 'CloudHeight_hsv.png          '
      KML_overlayX(ivar)      = '1.0'
      KML_overlayY(ivar)      = '0.0'
      KML_screenX(ivar)       = '1.0'
      KML_screenY(ivar)       = '0.4'
      KML_sizeX(ivar)         = '150'
      KML_sizeY(ivar)         = '305'
      KML_AltMode(ivar)       = 'absolute'

      ivar = 3 ! cloud height (bot)
      KMZ_filename(ivar)      = 'CloudBottom.kmz              '
      KML_filename(ivar)      = 'CloudBottom.kml              '
      KML_units(ivar)         = ' km  '
      KML_fid(ivar)           = fid_kmlbase + ivar
      KML_n_clrmp(ivar)       = 9
      KML_color_map(ivar,:) = (/0.24_ip,  3.0_ip,  6.0_ip, 10.0_ip, 13.0_ip, &
                             16.0_ip, 20.0_ip, 25.0_ip, 30.0_ip, 0.0_ip, 0.0_ip/)
      KML_Styles(ivar,1) = '00.-03.km';  KML_Colors(ivar,1) = '800080'
      KML_Styles(ivar,2) = '03.-06.km';  KML_Colors(ivar,2) = 'ff0000'
      KML_Styles(ivar,3) = '06.-10.km';  KML_Colors(ivar,3) = 'ff8000'
      KML_Styles(ivar,4) = '10.-13.km';  KML_Colors(ivar,4) = 'ffff00'
      KML_Styles(ivar,5) = '13.-16.km';  KML_Colors(ivar,5) = '80ff80'
      KML_Styles(ivar,6) = '16.-20.km';  KML_Colors(ivar,6) = '00ffff'
      KML_Styles(ivar,7) = '20.-25.km';  KML_Colors(ivar,7) = '0080ff'
      KML_Styles(ivar,8) = '25.-30.km';  KML_Colors(ivar,8) = '0000ff'
      KML_Styles(ivar,9) = '>>>30.0km';  KML_Colors(ivar,9) = '000080'
      KML_Styles(ivar,10)= '         ';  KML_Colors(ivar,10)= '      '
      KML_Styles(ivar,11)= '         ';  KML_Colors(ivar,11)= '      '
      KML_description(ivar)   = 'cloud height shading         '
      KML_legend(ivar)        = 'CloudHeight_hsv.png          '
      KML_overlayX(ivar)      = '1.0'
      KML_overlayY(ivar)      = '0.0'
      KML_screenX(ivar)       = '1.0'
      KML_screenY(ivar)       = '0.4'
      KML_sizeX(ivar)         = '150'
      KML_sizeY(ivar)         = '305'
      KML_AltMode(ivar)       = 'absolute'

      ivar = 4 ! cloud load
      KMZ_filename(ivar)      = 'CloudLoad.kmz                '
      KML_filename(ivar)      = 'CloudLoad.kml                '
      KML_units(ivar)         = 'T/km2'
      KML_fid(ivar)           = fid_kmlbase + ivar
      KML_n_clrmp(ivar)       = 9
      KML_color_map(ivar,:) = (/ 0.2_ip,   1.0_ip,   2.0_ip,    5.0_ip, 10.0_ip, &
                             30.0_ip, 100.0_ip, 300.0_ip, 1000.0_ip, 0.0_ip, 0.0_ip/)
      KML_Styles(ivar,1) = '0.20-1.00';  KML_Colors(ivar,1) = '800080'
      KML_Styles(ivar,2) = '1.00-2.00';  KML_Colors(ivar,2) = 'ff0000'
      KML_Styles(ivar,3) = '2.00-5.00';  KML_Colors(ivar,3) = 'ff8000'
      KML_Styles(ivar,4) = '5.00-10.0';  KML_Colors(ivar,4) = 'ffff00'
      KML_Styles(ivar,5) = '10.0-30.0';  KML_Colors(ivar,5) = '80ff80'
      KML_Styles(ivar,6) = '30.0-100.';  KML_Colors(ivar,6) = '00ffff'
      KML_Styles(ivar,7) = '100.-300.';  KML_Colors(ivar,7) = '0080ff'
      KML_Styles(ivar,8) = '300.0--1k';  KML_Colors(ivar,8) = '0000ff'
      KML_Styles(ivar,9) = '>>>1000.0';  KML_Colors(ivar,9) = '000080'
      KML_Styles(ivar,10)= '         ';  KML_Colors(ivar,10)= '      '
      KML_Styles(ivar,11)= '         ';  KML_Colors(ivar,11)= '      '
      KML_description(ivar)   = 'cloud load shading           '
      KML_legend(ivar)        = 'CloudLoad_hsv.png            '
      KML_overlayX(ivar)      = '1.0'
      KML_overlayY(ivar)      = '0.0'
      KML_screenX(ivar)       = '1.0'
      KML_screenY(ivar)       = '0.1'
      KML_sizeX(ivar)         = '150'
      KML_sizeY(ivar)         = '305'
      KML_AltMode(ivar)       = 'absolute'

      ivar = 5 ! cloud arrival time
      KMZ_filename(ivar)      = 'cloud_arrivaltimes_hours.kmz'
      KML_filename(ivar)      = 'cloud_arrivaltimes_hours.kml'
      KML_units(ivar)         = ' hrs '
      KML_fid(ivar)           = fid_kmlbase + ivar
      KML_n_clrmp(ivar)       = 9
      KML_color_map(ivar,:) = (/ 0.0_ip,  3.0_ip,  6.0_ip,  9.0_ip, 12.0_ip, &
                             15.0_ip, 18.0_ip, 24.0_ip, 36.0_ip, 0.0_ip, 0.0_ip/)
      KML_Styles(ivar,1) = '00-03_hrs';  KML_Colors(ivar,1) = '800080'
      KML_Styles(ivar,2) = '03-06_hrs';  KML_Colors(ivar,2) = 'ff0000'
      KML_Styles(ivar,3) = '06-09_hrs';  KML_Colors(ivar,3) = 'ff8000'
      KML_Styles(ivar,4) = '09-12_hrs';  KML_Colors(ivar,4) = 'ffff00'
      KML_Styles(ivar,5) = '12-15_hrs';  KML_Colors(ivar,5) = '80ff80'
      KML_Styles(ivar,6) = '15-18_hrs';  KML_Colors(ivar,6) = '00ffff'
      KML_Styles(ivar,7) = '18-24_hrs';  KML_Colors(ivar,7) = '0080ff'
      KML_Styles(ivar,8) = '24-36_hrs';  KML_Colors(ivar,8) = '0000ff'
      KML_Styles(ivar,9) = '>>>36_hrs';  KML_Colors(ivar,9) = '000080'
      KML_Styles(ivar,10)= '         ';  KML_Colors(ivar,10)= '      '
      KML_Styles(ivar,11)= '         ';  KML_Colors(ivar,11)= '      '
      KML_description(ivar)   = 'shades for cloud arrival time'
      KML_legend(ivar)        = 'CloudLoad_hsv.png            '
      KML_legend(ivar)        = 'cloud_arrival_time.png       '
      KML_overlayX(ivar)      = '1.0'
      KML_overlayY(ivar)      = '0.0'
      KML_screenX(ivar)       = '1.0'
      KML_screenY(ivar)       = '0.3'
      KML_sizeX(ivar)         = '150'
      KML_sizeY(ivar)         = '305'
      KML_AltMode(ivar)       = 'clampToGround'

      ivar = 6 ! cloud reflectivity
      KMZ_filename(ivar)      = 'reflectivity.kmz             '
      KML_filename(ivar)      = 'reflectivity.kml             '
      KML_units(ivar)         = ' dBZ '
      KML_fid(ivar)           = fid_kmlbase + ivar
      KML_n_clrmp(ivar)       = 9
      KML_color_map(ivar,:) = (/-20.0_ip, -10.0_ip,  0.0_ip, 10.0_ip, 20.0_ip, &
                           30.0_ip,  40.0_ip, 50.0_ip, 60.0_ip, 0.0_ip, 0.0_ip/)
      KML_Styles(ivar,1) = '-20_to-10';  KML_Colors(ivar,1) = '800080'
      KML_Styles(ivar,2) = '-10_to000';  KML_Colors(ivar,2) = 'ff0000'
      KML_Styles(ivar,3) = '000_to+10';  KML_Colors(ivar,3) = 'ff8000'
      KML_Styles(ivar,4) = '+10_to+20';  KML_Colors(ivar,4) = 'ffff00'
      KML_Styles(ivar,5) = '+20_to+30';  KML_Colors(ivar,5) = '80ff80'
      KML_Styles(ivar,6) = '+30_to+40';  KML_Colors(ivar,6) = '00ffff'
      KML_Styles(ivar,7) = '+40_to+50';  KML_Colors(ivar,7) = '0080ff'
      KML_Styles(ivar,8) = '+50_to+60';  KML_Colors(ivar,8) = '0000ff'
      KML_Styles(ivar,9) = '>>>>>>+60';  KML_Colors(ivar,9) = '000080'
      KML_Styles(ivar,10)= '         ';  KML_Colors(ivar,10)= '      '
      KML_Styles(ivar,11)= '         ';  KML_Colors(ivar,11)= '      '
      KML_description(ivar)   = 'radar reflectivity shading   '
      KML_legend(ivar)        = 'cloud_dbZ_hsv.png            '
      KML_overlayX(ivar)      = '1.0'
      KML_overlayY(ivar)      = '0.0'
      KML_screenX(ivar)       = '1.0'
      KML_screenY(ivar)       = '0.1'
      KML_sizeX(ivar)         = '150'
      KML_sizeY(ivar)         = '305'
      KML_AltMode(ivar)       = 'absolute'

      ivar = 7 ! deposit
      KMZ_filename(ivar)      = 'deposit_thickness_mm.kmz     '
      KML_filename(ivar)      = 'deposit_thickness_mm.kml     '
      KML_units(ivar)         = '  mm '
      KML_fid(ivar)           = fid_kmlbase + ivar
      KML_n_clrmp(ivar)       = 11
      KML_color_map(ivar,:) = (/ 0.01_ip, 0.03_ip, 0.1_ip,   0.3_ip,   1.0_ip,    3.0_ip, 10.0_ip,&
                             30.0_ip, 100.0_ip, 300.0_ip, 1000.0_ip/)
      KML_Styles(ivar,1) = '.01-.03mm';  KML_Colors(ivar,1) = '69ded6'
      KML_Styles(ivar,2) = '.03-0.1mm';  KML_Colors(ivar,2) = '71a7f9'
      KML_Styles(ivar,3) = '0.1-0.3mm';  KML_Colors(ivar,3) = '800080'
      KML_Styles(ivar,4) = '0.3-1.0mm';  KML_Colors(ivar,4) = 'ff0000'
      KML_Styles(ivar,5) = '1.0-3.0mm';  KML_Colors(ivar,5) = 'ff8000'
      KML_Styles(ivar,6) = '3.0-10.mm';  KML_Colors(ivar,6) = 'ffff00'
      KML_Styles(ivar,7) = '10.-30.mm';  KML_Colors(ivar,7) = '80ff80'
      KML_Styles(ivar,8) = '30.-100mm';  KML_Colors(ivar,8) = '00ffff'
      KML_Styles(ivar,9) = '100-300mm';  KML_Colors(ivar,9) = '0080ff'
      KML_Styles(ivar,10) = '300-1k_mm';  KML_Colors(ivar,10) = '0000ff'
      KML_Styles(ivar,11) = '>>>>1k_mm';  KML_Colors(ivar,11) = '000080'
      KML_description(ivar)   = 'Ash thickness shades         '
      KML_legend(ivar)        = 'deposit_thickness_hsv.png    '
      KML_overlayX(ivar)      = '1.0'
      KML_overlayY(ivar)      = '0.0'
      KML_screenX(ivar)       = '1.0'
      KML_screenY(ivar)       = '0.1'
      KML_sizeX(ivar)         = '150'
      KML_sizeY(ivar)         = '357'
      KML_AltMode(ivar)       = 'clampToGround'

      ivar = 8 ! deposit (NWS)
      KMZ_filename(ivar)      = 'deposit_thickness_inches.kmz  '
      KML_filename(ivar)      = 'deposit_thickness_inches.kml  '
      KML_units(ivar)         = '  in.'
      KML_fid(ivar)           = fid_kmlbase + ivar
      KML_n_clrmp(ivar)       = 3
      KML_color_map(ivar,:) = (/ 0.00394_ip, 0.0315_ip, 0.236_ip, 0.0_ip, 0.0_ip,&
                              0.0_ip, 0.0_ip, 0.0_ip, 0.0_ip, 0.0_ip, 0.0_ip/)
      KML_Styles(ivar,1) = '0.1-0.8mm';  KML_Colors(ivar,1) = 'ffff00'
      KML_Styles(ivar,2) = '0.8-6.0mm';  KML_Colors(ivar,2) = '00ffff'
      KML_Styles(ivar,3) = '>>>>6.0mm';  KML_Colors(ivar,3) = '0000ff'
      KML_Styles(ivar,4) = '         ';  KML_Colors(ivar,4) = '      '
      KML_Styles(ivar,5) = '         ';  KML_Colors(ivar,5) = '      '
      KML_Styles(ivar,6) = '         ';  KML_Colors(ivar,6) = '      '
      KML_Styles(ivar,7) = '         ';  KML_Colors(ivar,7) = '      '
      KML_Styles(ivar,8) = '         ';  KML_Colors(ivar,8) = '      '
      KML_Styles(ivar,9) = '         ';  KML_Colors(ivar,9) = '      '
      KML_Styles(ivar,10)= '         ';  KML_Colors(ivar,10)= '      '
      KML_Styles(ivar,11)= '         ';  KML_Colors(ivar,11)= '      '
      KML_description(ivar)   = 'Ash thickness shades         '
      KML_legend(ivar)        = 'GE_legend_dep_nws.png        '
      KML_overlayX(ivar)      = '1.0'
      KML_overlayY(ivar)      = '0.0'
      KML_screenX(ivar)       = '1.0'
      KML_screenY(ivar)       = '0.1'
      KML_sizeX(ivar)         = '150'
      KML_sizeY(ivar)         = '316'
      KML_AltMode(ivar)       = 'clampToGround'

      ivar = 9 ! deposit time
      KMZ_filename(ivar)      = 'ashfall_arrivaltimes_hours.kmz'
      KML_filename(ivar)      = 'ashfall_arrivaltimes_hours.kml'
      KML_units(ivar)         = ' hrs '
      KML_fid(ivar)           = fid_kmlbase + ivar
      KML_n_clrmp(ivar)       = 9
      KML_color_map(ivar,:) = (/ 0.0_ip,  3.0_ip,  6.0_ip,  9.0_ip, 12.0_ip, &
                             15.0_ip, 18.0_ip, 24.0_ip, 36.0_ip, 0.0_ip, 0.0_ip/)
      KML_Styles(ivar,1) = '00-03_hrs';  KML_Colors(ivar,1) = '800080'
      KML_Styles(ivar,2) = '03-06_hrs';  KML_Colors(ivar,2) = 'ff0000'
      KML_Styles(ivar,3) = '06-09_hrs';  KML_Colors(ivar,3) = 'ff8000'
      KML_Styles(ivar,4) = '09-12_hrs';  KML_Colors(ivar,4) = 'ffff00'
      KML_Styles(ivar,5) = '12-15_hrs';  KML_Colors(ivar,5) = '80ff80'
      KML_Styles(ivar,6) = '15-18_hrs';  KML_Colors(ivar,6) = '00ffff'
      KML_Styles(ivar,7) = '18-24_hrs';  KML_Colors(ivar,7) = '0080ff'
      KML_Styles(ivar,8) = '24-36_hrs';  KML_Colors(ivar,8) = '0000ff'
      KML_Styles(ivar,9) = '>>>36_hrs';  KML_Colors(ivar,9) = '000080'
      KML_Styles(ivar,10)= '         ';  KML_Colors(ivar,10)= '      '
      KML_Styles(ivar,11)= '         ';  KML_Colors(ivar,11)= '      '
      KML_description(ivar)   = 'shades for dep arrival time  '
      KML_legend(ivar)        = 'deposit_arrival_time.png     '
      KML_overlayX(ivar)      = '1.0'
      KML_overlayY(ivar)      = '0.0'
      KML_screenX(ivar)       = '1.0'
      KML_screenY(ivar)       = '0.3'
      KML_sizeX(ivar)         = '150'
      KML_sizeY(ivar)         = '316'
      KML_AltMode(ivar)       = 'clampToGround'

      ivar = 10 ! topography
      KMZ_filename(ivar)      = 'Topography.kmz                '
      KML_filename(ivar)      = 'Topography.kml                '
      KML_units(ivar)         = '  km '
      KML_fid(ivar)           = fid_kmlbase + ivar
      KML_n_clrmp(ivar)       = 8
      KML_color_map(ivar,:) = (/ 1.0_ip,   2.0_ip,   3.0_ip,    4.0_ip, 5.0_ip,&
                                 6.0_ip, 7.0_ip, 8.0_ip, 0.0_ip, 0.0_ip, 0.0_ip/)
      KML_Styles(ivar,1) = 'pale_blue';  KML_Colors(ivar,1) = 'ffe5e5'
      KML_Styles(ivar,2) = 'lite_blue';  KML_Colors(ivar,2) = 'ffcccc'
      KML_Styles(ivar,3) = 'med__blue';  KML_Colors(ivar,3) = 'ffb2b2'
      KML_Styles(ivar,4) = 'pale_pink';  KML_Colors(ivar,4) = 'ff99ff'
      KML_Styles(ivar,5) = 'lite_pink';  KML_Colors(ivar,5) = 'ff7fff'
      KML_Styles(ivar,6) = 'med__pink';  KML_Colors(ivar,6) = 'ff66ff'
      KML_Styles(ivar,7) = 'strg_pink';  KML_Colors(ivar,7) = 'ff4cff'
      KML_Styles(ivar,8) = 'dark_pink';  KML_Colors(ivar,8) = 'ff33ff'
      KML_Styles(ivar,9) = '         ';  KML_Colors(ivar,9) = '      '
      KML_Styles(ivar,10)= '         ';  KML_Colors(ivar,10)= '      '
      KML_Styles(ivar,11)= '         ';  KML_Colors(ivar,11)= '      '
      KML_description(ivar)   = 'Topography shades         '
      KML_legend(ivar)        = 'deposit_thickness_hsv.png    '
      KML_overlayX(ivar)      = '1.0'
      KML_overlayY(ivar)      = '0.0'
      KML_screenX(ivar)       = '1.0'
      KML_screenY(ivar)       = '0.1'
      KML_sizeX(ivar)         = '150'
      KML_sizeY(ivar)         = '305'
      KML_AltMode(ivar)       = 'clampToGround'

      end subroutine Set_OutVar_Specs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  OpenFile_KML
!
!  Called from: output_results and Ash3d_PostProc.F90
!  Arguments:
!    ivar = ID of variable to process
!
!  This subroutine opens and initializes the KML file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine OpenFile_KML(ivar)

      use io_data,       only : &
         VolcanoName,WriteDepositTS_KML

      use mesh,          only : &
         A3d_iprojflag,A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
         A3d_k0_scale,A3d_Re,IsLatLon, &
         latLL,lonLL,latUR,lonUR,xLL,yLL,xUR,yUR !,de,dn,dx,dy

      use time_data,     only : &
         BaseYear,useLeap,SimStartHour,OutputOffset

      use Source,        only : &
         neruptions, e_StartTime,e_Duration,e_Volume,e_PlumeHeight, &
         lon_volcano,lat_volcano,x_volcano,y_volcano

      use projection,    only : &
           PJ_proj_inv

      integer,intent(in) :: ivar

      character (len=13) :: yyyymmddhh
      character (len=2)  :: opacity
      real(kind=ip)      :: xleft, xright, ybottom, ytop
      real(kind=ip)      :: longLL,longUR
      real(kind=ip)      :: lattLL,lattUR

      real(kind=dp)      :: olam,ophi ! using precision needed by libprojection

      integer,       dimension(:), allocatable :: iyear, imonth, iday
      real(kind=ip), dimension(:), allocatable :: StartHour

      integer            :: ierup
      character(len=30)  :: filename
      integer            :: fid
      integer            :: n_clrmp,icmp
      character(len=9),dimension(11)   :: Styles
      character(len=6),dimension(11)   :: Colors
      character(len=30)  :: description
      character(len=30)  :: legend
      character(len=3)   :: overlayX
      character(len=3)   :: overlayY
      character(len=3)   :: screenX
      character(len=3)   :: screenY
      character(len=3)   :: sizeX
      character(len=3)   :: sizeY
      integer            :: iostatus
      character(len=120) :: iomessage
      character(len= 50) :: linebuffer050 
      character(len= 80) :: linebuffer080

      INTERFACE
        character (len=13) function HS_yyyymmddhh_since(HoursSince,byear,useLeaps)
          real(kind=8)               ::  HoursSince
          integer                    ::  byear
          logical                    ::  useLeaps
        end function HS_yyyymmddhh_since
      END INTERFACE

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine OpenFile_KML"
      endif;enddo

      filename    = KML_filename(ivar)
      fid         = KML_fid(ivar)
      n_clrmp     = KML_n_clrmp(ivar)
      Styles(:)   = KML_Styles(ivar,:)
      Colors(:)   = KML_Colors(ivar,:)
      description = KML_description(ivar)
      legend      = KML_legend(ivar)
      overlayX    = KML_overlayX(ivar)
      overlayY    = KML_overlayY(ivar)
      screenX     = KML_screenX(ivar)
      screenY     = KML_screenY(ivar)
      sizeX       = KML_sizeX(ivar)
      sizeY       = KML_sizeY(ivar)

      opacity = '80'
      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Opening KML file ",trim(adjustl(filename))
      endif;enddo
      open(fid,file=trim(adjustl(filename)),status='replace',action='write',err=2500)

      write(fid,1)                 ! write file header (35 lines)

      ! write StyleMap entries (each call is 10 lines)
      do icmp = 1,n_clrmp
        write(fid,2) Styles(icmp), Styles(icmp), Styles(icmp)
      enddo

      ! write highlighted styles (each call is 19 lines)
      do icmp = 1,n_clrmp
        write(fid,3) Styles(icmp), opacity, Colors(icmp)
      enddo

      ! write normal styles (each call is 19 lines)
      do icmp = 1,n_clrmp
        write(fid,4) Styles(icmp), opacity, Colors(icmp)
      enddo

      ! write legend (20 lines starting with ScreenOverlay)
      write(fid,5)description,legend,overlayX,overlayY,screenX,screenY,sizeX,sizeY

      ! Plot model region
      if (IsLatLon) then
        longLL  = lonLL
        longUR  = lonUR
        lattLL  = latLL
        lattUR  = latUR
        call PlotModelBoundary(longLL,longUR,lattLL,lattUR,fid)
      else
        xleft   = xLL
        xright  = xUR
        ybottom = yLL
        ytop    = yUR
        call PlotModelBoundary(xleft,xright,ybottom,ytop,fid)
      endif ! IsLatLon

      allocate(iyear(neruptions))
      allocate(imonth(neruptions))
      allocate(iday(neruptions))
      allocate(StartHour(neruptions))

      write(fid,7) VolcanoName, VolcanoName

      ! write table of ESP's to the volcano placemark
      do ierup=1,neruptions
        yyyymmddhh = HS_yyyymmddhh_since(e_StartTime(ierup)+SimStartHour+OutputOffset,&
                                         BaseYear,useLeap)
        read(yyyymmddhh,100,iostat=iostatus,iomsg=iomessage) &
                             iyear(ierup),imonth(ierup),iday(ierup), &
                             StartHour(ierup)
        linebuffer080 = yyyymmddhh
        linebuffer050 = "Reading date from date_string (write_KML)"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        write(fid,8) ierup,iyear(ierup),imonth(ierup),iday(ierup), &
                             StartHour(ierup), &
                             e_PlumeHeight(ierup),e_Duration(ierup),e_Volume(ierup)
100     format(i4,i2,i2,f5.2)
      enddo

      ! Plot volcano
      if (.not.IsLatLon) then                        !get lon_volcano and lat_volcano
        call PJ_proj_inv(real(x_volcano,kind=dp),real(y_volcano,kind=dp), &
                   A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
                   A3d_k0_scale,A3d_Re, &
                   olam,ophi)
        lon_volcano = real(olam,kind=ip)
        lat_volcano = real(ophi,kind=ip)
      endif
      if (lon_volcano.gt.180.0_ip) then
        write(fid,9) lon_volcano-360.0_ip, lat_volcano
      else
        write(fid,9) lon_volcano, lat_volcano
      endif

      ! Create forecast folder only for the files that have time steps
      if ((ivar.eq.5).or. &                                    ! cloud arrival time
         ((ivar.eq.7).and.(.not.WriteDepositTS_KML)).or. &     ! deposit
         ((ivar.eq.8).and.(.not.WriteDepositTS_KML)).or. &     ! deposit_NWS
          (ivar.eq.9)) then                                    ! deposit arrival time
        continue
      else
        write(fid,6)                                       ! create folder of forecasts
      endif

      deallocate(iyear,imonth,iday,StartHour)

      return

      ! Error traps
2500  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),20)
      endif;enddo
      stop 1

      ! Format statements
1     format('<?xml version="1.0" encoding="UTF-8"?>',/, &
             '<kml xmlns="http://www.opengis.net/kml/2.2">',/, &
             '<Document>',/, &
             '    <StyleMap id="VolcanoMarker">',/, &
             '      <Pair>',/, &
             '             <key>normal</key>',/, &
             '             <styleUrl>#VolcanoMarkerNormal</styleUrl>',/, &
             '      </Pair>',/, &
             '      <Pair>',/, &
             '             <key>highlight</key>',/, &
             '             <styleUrl>#VolcanoMarkerHighlight</styleUrl>',/, &
             '      </Pair>',/, &
             '    </StyleMap>',/, &
             '    <Style id="VolcanoMarkerNormal">',/, &
             '      <IconStyle>',/, &
             '             <scale>0.8</scale>',/, &
             '             <Icon>',/, &
             '                   <href>http://maps.google.com/mapfiles/kml/pal4/icon52.png</href>',/, &
             '             </Icon>',/, &
             '      </IconStyle>',/, &
             '      <LabelStyle>',/, &
             '             <scale>0</scale>',/, &
             '      </LabelStyle>',/, &
             '    </Style>',/, &
             '    <Style id="VolcanoMarkerHighlight">',/, &
             '      <IconStyle>',/, &
             '             <scale>1</scale>',/, &
             '             <Icon>',/, &
             '                   <href>http://maps.google.com/mapfiles/kml/pal4/icon52.png</href>',/, &
             '             </Icon>',/, &
             '      </IconStyle>',/, &
             '      <LabelStyle>',/, &
             '             <scale>1</scale>',/, &
             '      </LabelStyle>',/, &
             '    </Style>')
2     format('    <StyleMap id="',a9,'">',/, &
             '      <Pair>',/, &
             '             <key>normal</key>',/, &
             '             <styleUrl>#',a9,'Normal</styleUrl>',/, &
             '      </Pair>',/, &
             '      <Pair>',/, &
             '             <key>highlight</key>',/, &
             '             <styleUrl>#',a9,'Highlight</styleUrl>',/, &
             '      </Pair>',/, &
             '    </StyleMap>')
3     format('    <Style id="',a9,'Normal">',/, &
             '      <IconStyle>',/, &
             '             <scale>0</scale>',/, &
             '             <Icon>',/, &
             '                   <href>http://maps.google.com/mapfiles/kml/pal4/icon56.png</href>',/, &
             '             </Icon>',/, &
             '      </IconStyle>',/, &
             '      <LabelStyle>',/, &
             '             <scale>0</scale>',/, &
             '      </LabelStyle>',/, &
             '      <LineStyle>',/, &
             '        <color>ffffffff</color>',/, &
             '        <width>0</width>',/, &
             '      </LineStyle>',/, &
             '      <PolyStyle>',/, &
             '        <color>',a2,a6,'</color>',/, &
             '        <fill>1</fill>',/, &
             '      </PolyStyle>',/, &
             '    </Style>')
4     format('    <Style id="',a9,'Highlight">',/, &
             '      <IconStyle>',/, &
             '             <scale>0</scale>',/, &
             '             <Icon>',/, &
             '                   <href>http://maps.google.com/mapfiles/kml/pal4/icon56.png</href>',/, &
             '             </Icon>',/, &
             '      </IconStyle>',/, &
             '      <LabelStyle>',/, &
             '             <scale>1</scale>',/, &
             '      </LabelStyle>',/, &
             '      <LineStyle>',/, &
             '        <color>ffffffff</color>',/, &
             '        <width>5</width>',/, &
             '      </LineStyle>',/, &
             '      <PolyStyle>',/, &
             '        <color>',a2,a6,'</color>',/, &
             '        <fill>1</fill>',/, &
             '      </PolyStyle>',/, &
             '    </Style>')
5     format('    <ScreenOverlay>',/, &
             '      <name>Legend</name>',/, &
             '      <description>',a30,'</description>',/, &
             '        <Icon>',/, &
             '          <href>https://vsc-ash.wr.usgs.gov/images/',a30,'</href>',/, &
             '        </Icon>',/, &
             '      <overlayXY x="',a3,'" y="',a3,'" xunits="fraction" yunits="fraction"/>',/, &
             '      <screenXY x="',a3,'" y="',a3,'" xunits="fraction" yunits="fraction"/>',/, &
             '      <size x="',a3,'" y="',a3,'" xunits="pixels" yunits="pixels"/>',/, &
             '    </ScreenOverlay>',/, &
             '    <ScreenOverlay>',/, &
             '      <name>Warning</name>',/, &
             '      <description>USGS disclaimer</description>',/, &
             '        <Icon>',/, &
             '          <href>https://vsc-ash.wr.usgs.gov/images/USGS_warning3.png</href>',/, &
             '        </Icon>',/, &
             '      <overlayXY x="0.5" y="1.0" xunits="fraction" yunits="fraction"/>',/, &
             '      <screenXY x="0.5" y="1.0" xunits="fraction" yunits="fraction"/>',/, &
             '      <size x="337" y="68" xunits="pixels" yunits="pixels"/>',/, &
             '    </ScreenOverlay>')
6    format('   <Folder id="forecasts">',/, &
            '       <name>Forecasts</name>',/, &
            '       <open>1</open>',/, &
            '       <visibility>1</visibility>',/,&
            '       <Style>',/, &
            '           <ListStyle>',/, &
            '               <bgColor>00ffffff</bgColor>',/, &
            '           </ListStyle>',/, &
            '       </Style>')
7     format('  <Placemark>',/, &
             '    <name>',a35,'</name>',/, &
             '    <styleUrl>#VolcanoMarker</styleUrl>',/, &
             '    <gx:balloonVisibility>1</gx:balloonVisibility>',/, &
             '    <description>',/, &
             '      <![CDATA[',/, &
             '        <font size=4>',/, &
             '          <p><b>Volcano name: </b>',a35,'</p>',/, &
             '        </font>',/, &
             '        <hr />',/, &
             '        <font size=4><b>Eruption source parameters:</b></font>',/, &
             '        <table border="1" padding="0" width="400">',/, &
             '          <tr><td>Pulse #</td><td>Start year</td><td>month</td><td>day</td>',/, &
             '            <td>hour (UTC)</td><td>plume ht (km asl)</td>',/, &
             '            <td>Duration (hrs)</td><td>Volume (km3)</td></tr>')

8     format('          <tr><td>',i2,'</td><td>',i4,'</td><td>',i2,/, &
             '            </td><td>',i2,'</td><td>',f5.2,'</td><td>',f4.1,'</td>',/, &
             '            <td>',f6.2,'</td><td>',e8.3,'</td></tr>')
9     format('         </table>',/, &
             '      ]]>',/, &
             '    </description>',/, &
             '    <Point>',/, &
             '      <coordinates>',e12.6,',',e12.6,',0</coordinates>',/, &
             '    </Point>',/, &
             '  </Placemark>')
20    format(/,4x,'Error: Can''t open kml file for output.  Program stopped.')

      end subroutine OpenFile_KML

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Write_2D_KML
!
!  Called from: output_results and Ash3d_PostProc.F90
!  Arguments:
!    ivar        = ID of variable to process
!    OutVar      = 2-d array of data to write to KML
!    height_flag = flag to indicate heigh of cell (<0 : min, =0 : ground, >0 : max)
!    TS_flag     = 0 = not a time-series, 1 = time-series
!
!  This subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Write_2D_KML(ivar,OutVar,height_flag,TS_flag)
      
      use mesh,          only : &
         nxmax,nymax,A3d_iprojflag,A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2,&
         A3d_k0_scale,A3d_Re,de,dn,dx,dy,IsLatLon,&
         lon_cc_pd,lat_cc_pd,x_cc_pd,y_cc_pd

      use time_data,     only : &
         time,BaseYear,useLeap,SimStartHour,OutputOffset, &
         xmlTimeSpanStart,xmlTimeSpanEnd

      use Output_Vars,   only : &
         MinHeight,MaxHeight

      use projection,    only : &
           PJ_proj_inv

      integer      ,intent(in)  :: ivar
      real(kind=ip),intent(in)  :: OutVar(nxmax,nymax)
      integer      ,intent(in)  :: height_flag          ! <0 : min, =0 : ground, >0 : max
      integer      ,intent(in)  :: TS_flag              ! 0 = not a time-series, 1 = time-series

      real(kind=dp)       :: olam,ophi ! using precision needed by libprojection

      integer             :: i,j
      character(len=9)    :: StyleNow3
      real(kind=ip)       :: xleft,xright,ybottom,ytop
      real(kind=ip)       :: longLL,longUL,longLR,longUR,longCC
      real(kind=ip)       :: lattLL,lattUL,lattLR,lattUR,lattCC
      real(kind=ip)       :: longLL1,longUL1,longLR1,longUR1    !for polygons that cross the
      real(kind=ip)       :: longLL2,longUL2,longLR2,longUR2    !antimeridian
      real(kind=ip)       :: lattLL1,lattUL1,lattLR1,lattUR1    !for polygons that cross the
      real(kind=ip)       :: lattLL2,lattUL2,lattLR2,lattUR2    !antimeridian
      real(kind=ip)       :: longLR3,longUR3
      character (len=20)  :: xmlArrivalTime
      integer             :: height
      integer             :: fid
      integer             :: n_clrmp,icmp
      character(len=5)    :: units
      character(len=13)   :: AltMode
      real(kind=ip)   ,dimension(11)   :: color_map
      character(len=9),dimension(11)   :: Styles
      logical             :: CrossAntiMeridian     !if the polygon crosses the antimeridian

      INTERFACE
        character(len=20) function HS_xmltime(HoursSince,byear,useLeaps)
          real(kind=8)              :: HoursSince
          integer                   :: byear
          logical                   :: useLeaps
        end function HS_xmltime
      END INTERFACE

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Write_2D_KML"
      endif;enddo

      fid       = KML_fid(ivar)
      n_clrmp   = KML_n_clrmp(ivar)
      color_map = KML_color_map(ivar,:)
      Styles(:) = KML_Styles(ivar,:)
      units     = KML_units(ivar)
      AltMode   = KML_AltMode(ivar)

      xmlArrivalTime = HS_xmltime(SimStartHour + time + OutputOffset,&
                                  BaseYear,useLeap)

      StyleNow3 = 'PureWhite'

      if(TS_flag.ne.0)then
        write(fid,1) xmlArrivalTime, xmlArrivalTime,  &
                 xmlTimeSpanStart, xmlTimeSpanEnd
      else
        write(fid,15) 
      endif
      ! close folder if this is the final deposit in a deposit file

      do i=1,nxmax
        do j=1,nymax
          if (OutVar(i,j).lt.color_map(1)) cycle
          StyleNow3 = Styles(n_clrmp)
          do icmp = 1,n_clrmp-1
            if (OutVar(i,j).gt.color_map(icmp).and.&
                OutVar(i,j).le.color_map(icmp+1)) StyleNow3 = Styles(icmp)
          enddo
          if (IsLatLon) then
            longLL  = lon_cc_pd(i) - de/2.0_ip
            lattLL  = lat_cc_pd(j) - dn/2.0_ip
            longUL  = longLL
            lattUL  = lat_cc_pd(j) + dn/2.0_ip
            longLR  = lon_cc_pd(i) + de/2.0_ip
            lattLR  = lat_cc_pd(j) - dn/2.0_ip
            longUR  = longLR
            lattUR  = lat_cc_pd(j) + dn/2.0_ip
            longCC  = lon_cc_pd(i)
            lattCC  = lat_cc_pd(j)
          else
            xleft   = x_cc_pd(i) - dx/2.0_ip
            xright  = x_cc_pd(i) + dx/2.0_ip
            ybottom = y_cc_pd(j) - dy/2.0_ip
            ytop    = y_cc_pd(j) + dy/2.0_ip
            call PJ_proj_inv(real(xleft,kind=dp), real(ybottom,kind=dp),  &
                        A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
                        A3d_k0_scale,A3d_Re, &
                           olam,ophi)
            longLL = real(olam,kind=ip)
            lattLL = real(ophi,kind=ip)
            call PJ_proj_inv(real(xleft,kind=dp),    real(ytop,kind=dp),  &
                        A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
                        A3d_k0_scale,A3d_Re, &
                           olam,ophi)
            longUL = real(olam,kind=ip)
            lattUL = real(ophi,kind=ip)
            call PJ_proj_inv(real(xright,kind=dp),   real(ytop,kind=dp),  &
                        A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
                        A3d_k0_scale,A3d_Re, &
                           olam,ophi)
            longUR = real(olam,kind=ip)
            lattUR = real(ophi,kind=ip)
            call PJ_proj_inv(real(xright,kind=dp),real(ybottom,kind=dp),  &
                        A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
                        A3d_k0_scale,A3d_Re, &
                           olam,ophi)
            longLR = real(olam,kind=ip)
            lattLR = real(ophi,kind=ip)
            call PJ_proj_inv(real(x_cc_pd(i),kind=dp),real(y_cc_pd(j),kind=dp),  &
                        A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
                        A3d_k0_scale,A3d_Re, &
                           olam,ophi)
            longCC = real(olam,kind=ip)
            lattCC = real(ophi,kind=ip)
          endif ! IsLatLon
          if (longLL.gt.180.0_ip) longLL = longLL-360.0_ip
          if (longUL.gt.180.0_ip) longUL = longUL-360.0_ip
          if (longLR.gt.180.0_ip) longLR = longLR-360.0_ip
          if (longUR.gt.180.0_ip) longUR = longUR-360.0_ip
          if (longCC.gt.180.0_ip) longCC = longCC-360.0_ip
          if(height_flag.gt.0)then
            height = int(MaxHeight(i,j)*1000.0_ip)
          elseif(height_flag.lt.0)then
            height = int(MinHeight(i,j)*1000.0_ip)
          else
            height = 0
          endif
          if (LongUR.lt.LongUL) then
            CrossAntiMeridian = .true.       ! polygon crosses the antimeridian
            ! establish two polygons.  The first is left of the AM
            longLL1=longLL
            longUL1=longUL
            longLR1=179.9999_ip
            longUR1=179.9999_ip
            lattLL1=lattLL
            lattUL1=lattUL
            ! The second polygon is right of the AM
            longLL2=-179.9999_ip
            longUL2=-179.9999_ip
            longLR2=longLR
            longUR2=longUR
            lattLR2=lattLR
            lattUR2=lattUR
            ! interpolate to find lattLR1,lattUR1,lattLL2,lattUL2
            longUR3=longUR+360.0_ip
            longLR3=longLR+360.0_ip
            lattLR1=lattLL+(lattLR-lattLL)*(179.9999_ip-longLL)/(longLR3-longLL)
            lattUR1=lattUL+(lattUR-lattUL)*(179.9999_ip-longUL)/(longUR3-longUL)
            lattLL2=lattLR1
            lattUL2=lattUR1
          else
            CrossAntiMeridian = .false.
          endif
          ! write out polygon to  kml file
          if (.not.CrossAntiMeridian) then
            write(fid,2) OutVar(i,j),units,StyleNow3, &
                         longCC,lattCC,real(height,kind=4),  &
                         AltMode,                      &
                         longLL,lattLL,real(height,kind=4), &
                         longLR,lattLR,real(height,kind=4), &
                         longUR,lattUR,real(height,kind=4), &
                         longUL,lattUL,real(height,kind=4), &
                         longLL,lattLL,real(height,kind=4)
          else
            ! This fixes the antimeridian problem
            write(fid,25) OutVar(i,j),units,StyleNow3, &
                         longCC,lattCC,real(height,kind=4),  &
                         AltMode,                      &
                         longLL1,lattLL1,real(height,kind=4), &
                         longLR1,lattLR1,real(height,kind=4), &
                         longUR1,lattUR1,real(height,kind=4), &
                         longUL1,lattUL1,real(height,kind=4), &
                         longLL1,lattLL1,real(height,kind=4), &
                         longLL2,lattLL2,real(height,kind=4), &
                         longLR2,lattLR2,real(height,kind=4), &
                         longUR2,lattUR2,real(height,kind=4), &
                         longUL2,lattUL2,real(height,kind=4), &
                         longLL2,lattLL2,real(height,kind=4)
          endif
        enddo
      enddo

      write(fid,3)   ! close folder

      return
      
      ! format statements
1     format('      <Folder>',/, &
             '        <name>',a20,'</name>',/, &
             '        <open>0</open>',/, &
             '        <visibility>1</visibility>',/, &
             '        <description>Folder containing values on ',a20,'.</description>',/, &
             '        <TimeSpan>',/, &
             '           <begin>',a20,'</begin>',/, &
             '           <end>',a20,'</end>',/, &
             '        </TimeSpan>',/, &
             '        <Style>',/, &
             '             <ListStyle>',/, &
             '                   <listItemType>checkHideChildren</listItemType>',/, &
             '             </ListStyle>',/, &
             '        </Style>')
15    format('      <Folder>',/, &
             '        <name>Final</name>',/, &
             '        <open>0</open>',/, &
             '        <visibility>1</visibility>',/, &
             '        <description>Folder containing final values</description>',/, &
             '        <Style>',/, &
             '             <ListStyle>',/, &
             '                   <listItemType>checkHideChildren</listItemType>',/, &
             '             </ListStyle>',/, &
             '        </Style>')
2    format('         <Placemark>',/, &
            '         <name>',f9.2,a5,'</name>',/, &
            '         <styleUrl>#',a9,'</styleUrl>',/, &
            '           <MultiGeometry>',/, &
            '             <Point>',/, &
            '                 <altitudeMode>absolute</altitudeMode>',/, &
            '                 <coordinates>',e12.6,',',e12.6,',',e12.6,'</coordinates>',/, &
            '             </Point>',/, &
            '             <Polygon>',/, &
            '               <extrude>0</extrude>',/, &
            '               <tessellate>1</tessellate>',/, &
            '               <altitudeMode>',a13,'</altitudeMode>',/, &
            '               <outerBoundaryIs>',/, &
            '                 <LinearRing>',/, &
            '                   <coordinates>',/, &
                               17x,e12.6,',',e12.6,',',e12.6,/, &
                               17x,e12.6,',',e12.6,',',e12.6,/, &
                               17x,e12.6,',',e12.6,',',e12.6,/, &
                               17x,e12.6,',',e12.6,',',e12.6,/, &
                               17x,e12.6,',',e12.6,',',e12.6,/, &
            '                   </coordinates>',/, &
            '                 </LinearRing>',/, &
            '               </outerBoundaryIs>',/, &
            '             </Polygon>',/, &
            '           </MultiGeometry>',/, &
            '         </Placemark>')
!           For polygons that cross the antimeridian . . .
25   format('         <Placemark>  <!-- crosses the antimeridian -->',/, &
            '         <name>',f9.1,a5,'</name>',/, &
            '         <styleUrl>#',a9,'</styleUrl>',/, &
            '           <MultiGeometry>',/, &
            '             <Point>',/, &
            '                 <coordinates>',e12.6,',',e12.6,',',e12.6,'</coordinates>',/, &
            '             </Point>',/, &
            '             <Polygon>',/, &
            '               <extrude>0</extrude>',/, &
            '               <tessellate>1</tessellate>',/, &
            '               <altitudeMode>',a13,',</altitudeMode>',/, &
            '               <outerBoundaryIs>',/, &
            '                 <LinearRing>',/, &
            '                   <coordinates>',/, &
                               17x,e12.6,',',e12.6,',',e12.6,/, &
                               17x,e12.6,',',e12.6,',',e12.6,/, &
                               17x,e12.6,',',e12.6,',',e12.6,/, &
                               17x,e12.6,',',e12.6,',',e12.6,/, &
                               17x,e12.6,',',e12.6,',',e12.6,/, &
            '                   </coordinates>',/, &
            '                 </LinearRing>',/, &
            '               </outerBoundaryIs>',/, &
            '             </Polygon>',/, &
            '             <Polygon>',/, &
            '               <extrude>0</extrude>',/, &
            '               <tessellate>1</tessellate>',/, &
            '               <altitudeMode>clampToGround</altitudeMode>',/, &
            '               <outerBoundaryIs>',/, &
            '                 <LinearRing>',/, &
            '                   <coordinates>',/, &
                               17x,e12.6,',',e12.6,',',e12.6,/, &
                               17x,e12.6,',',e12.6,',',e12.6,/, &
                               17x,e12.6,',',e12.6,',',e12.6,/, &
                               17x,e12.6,',',e12.6,',',e12.6,/, &
                               17x,e12.6,',',e12.6,',',e12.6,/, &
            '                   </coordinates>',/, &
            '                 </LinearRing>',/, &
            '               </outerBoundaryIs>',/, &
            '             </Polygon>',/, &
            '           </MultiGeometry>',/, &
            '         </Placemark>')
3     format('      </Folder>')

      end subroutine Write_2D_KML

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Write_PointData_Airports_KML
!
!  Called from: output_results (in the isFinal_TS branch) and Ash3d_PostProc.F90
!  Arguments:
!    none
!
!  This subroutine writes all airport data to ash_arrivaltimes_airports.kml.
!  For each airport/POI, a short gnuplot script is written and executed which
!  creates a time-series plot of the deposit at accumulation at that point.
!  These are linked to the kml file.  If zip is installed on the system, the
!  kml file is compressed.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Write_PointData_Airports_KML

      use global_param,  only : &
         IsLinux,IsWindows,IsMacOS,usezip,zippath,&
         usegnuplot,gnuplotpath

      use io_data,       only : &
         nWriteTimes,VolcanoName,WriteTimes

      use mesh,          only : &
         A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2,&
         A3d_k0_scale,A3d_Re,IsLatLon

      use time_data,     only : &
         time,BaseYear,useLeap,SimStartHour,Simtime_in_hours,OutputOffset

      use Output_Vars,   only : &
         CloudLoad,CLOUDLOAD_THRESH,DEPRATE_THRESH,THICKNESS_THRESH

      use Airports,      only : &
         Airport_TS_plotindex,Airport_AshDuration,Airport_AshArrivalTime,&
         Airport_CloudArrivalTime,Airport_CloudDuration,nairports, &
         Airport_Thickness_TS,Airport_Name,Airport_AshArrived, &
         Airport_Longitude,Airport_Latitude,Airport_deprate, &
         Airport_i,Airport_j,Airport_Thickness,Airport_CloudArrived

      use Tephra,          only : &
         n_gs_max

      use Source,        only : &
         neruptions,e_StartTime,e_Duration,e_PlumeHeight,e_Volume,&
         lat_volcano,lon_volcano,x_volcano,y_volcano

      use projection,    only : &
           PJ_proj_inv

      integer             :: i
      integer             :: nWrittenOut
      character (len=13)  :: yyyymmddhh
      character (len=1)   :: cloud_morethan, deposit_morethan      !equals ">" if cloud is still overhead or ash is still falling
      real(kind=dp)       :: CloudTime
      character (len=20)  :: xmlTimeStart, xmlTimeEnd
      real(kind=ip)       :: airlon,airlat

      integer,       dimension(:), allocatable :: iyear, imonth, iday
      real(kind=ip), dimension(:), allocatable :: StartHour
      integer :: ierup
      integer :: ai
      integer :: plt_indx
      character(len=14) :: dp_outfile
      character(len=14) :: dp_gnufile
      character(len=14) :: dp_pngfile
      character(len=80) :: gnucom
      character(len=77) :: zipcom

      real(kind=ip)      :: ymaxpl
      logical            :: IsThere
      integer            :: stat
      real(kind=dp)      :: olam,ophi ! using precision needed by libprojection
      integer            :: iostatus
      character(len=120) :: iomessage
      character(len= 50) :: linebuffer050 
      character(len= 80) :: linebuffer080

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

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Write_PointData_Airports_KML"
      endif;enddo

      ! Loop of all airports in the computational domain and build list of
      ! impacted airports, incrementing a plot index and logging those airport
      ! in Airport_TS_plotindex by noting the polot index.
      plt_indx = 0
      Airport_TS_plotindex(:) = 0
      do ai=1,nairports
        ! First check if ash has accumulated there by more than 0.01 mm
        if(Airport_Thickness_TS(ai,nWriteTimes).lt.THICKNESS_THRESH)then
          cycle
        elseif(Airport_Thickness_TS(ai,nWriteTimes).lt.1.0_ip)then
          ymaxpl = 1.0_ip
        elseif(Airport_Thickness_TS(ai,nWriteTimes).lt.5.0_ip)then
          ymaxpl = 5.0_ip
        elseif(Airport_Thickness_TS(ai,nWriteTimes).lt.25.0_ip)then
          ymaxpl = 25.0_ip
        else
          ymaxpl = 100.0_ip
        endif

        plt_indx = plt_indx +1
        Airport_TS_plotindex(ai) = plt_indx

        ! Writing a gnuplot script for this airport
        write(dp_outfile,53) plt_indx,".dat"
        write(dp_gnufile,53) plt_indx,".gnu"
        write(dp_pngfile,53) plt_indx,".png"
 53     format('depTS_',i4.4,a4)

        open(fid_kmlgnuscr,file=dp_gnufile,status='replace',action='write')
        write(fid_kmlgnuscr,*)"set terminal png size 400,300"
        write(fid_kmlgnuscr,*)"set key bmargin left horizontal Right noreverse enhanced ",&
                              "autotitles box linetype -1 linewidth 1.000"
        write(fid_kmlgnuscr,*)"set border 31 lw 2.0 lc rgb '#000000'"
        write(fid_kmlgnuscr,*)"set style line 1 linecolor rgbcolor '#888888' linewidth 2.0 pt 7"
        write(fid_kmlgnuscr,*)"set ylabel 'Deposit Thickeness (mm)'"
        write(fid_kmlgnuscr,*)"set xlabel 'Time (hours after eruption)'"
        write(fid_kmlgnuscr,*)"set nokey"
        write(fid_kmlgnuscr,*)"set output '",dp_pngfile,"'"
        write(fid_kmlgnuscr,*)"set title '",Airport_Name(ai),"'"
        write(fid_kmlgnuscr,*)"plot [0:",ceiling(Simtime_in_hours),"][0:",&
                              nint(ymaxpl),"] '",dp_outfile,"' with filledcurve x1 ls 1"
        close(fid_kmlgnuscr)
        ! Writing the data file the gnuplot script will plot
        open(fid_kmlgnudat,file=dp_outfile,status='replace',action='write')
        do i = 1,nWriteTimes
           write(fid_kmlgnudat,*)WriteTimes(i),Airport_Thickness_TS(ai,i)
        enddo
        close(fid_kmlgnudat)
        ! Test if gnuplot is installed
        if(usegnuplot)then
          ! if we have gnuplot installed, just create the plots now
          write(gnucom,*)trim(adjustl(gnuplotpath)),' -p ',dp_gnufile
          call execute_command_line(gnucom)
          ! Now delete the script and data files
          open(unit=fid_kmlgnudat, iostat=stat, file=dp_outfile, status='old',action='write')
          if (stat.eq.0) close(fid_kmlgnudat, status='delete')
          open(unit=fid_kmlgnuscr, iostat=stat, file=dp_gnufile, status='old',action='write')
          if (stat.eq.0) close(fid_kmlgnuscr, status='delete')
        endif
      enddo

      ! Now starting the kml file
      open(unit=fid_kmlPOI,file='ash_arrivaltimes_airports.kml',status='replace',action='write',err=2001)
      write(fid_kmlPOI,5)                           ! write out file header info
      nWrittenOut = 0
      ! Loop through all the airports again so we can write out the airports that are hit.
      do ai=1,nairports
        ! First separate the airports into two categories: deposit cases or cloud cases
        if ((n_gs_max.gt.1).and.&         ! Check if this is a web-cloud run by testing the number of grain sizes
            Airport_AshArrived(ai)) then  ! Check if ash has actually arrived here
          ! Deposit case
          deposit_morethan = ' '
          cloud_morethan   = ' '
          airlon = Airport_Longitude(ai)
          airlat = Airport_Latitude(ai)
          xmlTimeStart = HS_xmltime(SimStartHour+Airport_AshArrivalTime(ai)+OutputOffset,&
                                    BaseYear,useLeap)
          !See whether cloud is still overhead, or whether ash is still falling
          if ((Airport_AshArrived(ai)).and.(Airport_deprate(ai).gt.DEPRATE_THRESH)) then
            Airport_AshDuration(ai) = time-Airport_AshArrivalTime(ai)
            deposit_morethan = '>'
          else
            deposit_morethan = ' '
          endif
          if (CloudLoad(Airport_i(ai),Airport_j(ai)).gt.CLOUDLOAD_THRESH) then
            Airport_CloudDuration(ai) = time-Airport_CloudArrivalTime(ai)
            cloud_morethan = '>'
          else
            cloud_morethan = ' '
          endif
          if (Airport_Longitude(ai).gt.180.0_ip) airlon=airlon-360.0_ip
          if (Airport_TS_plotindex(ai).gt.0)then
            write(dp_pngfile,53) Airport_TS_plotindex(ai),".png"

            ! A cumulative deposit plot exists for this point since it has a plot index
            ! Write out a placemark which includes the png of the deposit time-series
            write(fid_kmlPOI,16)Airport_Name(ai),      &
                        Airport_CloudArrivalTime(ai),  &
                        cloud_morethan,                &
                        Airport_CloudDuration(ai),     &
                        Airport_AshArrivalTime(ai),    &
                        deposit_morethan,              &
                        Airport_AshDuration(ai),       &
                        Airport_Thickness(ai),         &
                        Airport_Thickness(ai)/25.4_ip, &
                        dp_pngfile,                    &
                        xmlTimeStart,                  &
                        airlon,                        &
                        airlat
          else
            ! No ash has fallen here, but a cloud have been overhead
            ! Write out a normal placemark without image
            write(fid_kmlPOI,6) Airport_Name(ai),      &
                        Airport_CloudArrivalTime(ai),  &
                        cloud_morethan,                &
                        Airport_CloudDuration(ai),     &
                        Airport_AshArrivalTime(ai),    &
                        deposit_morethan,              &
                        Airport_AshDuration(ai),       &
                        Airport_Thickness(ai),         &
                        Airport_Thickness(ai)/25.4_ip, &
                        xmlTimeStart,                  &
                        airlon,                        &
                        airlat
          endif
          nWrittenOut = nWrittenOut + 1
        elseif (Airport_CloudArrived(ai)) then       ! If a cloud arrived but no deposit
          ! Cloud case
          cloud_morethan   = ' '
          airlon = Airport_Longitude(ai)
          airlat = Airport_Latitude(ai)
          CloudTime = SimStartHour+Airport_CloudArrivalTime(ai)
          xmlTimeStart = HS_xmltime(CloudTime+OutputOffset,&
                                    BaseYear,useLeap)
          ! See whether cloud is still overhead, or whether ash is still falling
          if (CloudLoad(Airport_i(ai),Airport_j(ai)).gt.CLOUDLOAD_THRESH) then
            Airport_CloudDuration(ai) = time-Airport_CloudArrivalTime(ai)
            cloud_morethan = '>'
          else
            cloud_morethan = ' '
          endif
          xmlTimeEnd   = HS_xmltime(CloudTime+Airport_CloudDuration(ai)+OutputOffset,&
                                    BaseYear,useLeap)
          if (Airport_Longitude(ai).gt.180.0_ip) airlon=airlon-360.0_ip
          ! Write out a normal placemark without image
          write(fid_kmlPOI,7) Airport_Name(ai),     &
                      Airport_CloudArrivalTime(ai), &
                      cloud_morethan,               &
                      Airport_CloudDuration(ai),    &
                      xmlTimeStart,                 &
                      xmlTimeEnd,                   &
                      airlon,                       &
                      airlat
          nWrittenOut = nWrittenOut + 1
        endif
      enddo ! ai=1,nairports
      if (IsLatLon.eqv..False.) then      ! Put a placemark at the location of the volcano
        call PJ_proj_inv(real(x_volcano,kind=dp), real(y_volcano,kind=dp),  &
                      A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
                      A3d_k0_scale,A3d_Re, &
                      olam,ophi)
        lon_volcano = real(olam,kind=ip)
        lat_volcano = real(ophi,kind=ip)
      endif

      allocate(iyear(neruptions))         ! Generate data needed to write ESP's
      allocate(imonth(neruptions))
      allocate(iday(neruptions))
      allocate(StartHour(neruptions))
      write(fid_kmlPOI,8) VolcanoName, VolcanoName
      !  Write table of ESP's
      do ierup=1,neruptions
        yyyymmddhh = HS_yyyymmddhh_since(e_StartTime(ierup)+SimStartHour+OutputOffset,&
                                         BaseYear,useLeap)
        read(yyyymmddhh,200,iostat=iostatus,iomsg=iomessage) &
                      iyear(ierup),imonth(ierup),iday(ierup), &
                      StartHour(ierup)
        linebuffer080 = yyyymmddhh
        linebuffer050 = "Reading date from date_string (write_KML)"
        if(iostatus.ne.0) call FileIO_Error_Handler(iostatus,linebuffer050,linebuffer080,iomessage)
        write(fid_kmlPOI,9) ierup,iyear(ierup),imonth(ierup),iday(ierup), &
                      StartHour(ierup), &
                      e_PlumeHeight(ierup),e_Duration(ierup),e_Volume(ierup)
200     format(i4,i2,i2,f5.2)
      enddo
      if (nWrittenOut.gt.0) then             ! If one or more airports were affected
        if (lon_volcano.gt.180.0_ip) then
          write(fid_kmlPOI,10) lon_volcano-360.0_ip, lat_volcano
        else
          write(fid_kmlPOI,10) lon_volcano, lat_volcano
        endif
      else                                   ! If no airports were affected
        if (lon_volcano.gt.180.0_ip) then
          write(fid_kmlPOI,11) lon_volcano-360.0_ip, lat_volcano
        else
          write(fid_kmlPOI,11) lon_volcano, lat_volcano
        endif
      endif
      write(fid_kmlPOI,12)                           ! write final lines of file
      close(fid_kmlPOI)                              ! close file

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),4) nWrittenOut      ! Write number of airports affected to log file & stdout
      endif;enddo

      ! Test if zip is installed
      IsThere = .false.
      if(usezip.or.IsLinux.or.IsMacOS)then
        ! double-check that zip is installed
        inquire( file=trim(adjustl(zippath)), exist=IsThere)
      elseif(IsWindows)then
        ! this is a placeholder for now
        IsThere = .false.
      else
        IsThere = .false.
      endif
      if(usezip)then
        write(zipcom,'(a77)')&
          'zip -r ash_arrivaltimes_airports.kmz ash_arrivaltimes_airports.kml depTS*.png'
        call execute_command_line(zipcom)
      endif

      return

      ! Error traps
2001  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)'Error opening ash_arrivaltimes_airports.kml.  Program stopped.'
      endif;enddo
      stop 1


4     format(/,4x,'Number of airports impacted by ash = ',i4)
5     format('<?xml version="1.0" encoding="UTF-8"?>',/, &
             '<kml xmlns="http://www.opengis.net/kml/2.2">',/, &
             '<Document>',/, &
             '    <StyleMap id="VolcanoMarker">',/, &
             '      <Pair>',/, &
             '             <key>normal</key>',/, &
             '             <styleUrl>#VolcanoMarkerNormal</styleUrl>',/, &
             '      </Pair>',/, &
             '      <Pair>',/, &
             '             <key>highlight</key>',/, &
             '             <styleUrl>#VolcanoMarkerHighlight</styleUrl>',/, &
             '      </Pair>',/, &
             '    </StyleMap>',/, &
             '    <Style id="VolcanoMarkerNormal">',/, &
             '      <IconStyle>',/, &
             '             <scale>0.8</scale>',/, &
             '             <Icon>',/, &
             '               <href>http://maps.google.com/mapfiles/kml/pal4/icon52.png</href>',/, &
             '             </Icon>',/, &
             '      </IconStyle>',/, &
             '      <LabelStyle>',/, &
             '             <scale>0</scale>',/, &
             '      </LabelStyle>',/, &
             '    </Style>',/, &
             '    <Style id="VolcanoMarkerHighlight">',/, &
             '      <IconStyle>',/, &
             '             <scale>1</scale>',/, &
             '             <Icon>',/, &
             '               <href>http://maps.google.com/mapfiles/kml/pal4/icon52.png</href>',/, &
             '             </Icon>',/, &
             '      </IconStyle>',/, &
             '      <LabelStyle>',/, &
             '             <scale>1</scale>',/, &
             '      </LabelStyle>',/, &
             '    </Style>',/, &
             '  <ScreenOverlay>',/, &
             '    <name>Warning</name>',/, &
             '    <description>USGS disclaimer</description>',/, &
             '      <Icon>',/, &
             '        <href>https://vsc-ash.wr.usgs.gov/images/USGS_warning3.png</href>',/, &
             '      </Icon>',/, &
             '    <overlayXY x="0.5" y="1.0" xunits="fraction" yunits="fraction"/>',/, &
             '    <screenXY x="0.5" y="1.0" xunits="fraction" yunits="fraction"/>',/, &
             '    <size x="337" y="68" xunits="pixels" yunits="pixels"/>',/, &
             '  </ScreenOverlay>',/, &
             '    <StyleMap id="AirportCldMarker">',/, &
             '      <Pair>',/, &
             '             <key>normal</key>',/, &
             '             <styleUrl>#AirportCldMarkerNormal</styleUrl>',/, &
             '      </Pair>',/, &
             '      <Pair>',/, &
             '             <key>highlight</key>',/, &
             '             <styleUrl>#AirportCldMarkerHighlight</styleUrl>',/, &
             '      </Pair>',/, &
             '    </StyleMap>',/, &
             '    <Style id="AirportCldMarkerNormal">',/, &
             '      <IconStyle>',/, &
             '             <scale>0.8</scale>',/, &
             '             <Icon>',/, &
             '               <href>http://maps.google.com/mapfiles/kml/pal4/icon48.png</href>',/, &
             '             </Icon>',/, &
             '      </IconStyle>',/, &
             '      <LabelStyle>',/, &
             '             <scale>0</scale>',/, &
             '      </LabelStyle>',/, &
             '    </Style>',/, &
             '    <Style id="AirportCldMarkerHighlight">',/, &
             '      <IconStyle>',/, &
             '             <scale>1</scale>',/, &
             '             <Icon>',/, &
             '               <href>http://maps.google.com/mapfiles/kml/pal4/icon48.png</href>',/, &
             '             </Icon>',/, &
             '      </IconStyle>',/, &
             '      <LabelStyle>',/, &
             '             <scale>1</scale>',/, &
             '      </LabelStyle>',/, &
             '    </Style>',/, &
             '    <StyleMap id="AirportAshMarker">',/, &
             '      <Pair>',/, &
             '             <key>normal</key>',/, &
             '             <styleUrl>#AirportAshMarkerNormal</styleUrl>',/, &
             '      </Pair>',/, &
             '      <Pair>',/, &
             '             <key>highlight</key>',/, &
             '             <styleUrl>#AirportAshMarkerHighlight</styleUrl>',/, &
             '      </Pair>',/, &
             '    </StyleMap>',/, &
             '    <Style id="AirportAshMarkerNormal">',/, &
             '      <IconStyle>',/, &
             '             <scale>0.8</scale>',/, &
             '             <Icon>',/, &
             '               <href>http://maps.google.com/mapfiles/kml/pal4/icon56.png</href>',/, &
             '             </Icon>',/, &
             '      </IconStyle>',/, &
             '      <LabelStyle>',/, &
             '             <scale>0</scale>',/, &
             '      </LabelStyle>',/, &
             '    </Style>',/, &
             '    <Style id="AirportAshMarkerHighlight">',/, &
             '      <IconStyle>',/, &
             '             <scale>1</scale>',/, &
             '             <Icon>',/, &
             '               <href>http://maps.google.com/mapfiles/kml/pal4/icon56.png</href>',/, &
             '             </Icon>',/, &
             '      </IconStyle>',/, &
             '      <LabelStyle>',/, &
             '             <scale>1</scale>',/, &
             '      </LabelStyle>',/, &
             '    </Style>')
6     format('  <Placemark>',/, &
             '    <name>',a35,'</name>',/, &
             '    <styleUrl>#AirportAshMarker</styleUrl>',/, &
             '    <description><![CDATA[',/, &
             '                <ul>',/, &
             '                  <li>Cloud arrival time: ',f10.2,' hrs after eruption start</li>',/, &
             '                  <li>Cloud will remain overhead for ',a1,f10.2,' hours**</li>',/, &
             '                  <li>Deposit start time: ',f10.2,' hrs* after eruption start</li>',/, &
             '                  <li>Deposit will fall for ',a1,f10.2,' hours**</li>',/, &
             '                  <li>Total ash thickness: ',f8.2,' mm (',f7.3,' in.)</li>',/, &
             '                </ul>',/, &
             '                **a duration preceded by a ">" sign indicates the cloud is still overhead ', &
                                'or the deposit is still falling at the end of the simulation.', &
             '                ]]>',/, &
             '    </description>',/, &
             '    <TimeSpan>',/, &
             '        <begin>',a20,'</begin>',/, &
             '    </TimeSpan>',/, &
             '    <Point>',/, &
             '      <coordinates>',e12.6,',',e12.6,',0</coordinates>',/, &
             '    </Point>',/, &
             '  </Placemark>')
16    format('  <Placemark>',/, &
             '    <name>',a35,'</name>',/, &
             '    <styleUrl>#AirportAshMarker</styleUrl>',/, &
             '    <description><![CDATA[',/, &
             '                <ul>',/, &
             '                  <li>Cloud arrival time: ',f10.2,' hrs after eruption start</li>',/, &
             '                  <li>Cloud will remain overhead for ',a1,f10.2,' hours**</li>',/, &
             '                  <li>Deposit start time: ',f10.2,' hrs* after eruption start</li>',/, &
             '                  <li>Deposit will fall for ',a1,f10.2,' hours**</li>',/, &
             '                  <li>Total ash thickness: ',f8.2,' mm (',f7.3,' in.)</li>',/, &
             '                </ul>',/, &
             '                **a duration preceded by a ">" sign indicates the cloud is still overhead ', &
                                'or the deposit is still falling at the end of the simulation.', &
             '                ]]>',/, &
             '                <img src="',a14,'" />', &
             '    </description>',/, &
             '    <TimeSpan>',/, &
             '        <begin>',a20,'</begin>',/, &
             '    </TimeSpan>',/, &
             '    <Point>',/, &
             '      <coordinates>',e12.6,',',e12.6,',0</coordinates>',/, &
             '    </Point>',/, &
             '  </Placemark>')
7     format('  <Placemark>',/, &
             '    <name>',a35,'</name>',/, &
             '    <styleUrl>#AirportCldMarker</styleUrl>',/, &
             '    <description><![CDATA[',/, &
             '                <ul>',/, &
             '                  <li>Cloud arrival time: ',f10.2,' hrs after eruption start</li>',/, &
             '                  <li>Cloud will remain overhead for ',a1,f10.2,' hours*</li>',/, &
             '                </ul>',/, &
             '                 *a duration preceded by a ">" sign indicates the cloud is still overhead ', &
                                'at the end of the simulation.', &
             '                ]]>',/, &
             '    </description>',/, &
             '    <TimeSpan>',/, &
             '        <begin>',a20,'</begin>',/, &
             '        <end>',a20,'</end>',/, &
             '    </TimeSpan>',/, &
             '    <Point>',/, &
             '      <coordinates>',e12.6,',',e12.6,',0</coordinates>',/, &
             '    </Point>',/, &
             '  </Placemark>')
8     format('  <Placemark>',/, &
             '    <name>',a35,'</name>',/, &
             '    <styleUrl>#VolcanoMarker</styleUrl>',/, &
             '    <gx:balloonVisibility>1</gx:balloonVisibility>',/, &
             '    <description>',/, &
             '      <![CDATA[',/, &
             '        <font size=4>',/, &
             '          <p><b>Volcano name: </b>',a35,'</p>',/, &
             '        </font>',/, &
             '        <hr />',/, &
             '        <font size=4><b>Eruption source parameters:</b></font>',/, &
             '        <table border="1" padding="0" width="400">',/, &
             '          <tr><td>Pulse #</td><td>Start year</td><td>month</td><td>day</td>',/, &
             '            <td>hour (UTC)</td><td>plume ht (km asl)</td>',/, &
             '            <td>Duration (hrs)</td><td>Volume (km3)</td></tr>')
9     format('          <tr><td>',i2,'</td><td>',i4,'</td><td>',i2,/, &
             ' </td><td>',i2,'</td><td>',f5.2,'</td><td>',f4.1,'</td>',/, &
             '            <td>',f6.2,'</td><td>',e8.3,'</td></tr>')
10    format('         </table>',/, &
             '      ]]>',/, &
             '    </description>',/, &
             '    <Point>',/, &
             '      <coordinates>',e12.6,',',e12.6,',0</coordinates>',/, &
             '    </Point>',/, &
             '  </Placemark>')
11    format('         </table>',/, &
             '      ]]>',/, &
             '      <p><b>No airports affected in this run</b></p>',/, &
             '    </description>',/, &
             '    <Point>',/, &
             '      <coordinates>',e12.6,',',e12.6,',0</coordinates>',/, &
             '    </Point>',/, &
             '  </Placemark>')
12    format('</Document>',/, &
             '</kml>')

      end subroutine Write_PointData_Airports_KML

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Close_KML
!
!  Called from: 
!  Arguments:
!    ivar        = ID of variable to process
!    TS_flag     = 0 = not a time-series, 1 = time-series
!
!  This subroutine closes the kml file associated with the variable ivar
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Close_KML(ivar,TS_flag)

      integer,intent(in) :: ivar
      integer,intent(in) :: TS_flag

      integer :: fid

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine Close_KML"
      endif;enddo

      fid = KML_fid(ivar)

      do io=1,2;if(VB(io).le.verbosity_info)then
        write(outlog(io),15)KML_filename(ivar)
      endif;enddo

      if(TS_flag.ne.0)then
        write(fid,11)
      else
        write(fid,12)
      endif
      close(fid)      !close kml file

15    format(/,5x,'Closing kml file ',a30)
11    format('   </Folder>',/,'</Document>',/,'</kml>')
12    format(/,'</Document>',/,'</kml>')

      end subroutine Close_KML

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  PlotModelBoundary
!
!  Called from: OpenFile_KML
!  Arguments:
!    xleft,xright,ybottom,ytop,fid
!
!  This subroutine finds the latitude & longitude of points in all sides of
!  the model region
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine PlotModelBoundary(xleft,xright,ybottom,ytop,fid)

      use mesh,          only : &
         A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2,&
         A3d_k0_scale,A3d_Re,IsLatLon

      use projection,    only : &
           PJ_proj_inv

      real(kind=ip) :: xplot(0:40),yplot(0:40),lonplot(0:40),latplot(0:40)
      real(kind=ip) :: xleft,xright,ybottom,ytop
      integer       :: ict, fid
      real(kind=dp)  :: olam,ophi ! using precision needed by libprojection

      do io=1,2;if(VB(io).le.verbosity_debug1)then
        write(outlog(io),*)"     Entered Subroutine PlotModelBoundary"
      endif;enddo

      write(fid,3) ! write style for model boundary

      ! For projected coordinates, plot ten points on each side of the polygon.
      do ict=0,10
        xplot(ict) = xleft+ict*(xright-xleft)/10.0_ip
        yplot(ict) = ybottom
      enddo
      do ict=11,20
        xplot(ict) = xright
        yplot(ict) = ybottom + (ict-10)*(ytop-ybottom)/10.0_ip
      enddo
      do ict=21,30
        xplot(ict) = xright - (ict-20)*(xright-xleft)/10.0_ip
        yplot(ict) = ytop
      enddo
      do ict=31,40
        xplot(ict) = xleft
        yplot(ict) = ytop - (ict-30)*(ytop-ybottom)/10.0_ip
      enddo

      ! Convert points to latitude & longitude
      if (IsLatLon) then
        lonplot = xplot
        latplot = yplot
       else
        do ict=0,40
          call PJ_proj_inv(real(xplot(ict),kind=dp),real(yplot(ict),kind=dp),  &
                        A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
                        A3d_k0_scale,A3d_Re, &
                        olam,ophi)
          lonplot(ict) = real(olam,kind=ip)
          latplot(ict) = real(ophi,kind=ip)
        enddo
      endif

      ! Make sure the longitude is between -180 and +180 degrees
      do ict=0,40
        if (lonplot(ict).gt.180.0_ip)  lonplot(ict) = lonplot(ict)-360.0_ip
        if (lonplot(ict).lt.-180.0_ip) lonplot(ict) = lonplot(ict)+360.0_ip
      enddo

      ! Write out the polygon
      write(fid,5) (lonplot(ict),latplot(ict), ict=0,40)

      ! Format statements
3     format('	  <Style id="boundary_style">',/, &             ! style for model boundary
             '		<PolyStyle>',/, &
             '			<color>33ffffff</color>',/, &
             '                  <fill>0</fill>',/, &
             '		</PolyStyle>',/, &
             '      <LineStyle>',/, &
             '          <width>1</width>',/, &
             '      </LineStyle>',/, &
             '	  </Style>')
5    format('	<Placemark>',/, &
            '     <name>Model boundary</name>',/,&
            '		<styleUrl>#boundary_style</styleUrl>',/, &
            '		<Polygon>',/, &
            '                   <altitudeMode>clampToGround</altitudeMode>',/, &
            '			<outerBoundaryIs>',/, &
            '				<LinearRing>',/, &
            '					<coordinates>',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
                                    20x,e12.6,',',e12.6,',0',/, &
            '                   </coordinates>',/, &
            '               </LinearRing>',/, &
            '           </outerBoundaryIs>',/, &
            '		</Polygon>',/, &
            '	</Placemark>')

      return

      end subroutine PlotModelBoundary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module Ash3d_KML_IO
!##############################################################################
