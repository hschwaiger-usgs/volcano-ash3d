      module Output_KML

      use precis_param

      use io_units

      integer, parameter :: nvars      = 10  ! Number of output variables with style profiles
      integer, parameter :: max_nclrmp = 11   ! Max number of colormap points

      !ivar = 1 :: cloud concentration
      !ivar = 2 :: cloud height (top)
      !ivar = 3 :: cloud height (bot)
      !ivar = 4 :: cloud load
      !ivar = 5 :: cloud arrival time
      !ivar = 6 :: cloud reflectivity
      !ivar = 7 :: deposit
      !ivar = 8 :: deposit (NWS)
      !ivar = 9 :: deposit time
      !ivar =10 :: topography

      character(len=30),dimension(nvars           ) :: KML_filename
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

!******************************************************************************

      subroutine Set_OutVar_Specs

      implicit none

      integer :: ivar

      ivar = 1 ! cloud concentration
      KML_filename(ivar)      = 'CloudConcentration.kml       '
      KML_units(ivar)         = 'mg/m3'
      KML_fid(ivar)           = 40
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
      KML_filename(ivar)      = 'CloudHeight.kml              '
      KML_units(ivar)         = ' km  '
      KML_fid(ivar)           = 170
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
      KML_filename(ivar)      = 'CloudBottom.kml              '
      KML_units(ivar)         = ' km  '
      KML_fid(ivar)           = 171
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
      KML_filename(ivar)      = 'CloudLoad.kml                '
      KML_units(ivar)         = 'T/km2'
      KML_fid(ivar)           = 160
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
      !KML_filename(ivar)      = 'CloudArrivalTime.kml         '
      KML_filename(ivar)      = 'cloud_arrivaltimes_hours.kml'
      KML_units(ivar)         = ' hrs '
      KML_fid(ivar)           = 390
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
      KML_filename(ivar)      = 'reflectivity.kml             '
      KML_units(ivar)         = ' dBZ '
      KML_fid(ivar)           = 206
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
      !KML_filename(ivar)      = 'Deposit.kml                  '
      KML_filename(ivar)      = 'deposit_thickness_mm.kml     '
      KML_units(ivar)         = '  mm '
      KML_fid(ivar)           = 50
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
      !KML_filename(ivar)      = 'Deposit_NWS.kml              '
      KML_filename(ivar)      = 'deposit_thickness_inches.kml  '
      KML_units(ivar)         = '  in.'
      KML_fid(ivar)           = 550
      KML_n_clrmp(ivar)       = 3
      KML_color_map(ivar,:) = (/ 0.003937_ip, 0.0315_ip, 0.2362_ip, 0.0_ip, 0.0_ip,&
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
      KML_sizeY(ivar)         = '319'
      KML_AltMode(ivar)       = 'clampToGround'

      ivar = 9 ! deposit time
      !KML_filename(ivar)      = 'DepositArrivalTime.kml       '
      KML_filename(ivar)      = 'ashfall_arrivaltimes_hours.kml'
      KML_units(ivar)         = ' hrs '
      KML_fid(ivar)           = 290
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
      KML_filename(ivar)      = 'Topography.kml                '
      KML_units(ivar)         = '  km '
      KML_fid(ivar)           = 540
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

!******************************************************************************

      subroutine OpenFile_KML(ivar)
      
      !Subroutine that opens and initializes the KML file

      use time_data,     only : &
         BaseYear,useLeap,SimStartHour,OutputOffset

      use mesh,          only : &
         A3d_iprojflag,A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
         A3d_k0_scale,A3d_radius_earth,de,dn,dx,dy,IsLatLon, &
         latLL,lonLL,latUR,lonUR,xLL,yLL,xUR,yUR

      use io_data,       only : &
         VolcanoName,WriteDepositTS_KML

      use Source,        only : &
         neruptions, e_StartTime,e_Duration,e_Volume,e_PlumeHeight, &
         lon_volcano,lat_volcano,x_volcano,y_volcano

      use projection,    only : &
           PJ_proj_inv

      implicit none

      integer             :: ivar
      !character (len=13)  :: HS_yyyymmddhh_since
      character (len=13)  :: yyyymmddhh
      character (len=2)   :: opacity
      real(kind=ip)       :: xleft, xright, ybottom, ytop
      real(kind=ip)       :: longLL,longUR
      real(kind=ip)       :: lattLL,lattUR

      integer,       dimension(:), allocatable :: iyear, imonth, iday
      real(kind=ip), dimension(:), allocatable :: StartHour

      integer             :: ierup
      character(len=30)   :: filename
      integer             :: fid
      integer             :: n_clrmp,icmp
      !real(kind=ip)   ,dimension(9)   :: color_map
      character(len=9),dimension(11)   :: Styles
      character(len=6),dimension(11)   :: Colors
      character(len=30)   :: description
      character(len=30)   :: legend
      character(len=3)    :: overlayX
      character(len=3)    :: overlayY
      character(len=3)    :: screenX
      character(len=3)    :: screenY
      character(len=3)    :: sizeX
      character(len=3)    :: sizeY

      INTERFACE
        character (len=13) function HS_yyyymmddhh_since(HoursSince,byear,useLeaps)
          real(kind=8)               ::  HoursSince
          integer                    ::  byear
          logical                    ::  useLeaps
        end function HS_yyyymmddhh_since
      END INTERFACE

      filename    = KML_filename(ivar)
      fid         = KML_fid(ivar)
      n_clrmp     = KML_n_clrmp(ivar)
      !color_map   = KML_color_map(ivar,:)
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
      
      open(fid,file=trim(adjustl(filename)),status='unknown',err=2500)

      write(fid,1)                 !write file header

      !write StyleMap entries
      do icmp = 1,n_clrmp
        write(fid,2) Styles(icmp), Styles(icmp), Styles(icmp)
      enddo

      !write highlighted styles
      !write(fid,3) 'PureWhite',   opacity, 'ffffff'
      !write(fid,3) 'pure__red',   opacity, '0000ff'
      do icmp = 1,n_clrmp
        write(fid,3) Styles(icmp), opacity, Colors(icmp)
      enddo

      !write normal styles
      do icmp = 1,n_clrmp
        write(fid,4) Styles(icmp), opacity, Colors(icmp)
      enddo

      !write legend
      write(fid,5)description,legend,overlayX,overlayY,screenX,screenY,sizeX,sizeY

      !PLOT MODEL REGION
      if (IsLatLon) then
        longLL  = lonLL - de/2.0_ip
        longUR  = lonUR + de/2.0_ip
        lattLL  = latLL - dn/2.0_ip
        lattUR  = latUR + dn/2.0_ip
        call PlotModelBoundary(longLL,longUR,lattLL,lattUR,fid)
      else
        xleft   = xLL - dx/2.0_ip
        xright  = xUR + dx/2.0_ip
        ybottom = yLL - dy/2.0_ip
        ytop    = yUR + dy/2.0_ip
        call PlotModelBoundary(xleft,xright,ybottom,ytop,fid)
      endif ! IsLatLon

      allocate(iyear(neruptions))
      allocate(imonth(neruptions))
      allocate(iday(neruptions))
      allocate(StartHour(neruptions))

      write(fid,7) VolcanoName, VolcanoName

!     WRITE TABLE OF ESP'S TO THE VOLCANO PLACEMARK
      do ierup=1,neruptions
        yyyymmddhh = HS_yyyymmddhh_since(e_StartTime(ierup)+SimStartHour+OutputOffset,&
                                         BaseYear,useLeap)
        read(yyyymmddhh,100) iyear(ierup),imonth(ierup),iday(ierup), &
                             StartHour(ierup)
        write(fid,8) ierup,iyear(ierup),imonth(ierup),iday(ierup), &
                             StartHour(ierup), &
                             e_PlumeHeight(ierup),e_Duration(ierup),e_Volume(ierup)
100     format(i4,i2,i2,f5.2)

      enddo
      !PLOT VOLCANO
      if (.not.IsLatLon) then                        !get lon_volcano and lat_volcano
        call PJ_proj_inv(x_volcano, y_volcano, &
                   A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
                   A3d_k0_scale,A3d_radius_earth, &
                   lon_volcano,lat_volcano)
      endif
      write(fid,9) lon_volcano, lat_volcano

      !CREATE FORECAST FOLDER ONLY FOR THE FILES THAT HAVE TIME STEPS
      if ((ivar.eq.5).or. &                                       !cloud arrival time
         ((ivar.eq.7).and.(.not.WriteDepositTS_KML)).or. & !deposit
         ((ivar.eq.8).and.(.not.WriteDepositTS_KML)).or. & !deposit_NWS
          (ivar.eq.9)) then                                       !deposit arrival time
        continue
      else
        write(fid,6)                                       !create folder of forecasts
      endif

      deallocate(iyear,imonth,iday,StartHour)

      return

      !Error traps
2500  write(6,20)
      write(9,20)
      stop 1

      !Format statements
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
!            '               <listItemType>radioFolder</listItemType>',/, &
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

!******************************************************************************
!******************************************************************************

      subroutine Write_2D_KML(ivar,OutVar,height_flag,TS_flag)
      
      use mesh,          only : &
         nxmax,nymax,A3d_iprojflag,A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2,&
         A3d_k0_scale,A3d_radius_earth,de,dn,dx,dy,IsLatLon,&
         lon_cc_pd,lat_cc_pd,x_cc_pd,y_cc_pd

      use time_data,     only : &
         time,BaseYear,useLeap,SimStartHour,OutputOffset, &
         xmlTimeSpanStart,xmlTimeSpanEnd

      use Output_Vars,   only : &
         MinHeight,MaxHeight

      use projection,    only : &
           PJ_proj_inv

      implicit none       

      integer        :: ivar
      real(kind=ip)  :: OutVar(nxmax,nymax)
      integer        :: height_flag          ! <0 : min, =0 : ground, >0 : max
      integer        :: TS_flag              ! 0 = not a time-series, 1 = time-series

      integer             :: i,j
      character (len=9)   :: StyleNow3
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
      !character(len=6),dimension(9)   :: Colors
      logical             :: CrossAntiMeridian     !if the polygon crosses the antimeridian
      !character(len=1)    :: answer                !for debugging

      INTERFACE
        character (len=20) function HS_xmltime(HoursSince,byear,useLeaps)
          real(kind=8)              :: HoursSince
          integer                   :: byear
          logical                   :: useLeaps
        end function HS_xmltime
      END INTERFACE


      fid       = KML_fid(ivar)
      n_clrmp   = KML_n_clrmp(ivar)
      color_map = KML_color_map(ivar,:)
      Styles(:) = KML_Styles(ivar,:)
      !Colors(:) = KML_Colors(ivar,:)
      units     = KML_units(ivar)
      AltMode   = KML_AltMode(ivar)

      xmlArrivalTime = HS_xmltime(SimStartHour + time + OutputOffset,&
                                  BaseYear,useLeap)

      StyleNow3 = 'PureWhite'

      if(TS_flag.ne.0)then
        !<<<< for debugging
        !if (ivar.eq.4) then
        !   write(6,*) 'in write_KML'
        !   write(6,*) 'time=',time
        !   write(6,*) 'eruption start time=',HS_xmltime(SimStartHour,BaseYear,useLeap)
        !   write(6,1) xmlArrivalTime, xmlArrivalTime,  &
        !         xmlTimeSpanStart, xmlTimeSpanEnd
        !   write(6,*) 'Continue?'
        !   read(5,'(a1)') answer
        !   if (answer.eq.'n') stop
        !end if
        write(fid,1) xmlArrivalTime, xmlArrivalTime,  &
                 xmlTimeSpanStart, xmlTimeSpanEnd
        !close forecast folder for deposit files on the last time step
        !elseif  &
        !  (((ivar.eq.7).and.WriteDepositTS_KML).or. & !deposit
        !   ((ivar.eq.8).and.WriteDepositTS_KML)) then  !deposit_NWS
        !    write(fid,3)
        !    write(fid,15) 
      else
        write(fid,15) 
      endif
  
      !close folder if this is the final deposit in a deposit file
      !if (final.and.WriteDepositFinal_KML)  write(fid,3)

      do i=1,nxmax
        do j=1,nymax
          !if (OutVar(i,j).lt.color_map(1)) goto 100
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
            call PJ_proj_inv(xleft, ybottom,  &
                        A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
                        A3d_k0_scale,A3d_radius_earth, &
                           longLL,lattLL)
            call PJ_proj_inv(xleft,    ytop,  &
                        A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
                        A3d_k0_scale,A3d_radius_earth, &
                           longUL,lattUL)
            call PJ_proj_inv(xright,   ytop,  &
                        A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
                        A3d_k0_scale,A3d_radius_earth, &
                           longUR,lattUR)
            call PJ_proj_inv(xright,ybottom,  &
                        A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
                        A3d_k0_scale,A3d_radius_earth, &
                           longLR,lattLR)
            call PJ_proj_inv(x_cc_pd(i),y_cc_pd(j),  &
                        A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
                        A3d_k0_scale,A3d_radius_earth, &
                           longCC,lattCC)
          endif ! IsLatLon
          if(height_flag>0)then
            height = int(MaxHeight(i,j)*1000.0_ip)
          elseif(height_flag.lt.0)then
            height = int(MinHeight(i,j)*1000.0_ip)
          else
            height = 0
          endif
          if (LongUR.lt.LongUL) then
            CrossAntiMeridian = .true.       !polygon crosses the antimeridian
            !establish two polygons.  The first is left of the AM
            longLL1=longLL
            longUL1=longUL
            longLR1=179.9999
            longUR1=179.9999
            lattLL1=lattLL
            lattUL1=lattUL
            !The second polygon is right of the AM
            longLL2=-179.9999
            longUL2=-179.9999
            longLR2=longLR
            longUR2=longUR
            lattLR2=lattLR
            lattUR2=lattUR
            !interpolate to find lattLR1,lattUR1,lattLL2,lattUL2
            longUR3=longUR+360.0_ip
            longLR3=longLR+360.0_ip
            lattUL1=lattUL
            !The second polygon is right of the AM
            longLL2=-179.9999
            longUL2=-179.9999
            longLR2=longLR
            longUR2=longUR
            lattLR2=lattLR
            lattUR2=lattUR
            !interpolate to find lattLR1,lattUR1,lattLL2,lattUL2
            longUR3=longUR+360.0_ip
            longLR3=longLR+360.0_ip
            lattLR1=lattLL+(lattLR-lattLL)*(179.9999-longLL)/(longLR3-longLL)
            lattUR1=lattUL+(lattUR-lattUL)*(179.9999-longUL)/(longUR3-longUL)
            lattLL2=lattLR1
            lattUL2=lattUR1
            !write(global_info,*) 'CrossAntiMeridian=',CrossAntiMeridian
            !write(global_info,*) 'longLL=',longLL,', longUR=',longUR
            !write(global_info,*) 'Continue?'
            !read(5,'(a1)') answer
            !if (answer.eq.'n') stop 1
          else
            CrossAntiMeridian = .false.
          endif
          !write out polygon to  kml file
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
             !This fixes the antimeridian problem
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
!100       continue
        enddo
      enddo

      write(fid,3)   !close folder

      return
      
      !format statements
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

!******************************************************************************

      subroutine Write_PointData_Airports_KML

      use Output_Vars,   only : &
         CloudLoad,CLOUDLOAD_THRESH,DEPRATE_THRESH

      use Airports,      only : &
         Airport_TS_plotindex,Airport_AshDuration,Airport_AshArrivalTime,&
         Airport_CloudArrivalTime,Airport_CloudDuration,nairports, &
         Airport_Thickness_TS,Airport_Name,Airport_AshArrived, &
         Airport_Longitude,Airport_Latitude,Airport_deprate, &
         Airport_i,Airport_j,Airport_Thickness,Airport_CloudArrived

      use time_data,     only : &
         time,BaseYear,useLeap,SimStartHour,Simtime_in_hours,OutputOffset

      use mesh,          only : &
         A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2,&
         A3d_k0_scale,A3d_radius_earth,IsLatLon

      use io_data,       only : &
         nWriteTimes,VolcanoName,WriteTimes

      use mesh,          only : &
         nsmax

      use Source,        only : &
         neruptions,e_StartTime,e_Duration,e_PlumeHeight,e_Volume,&
         lat_volcano,lon_volcano,x_volcano,y_volcano

      use projection,    only : &
           PJ_proj_inv

      implicit none

      integer             :: i
      integer             :: nWrittenOut
      character (len=13)  :: yyyymmddhh
      character (len=1)   :: cloud_morethan, deposit_morethan      !equals ">" if cloud is still overhead or ash is still falling
      real(kind=ip)       :: CloudTime
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
      real(kind=ip) :: ymaxpl

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

      plt_indx = 0
      Airport_TS_plotindex = 0
      do ai=1,nairports
        if(Airport_Thickness_TS(ai,nWriteTimes).lt.0.01_ip)then
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
        write(dp_outfile,53) plt_indx,".dat"
        write(dp_gnufile,53) plt_indx,".gnu"
        write(dp_pngfile,53) plt_indx,".png"
 53     format('depTS_',i4.4,a4)

        open(55,file=dp_gnufile,status='replace')
        write(55,*)"set terminal png size 400,300"
        write(55,*)"set key bmargin left horizontal Right noreverse enhanced ",&
                   "autotitles box linetype -1 linewidth 1.000"
        write(55,*)"set border 31 lw 2.0 lc rgb '#000000'"
        write(55,*)"set style line 1 linecolor rgbcolor '#888888' linewidth 2.0 pt 7"
        write(55,*)"set ylabel 'Deposit Thickeness (mm)'"
        write(55,*)"set xlabel 'Time (hours after eruption)'"
        write(55,*)"set nokey"
        write(55,*)"set output '",dp_pngfile,"'"
        write(55,*)"set title '",Airport_Name(ai),"'"
        write(55,*)"plot [0:",ceiling(Simtime_in_hours),"][0:",&
                   nint(ymaxpl),"] '",dp_outfile,"' with filledcurve x1 ls 1"
        close(55)

        open(54,file=dp_outfile,status='replace')
        do i = 1,nWriteTimes
           write(54,*)WriteTimes(i),Airport_Thickness_TS(ai,i)
        enddo
        close(54)
      enddo

      open(unit=60,file='ash_arrivaltimes_airports.kml',err=2001)
      write(60,5)                           !write out file header info
      nWrittenOut = 0
      do i=1,nairports                      !write out the airports that are hit.

        !if (Airport_AshArrived(i)) then      !if this is not a cloud run, and
        !a deposit has arrived at this location
        if ((nsmax.gt.1).and.Airport_AshArrived(i)) then      !if this is not a cloud run, and a deposit has arrived at this location
          deposit_morethan = ' '
          cloud_morethan   = ' '
          airlon = Airport_Longitude(i)
          airlat = Airport_Latitude(i)
          xmlTimeStart = HS_xmltime(SimStartHour+Airport_AshArrivalTime(i)+OutputOffset,&
                                    BaseYear,useLeap)
          !See whether cloud is still overhead, or whether ash is still
          !falling
          if ((Airport_AshArrived(i)).and.(Airport_deprate(i).gt.DEPRATE_THRESH)) then
            Airport_AshDuration(i) = time-Airport_AshArrivalTime(i)
            deposit_morethan = '>'
          else
            deposit_morethan = ' '
          endif
          if (CloudLoad(Airport_i(i),Airport_j(i)).gt.CLOUDLOAD_THRESH) then
            Airport_CloudDuration(i) = time-Airport_CloudArrivalTime(i)
            cloud_morethan = '>'
          else
            cloud_morethan = ' '
          endif
          if(ai.gt.0)then
            ! A cumulative deposit plot exists for this point
            write(60,16)Airport_Name(i), &
                        Airport_CloudArrivalTime(i), &
                        cloud_morethan, Airport_CloudDuration(i), &
                        Airport_AshArrivalTime(i), &
                        deposit_morethan, Airport_AshDuration(i), &
                        Airport_Thickness(i), Airport_Thickness(i)/25.4, &
                        dp_pngfile,&
                        xmlTimeStart, airlon, airlat
          else
            write(60,6) Airport_Name(i), &
                        Airport_CloudArrivalTime(i), &
                        cloud_morethan, Airport_CloudDuration(i), &
                        Airport_AshArrivalTime(i), &
                        deposit_morethan, Airport_AshDuration(i), &
                        Airport_Thickness(i), Airport_Thickness(i)/25.4, &
                        xmlTimeStart, airlon, airlat
          endif
          nWrittenOut = nWrittenOut + 1
        elseif (Airport_CloudArrived(i)) then       !if a cloud arrived but no deposit
          cloud_morethan   = ' '
          airlon = Airport_Longitude(i)
          airlat = Airport_Latitude(i)
          CloudTime = SimStartHour+Airport_CloudArrivalTime(i)
          xmlTimeStart = HS_xmltime(CloudTime+OutputOffset,&
                                    BaseYear,useLeap)
          !See whether cloud is still overhead, or whether ash is still
          !falling
          if (CloudLoad(Airport_i(i),Airport_j(i)).gt.CLOUDLOAD_THRESH) then
            Airport_CloudDuration(i) = time-Airport_CloudArrivalTime(i)
            cloud_morethan = '>'
          else
            cloud_morethan = ' '
          endif
          xmlTimeEnd   = HS_xmltime(CloudTime+Airport_CloudDuration(i)+OutputOffset,&
                                    BaseYear,useLeap)
          write(60,7) Airport_Name(i), &
                      Airport_CloudArrivalTime(i), &
                      cloud_morethan, Airport_CloudDuration(i), &
                      xmlTimeStart, xmlTimeEnd, airlon, airlat
          nWrittenOut = nWrittenOut + 1
        endif
      enddo
      if (IsLatLon.eqv..False.) then      !Put a placemark at the location of the volcano
        call PJ_proj_inv(x_volcano, y_volcano,  &
                      A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
                      A3d_k0_scale,A3d_radius_earth, &
                      lon_volcano,lat_volcano)
      endif

      allocate(iyear(neruptions))         !Generate data needed to write ESP's
      allocate(imonth(neruptions))
      allocate(iday(neruptions))
      allocate(StartHour(neruptions))
      write(60,8) VolcanoName, VolcanoName
      !  Write table of ESP's
      do ierup=1,neruptions
        yyyymmddhh = HS_yyyymmddhh_since(e_StartTime(ierup)+SimStartHour+OutputOffset,&
                                         BaseYear,useLeap)
        read(yyyymmddhh,200) iyear(ierup),imonth(ierup),iday(ierup), &
                      StartHour(ierup)
        write(60,9) ierup,iyear(ierup),imonth(ierup),iday(ierup), &
                      StartHour(ierup), &
                      e_PlumeHeight(ierup),e_Duration(ierup),e_Volume(ierup)
200     format(i4,i2,i2,f5.2)
      enddo
      if (nWrittenOut.gt.0) then              !If one or more airports were affected
        write(60,10) lon_volcano, lat_volcano
      else                                   !If no airports were affected
        write(60,11) lon_volcano, lat_volcano
      endif
      write(60,12)                           !write final lines of file
      close(60)                             !close file

      write(global_info,4) nWrittenOut               !Write number of airports affected to log file & stdout
      write(global_log ,4) nWrittenOut

      return

!     Error traps
2001  write(global_info,*)  'Error opening ash_arrivaltimes_airports.kml.  Program stopped.'
      write(global_info,*)  'Error opening ash_arrivaltimes_airports.kml.  Program stopped.'
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

!******************************************************************************

      subroutine Close_KML(ivar,TS_flag)

      implicit none

      integer :: ivar
      integer :: TS_flag

      integer :: fid

      fid = KML_fid(ivar)

      write(global_info,15)KML_filename(ivar)
      write(global_log ,15)KML_filename(ivar)

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

!******************************************************************************

      subroutine PlotModelBoundary(xleft,xright,ybottom,ytop,ifile)

      !subroutine that finds the latitude & longitude of points in all sides of the model region

      use mesh,          only : &
         A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2,&
         A3d_k0_scale,A3d_radius_earth,IsLatLon

      use projection,    only : &
           PJ_proj_inv

      use Source,        only : &
           lon_volcano

      implicit none

      real(kind=ip)  :: xplot(0:40),yplot(0:40),lonplot(0:40),latplot(0:40)
      real(kind=ip)  :: xleft,xright,ybottom,ytop
      integer        :: ict, ifile
      !character(len=1) :: answer                     !for debugging

      write(ifile,3) ! write style for model boundary

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

      !Convert points to latitude & longitude
      if (IsLatLon) then
        lonplot = xplot
        latplot = yplot
       else
        do ict=0,40
          call PJ_proj_inv(xplot(ict),yplot(ict),  &
                        A3d_iprojflag, A3d_lam0,A3d_phi0,A3d_phi1,A3d_phi2, &
                        A3d_k0_scale,A3d_radius_earth, &
                        lonplot(ict),latplot(ict))
        enddo
      endif

      !Make sure the longitude is within 180 degrees of lon_volcano
      do ict=0,40
        if ((lon_volcano-lonplot(ict)).gt.180.0_ip)  lonplot(ict) = lonplot(ict)-360.0_ip
        if ((lon_volcano-lonplot(ict)).lt.-180.0_ip)  lonplot(ict) = lonplot(ict)+360.0_ip
        !write(6,*) 'ict=',ict,', lonplot(ict)=',lonplot(ict)
      enddo

      !write out the polygon
      write(ifile,5) (lonplot(ict),latplot(ict), ict=0,40)
      !>>>> for debugging
      !write(6,*) 'Plotting model boundary.  Continue?'
      !read(5,'(a1)') answer
      !if (answer.eq.'n') stop
      !<<<<<<

      !Format statements
3     format('	  <Style id="boundary_style">',/, &             !style for model boundary
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
!11    format(' </Document>',/,'</kml>')

      return

      end subroutine PlotModelBoundary

!******************************************************************************
      function month(imonth)

!     function that returns the name of the month given the month number

      implicit none
      character (len=2) :: imonth
      character (len=3) :: month

      if (imonth.eq.'01') then
        month='jan'
      else if (imonth.eq.'02') then
        month='feb'
      else if (imonth.eq.'03') then
        month='mar'
      else if (imonth.eq.'04') then
        month='apr'
      else if (imonth.eq.'05') then
        month='may'
      else if (imonth.eq.'06') then
        month='jun'
      else if (imonth.eq.'07') then
        month='jul'
      else if (imonth.eq.'08') then
        month='aug'
      else if (imonth.eq.'09') then
        month='sep'
      else if (imonth.eq.'10') then
        month='oct'
      else if (imonth.eq.'11') then
        month='nov'
      else if (imonth.eq.'12') then
        month='dec'
      else
        month='xxx'
      endif

      return

      end function month

      end module Output_KML
