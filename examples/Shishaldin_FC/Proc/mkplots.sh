#!/bin/bash

WF=("091" "216" "221" "242" "ecmwf" "ERA5" "gfs" "MERRA" "NARR" "NASA" "NCEP")
RF=("FC"  "FC"  "FC"  "FC"  "FC"    "RN"   "FC"  "RN"    "RN"   "FC"   "RN")

# Full VD + topo
run="Full"
wc=(" " " " " " " " " " " " " " " " " " " " " ")    # wall clock (in seconds) taken from the *lst files
#wc=("816" "748" "745" "737" "631" "620" "717" "795" "615" "676" "801")
## Kh
#run="Kh"
#wc=("760" "693" "688" "695" "601" "533" "623" "729" "590" "621" "724")
## No VD
#run="None"
#wc=("429" "327" "327" "326" "272" "243" "301" "342" "278" "299" "344")
## No VD, No Topo
#run="None_NoTopo"
#wc=("433" "274" "273" "274" "249" "238" "332" "289" "229" "263" "292")

## Processing cloud load for hour=6
psz=35
tx=150
ty=400
ln -s ../3d_tephra_fall_LL_*.nc .
for (( i=0;i<11;i++)); do
 /opt/USGS/Ash3d/bin/Ash3d_PostProc 3d_tephra_fall_LL_${WF[i]}.nc 12 3 -1
 convert -pointsize ${psz} -fill blue -draw "text ${tx},${ty} \"VD-${WF[i]}:${run}: ${wc[i]} sec\"" Ash3d_CloudLoad_t006.00hrs.png Ash3d_CloudLoad_t006.00hrs_LL_${WF[i]}.png
done

montage -tile 3x4 -geometry 600x450 \
Ash3d_CloudLoad_t006.00hrs_LL_091.png \
Ash3d_CloudLoad_t006.00hrs_LL_216.png \
Ash3d_CloudLoad_t006.00hrs_LL_221.png \
Ash3d_CloudLoad_t006.00hrs_LL_242.png \
Ash3d_CloudLoad_t006.00hrs_LL_ecmwf.png \
Ash3d_CloudLoad_t006.00hrs_LL_ERA5.png \
Ash3d_CloudLoad_t006.00hrs_LL_gfs.png \
Ash3d_CloudLoad_t006.00hrs_LL_MERRA.png \
Ash3d_CloudLoad_t006.00hrs_LL_NARR.png \
Ash3d_CloudLoad_t006.00hrs_LL_NASA.png \
Ash3d_CloudLoad_t006.00hrs_LL_NCEP.png \
Shishaldin_test_AirCL_$VarDiff_${run}.png
rm Ash3d_CloudLoad*png

## Processing profiles
psz=30
tx=100
ty=175
ip=1
for (( i=0;i<11;i++)); do
  ASH3DPLOT=3 /opt/USGS/Ash3d/bin//Ash3d_PostProc 3d_tephra_fall_LL_${WF[i]}.nc 16 3
  convert -pointsize ${psz} -fill blue -draw "text ${tx},${ty} \"${WF[i]}:${run}: ${wc[i]} sec\"" gnupl_000${ip}.png Ash3d_profile_${WF[i]}.png
done

montage -tile 3x4 -geometry 600x450 \
Ash3d_profile_091.png \
Ash3d_profile_216.png \
Ash3d_profile_221.png \
Ash3d_profile_242.png \
Ash3d_profile_ecmwf.png \
Ash3d_profile_ERA5.png \
Ash3d_profile_gfs.png \
Ash3d_profile_MERRA.png \
Ash3d_profile_NARR.png \
Ash3d_profile_NASA.png \
Ash3d_profile_NCEP.png \
Shishaldin_test_Prof_$VarDiff_${run}.png
rm Ash3d_profile*png gnupl_000*png

## Processing Deposit thickness
psz=35
tx=150
ty=400
for (( i=0;i<11;i++)); do
 ASH3DPLOT=3 /opt/USGS/Ash3d/bin/Ash3d_PostProc 3d_tephra_fall_LL_${WF[i]}.nc 5 3
 convert -pointsize ${psz} -fill blue -draw "text ${tx},${ty} \"${WF[i]}:${run}: ${wc[i]} sec\"" Ash3d_Deposit____final.png Ash3d_Deposit____final_${WF[i]}.png
done

montage -tile 3x4 -geometry 600x450 \
Ash3d_Deposit____final_091.png \
Ash3d_Deposit____final_216.png \
Ash3d_Deposit____final_221.png \
Ash3d_Deposit____final_242.png \
Ash3d_Deposit____final_ecmwf.png \
Ash3d_Deposit____final_ERA5.png \
Ash3d_Deposit____final_gfs.png \
Ash3d_Deposit____final_MERRA.png \
Ash3d_Deposit____final_NARR.png \
Ash3d_Deposit____final_NASA.png \
Ash3d_Deposit____final_NCEP.png \
Shishaldin_test_Dep_VarDiff_${run}.png
rm Ash3d_Deposit_*png

rm Ash3d_pp.log
