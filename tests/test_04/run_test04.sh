#!/bin/bash
echo     "-----------------------------------------------------------"
echo "RUNNING TEST CASE 4: MSH with NCEP winds"
echo     "-----------------------------------------------------------"
Ash3d="../../bin/Ash3d"
Ash3d_ASCII_check="../../bin/tools/Ash3d_ASCII_check"
tol=0.01
n2Dfiles=11
ascii2Doutfiles=("CloudHeight_008.00hrs.dat" "CloudHeight_016.00hrs.dat" "CloudLoad_008.00hrs.dat" "CloudLoad_016.00hrs.dat" "CloudConcentration_008.00hrs.dat" "CloudConcentration_016.00hrs.dat" "CloudArrivalTime.dat" "DepositFile_008.00hrs.dat" "DepositFile_016.00hrs.dat" "DepositFile_____final.dat" "DepositArrivalTime.dat")

nSubCases=4   # 0       1              2         3                  4                5           6           7
SubCaseLabels=("Cloud" "umbrella_air" "deposit" "deposit umbrella" "cloud (global)" "Topo ID=0" "Topo ID=1" "Topo ID=2")

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

rc=0
rm -f Wind_nc
./clean.sh

if [ -z ${WINDROOT} ];then
 # Standard Linux location
 WINDROOT="/data/WindFiles"
 # Mac
 #WINDROOT="/opt/data/WindFiles"
fi
# Check to see if the NCEP data for 1980 is present
ls -1r ${WINDROOT}/NCEP/1980/air.1980.nc
rc=$((rc + $?))
if [[ "$rc" -gt 0 ]] ; then
  echo "Error: Could not find NCEP data for 1980"
  echo "To download the NCEP data, run:"
  echo "/opt/USGS/bin/autorun_scripts/autorun_scripts/get_NCEP_50YearReanalysis.sh 1980"
  exit
fi

for (( s=0;s<nSubCases;s++))
do
  echo     "-----------------------------------------------------------"
  echo "   Sub-case ${s} : ${SubCaseLabels[s]}"
  outdir="output${s}"
  ln -s ${WINDROOT}/NCEP Wind_nc
  ASH3DHOME=../../ ${Ash3d} TC4_LL_MSH_SC${s}.inp > /dev/null 2>&1
  rc=$((rc + $?))
  if [[ "$rc" -gt 0 ]] ; then
    echo "Error: Ash3d returned error code"
    exit 1
  fi
  grep "useVz_rhoG=.true." Ash3d.lst > /dev/null
  rc=$((rc + $?))
  if [[ "$rc" -gt 0 ]] ; then
    echo "Error: Ash3d was not compiled with useVz_rhoG=.true."
    echo "       Vz is calculated via finite-differncing dp/dz."
    echo "       Results may still be valid, but this script with report failures."
    exit 1
  fi
  for (( i=0;i<n2Dfiles;i++))
  do
    echo Checking 2d ASCII file "${ascii2Doutfiles[i]}"
    stat=`${Ash3d_ASCII_check} ${ascii2Doutfiles[i]} ${outdir}/${ascii2Doutfiles[i]} ${tol} | cut -f1 -d':'`\

    if [[ $stat == *PASS* ]];
    then
      printf " ---> ${GREEN}${stat}${NC}\n"
    else
      printf " ---> ${RED}${stat}${NC}\n"
    fi

  done
  ./clean.sh
done

echo     "-----------------------------------------------------------"
# Last NCEP sub-case is extra long so evaluate that separately since the hour-stamp is different than above
s=4
echo "   Sub-case ${s} : ${SubCaseLabels[s]}"
outdir="output${s}"
ascii2Doutfiles2=("CloudHeight_120.00hrs.dat" "CloudHeight_240.00hrs.dat" "CloudLoad_120.00hrs.dat" "CloudLoad_240.00hrs.dat" "CloudConcentration_120.00hrs.dat" "CloudConcentration_240.00hrs.dat" "CloudArrivalTime.dat" "DepositFile_120.00hrs.dat" "DepositFile_240.00hrs.dat" "DepositFile_____final.dat" "DepositArrivalTime.dat")

ln -s ${WINDROOT}/NCEP Wind_nc

ASH3DHOME=../../ ${Ash3d} TC4_LL_MSH_SC${s}.inp > /dev/null 2>&1
for (( i=0;i<n2Dfiles;i++))
do
  echo Checking 2d ASCII file "${ascii2Doutfiles2[i]}"
  stat=`${Ash3d_ASCII_check} ${ascii2Doutfiles2[i]} ${outdir}/${ascii2Doutfiles2[i]} ${tol} | cut -f1 -d':'`\

  if [[ $stat == *PASS* ]];
  then
    printf " ---> ${GREEN}${stat}${NC}\n"
  else
    printf " ---> ${RED}${stat}${NC}\n"
  fi
done
./clean.sh

echo     "-----------------------------------------------------------"
# We have 3 topography test cases at the end
#  5: ZScaling_ID = 0 (No topo)
#  6: ZScaling_ID = 1 (shifted topo) should be same model as for 5, but at a plateau of z=1.0 km
#  7: ZScaling_ID = 2 (scaled topo) same source as 5/6 at plateau of z=1.0, but different grid
for (( s=5;s<8;s++))
do
  echo     "-----------------------------------------------------------"
  echo "   Sub-case ${s} : ${SubCaseLabels[s]}"
  outdir="output${s}"
  ASH3DHOME=../../ ${Ash3d} TC4_LL_MSH_SC${s}.inp > /dev/null 2>&1
  rc=$((rc + $?))
  if [[ "$rc" -gt 0 ]] ; then
    echo "Error: Ash3d returned error code"
    exit 1
  fi
  grep "useVz_rhoG=.true." Ash3d.lst > /dev/null
  rc=$((rc + $?))
  if [[ "$rc" -gt 0 ]] ; then
    echo "Error: Ash3d was not compiled with useVz_rhoG=.true."
    echo "       Vz is calculated via finite-differncing dp/dz."
    echo "       Results may still be valid, but this script with report failures."
    exit 1
  fi
  echo "Checking 2d ASCII file DepositFile_____final.dat"
  stat=`${Ash3d_ASCII_check} DepositFile_____final.dat ${outdir}/DepositFile_____final.dat ${tol} | cut -f1 -d':'`\

  if [[ $stat == *PASS* ]];
  then
    printf " ---> ${GREEN}${stat}${NC}\n"
  else
    printf " ---> ${RED}${stat}${NC}\n"
  fi

  ./clean.sh
done


