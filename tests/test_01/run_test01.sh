#!/bin/bash
echo     "-----------------------------------------------------------"
echo "RUNNING TEST CASE 1: 2D-ADVECTION"
echo     "-----------------------------------------------------------"
Ash3d="../../bin/Ash3d"
Ash3d_ASCII_check="../../bin/tools/Ash3d_ASCII_check"
Ash3d_PP="../../bin/Ash3d_PostProc"
GenPlots=1
tol=0.01
n2Dfiles=7
ascii2Doutfiles=("CloudHeight_002.00hrs.dat" "CloudHeight_004.00hrs.dat" "CloudLoad_002.00hrs.dat" "CloudLoad_004.00hrs.dat" "CloudConcentration_002.00hrs.dat" "CloudConcentration_004.00hrs.dat" "CloudArrivalTime.dat")

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

outdir="output"

ASH3DHOME=../../ ${Ash3d} TC1_XY_MSH.inp > /dev/null 2>&1
rc=$((rc + $?))
if [[ "$rc" -gt 0 ]] ; then
  echo "Error: Ash3d returned error code"
  exit 1
fi

if [[ "$GenPlots" -gt 0 ]] ; then
  # Post-processing
  ASH3DHOME=../../ ${Ash3d_PP} TC1_pp_CL2.ctr > /dev/null 2>&1
  mv Ash3d_CloudLoad_t002.00hrs.png TC1_Ash3d_CloudLoad_t002.00hrs.png
  ASH3DHOME=../../ ${Ash3d_PP} TC1_pp_CL4.ctr > /dev/null 2>&1
  mv Ash3d_CloudLoad_t004.00hrs.png TC1_Ash3d_CloudLoad_t004.00hrs.png
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
