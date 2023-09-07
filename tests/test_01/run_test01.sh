#!/bin/bash
echo     "-----------------------------------------------------------"
echo "RUNNING TEST CASE 1: 2D-ADVECTION"
echo     "-----------------------------------------------------------"
Ash3d="../../bin/Ash3d"
Ash3d_ASCII_check="../../bin/tools/Ash3d_ASCII_check"
n2Dfiles=7
ascii2Doutfiles=("CloudHeight_002.00hrs.dat" "CloudHeight_004.00hrs.dat" "CloudLoad_002.00hrs.dat" "CloudLoad_004.00hrs.dat" "CloudConcentration_002.00hrs.dat" "CloudConcentration_004.00hrs.dat" "CloudArrivalTime.dat")

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

outdir="output"

${Ash3d} TC1_XY_MSH.inp > /dev/null 2>&1
rc=$((rc + $?))
if [[ "$rc" -gt 0 ]] ; then
  echo "Error: Ash3d returned error code"
  exit 1
fi

for (( i=0;i<n2Dfiles;i++))
do
  echo Checking 2d ASCII file "${ascii2Doutfiles[i]}"
  stat=`${Ash3d_ASCII_check} ${ascii2Doutfiles[i]} ${outdir}/${ascii2Doutfiles[i]} | cut -f1 -d':'`\

  if [[ $stat == *PASS* ]];
  then
    printf " ---> ${GREEN}${stat}${NC}\n"
  else
    printf " ---> ${RED}${stat}${NC}\n"
  fi

done
./clean.sh
