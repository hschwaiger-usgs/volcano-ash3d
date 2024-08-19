#!/bin/bash
echo     "-----------------------------------------------------------"
echo "RUNNING TEST CASE 3: Simple sources and fall models"
echo     "-----------------------------------------------------------"
Ash3d="../../bin/Ash3d"
Ash3d_ASCII_check="../../bin/tools/Ash3d_ASCII_check"
tol=0.01
n2Dfiles=8
ascii2Doutfiles=("CloudHeight_005.00hrs.dat" "CloudHeight_010.00hrs.dat" "CloudLoad_005.00hrs.dat" "CloudLoad_010.00hrs.dat" "CloudConcentration_005.00hrs.dat" "CloudConcentration_010.00hrs.dat" "CloudArrivalTime.dat" "DepositFile_005.00hrs.dat" "DepositFile_010.00hrs.dat" "DepositFile_____final.dat" "DepositArrivalTime.dat")

nSubCases=11  # 0            1       2                3           4                5                6         7            8                     9                    10
SubCaseLabels=("Suzuki (4)" "line" "point (tracer)" "point (WH)" "point (Ganser)" "point (Stokes)" "profile" "line (topo)" "point (WH) dz_plin" "point (WH) dz_clog" "point (WH) dz_cust")

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

./clean.sh
for (( s=0;s<nSubCases;s++))
do
  echo     "-----------------------------------------------------------"
  echo "   Sub-case ${s} : ${SubCaseLabels[s]}"
  outdir="output${s}"

  ASH3DHOME=../../ ${Ash3d} TC3_XY_MSH_SC${s}.inp > /dev/null 2>&1
  rc=$((rc + $?))
  if [[ "$rc" -gt 0 ]] ; then
    echo "Error: Ash3d returned error code"
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
