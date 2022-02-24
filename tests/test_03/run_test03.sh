#!/bin/bash
echo     "-----------------------------------------------------------"
echo "RUNNING TEST CASE 3: Simple sources and fall models"
echo     "-----------------------------------------------------------"
Ash3d="../../bin/Ash3d"
Ash3d_ASCII_check="../../bin/tools/Ash3d_ASCII_check"
n2Dfiles=8
ascii2Doutfiles=("CloudHeight_005.00hrs.dat" "CloudHeight_010.00hrs.dat" "CloudLoad_005.00hrs.dat" "CloudLoad_010.00hrs.dat" "CloudConcentration_005.00hrs.dat" "CloudConcentration_010.00hrs.dat" "CloudArrivalTime.dat" "DepositFile_005.00hrs.dat" "DepositFile_010.00hrs.dat" "DepositFile_____final.dat" "DepositArrivalTime.dat")

nSubCases=8   # 0            1       2                3           4                5                6         7
SubCaseLabels=("Suzuki (4)" "line" "point (tracer)" "point (WH)" "point (Ganser)" "point (Stokes)" "profile" "line (topo)")

./clean.sh
for (( s=0;s<nSubCases;s++))
do
  echo "   Sub-case ${s} : ${SubCaseLabels[s]}"
  outdir="output${s}"

  ${Ash3d} TC3_XY_MSH_SC${s}.inp > /dev/null 2>&1
  for (( i=0;i<n2Dfiles;i++))
  do
    echo Checking 2d ASCII file "${ascii2Doutfiles[i]}"
    ${Ash3d_ASCII_check} ${ascii2Doutfiles[i]} ${outdir}/${ascii2Doutfiles[i]}
  done
  ./clean.sh
done
