#!/bin/bash
echo     "-----------------------------------------------------------"
echo "RUNNING TEST CASE 4: MSH with NCEP winds"
echo     "-----------------------------------------------------------"
Ash3d="../../bin/Ash3d"
Ash3d_ASCII_check="../../bin/tools/Ash3d_ASCII_check"
n2Dfiles=8
ascii2Doutfiles=("CloudHeight_008.00hrs.dat" "CloudHeight_016.00hrs.dat" "CloudLoad_008.00hrs.dat" "CloudLoad_016.00hrs.dat" "CloudConcentration_008.00hrs.dat" "CloudConcentration_016.00hrs.dat" "CloudArrivalTime.dat" "DepositFile_008.00hrs.dat" "DepositFile_017.00hrs.dat" "DepositFile_____final.dat" "DepositArrivalTime.dat")

ln -s /data/WindFiles/NCEP Wind_nc

nSubCases=4   # 0       1              2         3                  4
SubCaseLabels=("Cloud" "umbrella_air" "deposit" "deposit umbrella" "cloud (global)")

./clean.sh
for (( s=0;s<nSubCases;s++))
do
  echo "   Sub-case ${s} : ${SubCaseLabels[s]}"
  outdir="output${s}"

  ${Ash3d} TC4_LL_MSH_SC${s}.inp > /dev/null 2>&1
  for (( i=0;i<n2Dfiles;i++))
  do
    echo Checking 2d ASCII file "${ascii2Doutfiles[i]}"
    ${Ash3d_ASCII_check} ${ascii2Doutfiles[i]} ${outdir}/${ascii2Doutfiles[i]}
  done
  ./clean.sh
done

# Last sub-case is extra long so evaluate that separately
s=4
echo "   Sub-case ${s} : ${SubCaseLabels[s]}"
outdir="output${s}"

${Ash3d} TC4_LL_MSH_SC${s}.inp > /dev/null 2>&1
${Ash3d_ASCII_check} CloudConcentration_120.00hrs.dat output4/CloudConcentration_120.00hrs.dat
${Ash3d_ASCII_check} CloudConcentration_240.00hrs.dat output4/CloudConcentration_240.00hrs.dat
${Ash3d_ASCII_check} CloudLoad_120.00hrs.dat output4/CloudLoad_120.00hrs.dat
${Ash3d_ASCII_check} CloudLoad_240.00hrs.dat output4/CloudLoad_240.00hrs.dat
${Ash3d_ASCII_check} CloudHeight_120.00hrs.dat output4/CloudHeight_120.00hrs.dat
${Ash3d_ASCII_check} CloudHeight_240.00hrs.dat output4/CloudHeight_240.00hrs.dat
${Ash3d_ASCII_check} DepositFile_120.00hrs.dat output4/DepositFile_120.00hrs.dat
${Ash3d_ASCII_check} DepositFile_240.00hrs.dat output4/DepositFile_240.00hrs.dat
./clean.sh

