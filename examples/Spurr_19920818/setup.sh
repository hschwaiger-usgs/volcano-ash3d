#!/bin/bash
if [ -z ${WINDROOT} ];then
 # Standard Linux location
 WINDROOT="/data/WindFiles"
 # Mac
 #WINDROOT="/opt/data/WindFiles"
fi
rc=0
ls -1r ${WINDROOT}/NCEP/1992/air.1992.nc
rc=$((rc + $?))
if [[ "$rc" -gt 0 ]] ; then
  echo "Error: Could not find NCEP data for 1992"
  echo "To download the NCEP data, run:"
  echo "/opt/USGS/bin/autorun_scripts/autorun_scripts/get_NCEP_50YearReanalysis.sh 1992"
  exit
fi
ln -s ${WINDROOT}/NCEP .
ln -s ../../bin/Ash3d Ash3d
