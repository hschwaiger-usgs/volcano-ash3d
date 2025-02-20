#!/bin/bash
if [ -z ${WINDROOT} ];then
 # default location
 WINDROOT="/data/WindFiles"
fi
rc=0
ls -1r ${WINDROOT}/NCEP/2008/air.2008.nc
rc=$((rc + $?))
if [[ "$rc" -gt 0 ]] ; then
  echo "Error: Could not find NCEP data for 2008"
  echo "To download the NCEP data, run:"
  echo "/opt/USGS/bin/autorun_scripts/autorun_scripts/get_NCEP_50YearReanalysis.sh 2008"
  exit
fi
ln -s ${WINDROOT}/NCEP .
ln -s ../../bin/Ash3d Ash3d
