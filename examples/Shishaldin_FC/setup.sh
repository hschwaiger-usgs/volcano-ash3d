#!/bin/bash
if [ -z ${WINDROOT} ];then
 # default location
 WINDROOT="/data/WindFiles"
fi
rc=0
ls -1r ${WINDROOT}/NWP_testfiles/lastdownload.txt
rc=$((rc + $?))
if [[ "$rc" -gt 0 ]] ; then
  echo "Error: Could not find the last NWP_testfile package (lastdownload.txt)"
  echo "To download the NWP Forecast data, run:"
  echo "/opt/USGS/bin/autorun_scripts/autorun_scripts/get_NCEP_50YearReanalysis.sh 1992"
  exit
fi
YYYYMMDD=`cat "${WINDROOT}/NWP_testfiles/lastdownload.txt"`
#ln -s ${WINDROOT}/NWP_testfiles/Forecasts/GRIB/${YYYYMMDD} WindGribFC
ln -s ${WINDROOT}/NWP_testfiles/Forecasts/NetCDF/${YYYYMMDD} WindNetcdfFC
ln -s ${WINDROOT}/NWP_testfiles/Reanalysis/NetCDF/${YYYYMMDD} WindNetcdfRN
ln -s /data/TOPO/GEBCO/GEBCO_23/GEBCO_2023.nc .
ln -s ../../bin/Ash3d Ash3d
