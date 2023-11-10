#!/bin/bash
if [ -z ${WINDROOT} ];then
 # Standard Linux location
 WINDROOT="/data/WindFiles"
 # Mac
 #WINDROOT="/opt/data/WindFiles"
fi
ln -s ../../bin/Ash3d Ash3d
ln -s ${WINDROOT}/NCEP .
