#!/bin/bash

# Set RunClass to one of:
# PJ_091 PJ_216 PJ_221 PJ_242
# LL_091 LL_216 LL_221 LL_242
# LL_gfs LL_ecmwf LL_NASA
# LL_NARR LL_ERA5 LL_MERRA LL_NCEP
# PJ_NARR
RunClass="PJ_NARR"
mkdir -p ${RunClass}

ln -s ../../../bin/Ash3d_PostProc .
ln -s ../3d_tephra_fall_${RunClass}.nc intemp.nc
infile=intemp.nc

# Generate vertical profile plots using plplot (change to ASH3DPLOT=3 for gnuplot)
ASH3DPLOT=3 ./Ash3d_PostProc ${infile} 16 3
mv gnupl_0001.png ${RunClass}/gnupl_0001.png
mv gnupl_0002.png ${RunClass}/gnupl_0002.png
mv gnupl_0003.png ${RunClass}/gnupl_0003.png

# Generate png output of deposit thickness in mm
ASH3DPLOT=4 ./Ash3d_PostProc ${infile} 5 3
mv Ash3d_Deposit____final.png ${RunClass}/Ash3d_Deposit____final.png
rm intemp.nc

