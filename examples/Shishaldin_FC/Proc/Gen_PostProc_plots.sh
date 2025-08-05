#!/bin/bash

# clean up
rm depothik.* Ash3d_CloudLoad_t0*png Ash3d_Deposit____final.png plplt_000*png

ln -s ../../../bin/Ash3d_PostProc .
ln -s ../3d_tephra_fall_LL_091.nc .
infile1=3d_tephra_fall_LL_091.nc
ln -s ../3d_tephra_fall_PJ_091.nc .
infile2=3d_tephra_fall_PJ_091.nc

# Generate vertical profile plots using plplot (change to ASH3DPLOT=3 for gnuplot)
ASH3DPLOT=2 ./Ash3d_PostProc ${infile1} 16 3
mv gnupl_0001.png gnupl_0001_LL.png
mv gnupl_0002.png gnupl_0002_LL.png
mv gnupl_0003.png gnupl_0003_LL.png

# Generate shapefile output of deposit thickness in mm
ASH3DPLOT=4 ./Ash3d_PostProc ${infile1} 5 3
mv Ash3d_Deposit____final.png Ash3d_Deposit____final_LL.png

# Generate vertical profile plots using plplot (change to ASH3DPLOT=3 for gnuplot)
ASH3DPLOT=2 ./Ash3d_PostProc ${infile2} 16 3
mv gnupl_0001.png gnupl_0001_PJ.png
mv gnupl_0002.png gnupl_0002_PJ.png
mv gnupl_0003.png gnupl_0003_PJ.png

# Generate shapefile output of deposit thickness in mm
ASH3DPLOT=4 ./Ash3d_PostProc ${infile2} 5 3
mv Ash3d_Deposit____final.png Ash3d_Deposit____final_PJ.png

# Generate map of deposit with custom contour levels (GMT) using a control file
#./Ash3d_PostProc pp.ctr

# Generate animation of cloud load using gnuplot
#tmax=`ncdump     -h ${infile} | grep "t = UNLIMITED" | grep -v pt | cut -c22-23` # maximum time dimension
#for t in `seq 0 $((tmax-1))`;
#do
#  ASH3DPLOT=3 ./Ash3d_PostProc ${infile} 12 3 $t
#done
#convert -delay 25 -loop 0 `ls -1tr Ash3d_CloudLoad_t*.png` CloudLoad_animation.gif

