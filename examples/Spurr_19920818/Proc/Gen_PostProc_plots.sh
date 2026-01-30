#!/bin/bash

# clean up
rm depothik.* Ash3d_CloudLoad_t0*png Ash3d_Deposit____final.png plplt_000*png

ln -s ../../../bin/Ash3d_PostProc .
ln -s ../3d_tephra_fall.nc .
infile=3d_tephra_fall.nc

# Generate vertical profile plots using plplot (change to ASH3DPLOT=3 for gnuplot)
ASH3DPLOT=2 ./Ash3d_PostProc ${infile} 16 3

# Generate shapefile output of deposit thickness in mm
ASH3DPLOT=3 ./Ash3d_PostProc ${infile} 5 5

# Generate map of deposit with custom contour levels (GMT) using a control file
./Ash3d_PostProc pp.ctr

# Generate animation of cloud load using gnuplot
tmax=`ncdump     -h ${infile} | grep "t = UNLIMITED" | grep -v pt | cut -c22-23` # maximum time dimension
for t in `seq 0 $((tmax-1))`;
do
  ASH3DPLOT=3 ./Ash3d_PostProc ${infile} 12 3 $t
done
convert -delay 25 -loop 0 `ls -1tr Ash3d_CloudLoad_t*.png` CloudLoad_animation.gif

