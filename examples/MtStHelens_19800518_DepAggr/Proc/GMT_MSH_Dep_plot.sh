#!/bin/bash
# GMT Post-processing script for the Ash3d Spurr deposit validation test.

# We need to know if we must prefix all gmt commands with 'gmt', as required by version 5/6
GMTv=5
type gmt >/dev/null 2>&1 || { echo >&2 "Command 'gmt' not found.  Assuming GMTv4."; GMTv=4;}
if [ $GMTv -eq 4 ] ; then
    echo "GMT 4 is no longer supported."
    echo "Please update to GMT 5 or 6"
    exit 1
 else
    GMTv=`gmt --version | cut -c1`
fi
GMTpen=("-" "-" "-" "-" "/" ",")
#echo "GMT version = ${GMTv}"

# MSH coordinates
vlt=46.20
vln=237.82

# Full domain of tephra samples
#lonw=235.0
#lone=255.0
#lats=45.0
#latn=51.0

# Near-field for comparison with historic Ash3d run
lonw=235.0
lone=243.0
lats=44.9
latn=48.0
DETAIL="-Dh"
BASE="-Bg1/g1 -P"
PROJ="-JM${vln}/${vlt}/8i"
AREA="-R${lonw}/${lats}/${lone}/${latn}r"
COAST="-G220/220/220 -W"

# Set up some default values
gmt gmtset PROJ_ELLIPSOID Sphere
gmt gmtset PAPER_MEDIA=Custom_720x510

#############################################################################
# Color map for the deposit sample thicknesses
CPT="dep.cpt"
echo "# file dep.cpt
#COLOR_MODEL = RGB
0.03    90      0       90      0.1     90      0       90
0.1	128	0	128	0.3	128	0	128
0.3	0	0	255	1.0	0	0	255
1.0	0	128	255	3.0	0	128	255
3.0	0	255	255	10.0	0	255	255
10.0	128	255	128	30.0	128	255	128
30.0	255	255	0	100.0	255	255	0
100.0	255	128	0	300.0	255	128	0
300.0	255	0	0	1000.0	255	0	0" > ${CPT}

#############################################################################
# Ash3d data file
# Preparing grid file
gmt grdconvert ../DepositFile_____final.dat=ef temp.grd

# Deposit data file
datafile="../Data/sample_mpua.xy"

# Ash3d deposit file (phi 1.0 rho2624)
gmt grdconvert ../Data/DepositFile_____final.dat=ef out.grd
#******************************************************************************

# Create Base Map
gmt pscoast $AREA $PROJ $BASE $DETAIL $COAST -S100/149/237 -K  > temp.ps

# Contour Ash3d output
echo "0.01   C" > dpm_0.01.lev   #deposit (0.01 mm)
echo "0.03   C" > dpm_0.03.lev   #deposit (0.03 mm)
echo "0.1    C"  > dpm_0.1.lev   #deposit (0.1 mm)
echo "0.3    C"  > dpm_0.3.lev   #deposit (0.3 mm)
echo "1.0    C"  >   dpm_1.lev   #deposit (1 mm)
echo "3.0    C"  >   dpm_3.lev   #deposit (3 mm)
echo "10.0   C"  >  dpm_10.lev   #deposit (1 cm)
echo "30.0   C"  >  dpm_30.lev   #deposit (3 cm)
echo "100.0  C"  > dpm_100.lev   #deposit (10cm)
echo "300.0  C"  > dpm_300.lev   #deposit (30cm)
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_0.03.lev -A- -W2,90/0/90,-   -O -K >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_0.1.lev  -A- -W2,128/0/128,- -O -K >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_0.3.lev  -A- -W2,0/0/255,-   -O -K >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_1.lev    -A- -W2,0/128/255,- -O -K >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_3.lev    -A- -W2,0/255/128,- -O -K >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_10.lev   -A- -W2,195/195/0,- -O -K >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_30.lev   -A- -W2,255/128/0,- -O -K >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_100.lev  -A- -W2,255/0/0,-   -O -K >> temp.ps

# Thin-line contours of the data from the Buckland paper
gmt grdcontour out.grd $AREA $PROJ $BASE -Cdpm_0.03.lev -A- -W1,0/0/0 -O -K >> temp.ps
gmt grdcontour out.grd $AREA $PROJ $BASE -Cdpm_0.1.lev  -A- -W1,0/0/0 -O -K >> temp.ps
gmt grdcontour out.grd $AREA $PROJ $BASE -Cdpm_0.3.lev  -A- -W1,0/0/0 -O -K >> temp.ps
gmt grdcontour out.grd $AREA $PROJ $BASE -Cdpm_1.lev    -A- -W1,0/0/0 -O -K >> temp.ps
gmt grdcontour out.grd $AREA $PROJ $BASE -Cdpm_3.lev    -A- -W1,0/0/0 -O -K >> temp.ps
gmt grdcontour out.grd $AREA $PROJ $BASE -Cdpm_10.lev   -A- -W1,0/0/0 -O -K >> temp.ps
gmt grdcontour out.grd $AREA $PROJ $BASE -Cdpm_30.lev   -A- -W1,0/0/0 -O -K >> temp.ps
gmt grdcontour out.grd $AREA $PROJ $BASE -Cdpm_100.lev  -A- -W1,0/0/0 -O -K >> temp.ps

# Plot legend
LEGLOC="-Dx5.9i/0.07i/2.0i/1.4i/BL"
gmt pslegend $AREA $PROJ $BASE -G255 $LEGLOC -K -O << EOF >> temp.ps
C black
H 14 1 Mt.St.Helens Dep.
D 1p
S 0.1i - 0.15i black  0.5p,black 0.3i Historic run
S 0.1i - 0.15i red    3.0p,red   0.3i Ash3d (all colors)
S 0.1i c 0.10i red  0.5p,black   0.3i measured thickness
EOF
gmt psscale -Dx6.75i/0.5i/1.75i/0.15ih -C$CPT -Q -B10f5/:"mm": -O -K >> temp.ps

# Plot the tephra site data
# First, reformat data file to something more easily ingested by psxy
cat ${datafile} | awk '{print $1,$2,$3*1.0}' > dep.dat
gmt psxy dep.dat $AREA $PROJ -Sc0.1i -C${CPT} -Wthinnest -O -K >> temp.ps

# Last gmt command is to plot the volcano and close out the ps file
echo $vln $vlt '1.0' | gmt psxy $AREA $PROJ -St0.1i -Gmagenta -Wthinnest -O >> temp.ps

# Save map
#ps2epsi temp.ps temp.eps
#convert temp.eps MSH_Deposit.png
#epstopdf temp.eps MSH_Deposit.pdf
ps2pdf temp.ps MSH_Deposit.pdf

# Clean up
rm -f temp.* dpm*lev gmt.history gmt.conf dep.cpt out.grd dep.dat

