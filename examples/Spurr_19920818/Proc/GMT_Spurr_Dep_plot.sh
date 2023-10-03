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
echo "GMT version = ${GMTv}"

# Spurr coordinates
vlt=61.299
vln=-152.251

lonw=-154.0
lone=-144.0
lats=60.0
latn=63.0
DETAIL="-Dh"
BASE="-Bg2/g1 -P"
PROJ=-JS-150.0/90/8i
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
infile="../3d_tephra_fall.nc"
if test -r ${infile} ; then
    echo "Preparing to read from ${infile} file"
  else
    echo "error: no ${infile} file. Exiting"
    exit 1
fi
# Preparing grid file
dep_grd=var_out_final.grd
gmt grdconvert "$infile?area" zero.grd
gmt grdmath 0.0 zero.grd MUL = zero.grd
# We need to convert the NaN's to zero to get the lowest contour
gmt grdconvert "$infile?depothickFin" temp.grd
gmt grdmath temp.grd zero.grd AND = temp.grd

# Deposit data file
datafile="../Data/Spurr_19920818_DepThick_mm.dat"
#******************************************************************************

# Create Base Map
gmt pscoast $AREA $PROJ $BASE $DETAIL $COAST -S255/255/255 -K  > temp.ps

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
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_0.03.lev -A- -W3,90/0/90   -O -K >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_0.1.lev  -A- -W3,128/0/128 -O -K >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_0.3.lev  -A- -W3,0/0/255   -O -K >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_1.lev    -A- -W3,0/128/255 -O -K >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_3.lev    -A- -W3,0/255/128 -O -K >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_10.lev   -A- -W3,195/195/0 -O -K >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_30.lev   -A- -W3,255/128/0 -O -K >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_100.lev  -A- -W3,255/0/0   -O -K >> temp.ps

# Plot legend
LEGLOC="-Dx0.1i/4.1i/3.0i/1.0i/BL"
gmt pslegend $AREA $PROJ $BASE -G255 $LEGLOC -K -O << EOF >> temp.ps
C black
H 14 1 Spurr August 18, 1992
D 1p
EOF
gmt psscale -Dx1.25i/4.6i/2i/0.15ih -C$CPT -Q -B10f5/:"mm": -O -K >> temp.ps
# Plot the tephra site data
gmt psxy ${datafile} $AREA $PROJ -Sc0.1i -C${CPT} -Wthinnest -O -K >> temp.ps

# Last gmt command is to plot the volcano and close out the ps file
echo $vln $vlt '1.0' | gmt psxy $AREA $PROJ -St0.1i -Gblack -Wthinnest -O >> temp.ps

# Save map
ps2epsi temp.ps temp.eps
convert temp.eps Spurr_Deposit.png
epstopdf temp.eps Spurr_Deposit.pdf

# Clean up
rm temp.* dpm*lev gmt.history gmt.conf dep.cpt zero.grd

