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

# Mazama coordinates
vlt=42.93
vln=-122.12

lonw=225.0
lone=269.0
lats=33.0
latn=57.0
DETAIL="-Dl"
BASE="-Bg5/g5 -P"
PROJ="-JM${vln}/${vlt}/7.5i"
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
#datafile="../Data/Mazama_DepThick_mm.dat"
datafile="../Data/thickness_expanded.txt"

# Ash3d deposit file (Run 22) from Buckland et al. 2022
gmt grdconvert ../Data/MZRun022.txt=ef out.grd
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
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_0.03.lev -A- -W5,90/0/90,-   -O -K >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_0.1.lev  -A- -W5,128/0/128,- -O -K >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_0.3.lev  -A- -W5,0/0/255,-   -O -K >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_1.lev    -A- -W5,0/128/255,- -O -K >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_3.lev    -A- -W5,0/255/128,- -O -K >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_10.lev   -A- -W5,195/195/0,- -O -K >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_30.lev   -A- -W5,255/128/0,- -O -K >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE -Cdpm_100.lev  -A- -W5,255/0/0,-   -O -K >> temp.ps

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
LEGLOC="-Dx5.4i/0.10i/2.0i/1.4i/BL"
gmt pslegend $AREA $PROJ $BASE -G255 $LEGLOC -K -O << EOF >> temp.ps
C black
H 14 1 Mazama Umb. Dep.
D 1p
S 0.1i - 0.15i black  0.5p,black 0.3i Run22 (Buckland,2022)
S 0.1i - 0.15i red    3.0p,red   0.3i Ash3d (all colors)
S 0.1i c 0.10i red  0.5p,black   0.3i measured thickness
EOF
gmt psscale -Dx6.25i/0.55i/1.75i/0.15ih -C$CPT -Q -B10f5/:"mm": -O -K >> temp.ps

# Plot the tephra site data
# First, reformat data file to something more easily ingested by psxy
cat ${datafile} | awk '{print $2,$1,$3*10.0}' | tail -n +2 > dep.dat
gmt psxy dep.dat $AREA $PROJ -Sc0.1i -C${CPT} -Wthinnest -O -K >> temp.ps

# Last gmt command is to plot the volcano and close out the ps file
echo $vln $vlt '1.0' | gmt psxy $AREA $PROJ -St0.1i -Gblack -Wthinnest -O >> temp.ps

# Save map
ps2epsi temp.ps temp.eps
convert temp.eps Mazama_UmbrellaDeposit.png
epstopdf temp.eps Mazama_UmbrellaDeposit.pdf

# Clean up
rm temp.* dpm*lev gmt.history gmt.conf dep.cpt zero.grd out.grd dep.dat

