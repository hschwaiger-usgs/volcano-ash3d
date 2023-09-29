#!/bin/bash
# GMT Post-processing script for the Ash3d Kasatochi validation test.
# You will need to unzip the compressed data files in Kasatochi_20080808_CldSat/Data

# Change the time step (i) to process here
# i=0  2008-08-08_13-40 ( 9.1 hours)
# i=1  2008-08-09_00-20 (20.2 hours)
# i=2  2008-08-09_12-45 (32.3 hours)
# i=3  2008-08-09_23-25 (43.3 hours)
# i=4  2008-08-10_11-49 (55.2 hours)
# i=5  2008-08-10_22-29 (64.6 hours)
# i=6  2008-08-11_10-53 (78.2 hours)
# i=7  2008-08-11_21-33 (89.5 hours)
#  2,3,4,6
i=2

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

# Kasatochi coordinates
vlt=52.177
vln=-175.508

lonw=179
lone=260
lats=29.0
latn=60.0
DETAIL="-Dl"
BASE1="-Bg10/g5 -P"
PROJ=-JS-150.0/90/8i
AREA="-R${lonw}/${lats}/${lone}/${latn}r"
COAST="-G220/220/220 -W"

# Set up some default values
gmt gmtset ELLIPSOID Sphere
gmt gmtset PAPER_MEDIA=Custom_720x510

#############################################################################
# Color map for the satellite ash retrievals
CPT="cload.cpt"
echo "# file cload.cpt
#COLOR_MODEL = RGB
#
0	255	234	255	1	255	234	255
1	206	192	255	2	206	192	255
2	0	255	255	3	0	255	255
3	0	255	0	4	0	255	0
4	251	145	0	5	251	145	0
B	255	255	255
F	0	0	0
N	128	128	128" > ${CPT}

#############################################################################
# Ash3d data file
infile="../3d_tephra_fall.nc"

# Satellite data file
if [ "$i" = "0" ]; then
  procfile="../Data/Kasatochi_SatCloudLoad_09hrs.nc"
  modis_file="2008-08-08_13-40 (9.1 hours)"
fi
if [ "$i" = "1" ]; then
  procfile="../Data/Kasatochi_SatCloudLoad_20hrs.nc"
  modis_file="2008-08-09_00-20 (20.2 hours)"
fi         
if [ "$i" = "2" ]; then
  procfile="../Data/Kasatochi_SatCloudLoad_32hrs.nc"
  modis_file="2008-08-09_12-45 (32.3 hours)"
fi         
if [ "$i" = "3" ]; then
  procfile="../Data/Kasatochi_SatCloudLoad_43hrs.nc"
  modis_file="2008-08-09_23-25 (43.3 hours)"
fi         
if [ "$i" = "4" ]; then
  procfile="../Data/Kasatochi_SatCloudLoad_55hrs.nc"
  modis_file="2008-08-10_11-49 (55.2 hours)"
fi         
if [ "$i" = "5" ]; then
  procfile="../Data/Kasatochi_SatCloudLoad_65hrs.nc"
  modis_file="2008-08-10_22-29 (64.6 hours)"
fi         
if [ "$i" = "6" ]; then
  procfile="../Data/Kasatochi_SatCloudLoad_78hrs.nc"
  modis_file="2008-08-11_10-53 (78.2 hours)"
fi         
if [ "$i" = "7" ]; then
  procfile="../Data/Kasatochi_SatCloudLoad_89hrs.nc"
  modis_file="2008-08-11_21-33 (89.5 hours)"
fi

BASE1="-Bg10/g5 -P"

# Create Base Map
gmt pscoast $AREA $PROJ $BASE1 $DETAIL $COAST -S100/149/237 -K  > temp.ps
# Plot satellite cloud load
gmt grdimage "${procfile}?Ash_Mass" -Q $AREA $PROJ $BASE1 -C${CPT} -K -O >> temp.ps

# Contour Ash3d output
echo "0.10 C" > ac0.1.lev
echo "1.00 C" > ac1.0.lev
echo "2.00 C" > ac2.0.lev
echo "3.00 C" > ac3.0.lev
gmt grdconvert "${infile}?cloud_load[$i]" temp.grd
gmt grdcontour temp.grd $AREA $PROJ $BASE1 -Cac0.1.lev  -A- -W1,0/0/0 -K -O >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE1 -Cac1.0.lev  -A- -W1,0/0/255 -K -O >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE1 -Cac2.0.lev  -A- -W1,255/0/255 -K -O >> temp.ps
gmt grdcontour temp.grd $AREA $PROJ $BASE1 -Cac3.0.lev  -A- -W1,255/0/0 -K -O >> temp.ps

# Plot legend
LEGLOC="-Dx0.1i/3.4i/4.0i/1.7i/BL"
gmt pslegend $AREA $PROJ $BASE -G255 $LEGLOC -K -O << EOF >> temp.ps
C black
H 14 1 Kasatochi ${modis_file}
D 1p
H 12 1 Ash3d
N 2
S 0.1i - 0.15i black   0.5p,black   0.3i 0.1 g/m^2
S 0.1i - 0.15i blue    0.5p,blue    0.3i 1.0 g/m^2
S 0.1i - 0.15i magenta 0.5p,magenta 0.3i 2.0 g/m^2
S 0.1i - 0.15i red     0.5p,red     0.3i 3.0 g/m^2
N 1
D 1p
H 12 1 Modis ash retrieval
EOF
gmt psscale -Dx1.25i/3.85i/2i/0.15ih -F+gwhite -C${CPT} -B1f1/:"g/m^2": -O -K >> temp.ps
echo $vln $vlt '1.0' | gmt psxy $AREA $PROJ -St0.1i -Gmagenta -Wthinnest -O >> temp.ps

# Save map
ps2pdf temp.ps
mv temp.pdf Kasatochi_cloudload_$i.pdf

# Clean up
rm temp.grd temp.ps ac*lev gmt.history gmt.conf cload.cpt

