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

# Kelud coordinates
vlt=-14.4
vln=99.0

lonw=100.0
lone=118.0
lats=-14.0
latn=0.0
DETAIL="-Dl"
BASE="-Bg5/g5 -P"
PROJ="-JM${vln}/${vlt}/7.0i"
AREA="-R${lonw}/${lats}/${lone}/${latn}r"
COAST="-G220/220/220 -W"

# Set up some default values
gmt gmtset PROJ_ELLIPSOID Sphere
gmt gmtset PAPER_MEDIA=Custom_720x510

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
for i in `seq 0 7`;
do

  gmt grdconvert "$infile?area" zero.grd
  gmt grdmath 0.0 zero.grd MUL = zero.grd
  # We need to convert the NaN's to zero to get the lowest contour
  gmt grdconvert "$infile?cloud_load[$i]" temp.grd
  gmt grdmath temp.grd zero.grd AND = temp.grd
  
  #******************************************************************************
  
  # Create Base Map
  gmt pscoast $AREA $PROJ $BASE $DETAIL $COAST -S100/149/237 -K  > temp.ps

  # Filled contour of Ash3d output
  #gmt grdimage "${infile}?cloud_load[$i]" -Q $AREA $PROJ $BASE -CGMT_hot.cpt -K -O >> temp.ps
  #gmt psscale -Dx1.25i/3.85i/2i/0.15ih -F+gwhite -CGMT_hot.cpt -B1f1/:"g/m^2": -O -K >> temp.ps

  # Contour Ash3d output
  echo "100.0 C" > ac100.0.lev
  gmt grdcontour temp.grd $AREA $PROJ $BASE -Cac100.0.lev  -A- -W2,255/0/0     -K -O >> temp.ps
  
  # Plot satellite cloud outline
  if [ $i -eq 0 ]; then
    cloud="../Data/m40_1619.csv"
    sat_file="2014-02-16:19 (0.11 hours)"
  elif [ $i -eq 1 ]; then
    cloud="../Data/m40_1659.csv"
    sat_file="2014-02-16:59 (0.77 hours)"
  elif [ $i -eq 2 ]; then
    cloud="../Data/m40_1709.csv"
    sat_file="2014-02-17:09 (0.94 hours)"
  elif [ $i -eq 3 ]; then
    cloud="../Data/m40_1739.csv"
    sat_file="2014-02-17:39 (1.44 hours)"
  elif [ $i -eq 4 ]; then
    cloud="../Data/m40_1809.csv"
    sat_file="2014-02-18:09 (1.94 hours)"
  elif [ $i -eq 5 ]; then
    cloud="../Data/m40_1839.csv"
    sat_file="2014-02-18:39 (2.44 hours)"
  elif [ $i -eq 6 ]; then
    cloud="../Data/m40_1909.csv"
    sat_file="2014-02-19:09 (2.94 hours)"
  elif [ $i -eq 7 ]; then
    cloud="../Data/m40_1939.csv"
    sat_file="2014-02-19:39 (3.44 hours)"
  fi
  
  cat ${cloud} | awk -F "," '{print $2, $1}' | gmt psxy $AREA $PROJ $BASE -W2,0/255/0 -V -K -O >> temp.ps
  
  # Plot legend
  LEGLOC="-Dx0.1i/4.0i/2.8i/1.0i/BL"
gmt pslegend $AREA $PROJ $BASE -G255 $LEGLOC -K -O << EOF >> temp.ps
C black
H 14 1 Kelud Umbrella Cloud
H 12 1 ${sat_file}
D 1p
S 0.1i - 0.15i red     0.5p,red     0.3i Ash3d 500 g/m^2
S 0.1i - 0.15i green   0.5p,green   0.3i Sat.Cloud
EOF
  
  # Last gmt command is to plot the volcano and close out the ps file
  echo $vln $vlt '1.0' | gmt psxy $AREA $PROJ -St0.1i -Gmagenta -Wthinnest -O >> temp.ps
 
  # Save map
  #ps2epsi temp.ps temp.eps
  #convert temp.eps Kelud_CloudOutline_$i.png
  #epstopdf temp.eps Kelud_CloudOutline_$i.pdf
  ps2pdf temp.ps Kelud_CloudOutline_$i.pdf
 
  # Clean up
  rm -f gmt.history gmt.conf zero.grd temp.* ac*lev 

done
