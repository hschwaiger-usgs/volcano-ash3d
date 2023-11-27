Ash3d Post-Processing
=====================

There are many output products that Ash3d can write at run-time, including
ESRI ASCII or kml/kmz files for each of final deposit thickness, transient
deposit thickness, ash-cloud concentration, ash-cloud load, ash-cloud height,
deposit arrival time, cloud arrival time, vertical profiles of ash concentration,
and a kml/kmz file showing the transient ash accumulation at airports or point-of-
interest.  Additionally, 3d transient data can be written in ASCII, binary
or NetCDF formats. NetCDF is the preferred output format as the full content of the
Ash3d run is written to this file, including all the information necessary to
recreate the input file, the state of all environment variables controlling the
run, as well as additional information about the run (date, time, user,
path, hostname, etc.).
This NetCDF file also contains all the data needed to recreate any of the
output products Ash3d could produce at run-time.

`Ash3d_PostProc` is a tool that can be used to post-process the Ash3d output
from the NetCDF output file.  Any of the standard output products (ESRI ASCII,
kml/kmz, binanry) can be recreated using this tool along with png images of
2d data in mapview (or contour plots for vertical prodiles) or shapefiles of 2d
data.  To produce image files, `Ash3d_PostProc` can use a variety of geographic
plotting software.  Both `dislin` and `plplot` are available as libraries and can
be linked to `Ash3d_PostProc` at compilation time.  `dislin` includes the feature
where contour data can be accessed by `Ash3d_PostProc` and used to generate
shapefiles.  This is the preferred means of generating shapefiles on a Microsoft
Windows system.  `Ash3d_PostProc` can also generate maps using `gnuplot` and
Generic Mapping Tools (`GMT`) via the writing and execution of temporary scripts.

## Running Ash3d_PostProc with command-line options

`Ash3d_PostProc` is designed to be a quick tool for converting the NetCDF
output of an Ash3d run to format that can be visualized. This tool can be
run interactively by just running the executable with no command-line
arguments.  If one argument is provided, it is interpreted to be the name
of a control file (described below).
Minimal instructions
for running this tool are available by running the program with a `-h` as
the only argument
`Ash3d_PostProc -h`
` Dislin   T`
` Plplot   T`
` Gnuplot  T`
` GMT      T`
`                                                                                `
` Ash3d post-processing tool: Ash3d_PostProc                                     `
`                                                                                `
`Usage: Ash3d_PostProc control_file [t_index]                                    `
`           or                                                                   `
`       Ash3d_PostProc infile output_product format                              `
`  where: infile   = the netcdf file written by Ash3d                            `
`   output_product = 1 full concentration array                                  `
`                    2 deposit granularity                                       `
`                    3 deposit thickness (mm time-series)                        `
`                    4 deposit thickness (inches time-series)                    `
`                    5 deposit thickness (mm final)                              `
`                    6 deposit thickness (inches final)                          `
`                    7 ashfall arrival time (hours)                              `
`                    8 ashfall arrival at airports/POI (mm)                      `
`                    9 ash-cloud concentration (mg/m3)                           `
`                   10 ash-cloud height (km)                                     `
`                   11 ash-cloud bottom (km)                                     `
`                   12 ash-cloud load (T/km2 or )                                `
`                   13 ash-cloud radar reflectivity (dBz)                        `
`                   14 ash-cloud arrival time (hours)                            `
`                   15 topography                                                `
`                   16 profile plots                                             `
`           format = 1 ASCII/ArcGIS                                              `
`                    2 KML/KMZ                                                   `
`                    3 image/png                                                 `
`                    4 binary                                                    `
`                    5 shape file                                                `
`         [t_index] = index of time slice to plot; -1 for final (optional)       `

For example, to generate a contour plot of the final deposit thickness in mm from
the output NetCDF file, enter:

`Ash3d_PostProc 3d_tephra_fall.nc 5 3`

The preferred graphics package is
set at the time of compilation depending on the system (Linux, Windows, MacOS)
and the availability of the libraries.  These can alway be overriden at
run-time with the environmet variable `ASH3DPLOT`, where: `1=dislin`, `2=plplot`,
`3=gnuplot`, and `4=GMT`.

`ASH3DPLOT=3 Ash3d_PostProc 3d_tephra_fall.nc 5 3`

To generate a shapefile of the cloud load at step 3 of the output file, enter:
`Ash3d_PostProc 3d_tephra_fall.nc 12 5 3`

## Running Ash3d_PostProc with a control file

Using a control file give much more flexibility and allows plotting maps
from 2d ASCII or binary files, or to customize the contour levels and colors.

Below is an example of a control file that can be used to convert a 2d ESRI ASCII
output file of the deposit thickness to a map image with custom contour levels.

`DepositFile_____final.dat       data_filename`
`1                       input format code (1=ascii, 2=binary, 3=netcdf)`
`test.inp                Ash3d_control_file_name   (skipped if datafile is netcdf)`
`5                       invar [outvar]    input variable and output variable, if different`
`2 1                     ndims Only needed for format 1 or 2, LatLonFlag`
`38 20                   nx ny [nz] Also only needed for format 1 or 2`
`1.0 1.0                 dx dy [dz] Needed if not a part of the data file`
`225.0 33.0              srtx srty         Start x and y`
`3                       output format     (1=ascii, 2=KML 3=image, 4=binary, 5=shapefile)`
`4                       plot_pref         (1=dislin, 2=plplot, 3=gnuplot, 4=GMT)`
`-1                      time_step         Only needed if input file is multi-timestep (eg netcdf) (-1 for final, -2 for all)`
`0                       Filled contour flag (1 for filled, 0 for lines)`
`1 5                     custom contour flag (1 for true, 0 for false), number of contours`
`1.0 3.0 10.0 50.0 100.0 lev(ncont)  : contour levels`
`100 100 100 100 255     R(ncont)    : Red channel of RGB`
`100 150 200 150   0     G(ncont)    : Green channel of RGB`
`200 150 100  50   0     B(ncont)    : Blue channel of RGB`

Line 1 contains the data file name, in this case `DepositFile_____final.dat`, but it could be
the output NetCDF file or any 2/3d binary or ASCII output file.
Line 2 specifies the format code for ASCII, binary or netcdf.
Line 3 is the name of the Ash3d control file for this run.  This is needed to know the geometry
of the simulation, the time steps, volcano name, eruption source parameters, and
other information needed for annotating the maps or other output products.  If the
input data file is the Ash3d NetCDF output file, then this Ash3d input control file
is ignored as all the information required is stored in the NetCDF file.
Line 4 is the variable code for the input file, which is assumed to be the output
variable code unless a second code for the output variable is provided. For example, if
the data file is the Ash3d NetCDF data file, then entering `5` in line 4 will result
in a plot of the final deposit thickness.  If the input data file is a 3d binary
file of the airborne ash concentration, but cloud load is needed for the output, line
4 would contain `1 12`.
Line 5 contains the number of dimensions of the data file (2 or 3) and a code indicating if
the data are projected or in lon/lat coordinates (0 or 1).  `Ash3d_PostProc` currently
only can plot lon/lat data.
Line 6 gives the number of nodes in the x,y, and possibly z directions for the data file.
Line 7 gives the grid spacing for all the coordinate directions of the file and line 8
gives the start coordinates.
Line 8 is the output format code where 1=ascii, 2=KML 3=image, 4=binary, 5=shapefile.
Line 9 is the plotting package preference code: 1=dislin, 2=plplot, 3=gnuplot, 4=GMT.
Line 10 is the time step to plot.  -1 can be entered to plot the last time step of the
file.
Line 11 is a flag to indicate if the plot should have contour lines (0) or filled contours (1).
This option is not yet fully functional and only contour lines are available.
Line 12 starts with a flag allowing custom contour levels (1) to override the default
levels for that variable.  If this flag is set, a second value is read, giving the
number of custom levels (`nlev`).
If custom contours are requested, then four additional lines are read in the control file.
Line 13 contains the `nlev` floating point values for the contour levels in the
default units of the output variable.
Lines 14, 15, and 16 are the `nlev` integer values of the RGB componets of the colors (0-255).


