Ash3d
==========

Ash3d is a 3-D Eulerian model to predict airborne volcanic ash concentration
and tephra deposition during volcanic eruptions. This model simulates downwind advection,
turbulent diffusion, and settling of ash injected into the atmosphere by a volcanic eruption
column. Ash advection is calculated using time-varying pre-existing wind data and a robust,
high-order, finite-volume method. The model can use either a spherical coordinate system
or a variety of projected coordinate systems, independent of the projection of the numerical
weather prediction (NWP) data used. Volcanic ash is specified with an arbitrary number of
grain sizes, which affects the fall velocity, distribution and duration of transport. Above
the source volcano, the vertical mass distribution with elevation can be specified using
a variety of source models, including a point source, line source, user-specified profile,
or with a Suzuki distribution for a given plume height, eruptive volume, and eruption
duration. Multiple eruptions separated in time may be included in a single simulation.
For larger events, an umbrella source can be used to account for the radial spreading
of the cloud.

The software is written in Fortran 2008 and is designed for a Linux system, although we
have had no trouble building the software on MacOS or Microsoft Windows.

For details on usage, please see the [User's Guide](doc/UsersGuide.md)
and look through the example programs.

## Building Ash3d

Ash3d is written in Fortran 2008 and is designed for a Linux system, although we
have had no trouble building the software on MacOS or Microsoft Windows.

### Prerequisite Software
This software relies on several auxiliary libraries for managing time calculations, projections,
and for interfacing with the variety of types of atmospheric data:  

- [HoursSince](https://code.usgs.gov/vsc/ash3d/volcano-ash3d-hourssince), [github mirror](https://github.com/DOI-USGS/volcano-ash3d-hourssince)
- [projection](https://code.usgs.gov/vsc/ash3d/volcano-ash3d-projection), [github mirror](https://github.com/DOI-USGS/volcano-ash3d-projection)
- [MetReader](https://code.usgs.gov/vsc/ash3d/volcano-ash3d-metreader), [github mirror](https://github.com/DOI-USGS/volcano-ash3d-metreader)

These libraries are currently available at the locations given above.
Installation instructions are given in the repositories for each of these libraries.  
The default installation
location for these libraries is `/opt/USGS`.  This could be changed to suit
your system by editing the makefile variable `INSTALLDIR` to another location,
but should be consistent since the Ash3d makefile expects
a single location. 

Additionally, several packages will significantly increase functionality and should be
installed, if possible:  

1. lapack:  This enables diffusion to be calculated via implicit algorithms, thereby significantly
            increasing the allowed time step.
2. netcdf:  Many NWP data products are provided in netcdf (typically v.4) format.  This library
            is needed to read these files. Moreover, Ash3d uses netcdf as the default format for
            storing output data.
3. ecCodes: This is used to directly read grib2 data, the format of choice for most forecast data.
4. zip    : This is needed to compress the kml output files into kmz files, as well as bundle
            shapefile output.
5. gnuplot: This is the default graphics package for creating plots directly from Ash3d.
6. dislin : This is an alternate graphics package for creating plots from Ash3d (preferred for Windows),
   available from `https://www.dislin.de`
7. plplot : Another alternate graphics package for creating plots from Ash3d (also works on Windows).

All of these packages (except dislin) are available on Red Hat and Ubuntu systems and can be installed
using the standard distribution software installer (yum/dnf for RedHat systems or apt for
Ubunto). For some of these packages and for some distributions, you might need to enable
extra repositories, such as epel and powertools/CRB (for Red Hat systems)

On RedHat-based systems, these can be installed with:  
`sudo yum install lapack lapack-devel blas blas-devel`  
`sudo yum install netcdf netcdf-devel netcdf-fortran netcdf-fortran-devel`  
`sudo yum install eccodes eccodes-devel`  
`sudo yum install zip`  
`sudo yum install gnuplot`  
`sudo yum install plplot`  

On Ubuntu-based systems:  
`sudo apt install liblapack64-3 liblapack64-dev libblas64-3 libblas64-dev`  
`sudo apt install libnetcdf-c++4 libnetcdf-dev libnetcdff-dev`  
`sudo apt libeccodes0 libeccodes-dev`  
`sudo apt install zip`  
`sudo apt install gnuplot`  
`sudo apt install plplot`  

Note that these particular package names are the distribution packages for the latest versions of RedHat and Ubuntu
at the time of writing. There may be slight variations to these names for older or newer
systems. If you prefer managing libraries either within a conda environment or via modules, you will likely have
to edit both `makefile` and `make_gfortran` (or whichever include file you need for the compiler
you are using), in order for the compiler to find the needed libraries and include files.
If you use modules, you will need to make sure that netcdf package you load (as well as eccodes, if desired)
is built with the compiler you plan to use.

### Compiling Ash3d with the default settings
Once the necessary USGS libraries and optional distribution packages are installed, Ash3d can
be built. The makefile in `volcano-ash3d/src` can be edited to suit your system.

To compile, verify that the makefile is consistent with the install directory and options
used in the installation of the preliminary software. Then simply type:  
`make`  
Alternatively, `make all` will also compile the `tools`, including `Ash3d_ASCII_check` which
compares two ASCII output files for run verification as well as `Ash3d_PostProc` which
can be used for post-processing Ash3d output. To test the installation, type:  
`make check`  
This will run all the tests in `tests/test_[1-4]` and run `Ash3d_ASCII_check` on the
resulting ASCII output files, reporting PASS or FAIL. These tests are described in
`tests/readme.txt` and are designed to be quick runs with output compared with expected
output. The test cases in 
`tests/test_04` require that the NCEP 50-year reanalysis data files for the year 1980
be stored in `/data/WindFiles/NCEP/1980`. Assuming MetReader was installed in the default
location, these reanayalsis files can be downloaded by running the following script
installed when MetReader was installed.  
`/opt/USGS/bin/autorun_scripts/get_NCEP_50YearReanalysis.sh 1980`  
This script expects that the directories `/data/WindFiles/NCEP` and
`/data/WindFiles/NCEP/dbuffer` exist and that the user
running the script has write permissions. If windfiles are stored elsewhere, please
make sure that the environment variable WINDROOT is set in your .bash\_profile or .bashrc
files.

To install the software, edit the `INSTALLDIR` variable of the makefile (the
default is `/opt/USGS`) and type:  
`make install`  
You will need to have write permission in `${INSTALLDIR}` or install as root.
This will install the following in `${INSTALLDIR}/bin/`:

Validation tests are found in the `examples` directory. These test cases compare
Ash3d results with both data and with past published results. Each of these test
cases require NCEP reanalysis data for the years of each eruption, which can be
downloaded as described above. The test cases are:

1. Mt. Spurr, Aug. 18, 1992  
This case models the fallout using a total grainsize distribution and fallout
measurements reported in [McGimsey et al, 2001](https://doi.org/10.3133/ofr01370).

2. Kasatochi, Aug. 8, 2008  
This case models the drifting ash cloud over several days and compares modeling
cloud load with satellite retrievals.

3. Mt. St. Helens, May 18, 1980  
This case does not model the full event, but just compares model output with a
published result from [Mastin et al, 2016](https://doi.org/10.5194/acp-16-9399-2016).

4. Mazama  
This case also does not model a particular event, but is intended to be compared
with published model results for ash fallout. The source term is an umbrella cloud.
Results are compared with those from [Buckland et al, 2022](https://doi.org/10.1007/s00445-022-01593-1).

5. Kelud, Feb. 13, 2014  
This case also uses an umbrella source, but just for tracking the ash cloud
(source is `umbrella_air`). Model results of the cloud development and advection
are compared with cloud outlines from satellite observations presented in
[Mastin and Van Eaton, 2020](https://doi.org/10.3390/atmos11101038).

The validation tests above demonstrate the utility of Ash3d in modeling volcanic
ash cloud transport and deposition for a variety of erpution scenarios, compared
with several types of observations. Verification tests, which demonstrate that
Ash3d solves the equations described in [Schwaiger et al.](https://doi.org/10.1029/2011JB008968),
can be found in the repository
[volcano-ash3d-optionalmodules](https://code.usgs.gov/vsc/ash3d/volcano-ash3d-optionalmodules).
This repository is structured as an example of how to construct an optional module
for user-created features and to build with the core Ash3d software. The example used in
this repository is a Testcases module that is composed of six sets of idealized tests
the compare Ash3d model results against analytic results for simplified cases of horizontal
advection, vertical advection, rigid rotation, diffusion, shear rotation and using the method
of manufactured solutions. Each test is set up to show the rate of convergence.

### Compiling Ash3d with the customized settings
To build Ash3d, we currently use a user-edited makefile. All the main
variables to edit are in the top block of the makefile up to the lines:  
`###############################################################################`  
`#####  END OF USER SPECIFIED FLAGS  ###########################################`  
`###############################################################################`  

The following are the variables available to edit:  

- `SYSTEM      = [gfortran] or ifort`  
              This controls which compiler to use and which libraries to link.
              If you use a different compiler than these two, you will have to
              edit the variables in the `SYSTEM` blocks to suit your system.  
- `RUN         = DEBUG, PROF, [OPT], or OMPOPT`  
              This variable sets which set of compiler flags will be used for
              the executable, either with debugging, profiling, optimization,
              or with `OMP` enabled.  
- `OS          = [LINUX], MACOS, WINDOWS`  
              This variable is used to indicate how Ash3d should build file
              paths when reading or writing files that require full paths.  
- `USGSROOT    = [/opt/USGS]`  
              This is the location of the USGS libraries and include files
              needed from `volcano-ash3d-hourssince`, `volcano-ash3d-projection`,
              and `volcano-ash3d-metreader`.  
- `ASH3DCCSRC  = [./]`  
              Location of the Ash3d core code source files. This is normally
              just the `cwd` unless you are building an optional module and
              need to link to the code Ash3d files.  
- `INSTALLDIR  = [/opt/USGS/Ash3d]`  
              Location of the installation directory.  
- `USENETCDF   = [T] or F`  
              Whether or not to build netCDF capabilities. If you only have
              netCDF v.3 available, then you will need to edit `ncFPPFLAG`
              to include `-DNC3`.  
- `USEGRIB     = [T] or F`  
              Whether or not to build GRIB2 capabilities via eccodes.  
- `USEPOINTERS = T or [F]`  
              The default is for arrays to be allocatable, but if you are
              building Ash3d to be called from C++ codes, such as forestclaw,
              then some arrays need to be defined as pointers. This does not
              work with older versions of gfortran.  
- `USEEXTDATA  = T or [F]`  
              Ash3d can be built with lists of Airports and global volcanoes
              as data variables set at compile-time.  Some low-memory systems
              might not be able to compile with these data variables (e.g. on
              Raspberry Pi systems).  Alternatively, these lists can be read
              at run-time from files installed on the system.  
- `FASTFPPFLAG = [] -DFAST_DT -DFAST_SUBGRID`  
              These are additional pre-processor flags that can be included
              to speed up the computation.  
           -  `-DFAST_DT` limits DT evaluations to just on the steps of the
               windfiles with linear interpolation of dt between steps. This
               cannot be used with sources that effect the winds (umbrella).  
          -  `-DFAST_SUBGRID` adjusts the computational domain to be just the
               min/max in x,y,z of the region where ash concentration exceeds
               some threshold.  
- `USEZIP      = T or [F]`  
              Zip is used to convert kml to kml and to bundle linked graphics,
              such as time-series ash accumulation plots, together with the kml
              file. Set this to F if zip is not installed on this system.
- `USEPLPLOT   = T or [F]`  
              If Ash3d_PostProc is to be built, this variable indicates if
              the plplot library is available on the system for plotting maps
              and/or vertical profiles.  
- `USEDISLIN   = T or [F]`  
              If Ash3d_PostProc is to be built, this variable indicates if
              the dislin library is available on the system for plotting maps
              and/or vertical profiles.  
- `USEGNUPLOT  = [T] or F`  
              Gnuplot is used by default to create time-series ash accumulation plots.
              If it is not installed, set this value to F.
- `USEGMT      = [F]`  
              If Ash3d_PostProc is to be built, this variable indicates if
              the GMT library is available on the system for plotting maps
              and/or vertical profiles. Currently, this is not enabled.
- `MFILE       = [makefile]`  
              The name of the makefile.  
- `LIMITER     = LIM_NONE, LIM_LAXWEN, LIM_BW, LIM_FROMM, LIM_MINMOD, [LIM_SUPERBEE], LIM_MC`  
              This variable specifies the limiter to use in the advection
              routines.  
- `DIFFMETH    = EXPLDIFF, [CRANKNIC]`  
              This variable specifies the diffusion scheme to be used.  
- `PII         = ON, [OFF]`  
              This variable allows information about a run to be included in
              the output netCDF file. This is set to OFF by default as the
              full path of the run directory, the hostname and the username
              might be considered sensitive information, but this information
              may be desired for logging of runs.

If you have a working executable, you can always see what settings were used to build it
by typing:  
`Ash3d -help make`  

### Ash3d Post-processing tool

Ash3d can write output products directly at run-time in a number of different formats, including
ERSI ASCII, binary, kml/kmz. Alternatively, it can also write a netcdf file which contains
all the information of these output products along with 3-D variables for transient
concentrations. If this netcdf file is written, the post-processing tool `Ash3d_PostProc`,
can be used to write ASCII, binary, kml, shapefiles, or png images of output variables.

To compile this tool, make sure the makefile variables `USEPLPLOT` and `USEDISLIN`
match your system and that the paths of the include directories and libraries are correct
for your system.
Then type:  
`make postproc`  

Usage
-----

Please see the [user's guide](doc/UsersGuide.md) for more information on using this software. A simplified
user's guide can be found in `volcano-ash3d/doc/UsersGuide.md`. A guide for using the
USGS web-based interface to Ash3d can be found at
[Mastin et al, 2013](https://doi.org/10.3133/ofr20131122).
Instructions on usage of the `Ash3d_PostProc` tool can be found [here](doc/PostProc.md)

References
----------

1. Buckland, et al, Modelling the transport and deposition of ash following a
magnitude 7 eruption: the distal Mazama tephra, Bulletin of Volcanology, v84,87,
2022,[doi:10.1007/s00445-022-01593-1](https://doi.org/10.1007/s00445-022-01593-1).

2. Mastin, L.G., M.J. Randall, H.F. Schwaiger, and R.P. Denlinger, 2013,
User’s Guide and Reference to Ash3d—A Three-Dimensional Model for Eulerian
Atmospheric Tephra Transport and Deposition, Open-File Report 2013-1122
[doi:10.3133/ofr20131122](https://doi.org/10.3133/ofr20131122).

3. Mastin, et al., 2016, Adjusting particle-size distribution to account for
aggregation in tephra-deposit model forecasts, Atmospheric Chemistry and
Physics, v16 n14, [doi:10.5194/acp-16-9399-2016](https://doi.org/10.5194/acp-16-9399-2016).

4. Mastin and Van Eaton, Comparing simulation of umbrella-cloud growth and ash
transport with observations from Pinatubo, Kelud, and Calbuco Volcanoes,
Atmosphere, v11,1038, 2020, [doi:10.3390/atmos11101038](https://doi.org/10.3390/atmos11101038).

5. McGimsey, R.G, Neal, C.A., and Riley, C.M., 2001, Areal distribution, thickness,
mass, volume, and grain size of tephra-fall deposits from the 1992 eruptions of
Crater Peak vent, Mt. Spurr Volcano, Alaska, U.S. Geological Survey Open-File
Report 01-370, [doi:10.3133/ofr01370](https://doi.org/10.3133/ofr01370).

6. Schwaiger, H.F., R.P. Denlinger, and L.G. Mastin, 2012, Ash3d, a finite-volume,
conservative numerical model for ash transport and tephra deposition,
Journal of Geophysical Research, 117, B04204,
[doi:10.1029/2011JB008968](https://doi.org/10.1029/2011JB008968).


Authors
-------

Hans F. Schwaiger <hschwaiger@usgs.gov>  
Larry G. Mastin <lgmastin@usgs.gov>  
Roger P. Denlinger <rdenlinger@usgs.gov>

Citation
--------

If you would like to reference Ash3d in your paper, please cite the software as follows.  

Schwaiger, H.F. et al. (2024) Ash3d (Version 1.0.0), U.S. Geological Survey Software Release,
[doi:10.5066/P1SJWAKZ](https://doi.org/10.5066/P1SJWAKZ).  

or you can cite the journal article describing the software given in reference 6 above.

License and Disclaimer
----------------------

[LICENSE](LICENSE.md): This project is in the public domain.  

[DISCLAIMER](DISCLAIMER.md): This software is preliminary or provisional and is subject to revision.

