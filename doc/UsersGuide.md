Ash3d User Guide
=================

Please see [README.md](../README.md) for instructions on installing Ash3d.

Once the Ash3d software is built and installed, the following executables are
copied to `$(INSTALLDIR)/bin/`:  
 `Ash3d`                   : main Ash3d program  
 `Ash3d_PostProc`          : Post-processing program for netcdf output.  
 `tools/Ash3d_ASCII_check` : program for calculating L2 norm of two ASCII output files.  

## Running Ash3d on the command line

Ash3d requires a control file to run. Example control files are given in this
repository in the examples folder. The control file format is described
[here](ControlFileFormat.md). To run Ash3d with a control file, simply
pass the control file name as a command-line argument to the executable:  
`./Ash3d Spurr_081992_ESP.inp`  
Alternatively, Ash3d can be run interactively, where the user will be prompted
for the control file, followed by a prompt asking if the run will be restarting
a previous run from a saved time step, if
desired. If a restart case is requested, the control file and the netcdf file
must be consistent.  

Ash3d can also provide basic help via the command:  
`Ash3d -h`  
or with more specific help requests, such as:  
`Ash3d -h make`  
`Ash3d -h run`  
`Ash3d -h input`  
`Ash3d -h postproc`  

### Run-time environment variables

There are several environment variables that can be used to modify the run
conditions of Ash3d, if desired.  

#### `ASH3DVERB`: 
The variable specifies the output verbosity level. The default
is level of 3 for output, but this variable can be used to override this
setting at run-time. Valid values are integers from 1-10. The lower the
level, the more verbose the output. For each level, information at that
level as well as all with higher levels are printed.  
1. = debug2     : Additional debugging information only written to stdout  
2. = debug1     : Debugging information only written to stdout  
3. = log        : Time step information (limit for writing to logfile)  
4. = info       : Additional information on run set up and shutdown  
5. = statistics : Details on health of run (timing, mass conservation)  
6. = production : Major program flow info only  
7. = essential  : Only start up and shutdown messages  
8. = error      : No logging to stdout, only stderr (and logfile)  
9. = silent     : No logging to stdout,stderr. Logfile written as normal  
10. = dark       : No logging to stdout,stderr or logfile  
To run Ash3d with no output to the terminal, but with a logfile, for example:  
`ASH3DVERB=9 ./Ash3d control.inp`  

#### `ASH3DHOME`: 
This is the path to the install folder.
If the executable is compiled such that it needs access to external files in
the installation directory, the path can be changed at run-time, if needed.  
`ASH3DHOME=/opt/USGS/Ash3d ./Ash3d control.inp`  
This path is also used with the post-processing tool.

#### `ASH3DCFL`: 
The Courant-Friedrichs-Lewy (CFL) factor can be changed at
run-time from the default of 0.8. The value must be greater than 0.0. The
advection routines become unstable if CFL>1.0.  
`ASH3DCFL=0.5 ./Ash3d control.inp`  

#### `OMP_NUM_THREADS`: 
This is not an Ash3d-specific environment variable,
but can be used if the executable was compiled with OpenMP. This specifies
the number of threads can be specified at run-time.  
`OMP_NUM_THREADS=4 ./Ash3d_opt control.in`  

### Run-directory

Ash3d reads files from and writes files to the run directory. Typically, the
control file governing the run will reside in this directory. This control
file will include paths to the windfiles that describe the transient atmospheric
conditions for the run. Also, an airport/point-of-interest file may be provided
with a path specified in the control file.  

For example, GFS forecast data might be stored in `/data/WindFiles/gfs/latest/`.
The run directory could have a link to this directory in the `cwd`.  
`ln -s /data/WindFiles/gfs/latest/ gfs`  
with the corresponding block of the control file having:  
`gfs/2015072100.f00.nc`  
`gfs/2015072100.f03.nc`  
`gfs/2015072100.f06.nc`  
`gfs/2015072100.f09.nc`  
`gfs/2015072100.f12.nc`  
`gfs/2015072100.f15.nc`  
`gfs/2015072100.f18.nc`  

Once the Ash3d program starts, the following files are opened:  
1. `Ash3d.lst`: output log file  
2. 3d output file (e.g. `3d_tephra_fall.nc`)  
3. `vprofile[##].txt` if vertical profile output is requested  
4. Any time-series kml file requested including:  
 - `CloudConcentration.kml`  
 - `CloudLoad.kml`  
 - `CloudHeight.kml`  
 - `deposit_thickness_mm.kml`  
 - `deposit_thickness_inches.kml`  
5. Any 2d ASCII files requested are written as the time-step is available:  
 - `CloudConcentration_006.00hrs.dat`  
 - `CloudHeight_006.00hrs.dat`  
 - `CloudLoad_006.00hrs.dat`  
 - `DepositFile_006.00hrs.dat`  
6. `progress.txt`: text file containing the fraction complete.  
7. `ash_arrivaltimes_airports.txt` if line 1 of block 6 is 'yes'  
8. `ash_arrivaltimes_airports.kml` if line 3 of block 6 is 'yes'  
9. `depTS_[####].png`: time-series plot of deposit accumulation at airports/POI.  

After completion, the KML/KMZ files can be viewed using your favorite KML-viewer,
such as GoogleEarth. The ESRI ASCII files can be imported into a GIS tool. The
output netcdf file can be viewed or processed with a variety of third-party tools.
Simple post-processing of the netcdf file (with limited mapping capibilities) can
be accomplished with the post-processing tool, `Ash3d_PostProc`. Instructions on
using this tool can be found [here](PostProc.md), along with other post-processing
recomendations.





