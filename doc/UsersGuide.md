Ash3d User Guide
=================

Please see README.md for instructions on installing Ash3d.

Once the Ash3d software is built and installed, the following executables are
copied to `$(INSTALLDIR)/bin/`:  
 `Ash3d`                   : main Ash3d program  
 `Ash3d_PostProc`          : Post-processing program for netcdf output.  
 `tools/Ash3d_ASCII_check` : program for calculating L2 norm of two ASCII output files.  

## Running Ash3d on the command line

Ash3d requires a control file to run.  Example control files are given in this
repository in the examples folder.  To run Ash3d with a control file, simply
pass the control file name as a command-line argument to the executable:  
`./Ash3d Spurr_081992_ESP.inp`  
Alternatively, Ash3d can be run interactively, where the user will be prompted
for the control file, followed by a prompt for restarting a previous run, if
desired.  If a restart case is requested, the control file and the netcdf file
must be consistent.  

### Run-time environment variables

There are several environment variables that can be used to modify the run
conditions of Ash3d, if desired.  
`ASH3DVERB = 1-10`  
   This overrides the default level of 3 for output.  
    1 = debug2     : Additional debugging information only written to stdout  
    2 = debug1     : Debugging information only written to stdout  
    3 = log        : Time step information (limit for writing to logfile)  
    4 = info       : Additional information on run set up and shutdown  
    5 = statistics : Details on health of run (timing, mass conservation)  
    6 = production : Major program flow info only  
    7 = essential  : Only start up and shutdown messages  
    8 = error      : No logging to stdout, only stderr (and logfile)  
    9 = silent     : No logging to stdout,stderr. Logfile written as normal  
   10 = dark       : No logging to stdout,stderr or logfile  
To run Ash3d with no output to the terminal, but with a logfile, for example:  
`ASH3DVERB=9 ./Ash3d control.inp`

`ASH3DHOME` = path to installation location

If the executable is compiled such that it needs access to external files in
the installation directory, the path can be changed at run-time, if needed.  
`ASH3DHOME=/opt/USGS/Ash3d ./Ash3d control.inp`

`ASH3DCFL = 0.0-1.0`

The Courant-Friedrichs-Lewy (CFL) factor can be changed at run-time from the
default of 0.8.  
`ASH3DCFL=0.5 ./Ash3d control.inp`

If the executable was compiled with OpenMP, then the number of threads can be
specified at run-time.  
`OMP_NUM_THREADS=4 ./Ash3d_opt control.in`

### Run-directory

Ash3d reads files from and writes files to the run directory.  Typically, the
control file governing the run will reside in this directory.  This control
file will include paths to the windfiles that describe the transient atmospheric
conditions for the run.  Also, an airport/point-of-interest file may be provided
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

### Control file

The control file is composed of at least 8 blocks in the following order:  
1. GRID INFO  
2. ERUPTION PARAMETERS  
3. WIND PARAMETERS  
4. OUTPUT OPTIONS  
5. INPUT WIND FILES  
6. AIRPORT FILE  
7. GRAIN-SIZE BINS, SETTLING VELOCITY  
8. VERTICAL PROFILES  
9. (Optional): NETCDF ANNOTATIONS  
10. OPTIONAL MODULES  
The blocks are marked with a line of `*` (or at least a line beginning with `*`) at
the start and end of the block.  Lines starting with `#` between blocks are treated
as comments.

#### BLOCK 1: Grid Info  

Line 1 of this block identifies the volcano by name.
If the volcano name begins with either 0 or 1, then the volcano
is assumed to be in the Smithonian database and default values for
Plume Height, Duration, Mass Flux Rate, Volume, and mass fraction of
fines are loaded.  These can be over-written by entering non-negative
values in the appropriate locations in this input file.  

##### Projection specification  
Line 2 of this block identifies the projection used and the form of
the input coordinates and is of the following format:  
   `latlonflag`, `projflag`,  followed by a variable list of projection parameters
`projflag` describes the projection used for the Ash3d run. The windfiles used to
populate the computational grid can have a different projection.  
For a particular `projflag`, additional values are read defining the projection.  
`latlonflag` = 0 if computational grid is projected or  
`latlonflag` = 1 if computational grid is lat/lon (all subsequent projection parameters ignored.)  
projflag   = 1 -- polar stereographic projection  
 lambda0 -- longitude of projection point  
 phi0    -- latitude of projection point  
 k0      -- scale factor at projection point  
 radius  -- earth radius for spherical earth  
e.g. for NAM 104,198, 216: 0 1 -105.0 90.0 0.933 6371.229  
 = 2 -- Alberts Equal Area ( not yet implemented)  
 = 3 -- UTM ( not yet implemented)  
 = 4 -- Lambert conformal conic  
          lambda0 -- longitude of origin  
             phi0 -- latitude of origin  
             phi1 -- latitude of secant1  
             phi2 -- latitude of secant2  
           radius -- earth radius for a spherical earth  
    e.g. for NAM 212: 0 4 265.0 25.0 25.0 25.0 6371.22  
              = 5 -- Mercator  
          lambda0 -- longitude of origin  
             phi0 -- latitude of origin  
           radius -- earth radius for a spherical earth  
    e.g. for NAM 196: 0 5 198.475 20.0 6371.229   

##### Vent and grid specification  
On line 3, the vent coordinates can optionally include a third value for elevation in km.
If the vent elevation is not given, 0 is used if topography is turned off.  

Line 4 is the width and height of the computational grid in km (if projected) or degrees.  
Line 5 is the vent x,y (or lon, lat) coordinates.  
Line 6, DX and DY resolution in km or degrees (for projected or lon/lat grid, respectively).  
Line 7, DZ can be given as a real number, indicating the vertical spacing in km.
Alternatively, it can be given as `dz_plin` (piece-wise linear), `dz_clog` (constant-
logarithmic), or `dz_cust` (custom specification).
If `dz_plin`, then a second line is read containing:
  number of line segments (N) followed by the steps and step-size of each segment  
  e.g. 4 6 0.25 5 0.5 5 1.0 10 2.0  
        This corresponds to 4 line segments with 6 cells of 0.25, then 5 cells of 0.5,
        5 cells of 1.0, and finally 10 cells of 2.0  
If `dz_clog`, then a second line is read containing:
  maximum z and number of steps of constant dlogz  
  e.g. 30.0 30  
        This corresponds to 30 steps from 0-30km with constant log-spacing  
If `dz_cust`, then a second line is read containing:  
  the number of dz values to read (ndz), followed by dz(1:ndz)  
  e.g. 20 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 5.5  
        This corresponds to 10 steps of 0.5, 9 steps of 1.5, followed by 1 step of 5.5.  

##### Eruption catagory and number
Line 8 is the the diffusivity (m2/s) followed by the eruption specifier.  The
eruption specifier can be a real number, in which case it is assumed to be the
positive constant specifying the Suzuki distribution.  Alternatively, it can be  
 `umbrella`: Suzuki (const. = 12) with radial spreading of the plume  
 `umbrella_air `: Suzuki (const. = 12) with radial spreading of the plume scaled to 5% of vol.  
 `point`: all mass inserted in cell containing PlmH  
 `linear`: mass uniformly distributed from z-vent to PlmH  
Line 9 : number of pulses to be read in BLOCK 2  

Example  
`******************* BLOCK 1 ***************************************************`  
`Eyjafjallajokull               # Volcano name (character*30)`  
`1 1 0.0 90.0 0.933 6367.470    # Proj flags and params; first term (LLflag) is 1, so all else ignored`  
`-25.0    45.0                  # x, y of LL corner of grid (km, or deg. if latlongflag=1)`  
`55.0     25.0                  # grid width and height (km, or deg. if latlonflag=1)`  
`-19.62   63.63                 # vent location         (km, or deg. if latlonflag=1)`  
`3.0      3.0                   # DX, DY of grid cells  (km, or deg. if latlonflag=1)`  
`1.0                            # DZ of grid cells      (always km)`  
`0.     4.                      # diffusion coefficient (m2/s), Suzuki constant`  
`9                              # neruptions, number of eruptions or pulses`  
`*******************************************************************************`  



#### BLOCK 2: Eruption Parameters
#### BLOCK 3: Wind Parameters
#### BLOCK 4: Output Options
#### BLOCK 5: Input Wind Files
#### BLOCK 6: Airpot File
#### BLOCK 7: Grain-size bins, Settling Velocity
#### BLOCK 8: Vertical Profiles
#### BLOCK 9 (Optional): netCDF Annotations
#### BLOCK 10: Optional Modules




