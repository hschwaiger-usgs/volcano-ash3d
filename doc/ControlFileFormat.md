Ash3d Control File Fomat
========================

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
is assumed to be in the Smithsonian database and default values for
Plume Height, Duration, Mass Flux Rate, Volume, and mass fraction of
fines are loaded. These can be over-written by entering non-negative
values in the appropriate locations in this input file.  

##### Projection specification  
Line 2 of this block identifies the projection used and the form of
the input coordinates and is of the following format:  
   `latlonflag`, `projflag`, followed by a variable list of projection parameters
`projflag` describes the projection used for the Ash3d run. The windfiles used to
populate the computational grid can have a different projection.  
For a particular `projflag`, additional values are read defining the projection.  
`latlonflag` = 0 if computational grid is projected  
`latlonflag` = 1 if computational grid is lat/lon (all subsequent projection parameters ignored.)  
1. projflag   = 1: polar stereographic projection  
 - lambda0 -- longitude of projection point  
 - phi0    -- latitude of projection point  
 - k0      -- scale factor at projection point  
 - radius  -- earth radius for spherical earth  
e.g. for NAM 104,198, 216: 0 1 -105.0 90.0 0.933 6371.229  
2. projflag   = 2: Albers Equal Area (not yet implemented)  
3. projflag   = 3: UTM (not yet implemented)  
4. projflag   = 4: Lambert conformal conic  
 - lambda0 -- longitude of origin  
 - phi0 -- latitude of origin  
 - phi1 -- latitude of secant1  
 - phi2 -- latitude of secant2  
 - radius -- earth radius for a spherical earth  
e.g. for NAM 212: 0 4 265.0 25.0 25.0 25.0 6371.22  
5. projflag   = 5: Mercator  
 - lambda0 -- longitude of origin  
 - phi0 -- latitude of origin  
 - radius -- earth radius for a spherical earth  
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

##### Eruption category and number
Line 8 is the the diffusivity (m2/s) followed by the eruption specifier. The
eruption specifier can be a real number, in which case it is assumed to be the
positive constant specifying the Suzuki distribution. Alternatively, it can be  
 `umbrella`: Suzuki (const. = 12) with radial spreading of the plume  
 `umbrella_air `: Suzuki (const. = 12) with radial spreading of the plume scaled to 5% of vol.  
 `point`: all mass inserted in cell containing PlmH  
 `linear`: mass uniformly distributed from z-vent to PlmH  
 `profile`: mass distributed with a user-specified vertical profile  
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
In the following line, each line represents one eruptive pulse.
Parameters are (1-4) start time (yyyy mm dd h.hh (UT)); (5) duration (hrs);
(6) plume height; (7) erupted volume (km3 DRE)
If neruptions=1 and the year is 0, then the model run in forecast mode where mm dd h.hh are
interpreted as the time after the start of the windfile. In this case, duration, plume
height and erupted volume are replaced with ESP if the values are negative.
This applies to source types: `suzuki`, `point`, `line`, `umbrella` and `umbrella_air`.
For profile sources, an additional two values are read: `dz` and `nz`  
`2010 04 14   0.00   1.0     18.0  0.16 1.0 18`  
`0.01 0.02 0.03 0.03 0.04 0.04 0.05 0.06 0.06 0.070 0.08 0.08 0.09 0.09 0.09 0.08 0.06 0.02`  
`******************* BLOCK 2 ***************************************************`  
`2010 4 14  9.0  3.0 7.4 8.46E-004`  
`2010 4 14 12.0  3.0 8.4 1.38E-003`  
`2010 4 14 15.0  3.0 6.8 6.11E-004`  
`2010 4 14 18.0  3.0 5.6 2.89E-004`  
`2010 4 14 21.0  3.0 5.1 2.01E-004`  
`2010 4 15  0.0  3.0 5.3 2.33E-004`  
`2010 4 15  3.0  3.0 5.2 2.17E-004`  
`2010 4 15  6.0  3.0 5.3 2.33E-004`  
`2010 4 15  9.0  3.0 5.7 3.09E-004`  
`*******************************************************************************`  

#### BLOCK 3: Wind Parameters
Ash3d will read from either a single 1-D wind sounding, or gridded, time-                            
dependent 3-D wind data, depending on the value of the parameter iwind.                              
`iwind =`  
1. read from a 1-D wind sounding                                                         
2. read from 3D gridded ASCII files                                                      
3/4. read directly from a single or multiple NetCDF files.                               
5. read directly from multiple multi-timestep NetCDF files.                              

The parameter iwindformat specifies the format of the wind files, as follows:  
`iwindformat =`  
0. User-defined via template                                                         
1. User-specified ASCII files                                                        
2. Global radiosonde data                                                            
3. NARR 221 Reanalysis (32 km)                                                       
4. NAM Regional North America 221 Forecast (32 km)                                   
5. NAM 216 Regional Alaska Forecast (45 km)                                          
6. NAM 104 Northern Hemisphere Forecast (90 km)                                      
7. NAM 212 40km Cont. US Forecast (40 km)                                            
8. NAM 218 12km Cont. US Forecast (12 km)                                            
9. NAM 227 Cont. US Forecast (5.08 km)                                               
10. NAM 242 11km Regional Alaska Forecast (11.25 km)                                  
11. NAM 196 Regional Hawaii Forecast (2.5 km)                                         
12. NAM 198 Regional Alaska Forecast (5.953 km)                                       
13. NAM 91 Regional Alaska Forecast (2.976 km)                                        
14. NAM Regional Cont. US Forecast (3.0 km)                                           
20. GFS 0.5 degree files Forecast                                                     
21. GFS 1.0 degree files Forecast                                                     
22. GFS 0.25 degree files Forecast                                                    
23. NCEP DOE Reanalysis 2.5 degree                                                    
24. NASA MERRA-2 Reanalysis                                                           
25. NCEP1 2.5 global Reanalysis (1948-pres)                                           
26. JRA-55 Reanalysis                                                                 
27. NOAA-CIRES II 2-deg global Reanalysis (1870-2010)                                 
28. ECMWF ERA-Interim Reanalysis                                                      
29. ECMWA ERA-5 Reanalysis                                                            
30. ECMWA ERA-20C Reanalysis                                                          
32. Air Force Weather Agency                                                          
33. CCSM 3.0 Community Atmospheric Model                                              
40. NASA GEOS-5 Cp                                                                    
41. NASA GEOS-5 Np                                                                    
50. Weather Research and Forecast (WRF) output                                        
  
`igrid` is the NCEP grid ID. If a NWP product is used, or the number of stations of
sonde data, if iwind = 1.  This is optional and defaults to that associated with 
`iwindformat`.  
`idata` is a flag for data type (1=ASCII, 2=netcdf, 3=grib).
This is also optional and defaults to 2 for netcdf.  
  
Many plumes extend higher than the maximum height of mesoscale models.  
Ash3d handles this as determined by the parameter iHeightHandler, as follows:  
`iHeightHandler =`  
1. stop the program if the plume height exceeds mesoscale height  
2. wind velocity at levels above the highest node
equal that of the highest node. Temperatures in the
upper nodes do not change between 11 and 20 km; above
20 km they increase by 2 C/km, as in the U.S. Standard
atmosphere. A warning is written to the log file.  

Simulation time in hours is the maximal length of the simulation.  
Ash3d can end the simulation early if desired, once 99% of the ash has deposited.
The last line of this block is the number of windfiles listed in block 5 below. If
iwind=5 and one of the NWP products is used that require a special file structure,
then nWindFiles should be set to 1 and only the root folder of the windfiles listed.  
`******************* BLOCK 3 ***************************************************`  
`4  20               # iwind, iwindformat, [igrid, idata]`  
`2                   # iHeightHandler`  
`60                  # Simulation time in hours`  
`no                  # stop computation when 99% of erupted mass has deposited?`  
`16                  # nWindFiles, number of gridded wind files (used if iwind>1)`  
`*******************************************************************************`  

#### BLOCK 4: Output Options
The list below allows users to specify the output options
All but the final deposit file can be written out at specified
times using the following parameters:
Line 15 asks for 3d output (yes/no) followed by an optional output format code;  
  1 = (default) output all the normal 2d products to the output file as well as the 3d concentrations  
  2 = only output the 2d products  
nWriteTimes  = if >0, number of times output are to be written. The following
line contains nWriteTimes numbers specifying the times of output
if =-1, it specifies that the following line gives a constant time
interval in hours between write times.
WriteTimes   = Hours between output (if nWritetimes=-1), or
Times (hours since start of first eruption) for each output
(if nWriteTimes >1).  
`******************* BLOCK 4 ***************************************************`  
`no      # Write out ESRI ASCII file of final deposit thickness?`  
`yes     # Write out        KML file of final deposit thickness?`  
`no      # Write out ESRI ASCII deposit files at specified times?`  
`no      # Write out        KML deposit files at specified times?`  
`no      # Write out ESRI ASCII files of ash-cloud concentration?`  
`no      # Write out        KML files of ash-cloud concentration?`  
`no      # Write out ESRI ASCII files of ash-cloud height?`  
`no      # Write out        KML files of ash-cloud height?`  
`yes     # Write out ESRI ASCII files of ash-cloud load (T/km2) at specified times?`  
`yes     # Write out        KML files of ash-cloud load (T/km2) at specified times?`  
`yes     # Write out ESRI ASCII file of deposit arrival times?`  
`yes     # Write out        KML file of deposit arrival times?`  
`yes     # write out ESRI ASCII file of cloud arrival times?`  
`yes     # Write out        KML file of cloud arrival times?`  
`yes 1   # Write out 3-D ash concentration at specified times? / [output code: 1=2d+concen,2=2d only]`  
`netcdf  # format of ash concentration files     (ascii, binary, or netcdf)`  
`-1      # nWriteTimes`  
`1       # WriteTimes (hours since eruption start)`  
`*******************************************************************************`  

#### BLOCK 5: Input Wind Files
The following block of data contains names of wind files.
If we are reading from a 1-D wind sounding (i.e. iwind=1) then there should
be only one wind file.
If we are reading gridded data there should be iWinNum wind files, each having
the format volcano_name_yyyymmddhh_FHhh.win
For iwindformat=25 (continuous wind data from 1948-pres), just give
the directory with the windfiles (e.g. Wind_nc).  
`******************* BLOCK 5 ***************************************************`  
`Wind_nc/gfs/gfs.2010041400/2010041400.f000.nc`  
`Wind_nc/gfs/gfs.2010041400/2010041400.f003.nc`  
`Wind_nc/gfs/gfs.2010041400/2010041400.f006.nc`  
`Wind_nc/gfs/gfs.2010041400/2010041400.f009.nc`  
`Wind_nc/gfs/gfs.2010041400/2010041400.f012.nc`  
`Wind_nc/gfs/gfs.2010041400/2010041400.f015.nc`  
`Wind_nc/gfs/gfs.2010041400/2010041400.f018.nc`  
`Wind_nc/gfs/gfs.2010041400/2010041400.f021.nc`  
`Wind_nc/gfs/gfs.2010041400/2010041400.f024.nc`  
`Wind_nc/gfs/gfs.2010041400/2010041400.f027.nc`  
`Wind_nc/gfs/gfs.2010041400/2010041400.f030.nc`  
`Wind_nc/gfs/gfs.2010041400/2010041400.f033.nc`  
`Wind_nc/gfs/gfs.2010041400/2010041400.f036.nc`  
`Wind_nc/gfs/gfs.2010041400/2010041400.f039.nc`  
`Wind_nc/gfs/gfs.2010041400/2010041400.f042.nc`  
`Wind_nc/gfs/gfs.2010041400/2010041400.f045.nc`  
`*******************************************************************************`  


#### BLOCK 6: Airport File
The following lines allow the user to specify whether times of ash arrival
at airports and other locations will be written out, and which file
to read for a list of airport locations.
Each line in the airport location file should contain the
airport latitude, longitude, projected x and y coordinates,
and airport name. If you are using a projected grid,
THE X AND Y MUST BE IN THE SAME PROJECTION as the computational grid.
Alternatively, if coordinates can be projected via libprojection
by typing "yes" to the last parameter.  
`******************* BLOCK 6 ***************************************************`  
`no                            # Write out ash arrival times at airports to ASCII FILE?`  
`no                            # Write out grain-size distribution to ASCII airports file?`  
`no                            # Write out ash arrival times to kml file?`  
`GlobalAirports.txt            # Name of file containing airport locations`  
`no                            # Have libprojection calculate projected coordinates?`  
`*******************************************************************************`  

#### BLOCK 7: Grain-size bins, Settling Velocity
The first line must contain the number of settling velocity groups, but
can optionally also include a flag for the fall velocity model to be used and a flag
specifying the shape parameterization used.  
Fall velocity model `FV_ID =`  
1. Wilson and Huang (default)  
2. Wilson and Huang + Cunningham slip  
3. Wilson and Huang + Mod by Pfeiffer Et al.  
4. Ganser   
5. Ganser + Cunningham slip  
6. Stokes flow for spherical particles + slip  
If no fall model is specified, FV_ID = 1, by default
The grain size bins can be enters with 2, 3, or 4 parameters.
If TWO are given, they are read as:   FallVel (in m/s), mass fraction
If THREE are given, they are read as: diameter (mm), mass fraction, density (kg/m3)
If FOUR are given, they are read as:  diameter (mm), mass fraction, density (kg/m3), Shape F
The shape factor is given as in Wilson and Huang: F=(b+c)/(2a), but converted
to sphericity (assuming b=c) for the Ganser model.
If a shape factor is not given, a default value of F=0.4 is used.
If FIVE are given, they are read as:  diameter (mm), mass fraction, density (kg/m3), Shape F, G
where G is an additional Ganser shape factor equal to c/b.  

Shape Parameterization `Shape_ID =`  
1. Shape factor from Wilson and Huang for column 4 (default)  
2. column 4 holds sphericity  


If the last grain size bin has a negative diameter, then the remaining mass fraction
will be distributed over the previous bins via a log-normal distribution in phi.
The last bin would be interpreted as:
`diam (neg value) , phi_mean, phi_stddev`  

`*******************************************************************************`  
`15                                 # Number of grain-size bins. FV_ID not given; defaults to 1`  
`2.000      0.0208  2003.   0.44    # grain size (mm), mass fraction, density (kg/m3), F=0.44`  
`1.414      0.0084  2350.   0.44`  
`1.000      0.0141  2005.   0.44`  
`0.707      0.0214  2248.   0.44`  
`0.500      0.0459  2624.   0.44`  
`0.354      0.0723  2644.   0.44`  
`0.250      0.0532  2639.   0.44`  
`0.177      0.0219  2690.   0.44`  
`0.125      0.0165  2691.   0.44`  
`0.088      0.0115  2730.   0.44`  
`0.2176     0.0714  600.    1.00    # Note that these are bigger again, low density and round`  
`0.2031     0.1428  600.    1.00    # These bottom five bins represent an aggregate distribution`  
`0.1895     0.2856  600.    1.00`  
`0.1768     0.1428  600.    1.00`  
`0.1649     0.0714  600.    1.00`  
`*******************************************************************************`  

#### BLOCK 8: Vertical Profiles
The first line below gives the number of locations (nlocs) where vertical
profiles are to be written.  That is followed by nlocs lines, each of which
contain the location, in the same coordinate as the computational grid.
Optionally, a site name can be provided in after the location.  
`******************* BLOCK 8 ***************************************************`  
`4                             # number of locations for vertical profiles (nlocs)`  
`12.4  51.4  Leipzig           # x,y (or lon/lat) [Site name]`  
`11.3  48.2  Munich            # Munich (Maisach)`  
`11.0  47.4  Schneefernerhaus  # Schneefernerhaus (Zugspitze)`  
`11.0  47.8  Hohenpeissenberg  # Hohenpeissenberg`  
`*******************************************************************************`  

#### BLOCK 9 (Optional): netCDF Annotations
This last block is optional.
The output file name can be give, but will default to 3d_tephra_fall.nc if absent
The title and comment lines are passed through to the netcdf header of the
output file.  
`******************* BLOCK 9 ***************************************************`  
`3d_tephra_fall.nc             # Name of output file`  
`Eyjafjallajokull              # Title of simulation`  
`no comment                    # Comment`  
`*******************************************************************************`  


#### BLOCK 10: Optional Modules
Optional Modules are identified by the text string at the top of the block
`OPTMOD=[module name]`.  These are not required to be in the control file, but
can be included to invoke features of the program other than the default
behavior. These additional input blocks allow a means for controlling user-provided
features such as non-standard source terms or physical processes.
There will need to be a custom block reader in the module to read this section
section of the input file.  There are two optional modules built in to the
main Ash3d code: RESTETPARAMS and TOPO.  


Below is the built-in example for resetting parameters.
Shown below are all the parameters available to be reset (along with the default
value). Only the parameters to be reset need to be listed.  
`*******************************************************************************`  
`OPTMOD=RESETPARAMS`  
`MagmaDensity         = 3500.0`  
`DepositDensity       = 1300.0`  
`LAM_GS_THRESH        = 250.0`  
`AIRBORNE_THRESH      = 1.0e-3`  
`GRAV                 = 9.81`  
`RAD_EARTH            = 6371.229`  
`CFL                  = 0.80`  
`DT_MIN               = 1.0e-5`  
`DT_MAX               = 1.0`  
`ZPADDING             = 1.3`  
`DEPO_THRESH          = 1.0e-1`  
`DEPRATE_THRESH       = 1.0e-2`  
`CLOUDCON_THRESH      = 2.0e-1`  
`CLOUDCON_GRID_THRESH = 2.0e-1`  
`CLOUDLOAD_THRESH     = 1.0e-2`  
`THICKNESS_THRESH     = 1.0e-3`  
`StopValue_FracAshDep = 0.99`  
`DBZ_THRESH           = -2.0e+1`  
`VelMod_umb           = 1`  
`lambda_umb           = 0.2`  
`N_BV_umb             = 0.02`  
`k_entrainment_umb    = 0.1`  
`SuzK_umb             = 12.0`  
`useMoistureVars      = F`  
`useVz_rhoG           = T`  
`useWindVars          = 0`  
`useOutprodVars       = 1`  
`useRestartVars       = 0`  
`cdf_institution      = USGS`  
`cdf_run_class        = Analysis`  
`cdf_url              = https://vsc-ash.wr.usgs.gov/ash3d-gui`  
`*******************************************************************************`  


Topography can be included in Ash3d by including the following optional block.  
`*******************************************************************************`  
`OPTMOD=TOPO`  
`yes 2                         # use topography?; z-mod (0=none,1=shift,2=sigma)`  
`1 1.0                         # Topofile format, smoothing radius`  
`GEBCO_08.nc                   # Topofile name`  
`*******************************************************************************`  

Line 1 indicates whether or not to use topography followed by the integer flag
describing how topography will modify the vertical grid.  
0. = no vertical modification; z-grid remains 0-> top throughout the domain  
1. = shifted; s = z-z_surf; computational grid is uniformly shifted upward everywhere by topography  
2. = sigma-altitude; s=(z-z_surf)/(z_top-z_surf); topography has decaying influence with height  
Line 2 indicates the topography data format followed by the smoothing radius in km.
Topofile format must be one of:  
1. Gridded lon/lat (netcdf)  
ETOPO : https://www.ncei.noaa.gov/products/etopo-global-relief-model  
GEBCO : https://www.gebco.net/  
2. Gridded Binary  
NOAA Globe (1-km/30 arcsec) https://www.ngdc.noaa.gov/mgg/topo/globe.html  
GTOPO30 (1-km/30 arcsec)  
3. ESRI ASCII  



