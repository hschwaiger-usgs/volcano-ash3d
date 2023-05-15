###############################################################################
#  Makefile for Ash3d
#
#    User-specified flags are in this top block
#
###############################################################################

#      This file is a component of the volcanic ash transport and dispersion model Ash3d,
#      written at the U.S. Geological Survey by Hans F. Schwaiger (hschwaiger@usgs.gov),
#      Larry G. Mastin (lgmastin@usgs.gov), and Roger P. Denlinger (roger@usgs.gov).

#      The model and its source code are products of the U.S. Federal Government and therefore
#      bear no copyright.  They may be copied, redistributed and freely incorporated 
#      into derivative products.  However as a matter of scientific courtesy we ask that
#      you credit the authors and cite published documentation of this model (below) when
#      publishing or distributing derivative products.

#      Schwaiger, H.F., Denlinger, R.P., and Mastin, L.G., 2012, Ash3d, a finite-
#         volume, conservative numerical model for ash transport and tephra deposition,
#         Journal of Geophysical Research, 117, B04204, doi:10.1029/2011JB008968. 

#      We make no guarantees, expressed or implied, as to the usefulness of the software
#      and its documentation for any purpose.  We assume no responsibility to provide
#      technical support to users of this software.

#      Sequence of commands:
#      "make"  compiles the Ash3d executable
#      "make all" builds the executables and copies to bin
#      "make install" copies the contents of branches/core_code/bin to the install location
#                        e.g. /opt/USGS/Ash3d
#
#  SYSTEM specifies which compiler to use
#    Current available options are:
#      gfortran , ifort
#    This variable cannot be left blank
#
SYSTEM = gfortran
#SYSTEM = ifort
#
#  RUN specifies which collection of compilation flags that should be run
#    Current available options are:
#      DEBUG : includes debugging info and issues warnings
#      PROF  : includes profiling flags with some optimization
#      OPT   : includes optimizations flags for fastest runtime
#      DEBUGOMP : includes debugging info and issues warnings and OpenMP directives
#      OMPOPT   : includes optimizations flags for fastest runtime and OpenMP directives
#                 To run, enter: env OMP_NUM_THREADS=4 Ash3d input_file.inp
#    This variable cannot be left blank
RUN = DEBUG
#RUN = PROF
#RUN = OPT
#RUN = OMPOPT

#    This variable cannot be left blank
OS = LINUX
#OS = MACOS
#OS = WINDOWS

# This is the location of the USGS libraries and include files
# This might be different from the INSTALLDIR path if you want to link to special MetReader
# projection and HoursSince packages.
USGSROOT=/opt/USGS

# This is the location of the Ash3d core code source files that we will build this
# research branch with.
ASH3dCCSRC=./

INSTALLDIR=/opt/USGS/Ash3d

# DATA FORMATS
#  For each data format you want to include in the library, set the corresponding
#  variable below to 'T'.  Set to 'F' any you do not want compiled or any unavailable
USENETCDF = T
USEGRIB   = T

#  VERB specifies the verbosity level for debugging purposes
#    Current available options are:
#      VERB=-1 : no output to either stdout or logfile
#      VERB=0  : no output to stdout; forces stop if interactive input needed
#      VERB=1  : standard level of logging to stdout and logfile
#      VERB=2  : normal output to stdout (VERB=1) with additional write statements
#      VERB=3  : 
#       -      : same as VERB=1
VERB = 0

# MEMORY
# If you need pointer arrays instead of allocatable arrays, set this to 'T'
USEPOINTERS = F

# EXTERNAL DATA FILES
#
USEEXTDATA = T

# SPEED-UP ALGORITHMS
#  FAST_DT
#  FAST_SUBGRID
FASTFPPFLAG = 
#FASTFPPFLAG = -DFAST_DT
#FASTFPPFLAG = -DFAST_SUBGRID
#FASTFPPFLAG = -DFAST_DT -DFAST_SUBGRID

# PLOTTING LIBRARIES
#  The default plotting utility for Ash3d_PostProc is gnuplot with plots generated
#  through system calls to that external program.  Alternatively, you can use either
#  the plplot or dislin plotting packages.  Currently, shapefile require contour to
#  be generated either with gnuplot or dislin.  Set the corresponding variable
#  below to 'T' if you want to build with these libraries.  Set to 'F' any you want
#  to use the default gnuplot.  
USEPLPLOT = T
USEDISLIN = T
USEGMT    = F

MFILE=makefile

###############################################################################
#  NOTE:  The following flags can either be set or commented out.  If they are
#         commented out, LIMITER defaults to LIM_SUPERBEE and DIFFMETH
#         defaults to CRANKNIC
#
#  LIMITER specifies which limiter to use
#    Current available options are:
#      LIM_NONE
#      LIM_LAXWEN
#      LIM_BW
#      LIM_FROMM
#      LIM_MINMOD
#      LIM_SUPERBEE
#      LIM_MC
LIMITER = LIM_SUPERBEE
#LIMITER = LIM_NONE
#
#  DIFFMETH specifies whether or not to use Crank-Nicolson
#
# Explicit calculation of diffusion term
#DIFFMETH=EXPLDIFF
# Implicit calculation of diffusion using Crank-Nicolson
#DIFFMETH=EXPLDIFF
DIFFMETH=CRANKNIC
#
OPTMODS_FPP= 
OPTMODS_OBJ=

# Ash3d can log some useful information about a run such as the username of 
# account that ran the case, then hostname of the system and run directory
# (or current working directory).  Although this is useful for record-keeping,
# this can be considered sensitive information.  To disable this logging,
# set PII=OFF
#PII=ON
PII=OFF

###############################################################################
#####  END OF USER SPECIFIED FLAGS  ###########################################
###############################################################################



###############################################################################
###############################################################################

OSFPPFLAG = -DLINUX
ifeq ($(OS), LINUX)
    OSFPPFLAG = -DLINUX
endif
ifeq ($(OS), MACOS)
    OSFPPFLAG = -DMACOS
endif
ifeq ($(OS), WINDOWS)
    OSFPPFLAG = -DWINDOWS
endif


LIMFPPFLAG = -DLIM_SUPERBEE
ifeq ($(LIMITER), LIM_NONE)
    LIMFPPFLAG = -DLIM_NONE
endif
ifeq ($(LIMITER), LIM_LAXWEN)
    LIMFPPFLAG = -DLIM_LAXWEN
endif
ifeq ($(LIMITER), LIM_BW)
    LIMFPPFLAG = -DLIM_BW
endif
ifeq ($(LIMITER), LIM_FROMM)
    LIMFPPFLAG = -DLIM_FROMM
endif
ifeq ($(LIMITER), LIM_MINMOD)
    LIMFPPFLAG = -DLIM_MINMOD
endif
ifeq ($(LIMITER), LIM_SUPERBEE)
    LIMFPPFLAG = -DLIM_SUPERBEE
endif
ifeq ($(LIMITER), LIM_MC)
    LIMFPPFLAG = -DLIM_MC
endif

DIFFFPPFLAG = -DCRANKNIC
ifeq ($(DIFFMETH), EXPLDIFF)
    DIFFFPPFLAG = -DEXPLDIFF
endif
ifeq ($(DIFFMETH), CRANKNIC)
    DIFFFPPFLAG = -DCRANKNIC
endif

VERBFPPFLAG = -DVERBOSE_L1
ifeq ($(VERB), 0)
    VERBFPPFLAG = -DVERBOSE_L0
endif
ifeq ($(VERB), 1)
    VERBFPPFLAG = -DVERBOSE_L1
endif
ifeq ($(VERB), 2)
    VERBFPPFLAG = -DVERBOSE_L2
endif
ifeq ($(VERB), 3)
    VERBFPPFLAG = -DVERBOSE_L3
endif

DATAFPPFLAG =
ifeq ($(USENETCDF), T)
 ncFPPFLAG = -DUSENETCDF
  # if you only have netcdf3, then add this flag to block out references to nc4 subroutines
 #ncFPPFLAG = -DUSENETCDF -DNC3
 ncLIBS = -lnetcdf -lnetcdff
endif
ifeq ($(USEGRIB), T)
 grbFPPFLAG = -DUSEGRIB
 # These are the libraries for grib_api
 #grbLIBS = -lgrib_api_f90 -lgrib_api
 # These are the libraries for ecCodes
 grbLIBS = -leccodes -leccodes_f90
endif

DATAFPPFLAG = $(ncFPPFLAG) $(grbFPPFLAG)
DATALIBS = $(ncLIBS) $(grbLIBS)

ifeq ($(USEPOINTERS), T)
 memFPPFLAG = -DUSEPOINTERS
endif

ifeq ($(USEEXTDATA), T)
 extdatFPPFLAG = -DUSEEXTDATA
endif

PLOTFPPFLAG =
ifeq ($(USEPLPLOT), T)
 plplotFPPFLAG = -DUSEPLPLOT
 #plplotLIBS = -lplplotf95cd -lplplotf95d
 #plplotINC  = -I/usr/include/plplot
 plplotLIBS = -L/usr/local/plplot/lib -lplplot -lplplotfortran
 plplotINC  = -I/usr/local/plplot/lib/fortran/modules/plplot
 plplotOBJ  = Ash3d_PostProc_plplot.o
endif
ifeq ($(USEDISLIN), T)
 dislinFPPFLAG = -DUSEDISLIN
 # Note that -ldislin if for single and -ldislin_d is double precision
 dislinLIBS = -L/usr/local/dislin -ldislin_d
 dislinINC  = -I/usr/local/dislin/gf
 dislinOBJ  = Ash3d_PostProc_dislin.o
endif
ifeq ($(USEGMT), T)
 gmtFPPFLAG = -DUSEGMT
 gmtLIBS = -L/usr/lib64 -lgmt -lpsl
 gmtINC  = -I/usr/include/gmt
 gmtOBJ  = gmt.o Ash3d_PostProc_GMT.o
 gmtfortloc = /home/hschwaiger/Programs/GIT/gmt-fortran
endif


# location of MetReader, HoursSince and projection
USGSLIBDIR = -L$(USGSROOT)/lib
USGSLIB = -lMetReader -lhourssince -lprojection
USGSINC = -I$(USGSROOT)/include

PLOTFPPFLAG = $(plplotFPPFLAG) $(dislinFPPFLAG) $(gmtFPPFLAG)
PLOTLIBS= $(plplotLIBS) $(dislinLIBS) $(gmtLIBS)
PLOTINC=$(plplotINC) $(dislinINC) $(gmtINC)

PIIFPP =
ifeq ($(PII), ON)
 PIIFPPFLAG = -DUSEPII
endif
###############################################################################
##########  GNU Fortran Compiler  #############################################
ifeq ($(SYSTEM), gfortran)
    FCHOME=/usr
    FC=/usr/bin/gfortran
    COMPINC = -I$(FCHOME)/include -I$(FCHOME)/lib64/gfortran/modules
    COMPLIBS = -L$(FCHOME)/lib64 -L$(FCHOME)/lib
    LIBS = $(COMPLIBS) $(USGSLIBDIR) $(USGSINC) $(COMPINC) -llapack -lblas $(USGSLIB) $(DATALIBS)
    # -lefence
# Debugging flags
ifeq ($(RUN), DEBUG)
    FFLAGS = -O0 -g3 -Wall -Wextra -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all -pedantic  -fbacktrace -Wunderflow -ffpe-trap=invalid,zero,overflow 
#-fdefault-real-8 -std=f2008 
    ASH3DEXEC=Ash3d_debug
endif
ifeq ($(RUN), DEBUGOMP)
    FFLAGS =  -g3 -pg -Wall -fbounds-check -pedantic -fimplicit-none -Wunderflow -Wuninitialized -Wmaybe-uninitialized -ffpe-trap=invalid,zero,overflow -fdefault-real-8 -fopenmp -lgomp
    ASH3DEXEC=Ash3d_debugomp
endif
# Profiling flags
ifeq ($(RUN), PROF)
    FFLAGS = -g -pg -w -fno-math-errno -funsafe-math-optimizations -fno-trapping-math -fno-signaling-nans -fcx-limited-range -fno-rounding-math -fdefault-real-8
    ASH3DEXEC=Ash3d_prof
endif
# Production run flags
ifeq ($(RUN), OPT)
    FFLAGS = -O3 -w -fno-math-errno -funsafe-math-optimizations -fno-trapping-math -fno-signaling-nans -fcx-limited-range -fno-rounding-math -fdefault-real-8
    ASH3DEXEC=Ash3d_opt
endif
ifeq ($(RUN), OMPOPT)
    FFLAGS = -O3 -w -fno-math-errno -funsafe-math-optimizations -fno-trapping-math -fno-signaling-nans -fcx-limited-range -fno-rounding-math -fdefault-real-8 -fopenmp -lgomp
    ASH3DEXEC=Ash3d_omp
endif
      # Preprocessing flags
    FPPFLAGS =  -x f95-cpp-input $(OSFPPFLAG) $(DATAFPPFLAG) $(VERBFPPFLAG) $(LIMFPPFLAG) $(DIFFFPPFLAG) $(memFPPFLAG) $(extdatFPPFLAG) $(FASTFPPFLAG) $(PLOTFPPFLAG) $(PIIFPPFLAG)
      # Extra flags
    #EXFLAGS = -xf95
    EXFLAGS =
endif
###############################################################################

OBJECTS = \
Ash3d_VariableModules.o  \
Ash3d_OS_info.o          \
Input_Data.o             \
Help.o                   \
VotW.o                   \
Airports.o               \
Calc_Mesh.o              \
Alloc_Arrays.o           \
Source.o                 \
Umbrella_Winds.o         \
Tephra.o                 \
Adjust_dt.o              \
TimeStepTotals.o         \
Atmosphere.o             \
MesoInterpolater.o       \
AdvectionHorz.o          \
AdvectionHorz_DCU.o      \
AdvectionVert_DCU.o      \
Diffusion.o              \
Set_BC.o                 \
Output_Vars.o            \
Output_Results.o         \
write_ASCII.o            \
write_BINARY.o           \
write_KML.o              \
ResetParams.o

PPOBJS = Ash3d_VariableModules.o        \
         Output_Vars.o Output_Results.o \
         Tephra.o                       \
         Source.o                       \
         Airports.o                     \
         Atmosphere.o                   \
         write_KML.o                    \
         write_ASCII.o                  \
         write_BINARY.o                 \
         write_Shapefile.o              \
         Ash3d_PostProc_cities.o        \
         Ash3d_OS_info.o

GRAFOBJ = Ash3d_PostProc_gnuplot.o $(plplotOBJ) $(dislinOBJ) $(gmtOBJ)

ifeq ($(USENETCDF), T)
  DATAOBJ = write_NETCDF.o
endif

Ash3d:  $(OBJECTS) $(DATAOBJ) Ash3d.F90 $(MFILE)
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) -c Ash3d.F90 
	$(FC) $(FFLAGS) $(EXFLAGS) $(OBJECTS) $(DATAOBJ) Ash3d.o $(LIBS) -o $(ASH3DEXEC)
	mkdir -p ../bin
	cp $(ASH3DEXEC) ../bin/Ash3d
	mv $(ASH3DEXEC) ../bin/$(ASH3DEXEC)


# Objects/Modules for initilization and model set-up
Ash3d_VariableModules.o: $(ASH3dCCSRC)/Ash3d_VariableModules.f90 $(MFILE)
	./get_version.sh
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/Ash3d_VariableModules.f90
Ash3d_OS_info.o: $(ASH3dCCSRC)/Ash3d_OS_info.F90 $(MFILE) $(ASH3dCCSRC)/Ash3d_VariableModules.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/Ash3d_OS_info.F90
Input_Data.o: $(ASH3dCCSRC)/Input_Data.f90 $(MFILE) Ash3d_VariableModules.o Tephra.o Help.o VotW.o Diffusion.o Airports.o $(DATAOBJ) Output_Vars.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/Input_Data.f90
ResetParams.o: $(ASH3dCCSRC)/ResetParams.f90 $(MFILE) Ash3d_VariableModules.o Tephra.o Output_Vars.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/ResetParams.f90
Help.o: $(ASH3dCCSRC)/Help.f90 $(MFILE) Ash3d_VariableModules.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) -c $(ASH3dCCSRC)/Help.f90
VotW.o: $(ASH3dCCSRC)/VotW.f90 $(MFILE) Ash3d_VariableModules.o Source.o
	echo "        Ash3dHome = '$(INSTALLDIR)'" > installpath.h
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/VotW.f90
Airports.o: $(ASH3dCCSRC)/Airports.f90 $(MFILE) Ash3d_VariableModules.o 
	echo "        Ash3dHome = '$(INSTALLDIR)'" > installpath.h
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/Airports.f90
Calc_Mesh.o: $(ASH3dCCSRC)/Calc_Mesh.f90 $(MFILE) Ash3d_VariableModules.o Tephra.o Atmosphere.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/Calc_Mesh.f90
Alloc_Arrays.o: $(ASH3dCCSRC)/Alloc_Arrays.f90 $(MFILE) Ash3d_VariableModules.o Output_Vars.o AdvectionHorz.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/Alloc_Arrays.f90
Source.o: $(ASH3dCCSRC)/Source.f90 $(MFILE) Ash3d_VariableModules.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/Source.f90
Umbrella_Winds.o: $(ASH3dCCSRC)/Umbrella_Winds.f90 $(MFILE) Source.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/Umbrella_Winds.f90
Tephra.o: $(ASH3dCCSRC)/Tephra.f90 $(MFILE) Atmosphere.o Ash3d_VariableModules.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/Tephra.f90

# Time
Adjust_dt.o: $(ASH3dCCSRC)/Adjust_dt.f90 $(MFILE) Ash3d_VariableModules.o Tephra.o Diffusion.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) -c $(ASH3dCCSRC)/Adjust_dt.f90
TimeStepTotals.o: $(ASH3dCCSRC)/TimeStepTotals.f90 $(MFILE) Ash3d_VariableModules.o Tephra.o Output_Vars.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) -c $(ASH3dCCSRC)/TimeStepTotals.f90

Atmosphere.o: $(ASH3dCCSRC)/Atmosphere.f90 $(MFILE) Ash3d_VariableModules.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/Atmosphere.f90
MesoInterpolater.o: $(ASH3dCCSRC)/MesoInterpolater.f90 $(MFILE) Atmosphere.o Tephra.o Ash3d_VariableModules.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/MesoInterpolater.f90

# PDE Solver (Advection, Diffusion, B.C.'s)
AdvectionHorz.o: $(ASH3dCCSRC)/AdvectionHorz.f90 $(MFILE) AdvectionHorz_DCU.o Ash3d_VariableModules.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) -c $(ASH3dCCSRC)/AdvectionHorz.f90
AdvectionHorz_DCU.o: $(ASH3dCCSRC)/AdvectionHorz_DCU.f90 $(MFILE) Ash3d_VariableModules.o Tephra.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) -c $(ASH3dCCSRC)/AdvectionHorz_DCU.f90
AdvectionVert_DCU.o: $(ASH3dCCSRC)/AdvectionVert_DCU.f90 $(MFILE) Ash3d_VariableModules.o Tephra.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) -c $(ASH3dCCSRC)/AdvectionVert_DCU.f90
Diffusion.o: $(ASH3dCCSRC)/Diffusion.f90 $(MFILE) Ash3d_VariableModules.o Tephra.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) -c $(ASH3dCCSRC)/Diffusion.f90
Set_BC.o: $(ASH3dCCSRC)/Set_BC.f90 $(MFILE) Ash3d_VariableModules.o Tephra.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) -c $(ASH3dCCSRC)/Set_BC.f90

# Output
Output_Vars.o: $(ASH3dCCSRC)/Output_Vars.f90 $(MFILE) Ash3d_VariableModules.o Tephra.o Airports.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/Output_Vars.f90
Output_Results.o: $(ASH3dCCSRC)/Output_Results.f90 $(MFILE) Ash3d_VariableModules.o Output_Vars.o Tephra.o write_ASCII.o write_KML.o write_NETCDF.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/Output_Results.f90
write_ASCII.o: $(ASH3dCCSRC)/write_ASCII.f90 $(MFILE) Ash3d_VariableModules.o Airports.o Output_Vars.o Tephra.o Source.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/write_ASCII.f90
write_BINARY.o: $(ASH3dCCSRC)/write_BINARY.f90 $(MFILE) Ash3d_VariableModules.o Airports.o Output_Vars.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/write_BINARY.f90
ifeq ($(USENETCDF), T)
write_NETCDF.o: $(ASH3dCCSRC)/write_NETCDF.f90 $(MFILE) Ash3d_VariableModules.o Output_Vars.o Tephra.o Source.o Airports.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/write_NETCDF.f90
endif
#ifeq ($(USEGRIB), T)
#write_GRIB.o: $(ASH3dCCSRC)/write_GRIB.f90 $(MFILE) Ash3d_VariableModules.o Output_Vars.o Tephra.o Source.o Airports.o
#	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/write_GRIB.f90
#endif
write_KML.o: $(ASH3dCCSRC)/write_KML.f90 $(MFILE) Ash3d_VariableModules.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/write_KML.f90

all: Ash3d tools


# Post-Processing tool
postproc: $(ASH3dCCSRC)/Ash3d_PostProc.f90 $(PPOBJS) $(GRAFOBJ) $(DATAOBJ) $(MFILE)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) $(PLOTINC) $(PLOTLIBS) -c $(ASH3dCCSRC)/Ash3d_PostProc.f90
	$(FC) $(FFLAGS) $(EXFLAGS) Ash3d_PostProc.o $(PPOBJS) $(GRAFOBJ) $(DATAOBJ) $(LIBS) $(PLOTINC) $(PLOTLIBS) -o ../bin/Ash3d_PostProc
# We always compile the gnuplot subroutine since these will just throw an error if gnuplot is not installed
Ash3d_PostProc_gnuplot.o: Ash3d_PostProc_gnuplot.f90 $(MFILE) Ash3d_PostProc_cities.o
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(PLOTINC) $(PLOTLIBS) -c $(ASH3dCCSRC)/Ash3d_PostProc_gnuplot.f90
ifeq ($(USEPLPLOT), T)
Ash3d_PostProc_plplot.o: Ash3d_PostProc_plplot.f90 $(MFILE) Ash3d_PostProc_cities.o
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(PLOTINC) $(PLOTLIBS) -c $(ASH3dCCSRC)/Ash3d_PostProc_plplot.f90
endif
ifeq ($(USEDISLIN), T)
Ash3d_PostProc_dislin.o: Ash3d_PostProc_dislin.f90 $(MFILE) Ash3d_PostProc_cities.o
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(PLOTINC) $(PLOTLIBS) -c $(ASH3dCCSRC)/Ash3d_PostProc_dislin.f90
endif
ifeq ($(USEGMT), T)
Ash3d_PostProc_GMT.o: Ash3d_PostProc_GMT.f90 $(MFILE) Ash3d_PostProc_cities.o gmt.o
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(PLOTINC) $(PLOTLIBS) -c $(ASH3dCCSRC)/Ash3d_PostProc_GMT.f90
gmt.o: $(gmtfortloc)/src/gmt.f90 $(MFILE)
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(PLOTINC) $(PLOTLIBS) -c $(gmtfortloc)/src/gmt.f90
endif
Ash3d_PostProc_cities.o: Ash3d_PostProc_cities.f90 $(MFILE) Ash3d_VariableModules.o
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(PLOTINC) $(PLOTLIBS) -c $(ASH3dCCSRC)/Ash3d_PostProc_cities.f90
write_Shapefile.o: write_Shapefile.f90 $(MFILE) Ash3d_VariableModules.o
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/write_Shapefile.f90

tools: ascii_check

ascii_check: $(ASH3dCCSRC)/tools/Ash3d_ASCII_check.f90 Ash3d_VariableModules.o write_ASCII.o $(MFILE)
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/tools/Ash3d_ASCII_check.f90
	$(FC) $(FFLAGS) $(EXFLAGS) $(OBJECTS) $(DATAOBJ) Ash3d_ASCII_check.o $(LIBS) -o Ash3d_ASCII_check
	mkdir -p ../bin/tools
	mv Ash3d_ASCII_check ../bin/tools/Ash3d_ASCII_check

ascii_stats: $(ASH3dCCSRC)/tools/Ash3d_ASCII_stats.f90 Ash3d_VariableModules.o write_ASCII.o $(MFILE)
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) -c $(ASH3dCCSRC)/tools/Ash3d_ASCII_stats.f90
	$(FC) $(FFLAGS) $(EXFLAGS) $(OBJECTS) $(DATAOBJ) Ash3d_ASCII_stats.o $(LIBS) -o Ash3d_ASCII_stats
	mkdir -p ../bin/tools
	mv Ash3d_ASCII_stats ../bin/tools/Ash3d_ASCII_stats

check: ascii_check Ash3d $(MFILE)
	./run_tests.sh

clean:
	rm -rf $(OBJECTS) $(PPOBJS) $(GRAFOBJ) Ash3d.o $(DATAOBJ) $(OPTMODS_OBJ) ../bin/$(ASH3DEXEC) ../bin/tools/Ash3d_ASCII_check \
        *.mod *__genmod.f90

install:
	install -d $(INSTALLDIR)/bin $(INSTALLDIR)/share/post_proc
	install -m 755 ../bin/$(ASH3DEXEC) $(INSTALLDIR)/bin/$(ASH3DEXEC)
	install -m 755 ../bin/Ash3d $(INSTALLDIR)/bin/Ash3d
	install -m 755 ../bin/Ash3d_PostProc $(INSTALLDIR)/bin/Ash3d_PostProc
	install -m 644 ../share/GlobalAirports_ewert.txt $(INSTALLDIR)/share/
	install -m 644 ../share/VotW_ESP_v12_csv.txt  $(INSTALLDIR)/share/
	install -m 644 ../share/post_proc/*png $(INSTALLDIR)/share/post_proc/
	install -m 644 ../share/post_proc/world_cities.txt $(INSTALLDIR)/share/post_proc/

uninstall:
	rm -f $(INSTALLDIR)/bin/$(ASH3DEXEC)
	rm -r $(INSTALLDIR)/share/GlobalAirports_ewert.txt
	rm -f $(INSTALLDIR)/share/VotW_ESP_v12_csv.txt
	rm -f $(INSTALLDIR)/share/post_proc/CloudLoad_hsv.png
	rm -f $(INSTALLDIR)/share/post_proc/USGS_warning3.png
	rmdir $(INSTALLDIR)/share/post_proc
	rmdir $(INSTALLDIR)/share
	rmdir $(INSTALLDIR)/bin
	rmdir $(INSTALLDIR)


