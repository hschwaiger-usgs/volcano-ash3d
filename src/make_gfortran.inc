###############################################################################
##########  GNU Fortran Compiler  #############################################
FCHOME=/usr
FC=/usr/bin/gfortran
#COMPINC = -I$(FCHOME)/include -I$(FCHOME)/lib64/gfortran/modules
#COMPLIBS = -L$(FCHOME)/lib64 -L$(FCHOME)/lib

#LAPACKLIB = -L$(FCHOME)/lib64 -llapack
LAPACKLIB = -L$(HOME)/Programs/Tarballs/lapack-3.11 -llapack

COMPINC = -I$(FCHOME)/lib64/gfortran/modules
COMPLIBS = -L$(FCHOME)/lib64 
LIBS = $(LAPACKLIB) $(COMPLIBS) $(USGSLIBDIR) $(USGSINC) $(COMPINC) $(USGSLIB) $(DATALIBS)
# -lefence

# Debugging flags
ifeq ($(RUN), DEBUG)
    FFLAGS = -O0 -g3 -Wall -Wextra -fimplicit-none -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fwhole-file -fcheck=all -pedantic -fbacktrace -Wunderflow -ffpe-trap=invalid,zero,overflow -Wuninitialized -Wmaybe-uninitialized
#-fdefault-real-8 -std=f2008 
endif
ifeq ($(RUN), DEBUGOMP)
    FFLAGS = -O0 -g3 -Wall -Wextra -fimplicit-none -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fwhole-file -fcheck=all -pedantic -fbacktrace -Wunderflow -ffpe-trap=invalid,zero,overflow -Wuninitialized -Wmaybe-uninitialized -fopenmp -lgomp
endif
# Profiling flags
ifeq ($(RUN), PROF)
    FFLAGS = -g -pg -w -fno-math-errno -funsafe-math-optimizations -fno-trapping-math -fno-signaling-nans -fcx-limited-range -fno-rounding-math -fdefault-real-8
endif
# Production run flags
ifeq ($(RUN), OPT)
    FFLAGS = -O3 -w -fno-math-errno -funsafe-math-optimizations -fno-trapping-math -fno-signaling-nans -fcx-limited-range -fno-rounding-math -fdefault-real-8
endif
ifeq ($(RUN), OMPOPT)
    FFLAGS = -O3 -w -fno-math-errno -funsafe-math-optimizations -fno-trapping-math -fno-signaling-nans -fcx-limited-range -fno-rounding-math -fdefault-real-8 -fopenmp -lgomp
endif
  # Preprocessing flags
FPPFLAGS =  -x f95-cpp-input $(OSFPPFLAG) $(DATAFPPFLAG) $(LIMFPPFLAG) $(DIFFFPPFLAG) $(memFPPFLAG) $(extdatFPPFLAG) $(FASTFPPFLAG) $(PIIFPPFLAG) $(SYSFPPFLAG) $(PLOTFPPFLAG)
  # Extra flags
EXFLAGS =
###############################################################################

