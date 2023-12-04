###############################################################################
##########  Intel Fortran Compiler  ###########################################
SYSFPPFLAG = -DIFORT

FCHOME = /opt/intel/oneapi/compiler/latest/linux/
FC = $(FCHOME)/bin/intel64/ifort

#LAPACKLIB = -L$(FCHOME)/lib64 -llapack
#LAPACKLIB = -L$(HOME)/Programs/Tarballs/lapack-3.11 -llapack
LAPACKLIB = -L/opt/intel/oneapi/mkl/latest/lib/intel64 -qmkl

COMPINC = -I./ -I$(FCHOME)/include
COMPLIBS = -L./ -L$(FCHOME)/lib -L/opt/intel/oneapi/mkl/latest/lib/intel64 -lirc -limf
LIBS = $(LAPACKLIB) $(COMPLIBS) $(USGSLIBDIR) $(USGSINC) $(COMPINC) $(USGSLIB) $(DATALIBS)

# Debugging flags
ifeq ($(RUN), DEBUG)
    FFLAGS = -g2 -pg -warn all -check all -real-size 64 -check uninit -traceback -ftrapuv -debug all 
endif
ifeq ($(RUN), DEBUGOMP)
    FFLAGS = -g2 -pg -warn all -check all -real-size 64 -check uninit -traceback -ftrapuv -debug all -openmp
endif
# Profiling flags
ifeq ($(RUN), PROF)
    FFLAGS = -g2 -pg
endif
# Production run flags
ifeq ($(RUN), OPT)
    FFLAGS = -O3 -ftz -w -ipo
endif
ifeq ($(RUN), OMPOPT)
    FFLAGS = -O3 -ftz -w -ipo -openmp
endif
  # Preprocessing flags
FPPFLAGS =  -fpp $(OSFPPFLAG) $(DATAFPPFLAG) $(LIMFPPFLAG) $(DIFFFPPFLAG) $(memFPPFLAG) $(extdatFPPFLAG) $(FASTFPPFLAG) $(PIIFPPFLAG) $(SYSFPPFLAG) $(PLOTFPPFLAG)
  # Extra flags
EXFLAGS =
###############################################################################
