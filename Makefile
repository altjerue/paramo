# -----  PHONY things  -----
.PHONY: clean

#-----  compiler level  -----
# COMPILER: 0 (GCC), 1 (INTEL)
COMPILER = 0

USEHDF5 = 0
USEMPI = 0

ifeq ($(USEHDF5),1)
DEFS=-DHDF5
ifeq ($(USEMPI),1)
FC=h5pfc

else
FC=h5fc

endif

else
DEFS=-UHDF5
ifeq ($(COMPILER),1)
FC=ifort

else
FC=gfortran

endif

endif

#-----  servers level  -----
# SERVER: 0 (UNIX PC) 1 (Brown@Purdue)
SERVER = 0
ifeq ($(SERVER),1)
DEFS+=-DBRWN
LIBS+=-L/usr/lib64

endif

#-----  optimization level  -----
# DEBUGGING: 0 (OFF), 1 (ON)
ifeq ($(DBG),1)
ifeq ($(COMPILER),1)
OPTIMIZATION=-g -debug all -check all -check nostack -warn all -fp-stack-check \
	-heap-arrays -ftrapuv -free

else
OPTIMIZATION=-g -Wall -ffree-form -ffree-line-length-none \
	-ffpe-trap=invalid,zero -fbacktrace -fcheck=all \
	-fbounds-check -fno-unsafe-math-optimizations -frounding-math \
	-fsignaling-nans

endif

else
ifeq ($(COMPILER),1)
OPTIMIZATION=-free -O3

ifeq ($(OPENMP),1)
OPTIMIZATION+=-qopenmp

endif

else
OPTIMIZATION=-O3 -ftree-vectorize -funroll-all-loops -ffree-form \
	-ffree-line-length-none -fbacktrace
ifeq ($(OPENMP),1)
OPTIMIZATION+=-fopenmp

endif

endif

endif

COPT=-c $(OPTIMIZATION) -cpp -dU $(DEFS) $(INCL)
LOPT=$(OPTIMIZATION) -cpp -dU $(DEFS) $(LIBS)


# -----  executables & dependencies  -----
# PROBLEM: 0 (tests), 1 (blazars), 2 (afterflow), 3 (turbulence)
PROBLEM=2
OBJECTS=misc.o params.o pwl_integ.o specialf.o SRtoolkit.o anaFormulae.o\
	radiation.o distribs.o transformers.o

ifeq ($(USEHDF5), 1)
OBJECTS+=h5_inout.o
endif

ifeq ($(PROBLEM),0)
OBJECTS+=tests.o
endif

ifeq ($(PROBLEM),1)
OBJECTS+=blazMag.o blazMag_main.o
endif

ifeq ($(PROBLEM),2) 
OBJECTS+=pairs.o blastwave.o afterglows.o afterglow_main.o
# WITHMEZCAL: 0 (1D blast-wave), 1 (Mezcal)
WITHMEZCAL=0
ifeq ($(WITHMEZCAL),1)
DEFS+=-DMEZCAL
else
DEFS+=-UMEZCAL
endif

endif

ifeq ($(PROBLEM),3)
OBJECTS+=turBlaz.o turBlaz_main.o
endif

# BENCH_OBJ=transformers.o misc.o pwl_integ.o SRtoolkit.o specialf.o \
# 	anaFormulae.o blastwave.o radiation.o distribs.o benchmarks.o \
# 	benchmarking.o


# -----  rules level  -----
all: Paramo

################################################################################
# WARNING: All below is DEPRECATED
# objects
# constants.o specialf.o pwl_integ.o misc.o h5_inout.o: data_types.o
# transformers.o SRtoolkit.o: data_types.o constants.o
# pairs.o: data_types.o constants.o misc.o
# params.o: data_types.o misc.o
# blastwave.o: data_types.o constants.o transformers.o SRtoolkit.o
# blazMag.o: data_types.o constants.o misc.o pwl_integ.o SRtoolkit.o \
# 	anaFormulae.o radiation.o distribs.o specialf.o
# blazMag_main.o: data_types.o misc.o blazMag.o
# afterglow_main.o: data_types.o misc.o afterglow.o
# turBlaz.o: data_types.o constants.o misc.o pwl_integ.o SRtoolkit.o \
# 	anaFormulae.o radiation.o distribs.o specialf.o
# turBlaz_main.o: data_types.o misc.o turBlaz.o
# anaFormulae.o: data_types.o constants.o misc.o pwl_integ.o
# radiation.o: data_types.o constants.o misc.o pwl_integ.o SRtoolkit.o \
# 	anaFormulae.o
# distribs.o: data_types.o constants.o misc.o pwl_integ.o SRtoolkit.o specialf.o
# tests.o: data_types.o constants.o misc.o pwl_integ.o SRtoolkit.o \
# 	anaFormulae.o radiation.o distribs.o specialf.o
# afterglows.o: data_types.o constants.o misc.o pwl_integ.o specialf.o\
# 	blastwave.o SRtoolkit.o anaFormulae.o radiation.o distribs.o
# benchmarks.o: data_types.o constants.o transformers.o misc.o pwl_integ.o \
# 	h5_inout.o specialf.o SRtoolkit.o blastwave.o anaFormulae.o distribs.o \
# 	radiation.o
# benchmarking.o: benchmarks.o
################################################################################

# -----  executables production  -----
Paramo: data_types.o constants.o $(OBJECTS)
	$(FC) $(LOPT) -o $@ $^

# -----  objects compilation -----
%.o: %.F90
	$(FC) $(COPT) $< -o $@

clean:
	rm -vf *.o *.mod *~
	rm -rvf x*
