################################################################################
#  PHONY things
.PHONY: clean


################################################################################
#  compiler level
ifndef COMPILER
# COMPILER: 0 (GCC), 1 (INTEL)
COMPILER = 0
endif

ifndef DEBUGGING
# DEBUGGING: 0 (OFF), 1 (ON)
DEBUGGING = 0
endif

ifndef USEHDF5
# USEHDF5: 0 (OFF), 1 (ON)
USEHDF5 = 0
endif

# ifndef USEMPI
# # USEMPI: 0 (OFF), 1 (ON)
# USEMPI = 0
# endif

### Choosing the compiler
ifeq ($(COMPILER),1)
FC=ifort
DEFS+=-DINTEL

else
DEFS+=-UINTEL
ifeq ($(USEHDF5),1)
DEFS+=-DHDF5
FC=h5fc

else
DEFS=-UHDF5
FC=gfortran

endif

endif


################################################################################
#  servers and libraries level
# SERVER: 0 (UNIX PC) 1 (Brown@Purdue)
SERVER = 0
ifeq ($(SERVER),1)
LIBS+=-L/usr/lib64

endif


################################################################################
#  optimization level
ifeq ($(DEBUGGING),1)
ifeq ($(COMPILER),1)
OPTIMIZATION=-g -debug all -check all -check nostack -warn all -fp-stack-check\
	-heap-arrays -ftrapuv -free

else
OPTIMIZATION=-g -Wall -ffree-form -ffree-line-length-none -fcheck=all\
	-ffpe-trap=invalid,zero -fbacktrace -fbounds-check -fsignaling-nans\
	-fno-unsafe-math-optimizations -frounding-math

endif

else
ifeq ($(COMPILER),1)
OPTIMIZATION=-free -O3

ifeq ($(OPENMP),1)
OPTIMIZATION+=-qopenmp

endif

else
OPTIMIZATION=-O3 -ftree-vectorize -funroll-all-loops -ffree-form -fbacktrace\
	-ffree-line-length-none
ifeq ($(OPENMP),1)
OPTIMIZATION+=-fopenmp

endif

endif

endif

COPT=-c $(OPTIMIZATION) -cpp -dU $(DEFS) $(INCL)
LOPT=$(OPTIMIZATION) -cpp -dU $(DEFS) $(LIBS)


################################################################################
#  executables & dependencies
OBJECTS=misc.o params.o pwl_integ.o specialf.o SRtoolkit.o anaFormulae.o\
	radiation.o distribs.o transformers.o

ifeq ($(USEHDF5), 1)
OBJECTS+=h5_inout.o
endif

ifndef PROBLEM
# PROBLEM: 0 (tests), 1 (blazars), 2 (afterflow 1D blast-wave), 3 (turbulence),
#          4 (afterglow wit Mezcal)
PROBLEM=0
endif

ifeq ($(PROBLEM),0)
OBJECTS+=tests.o
endif

ifeq ($(PROBLEM),1)
OBJECTS+=blazMag.o blazMag_main.o
endif

ifeq ($(PROBLEM),2) 
OBJECTS+=pairs.o blastwave.o afterglows.o afterglow_main.o
DEFS+=-UMEZCAL
endif

endif

ifeq ($(PROBLEM),3)
OBJECTS+=turBlaz.o turBlaz_main.o
endif

ifeq ($(PROBLEM),4)
OBJECTS+=pairs.o blastwave.o afterglows.o afterglow_main.o
DEFS+=-DMEZCAL
else

# BENCH_OBJ=transformers.o misc.o pwl_integ.o SRtoolkit.o specialf.o\
# 	anaFormulae.o blastwave.o radiation.o distribs.o benchmarks.o\
# 	benchmarking.o


################################################################################
#  rules level
all: Paramo

Paramo: data_types.o constants.o $(OBJECTS)
	$(FC) $(LOPT) -o $@ $^

%.o: %.F90
	$(FC) $(COPT) $< -o $@

clean:
	rm -vf *.o *.mod *~ Paramo
