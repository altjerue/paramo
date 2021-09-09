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
# SERVER: 0 (UNIX PC) 1 (Brown@Purdue) 2 (RC@RIT)
SERVER = 0
ifeq ($(SERVER),1)
LIBS:=-L/usr/lib64
endif

ifeq ($(SERVER),2)
ifeq ($(USEHDF5),1)
LIBS:=-L/.autofs/tools/spack/opt/spack/linux-rhel7-skylake_avx512/gcc-7.4.0/hdf5-1.10.5-cplfsiocth2wt5a24qp7uasw4tz5vyjm/lib
INCL:=-I/.autofs/tools/spack/opt/spack/linux-rhel7-skylake_avx512/gcc-7.4.0/hdf5-1.10.5-cplfsiocth2wt5a24qp7uasw4tz5vyjm/include
FC=h5pfc
endif

endif


################################################################################
#  optimization level
ifeq ($(DEBUGGING),1)
ifeq ($(COMPILER),1)
OPTIMIZATION=-g -debug all -check all -check nostack -warn all -fp-stack-check\
	-heap-arrays -ftrapuv -free

else #COMPILER
OPTIMIZATION=-g -Wall -ffree-form -ffree-line-length-none -fcheck=all\
	-ffpe-trap=invalid,zero -fbacktrace -fbounds-check -fsignaling-nans\
	-fno-unsafe-math-optimizations -frounding-math

endif #COMPILER

else #DEBUGGING
ifeq ($(COMPILER),1)
OPTIMIZATION=-free -O3

ifeq ($(OPENMP),1)
OPTIMIZATION+=-qopenmp

endif #OPENMP

else #COMPILER
OPTIMIZATION=-O3 -ftree-vectorize -funroll-all-loops -ffree-form -fbacktrace\
	-ffree-line-length-none
ifeq ($(OPENMP),1)
OPTIMIZATION+=-fopenmp

endif #OPENMP

endif #COMPILER

endif #DEBUGGING

COPT=-c $(OPTIMIZATION) -cpp -dU $(DEFS) $(INCL)
LOPT=$(OPTIMIZATION) -cpp -dU $(DEFS) $(LIBS)


################################################################################
#  executables & dependencies
OBJECTS=misc.o params.o pwl_integ.o specialf.o SRtoolkit.o radiation.o\
	distribs.o transformers.o

ifeq ($(USEHDF5), 1)
OBJECTS+=h5_inout.o
endif

ifndef CONFIG
# CONFIG: 0 (tests), 1 (blazars), 2 (afterflow 1D blast-wave), 3 (turbulence),
#         4 (afterglow wit Mezcal)
CONFIG=0
endif

ifeq ($(CONFIG),0)
#  The test is specified in the macro TEST_CHOICE in the file tests.F90
OBJECTS+=blastwave.o tests.o
DEFS+=-DTEST
endif

ifeq ($(CONFIG),1)
OBJECTS+=blazMag.o main.o #blazMag_main.o
DEFS+=-DBLAZ
endif

ifeq ($(CONFIG),2)
OBJECTS+=pairs.o blastwave.o afterglows.o main.o #afterglow_main.o
DEFS+=-DAGLOW
endif

ifeq ($(CONFIG),3)
OBJECTS+=turBlaz.o main.o #turBlaz_main.o
DEFS+=-DTURB
endif

ifeq ($(CONFIG),4)
OBJECTS+=pairs.o blastwave.o afterglows.o main.o #afterglow_main.o
DEFS+=-DMEZCAL
endif

# BENCH_OBJ=transformers.o misc.o pwl_integ.o SRtoolkit.o specialf.o\
# 	blastwave.o radiation.o distribs.o benchmarks.o benchmarking.o


################################################################################
#  rules level
all: Paramo

Paramo: data_types.o constants.o $(OBJECTS)
	$(FC) $(LOPT) -o $@ $^

%.o: %.F90
	$(FC) $(COPT) $< -o $@

clean:
	rm -vf *.o *.mod *~ Paramo
