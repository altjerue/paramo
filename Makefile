ifeq ($(IFC),1)

ifeq ($(MPI),1)
FC=ih5pfc
else  # MPI
FC=ih5fc
endif # MPI

else  # IFORT

ifeq ($(MPI),1)
FC=h5pfc
else  # MPI
FC=h5fc
endif # MPI

endif # IFORT

# definitions

ifeq ($(BROWN),1)

endif





# optimization level
ifeq ($(DBG),1)

ifeq ($(IFORT),1)
OPTIMIZATION=-g -debug all -check all -check nostack -warn all -fp-stack-check -heap-arrays \
	-ftrapuv -free
else  # IFORT
OPTIMIZATION=-g -Wall -ffree-form -ffree-line-length-none \
	-ffpe-trap=invalid,zero -fbacktrace -fcheck=all \
	-fbounds-check -fno-unsafe-math-optimizations -frounding-math \
	-fsignaling-nans
endif # IFORT

else # DBG

ifeq ($(IFORT),1)

ifeq ($(IFAST),1)
FASTI=-fast
else  # IFAST
FASTI=-O5
endif # IFAST

ifeq ($(IPAR),1)
PARI=-parallel
endif #IPAR

ifeq ($(OPENMP),1)
OMP=-qopenmp
endif #OPENMP

OPTIMIZATION=-mssse3 -xssse3 -free $(FASTI) $(PARI) $(OMP)

else # IFORT

ifeq ($(OPENMP),1)
OMP=-fopenmp
endif # OPENMP

OPTIMIZATION=-O5 -ftree-vectorize -funroll-all-loops -ffree-form \
	-ffree-line-length-none -fbacktrace $(OMP)

endif # IFORT

ifeq ($(COREI7),1)
OPTIMIZATION+=-march=corei7 -mtune=corei7
endif

ifeq ($(NATIVE),1)
OPTIMIZATION+=-march=native -mtune=native
endif

ifeq ($(BROWN),1)
OPTIMIZATION+=-arch ssse3 -mtune=corei7-avx
DEFS+=-DBRWN
endif

endif # DBG

COPT=-c $(OPTIMIZATION) $(DEFS) $(INCL)
LOPT=$(OPTIMIZATION) $(DEFS) $(LIBS)

# -----  executables  -----
PARAMO=xParamo
ITOBS=xITobs
TESTS=xTests
EBL=xEBL

# -----  dependencies  -----
PARAMO_OBJ = params.o misc.o pwl_integ.o h5_inout.o K1.o K2.o SRtoolkit.o anaFormulae.o \
	radiation.o dist_evol.o Paramo.o paramo_main.o
ITOBS_OBJ = misc.o h5_inout.o K2.o SRtoolkit.o pwl_integ.o anaFormulae.o \
	radiation.o IofTobs.o
TESTS_OBJ = params.o misc.o pwl_integ.o h5_inout.o K1.o K2.o SRtoolkit.o anaFormulae.o \
	radiation.o dist_evol.o tests.o

# -----  rules  -----
all: $(PARAMO) $(ITOBS) $(TESTS)

# objects
constants.o K2.o K1.o pwl_integ.o misc.o h5_inout.o: data_types.o
params.o: data_types.o misc.o
SRtoolkit.o: data_types.o constants.o K2.o
magnetobrem.o: data_types.o constants.o h5_inout.o misc.o anaFormulae.o \
	pwl_integ.o
IofTobs.o: data_types.o h5_inout.o SRtoolkit.o pwl_integ.o radiation.o
Paramo.o: data_types.o constants.o misc.o pwl_integ.o h5_inout.o SRtoolkit.o \
	anaFormulae.o magnetobrem.o radiation.o dist_evol.o K1.o K2.o
paramo_main.o: data_types.o misc.o magnetobrem.o Paramo.o
anaFormulae.o: data_types.o constants.o misc.o pwl_integ.o
radiation.o: data_types.o constants.o misc.o pwl_integ.o SRtoolkit.o \
	anaFormulae.o magnetobrem.o
dist_evol.o: data_types.o constants.o misc.o pwl_integ.o SRtoolkit.o
tests.o: data_types.o constants.o misc.o pwl_integ.o h5_inout.o SRtoolkit.o \
	anaFormulae.o magnetobrem.o radiation.o dist_evol.o K1.o K2.o

# executables
$(PARAMO): data_types.o constants.o $(PARAMO_OBJ)
	$(FC) $(LOPT) -o $@ $^

$(ITOBS): data_types.o constants.o $(ITOBS_OBJ)
	$(FC) $(LOPT) -o $@ $^

$(TESTS): data_types.o constants.o $(TESTS_OBJ)
	$(FC) $(LOPT) -o $@ $^

%.o: %.F90
	$(FC) $(COPT) $(DEFS) $< -o $@

# -----  PHONY things  -----
.PHONY: clean

clean:
	rm -vf *.o *.mod *~

clean_all:
	rm -vf *.o *.mod *~
	rm -rvf x*
