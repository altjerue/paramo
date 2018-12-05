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


ifeq ($(BROWN),1)
OBWN=-DBRWN
# else
# OBWN=-UBRWN
endif


ifeq ($(MBS),1)

ifeq ($(HYB),1)
OMBS=-DMBS -DHYB
else  # HYB
OMBS=-DMBS #-UHYB
endif # HYB

else  # MBS

ifeq ($(HYB),1)
OMBS=-DHYB #-UMBS
#else  # HYB
#OMBS=-UMBS -UHYB
endif # HYB

endif # MBS

ifeq ($(FBAR),1)
OFB=-DFBAR
INCL=-I/Users/jesus/lib/Fortran/forbear/static/mod/
LIBS=-L/Users/jesus/lib/Fortran/forbear/static/ -lforbear
#else
#OFB=-UFBAR
endif


# optimization level
ifeq ($(DBG),1)

ifeq ($(IFORT),1)
OPTIMIZATION=-g -debug all -check all -check nostack -warn all -fp-stack-check -heap-arrays \
	-ftrapuv -free $(OMBS) $(OBWN) $(OFB)
else  # IFORT
OPTIMIZATION=-g -Wall -ffree-form -ffree-line-length-none -DNONSTCPP \
	-mieee-fp -ffpe-trap=invalid,zero,overflow -fbacktrace -fcheck=all \
	-fbounds-check -fno-unsafe-math-optimizations -frounding-math \
	-fsignaling-nans $(OMBS) $(OBWN) $(OFB)
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

OPTIMIZATION=-mssse3 -xssse3 -free $(FASTI) $(PARI) $(OMP) $(OMBS) $(OFB) $(OBWN)

else # IFORT

ifeq ($(OPENMP),1)
OMP=-fopenmp
endif # OPENMP

OPTIMIZATION=-O5 -ftree-vectorize -funroll-all-loops -ffree-form \
	-ffree-line-length-none $(OMP) -DNONSTCPP $(OMBS) $(OFB) $(OBWN)

endif # IFORT

ifeq ($(COREI7),1)
OPTIMIZATION+=-march=corei7 -mtune=corei7
endif

ifeq ($(NATIVE),1)
OPTIMIZATION+=-march=native -mtune=native
endif

ifeq ($(BROWN),1)
OPTIMIZATION+=-arch ssse3 -mtune=corei7-avx
endif

endif # DBG



COPT=-c $(OPTIMIZATION) $(INCL)
LOPT=$(OPTIMIZATION) $(LIBS)

# -----  executables  -----
PARAMO=xParamo
ITOBS=xITobs

# -----  dependencies  -----
PARAMO_OBJ = misc.o pwl_integ.o h5_inout.o K2.o SRtoolkit.o anaFormulae.o \
	magnetobrem.o radiation.o dist_evol.o Paramo.o paramo_main.o
ITOBS_OBJ = misc.o h5_inout.o K2.o SRtoolkit.o pwl_integ.o IofTobs.o

# -----  rules  -----
all: $(PARAMO) $(ITOBS)

# objects
constants.o K2.o pwl_integ.o misc.o h5_inout.o: data_types.o
SRtoolkit.o: data_types.o constants.o K2.o
magnetobrem.o: data_types.o constants.o h5_inout.o misc.o anaFormulae.o \
	pwl_integ.o
IofTobs.o: data_types.o h5_inout.o SRtoolkit.o pwl_integ.o
Paramo.o: data_types.o constants.o misc.o pwl_integ.o h5_inout.o SRtoolkit.o \
	anaFormulae.o magnetobrem.o radiation.o dist_evol.o
paramo_main.o: data_types.o misc.o magnetobrem.o Paramo.o
anaFormulae.o: data_types.o constants.o misc.o pwl_integ.o
radiation.o: data_types.o constants.o misc.o pwl_integ.o SRtoolkit.o \
	anaFormulae.o magnetobrem.o
dist_evol.o: data_types.o constants.o misc.o pwl_integ.o SRtoolkit.o

# executables
$(PARAMO): data_types.o constants.o $(PARAMO_OBJ)
	$(FC) $(LOPT) -o $@ $^

$(ITOBS): data_types.o constants.o $(ITOBS_OBJ)
	$(FC) $(LOPT) -o $@ $^


%.o: %.F90
	$(FC) $(COPT) $(DEFS) $< -o $@

# -----  PHONY things  -----
.PHONY: clean

clean:
	rm -vf *.o *.mod *~ x*
	rm -rvf *.dSYM
