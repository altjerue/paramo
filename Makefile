ifeq ($(ARC),1)
	FC=h5pfc
else
	FC=h5fc
endif

ifeq ($(MPI),1)
	FC=h5pfc
endif

ifeq ($(IFORT),1)
	ifeq ($(MPI),1)
		FC=ih5pfc
	else
		FC=ih5fc
	endif
endif

ifeq ($(MBS),1)
	ifeq ($(HYB),1)
		OMBS=-DMBS -USHYB -DHYB
	else
		OMBS=-DMBS -USHYB -UHYB
	endif
else
	ifeq ($(HYB),1)
		OMBS=-UMBS -DSHYB -DHYB
	else
		OMBS=-UMBS -USHYB -UHYB
	endif
endif


# optimization level
ifeq ($(DBG),1)
	ifeq ($(IFORT),1)
		OPTIMIZATION=-m64 -g -debug all -check all -implicitnone -warn unused \
		-fp-stack-check -heap-arrays -ftrapuv -check pointers -check bounds -free
	else
		OPTIMIZATION=-g -Wall -ffree-form -ffree-line-length-none -DNONSTCPP \
		-mieee-fp -ffpe-trap=invalid,zero,overflow -fbacktrace -fcheck=all \
		-fbounds-check -fno-unsafe-math-optimizations -frounding-math \
		-fsignaling-nans $(OMBS)
	endif
else
	ifeq ($(IFORT),1)
		ifeq ($(IFAST),1)
			FASTI=-fast
		else
			FASTI=-O5
		endif

		ifeq ($(IPAR),1)
			PARI=-parallel
		endif

		ifeq ($(OPENMP),1)
			OMP=-openmp
		endif

		OPTIMIZATION=-mssse3 -xssse3 $(FASTI) $(PARI) $(OMP) -free $(OMBS)
	else
		ifeq ($(OPENMP),1)
			OMP=-fopenmp -fcheck=all
		endif

		OPTIMIZATION=-O5 -ftree-vectorize -funroll-all-loops -ffree-form \
		-ffree-line-length-none $(OMP) -DNONSTCPP $(OMBS)
	endif
endif

ifeq ($(COREI7),1)
	OPTIMIZATION+=-march=corei7 -mtune=corei7
endif

ifeq ($(NATIVE),1)
	OPTIMIZATION+=-march=native -mtune=native
endif

COPT = -c $(OPTIMIZATION)
LOPT = $(OPTIMIZATION)

# executables
PARAMO = xParamo
ITOBS = xITobs

# dependencies
PARAMO_OBJ = misc.o pwl_integ.o h5_inout.o K2.o SRtoolkit.o anaFormulae.o \
	magnetobrem.o radiation.o Paramo.o paramo_main.o
ITOBS_OBJ = misc.o h5_inout.o K2.o SRtoolkit.o pwl_integ.o IofTobs.o

# rules
all: $(PARAMO) $(ITOBS)

constants.o K2.o pwl_integ.o misc.o h5_inout.o: data_types.o
SRtoolkit.o: data_types.o constants.o K2.o
magnetobrem.o: data_types.o constants.o h5_inout.o misc.o anaFormulae.o pwl_integ.o
IofTobs.o: data_types.o h5_inout.o SRtoolkit.o pwl_integ.o
Paramo.o: data_types.o constants.o misc.o pwl_integ.o h5_inout.o SRtoolkit.o \
	anaFormulae.o magnetobrem.o radiation.o
paramo_main.o: data_types.o misc.o magnetobrem.o
anaFormulae.o: data_types.o constants.o misc.o pwl_integ.o
radiation.o: data_types.o constants.o misc.o pwl_integ.o SRtoolkit.o anaFormulae.o magnetobrem.o

$(PARAMO): data_types.o constants.o $(PARAMO_OBJ)
	$(FC) $(LOPT) -o $@ $^

$(ITOBS): data_types.o constants.o $(ITOBS_OBJ)
	$(FC) $(LOPT) -o $@ $^


%.o: %.F90
	$(FC) $(COPT) $(DEFS) $< -o $@

.PHONY: clean

clean:
	rm -vf *.o *.mod *~ x*
	rm -rvf *.dSYM
