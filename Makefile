ifeq ($(ARC),1)
	FC = h5pfc
else
	FC = h5fc
endif

ifeq ($(MPI),1)
	FC = h5pfc
endif

ifeq ($(IFORT),1)
	ifeq ($(MPI),1)
		FC = ih5pfc
	else
		FC = ih5fc
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
		OPTIMIZATION = -m64 -g -debug all -check all -implicitnone -warn unused\
		-fp-stack-check -heap-arrays -ftrapuv -check pointers\
		-check bounds -free
	else
		OPTIMIZATION = -g -Wall -ffree-form -ffree-line-length-none -DNONSTCPP \
		-mieee-fp -ffpe-trap=invalid,zero,overflow \
		-fbacktrace -fcheck=all -fbounds-check -fno-unsafe-math-optimizations \
		-frounding-math -fsignaling-nans $(OMBS)
	endif
else
	ifeq ($(IFORT),1)
		ifeq ($(IFAST),1)
			FASTI = -fast
		else
			FASTI = -O5
		endif

		ifeq ($(IPAR),1)
			PARI = -parallel
		endif

		ifeq ($(OMP),1)
			OMPI = -openmp
		endif

		OPTIMIZATION = -mssse3 -xssse3 $(FASTI) $(PARI) $(OMPI) -free $(OMBS)
	else
		ifeq ($(OMP),1)
			OMPI = -fopenmp
		endif

		OPTIMIZATION = -O5 -ftree-vectorize \
		-funroll-all-loops -ffree-form -ffree-line-length-none $(OMPI) \
		-DNONSTCPP $(OMBS)
	endif
endif

ifeq ($(COREI7),1)
        OPTIMIZATION+= -march=corei7 -mtune=corei7
	ifeq ($(OPENMP),1)
		OPTIMIZATION+= -fopenmp
	endif
endif

ifeq ($(NATIVE),1)
        OPTIMIZATION+= -march=native -mtune=native
	ifeq ($(OPENMP),1)
		OPTIMIZATION+= -fopenmp
	endif
endif

COPT = -c $(OPTIMIZATION)
LOPT = $(OPTIMIZATION)

# executables
ITOBS = xITobs

# dependencies
ITOBS_OBJ = misc.o h5_inout.o K2.o SRtoolkit.o pwl_integ.o IofTobs.o

# rules
all: $(ITOBS)
	@echo " ---------- " >> make.log
	@echo Compilation on: >> make.log
	@date >> make.log
	@echo Compiled with: >> make.log
	@echo $(MAKEFLAGS) >> make.log
	@echo " ---------- " >> make.log


SRtoolkit.o h5_inout.o K2.o IofTobs.o pwl_integ.o constants.o: data_types.o
SRtoolkit.o IofTobs.o: constants.o
SRtoolkit.o: K2.o
IofTobs.o: SRtoolkit.o h5_inout.o pwl_integ.o

$(ITOBS): data_types.o constants.o $(ITOBS_OBJ)
	$(FC) $(LOPT) -o $@ $^


%.o: %.F90
	$(FC) $(COPT) $(DEFS) $< -o $@

.PHONY: clean

clean:
	rm -vf *.o *.mod *~ x*
	rm -rvf *.dSYM
