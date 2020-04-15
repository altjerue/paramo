# definitions Is this a compiler?
FC=h5fc #if not using h5c this would be gfortran
# -----  executables  -----this creates the files names for the eventual compiled programs

TESTS=xturbulentemission

OPTIMIZATION=-O3 -ftree-vectorize -funroll-all-loops -ffree-form \
	-ffree-line-length-none -fbacktrace

COPT=-c $(OPTIMIZATION)

# -----  dependencies  ----- gets everything from use in the F90.slash represents inputs and output files. idk why hdf5 isnt here maybe it apart of fortran
TESTS_OBJ =misc.o data_types.o h5_inout.o dist_evol.o constants.o\
	pwl_integ.o SRtoolkit.o K2.o anaFormulae.o radiation.o turbulentemission.o
# -----  rules  ----- im not sure what this is yet ... this might be way of introducing short hand such as make all
all: $(TESTS)
# objects this seems to further break down dependencies for other scripts you have
constants.o K2.o K1.o pwl_integ.o misc.o h5_inout.o: data_types.o
dist_evol.o: data_types.o constants.o misc.o pwl_integ.o SRtoolkit.o K2.o
anaFormulae.o: data_types.o constants.o misc.o pwl_integ.o
radiation.o: data_types.o constants.o misc.o pwl_integ.o SRtoolkit.o \
	anaFormulae.o

# executables creating the programs from above
$(TESTS): $(TESTS_OBJ)
	echo "compiling $@ with these objects: $^"
	$(FC) -o $@ $^


#-o $@ $^ saying create all .o files listed in the varible name I think?? ask jesus
#$(FC) something to with hdf5
#not sure what this is jesus $(LOPT)

%.o: %.F90 #maybe all .0 outputs come from a corresponding .F90 input?
	$(FC) $(COPT) $< -o $@




# -----  PHONY things  -----
.PHONY: clean

clean:
	rm -vf *.o *.mod *~
	rm -rvf x*
