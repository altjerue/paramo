# definitions Is this a compiler?
FC=h5fc #if not using h5c this would be gfortran
# -----  executables  -----this creates the files names for the eventual compiled programs

TESTS=xTestszkd


# -----  dependencies  ----- gets everything from use in the F90.slash represents inputs and output files. idk why hdf5 isnt here maybe it apart of fortran
TESTS_OBJ = data_types.o h5_inout.o dist_evol.o zkdtests0.o constants.o misc.o \
	pwl_integ.o SRtoolkit.o K2.o radiation.o anaFormulae.o
# -----  rules  ----- im not sure what this is yet ... this might be way of introducing short hand such as make all

# objects this seems to further break down dependencies for other scripts you have
dist_evol.o: data_types.o constants.o misc.o pwl_integ.o SRtoolkit.o K2.o
# executables creating the programs from above
$(TESTS): $(TESTS_OBJ)
	echo "compiling $@ with these objects: $^"
	$(FC) -o $@ $^


#-o $@ $^ saying create all .o files listed in the varible name I think?? ask jesus
#$(FC) something to with hdf5
#not sure what this is jesus $(LOPT)

%.o: %.F90 #maybe all .0 outputs come from a corresponding .F90 input?
	$(FC) -c $< -o $@




# -----  PHONY things  -----
.PHONY: clean

clean:
	rm -vf *.o *.mod *~
	rm -rvf x*
