<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/afterglows.F90">afterglows.F90</a>

### Description: 

Simulates afterglows by, solving the blast-waves' equation of motion, evolving particles inside the emission region, and calculating the emission from the particles.


### Functions/Subroutines:

1. bw1D_afterglow
    
    Simulates afterglows by assuming the blast wave's equation is one dimensional and assumes an external density $$n(r)=A r^{-s}$$. s=0 is a homogenous cicumsteller medium and s=2 describes a constant peed wind from the progenitor.

    For s=0(with_wind=False) see [CD99] for model details. 

    For s>0 see  [PK00]

 
   __parameter__ : __in/out__: __type__

   * *params_file*:in
     * Input file that specifies the computational and physical parameters in the model
     * see https://bitbucket.org/comala/paramo/src/master/Examples/example_input.par
   * *output_file*:out
     * filename of the output hdf5 file
     * see https://bitbucket.org/comala/paramo/src/master/Examples/afterglow_examples.py for a list of outputs in the result classes
   * *with_wind*:in:logical
    
   * *cool_withKN*:in:logical
   
   * *blob*:in:logical


####
 
<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/Arriero.py">Arriero.py</a>


 
 
<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/benchmarking.F90">benchmarking.F90</a>


 
 
<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/benchmarks.F90">benchmarks.F90</a>


 
 
<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/blastwave.F90">blastwave.F90</a>


 
 
<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/blazMag.F90">blazMag.F90</a>


 
 
<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/constants.F90">constants.F90</a>


 
 
<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/data_types.F90">data_types.F90</a>


 
 
<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/distribs.F90">distribs.F90</a>


 
 
<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/h5_inout.F90">h5_inout.F90</a>


 
 
<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/interJets.F90">interJets.F90</a>


 
 
<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/main.F90">main.F90</a>


 
 
<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/misc.F90">misc.F90</a>


 
 
<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/pairs.F90">pairs.F90</a>


 
 
<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/params.F90">params.F90</a>


 
 
<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/pwl_integ.F90">pwl_integ.F90</a>


 
 
<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/radiation.F90">radiation.F90</a>


 
 
<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/specialf.F90">specialf.F90</a>


 
 
<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/SRtoolkit.F90">SRtoolkit.F90</a>


 
 
<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/tests.F90">tests.F90</a>


 
 
<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/Tests_Plots.py">Tests_Plots.py</a>


 
 
<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/transformers.F90">transformers.F90</a>


 
 
<a style="font-size:14px" href="https://bitbucket.org/comala/paramo/src/master/turBlaz.F90">turBlaz.F90</a>


 
 
