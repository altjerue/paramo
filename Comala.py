import os
from time import strftime, localtime
from SAPytho.misc import fortran_double

# TODO: get rid of b_index

#
#  #####    ##   #####    ##   #    #  ####
#  #    #  #  #  #    #  #  #  ##  ## #
#  #    # #    # #    # #    # # ## #  ####
#  #####  ###### #####  ###### #    #      #
#  #      #    # #   #  #    # #    # #    #
#  #      #    # #    # #    # #    #  ####
class parameters(object):
    '''This is the parameters class
    '''

    # -----  PARAMETERS  -----
    def rParams(self):
        self.R = 1e15                   # radius of emitting region (assuming spherical)
        self.R0 = 1e14                  # distance from central engine
        self.dLum = 4.0793e26           # luminosity distance (default Mrk 421)
        self.z = 0.03                   # redshift (default Mrk 421)
        self.theta_obs = 5.0            # observer viewing angle
        self.gamma_bulk = 1e2           # emitting region bulk Lorentz factor
        self.mu_mag = 1.0               # (1 + sigma) Gamma_bulk
        self.sigma = 1.0                # magnetization (sigma)
        self.f_rec = 1.0                # magnetic reconection dissipative efficiency
        self.b_index = 0.0              # magnetic field decay index
        self.Bfield = 1.0               # magnetic field decay index
        self.eps_B = 0.03               # epsilon_B
        self.eps_e = 0.1                # epsilon_e
        self.theta_e = 10.0             # electrons temperature
        self.zeta_e = 0.99              # fraction of non-thermal particles
        self.tstep = 1e-2               # time step factor
        self.tmax = 1e5                 # maximum time
        self.tmin = 0e0                 # minimum time
        self.tvar = 2e0                 # variability time scale
        self.L_jet = 1e45               # jet luminosity
        self.E0 = 1e50                  # energy of the blast wave
        self.n_ext = 1.0                # number density of the external medium
        self.g1 = 1e2                   # power-law min Lorentz factor
        self.g2 = 1e4                   # power-law max Lorentz factor
        self.gmin = 1.01                # EED minimum Lorentz factor
        self.gmax = 2e4                 # EED maximum Lorentz factor
        self.pind = 2.5                 # EED power-law index
        self.nu_ext = 1e14              # external radiation field freq.
        self.u_ext = 1e-4               # external radiation field ener. dens.
        self.numin = 1e7                # minimum frequency
        self.numax = 1e15               # maximum frequency
        self.numbins = 128              # number of EED bins
        self.numdt = 300                # number of time steps
        self.numdf = 256                # number of frequencies
        self.time_grid = 1              # kind of cooling
        self.params_file = 'input.par'  # name of the parameters file

    def __init__(self, **kwargs):
        self.rParams()
        self.__dict__.update(kwargs)

    def wParams(self):
        with open(self.params_file, 'w') as f:
            print(fortran_double(self.R), ' ! Radius', file=f)
            print(fortran_double(self.R0), ' ! Initial radius', file=f)
            print(fortran_double(self.dLum), ' ! luminosity distance', file=f)
            print(fortran_double(self.z), ' ! redshift', file=f)
            print(fortran_double(self.theta_obs), ' ! viewing angle', file=f)
            print(fortran_double(self.gamma_bulk), ' ! bulk Lorentz factor', file=f)
            print(fortran_double(self.mu_mag), ' ! (1 + sigma) Gamma', file=f)
            print(fortran_double(self.sigma), ' ! magnetization', file=f)
            print(fortran_double(self.f_rec), ' ! dissipative efficiency of magnetic reconection', file=f)
            print(fortran_double(self.b_index), ' ! magnetic field decay index', file=f)
            print(fortran_double(self.Bfield), ' ! magnetic field strength', file=f)
            print(fortran_double(self.theta_e), ' ! electrons temperature', file=f)
            print(fortran_double(self.zeta_e), ' ! fraction of nonthermal electrons', file=f)
            print(fortran_double(self.tstep), ' ! time step factor', file=f)
            print(fortran_double(self.tmax), ' ! maximum time', file=f)
            print(fortran_double(self.tmin), ' ! minimum time', file=f)
            print(fortran_double(self.tvar), ' ! variability time scale', file=f)
            print(fortran_double(self.L_jet), ' ! jet luminosity', file=f)
            print(fortran_double(self.E0), ' ! energy of the blast wave', file=f)
            print(fortran_double(self.n_ext), ' ! number density of the external medium', file=f)
            print(fortran_double(self.eps_e), ' ! epsilon_e', file=f)
            print(fortran_double(self.eps_B), ' ! epsilon_B', file=f)
            print(fortran_double(self.g1), ' ! power-law min Lorentz factor', file=f)
            print(fortran_double(self.g2), ' ! power-law max Lorentz factor', file=f)
            print(fortran_double(self.gmin), ' ! EED min Lorentz factor', file=f)
            print(fortran_double(self.gmax), ' ! EED max Lorentz factor', file=f)
            print(fortran_double(self.pind), ' ! EED power-law index', file=f)
            print(fortran_double(self.nu_ext), ' ! external rad. field frequency', file=f)
            print(fortran_double(self.u_ext), ' ! external rad. field ener. density', file=f)
            print(fortran_double(self.numin), ' ! min frequency', file=f)
            print(fortran_double(self.numax), ' ! max frequency', file=f)
            print(self.numbins, ' ! number of EED bins', file=f)
            print(self.numdt, ' ! number of time steps', file=f)
            print(self.numdf, ' ! number of frequencies', file=f)
            print(self.time_grid, ' ! kind of time grid', file=f)
            # print(self.file_label, ' ! label to identify each output', file=f)
        print("--> Parameters file: ", self.params_file)


#   ####   ####  #    # #####  # #      ######
#  #    # #    # ##  ## #    # # #      #
#  #      #    # # ## # #    # # #      #####
#  #      #    # #    # #####  # #      #
#  #    # #    # #    # #      # #      #
#   ####   ####  #    # #      # ###### ######
class compiler(object):
    '''This is the compilation class
    '''
    # -----  COMPILER FLAGS & RULES -----

    def flags(self):
        self.BROWN = False       # compile with BROWN=1 flag
        self.INTEL = False       # compile with IFORT=1 flag
        self.arch = ''           # compile with specific arch flag
        self.OMP = False         # compile with OpenMP
        self.DBG = False         # compile for debugging
        self.rules = 'all'       # rule to compile
        self.compile_dir = './'  # the path to Paramo... must end with '/'

    def __init__(self, **kwargs):
        self.flags()
        self.__dict__.update(kwargs)
        self.cwd = os.getcwd()

    def compile(self):
        make = 'make ' + self.rules + ' -j4'
        if self.arch in ['i7', 'corei7', 'I7', 'COREI7']:
            make += ' COREI7=1'
        else:
            make += ' NATIVE=1'

        if self.OMP:
            make += ' OPENMP=1'

        if self.DBG:
            make += ' DBG=1'

        if self.BROWN:
            make += ' BROWN=1'

        if self.INTEL:
            make += ' IFORT=1'

        os.chdir(self.compile_dir)
        print(os.getcwd())
        print("--> Running Makefile:\n ", make, "\n")
        log = strftime("%a, %d %b %Y %H:%M:%S %Z", localtime())
        os.system(make)
        with open("make.log", "a") as logfile:
            logfile.write(log + "\n" + make + "\n\n")
        logfile.close()
        os.chdir(self.cwd)
        print(os.getcwd())

    def cleanup(self):
        os.chdir(self.compile_dir)
        os.system("make clean")
        os.chdir(self.cwd)


#  #####  #####   ####     ###### # #      ######
#  #    # #    # #         #      # #      #
#  #    # #####   ####     #####  # #      #####
#  #####  #    #      #    #      # #      #
#  #      #    # #    #    #      # #      #
#  #      #####   ####     #      # ###### ######
def PBSfile(jname, qname, xcmd, depen=None, nodes=None, cores=None, mail=None, htime=None):
    '''This function generates the PBS file to queue a simulation
    '''
    from datetime import timedelta as td

    if htime is None:
        t = str(td(hours=2.0))
    else:
        t = str(td(hours=htime))
    sname = "{0}.sub".format(jname)

    if nodes is None:
        n = 1
    else:
        n = nodes

    if cores is None or qname is 'debug':
        c = 24
    else:
        c = cores

    with open(sname, 'w') as f:
        print("#!/bin/sh -l\n", file=f)
        print("# FILENAME: {0}\n".format(sname), file=f)
        print("#PBS -q " + qname, file=f)
        print("#PBS -l nodes={0:d}:ppn={1:d},naccesspolicy=singleuser".format(n, c), file=f)
        print("#PBS -l walltime={0}".format(t), file=f)
        print("#PBS -N " + jname, file=f)
        print("#PBS -o /scratch/brown/jruedabe/joboutput/{0}.out".format(jname), file=f)
        print("#PBS -e /scratch/brown/jruedabe/joboutput/{0}.err".format(jname), file=f)
        if depen is not None:
            print("#PBS -W depend=afterok:{0}".format(depen), file=f)
        if mail is not None:
            print("#PBS -M {0}".format(mail), file=f)
            print("#PBS -m bae", file=f)
        print("Working at: {0}".format(os.getcwd()))
        print("\ncd {0}".format(os.getcwd()), file=f)
        if (c > 1) or not (qname is 'debug'):
            print("export OMP_NUM_THREADS={0}".format(c), file=f)
        print("\n# RUN", file=f)
        print(xcmd, file=f)
    return sname


#  #####  #    # #    #
#  #    # #    # ##   #
#  #    # #    # # #  #
#  #####  #    # #  # #
#  #   #  #    # #   ##
#  #    #  ####  #    #
class Runner(object):
    '''This class sets up the exectuable instructions.
    '''

    def __init__(self, wCool=False, wAbs=False, wIC=False, flabel='DriverTest', par_kw={}, comp_kw={}):
        self.par = parameters(**par_kw)
        self.par.wParams()
        self.cwd = os.getcwd()
        self.comp_kw = comp_kw
        # -----  ARGS OF THE EXECUTABLE  -----
        if wCool:
            self.wCool = 'T'
        else:
            self.wCool = 'F'
        if wAbs:
            self.wAbs = 'T'
        else:
            self.wAbs = 'F'
        if wIC:
            self.wIC = 'T'
        else:
            self.wIC = 'F'
        self.flabel = flabel  # a label to identify each output

    def run_test(self, clean=False):
        comp = compiler(rules='xTests', **self.comp_kw)
        if clean:
            comp.cleanup()
        comp.compile()
        run_cmd = '{0}xTests'.format(comp.compile_dir)
        print("\n--> Running:\n  ", run_cmd, "\n")
        os.system(run_cmd)
        print("\n--> Finished")

    def run_blazMag(self, pream=None, clean=False, cl=False):
        comp = compiler(rules='xBlazMag', **self.comp_kw)
        if clean:
            comp.cleanup()
        comp.compile()
        outfile = self.flabel + '.jp.h5'
        if pream is None:
            run_cmd = '{0}xBlazMag {1} {2} {3} {4} {5}'.format(comp.compile_dir, self.par.params_file, self.wAbs, self.wIC, self.wCool, outfile)
        else:
            run_cmd = '{0} {1}xBlazMag {2} {3} {4} {5} {6}'.format(pream, comp.compile_dir, self.par.params_file, self.wAbs, self.wIC, self.wCool, outfile)
        print("\n--> Parameters:")
        os.system("cat -n " + self.par.params_file)
        if cl:
            print("\n--> Running:\n  ", run_cmd, "\n")
            return run_cmd
        else:
            print("\n--> Running:\n  ", run_cmd, "\n")
            os.system(run_cmd)
            print("\n--> Finished")

    def run_afterglow(self, pream=None, clean=False, cl=False):
        if clean:
            self.comp.cleanup()
        comp = compiler(rules='xAglow', **self.comp_kw)
        comp.compile()
        outfile = self.flabel + '.jp.h5'
        if pream is None:
            run_cmd = '{0}xAglow {1} {2} {3} {4} {5}'.format(comp.compile_dir, self.wAbs, self.wCool, self.wIC, self.par.params_file, outfile)
        else:
            run_cmd = '{0} {1}xAglow {2} {3} {4} {5} {6}'.format(pream, comp.compile_dir, self.wAbs, self.wCool, self.wIC, self.par.params_file, outfile)
            print("\n--> Parameters:")
        os.system("cat -n " + self.par.params_file)
        if cl:
            return run_cmd
        else:
            print("\n--> Running:\n  ", run_cmd, "\n")
            os.system(run_cmd)
            print("\n--> Finished")
