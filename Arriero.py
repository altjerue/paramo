# –¿Qué es lo que pasa, doña Eduviges?
# Ella sacudieo la cabeza como si despertara de un sueño
# –Es el caballo de Miguel Páramo, que galopa por el camino de la Media Luna.

import os
import subprocess
import sys
from time import strftime, localtime
from Eduviges.misc import fortran_double


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
        self.eps_acc = 1e-2             # epsilon_acc
        self.theta_e = 10.0             # electrons temperature
        self.zeta_e = 0.99              # fraction of non-thermal particles
        self.f_esc = 1.0                # electrons escape time factor
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
        self.gmax = self.g2*5e3                 # EED maximum Lorentz factor
        self.pind = 2.5                 # EED power-law index
        self.nu_ext = 1e14              # external radiation field freq.
        self.u_ext = 1e-4               # external radiation field ener. dens.
        self.numin = 1e7                # minimum frequency
        self.numax = 1e15               # maximum frequency
        self.NG = 128                   # number of EED bins
        self.NT = 300                   # number of time steps
        self.NF = 256                   # number of frequencies
        self.time_grid = 1              # kind of cooling
        self.lg1 = "d"              # logical parameter 1 with default = 'd'
        self.lg2 = "d"                # logical parameter 2 with default = 'd'
        self.lg3 = "d"                # logical parameter 3 with default = 'd'
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
            print(fortran_double(self.f_esc), ' ! electrons escape time', file=f)
            print(fortran_double(self.tstep), ' ! time step factor', file=f)
            print(fortran_double(self.tmax), ' ! maximum time', file=f)
            print(fortran_double(self.tmin), ' ! minimum time', file=f)
            print(fortran_double(self.tvar), ' ! variability time scale', file=f)
            print(fortran_double(self.L_jet), ' ! jet luminosity', file=f)
            print(fortran_double(self.E0), ' ! energy of the blast wave', file=f)
            print(fortran_double(self.n_ext), ' ! number density of the external medium', file=f)
            print(fortran_double(self.eps_e), ' ! epsilon_e', file=f)
            print(fortran_double(self.eps_B), ' ! epsilon_B', file=f)
            print(fortran_double(self.eps_acc), ' ! epsilon_acc', file=f)
            print(fortran_double(self.g1), ' ! power-law min Lorentz factor', file=f)
            print(fortran_double(self.g2), ' ! power-law max Lorentz factor', file=f)
            print(fortran_double(self.gmin), ' ! EED min Lorentz factor', file=f)
            print(fortran_double(self.gmax), ' ! EED max Lorentz factor', file=f)
            print(fortran_double(self.pind), ' ! EED power-law index', file=f)
            print(fortran_double(self.nu_ext), ' ! external rad. field frequency', file=f)
            print(fortran_double(self.u_ext), ' ! external rad. field ener. density', file=f)
            print(fortran_double(self.numin), ' ! min frequency', file=f)
            print(fortran_double(self.numax), ' ! max frequency', file=f)
            print(self.NG, ' ! number of EED bins', file=f)
            print(self.NT, ' ! number of time steps', file=f)
            print(self.NF, ' ! number of frequencies', file=f)
            print(self.time_grid, ' ! kind of time grid', file=f)
            print( self.lg1, ' !  logical parameter 1 with default=d', file=f)
            print( self.lg2, ' !  logical parameter 2 with default=d', file=f)
            print( self.lg3, ' !  logical parameter 3 with default=d', file=f)
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
    This function sets the value of the compilation options. The value of these
    options can be changed as argumens or kwargs of the compiler class.

    Parameters
    ----------
    COMP : int, optional
        Compiler to be used. 0 for GCC (default), 1 for Intel.
    OMP : bool, optional
        If True, compilation is done with OpenMP flag. Default False.
    DBG : bool, optional
        If True compilation is done with debugging flags. Default False.
    HDF5 : bool, optional
        If True (default), data saved in HDF5 data files. Default False.
        Note: other formats need to be included
    CONFIG : int, optional
        Program to be compiled. 0 for tests, 1 for blazars, 2 for afterflow, 3 for turbulence, 4 for Mezcal
    server : int, optional
        Computer where code whill be compiled. 0 for unix PC, 1 for Brown (Purdue) server, 2 for SPORC (RC-RIT)
    compileDir : str, optional
        Full path where Paramo is located. Must end with '/'. Default is './'
    '''

    # -----  COMPILER FLAGS & RULES -----
    def flags(self):
        self.COMP = 0           # 0 (GCC), 1 (INTEL)
        self.OMP = False        # compile with OpenMP
        self.DBG = False        # compile for debugging
        self.HDF5 = True        # save data with HDF5
        self.CONFIG = 0         # 0 (tests), 1 (blazars), 2 (afterflow), 3 (turbulence), 4 (Mezcal), 5 (interJets)
        self.server = 0         # 0 (UNIX PC) 1 (Brown@Purdue) 2 (RC@RIT)
        self.compileDir = './'  # the path to Paramo... must end with '/'

    def __init__(self, **kwargs):
        self.flags()
        self.__dict__.update(kwargs)
        self.cwd = os.getcwd()

    def compile(self):
        print(strftime("\n\n%a, %d %b %Y %H:%M:%S %Z", localtime()))
        # make = ["make", "Paramo", "CONFIG="+str(self.problem), "COMPILER="+str(self.COMP)]
        make = ["make", "Paramo", "CONFIG="+str(self.CONFIG), "COMPILER="+str(self.COMP)]
        if self.OMP:
            make.append("OPENMP=1")
        else:
            make.append("OPENMP=0")
        if self.DBG:
            make.append("DEBUGGING=1")
        else:
            make.append("DEBUGGING=0")
        if self.HDF5:
            make.append("USEHDF5=1")
        else:
            make.append("USEHDF5=0")
        os.chdir(self.compileDir)
        print("\nCompile dir: " + os.getcwd())
        print("--> Running Makefile\n ", " ".join(make), "\n")
        try:
            make_out = subprocess.run(make, capture_output=True, text=True, check=True)
            with open("make.out", "w") as logfile:
                logfile.write(make_out.stdout)
        except subprocess.CalledProcessError as err:
            print(err.stdout)
            print("Compilation Error")
            print(err.stderr)
            sys.exit(err.returncode)
        os.chdir(self.cwd)
        print("Working dir: " + os.getcwd())

    def cleanup(self):
        os.chdir(self.compileDir)
        os.system("make clean")
        os.chdir(self.cwd)


#  #####  #    # #    #
#  #    # #    # ##   #
#  #    # #    # # #  #
#  #####  #    # #  # #
#  #   #  #    # #   ##
#  #    #  ####  #    #
class Runner(object):
    '''
    Description
    -----------
    This class sets up the exectuable instructions.
    '''

    def __init__(self, flabel='DriverTest', par_kw={}, comp_kw={}):
        self.par = parameters(**par_kw)
        self.par.wParams()
        self.cwd = os.getcwd()
        self.comp_kw = comp_kw
        self.flabel = flabel  # a label to identify each output

    def run_test(self, clean=False, test_choice=0):
        if(test_choice==0):
            print("\n --> need to input test_choice")
            print("\n --> test_choice : 1 (steady_state), 2 (rad_procs), "
                  "3 (BlackBody_tests), 4 (syn_afterglow), 5 (ode_solver_test), time_dependence_test  (6)")
            exit(0)
        comp = compiler(rules='xTests', **self.comp_kw)
        if clean:
            comp.cleanup()
        comp.compile()
        outfile = self.flabel + '.jp.h5'
        run_cmd = '{0}Paramo {1} {2} {3}'.format(comp.compileDir,test_choice,self.par.params_file, outfile)
        print("\n--> Running:\n  ", run_cmd, "\n")
        os.system(run_cmd)
        print("\n--> Finished")

    # ----->   BlazMag compilation and run
    #cmd args
    #1. cool_withKN
    #2. with_abs
    def run_blazMag(self, cmd_args=(None, None), pream=None, clean=False, cl=False):
        if cmd_args[0] is None or cmd_args[0] is False:
            in1 = 'F'
        else:
            in1 = 'T'
        if cmd_args[1] is None or cmd_args[1] is False:
            in2 = 'F'
        else:
            in2 = 'T'
        comp = compiler(CONFIG=1, **self.comp_kw)
        if clean:
            comp.cleanup()
        comp.compile()
        outfile = self.flabel + '.jp.h5'
        if pream is None:
            run_cmd = '{0}Paramo {1} {2} {3} {4}'.format(comp.compileDir, self.par.params_file, outfile, in1, in2)
        else:
            run_cmd = '{0} {1}Paramo {2} {3} {4} {5}'.format(pream, comp.compileDir, self.par.params_file, outfile, in1, in2)
        print("\n--> Parameters:")
        os.system("cat -n " + self.par.params_file)
        if cl:
            print("\n--> Running:\n  ", run_cmd, "\n")
            return run_cmd
        else:
            print("\n--> Running:\n  ", run_cmd, "\n")
            os.system(run_cmd)
            print("\n--> Finished")

        # ----->   BlazMag compilation and run

    def run_interJets(self, pream=None, clean=False, cl=False):

        comp = compiler(CONFIG=5, **self.comp_kw)
        if clean:
            comp.cleanup()
        comp.compile()
        outfile = self.flabel + '.jp.h5'
        if pream is None:
            run_cmd = '{0}Paramo {1} {2}'.format(comp.compileDir, self.par.params_file, outfile)
        else:
            run_cmd = '{0} {1}Paramo {2} {3}'.format(pream, comp.compileDir, self.par.params_file, outfile)
        print("\n--> Parameters:")
        os.system("cat -n " + self.par.params_file)
        if cl:
            print("\n--> Running:\n  ", run_cmd, "\n")
            return run_cmd
        else:
            print("\n--> Running:\n  ", run_cmd, "\n")
            os.system(run_cmd)
            print("\n--> Finished")

    #####
    def run_Aglow(self, cmd_args=(False, False, False), pream=None, clean=False, cl=False, wMezcal=False):
        '''
        Description
        -----------
        Aglow compilation and run.

        Parameters
        ----------
        cmd_arg
        pream
        clean
        cl
        wMezcal
        '''
        if cmd_args[0] is None or cmd_args[0] is False:
            in1 = 'F'
        else:
            in1 = 'T'
        if cmd_args[1] is None or cmd_args[1] is False:
            in2 = 'F'
        else:
            in2 = 'T'
        if cmd_args[2] is None or cmd_args[2] is False:
            in3 = 'F'
        else:
            in3 = 'T'
        if wMezcal:
            comp = compiler(CONFIG=4, **self.comp_kw)
        else:
            comp = compiler(CONFIG=2, **self.comp_kw)
        if clean:
            comp.cleanup()
        comp.compile()
        outfile = self.flabel + '.jp.h5'
        if pream is None:
            run_cmd = '{0}Paramo {1} {2} {3} {4} {5}'.format(comp.compileDir, self.par.params_file, outfile, in1, in2, in3)
        else:
            run_cmd = '{0} {1}Paramo {2} {3} {4} {5} {6}'.format(pream, comp.compileDir, self.par.params_file, outfile, in1, in2, in3)
        print("\n--> Parameters:")
        os.system("cat -n " + self.par.params_file)
        if cl:
            print("\n--> Running:\n  ", run_cmd, "\n")
            return run_cmd
        else:
            print("\n--> Running:\n  ", run_cmd, "\n")
            os.system(run_cmd)
            print("\n--> Finished")

    # ----->   turbulence compilation and run
    def run_turb(self, cmd_args=(None, None), pream=None, clean=False, cl=False):
        if cmd_args[0] is None or cmd_args[0] is False:
            wCool = False
        else:
            wCool = True
        if cmd_args[1] is None or cmd_args[1] is False:
            wAbs = False
        else:
            wAbs = True

        comp = compiler(rules='xTurBlaz', **self.comp_kw)

        if clean:
            comp.cleanup()
        comp.compile()
        outfile = self.flabel + '.jp.h5'

        if pream is None:
            run_cmd = '{0}xTurBlaz {1} {2} {3} {4}'.format(comp.compileDir, self.par.params_file, wCool, wAbs, outfile)
        else:
            run_cmd = '{0} {1}xTurBlaz {2} {3} {4} {5}'.format(pream, comp.compileDir, self.par.params_file, wCool, wAbs, outfile)
        print("\n--> Parameters:")
        os.system("cat -n " + self.par.params_file)
        if cl:
            print("\n--> Running:\n  ", run_cmd, "\n")
            return run_cmd
        else:
            print("\n--> Running:\n  ", run_cmd, "\n")
            os.system(run_cmd)
            print("\n--> Finished")


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

    if (cores is None) or (qname == 'debug'):
        c = 24
    else:
        c = cores

    with open(sname, 'w') as f:
        print("#!/bin/bash -l\n", file=f)
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
        if (c > 1) or not (qname == 'debug'):
            print("export OMP_NUM_THREADS={0}".format(c), file=f)
        print("\n# RUN", file=f)
        print(xcmd, file=f)
    return sname


#
#  SLURM
#
def SlurmFile(jname, qname, xcmd, depen=None, nodes=None, cores=None, mail=None, htime=None, box=None):
    '''
    Description
    -----------
    This function generates the Slurm file to queue a simulation
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

    if (cores is None) or (qname == 'debug'):
        c = 24
    else:
        c = cores

    with open(sname, 'w') as f:
        print("#!/bin/sh -l\n", file=f)
        print("# FILENAME: {0}\n".format(sname), file=f)
        print("#SBATCH -A {0}".format(qname), file=f)
        print("#SBATCH -N {0:d}".format(n), file=f)
        print("#SBATCH -n {0:d}".format(c), file=f)
        print("#SBATCH --exclusive", file=f)
        print("#SBATCH -t {0}".format(t), file=f)
        print("#SBATCH --job-name={0}".format(jname), file=f)
        print("#SBATCH -o /scratch/brown/jruedabe/joboutput/{0}.out".format(jname), file=f)
        print("#SBATCH -e /scratch/brown/jruedabe/joboutput/{0}.err".format(jname), file=f)
        if depen is not None:
            print("#SBATCH --depend=afterok:{0}".format(depen), file=f)
        if mail is not None:
            print("#SBATCH --mail-user={0}".format(mail), file=f)
            print("#SBATCH --mail-type=FAIL", file=f)
        else:
            print("#SBATCH --mail-type=NONE", file=f)
        print("Working at: {0}".format(os.getcwd()))
        print("\ncd {0}".format(os.getcwd()), file=f)
        if (c > 1) or not (qname == 'debug'):
            print("export OMP_NUM_THREADS={0}".format(c), file=f)
        if box == 'brown':
            print("module load intel impi hdf5", file=f)
        print("\n# RUN", file=f)
        print(xcmd, file=f)
    return sname
