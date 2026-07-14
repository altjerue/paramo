# –¿Qué es lo que pasa, doña Eduviges?
# Ella sacudieo la cabeza como si despertara de un sueño
# –Es el caballo de Miguel Páramo, que galopa por el camino de la Media Luna.

import os
import subprocess
import sys
from time import strftime, localtime


def exp10(decimal):
    """Function that gets the exponent of a double number."""
    string = "{:.20e}".format(decimal)
    parts = string.split("e")
    return int(parts[-1])


# ------------------------------------------------------------------------------
# Fortran double precision functions
# ------------------------------------------------------------------------------
def fortran_double(number, dble=True):
    """Function that returns a floating point (as a string) in the FORTRAN
    notation
    """
    string = f"{number:.20e}"
    parts = string.split("e")

    ss = parts[0].split(".")
    s = ss[1]
    while s.endswith("0"):
        s = s[:-1]
    if len(s) > 0:
        ss = ss[0] + "." + s
    else:
        ss = ss[0]
    ee = parts[-1]
    if ee.startswith("-"):
        e = ee.split("-")[-1]
        while e.startswith("0"):
            e = e[1:]
        ee = "-" + e
    else:
        while ee.startswith("+") or ee.startswith("0"):
            ee = ee[1:]
    if len(ee) == 0:
        ee = "0"
    if dble:
        return f"{ss}d{ee}"
    else:
        return f"{ss}e{ee}"


# ------------------------------------------------------------------------------
# Parameters class
# ------------------------------------------------------------------------------
class parameters(object):
    '''This is the parameters class
    '''

    # -----  PARAMETERS  -----
    def rParams(self):
        self.R = 1e15                 # radius of emitting region (assuming spherical)
        self.R0 = 1e14                # distance from central engine
        self.dLum = 4.0793e26         # luminosity distance (default Mrk 421)
        self.z = 0.03                 # redshift (default Mrk 421)
        self.theta_obs = 5.0          # observer viewing angle
        self.gamma_bulk = 1e2         # emitting region bulk Lorentz factor
        self.mu_mag = 1.0             # (1 + sigma) Gamma_bulk
        self.sigma = 1.0              # magnetization (sigma)
        self.f_rec = 1.0              # magnetic reconnection dissipative efficiency
        self.b_index = 0.0            # magnetic field decay index
        self.Bfield = 1.0              # magnetic field decay index
        self.eps_B = 0.03             # epsilon_B
        self.eps_e = 0.1              # epsilon_e
        self.eps_acc = 1e-2           # epsilon_acc
        self.theta_e = 10.0           # electrons temperature
        self.zeta_e = 0.99            # fraction of non-thermal particles
        self.f_esc = 1.0              # electrons escape time factor
        self.tstep = 1e-2             # time step factor
        self.tmax = 1e5               # maximum time
        self.tmin = 0e0               # minimum time
        self.tvar = 2e0               # variability time scale
        self.L_jet = 1e45             # jet luminosity
        self.E0 = 1e50                # energy of the blast wave
        self.n_ext = 1.0              # number density of the external medium
        self.g1 = 1e2                 # power-law min Lorentz factor
        self.g2 = 1e4                 # power-law max Lorentz factor
        self.gmin = 1.01              # EED minimum Lorentz factor
        self.gmax = 2e4               # EED maximum Lorentz factor
        self.pind = 2.5               # EED power-law index
        self.nu_ext = 1e14            # external radiation field frequency
        self.u_ext = 1e-4             # external radiation field energy density
        self.numin = 1e7              # minimum frequency
        self.numax = 1e15             # maximum frequency
        self.NG = 128                 # number of EED bins
        self.NT = 300                 # number of time steps
        self.NF = 256                 # number of frequencies
        self.time_grid = 1            # kind of cooling
        self.params_file = 'input.par' # name of the parameters file

    def __init__(self, **kwargs):
        self.rParams()
        self.__dict__.update(kwargs)

    def wParams(self):
        with open(self.params_file, 'w') as f:
            print(f"{fortran_double(self.R):<26s} ! Radius", file=f)
            print(f"{fortran_double(self.R0):<26s} ! Initial radius", file=f)
            print(f"{fortran_double(self.dLum):<26s} ! luminosity distance", file=f)
            print(f"{fortran_double(self.z):<26s} ! redshift", file=f)
            print(f"{fortran_double(self.theta_obs):<26s} ! viewing angle", file=f)
            print(f"{fortran_double(self.gamma_bulk):<26s} ! bulk Lorentz factor", file=f)
            print(f"{fortran_double(self.mu_mag):<26s} ! (1 + sigma) Gamma", file=f)
            print(f"{fortran_double(self.sigma):<26s} ! magnetization", file=f)
            print(f"{fortran_double(self.f_rec):<26s} ! dissipative efficiency of magnetic reconection", file=f)
            print(f"{fortran_double(self.b_index):<26s} ! magnetic field decay index", file=f)
            print(f"{fortran_double(self.Bfield):<26s} ! magnetic field strength", file=f)
            print(f"{fortran_double(self.theta_e):<26s} ! electrons temperature", file=f)
            print(f"{fortran_double(self.zeta_e):<26s} ! fraction of nonthermal electrons", file=f)
            print(f"{fortran_double(self.f_esc):<26s} ! electrons escape time", file=f)
            print(f"{fortran_double(self.tstep):<26s} ! time step factor", file=f)
            print(f"{fortran_double(self.tmax):<26s} ! maximum time", file=f)
            print(f"{fortran_double(self.tmin):<26s} ! minimum time", file=f)
            print(f"{fortran_double(self.tvar):<26s} ! variability time scale", file=f)
            print(f"{fortran_double(self.L_jet):<26s} ! jet luminosity", file=f)
            print(f"{fortran_double(self.E0):<26s} ! energy of the blast wave", file=f)
            print(f"{fortran_double(self.n_ext):<26s} ! number density of the external medium", file=f)
            print(f"{fortran_double(self.eps_e):<26s} ! epsilon_e", file=f)
            print(f"{fortran_double(self.eps_B):<26s} ! epsilon_B", file=f)
            print(f"{fortran_double(self.eps_acc):<26s} ! epsilon_acc", file=f)
            print(f"{fortran_double(self.g1):<26s} ! power-law min Lorentz factor", file=f)
            print(f"{fortran_double(self.g2):<26s} ! power-law max Lorentz factor", file=f)
            print(f"{fortran_double(self.gmin):<26s} ! EED min Lorentz factor", file=f)
            print(f"{fortran_double(self.gmax):<26s} ! EED max Lorentz factor", file=f)
            print(f"{fortran_double(self.pind):<26s} ! EED power-law index", file=f)
            print(f"{fortran_double(self.nu_ext):<26s} ! external rad. field frequency", file=f)
            print(f"{fortran_double(self.u_ext):<26s} ! external rad. field ener. density", file=f)
            print(f"{fortran_double(self.numin):<26s} ! min frequency", file=f)
            print(f"{fortran_double(self.numax):<26s} ! max frequency", file=f)
            print(f"{self.NG:<26d} ! number of EED bins", file=f)
            print(f"{self.NT:<26d} ! number of time steps", file=f)
            print(f"{self.NF:<26d} ! number of frequencies", file=f)
            print(f"{self.time_grid:<26d} ! kind of time grid", file=f)
            # print(f"{self.file_label:<26s} ! label to identify each output", file=f)
        print(f"--> Parameters file: {self.params_file}")

# ------------------------------------------------------------------------------
# Compiler class
# ------------------------------------------------------------------------------
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
        self.CONFIG = 0         # 0 (tests), 1 (blazars), 2 (afterflow), 3 (turbulence), 4 (Mezcal)
        self.server = 0         # 0 (UNIX PC) 1 (Brown@Purdue) 2 (RC@RIT)
        self.compileDir = './'  # the path to Paramo... must end with '/'

    def __init__(self, **kwargs):
        self.flags()
        self.__dict__.update(kwargs)
        self.cwd = os.getcwd()

    def compile(self):
        print(strftime("\n\n%a, %d %b %Y %H:%M:%S %Z", localtime()))
        make = ["make", "Paramo", "CONFIG="+str(self.problem), "COMPILER="+str(self.COMP)]
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


# ------------------------------------------------------------------------------
# Runner class
# ------------------------------------------------------------------------------
class Runner(object):
    '''
    Description
    -----------
    This class sets up the executable instructions.
    '''

    def __init__(self, flabel='DriverTest', par_kw={}, comp_kw={}):
        self.par = parameters(**par_kw)
        self.par.wParams()
        self.cwd = os.getcwd()
        self.comp_kw = comp_kw
        self.flabel = flabel  # a label to identify each output

    def run_test(self, clean=False):
        comp = compiler(rules='xTests', **self.comp_kw)
        if clean:
            comp.cleanup()
        comp.compile()
        run_cmd = '{0}xTests'.format(comp.compileDir)
        print("\n--> Running:\n  ", run_cmd, "\n")
        os.system(run_cmd)
        print("\n--> Finished")

    # ----->   BlazMag compilation and run
    def run_blazMag(self, cmd_args=(None, None), pream=None, clean=False, cl=False):
        if cmd_args[0] is None or cmd_args[0] is False:
            in1 = 'F'
        else:
            in1 = 'T'
        if cmd_args[1] is None or cmd_args[1] is False:
            in2 = 'F'
        else:
            in2 = 'T'
        comp = compiler(problem=1, **self.comp_kw)
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
            comp = compiler(problem=4, **self.comp_kw)
        else:
            comp = compiler(problem=2, **self.comp_kw)
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


# ------------------------------------------------------------------------------
# PBS
# ------------------------------------------------------------------------------
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


# ------------------------------------------------------------------------------
#  SLURM
# ------------------------------------------------------------------------------
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
