import os
from time import strftime, localtime
from SAPyto.misc import fortran_double


#
#  #####    ##   #####    ##   #    #  ####
#  #    #  #  #  #    #  #  #  ##  ## #
#  #    # #    # #    # #    # # ## #  ####
#  #####  ###### #####  ###### #    #      #
#  #      #    # #   #  #    # #    # #    #
#  #      #    # #    # #    # #    #  ####
class parameters(object):

    def rParams(self):
        # -----  PARAMETERS  -----
        self.radius = 1e16                   # radius of emitting region (assuming spherical)
        self.pos_init = 1e15            # Distance from central engine
        self.dLum = 4.0793e26           # luminosity distance (default Mrk 421)
        self.z = 0.03                   # redshift (default Mrk 421)
        self.gamma_bulk = 1e2           # emitting region bulk Lorentz factor
        self.theta_obs = 5.0            # observer viewing angle
        self.sigma = 1.0                # magnetization (sigma)
        self.b_index = 0.0              # magnetic field decay index
        self.theta_e = 10.0             # electrons temperature
        self.zeta_e = 0.99
        self.tstep = 1e-2               # time step factor
        self.tmax = 1e5                 # maximum time
        self.L_j = 1e45                # num. dens. of particles injected per second
        self.eps_e = 0.1                # epsilon_e
        self.eps_B = 0.03                # epsilon_B
        self.g1 = 1e2                   # power-law min Lorentz factor
        self.g2 = 1e4                   # power-law max Lorentz factor
        self.gmin = 1.01                # EED minimum Lorentz factor
        self.gmax = 2e4                 # EED maximum Lorentz factor
        self.qind = 2.5                 # EED power-law index
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
            print(fortran_double(self.radius), ' ! Radius', file=f)
            print(fortran_double(self.pos_init), ' ! Initial radius', file=f)
            print(fortran_double(self.dLum), ' ! luminosity distance', file=f)
            print(fortran_double(self.z), ' ! redshift', file=f)
            print(fortran_double(self.gamma_bulk), ' ! bulk Lorentz factor', file=f)
            print(fortran_double(self.theta_obs), ' ! viewing angle', file=f)
            print(fortran_double(self.sigma), ' ! magnetization', file=f)
            print(fortran_double(self.b_index), ' ! magnetic field decay index', file=f)
            print(fortran_double(self.theta_e), ' ! electrons temperature', file=f)
            print(fortran_double(self.zeta_e), ' ! fraction of nonthermal electrons', file=f)
            print(fortran_double(self.tstep), ' ! time step factor', file=f)
            print(fortran_double(self.tmax), ' ! maximum time', file=f)
            print(fortran_double(self.L_j), ' ! injection luminosity', file=f)
            print(fortran_double(self.eps_e), ' ! epsilon_e', file=f)
            print(fortran_double(self.eps_B), ' ! epsilon_B', file=f)
            print(fortran_double(self.g1), ' ! power-law min Lorentz factor', file=f)
            print(fortran_double(self.g2), ' ! power-law max Lorentz factor', file=f)
            print(fortran_double(self.gmin), ' ! EED min Lorentz factor', file=f)
            print(fortran_double(self.gmax), ' ! EED max Lorentz factor', file=f)
            print(fortran_double(self.qind), ' ! EED power-law index', file=f)
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

#
#   ####   ####  #    # #####  # #      ######
#  #    # #    # ##  ## #    # # #      #
#  #      #    # # ## # #    # # #      #####
#  #      #    # #    # #####  # #      #
#  #    # #    # #    # #      # #      #
#   ####   ####  #    # #      # ###### ######


class compiler(object):

    # -----  COMPILER FLAGS & RULES -----
    def flags(self):
        self.HYB = False         # compile with HYB=1 flag
        self.MBS = False         # compile with MBS=1 flag
        self.arch = ''           # compile with specific arch flag
        self.OMP = False         # compile with OpenMP
        self.DBG = False         # compile for debugging
        self.rules = 'all'       # rule to compile
        self.compile_dir = './'  # address to Paramo, must end with '/'

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

        if self.HYB:
            make += ' HYB=1'

        if self.MBS:
            make += ' MBS=1'

        os.chdir(self.compile_dir)
        print("--> Running Makefile:\n ", make, "\n")
        log = strftime("%a, %d %b %Y %H:%M:%S %Z", localtime())
        os.system(make)
        with open("make.log", "a") as logfile:
            logfile.write(log + "\n" + make + "\n\n")
        logfile.close()
        os.chdir(self.cwd)

    def cleanup(self):
        os.chdir(self.compile_dir)
        os.system("make clean_all")
        os.chdir(self.cwd)


#
#  #####  #####   ####     ###### # #      ######
#  #    # #    # #         #      # #      #
#  #    # #####   ####     #####  # #      #####
#  #####  #    #      #    #      # #      #
#  #      #    # #    #    #      # #      #
#  #      #####   ####     #      # ###### ######
def PBSfile(jname, qname, xname, args=None, pream=None, depen=None, nodes=None, cores=None, mail=None, htime=None):
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
    if cores is None:
        c = 1
    else:
        c = cores
    with open(sname, 'w') as f:
        print("#!/bin/sh -l\n", file=f)
        print("# FILENAME: {0}\n".format(sname), file=f)
        print("#PBS -q " + qname, file=f)
        print("#PBS -l nodes={0:d}:ppn={1:d},nacesspolicy=singleuser".format(n, c), file=f)
        print("#PBS -l walltime={0}".format(t), file=f)
        print("#PBS -N " + jname, file=f)
        print("#PBS -o /home/jruedabe/joboutput/{0}.out".format(jname), file=f)
        print("#PBS -e /home/jruedabe/joboutput/{0}.err".format(jname), file=f)
        if depen is not None:
            print("#PBS -W depend=afterok:{0}".format(depen), file=f)
        if mail is not None:
            print("#PBS -M {0}".format(mail), file=f)
            print("#PBS -m bae", file=f)
        print("\n# Run command:", file=f)
        if pream is not None:
            xname = "{0} {1}".format(pream, xname)

        if args is None:
            print(xname, file=f)
        else:
            print("{0} {1}".format(xname, " ".join(args)), file=f)
    return sname


#
#  #####  #    # #    #
#  #    # #    # ##   #
#  #    # #    # # #  #
#  #####  #    # #  # #
#  #   #  #    # #   ##
#  #    #  ####  #    #
class Paramo(object):

    def __init__(self, wCool=False, wAbs=False, wSSC=False, Hyb=False, flabel='DriverTest', par_kw={}, comp_kw={}):
        self.par = parameters(**par_kw)
        self.comp = compiler(rules='xParamo', **comp_kw)
        self.par.wParams()
        self.cwd = os.getcwd()
        # -----  ARGS OF THE EXECUTABLE  -----
        self.wCool = wCool    # variable cooling
        self.wAbs = wAbs      # compute MBS self-absorption
        self.wSSC = wSSC      # compute SSC emissivity
        self.Hyb = Hyb
        self.flabel = flabel  # a label to identify each output

    def output_file(self):
        outf = ''
        argv = []

        if self.Hyb:
            outf += 'H'
            argv.append('T')
        else:
            outf += 'P'
            argv.append('F')

        if self.wCool:
            outf += 'V'
            argv.append('T')
        else:
            outf += 'C'
            argv.append('F')

        if self.wAbs:
            outf += 'O'
            argv.append('T')
        else:
            outf += 'T'
            argv.append('F')

        if self.wSSC:
            outf += 'wS'
            argv.append('T')
        else:
            outf += 'oS'
            argv.append('F')

        return outf + '-' + self.flabel + '.jp.h5', argv

    def run_Paramo(self, pream=None):
        self.comp.cleanup()
        self.comp.compile()
        outfile, argv = self.output_file()
        if pream is None:
            run_cmd = '{0}xParamo {1} {2} {3}'.format(self.comp.compile_dir, self.par.params_file, " ".join(argv), outfile)
        else:
            run_cmd = '{0} {1}xParamo {2} {3} {4}'.format(pream, self.comp.compile_dir, self.par.params_file, " ".join(argv), outfile)
        print("\n--> Parameters:")
        os.system("cat -n " + self.par.params_file)
        print("\n--> Running:\n  ", run_cmd, "\n")
        os.system(run_cmd)
        print("\n--> Paramo finished")
        return outfile

    def cl_Paramo(self, pream=None):
        self.comp.compile()
        self.outfile, self.argv = self.output_file()
        if pream is None:
            run_cmd = '{0}xParamo'.format(self.comp.compile_dir)
        else:
            run_cmd = '{0} {1}xParamo'.format(pream, self.comp.compile_dir)
        args = []
        args.append(self.par.params_file)
        args += self.argv
        args.append(self.outfile)
        return run_cmd, args


class ITobs(object):

    def __init__(self, paramo_file, comp_kw={}):
        self.comp = compiler(rules='xITobs', **comp_kw)
        self.Pfile = paramo_file
        self.cwd = os.getcwd()

    def run_ITobs(self, pream=None):
        self.comp.compile()
        if pream is None:
            run_cmd = '{0}xITobs {1}'.format(self.comp.compile_dir, self.Pfile)
        else:
            run_cmd = '{0} {1}xITobs {2}'.format(pream, self.comp.compile_dir, self.Pfile)
        print("\n--> Running:\n  ", run_cmd, "\n")
        os.system(run_cmd)
        print("\n--> Paramo finished")

    def cl_ITobs(self, pream=None):
        self.comp.compile()
        if pream is None:
            run_cmd = '{0}xITobs'.format(self.comp.compile_dir)
        else:
            run_cmd = '{0} {1}xITobs'.format(pream, self.comp.compile_dir)
        args = []
        args.append(self.Pfile)
        return run_cmd, args
