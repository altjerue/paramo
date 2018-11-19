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
        self.R = 1e18                   # radius of emitting region (assuming spherical)
        self.R0 = 1e15                  # radius of emitting region (assuming spherical)
        self.Rinit = 1e15               # radius of emitting region (assuming spherical)
        self.dLum = 4.0793e26           # luminosity distance (default Mrk 421)
        self.z = 0.03                   # redshift (default Mrk 421)
        self.gamma_bulk = 1e2           # emitting region bulk Lorentz factor
        self.theta_obs = 5.0            # observer viewing angle
        self.B = 1.0                    # magnetic field magnitude
        self.b_index = 0.0              # magnetic field decay index
        self.theta_e = 10.0             # electrons temperature
        self.dtacc = 1e2                # injection period
        self.tstep = 1e-2               # time step factor
        self.tmax = 1e5                 # maximum time
        self.Q0 = 1.0                   # num. dens. of particles injected per second
        self.g1 = 1e2                   # power-law min Lorentz factor
        self.g2 = 1e4                   # power-law max Lorentz factor
        self.gmin = 1.01                # EED minimum Lorentz factor
        self.gmax = 2e4                 # EED maximum Lorentz factor
        self.qind = 2.5                 # EED power-law index
        self.numin = 1e7                # minimum frequency
        self.numax = 1e15               # maximum frequency
        self.numbins = 128              # number of EED bins
        self.numdt = 300                # number of time steps
        self.numdf = 256                # number of frequencies
        self.cool_kind = 1              # kind of cooling
        self.time_grid = 1              # kind of cooling

        self.params_file = 'input.par'  # name of the parameters file

    def __init__(self, **kwargs):
        self.rParams()
        self.__dict__.update(kwargs)

    def wParams(self):
        with open(self.params_file, 'w') as f:
            print(fortran_double(self.R), '! Radius', file=f)
            print(fortran_double(self.R0), '! Normalization radius', file=f)
            print(fortran_double(self.Rinit), '! Initial radius', file=f)
            print(fortran_double(self.dLum), '! luminosity distance', file=f)
            print(fortran_double(self.z), '! redshift', file=f)
            print(fortran_double(self.gamma_bulk), '! bulk Lorentz factor', file=f)
            print(fortran_double(self.theta_obs), '! viewing angle', file=f)
            print(fortran_double(self.B), '! magnetic field', file=f)
            print(fortran_double(self.b_index), '! magnetic field decay index', file=f)
            print(fortran_double(self.theta_e), '! electrons temperature', file=f)
            print(fortran_double(self.dtacc), '! injection period',  file=f)
            print(fortran_double(self.tstep), '! time step factor', file=f)
            print(fortran_double(self.tmax), '! maximum time', file=f)
            print(fortran_double(self.Q0), '! num. dens. of particles injected', file=f)
            print(fortran_double(self.g1), '! power-law min Lorentz factor', file=f)
            print(fortran_double(self.g2), '! power-law max Lorentz factor', file=f)
            print(fortran_double(self.gmin), '! EED min Lorentz factor', file=f)
            print(fortran_double(self.gmax), '! EED max Lorentz factor', file=f)
            print(fortran_double(self.qind), '! EED power-law index', file=f)
            print(fortran_double(self.numin), '! min frequency', file=f)
            print(fortran_double(self.numax), '! max frequency', file=f)
            print(self.numbins, '! number of EED bins', file=f)
            print(self.numdt, '! number of time steps', file=f)
            print(self.numdf, '! number of frequencies', file=f)
            print(self.cool_kind, '! kind of cooling', file=f)
            print(self.time_grid, '! kind of time grid', file=f)
            # print(self.file_label, '! label to identify each output', file=f)


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
        print("--> Running Makefile:\n   ", make, "\n")
        log = strftime("%a, %d %b %Y %H:%M:%S %Z", localtime())
        os.system(make)
        with open("make.log", "a") as logfile:
            logfile.write(log + "\n" + make + "\n\n")
        logfile.close()
        os.chdir(self.cwd)

    def cleanup(self):
        os.chdir(self.compile_dir)
        os.system("make clean")
        os.chdir(self.cwd)

#
#  #####  #    # #    #
#  #    # #    # ##   #
#  #    # #    # # #  #
#  #####  #    # #  # #
#  #   #  #    # #   ##
#  #    #  ####  #    #


class Paramo(object):

    def __init__(self, wCool=False, wAbs=False, wSSC=False, flabel='DriverTest', par_kw={}, comp_kw={}):
        self.par = parameters(**par_kw)
        self.comp = compiler(rules='xParamo', **comp_kw)
        self.par.wParams()
        self.cwd = os.getcwd()
        # -----  ARGS OF THE EXECUTABLE  -----
        self.wCool = wCool    # variable cooling
        self.wAbs = wAbs      # compute MBS self-absorption
        self.wSSC = wSSC      # compute SSC emissivity
        self.flabel = flabel  # a label to identify each output

    def output_file(self):
        outf = ''
        argv = ''
        if self.comp.HYB:
            outf += 'H'
        else:
            outf += 'P'

        if self.comp.MBS:
            outf += 'M'
        else:
            outf += 'S'

        if self.wCool:
            outf += 'V'
            argv += ' T'
        else:
            outf += 'C'
            argv += ' F'

        if self.wAbs:
            outf += 'O'
            argv += ' T'
        else:
            outf += 'T'
            argv += ' F'

        if self.wSSC:
            outf += 'wSSC'
            argv += ' T'
        else:
            outf += 'oSSC'
            argv += ' F'

        return outf + '-' + self.flabel + '.jp.h5', argv

    def run_Paramo(self, pream=None):
        self.comp.compile()
        self.outfile, self.argv = self.output_file()
        if pream is None:
            run_cmd = '{0}xParamo {1}{2} {3}'.format(self.comp.compile_dir, self.par.params_file, self.argv, self.outfile)
        else:
            run_cmd = '{0} {1}xParamo {2}{3} {4}'.format(pream, self.comp.compile_dir, self.par.params_file, self.argv, self.outfile)
        print("\n--> Running:\n  ", run_cmd, "\n")
        os.system(run_cmd)
        print("\n--> Paramo finished")
        return self.outfile


class ITobs(object):

    def __init__(self, paramo_file, comp_kw={}):
        self.comp = compiler(rules='xITobs', **comp_kw)
        self.Pfile = paramo_file
        self.cwd = os.getcwd()

    def run_ITobs(self):
        self.comp.compile()
        run_cmd = '{0}xITobs {1}'.format(self.comp.compile_dir, self.Pfile)
        print("\n--> Running:\n  ", run_cmd, "\n")
        os.system(run_cmd)
        print("\n--> Paramo finished")
