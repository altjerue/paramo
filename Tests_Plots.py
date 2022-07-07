import matplotlib.pyplot as plt
import numpy as np
import Arriero as ar
import Eduviges.extractor as extr
import matplotlib.cm as cm
import matplotlib

class SS_results:

    def __init__(self, outfile):
        ex = extr.fromHDF5(h5fname=outfile)
        self.numt = ex.hdf5ExtractScalar('numt', group='Parameters')
        self.numf = ex.hdf5ExtractScalar('numf', group='Parameters')
        self.numg = ex.hdf5ExtractScalar('numg', group='Parameters')
        self.tmax = ex.hdf5ExtractScalar('t_max', group='Parameters')
        self.tstep = ex.hdf5ExtractScalar('tstep', group='Parameters')
        self.gmin = ex.hdf5ExtractScalar('gamma_min', group='Parameters')
        self.gmax = ex.hdf5ExtractScalar('gamma_max', group='Parameters')
        self.g1 = ex.hdf5ExtractScalar('gamma_1', group='Parameters')
        self.g2 = ex.hdf5ExtractScalar('gamma_2', group='Parameters')
        self.qind = ex.hdf5ExtractScalar('pwl-index', group='Parameters')
        self.numin = ex.hdf5ExtractScalar('nu_min', group='Parameters')
        self.numax = ex.hdf5ExtractScalar('nu_max', group='Parameters')
        self.R = ex.hdf5ExtractScalar('emission_R', group='Parameters')
        self.B = ex.hdf5ExtractScalar('B_0', group='Parameters')
        self.uB = ex.hdf5ExtractScalar('uB', group='Parameters')
        self.C0 = ex.hdf5ExtractScalar('C0', group='Parameters')
        self.tacc = ex.hdf5ExtractScalar('tacc', group='Parameters')
        self.tesc = ex.hdf5ExtractScalar('tesc', group='Parameters')
        self.D0 = ex.hdf5Extract1D('D0', group='Parameters')

        self.nu = ex.hdf5Extract1D('freqs', group='Numeric')
        self.g = ex.hdf5Extract1D('gamma', group='Numeric')
        self.dg = ex.hdf5Extract1D('dg', group='Numeric')
        self.t = ex.hdf5Extract1D('time', group='Numeric')
        self.dt = ex.hdf5Extract1D('dt', group='Numeric')
        self.Inu1 = ex.hdf5Extract1D('Inu1', group='Numeric')
        self.Inu4 = ex.hdf5Extract1D('Inu4', group='Numeric')
        self.Inu5 = ex.hdf5Extract1D('Inu5', group='Numeric')
        self.Inu6 = ex.hdf5Extract1D('Inu6', group='Numeric')
        self.Ntot1 = ex.hdf5Extract1D('Ntot1', group='Numeric')
        self.Ntot2 = ex.hdf5Extract1D('Ntot2', group='Numeric')
        self.Ntot3 = ex.hdf5Extract1D('Ntot3', group='Numeric')
        self.Ntot4 = ex.hdf5Extract1D('Ntot4', group='Numeric')
        self.Ntot5 = ex.hdf5Extract1D('Ntot5', group='Numeric')
        self.Ntot6 = ex.hdf5Extract1D('Ntot6', group='Numeric')
        self.n1 = ex.hdf5Extract2D('n1', group='Numeric')
        self.n2 = ex.hdf5Extract2D('n2', group='Numeric')
        self.n3 = ex.hdf5Extract2D('n3', group='Numeric')
        self.n4 = ex.hdf5Extract2D('n4', group='Numeric')
        self.n5 = ex.hdf5Extract2D('n5', group='Numeric')
        self.n6 = ex.hdf5Extract2D('n6', group='Numeric')
        self.jmbs1 = ex.hdf5Extract2D('jmbs1', group='Numeric')
        self.jmbs4 = ex.hdf5Extract2D('jmbs4', group='Numeric')
        self.jmbs5 = ex.hdf5Extract2D('jmbs5', group='Numeric')
        self.jmbs6 = ex.hdf5Extract2D('jmbs6', group='Numeric')
        self.jssc1 = ex.hdf5Extract2D('jssc1', group='Numeric')
        self.jssc4 = ex.hdf5Extract2D('jssc4', group='Numeric')
        self.jssc5 = ex.hdf5Extract2D('jssc5', group='Numeric')
        self.jssc6 = ex.hdf5Extract2D('jssc6', group='Numeric')
        self.ambs1 = ex.hdf5Extract2D('ambs1', group='Numeric')
        self.ambs4 = ex.hdf5Extract2D('ambs4', group='Numeric')
        self.ambs5 = ex.hdf5Extract2D('ambs5', group='Numeric')
        self.ambs6 = ex.hdf5Extract2D('ambs6', group='Numeric')
        self.ns = [self.n1,self.n2,self.n3,self.n4,self.n5,self.n6]
        self.jmbss=[self.jmbs1,self.jmbs4,self.jmbs5,self.jmbs6]
        self.ambss = [self.ambs1,self.ambs4,self.ambs5,self.ambs6]
        self.jsscs = [self.jssc1,self.jssc4,self.jssc5,self.jssc6]

outfile = './steady_state_test.h5'

plots_folder = '/home/zach/Documents/Code_Projects/paramo/Plots/'

def run_steady_state_test():
    rr = ar.Runner(comp_kw={'OMP': True, 'HDF5': True})
    rr.run_test(output_file=outfile, clean=False, test_choice=1)


def get_steady_state_results():
    ssr = SS_results(outfile=outfile)
    return ssr


def ssr_n_plots():
    ssr = get_steady_state_results()

    for j in ssr.ns:
        fig, ax = plt.subplots()
        ax.set_yscale('log')
        ax.set_xscale('log')



        g = ssr.g
        t = ssr.t
        n = j

        cmap = cm.rainbow
        sm = plt.cm.ScalarMappable(cmap=cmap,  norm=matplotlib.colors.LogNorm(vmin=t[0], vmax=t[-1]))

        pls = []
        my = None

        for i in range(len(t)):
            if(i%25 !=0 and i != len(t) - 1 ):
                continue
            y = n[:,i]
            if(my == None):
                my = max(y)
            if(max(y)> my):
                my= max(y)
            pl = ax.plot(g,y,color=cmap(i))
            pls.append(pl)

        ax.set_ylim(1e-4, 2*my)
        ax.set_xlim(1e0, 2e6)

        cbticks = []
        for i in range(8):
            cbticks.append(r"$10^{{{0}}}$".format(i))
        cbar = plt.colorbar(sm,anchor=(-0.6,0.0),ticks=np.logspace(0,7,8))
        cbar.ax.set_yticklabels(cbticks)
        cbar.ax.minorticks_off()
        cbar.ax.set_ylabel(r"t [s]",fontsize=18)
        cbar.ax.tick_params(
            labelsize=15
        )

        plt.tick_params(
            axis='x',
            which='minor',
            bottom = False,
            top = False,
            labelbottom = False
        )

        plt.tick_params(
            axis='both',
            which='both',
            labelsize=15
        )
        plt.tick_params(
            axis='both',
            which='major',
            size=10
        )
        plt.tick_params(
            axis='y',
            which='minor',
            size=5
        )

        ax.set_xlabel(r"$\gamma$",fontsize=18)
        ax.set_ylabel(r"n $[cm^{-3}]$",fontsize=18)

        plt.tight_layout()

        # plt.savefig(plots_folder+"n1vsg.png")
        plt.show()

def ssr_j_plots():
    ssr = get_steady_state_results()

    for j in range(len(ssr.jmbss)):
        fig, ax = plt.subplots()
        ax.set_yscale('log')
        ax.set_xscale('log')



        nu = ssr.nu
        t = ssr.t
        jm = ssr.jmbss[j]
        jssc = ssr.jsscs[j]

        cmap = cm.rainbow
        sm = plt.cm.ScalarMappable(cmap=cmap,  norm=matplotlib.colors.LogNorm(vmin=t[0], vmax=t[-1]))

        pls = []
        my = None

        for i in range(len(t)):
            if(i%25 !=0 and i != len(t) - 1 ):
                continue
            y = jm[i,:] + jssc[i,:]
            if(my == None):
                my = max(y)
            if(max(y)> my):
                my= max(y)
            pl = ax.plot(nu,y,color=cmap(i))
            pls.append(pl)

        ax.set_ylim(my/1e10, 2*my)
        ax.set_xlim(ssr.numin, ssr.numax)

        cbticks = []
        for i in range(8):
            cbticks.append(r"$10^{{{0}}}$".format(i))
        cbar = plt.colorbar(sm,anchor=(-0.6,0.0),ticks=np.logspace(0,7,8))
        cbar.ax.set_yticklabels(cbticks)
        cbar.ax.minorticks_off()
        cbar.ax.set_ylabel(r"t [s]",fontsize=18)
        cbar.ax.tick_params(
            labelsize=15
        )

        plt.tick_params(
            axis='x',
            which='minor',
            bottom = False,
            top = False,
            labelbottom = False
        )

        plt.tick_params(
            axis='both',
            which='both',
            labelsize=15
        )
        plt.tick_params(
            axis='both',
            which='major',
            size=10
        )
        plt.tick_params(
            axis='y',
            which='minor',
            size=5
        )

        ax.set_xlabel(r"$\nu$ [Hz]",fontsize=18)
        ax.set_ylabel(r"j $[\frac{erg}{s cm^{3}}]$",fontsize=18)

        plt.tight_layout()

        # plt.savefig(plots_folder+"n1vsg.png")
        plt.show()

run_steady_state_test()
ssr_n_plots()
# ssr_j_plots()
