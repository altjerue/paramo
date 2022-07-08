import matplotlib.pyplot as plt
import numpy as np
import Arriero as ar
import Eduviges.extractor as extr
import matplotlib.cm as cm
import matplotlib

outfile = 'interJets_test'

plots_folder = '/home/zach/Documents/Code_Projects/paramo/Plots/'


class ij_results:

    def __init__(self, outfile):
        ex = extr.fromHDF5(h5fname=outfile)
        self.numt = ex.hdf5ExtractScalar('numt', group='Parameters')
        self.numf = ex.hdf5ExtractScalar('numdf', group='Parameters')
        self.numg = ex.hdf5ExtractScalar('numbins', group='Parameters')
        self.time_grid = ex.hdf5ExtractScalar('time-grid', group='Parameters')
        self.t_max = ex.hdf5ExtractScalar('t_max', group='Parameters')
        self.t_min = ex.hdf5ExtractScalar('t_min', group='Parameters')
        self.tstep = ex.hdf5ExtractScalar('tstep', group='Parameters')
        self.d_lum = ex.hdf5ExtractScalar('d_lum', group='Parameters')
        self.redshift = ex.hdf5ExtractScalar('redshift', group='Parameters')
        self.Gamma_bulk0 = ex.hdf5ExtractScalar('Gamma_bulk0', group='Parameters')
        self.view_angle = ex.hdf5ExtractScalar('view-angle', group='Parameters')
        self.gmin = ex.hdf5ExtractScalar('gamma_min', group='Parameters')
        self.gmax = ex.hdf5ExtractScalar('gamma_max', group='Parameters')
        self.g1 = ex.hdf5ExtractScalar('gamma_1', group='Parameters')
        self.g2 = ex.hdf5ExtractScalar('gamma_2', group='Parameters')
        self.qind = ex.hdf5ExtractScalar('pwl-index', group='Parameters')
        self.L_j = ex.hdf5ExtractScalar('L_j', group='Parameters')
        self.epsilon_e = ex.hdf5ExtractScalar('epsilon_e', group='Parameters')
        self.epsilon_B = ex.hdf5ExtractScalar('epsilon_B', group='Parameters')
        self.numin = ex.hdf5ExtractScalar('nu_min', group='Parameters')
        self.numax = ex.hdf5ExtractScalar('nu_max', group='Parameters')
        self.E0 = ex.hdf5ExtractScalar('E0', group='Parameters')
        self.Ejet = ex.hdf5ExtractScalar('Ejet', group='Parameters')
        self.R0 = ex.hdf5ExtractScalar('R0', group='Parameters')
        self.n_ext = ex.hdf5ExtractScalar('n_ext', group='Parameters')

        self.Rd = ex.hdf5ExtractScalar('Rd')
        self.td = ex.hdf5ExtractScalar('td')
        self.t = ex.hdf5Extract1D('time')
        self.t_obs = ex.hdf5Extract1D('t_obs')
        self.Rb = ex.hdf5Extract1D('Rb')
        self.Rbw = ex.hdf5Extract1D('Rbw')
        self.volume = ex.hdf5Extract1D('volume')
        self.Gamma_bulk = ex.hdf5Extract1D('Gamma_bulk')
        self.Doppler = ex.hdf5Extract1D('Doppler')
        self.nu = ex.hdf5Extract1D('nu')
        self.nu_obs = ex.hdf5Extract1D('nu_obs')
        self.g = ex.hdf5Extract1D('gamma')
        self.jnut = ex.hdf5Extract2D('jnut')
        self.jsyn = ex.hdf5Extract2D('jmbs')
        self.jssc = ex.hdf5Extract2D('jssc')
        self.jeic = ex.hdf5Extract2D('jeic')
        self.anut = ex.hdf5Extract2D('anut')
        self.ambs = ex.hdf5Extract2D('ambs')
        self.Qinj = ex.hdf5Extract2D('Qinj')
        self.n = ex.hdf5Extract2D('n_e')
        self.gdotty = ex.hdf5Extract2D('cool-coef')
        self.Ddiff = ex.hdf5Extract2D('diffusion')





def run_interJets():
    rr = ar.Runner(flabel=outfile,comp_kw={'OMP': True, 'HDF5': True,'compileDir':'../'})
    ##adjust parameters
    rr.par.numax = 1e28
    rr.par.wParams()
    ###
    rr.run_interJets(clean=True)

def get_ij_results():
    ijr = ij_results(outfile=outfile+'.jp.h5')
    return ijr

def ijr_n_plot():
    ijr = get_ij_results()

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')



    g = ijr.g
    t = ijr.t
    n = ijr.n

    cmap = cm.rainbow
    sm = plt.cm.ScalarMappable(cmap=cmap,  norm=matplotlib.colors.LogNorm(vmin=t[0], vmax=t[-1]))

    pls = []
    my = None

    for i in range(len(t)):
        if(i%25 !=0 and i != len(t) - 1 ):
            continue
        y = n[i,:]
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


def ijr_j_plots():
    ijr = get_ij_results()

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')



    nu = ijr.nu
    t = ijr.t
    jm = ijr.jsyn
    jssc = ijr.jssc
    jeic = ijr.jeic

    cmap = cm.rainbow
    sm = plt.cm.ScalarMappable(cmap=cmap,  norm=matplotlib.colors.LogNorm(vmin=t[0], vmax=t[-1]))

    pls = []
    my = None

    for i in range(len(t)):
        if(i%25 !=0 and i != len(t) - 1 ):
            continue
        y = jm[i,:] + jssc[i,:] + jeic[i,:]
        if(my == None):
            my = max(y)
        if(max(y)> my):
            my= max(y)
        pl = ax.plot(nu,y,color=cmap(i))
        pls.append(pl)

    ax.set_ylim(my/1e10, 2*my)
    ax.set_xlim(ijr.numin, ijr.numax)

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

run_interJets()
# get_ij_results()
ijr_n_plot()
ijr_j_plots()