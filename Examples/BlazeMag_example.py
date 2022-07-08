import matplotlib.pyplot as plt
import numpy as np
import Arriero as ar
import Eduviges.extractor as extr
import matplotlib.cm as cm
import matplotlib

outfile = 'blazeMag_test'
outdir='./'
plots_folder = '/home/zach/Documents/Code_Projects/paramo/Plots/'


class bm_results:

    def __init__(self, outfile):
        ex = extr.fromHDF5(h5fname=outdir+outfile)
        self.numt = ex.hdf5ExtractScalar('numdt', group='Parameters')
        self.numf = ex.hdf5ExtractScalar('numdf', group='Parameters')
        self.numg = ex.hdf5ExtractScalar('numbins', group='Parameters')
        self.t_max = ex.hdf5ExtractScalar('t_max', group='Parameters')
        self.tstep = ex.hdf5ExtractScalar('tstep', group='Parameters')
        self.R_b = ex.hdf5ExtractScalar('R_b', group='Parameters')
        self.R_em = ex.hdf5ExtractScalar('R_em', group='Parameters')
        self.d_lum = ex.hdf5ExtractScalar('d_lum', group='Parameters')
        self.redshift = ex.hdf5ExtractScalar('redshift', group='Parameters')
        self.Gamma_bulk = ex.hdf5ExtractScalar('Gamma_bulk', group='Parameters')
        self.sigma = ex.hdf5ExtractScalar('sigma', group='Parameters')
        self.theta_obs_deg = ex.hdf5ExtractScalar('theta_obs_deg', group='Parameters')
        self.gmin = ex.hdf5ExtractScalar('gamma_min', group='Parameters')
        self.gmax = ex.hdf5ExtractScalar('gamma_max', group='Parameters')
        self.g1 = ex.hdf5ExtractScalar('gamma_1', group='Parameters')
        self.g2 = ex.hdf5ExtractScalar('gamma_2', group='Parameters')
        self.qind = ex.hdf5ExtractScalar('pwl-index', group='Parameters')
        self.u_ext = ex.hdf5ExtractScalar('u_ext', group='Parameters')
        self.nu_ext = ex.hdf5ExtractScalar('nu_ext', group='Parameters')
        self.L_j = ex.hdf5ExtractScalar('L_jet', group='Parameters')
        self.numin = ex.hdf5ExtractScalar('nu_min', group='Parameters')
        self.numax = ex.hdf5ExtractScalar('nu_max', group='Parameters')
        self.numag = ex.hdf5ExtractScalar('mu_mag', group='Parameters')

        self.t_inj = ex.hdf5ExtractScalar('t_inj')
        self.t_esc = ex.hdf5ExtractScalar('t_esc')

        self.Bfield = ex.hdf5ExtractScalar('Bfield', group='electrons-energy')
        self.L_B = ex.hdf5ExtractScalar('L_B', group='electrons-energy')
        self.uB = ex.hdf5ExtractScalar('uB', group='electrons-energy')
        self.L_e = ex.hdf5ExtractScalar('L_e', group='electrons-energy')
        self.L_e2 = ex.hdf5ExtractScalar('L_e2', group='electrons-energy')
        self.Q_nth = ex.hdf5ExtractScalar('Q_nth', group='electrons-energy')
        self.Q_nth2 = ex.hdf5ExtractScalar('Q_nth2', group='electrons-energy')
        self.Q_nth3 = ex.hdf5ExtractScalar('Q_nth3', group='electrons-energy')


        self.t = ex.hdf5Extract1D('time')
        self.t_obs = ex.hdf5Extract1D('t_obs')
        self.nu = ex.hdf5Extract1D('nu')
        self.nu_obs = ex.hdf5Extract1D('nu_obs')
        self.g = ex.hdf5Extract1D('gamma')
        self.jnut = ex.hdf5Extract2D('jnut')
        self.jsyn = ex.hdf5Extract2D('jmbs')
        self.jssc = ex.hdf5Extract2D('jssc')
        self.jeic = ex.hdf5Extract2D('jeic')
        self.anut = ex.hdf5Extract2D('anut')
        self.ambs = ex.hdf5Extract2D('ambs')
        self.n = ex.hdf5Extract2D('n_e')
        self.gdotty = ex.hdf5Extract2D('dgdt')






def run_BlazeMag():
    rr = ar.Runner(flabel=outfile,comp_kw={'OMP': True, 'HDF5': True,'compileDir':'/home/zach/Documents/Code_Projects/paramo/'})
    ##adjust parameters
    rr.par.numax = 1e28
    rr.par.wParams()
    ###
    rr.run_blazMag(cmd_args=(True, True),clean=True)

def get_bm_results():
    bmr = bm_results(outfile=outfile+'.jp.h5')
    return bmr

def bmr_n_plot():
    bmr = get_bm_results()

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')



    g = bmr.g
    t = bmr.t
    n = bmr.n

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


def bmr_j_plots():
    bmr = get_bm_results()

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')



    nu = bmr.nu
    t = bmr.t
    jm = bmr.jsyn
    jssc = bmr.jssc
    jeic = bmr.jeic

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
    ax.set_xlim(bmr.numin, bmr.numax)

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

# run_BlazeMag()
# get_bm_results()
bmr_n_plot()
bmr_j_plots()