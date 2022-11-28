import Eduviges.constants
import matplotlib.pyplot as plt
import numpy as np
import Arriero as ar
import Eduviges.extractor as extr
import matplotlib.cm as cm
import matplotlib
import Eduviges.constants as acons

outfile = 'turblaz_test'
outdir='./'
plots_folder = '/home/zach/Documents/Code_Projects/paramo/Example/Plots/'


class tb_results:

    def __init__(self, outfile):
        ex = extr.fromHDF5(h5fname=outdir+outfile)
        self.numt = ex.hdf5ExtractScalar('numdt', group='Parameters')
        self.numf = ex.hdf5ExtractScalar('numdf', group='Parameters')
        self.numg = ex.hdf5ExtractScalar('numbins', group='Parameters')
        self.t_max = ex.hdf5ExtractScalar('t_max', group='Parameters')
        self.tstep = ex.hdf5ExtractScalar('tstep', group='Parameters')
        self.Gamma_bulk = ex.hdf5ExtractScalar('Gamma_bulk', group='Parameters')
        self.sigma = ex.hdf5ExtractScalar('sigma', group='Parameters')
        self.gmin = ex.hdf5ExtractScalar('gamma_min', group='Parameters')
        self.gmax = ex.hdf5ExtractScalar('gamma_max', group='Parameters')
        self.u_ext = ex.hdf5ExtractScalar('u_ext', group='Parameters')
        self.nu_ext = ex.hdf5ExtractScalar('nu_ext', group='Parameters')

        self.numin = ex.hdf5ExtractScalar('nu_min', group='Parameters')
        self.numax = ex.hdf5ExtractScalar('nu_max', group='Parameters')

        self.R = ex.hdf5ExtractScalar('R')
        self.L_j = ex.hdf5ExtractScalar('Lj')
        self.gam0 = ex.hdf5ExtractScalar('gam0')
        self.tc = ex.hdf5ExtractScalar('tc')
        self.tesc = ex.hdf5ExtractScalar('t_esc')
        self.B_0 = ex.hdf5ExtractScalar('Bfield')
        self.n0 = ex.hdf5ExtractScalar('n0')
        self.uB = ex.hdf5ExtractScalar('uB')
        self.uph = ex.hdf5ExtractScalar('uph')
        self.t = ex.hdf5Extract1D('time')
        self.nu = ex.hdf5Extract1D('nu')
        self.ubol = ex.hdf5Extract1D('ubol')
        self.g = ex.hdf5Extract1D('gamma')
        self.dotgkn = ex.hdf5Extract1D('dotgkn')
        self.jnut = ex.hdf5Extract2D('jnut')
        self.jsyn = ex.hdf5Extract2D('jmbs')
        self.jssc = ex.hdf5Extract2D('jssc')
        self.jeic = ex.hdf5Extract2D('jeic')
        self.anut = ex.hdf5Extract2D('anut')
        self.ambs = ex.hdf5Extract2D('ambs')
        self.n = ex.hdf5Extract2D('n_e')
        self.gdotty = ex.hdf5Extract2D('dgdt')






def run_Turblaz():
    rr = ar.Runner(flabel=outfile,comp_kw={'OMP': True, 'HDF5': True,'compileDir':'/home/zach/Documents/Code_Projects/paramo/'})
    ##adjust parameters
    rr.par.lg1 = 'F' #only calculate radiation at the end
    rr.par.ed1 = 0.612 #del_ph
    rr.par.ed2 = 0.212 # ninj=par_ed3 !efficienty of driving frequency to plasma
    rr.par.R0 = 0.03 #rtm
    rr.par.R = 1e18
    rr.par.sigma = 3.82 #sigma
    rr.par.gmax = 1e8
    rr.par.NG =150
    rr.par.NT = 150
    rr.par.gmin = 1e0 + 1e-2
    rr.par.gamma_bulk=17.57
    rr.par.L_jet = 3.295e46
    rr.par.nu_ext = 1e15*rr.par.gamma_bulk
    rr.par.numin = 1e10
    rr.par.numax = 1e28
    rr.par.nu_ext = 1e18
    R_BLR = (1e17 / np.sqrt(1e45)) * np.sqrt(0.1 * rr.par.L_jet )
    if (rr.par.R > R_BLR):
        rr.par.ed1 = rr.par.ed1 * (rr.par.R / R_BLR) ** (-3)

    rr.par.wParams()
    ###
    rr.run_turb(cmd_args=(False, True),clean=True) #cmd_arges cool_withKN,with_abs

def get_tb_results():
    tbr = tb_results(outfile=outfile+'.jp.h5')
    return tbr

def tb_n_plot():
    tbr = get_tb_results()

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')



    g = tbr.g
    t = tbr.t
    tempdgdt = tbr.gdotty
    n = np.zeros((len(t),len(g)))
    # print(tempdgdt[1,:])
    # print(tempdgdt[0,:])
    # for i in range(len(t)):
    #     for j in range(len(g)):
    #         tempdgdt[i,j] =  tempdgdt[i,j]/ ((tempdgdt[0,0]/(g[0]**2))*(g[j]**2))
    # for j  in range(len(t)):
    #     n[j,:] = tbr.dotgkn
    # n = tempdgdt# tbr.n
    n =tbr.n
    # tempdgdt = tempdgdt[:,-1]/((4/3)*acons.sigmaT*acons.cLight*(tbr.uph * tbr.uB)*(tbr.g**2))
    cmap = cm.rainbow
    sm = plt.cm.ScalarMappable(cmap=cmap,  norm=matplotlib.colors.LogNorm(vmin=t[0], vmax=t[-1]))

    pls = []
    my = None

    for i in range(len(t)):
        if(i%10 !=0 and i != len(t) - 1 ):
            continue
        y = n[i,:]
        if(my == None):
            my = max(y)
        if(max(y)> my):
            my= max(y)
        pl = ax.plot(g,y,color=cmap(i))

        pls.append(pl)
    # pl = ax.plot(g, y[0] * (g / g[0]) ** 2)
    # pl = ax.plot(g, (4/3)*Eduviges.constants.cLight*Eduviges.constants.sigmaT*tbr.uB*(g**2))
    # ax.set_ylim(1e-4, 2*my)
    # ax.set_xlim(1e0, 2e13)

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


def tb_j_plots():
    bmr = get_tb_results()

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
        if(i%10 !=0 and i != len(t) - 1 ):
            continue
        y = jm[i,:] + jssc[i,:] + jeic[i,:]
        y = bmr.nu*y
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

run_Turblaz()



# get_bm_results()
tb_n_plot()
tb_j_plots()
