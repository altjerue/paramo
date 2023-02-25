import os
import random
import Eduviges.SRtoolkit
import matplotlib.pyplot as plt
import numpy as np
import Arriero as ar
import Eduviges.extractor as extr
import Eduviges.magnetobrem as mb
import matplotlib.cm as cm
import scipy.integrate as intergrate
import scipy.special as scisp
import Eduviges.constants as aCons
import matplotlib
import analytical_solutions as ansol
import Plotter as PL
import ModelFitter.PlotExplorer as PE
from Eduviges import spectra as spec

outfile = 'katarzynski_2006_test'

plots_folder = '/home/zach/Documents/Code_Projects/paramo/Plots/'


def get_data(file):
    x = []
    y = []
    fn = file
    f = open(fn, 'r')
    ls = f.readlines()
    for l in ls:
        l = l.split(", ")
        x.append(float(l[0]))
        y.append(float(l[1]))
    return x, y

class katarzynski_2006_results:

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
        self.Inu7 = ex.hdf5Extract2D('Inu7', group='Numeric')
        self.Inu8 = ex.hdf5Extract2D('Inu8', group='Numeric')
        self.dotgKN7 = ex.hdf5Extract1D('dotgKN7', group='Numeric')
        self.dotgKN8 = ex.hdf5Extract1D('dotgKN8', group='Numeric')
        self.Ntot1 = ex.hdf5Extract1D('Ntot1', group='Numeric')
        self.Ntot2 = ex.hdf5Extract1D('Ntot2', group='Numeric')
        self.Ntot3 = ex.hdf5Extract1D('Ntot3', group='Numeric')
        self.Ntot4 = ex.hdf5Extract1D('Ntot4', group='Numeric')
        self.Ntot5 = ex.hdf5Extract1D('Ntot5', group='Numeric')
        self.Ntot6 = ex.hdf5Extract1D('Ntot6', group='Numeric')
        self.Ntot7 = ex.hdf5Extract1D('Ntot7', group='Numeric')
        self.Ntot8 = ex.hdf5Extract1D('Ntot8', group='Numeric')
        self.n1 = ex.hdf5Extract2D('n1', group='Numeric')
        self.n2 = ex.hdf5Extract2D('n2', group='Numeric')
        self.n3 = ex.hdf5Extract2D('n3', group='Numeric')
        self.n4 = ex.hdf5Extract2D('n4', group='Numeric')
        self.n5 = ex.hdf5Extract2D('n5', group='Numeric')
        self.n6 = ex.hdf5Extract2D('n6', group='Numeric')
        self.n7 = ex.hdf5Extract2D('n7', group='Numeric')
        self.n8 = ex.hdf5Extract2D('n8', group='Numeric')
        self.gdot7 = ex.hdf5Extract2D('gdot7', group='Numeric')
        self.gdot8 = ex.hdf5Extract2D('gdot8', group='Numeric')
        self.jmbs7 = ex.hdf5Extract2D('jmbs7', group='Numeric')
        self.jmbs8 = ex.hdf5Extract2D('jmbs8', group='Numeric')
        self.jssc7 = ex.hdf5Extract2D('jssc7', group='Numeric')
        self.jssc8 = ex.hdf5Extract2D('jssc8', group='Numeric')
        self.ambs7 = ex.hdf5Extract2D('ambs7', group='Numeric')
        self.agg7 = ex.hdf5Extract2D('agg7', group='Numeric')
        self.ambs8 = ex.hdf5Extract2D('ambs8', group='Numeric')
        self.ns = [self.n1, self.n2, self.n3, self.n4, self.n5, self.n6,self.n7,self.n8]
        self.jmbss = [self.jmbs7, self.jmbs8]
        self.ambss = [self.ambs7, self.ambs8]
        self.jsscs = [self.jssc7, self.jssc8]

def run_katarzynski_2006():
    rr = ar.Runner(flabel=outfile, comp_kw={'OMP': True, 'HDF5': True, 'compileDir': '../'})
    ##adjust parameters
    rr.par.NG = 100#1000
    rr.par.NT = 100 #10000
    rr.par.NF = 100
    rr.par.g1 = 1e4
    rr.par.g2 = 1e6
    rr.par.gmin = 1.01e0
    rr.par.gmax = 1.5e2 * rr.par.g2
    rr.par.numin = 1e9
    rr.par.numax = 1e29
    rr.par.pind = 0e0
    rr.par.lg1 = 'T'
    rr.par.wParams()
    ###
    rr.run_test(clean=True, test_choice=7)


def get_katarzynski_2006():
    kr = katarzynski_2006_results(outfile=outfile + '.jp.h5')
    return kr


def data_compare(x,y,data_x,data_y):
    lx = np.log10(x)
    ly = np.log10(y)
    ldata_x = np.log10(data_x)
    ldata_y = np.log10(data_y)
    r = lx**2 + ly**2
    r = np.sqrt(r)
    rdata_x =[]
    rdata_y =[]
    for i in range(len(data_x)):
        rdata = ldata_x[i]**2 + ldata_y[i]**2
        rdata = np.sqrt(rdata)
        minv = np.min(np.abs(r-rdata))
        if(minv<5e-3):
            rdata_x.append(data_x[i])
            rdata_y.append(data_y[i])
    return np.array(rdata_x),np.array(rdata_y)


def n_plot(n,g,t,results,compare_fig):

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')

    no_data=False
    if(compare_fig=='fig1_a'):
        # data_x , data_y = get_data("./katarzynski_2006_dig_extr_data/fig1_a")
        ax.set_ylim(1e-12, 1e1)
        ax.set_xlim(1e0, 1e7)
    elif(compare_fig=='fig1_b'):
        # data_x, data_y = get_data("./katarzynski_2006_dig_extr_data/fig1_b")
        ax.set_ylim(1e-10, 1e3)
        ax.set_xlim(1e0, 1e7)
    elif(compare_fig == 'fig1_c'):
        # data_x, data_y = get_data("./katarzynski_2006_dig_extr_data/fig1_c")
        ax.set_ylim(1e-9, 1e4)
        ax.set_xlim(1e0, 1e7)
    elif (compare_fig == 'fig2_a'):
        # data_x, data_y = get_data("./katarzynski_2006_dig_extr_data/fig1_c")
        ax.set_ylim(1e-10, 1e5)
        ax.set_xlim(1e0, 1e7)
    elif (compare_fig == 'fig2_b'):
        # data_x, data_y = get_data("./katarzynski_2006_dig_extr_data/fig1_c")
        ax.set_ylim(1e-8, 1e7)
        ax.set_xlim(1e0, 1e7)
    elif (compare_fig == 'fig2_c'):
        # data_x, data_y = get_data("./katarzynski_2006_dig_extr_data/fig1_c")
        ax.set_ylim(1e-8, 1e7)
        ax.set_xlim(1e0, 1e7)
    elif (compare_fig == 'fig3_a'):
        data_x, data_y = get_data("./katarzynski_2006_dig_extr_data/fig3_a")
        ax.set_ylim(1e-12, 1e2)
        ax.set_xlim(1e0, 1e8)
    else:
        no_data=True
        ax.set_ylim(1e-9, 1e4)
        ax.set_xlim(1e0, 1e7)
    


    cmap = cm.rainbow
    sm = plt.cm.ScalarMappable(cmap=cmap,  norm=matplotlib.colors.LogNorm(vmin=t[0], vmax=t[-1]))

    pls = []
    ts = [0.01,0.1,0.5,1,2,3,5,10,15,20,40,60,80]
    t_tac = t / results.tacc
    i_list=[]
    t_tac_2 = []
    for i in range(len(ts)):
        i_list.append(np.argmin(np.abs(t_tac-ts[i])))
    for i in range(len(t)):
        if(i not in i_list and i != len(t) - 1  and i!=0):
            if(i%10 !=0):
                continue
        t_tac_2.append(t_tac[i])
        y = n[:,i]
        if(no_data != True):
            # xdata, ydata = data_compare(g,y,data_x,data_y)
            pl2 = ax.scatter(data_x,data_y,s=4,marker='d',color='black')
        pl = ax.plot(g,y,color=cmap(i))
        pls.append(pl)



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


def nuFnu_plots(jm,nu,t,results,compare_fig):


    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')


    cmap = cm.rainbow
    sm = plt.cm.ScalarMappable(cmap=cmap,  norm=matplotlib.colors.LogNorm(vmin=t[0], vmax=t[-1]))

    pls = []
    my = None

    ts = [0.01, 0.1, 0.5, 1, 2, 3, 5, 10, 15,20,40,60,80]
    t_tac = t / results.tacc
    i_list = []
    t_tac_2 = []
    for i in range(len(ts)):
        i_list.append(np.argmin(np.abs(t_tac - ts[i])))
    for i in range(len(t)):
        if (i not in i_list and i != len(t) - 1 and i != 0):
            if (i % 10 != 0):
                continue
        t_tac_2.append(t_tac[i])
        y = spec.Luminosity(nu,jm[i,:],results.ambs7[i,:],21,results.R,(4/3)*np.pi*(results.R**3))/(4*np.pi*((4.32e26)**2))
        # y=jm[i,:]
        # xdata, ydata = data_compare(g,y,data_x,data_y)
        # pl2 = ax.scatter(data_x,data_y,s=8,marker='d',color='black')
        nu = nu
        pl = ax.plot(nu, y, color=cmap(i))
        pls.append(pl)
        if(my is None):
            my = max(y)
        else:
            if(my<max(y)):
                my = max(y)



    if (compare_fig == 'fig3_b'):
        # data_x, data_y = get_data("./katarzynski_2006_dig_extr_data/fig1_c")
        ax.set_ylim(1e-16, 1e-5)
        ax.set_xlim(1e10, 1e22)
    else:
        # ax.set_ylim(my / 1e10, 2 * my)
        ax.set_xlim(nu[0], nu[-1])

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
    ax.set_ylabel(r"$\nu F_{\nu}$ $[\frac{erg}{s cm^{2}}]$",fontsize=18)

    plt.tight_layout()

    # plt.savefig(plots_folder+"n1vsg.png")
    plt.show()

def abs_plots(am,nu,t,results,compare_fig):


    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')


    cmap = cm.rainbow
    sm = plt.cm.ScalarMappable(cmap=cmap,  norm=matplotlib.colors.LogNorm(vmin=t[0], vmax=t[-1]))

    pls = []
    my = None

    ts = [0.01, 0.1, 0.5, 1, 2, 3, 5, 10, 15,20,40,60,80]
    t_tac = t / results.tacc
    i_list = []
    t_tac_2 = []
    for i in range(len(ts)):
        i_list.append(np.argmin(np.abs(t_tac - ts[i])))
    for i in range(len(t)):
        if (i not in i_list and i != len(t) - 1 and i != 0):
            if (i % 10 != 0):
                continue
        t_tac_2.append(t_tac[i])
        y = am[i,:]
        nu = nu
        pl = ax.plot(nu, y, color=cmap(i))
        pls.append(pl)
        if(my is None):
            my = max(y)
        else:
            if(my<max(y)):
                my = max(y)



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
    ax.set_ylabel(r"$\frac{d\tau}{dx}$",fontsize=18)

    plt.tight_layout()

    # plt.savefig(plots_folder+"n1vsg.png")
    plt.show()

run_katarzynski_2006()
k_results = get_katarzynski_2006()
compare_figs =['fig1_a','fig1_b','fig1_c','fig2_a','fig2_b','fig2_c','fig3_a','fig4_a']
compare_synfigs=['fig3_b','fig4_b']
compare_sscfigs=['fig3_c','fig4_c']
for i in range(6,len(k_results.ns)):
    n_plot(k_results.ns[i],k_results.g,k_results.t,results=k_results,compare_fig=compare_figs[i])
    if(i>5):
        nuFnu_plots(k_results.jmbss[i-6],k_results.nu,k_results.t,results=k_results,compare_fig=compare_synfigs[i-6])
        print("ssc1")
        nuFnu_plots(k_results.jssc7,k_results.nu,k_results.t,results=k_results,compare_fig=compare_sscfigs[i-6])
        abs_plots(k_results.agg7,k_results.nu,k_results.t,results=k_results,compare_fig='asldfkj')

test = k_results.jssc7
test2 = k_results.jssc8
nuFnu_plots(k_results.Inu7,k_results.nu,k_results.t,results=k_results,compare_fig=compare_sscfigs[0])
n_plot(k_results.gdot7,k_results.g,k_results.t,results=k_results,compare_fig='')
n_plot(k_results.gdot7[:,0:1],k_results.g,k_results.t[0:1],results=k_results,compare_fig='')
n_plot(k_results.gdot8,k_results.g,k_results.t,results=k_results,compare_fig='')
n_plot(k_results.gdot8[:,0:1],k_results.g,k_results.t[0:1],results=k_results,compare_fig='')

# def test_sum_f(n,x):
#     f = (-1)**(n-1)
#     f = f*(n**(-2))
#     f = f*(x**(-n))
#     return f
#
# def test_sum(nmax,xrange):
#     if(len(xrange)>1):
#         f=np.zeros(len(xrange))
#         for i in range(len(xrange)):
#             for n in range(1,int(nmax)):
#                 f[i]=f[i] + test_sum_f(n,xrange[i])
#     else:
#         f=0
#         for n in range(1, int(nmax)):
#             f = f + test_sum_f(n, xrange[0])
#     return f
#
# fig, ax = plt.subplots()
# # ax.set_yscale('log')
# ax.set_xscale('log')
# cmap = cm.rainbow
# x = np.linspace(1,10)
# nrange = np.logspace(0,4,50)
# # sm = plt.cm.ScalarMappable(cmap=cmap, norm=matplotlib.colors.Normalize(vmin=nrange[0], vmax=nrange[-1]))
# # for i in range(len(nrange)):
# #     pl = ax.plot(x, test_sum(i,x), color=cmap(i))
# #
# # cbar = plt.colorbar(sm, anchor=(-0.6, 0.0), ticks=np.logspace(0, 7, 8))
# # cbar.ax.minorticks_off()
# # cbar.ax.set_ylabel(r"nrange", fontsize=18)
# y = []
# for n in nrange:
#     y.append(test_sum(n,[10]))
# ax.plot(nrange,y)
# ax.legend()
# plt.show()