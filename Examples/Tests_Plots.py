import os

import Eduviges.SRtoolkit
import matplotlib.pyplot as plt
import numpy as np
import Arriero as ar
import Eduviges.extractor as extr
import matplotlib.cm as cm
import scipy.integrate as intergrate
import Eduviges.constants as aCons
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
        self.ns = [self.n1, self.n2, self.n3, self.n4, self.n5, self.n6]
        self.jmbss = [self.jmbs1, self.jmbs4, self.jmbs5, self.jmbs6]
        self.ambss = [self.ambs1, self.ambs4, self.ambs5, self.ambs6]
        self.jsscs = [self.jssc1, self.jssc4, self.jssc5, self.jssc6]


outfile = __file__.split('Tests_Plots.py')[0] + 'steady_state_test'

plots_folder = '/home/zach/Documents/Code_Projects/paramo/Plots/'


def run_steady_state_test():
    rr = ar.Runner(flabel=outfile, comp_kw={'OMP': True, 'HDF5': True, 'compileDir': './'})
    ##adjust parameters
    rr.par.numax = 1e28
    rr.par.NG = 12800
    rr.par.NT = 300
    rr.par.NF = 192
    rr.par.g1 = 1e4
    rr.par.g2 = 1e6
    rr.par.gmin = 1.01e0
    rr.par.gmax = 1.5e0 * rr.par.g2
    rr.par.numin = 1e10
    rr.par.numax = 1e27
    rr.par.pind = 0e0
    rr.par.lg1 = 'F'
    rr.par.wParams()
    ###
    rr.run_test(clean=True, test_choice=1)


def get_steady_state_results():
    ssr = SS_results(outfile=outfile + '.jp.h5')
    return ssr


def run_convergence_test(numts, numgs):
    check_exist=True

    for numt in numts:
        numt = int(numt)
        for numg in numgs:
            if (os.path.exists(outfile + "_" + str(int(numg)) + "_" + str(int(numt)) + "_.jp.h5") and check_exist):
                continue
            numg = int(numg)
            rr = ar.Runner(flabel=outfile + "_" + str(int(numg)) + "_" + str(int(numt)) + "_",
                           comp_kw={'OMP': True, 'HDF5': True, 'compileDir': './'})
            ##adjust parameters
            rr.par.numax = 1e28
            rr.par.NG = numg
            rr.par.NT = numt
            rr.par.NF = 192
            # rr.par.g1 = 1e4
            # rr.par.g2 = 1e6
            # rr.par.gmin = 1e0 + 1e-10
            # rr.par.gmax = 1.5e2 * rr.par.g2
            g1 = 1e4
            g2=1e6
            rr.par.g1 =g1#Eduviges.SRtoolkit.pofg(g1)*aCons.me*aCons.cLight
            rr.par.g2 =g2 #Eduviges.SRtoolkit.pofg(g2)*aCons.me*aCons.cLight
            rr.par.gmin =1e0  #+ 1e-90 #Eduviges.SRtoolkit.pofg(1 + 1e-15)*aCons.me*aCons.cLight
            rr.par.gmax =1.5e2 * g2 #Eduviges.SRtoolkit.pofg(1.5e2 * g2)*aCons.me*aCons.cLight
            rr.par.numin = 1e10
            rr.par.numax = 1e27
            rr.par.pind = 0e0
            rr.par.lg1 = 'F'  # no radiation
            rr.par.wParams()
            ###
            rr.run_test(clean=True, test_choice=1)


def get_convergence_results(numt, numg):
    ssr = SS_results(outfile=outfile + "_" + str(int(numg)) + "_" + str(int(numt)) + "_" + '.jp.h5')
    return ssr


def get_error(ef, efg, ei, eig):
    er = []
    for i in range(len(eig)):
        closest=np.argsort(np.abs(efg - eig[i]))
        if(efg[closest[0]] - eig[i] == 0 ):
            ci = np.argmin(np.abs(efg - eig[i]))
            e = (ef[ci] - ei[i]) ** 2
            er.append(e)
        else:
            eic = ei[i]
            c0 = closest[0]
            c1 = closest[1]
            gf0 = efg[c0]
            gf1 = efg[c1]
            nf0 = ef[c0]
            nf1 = ef[c1]
            ncomp = nf0
            m= (nf1-nf0)/(gf1-gf0)
            ncomp += m*np.abs(eig[i] - gf0)
            e = ((ncomp - eic)) ** 2
            er.append(e)


    er = np.sqrt(sum(er) )/ len(eig)
    # er = np.sqrt(sum(er))
    return er


def get_error2(ef, efg, ei, eig):
    er = []
    # g1cut = 1e1
    # g2cut = 1e5
    # gcmin,gcmax=np.argmin(np.abs(efg-g1cut)),np.argmin(np.abs(efg-g2cut))
    # gcimin,gcimax = np.argmin(np.abs(eig-g1cut)),np.argmin(np.abs(eig-g2cut))
    # ef = ef[gcmin:gcmax]
    # efg = efg[gcmin:gcmax]
    # ei = ei[gcimin:gcimax]
    # eig = eig[gcimin:gcimax]
    for i in range(len(efg)):
        closest=np.argsort(np.abs(eig - efg[i]))
        fucku = True
        if(eig[closest[0]] - efg[i] == 0 or fucku):
            ci = np.argmin(np.abs(eig - efg[i]))
            e = (ei[ci] - ef[i]) ** 2
            er.append(e)
        else:
            efc = ef[i]
            c0 = closest[0]
            c1 = closest[1]
            gi0 = eig[c0]
            gi1 = eig[c1]
            ni0 = ei[c0]
            ni1 = ei[c1]
            ncomp = ni0
            m= np.abs((ni1-ni0)/(gi1-gi0))
            # if(efg[i]<gi0):
            #     m=-m
            ncomp += (m*(efg[i] - gi0))
            e = ((ncomp - efc)/ncomp) ** 2
            # if(efg[i]>= 3e4):
            #     print("afd")
            er.append(e)
    # fig,ax = plt.subplots()
    # # ax.plot(range(len(efg)),er)
    # ax.plot(efg,er)
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    # plt.show()
    eers = np.sqrt(er)/len(eig)
    er = np.sqrt(sum(er) / len(efg))



    return er,eers


def get_error3(ef, efg, ei, eig):
    er = []
    eers =[]
    # g1cut = 1.001e1
    # g2cut = 1.8e7
    # gcmin,gcmax=np.argmin(np.abs(efg-g1cut)),np.argmin(np.abs(efg-g2cut))
    # gcimin,gcimax = np.argmin(np.abs(eig-g1cut)),np.argmin(np.abs(eig-g2cut))
    # ef = ef[gcmin:gcmax]
    # efg = efg[gcmin:gcmax]
    # ei = ei[gcimin:gcimax]
    # eig = eig[gcimin:gcimax]
    # ef = np.log10(ef)
    # efg = np.log10(efg)
    # ei = np.log10(ei)
    # eig = np.log10(eig)
    for i in range(len(eig)):
        closest=np.argsort(np.abs(efg - eig[i]))
        fucku = True
        if(np.abs(efg[closest[0]] - eig[i]) == 0 or fucku):
            ci = np.argmin(np.abs(efg - eig[i]))
            e = np.abs((ef[ci] - ei[i])/ef[ci])
            er.append(e)
        else:
            eic = ei[i]
            c0 = closest[0]
            c1 = closest[1]
            gf0 = efg[c0]
            gf1 = efg[c1]
            nf0 = ef[c0]
            nf1 = ef[c1]
            ncomp = nf0
            m= np.abs((nf1-nf0)/(gf1-gf0))
            # if (efg[i] >= 3e4):
            #     m=-m
            if(efg[i]<gf0):
                m=-m
            ncomp -= (m*np.abs(eig[i] - gf0))
            e = ((ncomp - eic)/ncomp) ** 2
            # if(efg[i]>= 3e4):
            #     print("afd")
            er.append(e)
    eers = np.array(er)/len(eig)
    # fig,ax = plt.subplots()
    # # ax.plot(range(len(efg)),er)
    # ax.plot(eig,er)
    # ax.set_xscale('log')
    # # ax.set_yscale('log')
    # plt.show()

    er = sum(er) / len(eig)



    return er,eers

def get_error4(ef, efg, ei, eig):
    g1cut = 1e1
    g2cut = 1e5
    # gcmin,gcmax=np.argmin(np.abs(efg-g1cut)),np.argmin(np.abs(efg-g2cut))
    # gcimin,gcimax = np.argmin(np.abs(eig-g1cut)),np.argmin(np.abs(eig-g2cut))
    # ef = ef[gcmin:gcmax]
    # efg = efg[gcmin:gcmax]
    # ei = ei[gcimin:gcimax]
    # eig = eig[gcimin:gcimax]
    ef = np.log10(ef)
    efg = np.log10(efg)
    ei = np.log10(ei)
    eig = np.log10(eig)
    eef = np.trapz(ef,efg)
    eei=np.trapz(ei,eig)
    er = ((eef-eei)**2)
    er = np.sqrt(er)/ len(eig)
    return er

def analytical_form(g):
    C0 = 3.48e-11
    tesc = 1e0 / (C0* ((10e0) ** (4.5e0)))
    return (g**2)*np.exp(-2*C0*tesc*(g-1))
def analytical_solution(n0,g):
    x = n0/(intergrate.quad(lambda x: analytical_form(x),1,np.inf)[0])
    return x*analytical_form(g)

def get_error_analytic(ei, eig,n0):

    er = []
    for i in range(len(eig)):
        ef = analytical_solution(n0,eig[i])
        if(ef>0):
            er.append(((ef-ei[i])/ef)**2)
        else:
            er.append(0)
        # er.append(((ef-ei[i]))**2)
    eers = np.sqrt(er)
    er = np.sqrt(sum(er)/len(eig))
    return er,eers


def convergence_plots_analytic(numt,numgs):
    xs = []
    ys = []
    gs = []
    ns = []
    nios=[]
    eers = []
    p = 0e0
    gmin = 1e4
    gmax = 1e6
    n0 = intergrate.quad(lambda x: x**p,gmin,gmax)[0]
    for i in range(len(numgs)):
        mg = numgs[i]
        eissr = get_convergence_results(numt, mg)
        ei = eissr.n4[:, -1]
        eig = eissr.g
        nio = np.trapz(ei,eig)
        nio0 = np.trapz(eissr.n4[:, 0],eig)
        nios.append(nio/nio0)
        # ei = n0*ei/nio

        gminci = np.argmin(np.abs(eig - 10e0))
        gmaxci = np.argmin(np.abs(eig - 1e7))
        eig=eig[gminci:gmaxci]
        ei=ei[gminci:gmaxci]
        er, eer = get_error_analytic(ei, eig,n0)
        gs.append(eig)
        ns.append(ei)
        eers.append(eer)
        ys.append(er)
        xs.append(len(eig))

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')
    fig1, ax1 = plt.subplots()
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    fig2, ax2 = plt.subplots()
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    fig3, ax3 = plt.subplots()
    # ax3.set_yscale('log')
    ax3.set_xscale('log')
    pl2 = ax2.plot()
    pl = ax.scatter(xs, ys)
    pl3 = ax3.scatter(xs,nios)
    xc = np.logspace(0.5, 1.7)
    xc2 = np.logspace(1.5, 2.9)
    pl = ax.plot(xc*xc[1], ys[1] * (xc / xc[1]) ** -2, label='p=-2')
    pl = ax.plot(xc2, ys[4] * (xc2 / xc2[0]) ** -1, label='p=-1')
    for i in range(len(gs)):
            pl2 = ax2.scatter(gs[i], eers[i], 4)
            pl2 = ax1.plot(gs[i], ns[i], label='line: ' + str(i))
    ax1.plot(gs[-1],analytical_solution(n0,gs[-1]),'--')
    ax1.legend()
    ax.legend()
    plt.show()

def convergence_plots_numg(numt, numgs):
    mgi = np.argmax(numgs)
    mg = int(numgs[mgi])
    numt = int(numt)
    efssr = get_convergence_results(numt, mg)
    ef = efssr.n4[:, -1]

    for ii in range(len(ef)):
        if ef[ii]<1e-150:
            ef[ii]=1e-150
    efg = efssr.g
    xs = []
    ys = []
    gs =[]
    ns =[]
    eers=[]
    g1cut = 5.001e0
    g2cut = 1.8e8
    gcmin,gcmax=np.argmin(np.abs(efg-g1cut)),np.argmin(np.abs(efg-g2cut))

    ef = ef[gcmin:gcmax]
    efg = efg[gcmin:gcmax]

    for numg in numgs:
        numg = int(numg)
        if(numg == mg ):
            continue
        eissr = get_convergence_results(numt, numg)
        ei = eissr.n4[:, -1]
        eig = eissr.g
        for ii in range(len(ei)):
            if ei[ii] < 1e-150:
                ei[ii] = 1e-150
        gcimin,gcimax = np.argmin(np.abs(eig-g1cut)),np.argmin(np.abs(eig-g2cut))
        ei = ei[gcimin:gcimax]
        eig = eig[gcimin:gcimax]
        er,eer= get_error3(ef, efg, ei, eig)
        gs.append(eig)
        ns.append(ei)
        # er = get_error(ef, efg, ei, eig)
        eers.append(eer)
        ys.append(er)
        xs.append(len(eig))
    gs.append(efg)
    ns.append(ef)
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')
    fig1, ax1 = plt.subplots()
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    fig2, ax2 = plt.subplots()
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    pl2 = ax2.plot()
    pl = ax.scatter(xs, ys)
    xc = np.logspace(0.5,1.7)
    xc2 = np.logspace(1.5,2.9)
    pl = ax.plot(xc,ys[0]*(xc/xc[0])**-3,label='p=-3')
    pl = ax.plot(xc2,ys[2]*(xc2/xc2[0])**-(1),label='p=-1')
    for i in range(len(gs)):

        if(i==len(gs)-1):
            pl2 = ax1.plot(gs[i],ns[i],'--',label='final')

        else:
            pl2 = ax2.scatter(gs[i], eers[i],4)
            pl2 = ax1.plot(gs[i], ns[i],label = 'line: ' + str(i))
    ax1.legend()
    ax.legend()
    plt.show()


def convergence_plots_numt():
    print("TBA")


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
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=matplotlib.colors.LogNorm(vmin=t[0], vmax=t[-1]))

        pls = []
        my = None

        for i in range(len(t)):
            if (i % 25 != 0 and i != len(t) - 1):
                continue
            y = n[:, i]
            if (my == None):
                my = max(y)
            if (max(y) > my):
                my = max(y)
            pl = ax.plot(g, y, color=cmap(i))
            pls.append(pl)

        ax.set_ylim(1e-4, 2 * my)
        ax.set_xlim(1e0, 2e6)

        cbticks = []
        for i in range(8):
            cbticks.append(r"$10^{{{0}}}$".format(i))
        cbar = plt.colorbar(sm, anchor=(-0.6, 0.0), ticks=np.logspace(0, 7, 8))
        cbar.ax.set_yticklabels(cbticks)
        cbar.ax.minorticks_off()
        cbar.ax.set_ylabel(r"t [s]", fontsize=18)
        cbar.ax.tick_params(
            labelsize=15
        )

        plt.tick_params(
            axis='x',
            which='minor',
            bottom=False,
            top=False,
            labelbottom=False
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

        ax.set_xlabel(r"$\gamma$", fontsize=18)
        ax.set_ylabel(r"n $[cm^{-3}]$", fontsize=18)

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
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=matplotlib.colors.LogNorm(vmin=t[0], vmax=t[-1]))

        pls = []
        my = None

        for i in range(len(t)):
            if (i % 25 != 0 and i != len(t) - 1):
                continue
            y = jm[i, :] + jssc[i, :]
            if (my == None):
                my = max(y)
            if (max(y) > my):
                my = max(y)
            pl = ax.plot(nu, y, color=cmap(i))
            pls.append(pl)

        ax.set_ylim(my / 1e10, 2 * my)
        ax.set_xlim(ssr.numin, ssr.numax)

        cbticks = []
        for i in range(8):
            cbticks.append(r"$10^{{{0}}}$".format(i))
        cbar = plt.colorbar(sm, anchor=(-0.6, 0.0), ticks=np.logspace(0, 7, 8))
        cbar.ax.set_yticklabels(cbticks)
        cbar.ax.minorticks_off()
        cbar.ax.set_ylabel(r"t [s]", fontsize=18)
        cbar.ax.tick_params(
            labelsize=15
        )

        plt.tick_params(
            axis='x',
            which='minor',
            bottom=False,
            top=False,
            labelbottom=False
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

        ax.set_xlabel(r"$\nu$ [Hz]", fontsize=18)
        ax.set_ylabel(r"j $[\frac{erg}{s cm^{3}}]$", fontsize=18)

        plt.tight_layout()

        # plt.savefig(plots_folder+"n1vsg.png")
        plt.show()

def build_g(gmax,gmin,numg):
    g=np.zeros(numg)
    g3=np.zeros(numg)
    xs = []
    for k in range(1,numg+1):
        g[k-1] = gmin * (gmax / gmin)**(float(k-1) / float(numg-1))
        g3[k-1] = (gmax - gmin)*k/numg
        xs.append(k)
    g2 = np.logspace(np.log10(gmin),np.log10(gmax),numg,endpoint=True)
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    # ax.set_xscale('log')
    ax.plot(xs,g)
    ax.plot(xs,g2,'--')
    ax.plot(xs,g3,'-.')
    plt.show()
    return g

# g1 = build_g(1.5e0 * 1e6,1.01e0,8)
# g2 = build_g(1.5e0 * 1e6,1.01e0,64)
# print("")
# run_steady_state_test()
# ssr_n_plots()
# ssr_j_plots()
# run_convergence_test(np.logspace(1,4,4),np.logspace(1,4,4))
# garr = [5,10,20,30,50,60,80,100,150,200,300,400,500,600,800,2000,3000,3500]
garr = [20,30,60,80,100,150,200,250,300,400,500,700,900,1000,3000,6000,10000]#,3500]
numtarr =[300]
run_convergence_test(numtarr,garr)
# convergence_plots_numg(300, garr)
# build_g(1.5e6,1.0001,9000)

def debug_distribsFP():
    numg = int(5000)
    numt = 2
    rr = ar.Runner(flabel=outfile + "_" + str(int(numg)) + "_" + str(int(numt)) + "_",
                   comp_kw={'OMP': True, 'HDF5': True, 'compileDir': './'})
    ##adjust parameters
    rr.par.numax = 1e28
    rr.par.NG = numg
    rr.par.NT = numt
    rr.par.NF = 192
    rr.par.g1 = 1e2
    rr.par.g2 = 1e6
    rr.par.gmin = 1.01e0
    rr.par.gmax = 1.5e1 * rr.par.g2
    rr.par.numin = 1e10
    rr.par.numax = 1e27
    rr.par.pind = 0e0
    rr.par.lg1 = 'F'  # no radiation
    rr.par.wParams()
    ###
    rr.run_test(clean=True, test_choice=1)
    eissr = get_convergence_results(numt, numg)
    print("afd")

# debug_distribsFP()
convergence_plots_analytic(numtarr[0],garr)