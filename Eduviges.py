#  ===================================================================
#   #####                          #######
#  #     # #####  ######  ####     #       #    # #    #  ####   ####
#  #       #    # #      #    #    #       #    # ##   # #    # #
#   #####  #    # #####  #         #####   #    # # #  # #       ####
#        # #####  #      #         #       #    # #  # # #           #
#  #     # #      #      #    #    #       #    # #   ## #    # #    #
#   #####  #      ######  ####     #        ####  #    #  ####   ####
#  ===================================================================
import SAPyto.spectra as spec
import SAPyto.SRtoolkit as SR
import extractor.fromHDF5 as extr
from Comala import Paramo


#  #                              #####
#  #       #  ####  #    # ##### #     # #    # #####  #    # ######  ####
#  #       # #    # #    #   #   #       #    # #    # #    # #      #
#  #       # #      ######   #   #       #    # #    # #    # #####   ####
#  #       # #  ### #    #   #   #       #    # #####  #    # #           #
#  #       # #    # #    #   #   #     # #    # #   #   #  #  #      #    #
#  ####### #  ####  #    #   #    #####   ####  #    #   ##   ######  ####
def build_LCs(nu_min, nu_max, only_load=True, Jy=False, eV=False, days=False, par_kw={}, comp_kw={}):
    run = Paramo(par_kw=par_kw, comp_kw=comp_kw)
    if not only_load:
        run.run_Paramo()
    file = run.run_Paramo()
    Gbulk = extr.hdf5ExtractScalar(file, 'gamma_bulk', group='Params')
    theta_obs = extr.hdf5ExtractScalar(file, 'view-angle', group='Params')
    D = SR.Doppler(Gbulk, theta_obs)
    nu = extr.hdf5Extract1D(file, 'nu_obs')
    if eV:
        nu = spec.Hz2eV(nu)
    t = extr.hdf5Extract1D(file, 't_obs')
    if days:
        t = spec.sec2dy(t)
    Inu = extr.hdf5Extract2D(file, 'Iobs')
    R = extr.hdf5ExtractScalar(file, 'R', group='Params')
    z = extr.hdf5ExtractScalar(file, 'redshift', group='Params')
    dL = extr.hdf5ExtractScalar(file, 'd_lum', group='Params')
    numdt = t.size
    Fnu = spec.flux_dens(Inu, dL, z, D, R)
    LC = spec.LightCurves()

    if Jy:
        return t, spec.conv2Jy(LC.integ(nu_min, nu_max, numdt, nu, Fnu))
    else:
        return t, LC.integ(nu_min, nu_max, numdt, nu, Fnu)


#     #                   #####
#    # #   #    #  ####  #     # #####  ######  ####
#   #   #  #    # #    # #       #    # #      #    #
#  #     # #    # #       #####  #    # #####  #
#  ####### #    # #  ###       # #####  #      #
#  #     #  #  #  #    # #     # #      #      #    #
#  #     #   ##    ####   #####  #      ######  ####
def build_avgSpec(t_min, t_max, dset='Iobs', only_load=True, Jy=False, eV=False, days=False, par_kw={}, comp_kw={}):
    run = Paramo(par_kw=par_kw, comp_kw=comp_kw)
    if not only_load:
        run.run_Paramo()
    file = run.run_Paramo()
    Gbulk = extr.hdf5ExtractScalar(file, 'gamma_bulk', group='Params')
    theta_obs = extr.hdf5ExtractScalar(file, 'view-angle', group='Params')
    D = SR.Doppler(Gbulk, theta_obs)
    nu = extr.hdf5Extract1D(file, 'nu_obs')
    if eV:
        nu = spec.Hz2eV(nu)
    t = extr.hdf5Extract1D(file, 't_obs')
    if days:
        t = spec.sec2dy(t)
    Inu = extr.hdf5Extract2D(file, 'Iobs')
    R = extr.hdf5ExtractScalar(file, 'R', group='Params')
    z = extr.hdf5ExtractScalar(file, 'redshift', group='Params')
    dL = extr.hdf5ExtractScalar(file, 'd_lum', group='Params')
    numdf = nu.size
    Fnu = spec.flux_dens(Inu, dL, z, D, R)
    sp = spec.spectrum()

    if Jy:
        return nu, spec.conv2Jy(sp.averaged(t_min, t_max, numdf, t, Fnu))
    else:
        return nu, sp.averaged(t_min, t_max, numdf, t, Fnu)
