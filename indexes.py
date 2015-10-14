from math import *
import numpy as np

from tools import integrate

# Index parameters
indparams = {'NaID':
             {'w1': {'ww0': 5895.92, 'dww': 0.5},  # Core D1
              'w2': {'ww0': 5889.95, 'dww': 0.5},  # Core D2
              'wr1': {'ww0': 5805.0, 'dww': 10.0},
              'wr2': {'ww0': 6090.0, 'dww': 20.0}
              },

             'Halpha_cincunegui':
                 {'w1': {'ww0': 6562.808, 'dww': 0.6},
                  'wr1': {'ww0': 6605.0, 'dww': 20.0},
                  'wr2': {'ww0': 6605.0, 'dww': 0.0}
                  },

             'Halpha_narrow':
                 {'w1': {'ww0': 6562.808, 'dww': 40.0},
                  'wr1': {'ww0': 6605.0, 'dww': 20.0},
                  'wr2': {'ww0': 6605.0, 'dww': 0.0}
                  },

             'Halpha_dasilva':
                 {'w1': {'ww0': 6562.808, 'dww': 1.6},
                  'wr1': {'ww0': 6550.87, 'dww': 10.75},
                  'wr2': {'ww0': 6580.31, 'dww': 8.75}
                  },

             'CaII':
                 {'w1': {'ww0': 3933.664, 'dww': 2.0},
                  'w2': {'ww0': 3968.470, 'dww': 2.0},
                  'wr1': {'ww0': 4001.07, 'dww': 20.0},
                  'wr2': {'ww0': 3901.07, 'dww': 20.0}
                  },

             }


def compute_flux(e2ds, blaze, weight, ww, wwmin, wwmax, noise=0.0):

    # Find the spectral order(s) where the relevant window is located
    orders = np.unique(np.where((ww > wwmin) * (ww < wwmax))[0])

    # Initialize arrays
    flux = np.zeros(len(orders))
    vflux = np.zeros(len(orders))

    for i, order in enumerate(orders):
        #
        e2ds_order = e2ds[order]
        blaze_order = blaze[order]
        ww_order = ww[order]
        weight_order = weight[order]

        # Normalize
        e2dsn = e2ds_order/blaze_order

        # Compute variance in e2dsn (error is sqrt(e2ds))
        var_e2dsn = (e2ds_order + noise**2)/blaze_order**2

        flux[i], vflux[i] = integrate(ww_order, e2dsn, weight_order,
                                      wwmin, wwmax, vary=var_e2dsn)

    # Obtain weighted average of fluxes, if present in more than one order
    ff = np.sum(flux/vflux)/np.sum(1/vflux)
    vv = 1.0/np.sum(1.0/vflux)

    return ff, vv


def CaII(ww, e2ds, blaze, noise=0.0, full_output=False):

    # Weight function in 1.0
    wfunc = ww*0.0 + 1.0

    # Select indices
    windows = indparams['CaII'].copy()

    # Compute wavelngth limits for computation
    for win in ('w1', 'w2', 'wr1', 'wr2'):
        windows[win]['wwmin'] = windows[win]['ww0'] - windows[win]['dww']*0.5
        windows[win]['wwmax'] = windows[win]['ww0'] + windows[win]['dww']*0.5

    # Compute fluxes for each window

    # Triangular weight function
    wfunc = np.where(ww >= windows['w1']['ww0'],
                     -(ww - windows['w1']['ww0'])/1.09 + 1.0,
                     (ww - windows['w1']['ww0'])/1.09 + 1.0)
    wfunc = np.where(wfunc > 0, wfunc, 0.0)

    flux_k, var_fluxk = compute_flux(e2ds, blaze, wfunc, ww,
                                     windows['w1']['wwmin'],
                                     windows['w1']['wwmax'], noise)

    # ## Triangular weight function
    wfunc = np.where(ww >= windows['w2']['ww0'],
                     -(ww - windows['w2']['ww0'])/1.09 + 1.0,
                     (ww - windows['w2']['ww0'])/1.09 + 1.0)
    wfunc = np.where(wfunc > 0, wfunc, 0.0)

    flux_h, var_fluxh = compute_flux(e2ds, blaze, wfunc, ww,
                                     windows['w2']['wwmin'],
                                     windows['w2']['wwmax'], noise)
    flux_r1, var_fluxr1 = compute_flux(e2ds, blaze, wfunc, ww,
                                       windows['wr1']['wwmin'],
                                       windows['wr1']['wwmax'], noise)
    flux_r2, var_fluxr2 = compute_flux(e2ds, blaze, wfunc, ww,
                                       windows['wr2']['wwmin'],
                                       windows['wr2']['wwmax'], noise)

    # Compute NaD index
    s = (flux_k + flux_h)/(flux_r1 + flux_r2)

    # Compute error
    var_s = (var_fluxk + var_fluxh)/(flux_r1 + flux_r2)**2 + \
        s**2/(flux_r1 + flux_r2)**2 * (var_fluxr1 + var_fluxr2)

    if full_output:
        return s, sqrt(var_s), np.array([[flux_k, var_fluxk],
                                        [flux_h, var_fluxh],
                                        [flux_r1, var_fluxr1],
                                        [flux_r2, var_fluxr2]])
    else:
        return s, sqrt(var_s)


def NaID(ww, e2ds, blaze, noise=0, full_output=False):

    # Weight function in 1.0
    wfunc = ww*0.0 + 1.0

    # Select indices
    windows = indparams['NaID'].copy()

    # Compute wavelngth limits for computation
    for win in ('w1', 'w2', 'wr1', 'wr2'):
        windows[win]['wwmin'] = windows[win]['ww0'] - windows[win]['dww']*0.5
        windows[win]['wwmax'] = windows[win]['ww0'] + windows[win]['dww']*0.5

    # Compute fluxes for each window
    flux_d1, var_fluxd1 = compute_flux(e2ds, blaze, wfunc, ww,
                                       windows['w1']['wwmin'],
                                       windows['w1']['wwmax'], noise)
    flux_d2, var_fluxd2 = compute_flux(e2ds, blaze, wfunc, ww,
                                       windows['w2']['wwmin'],
                                       windows['w2']['wwmax'], noise)
    flux_r1, var_fluxr1 = compute_flux(e2ds, blaze, wfunc, ww,
                                       windows['wr1']['wwmin'],
                                       windows['wr1']['wwmax'], noise)
    flux_r2, var_fluxr2 = compute_flux(e2ds, blaze, wfunc, ww,
                                       windows['wr1']['wwmin'],
                                       windows['wr2']['wwmax'], noise)

    # Compute NaD index
    d = (flux_d1 + flux_d2)/(flux_r1 + flux_r2)

    # Compute error
    var_d = (var_fluxd1 + var_fluxd2)/(flux_r1 + flux_r2)**2 + \
        d**2/(flux_r1 + flux_r2)**2 * (var_fluxr1 + var_fluxr2)

    if full_output:
        return d, sqrt(var_d), np.array([[flux_d1, var_fluxd1],
                                        [flux_d2, var_fluxd2],
                                        [flux_r1, var_fluxr1],
                                        [flux_r2, var_fluxr2]])
    else:
        return d, sqrt(var_d)


def Halpha(ww, e2ds, blaze, noise=0, full_output=False, version='cincunegui'):

    # Weight function in 1.0
    wfunc = ww*0.0 + 1.0

    # Select indices
    windows = indparams['Halpha_'+version].copy()

    # Compute wavelngth limits for computation
    for win in ('w1', 'wr1', 'wr2'):
        windows[win]['wwmin'] = windows[win]['ww0'] - windows[win]['dww']*0.5
        windows[win]['wwmax'] = windows[win]['ww0'] + windows[win]['dww']*0.5

    # Compute fluxes for each window
    flux_h, var_fluxh = compute_flux(e2ds, blaze, wfunc, ww,
                                     windows['w1']['wwmin'],
                                     windows['w1']['wwmax'], noise)
    flux_r1, var_fluxr1 = compute_flux(e2ds, blaze, wfunc, ww,
                                       windows['wr1']['wwmin'],
                                       windows['wr1']['wwmax'], noise)
    if version == 'dasilva':
        flux_r2, var_fluxr2 = compute_flux(e2ds, blaze, wfunc, ww,
                                           windows['wr2']['wwmin'],
                                           windows['wr2']['wwmax'], noise)
    else:
        flux_r2, var_fluxr2 = flux_r1, 0.0

    # Compute Halpha index
    h = 2*flux_h/(flux_r1 + flux_r2)

    # Compute error
    var_h = 4*var_fluxh/(flux_r1 + flux_r2)**2 + \
        h**2/(flux_r1 + flux_r2)**2 * (var_fluxr1 + var_fluxr2)

    if full_output:
        return h, sqrt(var_h), np.array([[flux_h, var_fluxh],
                                        [flux_r1, var_fluxr1],
                                        [flux_r2, var_fluxr2]])
    else:
        return h, sqrt(var_h)


__author__ = 'Rodrigo F. Diaz'
