from math import *

import numpy as np

from indexes import indparams
from tools import integrate


def compute_flux(e2ds, blaze, weight, ww, wwmin, wwmax, noise=0.0):
    """
    Compute mean flux in e2ds spectrum between a given wavelngth range.

    :param np.array e2ds: 2-D spectrum; must have shape (norders, npixels).
    :param np.array blaze: 2-D blaze function; same shape as e2ds.
    :param np.array weight: a 2-D weigth function; same shape as e2ds.
    :param np.array ww: 2-D wavelength array; same shape as e2ds.
    :param float wwmin: lowelimit of integration window.
    :param float wwmax: upper limit of integration window.
    :param float or np.array noise: additional noise in e-. Must be a float, in
    which case the same noise is added to all pixels, or have the same shape as
    e2ds.

    :return: average flux and variance in window.
    """

    # Find the spectral order(s) where the relevant window is located
    orders = np.unique(np.where(np.logical_and(ww > wwmin, ww < wwmax))[0])

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
    ff = np.sum(flux/vflux)/np.sum([1.0/vflux])
    vv = 1.0/np.sum([1.0/vflux])

    return ff, vv


def CaII(ww, e2ds, blaze, noise=0.0, full_output=False):
    """
    Compute the Ca II H and K lines activity index (S) using the blaze function
    as weight and triangular windows around the line cores.

    :param np.array ww: wavelength array. Shape is (k x n), where k is number of
    echelle orders and n is the pixel size of the orders.
    :param np.array e2ds: flux array. Shape is (k x n).
    :param np.array blaze: blaze function of spectrum. Shape is (k x n). This is
    used to flatten the input spectrum and as weight for the mean flux
    computation.
    :param float noise: a level of extra noise, in e-
    :param bool full_output: if True, return also the measurements of the
    mean flux on each window.
    :return:
    """

    # Use blaze as weight function.
    wfunc = blaze

    # Select indices
    windows = indparams['CaII'].copy()

    # Compute wavelngth limits for computation
    for win in ('w1', 'w2', 'wr1', 'wr2'):
        windows[win]['wwmin'] = windows[win]['ww0'] - windows[win]['dww']*0.5
        windows[win]['wwmax'] = windows[win]['ww0'] + windows[win]['dww']*0.5

    # Compute fluxes for each reference window
    flux_r1, var_fluxr1 = compute_flux(e2ds, blaze, np.ones_like(ww), ww,
                                       windows['wr1']['wwmin'],
                                       windows['wr1']['wwmax'], noise)

    flux_r2, var_fluxr2 = compute_flux(e2ds, blaze, np.ones_like(ww), ww,
                                       windows['wr2']['wwmin'],
                                       windows['wr2']['wwmax'], noise)

    # Compute line fluxes

    #  Add triangular weight function for K line
    triang1 = np.where(ww >= windows['w1']['ww0'],
                       -(ww - windows['w1']['ww0'])/1.09 + 1.0,
                       (ww - windows['w1']['ww0'])/1.09 + 1.0)

    triang1 = np.where(triang1 > 0, triang1, triang1*0.0)

    flux_k, var_fluxk = compute_flux(e2ds, blaze, wfunc * triang1, ww,
                                     windows['w1']['wwmin'],
                                     windows['w1']['wwmax'], noise)

    # ## Add triangular weight function for H line
    triang2 = np.where(ww >= windows['w2']['ww0'],
                       -(ww - windows['w2']['ww0'])/1.09 + 1.0,
                       (ww - windows['w2']['ww0'])/1.09 + 1.0)
    triang2 = np.where(triang2 > 0, wfunc, triang2*0.0)

    flux_h, var_fluxh = compute_flux(e2ds, blaze, wfunc * triang2, ww,
                                     windows['w2']['wwmin'],
                                     windows['w2']['wwmax'], noise)

    # Compute Ca II index
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


def NaID(ww, e2ds, blaze, noise=0, full_output=False, version='classic'):

    # Use blaze as weight function.
    wfunc = blaze

    # Select indices
    windows = indparams['NaID_{}'.format(version)].copy()

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
    flux_r2, var_fluxr2 = compute_flux(e2ds, blaze, wfunc, ww,
                                       windows['wr2']['wwmin'],
                                       windows['wr2']['wwmax'], noise)

    if version == 'classic':
        flux_r1, var_fluxr1 = compute_flux(e2ds, blaze, wfunc, ww,
                                           windows['wr1']['wwmin'],
                                           windows['wr1']['wwmax'], noise)
    elif version == 'telluric':
        flux_r1, var_fluxr1 = flux_r2, 0.0

    else:
        raise NameError('version not recognised.')

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


def Halpha(ww, e2ds, blaze, noise=0, full_output=False, version='narrow'):

    # Use blaze as weight function.
    wfunc = blaze

    # Select indices
    windows = indparams['Halpha_{}'.format(version)].copy()

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
    elif version == 'narrow' or version == 'broad':
        flux_r2, var_fluxr2 = flux_r1, 0.0

    else:
        raise NameError('version not recognised.')

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
