import os
import glob

import numpy as np
import pyfits

# DACEi3
datadir = {'harps': '/dace/data/harps/DRS-3.5/reduced',
           'sophie': '/dace/data/sophie/DRS-3.1/reduced',
           'coralie98': '/dace/data/coralie98/DRS-3.3/reduced',
           'coralie07': '/dace/data/coralie07/DRS-3.4/reduced',
           'coralie14': '/dace/data/coralie14/DRS-3.8/reduced'
           }

# Blaze dir (for tests)
# blazedir = '/home/spectro/diazr/DACE/blazefiles'
blazedir = datadir


# Relevant physical constants
lightspeed = 299792458.0*1e-3  # km/s


def read_e2ds(rootname, night, instrument, header=False):
    """
    Read e2ds_A.fits file.

    :param str rootname: name of file without e2ds termination.
    :param str night: night of observation. Format YYYY-MM-DD
    :param str instrument: instrument must exist in datadir dictionary.
    :param bool header: determines if header instance is returned.
    :return: e2ds data array (in e-), format (k x n), where k is number of
     orders and n is the number of pixels per order. If header is True,
     then the header instance is also returned.
    """

    e2dsname = os.path.join(night, rootname+'_e2ds_A.fits')
    try:
        dd = datadir[instrument]
    except KeyError:
        print('Instrument not recognised.')
        raise
    try:
        ar = pyfits.open(os.path.join(dd, e2dsname))
    except IOError:
        print('File not found.')
        raise

    if header:
        return ar[0].data, ar[0].header
    else:
        return ar[0].data


def read_s1d(rootname, night, instrument, header=False):
    """
    Read e2ds_A.fits file given a rootname (without the e2ds termination).
    """
    s1dname = os.path.join(night, rootname+'_s1d_A.fits')
    ar = pyfits.open(os.path.join(datadir[instrument], s1dname))

    if header:
        return ar[0].data, ar[0].header
    else:
        return ar[0].data


def read_blaze(e2dshead, night, instrument):
    """
    Read the blaze function for a given rootname.
    """
    if instrument == 'harps' or 'coralie' in instrument:
        obs = 'ESO'
        # Get blaze name from e2ds header
        blazename = e2dshead['HIERARCH {} DRS BLAZE FILE'.format(obs)]

    elif instrument == 'sophie':
        raise NotImplemented('SOPHIE not yet implemented.')

    else:
        raise NameError('Unrecognized instrument.')

    blaze = pyfits.getdata(os.path.join(datadir[instrument], night, blazename))
    return blaze


def wavelength(e2ds, header, instrument):
    """
    Compute the wavelength solution based on the CCF keywords.
    """
    if instrument == 'harps' or 'coralie' in instrument:
        obs = 'ESO'
    elif instrument == 'sophie':
        obs = 'OHP'
    else:
        raise NameError('Unrecognized instrument.')

    berv = header['HIERARCH {} DRS BERV'.format(obs)]
    deg = header['HIERARCH {} DRS CAL TH DEG LL'.format(obs)]

    ll_coeff = np.zeros((len(e2ds), deg + 1))

    # Read coefficients
    for i in range(len(e2ds)):
        for j in range(deg + 1):
            ll_coeff[i, j] = header['HIERARCH {} DRS CAL TH COEFF '
                                    'LL{}'.format(obs, (j + (deg + 1)*i))]

    # Evaluate polynomials
    x = np.arange(e2ds.shape[1])  # Pixel array
    ww = np.zeros(e2ds.shape)  # Wavelength 2D array
    ##
    for i in range(len(ww)):
        ww[i] = np.poly1d(ll_coeff[i][::-1])(x)

    # Correct for BERV
    ww /= 1 - berv / lightspeed

    return ww


def wavelength_s1d(header, instrument):
    if instrument == 'harps' or 'coralie' in instrument:
        obs = 'ESO'
    elif instrument == 'sophie':
        obs = 'OHP'
    else:
        raise NameError('Unrecognized instrument.')

    berv = header['HIERARCH {} DRS BERV'.format(obs)]

    ww = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    return ww/(1 - berv/lightspeed)


def ron(header, instrument):
    if instrument == 'harps' or 'coralie' in instrument:
        obs = 'ESO'
    elif instrument == 'sophie':
        obs = 'OHP'
    else:
        raise NameError('Unrecognized instrument.')

    return header['HIERARCH {} DRS CCD SIGDET'.format(obs)]


def get_gain(header, instrument):
    if instrument == 'harps' or 'coralie' in instrument:
        obs = 'ESO'
    elif instrument == 'sophie':
        obs = 'OHP'
    else:
        raise NameError('Unrecognized instrument.')

    return header['HIERARCH {} DRS CCD CONAD'.format(obs)]


def velocity(rootname, night, instrument):
    """
    Get the measured radial velocity for a given observation.
    """

    # Get instrument
    if instrument is None:
        if 'HARPS' in rootname:
            obs = 'ESO'
            instrument = 'harps'
        elif 'SOPHIE' in rootname:
            obs = 'SOPHIE'
            instrument = 'sophie'
        elif 'CORALIE' in rootname:
            obs = 'ESO'
            instrument = 'coralie'

    elif instrument in datadir:
        pass

    else:
        raise NameError('Unrecognized instrument.')

    if instrument == 'harps' or 'coralie' in instrument:
        obs = 'ESO'
    elif instrument == 'sophie':
        obs = 'OHP'
    else:
        raise NameError('Unrecognized instrument.')

    # List all CCF file for that night and rootname
    ccfl = glob.glob(os.path.join(datadir[instrument], night,
                                  rootname+'*ccf*A.fits'))
    return pyfits.getheader(ccfl[0])['HIERARCH {} DRS CCF RVC'.format(obs)]


def read_data(rootname, night, instrument=None, s1d=False):

    # Get instrument
    if instrument is None:
        if 'HARPS' in rootname:
            instrument = 'harps'
        elif 'SOPHIE' in rootname:
            instrument = 'sophie'
        elif 'CORALIE' in rootname:
            instrument = 'coralie07'
    elif instrument in datadir:
        pass
    else:
        raise NameError('Unrecognized instrument.')
            
    if s1d:
        s1d, header = read_s1d(rootname, night, instrument, header=True)

        # Get wavelength solution for s1d (already corrected for BERV)
        ww = wavelength_s1d(header, instrument)

    else:
        e2ds, header = read_e2ds(rootname, night, instrument, header=True)
        blaze = read_blaze(header, night, instrument)

        # Get wavelength solution for e2ds (already corrected for BERV)
        ww = wavelength(e2ds, header, instrument)

    # Go to rest frame
    rv = velocity(rootname, night, instrument)
    ww /= 1 + rv / lightspeed

    # Get Read out noise (in e-)
    readout_noise = ron(header, instrument)

    if s1d:
        return ww, s1d
    else:
        return ww, e2ds, blaze, readout_noise
