import os
import numpy as np

import activist.data
from activist.tools import load_rdbfile
import activist.indexes as indexes

__author__ = 'Rodrigo F. Diaz'


def get_indexes(target, instrument, innights, outdir=os.getenv('HOME')):
    """
    Compute activity indexes for a given star and target.

    :param str target:
    :param str instrument:
    :param iterable innights: nights where to compute indexes
    :param str outdir: directory where output files are to be saved
    :return:
    """
    # Read data
    try:
        basedir = activist.data.datadir[instrument].replace('reduced', 'dbase')
    except KeyError:
        raise NameError('Instrument not recognized.')

    if 'coralie' in instrument:
        inststr = 'coralie'
    else:
        inststr = instrument

    try:
        data = load_rdbfile(os.path.join(basedir,
                                         '{}_{}_extract.rdb'.format(target,
                                                                    inststr)))
    except IOError:
        print('No {} data found for star {}.'.format(instrument, target))
        return 0

    nights = data['night']
    rootnames = data['file_root']

    # Select nights
    if innights != 'all':
        nights = nights[innights]
        rootnames = rootnames[innights]

    # Initialize arrays
    halpha = np.empty(len(nights))
    sigma_halpha = np.empty(len(nights))
    fluxes_halpha = np.empty((len(nights), 3, 2))

    sodium = np.empty(len(nights))
    sigma_sodium = np.empty(len(nights))
    fluxes_sodium = np.empty((len(nights), 4, 2))

    calcium = np.empty(len(nights))
    sigma_calcium = np.empty(len(nights))
    fluxes_calcium = np.empty((len(nights), 4, 2))

    noncomputed = []

    for i in range(len(nights)):
        if i % 10 == 0:
            print('%.2f%%' % (100*float(i)/len(nights)))

        try:
            ww, e2ds, blaze, noise = activist.data.read_data(rootnames[i],
                                                             nights[i])

            halpha[i], sigma_halpha[i], \
                fluxes_halpha[i] = indexes.Halpha(ww, e2ds, blaze, noise,
                                                  full_output=True,
                                                  version='cincunegui')

            sodium[i], sigma_sodium[i], \
                fluxes_sodium[i] = indexes.NaID(ww, e2ds, blaze, noise,
                                                full_output=True)

            calcium[i], sigma_calcium[i], \
                fluxes_calcium[i] = indexes.CaII(ww, e2ds, blaze, noise,
                                                 full_output=True)

        except:
            halpha[i], sigma_halpha[i] = 999.9, 999.9
            sodium[i], sigma_sodium[i] = 999.9, 999.9
            calcium[i], sigma_calcium[i] = 999.9, 999.9

            noncomputed.append(rootnames[i])

    f = open(os.path.join(outdir, '{}_indexes.txt'.format(target)), 'w')
    f.write('bjd\tSindex\tsig_Sindex\tHindex\tsig_Hindex\tDindex\tsig_Dindex\n')
    f.write('---\t------\t----------\t------\t----------\t------\t'
            '-----------\n')

    for i in range(len(sodium)):
        f.write('{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\t'
                '{:.6f}\n'.format(data['jdb'][i], calcium[i], sigma_calcium[i],
                                  halpha[i], sigma_halpha[i], sodium[i],
                                  sigma_sodium[i]))
    f.close()

    f = open(os.path.join(outdir, '{}_noncomputed.txt'.format(target)), 'w')
    for i in range(len(noncomputed)):
        f.write(noncomputed+'\n')
    f.close()