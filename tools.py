import numpy as np
import pylab as plt


def integrate(x, y, weight, xmin, xmax, vary=None):
    """
    Return integrated flux of array y (already integrated in
    each pixel) with respecto to x, with variance vary and
    weight function weigth between limits xmin and xmax.
    """

    # Keep only elements of x, y, and weight within the interval,
    # include fractional pixels
    deltax = np.diff(x)
    # Reduce x of first element for compatibility with deltax
    x = x[1:]

    cond = (x > xmin - deltax/2.0) * (x <= xmax + deltax/2.0)

    # Compute fraction of pixel within interval
    fraction_left = 0.5 - (xmin - x[cond])/deltax[cond]
    fraction_right = 0.5 - (x[cond] - xmax)/deltax[cond]

    fraction_left = np.where(fraction_left > 1, 1.0, fraction_left)
    fraction_right = np.where(fraction_right > 1, 1.0, fraction_right)

    fraction = np.minimum(fraction_left, fraction_right)

    # Sum contributions of pixels inside interval, considering fractions
    # when necessary
    summed_y = np.sum(y[1:][cond] * weight[1:][cond] * fraction)
    summed_weight = np.sum(weight[1:][cond] * fraction)

    integral = np.divide(summed_y, summed_weight)

    if vary is not None:
        # Also compute error
        summed_var = np.sum(vary[1:][cond] * (weight[1:][cond] * fraction)**2)
        var_integral = np.divide(summed_var, summed_weight**2)

    else:
        var_integral = None

    return integral, var_integral


def plot_windows(ww, y, windows):

    ax = plt.subplots(1, 1)[1]

    for win in windows:
        indices = np.where(np.logical_and(ww > windows[win]['wwmin'],
                                          ww < windows[win]['wwmax']))

        for order in np.unique(indices[0]):

            ax.plot(ww[order], y[order], 'k')

            for ind in enumerate(indices):
                ax.plot(ww[ind], y[ind], 'r', lw=2)
    return

__author__ = 'Rodrigo F. Diaz'
