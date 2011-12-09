from matplotlib import pyplot as mpl
from matplotlib import cm
from collections import defaultdict
from scipy.stats import spearmanr, pearsonr, gaussian_kde
from numpy import log, array, Inf, median, exp


def loglog_heat(x,y, **kwargs):
    if 's' not in kwargs:
        kwargs['s'] = 10
    if 'edgecolors' not in kwargs:
        kwargs['edgecolors'] = 'none'
    if 'cmap' not in kwargs:
        kwargs['cmap'] = cm.jet
    logx = log(array(x))
    logy = log(array(y))
    estimator = gaussian_kde([logx, logy])
    density = estimator.evaluate([logx, logy])

    normdensity = exp(density.clip(median(density), Inf))

    xlim = kwargs.pop('xlim', (min(x), max(x)))
    ylim = kwargs.pop('ylim', (min(y), max(y)))
    retval = mpl.scatter(x, y, c=normdensity, **kwargs)
    ax = mpl.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    return retval



