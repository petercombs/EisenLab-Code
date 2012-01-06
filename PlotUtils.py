from matplotlib import pyplot as mpl
from matplotlib import cm
from collections import defaultdict
from scipy.stats import spearmanr, pearsonr, gaussian_kde
from numpy import log, array, Inf, median, exp, argsort


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


def hist_sorted(*args, labels=None, **kwargs):
    all_ns = []
    all_patches = []

    if not labels:
        labels = ['data %d' % (i+1) for i in range(len(args))]
    elif len(labels) != len(args):
        raise ValueError('length of labels not equal to length of data')

    for data, label in zip(args, labels):
        ns, bins, patches = mpl.hist(data, label=label, **kwargs)
        all_ns.append(ns)
        all_patches.append(patches)
    z_orders = argsort(all_ns, axis=0)

    for zrow, patchrow in zip(z_orders, all_patches):
        assert len(zrow) == len(patchrow)
        for z_val, patch in zip(zrow, patchrow):
            patch.set_zorder(z_val)

    return all_ns, bins, all_patches





