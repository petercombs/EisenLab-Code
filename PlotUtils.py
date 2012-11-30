from matplotlib import pyplot as mpl
from matplotlib.colors import hsv_to_rgb
from matplotlib import cm
from scipy.stats import gaussian_kde
from numpy import log, array, Inf, median, exp, argsort, linspace
import numpy as np


def scatter_heat(x, y, **kwargs):
    if 's' not in kwargs:
        kwargs['s'] = 10
    if 'edgecolors' not in kwargs:
        kwargs['edgecolors'] = 'none'
    if 'cmap' not in kwargs:
        kwargs['cmap'] = cm.jet
    estimator = gaussian_kde([x, y])
    density = estimator.evaluate([x, y])

    normdensity = exp(density.clip(median(density), Inf))

    xlim = kwargs.pop('xlim', (min(x), max(x)))
    ylim = kwargs.pop('ylim', (min(y), max(y)))
    retval = mpl.scatter(x, y, c=normdensity, **kwargs)
    ax = mpl.gca()
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    return retval

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


def hist_sorted(*args, **kwargs):
    all_ns = []
    all_patches = []

    labels = kwargs.pop('labels', None)
    if not labels:
        labels = ['data %d' % (i+1) for i in range(len(args))]
    elif len(labels) != len(args):
        raise ValueError('length of labels not equal to length of data')

    bins = kwargs.pop('bins', linspace(min(min(a) for a in args),
                                       max(max(a) for a in args),
                                       num = 11))

    for data, label in zip(args, labels):
        ns, bins, patches = mpl.hist(data, bins=bins, label=label, **kwargs)
        all_ns.append(ns)
        all_patches.append(patches)
    z_orders = -argsort(all_ns, axis=0)

    for zrow, patchrow in zip(z_orders, all_patches):
        assert len(zrow) == len(patchrow)
        for z_val, patch in zip(zrow, patchrow):
            patch.set_zorder(z_val)

    return all_ns, bins, all_patches




def plot_likelihoods(likelihoods, starts, column_headers):
    n_samples = len(column_headers)
    max_val = np.argmax(starts>150)
    print max_val
    plots = []
    for i in range(n_samples):
        hsv = np.array([0.7*i/n_samples, 1, 1])
        color = tuple(hsv_to_rgb(np.reshape(hsv, (1,1,3))))[0].flatten()
        print color
        plots.append(mpl.plot(starts[:max_val], likelihoods[i,:max_val],
                              label=column_headers[i], color=color))
        best = np.argmax(likelihoods[i,:])
        print best
        plots.append(mpl.plot(starts[best], likelihoods[i,best], '*',
                          color=color))
    return plots


