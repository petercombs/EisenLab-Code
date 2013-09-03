from matplotlib import pyplot as mpl
from matplotlib.colors import hsv_to_rgb
from matplotlib import cm
from scipy.stats import gaussian_kde
from numpy import log, array, Inf, median, exp, argsort, linspace
import numpy as np

import urllib, time
from os import path
def imget(imname):
    """ Use cached, or fetch an image from FlyExpress

Assumes that the image name is one from BDGP, in which case
it's pretty easy to look at the source of the FlyExpress
report pages and see what the format is.

"""
    im_basename = path.splitext(path.basename(imname))[0]
    filename = path.join('figures', 'BDGP', im_basename+'.bmp')
    if not path.exists(filename):
        base_web = "http://www.flyexpress.net/ZOOX4_DBImages/BDGP/thumbnails/%s_s.bmp"
        print "1 second delay to avoid spamming server"
        time.sleep(1)
        urllib.urlretrieve(base_web % im_basename, filename)
    return mpl.imread(filename)

def scatter_heat(x, y, **kwargs):
    if 's' not in kwargs:
        kwargs['s'] = 10
    if 'edgecolors' not in kwargs:
        kwargs['edgecolors'] = 'none'
    if 'cmap' not in kwargs:
        kwargs['cmap'] = cm.jet
    if 'density' not in kwargs:
        estimator = gaussian_kde([x, y])
        density = estimator.evaluate([x, y])
    else:
        density = kwargs['density']

    normdensity = exp(density.clip(median(density), Inf))

    xlim = kwargs.pop('xlim', (min(x), max(x)))
    ylim = kwargs.pop('ylim', (min(y), max(y)))
    retval = mpl.scatter(x, y, c=normdensity, **kwargs)
    ax = mpl.gca()
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    return retval, density

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


def svg_heatmap(data, filename, row_labels=None, boxsize=4,
               cmap=cm.Blues, norm_rows_by = None, draw_row_labels=False,
               col_sep = ''):
    """
    """
    import svgwrite as svg
    import pandas as pd

    dwg = svg.Drawing(filename)
    if not isinstance(data, tuple):
        data = (data,)

    rows, cols = np.shape(data[0])
    if row_labels is None:
        if hasattr(data[0], 'index'):
            row_labels = data[0].index
        else:
            row_labels = ['' for row in range(rows)]



    if not hasattr(cmap, "__len__"):
        cmap = [cmap for frame in data]

    x_start = 0
    for frame, c_cmap in zip(data, cmap):
        frame = pd.DataFrame(frame)
        if norm_rows_by is None:
            norm_data = frame.copy()
        elif norm_rows_by is 'mean':
            norm_data = frame.divide(frame.mean(axis=1), axis=0)
        elif norm_rows_by is 'max':
            norm_data = frame.divide(frame.max(axis=1), axis=0)
        elif hasattr(norm_rows_by, "__len__") and len(norm_rows_by) == rows:
            norm_data = frame.divide(norm_rows_by, axis=0)

        else:
            raise TypeError("norm_rows_by should be the same shape "
                            "as the number of rows")
        new_rows, new_cols = np.shape(frame)
        if hasattr(frame, 'index'):
            col_labels = frame.columns
        else:
            col_labels = ['' for col in range(new_cols)]
        if new_rows != rows:
            raise ValueError("All input elements must have the same number of"
                             " rows (and same row meanings --unchecked)")

        for i in range(rows):
            prefix = col_labels[0][:col_labels[0].find(col_sep)]
            for j in range(new_cols):
                g = dwg.g()
                g.add(svg.base.Title("{}, {}: {:.2f}".format(row_labels[i],
                                                             col_labels[j],
                                                             frame.ix[i,j])))
                g.add(dwg.rect((x_start + boxsize*j, i*boxsize),
                               (boxsize, boxsize),
                               style="fill:#{:02x}{:02x}{:02x}"
                                .format(*[int(255*x) for x in
                                          c_cmap(norm_data.ix[i,j])])))
                dwg.add(g)
                col_base = col_labels[j][:col_labels[j].find(col_sep)] 
                if col_base != prefix:
                    prefix = col_base
                    g.add(dwg.line((x_start+boxsize*j, i*boxsize),
                                   (x_start+boxsize*j, (i+1)*boxsize),
                                   style="stroke-width:{}; stroke:#000000"
                                   .format(.1 * boxsize)))
        x_start += new_cols * boxsize + boxsize


    if draw_row_labels:
        for i in range(rows):
            dwg.add(dwg.text(row_labels[i], (x_start, i*boxsize+boxsize),))
    dwg.saveas(filename)
