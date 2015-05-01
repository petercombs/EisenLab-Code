from __future__ import division, print_function
from matplotlib import pyplot as mpl
from matplotlib.colors import hsv_to_rgb, LinearSegmentedColormap
from matplotlib import cm
from scipy.stats import gaussian_kde
from numpy import log, array, Inf, median, exp, argsort, linspace, isfinite
from itertools import repeat
import numpy as np
import subprocess

import urllib
import time
from os import path

ISH_ROT_4 = hsv_to_rgb(array(
    [[[(0.65+offset)%1, 0.00, 1.00],
      [(0.65+offset)%1, 0.53, 1.00],
      [(0.65+offset)%1, 0.53, 0.38],]
     for offset in linspace(0, 1, 4, endpoint=False)
    ]))
ISH_ROT_5 = hsv_to_rgb(array(
    [[[(0.65+offset)%1, 0.00, 1.00],
      [(0.65+offset)%1, 0.53, 1.00],
      [(0.65+offset)%1, 0.53, 0.38],]
     for offset in linspace(0, 1, 5, endpoint=False)
    ]))
ISH_ROT_6 = hsv_to_rgb(array(
    [[[(0.65+offset)%1, 0.00, 1.00],
      [(0.65+offset)%1, 0.53, 1.00],
      [(0.65+offset)%1, 0.53, 0.38],]
     for offset in linspace(0, 1, 6, endpoint=False)
    ]))

ISH_CMS_4 = []
ISH_CMS_5 = []
ISH_CMS_6 = []

for CMS, ROT in [(ISH_CMS_4, ISH_ROT_4),
                 (ISH_CMS_5, ISH_ROT_5),
                 (ISH_CMS_6, ISH_ROT_6)]:
    for I, ARR in enumerate(ROT):
        CMS.append(
            LinearSegmentedColormap('ish{}'.format(I),
                                    dict(red=((0.0, ARR[0, 0], ARR[0, 0]),
                                              (0.7, ARR[1, 0], ARR[1, 0]),
                                              (1.0, ARR[2, 0], ARR[2, 0])),
                                         green=((0.0, ARR[0, 1], ARR[0, 1]),
                                                (0.7, ARR[1, 1], ARR[1, 1]),
                                                (1.0, ARR[2, 1], ARR[2, 1])),
                                         blue=((0.0, ARR[0, 2], ARR[0, 2]),
                                               (0.7, ARR[1, 2], ARR[1, 2]),
                                               (1.0, ARR[2, 2], ARR[2, 2])),
                                        )))


ISH = LinearSegmentedColormap('ish',
                              dict(red=((0, 1, 1),
                                        (.7, 120/255, 120/255),
                                        (1, 46/255, 46/255)),
                                   green=((0, 1, 1),
                                          (.7, 129/255, 129/255),
                                          (1, 46/255, 46/255)),
                                   blue=((0, 1, 1),
                                         (.7, 1, 1),
                                         (1, 98/255, 98/255))))


def imget(imname):
    """ Use cached, or fetch an image from FlyExpress

Assumes that the image name is one from BDGP, in which case
it's pretty easy to look at the source of the FlyExpress
report pages and see what the format is.

"""
    im_basename = path.splitext(path.basename(imname))[0]
    filename = path.join('figures', 'BDGP', im_basename+'.bmp')
    if not path.exists(filename):
        base_web = ("http://www.flyexpress.net/"
                    "ZOOX4_DBImages/BDGP/thumbnails/%s_s.bmp")
        print("1 second delay to avoid spamming server")
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


def loglog_heat(x, y, **kwargs):
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
                                       num=11))

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
    max_val = np.argmax(starts > 150)
    print(max_val)
    plots = []
    for i in range(n_samples):
        hsv = np.array([0.7*i/n_samples, 1, 1])
        color = tuple(hsv_to_rgb(np.reshape(hsv, (1, 1, 3))))[0].flatten()
        print(color)
        plots.append(mpl.plot(starts[:max_val], likelihoods[i, :max_val],
                              label=column_headers[i], color=color))
        best = np.argmax(likelihoods[i, :])
        print(best)
        plots.append(mpl.plot(starts[best], likelihoods[i, best], '*',
                              color=color))
    return plots


def svg_heatmap(data, filename, row_labels=None, box_size=4,
                index=None,
                cmap=ISH, norm_rows_by=None, draw_row_labels=False,
                color_row_labels=False,
                col_sep='', box_height=None, total_width=None,
                draw_box=False, draw_name=False, data_names=None,
                make_hyperlinks = False,
                progress_bar = False,
                max_width=np.inf,
                spacers=None,
                convert=False,
                cmap_by_prefix=None,
                draw_average=False,
                draw_average_only=False,
                average_scale=1,
                split_columns=False,
                vspacer=30,
                hatch_nan=True, hatch_size=20,
                first_col='', last_col=''):
    """
    Draw heatmap as an SVG file stored in filename

    *data* can be either a 2D array-like type (list of lists, numpy array,
    pandas DataFrame, etc), or a tuple of 2D array-likes, in which case a
    separator will be added between each one in the output

    *cmap* is a matplotlib-like colormap (i.e. a callable that expects floats
    in the range 0.0-1.0.), or an iterable of the same length as the tuple
    *data* containing colormaps

    *row_labels* can be supplied, otherwise they will detected from the first
    item in *data*, if available, and if not they will be blank.

    If *total_width* is supplied, width of each dataset in *data* will be
    scaled to that constant. If *box_height* is supplied, the height of each
    row will be *box_height*, otherwise it will be equal to the width of each
    element. If neither are supplied, elements will be squares equal to
    *box_size*. IT IS STRONGLY RECOMMENDED that if if supplying *total_width*,
    *box_height* also be specified, but this is not enforced.

    *draw_row_labels*, if True, will label the rows on the right hand side. As
    of 2013/09/03, this won't scale the SVG properly, so including the
    resulting file in an html element won't display properly.

    *spacers* is the distance between adjacent datasets.  Can either be a
    number, in which case it will apply to all datasets, or an interable for
    different distances. If the iterable is shorter than the number of
    datasets, the last value will be repeated.

    """
    import svgwrite as svg
    import pandas as pd

    if split_columns and isinstance(data, pd.DataFrame):
        from Utils import sel_startswith
        colnames = list(sorted(
            {col.split(col_sep)[0] for col in data.columns}))
        data = tuple(
            data.select(**sel_startswith(colname)) for colname in colnames
        )
    elif not isinstance(data, tuple):
        data = (data,)

    rows, cols = np.shape(data[0])
    if index is not None:
        rows = len(index)
    if box_height is None:
        box_height = box_size

    if total_width is not None and max_width is not np.inf:
        dwg = svg.Drawing(filename,
                          size=(max_width,
                                np.ceil((len(data) * total_width)/max_width +
                                        (draw_average or draw_average_only))
                                * (box_height+vspacer)))
    else:
        dwg = svg.Drawing(filename)
    dwg.add(svg.base.Title(path.basename(filename)))

    pat = dwg.pattern(id='hatch', insert=(0, 0), size=(hatch_size, hatch_size),
                      patternUnits='userSpaceOnUse')
    g = pat.add(dwg.g(style="fill:none; stroke:#B0B0B0; stroke-width:1"))
    g.add(dwg.path(('M0,0', 'l{hatch},{hatch}'.format(hatch=hatch_size))))
    g.add(dwg.path(('M{hatch2},0 l{hatch2},{hatch2}'.format(hatch2=hatch_size/2).split())))
    g.add(dwg.path(('M0,{hatch2} l{hatch2},{hatch2}'.format(hatch2=hatch_size/2).split())))

    dwg.add(pat)

    if row_labels is None:
        if index is not None:
            row_labels = index
        elif hasattr(data[0], 'index'):
            row_labels = data[0].index
        else:
            row_labels = ['' for row in range(rows)]

    if box_height is None:
        box_height = box_size

    if not hasattr(cmap, "__len__"):
        cmap = [cmap for frame in data]

    if data_names is None:
        data_names = ["" for frame in data]

    if len(cmap) != len(data):
        raise ValueError("cmap and data should be the same length")

    if not hasattr(spacers, "__len__"):
        spacers = [spacers]
    else:
        spacers = list(spacers)
    while len(spacers) < len(data):
        spacers.append(spacers[-1])

    if not isinstance(norm_rows_by, tuple):
        norm_rows_by = repeat(norm_rows_by)

    if 'center0all' in norm_rows_by:
        all_data = pd.concat(data, axis=1)

    x_start = 0
    y_start = 0
    y_diff = 0
    if progress_bar:
        from progressbar import ProgressBar
        iterator = zip(data, cmap, data_names, norm_rows_by, spacers)
        pbar = ProgressBar(maxval=len(iterator)*rows).start()
        pbar_val = 0
    else:
        iterator = zip(data, cmap, data_names, norm_rows_by, spacers)

    for frame, c_cmap, name, normer, spacer in iterator:
        if frame is None:
            if total_width is not None:
                if spacer is None:
                    x_start += total_width * 1.1
                else:
                    x_start += total_width + spacer
            else:
                if spacer is None:
                    x_start += box_size
                else:
                    x_start += spacer
            if x_start > max_width:
                x_start = 0
                y_start += y_diff
            continue
        frame = pd.DataFrame(frame)
        if index is not None:
            frame = frame.ix[index]
        if normer is None:
            norm_data = frame.copy()
        elif normer is 'mean':
            norm_data = frame.divide(frame.dropna(axis=1, how='all').mean(axis=1)+10, axis=0)
        elif normer is 'max':
            norm_data = frame.divide(frame.dropna(axis=1, how='all').max(axis=1)+10, axis=0)
        elif normer is 'center0':
            norm_data = (0.5 +
                         0.5 * frame.divide(frame.dropna(axis=1).abs().max(axis=1),
                                      axis=0)
                        )
        elif normer is 'center0all':
            norm_data = (0.5 +
                         0.5 *
                         frame.divide(all_data.dropna(how='all', axis=1).abs().max(axis=1),
                                      axis=0)
                        )
        elif index is not None and hasattr(normer, "ix"):
            norm_data = frame.divide(normer.ix[index], axis=0)
        elif hasattr(normer, "__len__") and len(normer) == rows:
            norm_data = frame.divide(normer, axis=0)

        elif hasattr(normer, "__len__"):
            raise TypeError("norm_rows_by should be the same shape "
                            "as the number of rows")
        else:
            norm_data = frame.divide(normer, axis=0)

        if not c_cmap or str(c_cmap).lower() == 'default':
            c_cmap = ISH

        new_rows, new_cols = np.shape(frame)
        if hasattr(frame, 'index'):
            col_labels = frame.columns
        else:
            col_labels = ['' for col in range(new_cols)]
        if new_rows != rows:
            raise ValueError("All input elements must have the same number of"
                             " rows (and same row meanings --unchecked)")

        if total_width is not None:
            box_size = total_width / float(new_cols)

        i = 0
        if not draw_average_only:
            for i in range(rows):
                if progress_bar:
                    pbar.update(pbar_val)
                    pbar_val += 1
                prefix = col_labels[0][:col_labels[0].find(col_sep)]
                if cmap_by_prefix:
                    c_cmap = cmap_by_prefix(prefix)
                for j in range(new_cols):
                    g = dwg.g()
                    g.add(svg.base.Title("{}, {}: {:.2f}".format(row_labels[i],
                                                                 col_labels[j],
                                                                 frame.ix[i, j])))
                    hatch = not isfinite(norm_data.ix[i, j])
                    if hatch:
                        n = 0
                        norm_data.ix[i, j] = 0
                        if j > 0 and isfinite(norm_data.ix[i,j-1]):
                            norm_data.ix[i, j] += norm_data.ix[i, j-1]
                            n += 1
                        if (j + 1 < len(norm_data.columns)
                            and isfinite(norm_data.ix[i, j+1])):
                            norm_data.ix[i, j] += norm_data.ix[i, j+1]
                            n += 1
                        norm_data.ix[i, j] /= n
                    g.add(dwg.rect((x_start + box_size*j, y_start + i*box_height),
                                   (box_size, box_height),
                                   style="fill:#{:02x}{:02x}{:02x}"
                                   .format(*[int(255*x) for x in
                                             c_cmap(norm_data.ix[i, j])])))
                    dwg.add(g)
                    if hatch_nan and hatch:
                        g.add(dwg.rect((x_start + box_size*j,
                                        y_start + i*box_height),
                                       (box_size, box_height),
                                       style="fill:url(#hatch)"
                                      )
                             )
                    col_base = col_labels[j][:col_labels[j].find(col_sep)]
                    if col_base != prefix:
                        prefix = col_base
                        if cmap_by_prefix:
                            c_cmap = cmap_by_prefix(prefix)
                        g.add(dwg.line((x_start + box_size * j,
                                        y_start + i * box_height),
                                       (x_start + box_size * j,
                                        y_start + (i + 1) * box_height),
                                       style="stroke-width:{}; stroke:#000000"
                                       .format(.1 * box_size)))
        else:
            for j in range(new_cols):
                hatch = not isfinite(norm_data.ix[0, j])
                if hatch:
                    n = 0
                    norm_data.ix[:, j] = 0
                    if j > 0 and isfinite(norm_data.ix[0,j-1]):
                        norm_data.ix[:, j] += norm_data.ix[:, j-1]
                        n += 1
                    if (j + 1 < len(norm_data.columns)
                        and isfinite(norm_data.ix[0, j+1])):
                        norm_data.ix[:, j] += norm_data.ix[:, j+1]
                        n += 1
                    norm_data.ix[:, j] /= n
        dwg.add(dwg.text(first_col, (x_start,
                                     y_start + (i + 1) * box_height)))
        dwg.add(dwg.text(last_col, (x_start + (new_cols - 1) * box_size,
                                    y_start + (i + 1) * box_height)))
        if draw_box and not draw_average_only:
            dwg.add(dwg.rect((x_start, y_start + 0),
                             (new_cols*box_size, rows*box_height),
                             style="stroke-width:1; "
                             "stroke:#000000; fill:none"))
        if draw_average or draw_average_only:
            avg_frame = norm_data.mean(axis=0)
            for j in range(new_cols):
                col_base = col_labels[j][:col_labels[j].find(col_sep)]
                prefix = col_base
                if cmap_by_prefix:
                    c_cmap = cmap_by_prefix(prefix)
                g = dwg.g()
                g.add(svg.base.Title("Average, {}: {:.2f}".format(col_labels[j],
                                                                  avg_frame.ix[j])))
                g.add(dwg.rect((x_start + box_size*j,
                                y_start + (i+(not draw_average_only))*box_height),
                               (box_size, box_height),
                               style="fill:#{:02x}{:02x}{:02x}"
                               .format(*[int(255*x) for x in
                                         c_cmap(average_scale*avg_frame.ix[j])])))
                if not isfinite(frame.ix[0, j]) and hatch_nan:
                    g.add(dwg.rect((x_start + box_size*j,
                                    y_start + (i+(not draw_average_only))*box_height),
                                   (box_size, box_height),
                                   style="fill:url(#hatch)"
                                  )
                         )

                dwg.add(g)
            dwg.add(dwg.rect((x_start,
                              y_start + (i+(not draw_average_only))*box_height),
                             (new_cols*box_size, 1*box_height),
                             style="stroke-width:1; stroke:#000000; fill:none"
                            ))


        if draw_name:
            if name == "" and split_columns:
                name = col_base
            xpos = x_start + box_size * new_cols / 2.0
            text = dwg.text('',
                             (xpos,
                              y_start
                              + box_height * (rows) * (1-draw_average_only)
                              + box_height * (draw_average or draw_average_only)
                              + 13),
                             style="text-anchor: middle;")
            text.add(dwg.tspan("", dy=["-1.5em"]))
            for line in name.split('_'):
                text.add(dwg.tspan(line,
                                   dy=["1.5em"],
                                   x=[xpos],
                                   style="text-anchor: middle;",
                                   ))
            dwg.add(text)

        if total_width is not None:
            if spacer is None:
                x_start += total_width * 1.1
            else:
                x_start += total_width + spacer
        else:
            if spacer is None:
                x_start += new_cols * box_size + box_size
            else:
                x_start += new_cols * box_size + spacer

        y_diff = new_rows * box_height + 30
        if x_start + total_width >= max_width:
            x_start = 0
            y_start += new_rows*box_height*(not draw_average_only) + vspacer
            y_start += box_height*(draw_average_only or draw_average)

    if draw_row_labels and not draw_average_only:
        for i in range(rows):
            if color_row_labels:
                style = "font-size: {size}; fill: {color};".format(
                    size=box_height,
                    color='red' if row_labels[i] in color_row_labels else 'black',
                )
            else:
                style = "font-size: {}".format(box_height)
            labeltext = (dwg.text(row_labels[i],
                             (x_start, y_start + i*box_height+box_height),
                             style=style,
                            ))
            if make_hyperlinks:
                link = dwg.a('http://insitu.fruitfly.org/cgi-bin/ex/report.pl?ftype=0&ftext={}'
                             .format(row_labels[i]),
                             target='_replace',
                            )
                link.add(labeltext)
                dwg.add(link)
            else:
                dwg.add(labeltext)
    if progress_bar:
        pbar.finish()
    dwg.saveas(filename)
    if convert:
        cmd = [
            'convert',
            filename,
            '-density', '300',
            '-background', 'none',
            filename.replace('svg', 'png'),
        ]
        subprocess.Popen(cmd)



def cmap_by_prefix(prefix):
    cms = dict(
        WT = ISH_CMS_5[0],
        bcd = ISH_CMS_5[1],
        zld = ISH_CMS_5[2],
        G20 = ISH_CMS_5[3],
        hb = ISH_CMS_5[4],
    )
    for p in cms:
        if prefix.startswith(p):
            return cms[p]
    return ISH
