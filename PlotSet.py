""" PlotSet.py

quick and dirty scratch function. If it proves useful, should probably be rolled
into PlotUtils.py

"""
from numpy import linspace, mean, std
from matplotlib.pyplot import (clf, subplot, errorbar, plot, legend, xticks, yticks,
                               xlabel, ylabel, ylim, savefig, figure)
import pandas as pd

try:
    len(has_zld)
except NameError:
    zld_bind_all = pd.read_table('prereqs/journal.pgen.1002266.s005.xls',
                                 header=1, keep_default_na=False)
    has_zld = {gene.strip() for gene in zld_bind_all.TSS_gene}

startswith = lambda y: lambda x: x.startswith(y)

def PlotSet(theset, wt, zld,
            filter_zyg=False,
            draw_errorbars=False,
            outfn=None, normstage='cyc14B'):
    """Plot a matched set of the genes provided, broken out by zld binding

    *outfn* The file to save to, otherwise no save

    *normstage*

    """

    w11 = wt.select(startswith('cyc11'), axis=1)
    w13 = wt.select(startswith('cyc13'), axis=1)
    w14A = wt.select(startswith('cyc14A'), axis=1)
    w14B = wt.select(startswith('cyc14B'), axis=1)

    z11 = zld.select(startswith('cyc11'), axis=1)
    z13 = zld.select(startswith('cyc13_rep1'), axis=1)
    z14A = zld.select(startswith('cyc14A'), axis=1)
    z14B = zld.select(startswith('cyc14B'), axis=1)

    myset = set(theset)

    if filter_zyg:
        is_zyg = mean(w14B, axis=1) < mean(w11, axis=1)
        myset.difference_update(is_zyg[is_zyg].index)

    normer = (wt.select(startswith(normstage), axis=1).mean(axis=1))

    max_val = 0
    figure(figsize=(16,12))
    clf()
    for i, (w, z) in enumerate(zip([w11, w13, w14A, w14B],
                                   [z11, z13, z14A, z14B])):
        subplot(1, 4, i+1)
        normw_y = (w.ix[myset.intersection(has_zld)]
                   .divide(normer.ix[myset.intersection(has_zld)].mean(axis=1),
                              axis=0))
        normw_n = (w.ix[myset.difference(has_zld)]
                   .divide(normer.ix[myset.intersection(has_zld)].mean(axis=1),
                              axis=0))
        normz_y = (z.ix[myset.intersection(has_zld)]
                      .divide(normer.ix[myset.intersection(has_zld)].mean(axis=1),
                              axis=0))
        normz_n = (z.ix[myset.difference(has_zld)]
                      .divide(normer.ix[myset.intersection(has_zld)].mean(axis=1),
                              axis=0))
        if draw_errorbars:
            errorbar(linspace(0, 1, len(w.columns), endpoint=True),
                     mean(normw_y), yerr=std(normw_y),
                     fmt='b-*', label='WT +ZLD'
                    )
            errorbar(linspace(0, 1, len(w.columns), endpoint=True),
                     mean(normw_n), yerr=std(normw_n),
                     fmt='b-*', label='WT -ZLD'
                    )
            errorbar(linspace(0, 1, len(z.columns), endpoint=True),
                     mean(normz_y), yerr=std(normz_y),
                     fmt='r-*',
                     label='zld- +ZLD')
            errorbar(linspace(0, 1, len(z.columns), endpoint=True),
                     mean(normz_n), yerr=std(normz_n),
                     fmt='r-o', markerfacecolor='w',
                     label='zld- -ZLD')

        else:
            plot(linspace(0, 1, len(w.columns), endpoint=True),
                 mean(normw_y),
                 'b-*',
                 label='WT +ZLD')
            plot(linspace(0, 1, len(w.columns), endpoint=True),
                 mean(normw_n),
                 'b-o', markerfacecolor='w',
                 label='WT -ZLD')
            plot(linspace(0, 1, len(z.columns), endpoint=True),
                 mean(normz_y),
                 'r-*',
                 label='zld- +ZLD')
            plot(linspace(0, 1, len(z.columns), endpoint=True),
                 mean(normz_n),
                 'r-o', markerfacecolor='w',
                 label='zld- -ZLD')
        xticks([0, 1], 'AP')
        max_val = max(max_val, ylim()[1])
        if i != 0:
            yticks([])
        xlabel(w.columns[0].split('_')[0])

    for i in range(4):
        subplot(1,4,i+1)
        ylim(0, max_val)
    subplot(1,4,1)
    legend(numpoints=1)
    if outfn:
        savefig(outfn)
