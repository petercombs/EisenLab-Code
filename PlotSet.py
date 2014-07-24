""" PlotSet.py

quick and dirty scratch function. If it proves useful, should probably be rolled
into PlotUtils.py

"""
from numpy import linspace, mean
from matplotlib.pyplot import (clf, subplot, plot, legend, xticks, yticks,
                               xlabel, ylim, savefig, figure)
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
        plot(linspace(0, 1, len(w.columns), endpoint=True),
             mean(w.ix[myset.intersection(has_zld)]
                  .divide(normer.ix[myset.intersection(has_zld)].mean(axis=1),
                          axis=0)),
             'b-*',
             label='WT +ZLD')
        plot(linspace(0, 1, len(w.columns), endpoint=True),
             mean(w.ix[myset.difference(has_zld)]
                  .divide(normer.ix[myset.intersection(has_zld)].mean(axis=1),
                          axis=0)),
             'b-o', markerfacecolor='w',
             label='WT -ZLD')
        plot(linspace(0, 1, len(z.columns), endpoint=True),
             mean(z.ix[myset.intersection(has_zld)]
                  .divide(normer.ix[myset.intersection(has_zld)].mean(axis=1),
                          axis=0)),
             'r-*',
             label='zld- +ZLD')
        plot(linspace(0, 1, len(z.columns), endpoint=True),
             mean(z.ix[myset.difference(has_zld)]
                  .divide(normer.ix[myset.intersection(has_zld)].mean(axis=1),
                          axis=0)),
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
