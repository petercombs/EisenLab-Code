from scipy.stats import scoreatpercentile
import pandas as pd
from emd import emd
import PlotUtils
from matplotlib.pyplot import clf, xlabel, hist, savefig
from numpy import median, linspace
import numpy as np
from itertools import combinations

expr = pd.read_table('analysis/summary_merged.tsv', keep_default_na=False, na_values='-',
                     index_col=0).dropna(how='any')
emds1_2 = pd.Series(index=expr.index, data=0)
n_emds = pd.Series(index=expr.index, data=0)

contains = lambda y: lambda x: y in x
rep1 = expr.select(contains('rep1'), axis=1)
rep2 = expr.select(contains('rep2'), axis=1)
rep3 = expr.select(contains('rep3'), axis=1)

combs = [rep1, rep2, rep3]
for a,b in combinations([rep1, rep3], 2):
    for i in emds1_2.index:
        if max(a.ix[i]) > 5 and max(b.ix[i])>5:
            n_emds.ix[i] += 1
            emds1_2.ix[i] += emd(linspace(0,1,len(a.columns), endpoint=True),
                                linspace(0,1,len(b.columns), endpoint=True),
                                a.ix[i]/sum(a.ix[i]+.001)+.001,
                                b.ix[i]/sum(b.ix[i]+.001)+.001,
                               )

max_combs = len(combs)*(len(combs)-1)/2
dists = emds1_2.ix[n_emds == max_combs]/max_combs
dists.sort()

print median(dists), scoreatpercentile(dists,95)

clf()
hist(dists, bins=100)
xlabel('Earth Mover Distance')
savefig('analysis/results/all_emds_1vs2.png', dpi=300)

c10 = dists.ix[dists>median(dists)].ix[:10]
pct95 = dists.ix[dists > scoreatpercentile(dists, 95)].ix[:10]
#pct95 = expr.ix[(emds1_2 > scoreatpercentile(dists,95))
                 #* (emds1_2 < scoreatpercentile(dists,95) + .01)].ix[:10]


PlotUtils.svg_heatmap((rep1.ix[c10.index,:], rep2.ix[c10.index,:],
                       rep3.ix[c10.index,:]),
                      'analysis/results/median_emd.svg', norm_rows_by='max', col_sep='rep',
                      draw_box=True,
                      total_width=100, draw_row_labels=False, box_height=25)
PlotUtils.svg_heatmap((rep1.ix[pct95.index,:], rep2.ix[pct95.index,:],
                       rep3.ix[pct95.index,:]),
                      'analysis/results/95th_emd.svg', norm_rows_by='max', col_sep='rep',
                      draw_box=True,
                      total_width=100, draw_row_labels=False, box_height=25)

