from scipy.stats import scoreatpercentile
import pandas as pd
from emd import emd
import PlotUtils
from matplotlib.pyplot import clf, xlabel
from numpy import median

expr = pd.read_table('../summary.tsv', keep_default_na=False, na_values='-',
                     index_col=0)
emds1_2 = pd.Series(index=expr.index, data=np.nan)
contains = lambda y: lambda x: y in x
rep1 = expr.select(contains('rep1'), axis=1)
rep2 = expr.select(contains('rep2'), axis=1)
rep3 = expr.select(contains('rep3'), axis=1)
for i in emds1_2.index:
    if max(rep1.ix[i]) > 1 or max(rep2.ix[i])>1:
        emds1_2.ix[i] = emd(linspace(0,1,6, endpoint=True),
                            linspace(0,1,6, endpoint=True),
                            rep1.ix[i]/sum(rep1.ix[i]+.001)+.001, 
                            rep2.ix[i]/sum(rep2.ix[i]+.001)+.001, 
                           )

dists = emds1_2.dropna()
dists.sort()

print median(dists), scoreatpercentile(dists,95)

clf()
hist(emds1_2.dropna(), bins=100)
xlabel('Earth Mover Distance')
savefig('all_emds_1vs2.png', dpi=300)

c10 = dists.ix[dists>median(dists)].ix[:10]
pct95 = dists.ix[dists > scoreatpercentile(dists, 95)].ix[:10]
#pct95 = expr.ix[(emds1_2 > scoreatpercentile(dists,95)) 
                 #* (emds1_2 < scoreatpercentile(dists,95) + .01)].ix[:10]


PlotUtils.svg_heatmap((rep1.ix[c10.index,:], rep2.ix[c10.index,:],
                       rep3.ix[c10.index,:]),
                      'median_emd.svg', norm_rows_by='max', col_sep='rep',
                      total_width=100, draw_row_labels=True, box_height=25)
PlotUtils.svg_heatmap((rep1.ix[pct95.index,:], rep2.ix[pct95.index,:],
                       rep3.ix[pct95.index,:]),
                      '95th_emd.svg', norm_rows_by='max', col_sep='rep',
                      total_width=100, draw_row_labels=True, box_height=25)

