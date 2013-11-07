import pandas as pd
from glob import glob
from collections import defaultdict
from PlotUtils import svg_heatmap
#from matplotlib.pyplot import subplot2grid, pcolormesh, figure, xlim, ylim, \
    #xticks, yticks, vlines, show, cm, title, axis, gca, savefig, pcolor
import matplotlib.pyplot as mpl

from numpy import shape, all

import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = ['Arial']

xlf = pd.ExcelFile('prereqs/ZLD_vs_expression-1203.xlsx')
zld_bind = xlf.parse(xlf.sheet_names[0])
zld_sorted = zld_bind.sort('ZLD Prom st3', ascending=False)
zld_sites = zld_sorted.Gene.dropna().unique()

wt_stages = ('cyc11', 'cyc14A', 'cyc14B')

zld_exp = pd.read_table('analysis/summary.tsv', index_col=0).ix[zld_sites]
#wt_exp = pd.read_table('prereqs/journal.pone.0071820.s008.txt', index_col=0)
wt_exp = pd.read_table('prereqs/WT5.53_summary.tsv', index_col=0).ix[zld_sites]
in_both = set(wt_exp.index).intersection(set(zld_exp.index))

wt_norm = wt_exp.select(lambda x: x.startswith('25u_emb1')
                                  or x.startswith('25u_emb2'),
                         axis=1)
zld_norm = zld_exp.divide(wt_norm.max(axis=1)+.01, axis=0)
wt_norm = wt_norm.divide(wt_norm.max(axis=1)+.01, axis=0)

clust_pairs = defaultdict(set)
clusts = open('prereqs/journal.pone.0071820.s010.kgg')
clusts.readline()
for line in clusts:
    clust_pairs[int(line.split()[1])].add(line.split()[0])

wt_cols = len(wt_norm.columns)
zld_cols = len(zld_norm.columns)

for clustfn in sorted(clust_pairs.keys()):
    print clustfn
    clust = set(clust_pairs[clustfn]).intersection(in_both)
    wt_tmp = (wt_exp.select(lambda x:x in clust)
              .select(lambda x: x.startswith(wt_stages), axis=1)
              #.sort_index()
             )
    zld_tmp = (zld_exp.select(lambda x:x in clust)
               .select(lambda x: x.startswith(wt_stages), axis=1)
               #.sort_index()
              )
    assert all(wt_tmp.index == zld_tmp.index)
    norm_col = wt_tmp.max(axis=1)
    for row in norm_col.index:
        norm_col[row] = max(norm_col[row], 10)
    svg_heatmap((wt_tmp, zld_tmp), 'analysis/results/withzld_{}.svg'.format(clustfn),
                norm_rows_by=norm_col,
                col_sep='_sl',
                total_width=300,
                cmap=(mpl.cm.Blues, mpl.cm.Blues),
                box_size=15,
                draw_row_labels=True)


mpl.draw_if_interactive()
