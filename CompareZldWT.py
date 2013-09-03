import pandas as pd
from glob import glob
from collections import defaultdict
from PlotUtils import svg_heatmap
#from matplotlib.pyplot import subplot2grid, pcolormesh, figure, xlim, ylim, \
    #xticks, yticks, vlines, show, cm, title, axis, gca, savefig, pcolor
import matplotlib.pyplot as mpl

from numpy import shape, arange

import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = ['Arial']



zld_exp = pd.read_table('analysis/summary.tsv', index_col=0)
wt_exp = pd.read_table('prereqs/journal.pone.0071820.s008.txt', index_col=0)
in_both = set(wt_exp.index).intersection(set(zld_exp.index))

wt_norm = wt_exp.select(lambda x: x.startswith('25u_emb1') or x.startswith('25u_emb2'),
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
    clust = set(clust_pairs[clustfn]).intersection(in_both)
    wt_tmp = wt_exp.select(lambda x:x in clust).select(lambda x: x in
                                                       wt_norm.columns, 
                                                       axis=1)
    zld_tmp = zld_exp.select(lambda x:x in clust).select(lambda x: x in
                                                         zld_norm.columns, 
                                                         axis=1)
    norm_col = wt_tmp.max(axis=1)
    svg_heatmap((wt_tmp, zld_tmp), 'analysis/results/{}.svg'.format(clustfn),
                norm_rows_by=norm_col,
                col_sep='_sl',
                cmap=(mpl.cm.Blues, mpl.cm.Reds), boxsize=10,
                draw_row_labels=True)


mpl.draw_if_interactive()
