import pandas as pd
from glob import glob
from collections import defaultdict
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
    fig = mpl.figure(figsize=((wt_cols+zld_cols)/4.0, len(clust)/4.0))

    #fig.add_subplot(1,1,1, axis_bgcolor=(0,0,0,0), clip_on=False)
    #title(clustfn)
    #xticks()

    #ax = subplot2grid((10,wt_cols + zld_cols), (0,0), colspan=wt_cols, rowspan=10)
    ax = mpl.subplot(1,2,1)
    wt_tmp = wt_norm.select(lambda x:x in clust).sort_index(axis=0)
    wt_shape = shape(wt_tmp)
    print wt_shape,
    mpl.pcolor(wt_tmp.as_matrix(), cmap=mpl.cm.Blues)
    mpl.xlim(0, wt_shape[1])
    mpl.ylim(0, wt_shape[0])
    mpl.xticks([])
    mpl.yticks(arange(0,wt_shape[0]) + 0.5, wt_tmp.index)
    mpl.vlines(14,0,wt_shape[0])
    ax.set_aspect('equal')

    #ax = subplot2grid((10,wt_cols + zld_cols), (0,wt_cols), colspan=zld_cols,
    #             rowspan=10)
    ax = mpl.subplot(1,2,2)
    zld_tmp = zld_norm.select(lambda x:x in clust).sort_index(axis=0)
    zld_shape = shape(zld_tmp)
    print zld_shape
    mpl.pcolor(zld_tmp.as_matrix(), cmap=mpl.cm.Reds)
    mpl.xlim(0, zld_shape[1])
    mpl.ylim(0, zld_shape[0])
    mpl.xticks([])
    mpl.yticks([])
    mpl.vlines(13,0,wt_shape[0])
    ax.set_aspect('equal')

    #subplot2grid((11,wt_cols + zld_cols), (0,0), colspan=(wt_cols + zld_cols))
    #title(clustfn)
    #axis('off')
    mpl.savefig('analysis/results/{}.pdf'.format(clustfn))

mpl.draw_if_interactive()
