import pandas as pd
from glob import glob
from collections import defaultdict

zld_exp = pd.read_table('zeldaSummary.tsv', index_col=0)
wt_exp = pd.read_table('wtSummary.tsv', index_col=0)
in_both = set(wt_exp.index).intersection(set(zld_exp.index))

wt_norm = wt_exp.select(lambda x: x.startswith('emb1') or x.startswith('emb2'),
                         axis=1)
wt_norm = wt_norm.divide(wt_norm.max(axis=1)+.01, axis=0)
zld_norm = zld_exp.divide(wt_norm.max(axis=1)+.01, axis=0)

clust_pairs = defaultdict(set)
clusts = open('../May2013PaperPush/60um_filtered_norm_K_G20.kgg')
clusts.readline()
for line in clusts:
    clust_pairs[line.split()[1]].add(line.split()[0])

wt_cols = len(wt_norm.columns)
zld_cols = len(zld_norm.columns)

for clustfn in clust_pairs:
    fig = figure()
    clust = set(clust_pairs[clustfn]).intersection(in_both)

    #fig.add_subplot(1,1,1, axis_bgcolor=(0,0,0,0), clip_on=False)
    #title(clustfn)
    #xticks()

    subplot2grid((1,wt_cols + zld_cols), (0,0), colspan=wt_cols)
    wt_tmp = wt_norm.select(lambda x:x in clust).sort_index(axis=0)
    wt_shape = shape(wt_tmp)
    print wt_shape,
    pcolormesh(wt_tmp.as_matrix(), cmap=cm.Blues)
    xlim(0, wt_shape[1])
    ylim(0, wt_shape[0])
    xticks([])
    yticks(arange(0,wt_shape[0]) + 0.5, wt_tmp.index)
    vlines(14,0,wt_shape[0])

    subplot2grid((1,wt_cols + zld_cols), (0,wt_cols), colspan=zld_cols)
    zld_tmp = zld_norm.select(lambda x:x in clust).sort_index(axis=0)
    zld_shape = shape(zld_tmp)
    print zld_shape
    pcolormesh(zld_tmp.as_matrix(), cmap=cm.Reds)
    xlim(0, zld_shape[1])
    ylim(0, zld_shape[0])
    xticks([])
    yticks([])
    vlines(13,0,wt_shape[0])
    savefig('analysis/results/clust_{}.png'.format(clust_fn))

show()
