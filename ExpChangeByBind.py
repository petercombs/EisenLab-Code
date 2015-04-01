from __future__ import division, print_function
import BindUtils as bu
import PlotUtils as pu
import pandas as pd
import DistributionDifference as dd
from scipy.stats import ks_2samp
from Utils import sel_contains, load_to_locals, center_of_mass
from collections import defaultdict
from progressbar import ProgressBar as pb
from matplotlib.pyplot import (hist, savefig, legend, clf,
                               title, xlabel, gca)
from numpy import arange
from multiprocessing import Pool
from itertools import repeat
from sys import argv
import setcolor

def get_dists_mp(args):
    gene, wts, mut = args
    return dd.earth_mover_multi_rep(
        wts.ix[gene] + eps,
        mut.ix[gene] + eps
    )



if __name__ == "__main__":
    screen = '--screen' in argv
    expr_min = 5
    eps = 1
    read_table_args = dict(index_col=0,
                           keep_default_na=False,
                           na_values=['---', ''])

    if 'all_expr' not in locals():
        all_expr, (wt, bcd, zld, g20, hb), (wts, bcds, zlds, g20s, hbs), by_cycle\
        = load_to_locals(locals(), expr_min)
    print("Read expression in")
    if 'keep_old' not in locals():
        keep_old = False


    cyc13 = [
        (g20, ('cyc13_rep1', 'cyc13_rep2')),
        (bcd, ('cyc13_rep1', 'cyc13_rep2')),
        (zld, ('cyc13_rep2', 'cyc13_rep3')),
    ]
    cyc14 = [
        (g20, ('cyc14D_rep1',)),
        (bcd, ('cyc14D_rep1', 'cyc14D_rep2')),
        (zld, ('cyc14B', 'cyc14D',)),
    ]

    #mut_st, wt_st = cyc13, 'cyc13'
    mut_st, wt_st = cyc14, 'cyc14D'

    p = Pool()
    if 'all_dists' not in locals() or not keep_old:
        all_dists = pd.Series(index=all_expr.index, data=0.0)
        for mut, cycs in pb()(mut_st):
            cyc_expr = mut.select(**sel_contains(cycs))
            all_dists += p.map(get_dists_mp,
                               zip(all_expr.index,
                                   repeat(wt.select(**sel_contains(wt_st))),
                                   repeat(cyc_expr),
                                  )
                              )


        keep_old=True

    del p


    if 'keep_old_bm' not in locals():
        keep_old_bm = False
    if 'bm' not in locals() or keep_old_bm == False:
        bm = bu.get_binding_matrix(all_expr.index)
        bm = bm.ix[:, bu.ap_early_zld]
        keep_old_bm = True
    print("Read binding in")
    bprofs = defaultdict(set)
    for gene in bm.index:
        bprofs[tuple(bm.ix[gene])].add(gene)
    print('Of {} possible patterns, {} are used {:%}'.format(
        2**len(bm.columns),
        len(bprofs),
        len(bprofs)/2**len(bm.columns)))
    gene_in_cluster_len = pd.Series(index=all_expr.index)
    for cluster in bprofs.itervalues():
        for gene in cluster:
            gene_in_cluster_len.ix[gene] = len(cluster)

    heatmap_data = []
    heatmap_data.append(wt.select(**sel_contains(wt_st)))
    for mut, cycs in mut_st:
        for cyc in cycs:
            cyc = mut.select(**sel_contains(cyc))
            heatmap_data.append(cyc)
    cdt = pd.read_table('analysis/results/all_earth_mover_multi.cdt',
                        index_col='NAME')

    good_clusters = [binds for binds in bprofs if len(bprofs[binds]) > 30]
    for binds in good_clusters:
        cluster = bprofs[binds]
        score, pv = ks_2samp(all_dists, all_dists.ix[cluster])
        if pv * len(good_clusters) < .05:
            clf()
            binds = bm.columns[[bool(b) for b in binds]]
            print(' '.join(binds), len(cluster))
            print('{:30} {:10.6g} {}'.format(' '.join(binds), pv, score))
            hist(all_dists.ix[cluster], bins=arange(0, 1, 0.05),
                 normed=True, label=' '.join(binds))
            hist(all_dists, bins=arange(0, 1, 0.05),
                 normed=True, histtype='step', color='r', label='All')
            legend()
            xlabel('Summed distance')
            title('{} genes'.format(len(cdt.select(cluster.__contains__))))
            if screen:
                setcolor.set_foregroundcolor(gca(), 'w')
                setcolor.set_backgroundcolor(gca(), 'k')
            savefig('analysis/results/bindpats/binds-'+'-'.join(binds)+'.png',
                   transparent=True, dpi=300)
            pu.svg_heatmap(
                data=tuple(heatmap_data),
                filename=('analysis/results/bindpats/binds-'
                          + '-'.join(binds)
                          +'-heatmap.svg'),
                data_names=[c.columns[0].split("_sl")[0] for c in heatmap_data],
                cmap_by_prefix=pu.cmap_by_prefix,
                spacers=10,
                draw_box=True,
                total_width=100,
                index=cdt.select(cluster.__contains__).index,
                progress_bar=True,
                norm_rows_by='max',
                draw_name=True,
                convert=True,
                col_sep='sl',
            )




