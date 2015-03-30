from __future__ import division, print_function
import BindUtils as bu
import pandas as pd
import DistributionDifference as dd
from scipy.stats import ks_2samp
from Utils import sel_contains, load_to_locals
from collections import defaultdict
from progressbar import ProgressBar as pb
from matplotlib.pyplot import (hist, savefig, legend, clf,
                               title, xlabel)
from numpy import arange
from multiprocessing import Pool
from itertools import repeat

def get_dists_mp(args):
    gene, wts, mut = args
    return dd.earth_mover_multi_rep(
        wts.ix[gene] + eps,
        mut.ix[gene] + eps
    )



if __name__ == "__main__":
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
        (bcd, ('cyc13_rep1', 'cyc13_rep2')),
        (zld, ('cyc13_rep1', 'cyc13_rep2', 'cyc13_rep3')),
        (g20, ('cyc13_rep1', 'cyc13_rep2'))
    ]
    cyc14 = [
        (bcd, ('cyc14D_rep1', 'cyc14D_rep2')),
        (zld, ('cyc14B')),
        (g20, ('cyc14D_rep1'))
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




    if 'keep_old_bm' not in locals():
        keep_old_bm = False
    if 'bm' not in locals() or keep_old_bm == False:
        bm = bu.get_binding_matrix(all_expr.index)
        bm = bm.ix[:, 'bcd cad gt dl gt hb hkb kni kr run shn zld twi tll '.split()]
        keep_old_bm = True
    print("Read binding in")
    bprofs = defaultdict(list)
    for gene in bm.index:
        bprofs[tuple(bm.ix[gene])].append(gene)
    gene_in_cluster_len = pd.Series(index=all_expr.index)
    for cluster in bprofs.itervalues():
        for gene in cluster:
            gene_in_cluster_len.ix[gene] = len(cluster)

    for binds in bprofs:
        cluster = bprofs[binds]
        if len(cluster) > 30:
            score, pv = ks_2samp(all_dists, all_dists.ix[cluster])
            if pv < .05:
                clf()
                binds = bm.columns[[bool(b) for b in binds]]
                print('{:30} {:10.6g} {}'.format(' '.join(binds), pv, score))
                hist(all_dists.ix[cluster], bins=arange(0, 1, 0.05),
                     normed=True, label=' '.join(binds))
                hist(all_dists, bins=arange(0, 1, 0.05),
                     normed=True, histtype='step', label='All')
                legend()
                xlabel('Summed distance')
                title('{} genes'.format(len(cluster)))
                savefig('analysis/results/bindpats/binds-'+'-'.join(binds)+'.png')


