from __future__ import division
import pandas as pd
from numpy import log2
from Utils import sel_startswith, sel_contains
from multiprocessing import Pool
import DistributionDifference as dd
from matplotlib.pyplot import (violinplot, ylabel, savefig, xticks,
                               close)

def get_dists_mp(args):
    gene, row = args
    temp = pd.Series(data=1, index=row.index)
    return dd.earth_mover_multi_rep(row+eps, temp)

if __name__ == "__main__":
    expr_min = 15
    eps = 1
    max_diff = .04
    read_table_args = dict(index_col=0,
                           keep_default_na=False,
                           na_values=['---', ''])

    if 'all_expr' not in locals():
        all_expr = (pd.read_table('analysis/summary.tsv', **read_table_args)
                    .sort_index())
        top_expr = all_expr.max(axis=1)
        all_expr = all_expr.ix[top_expr > expr_min]
        wt  = all_expr.select(**sel_startswith('WT'))
        bcd = all_expr.select(**sel_startswith('bcd'))
        zld = all_expr.select(**sel_startswith('zld'))
        g20 = all_expr.select(**sel_startswith('G20'))
        hb  = all_expr.select(**sel_startswith('hb'))

        wts = bcds = zlds = g20s = hbs = 0
        by_cycle = {}
        for sub_df_name in 'wt bcd zld g20 hb'.split():
            sub_df = locals()[sub_df_name]
            cycs = {col.split('_sl')[0].split('_',1)[1] for col in sub_df.columns}
            cycs.update({col.split('_')[1] for col in sub_df.columns})
            cyc_embs = {}
            by_cycle[sub_df_name] = cyc_embs
            for cyc in cycs:
                cyc_embs[cyc] = sub_df.select(**sel_contains(cyc))
            locals()[sub_df_name+'s'] = cyc_embs
    print("Read expression in")


    p = Pool()
    diff_from_ubiq_wt = pd.Series(data=p.map(get_dists_mp,
                                             wts['cyc14D'].iterrows()),
                                  index=all_expr.index)
    diff_from_ubiq_bcd = pd.Series(data=p.map(get_dists_mp,
                                              bcds['cyc14D'].iterrows()),
                                  index=all_expr.index)
    diff_from_ubiq_g20 = pd.Series(data=p.map(get_dists_mp,
                                              g20s['cyc14D'].iterrows()),
                                  index=all_expr.index)

    ubiq_all = ((diff_from_ubiq_wt < max_diff)
                & (diff_from_ubiq_bcd < max_diff)
                & (diff_from_ubiq_g20 < max_diff))

    up_both = []
    down_both = []
    thresh = 1.5
    wt_bcd_avg = wts['cyc14D'].mean(axis=1)/(eps + bcds['cyc14D'].mean(axis=1))
    wt_g20_avg = wts['cyc14D'].mean(axis=1)/(eps + g20s['cyc14D'].mean(axis=1))
    bcd_g20_avg = ((eps+bcds['cyc14D'].mean(axis=1))
                   /(eps+g20s['cyc14D'].mean(axis=1)))
    for gene in ubiq_all.ix[ubiq_all].index:
        if (wt_bcd_avg.ix[gene] < 1/thresh) and (wt_g20_avg.ix[gene] < 1/thresh):
            up_both.append(gene)
        elif ((wt_bcd_avg.ix[gene] > thresh)
              and (wt_g20_avg.ix[gene] > thresh)):
            down_both.append(gene)


    print("Up in both: {:>12} \nDown in both: {:>10}"
          .format(len(up_both), len(down_both)))

    close('all')
    violinplot([log2(bcd_g20_avg.ix[up_both]),
                log2(bcd_g20_avg.ix[down_both])],
               showmedians=True, showextrema=False)
    xticks([1, 2], ['Up in both', 'Down in both'])
    ylabel('$\\log_2 bcd^-/bcd^{+++}$')
    savefig('analysis/results/updownboth_log2')




