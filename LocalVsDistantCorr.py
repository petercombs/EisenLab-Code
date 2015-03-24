from __future__ import division
import pandas as pd
from Utils import sel_startswith, sel_contains
from scipy.stats import spearmanr
from matplotlib import pyplot as mpl

if __name__ == "__main__":
    expr_min = 5
    eps = 1
    sep = 1/3.
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
        for sub_df_name in 'wt bcd zld g20 hb'.split():
            sub_df = locals()[sub_df_name]
            cycs = {col.split('_sl')[0].split('_',1)[1] for col in sub_df.columns}
            cycs.update({col.split('_')[1] for col in sub_df.columns})
            cyc_embs = {}
            for cyc in cycs:
                cyc_embs[cyc] = sub_df.select(**sel_contains(cyc))
            locals()[sub_df_name+'s'] = cyc_embs
    print("Read expression in")

    rhos_close = []
    rhos_far = []
    for cycs in [wts, bcds, zlds, g20s, hbs]:
        for cyc in cycs:
            cyc = cycs[cyc]
            n_cols = len(cyc.columns)
            for col1 in range(n_cols):
                for col2 in range(n_cols):
                    if col2 > col1:
                        break
                    if col1 == col2:
                        continue
                    if ((len(cyc.ix[:, col1].dropna()) == 0)
                        or (len(cyc.ix[:, col2].dropna()) == 0)):
                        continue
                    if col1 - col2 == 1:
                        rhos_close.append(
                            spearmanr(cyc.ix[:, col1], cyc.ix[:, col2])[0]
                        )
                    elif col1 - col2 > (n_cols*sep):
                        rhos_far.append(
                            spearmanr(cyc.ix[:, col1], cyc.ix[:, col2])[0]
                        )

    mpl.close('all')
    mpl.violinplot([rhos_close, rhos_far], showmedians=True)
    mpl.ylabel('Spearman Rank Correlation')
    mpl.xticks([1, 2],
               ['Adjacent\nN={}'.format(len(rhos_close)),
                'Distant\n>{:0.1%} embryo length\nN={}'.format(sep, len(rhos_far))],
               rotation=0,
              )
    mpl.tight_layout()
    mpl.savefig('analysis/results/close_vs_far.png')


