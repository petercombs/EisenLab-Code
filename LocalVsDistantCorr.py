from __future__ import division
from Utils import load_to_locals
from scipy.stats import spearmanr
from matplotlib import pyplot as mpl
import setcolor
from sys import argv

if __name__ == "__main__":
    screen = '--screen' in argv
    expr_min = 5
    eps = 1
    sep = 1/3.
    read_table_args = dict(index_col=0,
                           keep_default_na=False,
                           na_values=['---', ''])

    if 'all_expr' not in locals():
        all_expr, (wt, bcd, zld, g20, hb), (wts, bcds, zlds, g20s, hbs), by_cycle\
        = load_to_locals(locals(), expr_min)
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
    if screen:
        setcolor.set_screen(mpl.gcf())
    mpl.savefig('analysis/results/close_vs_far.png', transparent=True)


