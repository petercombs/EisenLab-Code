from __future__ import print_function, division
import pandas as pd
import matplotlib.pyplot as mpl
from numpy import arange, median, mean
import DistributionDifference as DD
from scipy.stats import scoreatpercentile, gaussian_kde
import setcolor
from Utils import sel_startswith, sel_contains
from multiprocessing import Pool
from sys import argv
from itertools import repeat
from progressbar import ProgressBar as pb


def get_dists_mp(args):
    set1, set2, gene = args
    return DD.earth_mover_multi_rep(set1.ix[gene], set2.ix[gene])

def get_dists(gene):
    return DD.earth_mover_multi_rep(set1.ix[gene], set2.ix[gene])

if __name__ == "__main__":
    screen = '--screen' in argv
    expr_min = 10
    eps = .1
    read_table_args = dict(keep_default_na=False, na_values=['---', ''], index_col=0)
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
            cyc_embs = {}
            for cyc in cycs:
                cyc_embs[cyc] = sub_df.select(**sel_contains(cyc))
            locals()[sub_df_name+'s'] = cyc_embs
    print("Read expression in")

    p = Pool()

    combos14D = [
        ('WT_cyc14D', 'WT_cyc14C'),
        ('WT_cyc14D', 'WT_cyc14E'),
        ('WT_cyc14D', 'bcd_cyc14D'),
        ('WT_cyc14D', 'G20_cyc14D'),
        ('WT_cyc14D', 'zld_cyc14D'),
        ('WT_cyc14D', 'zld_cyc14B'),
        ('WT_cyc14D', 'hb_cyc14D'),
        ('bcd_cyc14D_rep1', 'bcd_cyc14D_rep2'),
        ('hb_cyc14D_rep1', 'hb_cyc14D_rep2'),
    ]

    combos13 = [
        ('WT_cyc13', 'WT_cyc11'),
        ('WT_cyc13', 'WT_cyc14A'),
        ('WT_cyc13', 'bcd_cyc13'),
        ('WT_cyc13', 'G20_cyc13'),
        ('WT_cyc13', 'zld_cyc13'),
        ('WT_cyc13', 'hb_cyc13_rep1'),
        ('bcd_cyc13_rep1', 'bcd_cyc13_rep2'),
        ('G20_cyc13_rep1', 'G20_cyc13_rep2'),
        ('zld_cyc13_rep2', 'zld_cyc13_rep3'),
    ]

    cyc='13'
    for cyc in ['13', '14D']:
        print('='*72)
        print(cyc)
        print('='*72)
        combos = locals()['combos'+cyc]
        mpl.close('all')
        mpl.figure()

        all_dists = {}
        for set1_name, set2_name in combos:
            print('-'*50)
            print(set1_name, set2_name)
            print('-'*50)
            set1 = all_expr.select(**sel_startswith(set1_name))
            set2 = all_expr.select(**sel_startswith(set2_name))

            max_1 = set1.max(axis=1)
            max_2 = set2.max(axis=1)
            keep = (max_1 > expr_min) | (max_2 > expr_min)

            set1 = set1.ix[keep] + eps
            set2 = set2.ix[keep] + eps


            #dists = map(get_dists, set1.index)
            dists = p.map(get_dists_mp, zip(repeat(set1), repeat(set2), set1.index))
            dist_12 = pd.Series(index=set1.index, data=dists)
            all_dists[set1_name + '-' + set2_name] = dist_12

            kernel = gaussian_kde(dist_12)
            if set2_name.startswith('WT'):
                linestyle = ':'
            elif 'rep' in set1_name and 'rep' in set2_name:
                linestyle = '-.'
            else :
                linestyle = '-'
            xvals = arange(0, .2, 1e-4)
            mpl.plot(xvals, kernel(xvals), linestyle,
                     label='{} vs {}'.format(set1_name, set2_name))
            print(mean(dist_12), median(dist_12))
            print(scoreatpercentile(dist_12, [5, 25, 75, 95]))
            print(sum(dist_12 > .15))

        mpl.xlabel('EMD')
        mpl.ylabel('Density')
        mpl.legend()
        if screen:
            setcolor.set_foregroundcolor(mpl.gca(), 'w')
            setcolor.set_backgroundcolor(mpl.gca(), 'b')
        mpl.savefig('analysis/results/pattern_dists_{}.png'.format(cyc), dpi=300)

        mpl.clf()
        dist_list = [all_dists[n1 + '-' + n2]
                     for n1, n2 in combos]
        result = mpl.violinplot(dist_list,
                                showextrema=False,
                                widths=0.9,
                                showmedians=True,
                                points=1000)
        for (n1, n2), body in zip(combos, result['bodies']):
            if 'rep' in n1 and 'rep' in n2:
                body.set_color('b')
            elif 'WT' in n2:
                body.set_color('g')

        mpl.xticks(range(1, len(result['bodies'])+1),
                   ['\n'.join(k) for k in combos],
                   rotation=90)
        mpl.ylim(0, 0.2)
        mpl.tight_layout()
        mpl.savefig('analysis/results/violin{}'.format(cyc), dpi=300)


