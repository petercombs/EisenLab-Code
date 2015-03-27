import pandas as pd
#import PlotUtils as pu
from Utils import sel_startswith, sel_contains
from PeakFinder import has_anterior_peak, has_posterior_peak, has_central_peak
from Utils import center_of_mass
from matplotlib.pyplot import (plot, title, savefig, clf, xlim, ylim, xlabel)
from scipy.stats import linregress
from numpy import array, mean
from progressbar import ProgressBar

if __name__ == "__main__":
    read_table_args = dict(index_col=0,
                           keep_default_na=False,
                           na_values=['---', ''])
    if 'all_expr' not in locals():
        all_expr = (pd.read_table('analysis/summary.tsv', **read_table_args)
                    .sort_index())
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

    wt_ants = has_anterior_peak(wts['cyc14D'], fold=2)
    wt_cents = has_central_peak(wts['cyc14D'], fold=2)
    wt_posts = has_posterior_peak(wts['cyc14D'], fold=2)
    xs = array([0, 0, 1, 2.4])
    flops = []
    for pos, genes in [('ants', wt_ants), ('cents', wt_cents)]:
        for gene in ProgressBar()(all_expr.index[genes & ~wt_posts]):
            ys = [
                center_of_mass(bcds['cyc14D_rep1'].ix[gene]),
                center_of_mass(bcds['cyc14D_rep2'].ix[gene]),
                center_of_mass(wts['cyc14D'].ix[gene]),
                center_of_mass(g20s['cyc14D_rep1'].ix[gene]),
            ]

            clf()
            plot(xs, ys, 'b.')
            m, b, r, p, e = linregress(xs, ys)
            if ((mean(ys[:2]) < ys[2] < ys[3])
                or (mean(ys[:2]) > ys[2] > ys[3])):
                plot(xs, m * xs + b, 'r:')
                print gene, mean(ys[:2]), ys[2], ys[3]
            title(gene)
            xlim(-.1, 2.5)
            ylim(0, 0.75)
            xlabel('BCD dosage')
            savefig('analysis/results/coms/{}/{}.png'.format(pos, gene),
                    dpi=300)
            clf()

            if (has_anterior_peak(wts['cyc14D'].ix[gene], fold=2)
                and (has_anterior_peak(bcds['cyc14D_rep1'].ix[gene], fold=-2)
                    and has_anterior_peak(bcds['cyc14D_rep2'].ix[gene], fold=-2))):

                print gene
                flops.append(gene)




