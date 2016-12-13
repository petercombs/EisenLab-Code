from Utils import load_to_locals
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
        expr_min = 15
        all_expr, (wt, bcd, zld, g20, hb), (wts, bcds, zlds, g20s, hbs), by_cycle\
        = load_to_locals(locals(), expr_min)

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
            savefig('analysis/results/coms/{}/{}.eps'.format(pos, gene),
                    dpi=300)
            clf()

            if (has_anterior_peak(wts['cyc14D'].ix[gene], fold=2)
                and (has_anterior_peak(bcds['cyc14D_rep1'].ix[gene], fold=-2)
                    and has_anterior_peak(bcds['cyc14D_rep2'].ix[gene], fold=-2))):

                print gene
                flops.append(gene)




