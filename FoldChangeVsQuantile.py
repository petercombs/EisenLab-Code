import pandas as pd
from pylab import figure, close, hist, subplot, ylabel, plot, vlines, \
        ylim
from numpy import arange, log2, mean

startswith = lambda y: lambda x: x.startswith(y)

zld_exp = pd.read_table('analysis/summary.tsv', index_col=0)
wt_exp = pd.read_table('prereqs/WT5.53_summary.tsv', index_col=0)

zlds = [zld_exp.select(startswith('cyc13_rep' + str(i)),
                       axis=1).mean(axis=1).sort_index()
        for i in range(1,4)]
wt13 = wt_exp.select(startswith('cyc13'), axis=1).mean(axis=1).sort_index()


close('all')
n = 10
wt_mean = wt13.copy()
wt_mean.sort()
wt_ix = wt_mean[(wt_mean > 3)].index
tot = len(wt_ix)

figure(figsize=(len(zlds), n))
quants = [None, None, None]

for col, zld in enumerate(zlds):

    fc = zld / wt13

    quants_i = []

    for i in range(n):
        subplot(n, len(zlds), i*len(zlds) + col + 1)
        dat = fc.ix[wt_ix[(i *tot//n):((i+1) * tot // n)]]
        quants_i.append(log2(mean(dat)))
        hist(log2(dat),
             arange(-4, 4, .1),
             alpha=.3,
             #normed=True,
             #label="Quantile: {}, n={}, avg FPKM={}".format(i+1,
                                               #sum(isfinite(dat)),
                                              #mean(wt_mean.ix[dat.index]))
            )

        #legend()
        vlines([0], *ylim(), color='r')
        #xlabel(r'$\log_2$ Fold Change zld-/WT')
    quants[col] = quants_i

figure()
plot(quants[0])
plot(quants[1])
plot(quants[2])
