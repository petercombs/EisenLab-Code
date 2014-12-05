from __future__ import print_function, division
import pandas as pd
import numpy as np
from matplotlib.pyplot import (savefig, tight_layout, subplot, figure, title,
                               xlim, legend, gca, hist)
from numpy.random import multinomial
from scipy.stats import linregress, scoreatpercentile
from os.path import join
from Utils import startswith
import progressbar as pb

amts = np.array([0, 5, 10, 20])
truth_amt = 20

expr_min = 20
expr_max = 1e6
bin_step = 1

truth = pd.read_table(join('analysis',
                           'Truseq_V{}'.format(truth_amt),
                           'subset_min',
                           'subset_genes.fpkm_tracking'),
                      sep='\t',
                      index_col=0,
                      header=None,
                      names=['expression'],
                      #skipfooter=5,
                      keep_default_na=False, na_values=['---'],
                     )

truth = truth.ix[truth.expression > expr_min]
n_reads = sum(truth.expression)

sims = pd.DataFrame(index=truth.index, columns=amts)
is_vir = truth.index.map(startswith('Dvir\\'))

n_sim = 1000

all_slope_vars = np.empty(n_sim)
all_intercept_vars = np.empty(n_sim)
all_corrs_means = np.empty(n_sim)
all_corr_vars = np.empty(n_sim)


pbar = pb.ProgressBar()
for sim_num in pbar(range(n_sim)):
    for amt in amts:
        sim_truth = truth.expression.copy()
        n_mel_old = sum(sim_truth.ix[~is_vir])
        n_mel_new = n_mel_old + (1 - amt/truth_amt) * sum(sim_truth.ix[is_vir])
        sim_truth.ix[is_vir] = sim_truth.ix[is_vir] * amt/truth_amt
        sim_truth.ix[~is_vir] = sim_truth.ix[~is_vir] * n_mel_new / n_mel_old

        assert abs(sum(sim_truth) - n_reads) < 5
        sims.ix[:, amt] = multinomial(n_reads,
                                      sim_truth/n_reads)


    normer = np.sum(amts * sims, axis=1)/sum(amts**2)
    simsN = sims.divide(normer, axis='index')

    slopes = pd.Series(index=truth.index[is_vir])
    intercepts = pd.Series(index=truth.index[is_vir])
    r_values = pd.Series(index=truth.index[is_vir])
    for i, gene in enumerate(truth.index[is_vir]):
        res = linregress(amts, simsN.ix[gene])
        slopes.ix[gene] = res[0]
        intercepts.ix[gene] = res[1]
        r_values.ix[gene] = res[2]
    all_slope_vars[sim_num] = (scoreatpercentile(slopes, 75)
                               - scoreatpercentile(slopes, 25))
    all_intercept_vars[sim_num] = (scoreatpercentile(intercepts, 75)
                                   - scoreatpercentile(intercepts, 25))
    all_corrs_means[sim_num] = np.mean(r_values)
    all_corr_vars[sim_num] = np.std(r_values)


print('-'*30, '\n', 'Simulation')
print("analyzed {} genes".format(i))
figure(figsize=(16,16))
subplot(3,1,1)
slopes = slopes.dropna()
#slopes_by_expr = [slopes.ix[(10**(i) < sims.ix[:,200]) *
                            #(sims.ix[:,200] < 10**(i+bin_step))]
                  #.dropna()
                  #for i in np.arange(np.floor(np.log10(expr_min)),
                                     #np.ceil(np.log10(expr_max)),
                                     #bin_step)]
#hist(slopes_by_expr,bins=np.linspace(0, 2, 100), range=(-0,2),
#stacked=True, histtype='bar',normed=True,
#label=['{} < FPKM < {}'.format(10**i, 10**(i+bin_step))
#for i in np.arange(np.floor(np.log10(expr_min)),
#np.ceil(np.log10(expr_max)),
#bin_step)])
hist(slopes, bins=np.linspace(0, 2, 100), range=(0,2))
title('Slopes')
#legend()
xlim(-0,2)
print("Slopes")
print(np.median(slopes), "+/-",)
print(scoreatpercentile(slopes, 75) - scoreatpercentile(slopes, 25))

        #setcolor.set_foregroundcolor(gca(), 'w')
        #setcolor.set_backgroundcolor(gca(), 'k')
subplot(3,1,2)
intercepts = intercepts.dropna()/(5*max(amts)/100)
hist(intercepts,bins=100, range=(-10,10))
xlim(-10,10)
title('Intercepts (% D. vir)')
print("Intercepts")
print(np.median(intercepts.dropna()), "+/-",)
print(scoreatpercentile(intercepts, 75)
      - scoreatpercentile(intercepts, 25))
#setcolor.set_foregroundcolor(gca(), 'w')
#setcolor.set_backgroundcolor(gca(), 'k')

subplot(3,1,3)
hist(r_values.dropna(),bins=np.arange(.5,1,.01))
xlim(.5,1)
title('R Values')
print("R Values")
print(np.mean(r_values.dropna()), "+/-",)
print(np.std(r_values.dropna()))

tight_layout()

savefig('{outdir}/{}_virslopes.png'.format('simulated',
                                           outdir='analysis/results_subset_min'),
        dpi=150,
        transparent=True)

