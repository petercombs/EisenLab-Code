#!/usr/bin/env python

""" BayesAssign.py

Use a bayesian approach to assign a set of RNAseq FPKMs to pre-existing
data-sets that are well-determined in either time (Lott et al., 2011) or space
(Fowlkes et al., 2008).

"""

import pandas as pd
import numpy as np
from glob import glob
from scipy import stats
from progressbar import ProgressBar

def prob(sample, reference):
    sample_mean = 0
    sample_var_lo = 1
    sample_var_hi = 1
    n_samples = 0
    for col in sample.index:
        if col.endswith('FPKM'):
            sample_mean += sample[col]
            n_samples += 1
            col_base = col.strip('_FPKM')
            if col_base + '_conf_lo' in sample.index:
                D_lo = sample[col] - sample[col_base + '_conf_lo']
                D_hi = sample[col_base + '_conf_hi'] - sample[col]
                sample_var_lo += D_lo**2
                sample_var_hi += D_hi**2
            else:
                sample_var_lo += sample[col]
                sample_var_hi += sample[col]

    sample_mean /= n_samples
    lo = sample_mean - np.sqrt(sample_var_lo)
    hi = sample_mean + np.sqrt(sample_var_hi)
    lo_prob = stats.zprob((lo - np.mean(reference,axis=1)) /
                          (np.std(reference,axis=1) + 10))
    hi_prob = stats.zprob((hi - np.mean(reference,axis=1)) /
                          (np.std(reference,axis=1) + 10))
    return float(hi_prob - lo_prob)

def bayes(priors, probabilities):
    # P(H|E) = P(E|H) * P(H) / P(E)
    #        = P(E|H) * P(H) / sum(P(E|H_i) * P(H_i))
    probabilities = np.array(probabilities)
    probabilities += .001
    posteriors = [P*H / sum(Pi*Hi for Pi, Hi in zip(probabilities, priors))
                  for P,H in zip(probabilities, priors)]

    # Divide to prevent slow divergence from sum(P_i) == 1
    posteriors = np.array(posteriors)
    posteriors += 1e-10
    return posteriors / sum(posteriors)



frame = pd.read_table('merged_genes.fpkm_tracking')
frame = frame.dropna(how='any')
frame.index = frame['gene_short_name']



cycles = [pd.read_table(f, index_col = 0) for f in glob('../susan/by_cycle/*')]

whole_frame = frame.select(lambda x: x in cycles[0].index)

priors = np.ones(len(cycles)) / len(cycles)

old_priors = np.zeros((len(frame.index), len(priors)))

for set in ['CaS1', 'CaS2', 'CaS3', 'Bcd1', 'Bcd2']:
    frame = whole_frame.select(lambda x: x.startswith(set), axis=1)
    progress = ProgressBar()
    for i, gene in enumerate(progress(frame.index)):
        if gene == 'nan':
            old_priors[i,:] = old_priors[i-1,:]
        all_probs = [prob(frame.ix[gene], cycle.ix[gene]) for cycle in cycles]
        if 0 not in all_probs and np.nan not in all_probs:
            posterior = bayes(priors, all_probs)
            assert not sum(np.isnan(posterior))
            assert all(posterior)
            old_priors[i,:] = priors
            priors = posterior
        else:
            old_priors[i,:] = old_priors[i-1,:]


    best_cycle = cycles[np.argmax(priors)]
    print "Best hit in ", glob('../susan/by_cycle/*')[np.argmax(priors)]


