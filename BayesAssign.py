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
import pickle as pkl
import PointClouds as pc

def prob(sample, reference):
    sample_mean = 0
    sample_var_lo = 1
    sample_var_hi = 1
    n_samples = 0
    for col in sample.index:
        if col.endswith('FPKM') or '_' not in col:
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

def bayes(priors, probabilities, prob_boost = .001, post_min = 1e-10):
    # P(H|E) = P(E|H) * P(H) / P(E)
    #        = P(E|H) * P(H) / sum(P(E|H_i) * P(H_i))
    probabilities = np.array(probabilities)
    probabilities += prob_boost/len(priors)
    posteriors = [P*H / sum(Pi*Hi for Pi, Hi in zip(probabilities, priors))
                  for P,H in zip(probabilities, priors)]

    # Divide to prevent slow divergence from sum(P_i) == 1
    posteriors = np.array(posteriors).clip(post_min, 1)
    return posteriors / sum(posteriors)



frame = pd.read_table('merged_genes.fpkm_tracking')
frame = frame.dropna(how='any')
frame.index = frame['gene_short_name']



cycles = [pd.read_table(f, index_col = 0) for f in glob('../susan/by_cycle/*')]

whole_frame = frame.select(lambda x: x in cycles[0].index)

mpl.ion()

for set in ['CaS1', 'CaS2', 'CaS3', 'Bcd1', 'Bcd2', 'Bcd3']:
    print '-'*60, '\n', set, '\n', '-'*60
    priors = np.ones(len(cycles)) / len(cycles)
    old_priors = np.zeros((len(frame.index), len(priors)))

    frame = whole_frame.select(lambda x: x.startswith(set), axis=1)

    progress = ProgressBar()
    for i, gene in enumerate(progress(frame.index)):
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

    pkl_file = open('../Slice60u-NaN-std.pkl')
    bdtnp_parser = pc.PointCloudReader(open('../D_mel_wt_atlas_r2.vpc'))

    starts = pkl.load(pkl_file)
    slices = pkl.load(pkl_file)
    n_pos, n_genes, n_times = np.shape(slices)

    FPKM_cols = [c for c in frame.columns if c.endswith('FPKM') or '_' not in c]
    slice_frames = [pd.DataFrame(slices[:,:,i].T) for i in range(n_times)]
    for slice_frame in slice_frames:
        slice_frame.index = bdtnp_parser.get_gene_names()

    for ts, slice in enumerate(slice_frames):
        priors = np.ones((n_pos, len(FPKM_cols))) / n_pos
        progress = ProgressBar()
        for gene in progress(slice.index):
            if gene not in frame.index: continue
            if sum(np.isnan(slice.ix[gene])): continue
            normed = (slice.ix[gene] / max(slice.ix[gene]) *
                      np.mean(best_cycle.ix[gene], axis=1))
            for i, col in enumerate(FPKM_cols):
                std = (frame[col.replace("FPKM","conf_hi")][gene]
                       - frame[col.replace("FPKM","conf_lo")][gene]) / 2
                if not std:
                    std = frame[col + "_conf_range"][gene]
                evidence = stats.zprob(-np.abs((normed -
                                                frame[col][gene])/(std+1)))


                updated = bayes(priors[:,i], evidence)
                assert not sum(np.isnan(updated))
                priors[:,i] = updated

        print "In time slice", ts
        print np.argmax(priors, axis=0)


