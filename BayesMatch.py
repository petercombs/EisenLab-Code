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
from os import path, makedirs
from progressbar import ProgressBar, Bar, ETA, Percentage
import pickle as pkl
import PointClouds as pc
import sys
import argparse

from matplotlib import pyplot as mpl
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = ['Arial']

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
    priors = np.array(priors)
    probabilities = np.array(probabilities)
    probabilities += prob_boost/len(priors)
    denom = sum(probabilities * priors)
    posteriors = probabilities * priors / denom
    # Divide to prevent slow divergence from sum(P_i) == 1
    posteriors = np.array(posteriors).clip(post_min, 1)
    return posteriors / sum(posteriors)

def get_std(column, data):
    if column.replace("FPKM", "conf_hi") in data.columns:
        std = (data[col.replace("FPKM","conf_hi")][gene]
               - data[col.replace("FPKM","conf_lo")][gene]) / 2
    elif column.replace("FPKM", "conf_range") in data.columns:
        std = data[col + "_conf_range"][gene]
    else:
        std = .3 * data[col][gene]
        sys.stderr.write("Warning: Can't find stddev for %s in"
                         "%s" % (gene, col))
    return std

def parse_args():
    description = ('Takes a set of FPKM values from sliced RNAseq data, and'
                   'matches those slices to known gold-standards.')
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('fname', type=open)
    parser.add_argument('--rnaseq-standard-dir', '-r',
                        default='../susan/by_cycle')
    parser.add_argument('--slice-pickle', '-p', type=open,
                        default=open('../Slice60u-NaN-std.pkl'))
    parser.add_argument('--atlas', '-a', type=open,
                        default=open('../D_mel_wt_atlas_r2.vpc'))
    parser.add_argument('--set', '-s', action='append',
                       help='Prefix of columns to use (May include a comma to '
                        'allow multiple prefixes; e.g. --set A,P)',
                       default=[])
    parser.add_argument('--figdir', '-f', default='figures',
                       help='Directory for writing figures into')
    parser.add_argument('--colormap', '-c', default='jet',
                        help='Colormap for plotting')
    parser.add_argument('--figwidth', '-W', default=4, type=float)
    parser.add_argument('--figheight', '-H', default=3, type=float)
    parser.add_argument('--dpi', '-d', default=300, type=int)

    args = parser.parse_args()
    print args.set
    for i, set in enumerate(args.set):
        args.set[i] = tuple(set.split(','))
    if args.set == []:
        args.set.append('')
    print args.set
    return args


args = parse_args()
frame = pd.read_table(args.fname)
frame = frame.dropna(how='any')
frame.index = frame['gene_short_name']


pkl_file = args.slice_pickle
bdtnp_parser = pc.PointCloudReader(args.atlas)

starts = pkl.load(pkl_file)
slices = pkl.load(pkl_file)
slices = slices[0:400,:,:]
n_pos, n_genes, n_times = np.shape(slices)

cycnames = glob(path.join(args.rnaseq_standard_dir, '*'))
cycles = [pd.read_table(f, index_col = 0) for f in cycnames]

whole_frame = frame.select(lambda x: x in cycles[0].index)

mpl.ion()

# Each set of slices should be treated independently.
for set in args.set:
    # Print Header line
    print '-'*60, '\n', set, '\n', '-'*60


    frame = whole_frame.select(lambda x: x.startswith(set), axis=1)

    priors = np.ones(len(cycles)) / len(cycles)
    old_priors = np.zeros((len(frame.index), len(priors)))
    FPKM_cols = [c for c in frame.columns if c.endswith('FPKM') ]
    FPKM_cols = sorted(FPKM_cols, key= lambda x: x.count('P') - x.count('A'))

    print "Found columns: ", FPKM_cols

    # ==============================================================
    # Match to the correct time-slice
    # ==============================================================
    widgets = ['Susan: ' + str(set) + ':', Percentage(), Bar(), ETA()]
    progress = ProgressBar(widgets=widgets)

    mpl.figure()
    for i, gene in enumerate(progress(frame.index)):
        all_probs = np.array([prob(frame.ix[gene], cycle.ix[gene])
                              for cycle in cycles])
        if 0 not in all_probs and np.nan not in all_probs:
            posterior = bayes(priors, all_probs, post_min=0)
            #assert not sum(np.isnan(posterior))
            #assert all(posterior)
            old_priors[i,:] = priors
            priors = posterior
        else:
            old_priors[i,:] = old_priors[i-1,:]
    mpl.plot(old_priors)

    best_cycle = cycles[np.argmax(priors)]
    print "Best hit in ", cycnames[np.argmax(priors)]
    sys.stdout.flush()

    #==============================================================
    # Match each slice to the correct positions
    #==============================================================

    slice_frames = [pd.DataFrame(slices[:,:,i].T) for i in range(n_times)]
    for slice_frame in slice_frames:
        slice_frame.index = bdtnp_parser.get_gene_names()

    for ts, slice in enumerate(slice_frames):
        mpl.figure(figsize=(args.figwidth, args.figheight))
        slice = slice.dropna(how='any')
        priors = np.ones((n_pos, len(FPKM_cols))) / n_pos
        widgets = ['Time %s:'%ts, Percentage(), Bar(), ETA()]
        progress = ProgressBar(widgets=widgets)
        for gene in progress(slice.index):
            if gene not in frame.index: continue
            #if sum(np.isnan(slice.ix[gene])):
            #    assert False
            #    continue
            normed = (slice.ix[gene] / max(slice.ix[gene]) *
                      np.mean(best_cycle.ix[gene], axis=1))
            for i, col in enumerate(FPKM_cols):
                std = get_std(col, frame)

                evidence = stats.zprob(-np.abs((normed -
                                                frame[col][gene])/(std+1)))

                updated = bayes(priors[:,i], evidence)
                #assert not sum(np.isnan(updated))
                priors[:,i] = updated
        my_cm = mpl.cm.__getattribute__(args.colormap)
        n_pos, n_samples = np.shape(priors)
        plots = []
        for i in range(n_samples):
            plots.extend(mpl.plot(priors[:,i],
                                  label=FPKM_cols[i].replace('_FPKM', ''),
                                  color = my_cm(i * 256 / (n_samples-1))))
        ax = mpl.gca()
        Y = priors.max()
        dY = 0.25 * Y
        Y += dY
        ests = np.argmax(priors, axis=0)
        for plot, x in zip(plots, ests):
            mpl.vlines(x, plot.get_data()[1][x], Y+dY,
                       edgecolor=plot.get_color(), linestyles='dotted')
            ax.add_artist(mpl.Rectangle((x, Y), 60, dY,
                                        facecolor=plot.get_color(), alpha=0.7))

        ax.set_ylim(0, Y+dY)
        ax.set_xlim(0, n_pos + 60)
        mpl.title("Slice Position estimates for %s" % str(set))
        mpl.xlabel("A/P position ($\\mu$m)")
        mpl.ylabel("P(start @ $x\\pm1\\mu$m)")
        mpl.tight_layout()
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
        mpl.draw()
        try:
            makedirs(args.figdir)
        except OSError:
            pass
        filename = "Sample_%s_T%s.pdf" % (str(set).translate(None,
                                                             ' \'(),"'),
                                          ts)
        print "Saving to: ", path.abspath(path.join(args.figdir, filename))
        mpl.savefig(path.join(args.figdir,filename),
                   dpi=args.dpi)


        print "In time slice", ts
        print "Mode position", np.argmax(priors, axis=0)
        print "Mean position", [sum(np.arange(n_pos) * priors[:,i])
                                for i in range(n_samples)]
        sys.stdout.flush()


