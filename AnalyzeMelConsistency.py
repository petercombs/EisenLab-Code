""" AnalyzeMelConsistency
* Map to mel only (discard vir): 

    * Within a given kit, do the different libraries agree with each other?  Mean
      log FC = 0, RMS(log FC) should be small. Spearman  and Pearson coefficients
      high > .95

    * Do they agree with the TruSeq? This is optional, but nice
      to know. Ideally both Spearman and Pearson coefficients are
      high.
"""

import sys
import pandas as pd
import numpy as np
from numpy import outer, triu_indices, mean, std, log, arange
from numpy import log10 as log
from scipy.stats import spearmanr, pearsonr
from scipy.stats.mstats import gmean
from matplotlib.pyplot import subplot, hist, figure, loglog, ylabel, xlabel,\
        savefig, xlim, ylim, xticks, yticks, close, text, gca
from os import path


startswith = lambda y: lambda x: x.startswith(y)
contains = lambda y: lambda x: y not in x
expr = pd.read_table(sys.argv[1], index_col=0, converters={'gene_short_name':str})
outdir = path.join(path.dirname(sys.argv[1]), 'results')

protocols = {c.split('_')[0] for c in expr.columns}
print protocols

for protocol in protocols:
    try:
        samples = expr.select(startswith(protocol+'_'), axis=1)
        samples = samples.select(contains('subset'), axis=1)
        #selector = lambda x: 3 < np.max(samples.ix[x]) < 1000
        samples = samples.ix[(3 < samples.max(axis=1)) 
                             * (samples.min(axis=1) < 1000),
                             :]
        #samples = samples.select(selector)

        fcs = pd.Series(index=samples.index)
        for i, gene in enumerate(samples.index):
            a = samples.ix[gene]
            # Geometric mean is proper way to weight all the fold changes
            # Outer product multiplies a by 1/a to give all by all
            # Triu_indices(..,1) only takes elements above the diagonal
            fcs.ix[gene] = gmean(outer(a, 1/a)[triu_indices(len(a),1)])


        print 'log(10) = ', log(10)
        print '-'*30, '\n', protocol
        print "analyzed {} genes".format(i)

        print "Mean Log FC", mean(log(fcs))
        print "Std  Log FC", std(log(fcs[np.abs(log(fcs)) < 2]))

        figure(figsize=(16,16))
        L = len(a)
        for row, samp1 in enumerate(samples.columns):
            for col, samp2 in enumerate(samples.columns):
                subplot(L, L, row*L + col + 1)
                if row == col:
                    text(.5, .5, samp1, horizontalalignment='center',
                         verticalalignment='center', transform=gca().transAxes,
                         fontsize=20)
                    xticks([])
                    yticks([])
                elif row > col:
                    #subplot(L-1, L-1, (row-1) * (L-1) + col + 1)
                    hist(log(fcs), arange(-1, 1, .01))
                else:
                    #subplot(L-1, L-1,  (row-1) * (L-1) + col + 1)
                    loglog(samples.ix[:, samp2],
                           samples.ix[:, samp1], 
                           'k.', alpha=0.5)
                    xlim(.1,1e5)
                    ylim(.1,1e5)
                    if col == 0:
                        ylabel(samp1)
                        yticks([1e-1, 1e1, 1e3, 1e5])
                    else:
                        yticks([1e-1, 1e1, 1e3, 1e5], ('',)*4)

                    if row + 1 == L:
                        xlabel(samp2)
                        xticks([1e-1, 1e1, 1e3, 1e5])
                    else:
                        xticks([1e-1, 1e1, 1e3, 1e5], ('',)*4)


                    print "{} vs {}:\tspearman: {}\t\tpearson: {}".format(
                        row, col, spearmanr(samples.ix[:,samp1],
                                            samples.ix[:,samp2]),
                        pearsonr(samples.ix[:,samp1],
                                            samples.ix[:,samp2])
                    )


        savefig('{outdir}/{}_logfc.png'.format(protocol, outdir=outdir), dpi=150)
        figure(figsize=(15,15))
        hist(log(fcs), arange(-1,1,.01))
        xlim(-1, 1)
        xlabel('$log_{10} FC$')
        savefig('{outdir}/{}_hist_logfc.png'.format(protocol, outdir=outdir), dpi=150)
        close('all')
    except Exception as error:
        if 'die' in sys.argv:
            raise(error)
        print("Skipping "+protocol)


