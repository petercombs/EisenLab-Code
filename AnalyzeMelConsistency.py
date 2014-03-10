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
from numpy import outer, triu_indices, mean, std, log
from scipy.stats import spearmanr, pearsonr
from scipy.stats.mstats import gmean
from matplotlib import subplot, hist, figure


startswith = lambda y: lambda x: x.startswith(y)
expr = pd.read_table(sys.argv[1], index_col=0, converters={'gene_short_name':str})

protocols = {c.split('_')[0] for c in expr.columns}

for protocol in protocols:
    samples = expr.select(startswith(protocol), axis=1)
    selector = lambda x: 3 < np.mean(samples.ix[x]) < 1000

    fcs = pd.Series(index=samplex.index)
    for i, gene in enumerate(samples.select(selectors).index):
        a = samplex.ix[gene]
        # Geometric mean is proper way to weight all the fold changes
        # Outer product multiplies a by 1/a to give all by all
        # Triu_indices(..,1) only takes elements above the diagonal
        fcs.ix[gene] = gmean(outer(a, 1/a)[triu_indices(len(a),1)])


    print '-'*30, '\n', protocol
    print "analyzed {} genes".format(i)

    print "Mean Log FC", mean(log(fcs))
    print "Std  Log FC", std(log(fcs))

    figure()
    L = len(a)
    for row, samp1 in enumerate(samples.columns):
        for col, samp2 in enumerate(samples.columns):
            if row == col: break
            subplot(L-1, L-1,  (row-1) * (L-1) + col + 1)
            loglog(samples.select(selectors).ix[:, samp2],
                   samples.select(selectors).ix[:, samp1], 
                   'k.', alpha=0.5)
            if col == 0:
                ylabel(row)
            if col + 1 == L:
                xlabel(col)

            print "{} vs {}:\tspearman: {}\t\tpearson: {}".format(
                row, col, spearmanr(samples.select(selectors).ix[:,samp1],
                                    samples.select(selectors).ix[:,samp2]),
                pearsonr(samples.select(selectors).ix[:,samp1],
                                    samples.select(selectors).ix[:,samp2])
            )


    savefig('analysis/results/{}.png'.format(protocol))


