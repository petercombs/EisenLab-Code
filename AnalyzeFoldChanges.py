import pandas as pd
import numpy as np
import sys
from scipy.stats import linregress
from matplotlib.pyplot import figure, subplot, hist, title, savefig, tight_layout

"""

* Map to mel + vir reference (especially the cufflinks step)

* Do we see the nice dilution series of the vir reads. Calculate a
  histogram of the slopes of (FPKM_i (yeast input) for each gene
  i).  Perhaps most interested to restrict this to FPKMs > 3 (and
  also potentially less than ~1000). 

"""

startswith = lambda y: lambda x: x.startswith(y)
expr = pd.read_table(sys.argv[1], converters={'gene_short_name':str})
expr.set_index('gene_short_name', inplace=True, verify_integrity=True)

protocols = {c.split('_')[0] for c in expr.columns}
print protocols

for protocol in protocols:
    samples = expr.select(startswith(protocol), axis=1)
    selector = lambda x: x.startswith('Dvir') and 3 < np.max(samples.ix[x]) < 1000
    x_values = [float(c.split('V')[-1].split('_FPKM')[0]) for c in samples.columns]
    samples = samples.select(selector)

    slopes = pd.Series(index=samples.index)
    intercepts = pd.Series(index=samples.index)
    r_values = pd.Series(index=samples.index)
    for i, gene in enumerate(samples.index):
        res = linregress(x_values, samples.ix[gene])
        slopes.ix[gene] = res[0]
        intercepts.ix[gene] = res[1]
        r_values.ix[gene] = res[2]

    print '-'*30, '\n', protocol
    print "analyzed {} genes".format(i)
    figure(figsize=(16,16))
    subplot(3,1,1)
    hist(slopes,bins=100)
    title('Slopes')
    print "Slopes"
    print np.mean(slopes), "+/-",
    print np.std(slopes)

    subplot(3,1,2)
    hist(intercepts,bins=100)
    title('Intercepts')
    print "Intercepts"
    print np.mean(intercepts), "+/-",
    print np.std(intercepts)

    subplot(3,1,3)
    hist(r_values,bins=100)
    title('R Values')
    print "R Values"
    print np.mean(r_values), "+/-",
    print np.std(r_values)

    tight_layout()
    savefig('analysis/results/{}_virslopes.png'.format(protocol), dpi=150)

