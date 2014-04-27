from __future__ import print_function
import pandas as pd
import numpy as np
import sys
from scipy.stats import linregress, scoreatpercentile
from matplotlib.pyplot import figure, subplot, hist, title, \
        savefig, tight_layout, close, xlim

"""

* Map to mel + vir reference (especially the cufflinks step)

* Do we see the nice dilution series of the vir reads. Calculate a
  histogram of the slopes of (FPKM_i (yeast input) for each gene
  i).  Perhaps most interested to restrict this to FPKMs > 3 (and
  also potentially less than ~1000).

"""

startswith = lambda y: lambda x: x.startswith(y)
contains = lambda y: lambda x: y in x
expr = pd.read_table(sys.argv[1], converters={'gene_short_name':str})
expr.set_index('gene_short_name', inplace=True, verify_integrity=True)

protocols = {c.split('_')[0] for c in expr.columns}
all_samples = {}
print(protocols)

for protocol in protocols:
    try:
        samples = expr.select(startswith(protocol+"_"), axis=1)
        samples = samples.select(contains('subset'), axis=1)
        selector = lambda x: x.startswith('Dvir') and 3 < np.max(samples.ix[x]) < 1000
        x_values = [float(c.split('V')[-1].split('_')[0]) for c in samples.columns]
        samples = samples.select(selector)
        samples = samples.divide(samples.ix[:,2], axis='index')
        all_samples[protocol] = samples

        slopes = pd.Series(index=samples.index)
        intercepts = pd.Series(index=samples.index)
        r_values = pd.Series(index=samples.index)
        for i, gene in enumerate(samples.index):
            res = linregress(x_values, samples.ix[gene])
            slopes.ix[gene] = res[0] * x_values[2]
            intercepts.ix[gene] = res[1]
            r_values.ix[gene] = res[2]

        print('-'*30, '\n', protocol)
        print("analyzed {} genes".format(i))
        figure(figsize=(16,16))
        subplot(3,1,1)
        slopes = slopes.dropna()
        hist(slopes.dropna(),bins=100, range=(-1,5))
        title('Slopes')
        xlim(-1,5)
        print("Slopes")
        print(np.median(slopes), "+/-",)
        print(scoreatpercentile(slopes, 75) - scoreatpercentile(slopes, 25))

        subplot(3,1,2)
        intercepts = intercepts.dropna()
        hist(intercepts,bins=100, range=(-2,2))
        xlim(-2,2)
        title('Intercepts')
        print("Intercepts")
        print(np.median(intercepts.dropna()), "+/-",)
        print(scoreatpercentile(intercepts, 75)
              - scoreatpercentile(intercepts, 25))

        subplot(3,1,3)
        hist(r_values.dropna(),bins=np.arange(.5,1,.01))
        xlim(.5,1)
        title('R Values')
        print("R Values")
        print(np.mean(r_values.dropna()), "+/-",)
        print(np.std(r_values.dropna()))

        tight_layout()
        savefig('analysis/results/{}_virslopes.png'.format(protocol), dpi=150)
    except Exception as error:
        if 'die' in sys.argv:
            raise(error)
        print("Skipping {}, because of a {}".format(protocol, error))

close('all')

