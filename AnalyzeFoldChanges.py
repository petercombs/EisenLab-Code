from __future__ import print_function
import pandas as pd
import numpy as np
import sys
from scipy.stats import linregress, scoreatpercentile
from matplotlib.pyplot import figure, subplot, hist, title, \
        savefig, tight_layout, close, xlim, legend, gca
import setcolor
from os import path

"""

* Map to mel + vir reference (especially the cufflinks step)

* Do we see the nice dilution series of the vir reads. Calculate a
  histogram of the slopes of (FPKM_i (yeast input) for each gene
  i).  Perhaps most interested to restrict this to FPKMs > 3 (and
  also potentially less than ~1000).

"""

startswith = lambda y: lambda x: x.startswith(y)
contains = lambda y: lambda x: y not in x
expfile = sys.argv[1]
if 'in' in expfile:
    outdir = path.join(path.dirname(expfile),
                       'results_{}'.format(expfile[expfile.find('in') +3:
                                                  expfile.find('.')]))
    print("Saving to "+outdir)
else:
    outdir = 'analysis/results'


expr = pd.read_table(expfile, converters={'gene_short_name':str})
expr.set_index('gene_short_name', inplace=True, verify_integrity=True)

protocols = {c.split('_')[0] for c in expr.columns}
all_samples = {}
all_samples_nonorm = {}
all_slopes = {}
all_rs = {}
expr_min = .1
expr_max = 1e3
bin_step = 1
print(protocols)

for protocol in protocols:
    try:
        samples = expr.select(startswith(protocol+"_"), axis=1)
        samples = samples.select(contains('subset'), axis=1)
        selector = lambda x: (x.startswith('Dvir')
                              and expr_min < samples.ix[x, -1]
                              and samples.ix[x, -1] < expr_max)

        x_values = np.array([float(c.split('V')[-1].split('_')[0])
                             for c in samples.columns])
        samples = samples.select(selector)
        normer = np.sum(x_values*samples, axis=1)/(sum(x_values))
        samplesN = samples.divide(normer, axis='index')
        samplesN /= (samplesN.ix[:,-1].mean() / x_values[-1])

        #samples = samplesN
        all_samples_nonorm[protocol] = samples
        all_samples[protocol] = samplesN
        #all_samples[protocol] = samples

        slopes = pd.Series(index=samplesN.index)
        intercepts = pd.Series(index=samplesN.index)
        r_values = pd.Series(index=samplesN.index)
        for i, gene in enumerate(samplesN.index):
            res = linregress(x_values, samplesN.ix[gene])
            slopes.ix[gene] = res[0]
            intercepts.ix[gene] = res[1]
            r_values.ix[gene] = res[2]

        print('-'*30, '\n', protocol)
        print("analyzed {} genes".format(i))
        figure(figsize=(16,16))
        subplot(3,1,1)
        slopes = slopes.dropna()
        slopes_by_expr = [slopes.ix[(10**(i) < samples.ix[:,-1]) *
                                    (samples.ix[:,-1] < 10**(i+bin_step))]
                          .dropna()
                          for i in np.arange(np.floor(np.log10(expr_min)),
                                             np.ceil(np.log10(expr_max)),
                                             bin_step)]
        hist(slopes_by_expr,bins=np.linspace(0, 2, 100), range=(-0,2),
             stacked=True, histtype='bar',normed=True,
             label=['{} < FPKM < {}'.format(10**i, 10**(i+bin_step))
                    for i in np.arange(np.floor(np.log10(expr_min)),
                                       np.ceil(np.log10(expr_max)),
                                       bin_step)])
        title('Slopes')
        legend()
        xlim(-0,2)
        print("Slopes")
        print(np.median(slopes), "+/-",)
        print(scoreatpercentile(slopes, 75) - scoreatpercentile(slopes, 25))

        setcolor.set_foregroundcolor(gca(), 'w')
        setcolor.set_backgroundcolor(gca(), 'k')
        subplot(3,1,2)
        intercepts = intercepts.dropna()
        hist(intercepts,bins=100, range=(-50,50))
        xlim(-50,50)
        title('Intercepts')
        print("Intercepts")
        print(np.median(intercepts.dropna()), "+/-",)
        print(scoreatpercentile(intercepts, 75)
              - scoreatpercentile(intercepts, 25))
        setcolor.set_foregroundcolor(gca(), 'w')
        setcolor.set_backgroundcolor(gca(), 'k')

        subplot(3,1,3)
        hist(r_values.dropna(),bins=np.arange(.5,1,.01))
        xlim(.5,1)
        title('R Values')
        print("R Values")
        print(np.mean(r_values.dropna()), "+/-",)
        print(np.std(r_values.dropna()))

        tight_layout()
        setcolor.set_foregroundcolor(gca(), 'w')
        setcolor.set_backgroundcolor(gca(), 'k')
        savefig('{outdir}/{}_virslopes.png'.format(protocol, outdir=outdir), dpi=150,
                transparent=True)
        all_slopes[protocol] = slopes
        all_rs[protocol] = r_values
    except Exception as error:
        if 'die' in sys.argv:
            raise(error)
        print("Skipping {}, because of a {}".format(protocol, error))

close('all')

