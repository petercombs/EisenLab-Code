from __future__ import print_function
import pandas as pd
import numpy as np
import sys
from scipy.stats import linregress, scoreatpercentile
from matplotlib.pyplot import figure, subplot, hist, title, \
        savefig, tight_layout, close, xlim, legend, gca
import setcolor
from os import path, makedirs

# Python 2/3 compatibility code
try:
    import builtins
except ImportError:
    import __builtin__ as builtins

if "FileExistsError" not in dir(builtins):
    FileExistsError = OSError

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
                                                  expfile.rfind('.')]))
    print("Saving to "+outdir)
    try:
        makedirs(outdir)
    except (FileExistsError, OSError):
        pass
else:
    outdir = path.join(path.dirname(expfile),
                       'results')

protocol_map = dict(
    Truseq = 'TruSeq',
    Clontech = 'Clontech',
    EpicentreOligodT = 'TotalScript',
    SMART2 = 'Smart-seq2',
    SMART2dil2 = 'Smart-seq2, 2.5 fold dilution',
    SMART2dil5 = 'Smart-seq2, 5 fold dilution',
)


outtex = open(path.join(outdir, 'foldchanges.tex'), 'w')
outtex.write(r'''
\begin{table}[htdp]

\caption{Distribution of fit parameters. A simple linear fit,
             $\hat{A}_{ij} = m \cdot Q_{j} + b$
             was computed for each gene $i$, and a correlation coefficent $r$
             calculated.  For brevity,
             $\bar{x}$ is the mean of some variable $x$, and $\sigma_x$ is its
             standard deviation.  }
             \begin{center}
             \begin{tabular}{|c|r@{$\pm$}l|r@{$\pm$}l|r@{$\pm$}l|}
             \hline Protocol & $\bar{m}$ & $\sigma_m$ & $\bar{b}$ & $\sigma_b$
                    & $\bar{r}$ & $\sigma_r$ \\\hline
             ''')

#expr = pd.read_table(expfile, converters={'gene_short_name':str})
#expr.set_index('gene_short_name', inplace=True, verify_integrity=True)
expr = pd.read_table(expfile, converters={"0":str})
expr.rename(columns={'0':'gene_short_name'}, inplace=True)
expr.set_index('gene_short_name', inplace=True, verify_integrity=True)

protocols = {c.split('_')[0] for c in expr.columns}
all_samples = {}
all_samples_nonorm = {}
all_slopes = {}
all_intercepts = {}
all_rs = {}
expr_min = 20
expr_max = 1e6
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
        if 0 in samples.shape:
            print("Skipping", protocol, "due to low reads")
            continue
        normer = np.sum(x_values*samples, axis=1)/(sum(x_values**2))
        samplesN = samples.divide(normer, axis='index')
        #samplesN /= (samplesN.ix[:,-1].mean() / x_values[-1])

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
        intercepts = intercepts.dropna()/(5*max(x_values)/100)
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

        outtex.write('{prot} & {slope:0.3} & {stdslope:0.3} & {b:0.3} & '
                     '{stdb:0.3} & {r:0.2} & {stdr:0.2} \\\\\n'
                     .format(prot=protocol_map[protocol],
                             slope=np.mean(slopes),
                             stdslope=np.std(slopes),
                             b=np.mean(intercepts),
                             stdb=np.std(intercepts),
                             r=np.mean(r_values.dropna()),
                             stdr=np.std(r_values.dropna())
                            )
                    )

        tight_layout()
        #setcolor.set_foregroundcolor(gca(), 'w')
        #setcolor.set_backgroundcolor(gca(), 'k')
        savefig('{outdir}/{}_virslopes.png'.format(protocol, outdir=outdir), dpi=150,
                transparent=True)
        all_slopes[protocol] = slopes
        all_intercepts[protocol] = intercepts
        all_rs[protocol] = r_values
    except Exception as error:
        if 'die' in sys.argv:
            raise(error)
        print("Skipping {}, because of a {}".format(protocol, error))

close('all')

outtex.write('''
\hline

\end{tabular}
\label{tab:fits}
\end{center}
\end{table}



              ''')
outtex.close()
