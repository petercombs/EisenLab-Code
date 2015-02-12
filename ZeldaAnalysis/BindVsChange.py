from __future__ import division
import pandas as pd
import DistributionDifference
import progressbar as pb
import matplotlib.pyplot as mpl
import numpy as np
import setcolor

zldbind = pd.read_table('prereqs/journal.pgen.1002266.s005.xls', skiprows=1)
has_zld = set(item.strip() for item in zldbind.TSS_gene)

for i, row in zldbind.iterrows():
    zldbind.TSS_gene[i] = row['TSS_gene'].strip()

assert not zldbind.TSS_gene[0].startswith(' ')

exp_min = 5

wt_exp = pd.read_table('prereqs/WT5.53_summary.tsv', index_col=0)
zld_exp = pd.read_table('analysis/summary_merged.tsv', index_col=0)
wt_exp.sort_index(inplace=True)
zld_exp.sort_index(inplace=True)
assert np.all(wt_exp.index == zld_exp.index)

wt_14 =  wt_exp .select(lambda x: x.startswith('cyc14A'), axis=1)+.1
zld_14 = zld_exp.select(lambda x: x.startswith('cyc14A'), axis=1)+.1

is_expr = (wt_14.max(axis=1) > 5) | (zld_14.max(axis=1) > 5)

pbar = pb.ProgressBar()
diffs = pd.Series(index=wt_exp.ix[is_expr].index)
for gene in pbar(diffs.index):
    metric = DistributionDifference.earth_mover
    diffs[gene] = metric(wt_14.ix[gene], zld_14.ix[gene])

diffs.sort()

binsize = 100
xs = []
ys = []
xs2 = []
ys2 = []
pbar = pb.ProgressBar()
for i in pbar(np.arange(0, len(diffs), binsize)):
    dat = diffs[i:i + binsize]
    n = len(dat)
    mean_diffval = dat.mean()
    p_bind = 0
    ix = set(dat.index)
    for gene in ix:
        if gene in has_zld:
            p_bind += 1

    binds = zldbind.select(lambda x: zldbind.ix[x]['TSS_gene'] in ix)
    quartile = '25%'
    upperQ = binds.Height.describe()[quartile]

    if not np.isnan(upperQ):
        xs2.append(mean_diffval)
        ys2.append(upperQ)

    p_bind /= n
    xs.append(mean_diffval)
    ys.append(p_bind)


mpl.figure(figsize=(6,4.5))
mpl.plot(xs, ys, 'r.')
mpl.xlabel('Difference metric: ' + metric.__name__)
mpl.ylabel('Frac Bound')
setcolor.set_backgroundcolor(mpl.gca(), 'k')
setcolor.set_foregroundcolor(mpl.gca(), 'w')
mpl.savefig('analysis/results/PbindVsEMD.png', background='None',
            dpi=600,
            transparent=True)

mpl.figure()
mpl.semilogx(xs2, ys2)
mpl.xlabel('Difference metric: ' + metric.__name__)
mpl.ylabel(quartile + 'quartile')
mpl.savefig('analysis/results/QuartileVsEMD.png')
