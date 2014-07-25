import DistributionDifference as DD
import pandas as pd
import numpy as np
from matplotlib.pyplot import (clf, plot, xlabel, ylabel, savefig)

zld = pd.read_table("analysis/summary_merged.tsv", index_col=0,
                                        na_values='-', keep_default_na=False)

wt = pd.read_table("prereqs/WT5.57_summary.tsv", index_col=0,
                                        na_values='-', keep_default_na=False)

contains = lambda y: lambda x: y in x
startswith = lambda y: lambda x: x.startswith(y)

'''
emds1_2 = pd.Series(index=zld.index, data=np.nan)
emds1_3 = pd.Series(index=zld.index, data=np.nan)
emds2_3 = pd.Series(index=zld.index, data=np.nan)

rep1 = zld.select(contains('rep1'), axis=1)
rep2 = zld.select(contains('rep2'), axis=1)
rep3 = zld.select(contains('rep3'), axis=1)


for i in emds1_2.index:
    if max(rep1.ix[i]) > 1 or max(rep2.ix[i]) > 1:
        emds1_2.ix[i] = emd(linspace(0,1,len(rep1.columns), endpoint=True),
                            linspace(0,1,len(rep2.columns), endpoint=True),
                            rep1.ix[i]/sum(rep1.ix[i]+.001)+.001,
                            rep2.ix[i]/sum(rep2.ix[i]+.001)+.001,
                           )

for i in emds2_3.index:
    if max(rep2.ix[i]) > 1 or max(rep3.ix[i]) > 1:
        emds2_3.ix[i] = emd(linspace(0,1,len(rep2.columns), endpoint=True),
                            linspace(0,1,len(rep3.columns), endpoint=True),
                            rep2.ix[i]/sum(rep2.ix[i]+.001)+.001,
                            rep3.ix[i]/sum(rep3.ix[i]+.001)+.001,
                           )
for i in emds1_2.index:
    if max(rep1.ix[i]) > 1 or max(rep3.ix[i]) > 1:
        emds1_3.ix[i] = emd(linspace(0,1,len(rep1.columns), endpoint=True),
                            linspace(0,1,len(rep3.columns), endpoint=True),
                            rep1.ix[i]/sum(rep1.ix[i]+.001)+.001,
                            rep3.ix[i]/sum(rep3.ix[i]+.001)+.001,
                           )
print median(emds1_2.dropna())
print median(emds1_3.dropna())
print median(emds2_3.dropna())

'''

wt_select = wt.select(startswith(tuple('cyc11 cyc13 cyc14A cyc14B'.split())), axis=1)
renamer = lambda x: x.replace('_rep1', '')
zld_select = zld.select(startswith(tuple('cyc11 cyc13_rep1 cyc14A cyc14B'.split())),
                        axis=1).rename(columns=renamer)

allEMDs = pd.Series(data=np.nan, index=zld_select.index)

for i in zld_select.index:
    if max(zld_select.ix[i] > 5) or max(wt_select.ix[i]) > 5:
        allEMDs.ix[i] = DD.earth_mover_multi(zld_select.ix[i],
                                              wt_select.ix[i])

allEMDsg = allEMDs.dropna()
allEMDsg.sort(ascending=False)
zld_bind_all = pd.read_table('prereqs/journal.pgen.1002266.s005.xls', header=1)
zld_bind = set(g.strip() for g in zld_bind_all.TSS_gene)

bin = 100
xs = []
ys = []
print "figuring out plot"
for i in range(0, len(allEMDsg), bin):
    xs.append( mean(allEMDsg.ix[i:i+bin]))
    ys.append(100.0 * sum(i in zld_bind
                          for i in allEMDsg.index[i:i+bin])
              / len(allEMDsg.index[i:i+bin]))


clf()
plot(xs, ys, '.-')
xlabel('sum EMD - WT vs zld-')
ylabel('% with ZLD ChIP-seq site')
savefig('analysis/results/ZldFrac_vs_EMD.png', dpi=300)
