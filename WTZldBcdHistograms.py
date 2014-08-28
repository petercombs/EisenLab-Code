from __future__ import print_function, division
import pandas as pd
import matplotlib.pyplot as mpl
from numpy import arange
import DistributionDifference as DD
from bisect import bisect


read_table_args = dict(keep_default_na=False, na_values='---', index_col=0)

bind_dist = 1000

try:
    wt = locals()['wt']
except KeyError:
    wt = pd.read_table('prereqs/WT6.01.summary.tsv', **read_table_args)
    zld = pd.read_table('prereqs/Zld6.01.summary_merged.tsv', **read_table_args)
    bcd = pd.read_table('analysis/summary.tsv', **read_table_args)
    bcd = bcd.dropna(how='all', axis=1)

try:
    tss_dict = locals()['tss_dict']
except KeyError:
    from collections import defaultdict
    tss_dict = defaultdict(list)
    for i, row in pd.read_table('Reference/tss').iterrows():
        tss_dict[row['chr'].replace('dmel_', '')].append((row['TSS_start'], row['gene_name']))

    for chrom in tss_dict.itervalues():
        chrom.sort()

def find_near(chrom, coord, dist):
    out = set()
    center = bisect([i[0] for i in chrom], coord)
    for (coord2, gene) in reversed(chrom[:center]):
        if coord - coord2 > dist: break
        out.add(gene)

    for (coord2, gene) in chrom[center:]:
        if coord2 - coord > dist: break
        out.add(gene)

    return out




zld_bind = pd.read_table('Reference/zld_peaks')
has_zld = set()
for coord in zld_bind.NewPeak:
    chr, coord = coord.split(':')
    has_zld.update(find_near(tss_dict[chr], int(coord), bind_dist))

bcd_bind = pd.read_table('Reference/bcd_peaks')
has_bcd = set()
for coord in bcd_bind.NewPeak:
    chr, coord = coord.split(':')
    has_bcd.update(find_near(tss_dict[chr], int(coord), bind_dist))



startswith = lambda x: lambda y: y.startswith(x)

try:
    assert locals()['keep_old']
    wt_zld = locals()['wt_zld']
    wt_bcd = locals()['wt_bcd']
except (KeyError, AssertionError):
    wt_zld = pd.Series(index=wt.index)
    wt_bcd = pd.Series(index=wt.index)



    print("Calculating Distances")
    for gene in wt.index.intersection(zld.index):
        wt_zld.ix[gene] = (DD.earth_mover(wt.ix[gene].select(startswith('cyc14D'))+.01,
                                          zld.ix[gene].select(startswith('cyc14D'))+.01)
                          )

    for gene in wt.index.intersection(bcd.index):
        wt_bcd.ix[gene] = (DD.earth_mover(wt.ix[gene].select(startswith('cyc14D'))+.01,
                                          bcd.ix[gene].select(startswith('cyc14D_rep1'))+.01)
                           +DD.earth_mover(wt.ix[gene].select(startswith('cyc14D'))+.01,
                                           bcd.ix[gene].select(startswith('cyc14D_rep2'))+.01))/2



both = wt_zld.dropna().index.intersection(wt_bcd.dropna().index)
wt14D_hi = wt.select(startswith('cyc14D'), axis=1).max(axis=1) > 3

print("Plotting")

mpl.clf()
mpl.hist(wt_zld.dropna(), bins=arange(0, 1, .05), normed=True, label='WT vs Zld', histtype='step', color='r')
mpl.hist(wt_bcd.dropna(), bins=arange(0, 1, .05), normed=True, label='WT vs Bcd', histtype='step', color='g')
mpl.legend()
mpl.savefig('analysis/results/WTZldBcdHist.png')

colors = pd.Series('y g b c'.split())
cmaps = pd.Series(index=wt.index, data=0, dtype=int)
for gene in cmaps.index:
    for i, has_bind in enumerate([has_zld, has_bcd]):
            cmaps.ix[gene] += 2**i * (gene in has_bind)

cmaps2 = colors[cmaps]
cmaps2.index = cmaps.index

mpl.clf()
mpl.plot([0,1], [0,1], 'r:')
mpl.scatter(x=wt_zld.ix[both].ix[wt14D_hi],
            y=wt_bcd.ix[both].ix[wt14D_hi],
            c=cmaps2.ix[both].ix[wt14D_hi],
            marker = '.',
            edgecolors='none',
            alpha=0.2)
mpl.xlabel('WT vs Zld'); mpl.ylabel('WT vs Bcd')
mpl.savefig('analysis/results/WTBcdWTZldCorr.png', dpi=600)
