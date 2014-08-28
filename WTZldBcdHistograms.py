from __future__ import print_function, division
import pandas as pd
import matplotlib.pyplot as mpl
from numpy import arange
import DistributionDifference as DD
from bisect import bisect

cyc_of_interest = 'cyc14D'

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
    wt_zld = pd.Series(index=wt.index
                       .intersection(zld.index)
                       .intersection(bcd.index),
                      data=0)
    wt_bcd = pd.Series(index=wt_zld.index,
                      data=0)

    print("Calculating Distances")
    zld_reps = {col.split('sl')[0]
                for col in zld.columns
                if col.startswith(cyc_of_interest)
               }
    for zld_rep in zld_reps:
        print('Calculating distance for ', zld_rep)
        for gene in wt_zld.index:
            wt_zld.ix[gene] += (
                DD.earth_mover(wt.ix[gene].select(startswith(cyc_of_interest))+.01,
                               zld.ix[gene].select(startswith(zld_rep))+.01)
            )
    wt_zld /= len(zld_reps)

    bcd_reps = {col.split('sl')[0]
                for col in bcd.columns
                if col.startswith(cyc_of_interest)
               }
    for bcd_rep in bcd_reps:
        print('Calculating distance for ', bcd_rep)
        for gene in wt_bcd.index:
            wt_bcd.ix[gene] += (
                DD.earth_mover(wt.ix[gene].select(startswith(cyc_of_interest))+.01,
                               bcd.ix[gene].select(startswith(bcd_rep))+.01)
            )
    wt_bcd /= len(bcd_reps)

    keep_old = True



both = wt_zld.dropna().index.intersection(wt_bcd.dropna().index)
wt_hi = wt.select(startswith(cyc_of_interest), axis=1).max(axis=1) > 3

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
mpl.scatter(x=wt_zld.ix[both].ix[wt_hi],
            y=wt_bcd.ix[both].ix[wt_hi],
            c=cmaps2.ix[both].ix[wt_hi],
            marker = '.',
            edgecolors='none',
            alpha=0.999999999)
mpl.xlabel('WT vs Zld')
mpl.ylabel('WT vs Bcd')
mpl.xlim(0, 1)
mpl.ylim(0, 1)
mpl.title('{cyc} - {dist:3.1f}kb'
          .format(cyc=cyc_of_interest,
                  dist=bind_dist/1000.))
mpl.savefig('analysis/results/WTBcdWTZldCorr-{}.png'.format(cyc_of_interest), dpi=600)
