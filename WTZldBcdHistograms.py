from __future__ import print_function, division
import pandas as pd
import matplotlib.pyplot as mpl
from numpy import arange, array, histogram, isfinite, all, abs, log10
import DistributionDifference as DD
from bisect import bisect
from scipy.stats import scoreatpercentile, chi2_contingency
import setcolor

cyc_of_interest = 'cyc13'
eps = .1

read_table_args = dict(keep_default_na=False, na_values='---', index_col=0)

bind_dist = 1000

exp_cutoff = 5 + eps
expr_diff_weight =.01

def earth_mover_multi(points1, points2):
    dist = 0.0
    cyc1 = {col.split('_')[0] for col in points1.index}
    cyc2 = {col.split('_')[0] for col in points2.index}
    sums = [[],[]]
    for cyc in cyc1.intersection(cyc2):
        reps1 = {col.split('sl')[0] for col in
                 points1.select(startswith(cyc)).index
                }
        reps2 = {col.split('sl')[0] for col in
                 points2.select(startswith(cyc)).index
                }
        for rep1 in reps1:
            for rep2 in reps2:

                dist += (DD.earth_mover(points1.select(startswith(rep1)),
                                       points2.select(startswith(rep2)))**2
                         / (len(reps1)*len(reps2)))
        sums[0].append(points1.select(startswith(cyc)).mean())
        sums[1].append(points2.select(startswith(cyc)).mean())
    dist += DD.earth_mover(sums[0], sums[1])
    return dist**.5

dist_func = earth_mover_multi


try:
    wt = locals()['wt']
    assert min(wt.min()) == eps
except (KeyError, AssertionError):
    print("Importing data")
    keep_old = False
    wt = pd.read_table('prereqs/WT6.01.summary.tsv', **read_table_args)
    zld = pd.read_table('prereqs/Zld6.01.summary_merged.tsv', **read_table_args)
    bcd = pd.read_table('analysis/summary.tsv', **read_table_args)
    for dataset in (wt, zld, bcd):
        dataset.ix[:,:] += eps
        for i, column in enumerate(dataset.columns):
            if not all(isfinite(dataset.ix[:, column])):
                try:
                    dataset.ix[:, column] = (dataset.ix[:, i-1] + dataset.ix[:, i+1])/2
                except:
                    print("Whoops! {}".format(column))


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
    bcd_zld= pd.Series(index=wt_zld.index,
                       data = 0)

    print("Calculating Distances")
    zld_reps = {col.split('sl')[0]
                for col in zld.columns
                if col.startswith(cyc_of_interest)
               }
    for zld_rep in zld_reps:
        print('Calculating distance for ', zld_rep)
        for gene in wt_zld.index:
            wt_zld.ix[gene] += (
                dist_func(wt.ix[gene].select(startswith(cyc_of_interest)),
                          zld.ix[gene].select(startswith(zld_rep)))
            )
            wt_zld.ix[gene] += abs(log10(
                wt.ix[gene].select(startswith(cyc_of_interest)).mean()/
                zld.ix[gene].select(startswith(zld_rep)).mean()
            )) * expr_diff_weight
    wt_zld /= len(zld_reps)

    bcd_reps = {col.split('sl')[0]
                for col in bcd.columns
                if col.startswith(cyc_of_interest)
               }
    for bcd_rep in bcd_reps:
        print('Calculating distance for ', bcd_rep)
        for gene in wt_bcd.index:
            wt_bcd.ix[gene] += (
                dist_func(wt.ix[gene].select(startswith(cyc_of_interest)),
                          bcd.ix[gene].select(startswith(bcd_rep)))
            )
            wt_bcd.ix[gene] += abs(log10(
                wt.ix[gene].select(startswith(cyc_of_interest)).mean()/
                bcd.ix[gene].select(startswith(bcd_rep)).mean()
            )) * expr_diff_weight
    wt_bcd /= len(bcd_reps)

    for bcd_rep in bcd_reps:
        for zld_rep in zld_reps:
            print('Calculating distance for ', bcd_rep, zld_rep)
            for gene in wt_bcd.index:
                bcd_zld.ix[gene] += (
                    dist_func(zld.ix[gene].select(startswith(zld_rep)),
                              bcd.ix[gene].select(startswith(bcd_rep)))
                )
                bcd_zld.ix[gene] += abs(log10(
                    zld.ix[gene].select(startswith(zld_rep)).mean()/
                    bcd.ix[gene].select(startswith(bcd_rep)).mean()
                )) * expr_diff_weight
    bcd_zld /= (len(bcd_reps)*len(zld_reps))

    keep_old = True



both = wt_zld.dropna().index.intersection(wt_bcd.dropna().index)
wt_hi = wt.select(startswith(cyc_of_interest), axis=1).max(axis=1) > exp_cutoff
bcd_hi = bcd.select(startswith(cyc_of_interest), axis=1).max(axis=1) > exp_cutoff
keep = (wt_hi + bcd_hi) > 0

print("Plotting")

mpl.clf()
mpl.hist(wt_zld.dropna(), bins=arange(0, 1, .01), normed=True, label='WT vs Zld', histtype='step', color='r')
mpl.hist(wt_bcd.dropna(), bins=arange(0, 1, .01), normed=True, label='WT vs Bcd', histtype='step', color='g')
mpl.legend()
setcolor.set_foregroundcolor(mpl.gca(), 'w')
setcolor.set_backgroundcolor(mpl.gca(), 'b')
mpl.savefig('analysis/results/WTZldBcdHist.png', dpi=600, transparent=True)
mpl.savefig('analysis/results/WTZldBcdHist.eps', transparent=True)

colors = pd.Series('y g b c'.split())
cmaps = pd.Series(index=wt.index, data=0, dtype=int)
for gene in cmaps.index:
    for i, has_bind in enumerate([has_zld, has_bcd]):
            cmaps.ix[gene] += 2**i * (gene in has_bind)

cmaps.sort(ascending=False)
cmaps2 = colors[cmaps]
cmaps2.index = cmaps.index

keep = keep.ix[cmaps.index]

yy = scoreatpercentile(wt_bcd, 80)
xx = scoreatpercentile(wt_zld, 80)
nochgx = scoreatpercentile(wt_zld, 50)
nochgy = scoreatpercentile(wt_bcd, 50)
mpl.clf()

import matplotlib.path as mpath
import matplotlib.patches as mpatches

gene_specific_point = 1/2
both_changed_point = 1/2
gene_specific_slope = 1/3
both_changed_slope = 1/2

# How far out to go in the distance distribution. 1 ought to be fine, but
# different distance metrics can have higher than this.
extent = 10
bcd_only = mpath.Path([[0, yy], [xx * gene_specific_point, yy],
                       [gene_specific_slope*(extent - yy) + gene_specific_point * xx, extent],
                       [0, extent], [0, yy]])
zld_only = mpath.Path([[xx, 0], [xx, yy * gene_specific_point],
                       [extent, gene_specific_slope * (extent - xx) + gene_specific_point * yy],
                       [extent, 0], [xx, 0]])
both_change = mpath.Path([[xx * both_changed_point, yy],
                          [both_changed_slope * (extent - yy) + both_changed_point * xx, extent],
                          [extent,extent],
                          [extent, both_changed_slope * (extent - xx) + both_changed_point * yy],
                          [xx, both_changed_point * yy], [xx, yy],
                          [xx * both_changed_point, yy]])
no_change = mpath.Path([[0, 0], [0, nochgy], [nochgx, nochgy], [nochgx, 0],
                        [0, 0]])

wt_zld_exp = wt_zld.ix[cmaps.index].ix[keep]
wt_bcd_exp = wt_bcd.ix[cmaps.index].ix[keep]
bcd_zld_exp = bcd_zld.ix[cmaps.index].ix[keep]
cmaps_exp =  cmaps .ix[cmaps.index].ix[keep]
dist2 = bcd_zld_exp - abs(wt_zld_exp - wt_bcd_exp)
dist2_sorted = dist2.copy()
dist2_sorted.sort()

ax = mpl.gca()

mpl.plot([0,1], [0,1], 'r:', zorder=5)
mpl.plot([0, xx, xx], [yy, yy, 0], 'k-', zorder=5, alpha=0.4)
for c in range(4):
    mpl.scatter(x=wt_zld_exp.ix[(cmaps_exp == c)],
                y=wt_bcd_exp.ix[(cmaps_exp == c)],
                c='w',
                marker = '.',
                edgecolors='none',
                zorder=c+1,
                alpha=0.999999999)
mpl.xlabel('WT vs Zld')
mpl.ylabel('WT vs Bcd')
mpl.xlim(0, max(1, max(wt_zld_exp)*1.1))
mpl.ylim(0, max(1, max(wt_bcd_exp)*1.1))
mpl.title('{cyc} - {dist:3.1f}kb'
          .format(cyc=cyc_of_interest,
                  dist=bind_dist/1000.))
ax.set_aspect(1)
setcolor.set_foregroundcolor(mpl.gca(), 'w')
setcolor.set_backgroundcolor(mpl.gca(), 'b')
mpl.gcf().tight_layout()
mpl.savefig('analysis/results/WTBcdWTZldCorr-{}.png'.format(cyc_of_interest),
            dpi=600, transparent=True)
ax.add_patch(mpatches.PathPatch(bcd_only, facecolor='b', alpha=1)).set_zorder(0)
ax.add_patch(mpatches.PathPatch(zld_only, facecolor='g', alpha=1)).set_zorder(0)
ax.add_patch(mpatches.PathPatch(both_change, facecolor='c', alpha=1)).set_zorder(0)
ax.add_patch(mpatches.PathPatch(no_change, facecolor='r', alpha=1)).set_zorder(0)
mpl.savefig('analysis/results/WTBcdWTZldCorr-{}-withpatches.png'.format(cyc_of_interest),
            dpi=600, transparent=True)

mpl.clf()
ax = mpl.gca()
mpl.scatter(x=wt_zld_exp.ix[dist2_sorted.index],
            y=wt_bcd_exp.ix[dist2_sorted.index],
            c=dist2.ix[dist2_sorted.index],
            edgecolors=(0., 0., 0., 0.))
mpl.xlabel('WT vs Zld')
mpl.ylabel('WT vs Bcd')
ax.add_patch(mpatches.PathPatch(bcd_only, facecolor=(0., 0., 1., 0.1),
                                edgecolor=(.5, .5, .5, 1))).set_zorder(5)
ax.add_patch(mpatches.PathPatch(zld_only, facecolor=(0., 1., 0., 0.1),edgecolor=(.5, .5, .5, 1) )).set_zorder(5)
ax.add_patch(mpatches.PathPatch(both_change, facecolor=(0., 1., 1., 0.1), edgecolor=(.5, .5, .5, 1))).set_zorder(5)
ax.add_patch(mpatches.PathPatch(no_change, facecolor=(1., 0., 0., 0.1), edgecolor=(.5, .5, .5, 1))).set_zorder(5)
mpl.xlim(0, max(1, max(wt_zld_exp)*1.1))
mpl.ylim(0, max(1, max(wt_bcd_exp)*1.1))
ax.set_aspect(1)
mpl.colorbar()
setcolor.set_foregroundcolor(mpl.gca(), 'w')
setcolor.set_backgroundcolor(mpl.gca(), 'b')
mpl.savefig('analysis/results/WT_Zld_Bcd_3way.png', dpi=600, transparent=True)


in_bcd = bcd_only.contains_points(zip(wt_zld_exp, wt_bcd_exp))
in_zld = zld_only.contains_points(zip(wt_zld_exp, wt_bcd_exp))
in_both = both_change.contains_points(zip(wt_zld_exp, wt_bcd_exp))
in_nochg =  no_change.contains_points(zip(wt_zld_exp, wt_bcd_exp))

in_category = pd.Series(index=wt_zld_exp.index, data='')
for cat, genes in zip(['Bcd Only', 'Zld Only', 'Both', 'No Change'],
                      [in_bcd, in_zld, in_both, in_nochg]):
    in_category.ix[genes] = cat

in_category.to_csv('analysis/results/change_cats.txt', sep='\t')

contingency = array(
    [histogram(cmaps_exp.ix[in_bcd],  bins=arange(5))[0],
     histogram(cmaps_exp.ix[in_zld],  bins=arange(5))[0],
     histogram(cmaps_exp.ix[in_both], bins=arange(5))[0]])

print(contingency)
print(chi2_contingency(contingency))

all_dists = pd.DataFrame(data={'WT_Zld': wt_zld_exp, 'WT_Bcd':wt_bcd_exp,
                               'Bcd_Zld': bcd_zld_exp, 'Dist2' : dist2},
                        columns=['WT_Zld', 'WT_Bcd', 'Bcd_Zld', 'Dist2'])

all_dists.ix[in_bcd].sort(columns='Dist2').to_csv('analysis/results/change_bcd.tsv', sep='\t')

all_dists.ix[in_zld].sort(columns='Dist2').to_csv('analysis/results/change_zld.tsv', sep='\t')

all_dists.ix[in_both].sort(columns='Dist2').to_csv('analysis/results/change_both.tsv', sep='\t')

