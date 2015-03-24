from __future__ import print_function, division
import pandas as pd
import matplotlib.pyplot as mpl
from numpy import arange, array, histogram, abs, median, mean
import DistributionDifference as DD
from bisect import bisect
from scipy.stats import scoreatpercentile, chi2_contingency, gaussian_kde
import setcolor
from collections import defaultdict
from Utils import sel_startswith, sel_contains, contains
from itertools import combinations
from sys import argv

cyc_of_interest = 'cyc14D'
eps = .1

read_table_args = dict(keep_default_na=False, na_values=['---', ''], index_col=0)

bind_dist = 1000

exp_cutoff = 5 + eps
expr_diff_weight =.01

def earth_mover_multi(points1, points2, abs_expr=False):
    dist = 0.0
    cyc1 = {col.split('_')[1] for col in points1.index}
    cyc2 = {col.split('_')[1] for col in points2.index}
    sums = [[],[]]
    for cyc in cyc1.intersection(cyc2):
        reps1 = {col.split('sl')[0] for col in
                 points1.select(contains(cyc)).index
                }
        reps2 = {col.split('sl')[0] for col in
                 points2.select(contains(cyc)).index
                }
        for rep1 in reps1:
            for rep2 in reps2:

                dist += (DD.earth_mover_interp(points1.select(contains(rep1)),
                                       points2.select(contains(rep2)))**2
                         / (len(reps1)*len(reps2)))
        if abs_expr:
            sums[0].append(points1.select(contains(cyc)).mean())
            sums[1].append(points2.select(contains(cyc)).mean())
    if abs_expr:
        dist += DD.earth_mover_interp(sums[0], sums[1])
    return dist**.5

dist_func = earth_mover_multi

try:
    tss_dict = locals()['tss_dict']
except KeyError:
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

def get_dists(gene):
    return (dist_func(set1.ix[gene], set2.ix[gene]),
            dist_func(set1.ix[gene], set3.ix[gene]),
            dist_func(set2.ix[gene], set3.ix[gene]),
           )

if __name__ == "__main__":

    screen = '--screen' in argv

    expr_min = 10
    eps = .1

    if 'all_expr' not in locals():
        all_expr = (pd.read_table('analysis/summary.tsv', **read_table_args)
                    .sort_index())
        top_expr = all_expr.max(axis=1)
        all_expr = all_expr.ix[top_expr > expr_min]
        wt  = all_expr.select(**sel_startswith('WT'))
        bcd = all_expr.select(**sel_startswith('bcd'))
        zld = all_expr.select(**sel_startswith('zld'))
        g20 = all_expr.select(**sel_startswith('G20'))
        hb  = all_expr.select(**sel_startswith('hb'))

        wts = bcds = zlds = g20s = hbs = 0
        for sub_df_name in 'wt bcd zld g20 hb'.split():
            sub_df = locals()[sub_df_name]
            cycs = {col.split('_sl')[0].split('_',1)[1] for col in sub_df.columns}
            cycs.update({col.split('_')[1] for col in sub_df.columns})
            cyc_embs = {}
            for cyc in cycs:
                cyc_embs[cyc] = sub_df.select(**sel_contains(cyc))
            locals()[sub_df_name+'s'] = cyc_embs
    print("Read expression in")

    combos = combinations( ['wt', 'bcd', 'zld', 'g20', 'hb'], 3)
    all_d2s = {}
    for set1_name, set2_name, set3_name in combos:
        print('-'*50)
        print(set1_name, set2_name, set3_name)
        print('-'*50)
        set1 = locals()[set1_name].select(**sel_contains(cyc_of_interest))
        set2 = locals()[set2_name].select(**sel_contains(cyc_of_interest))
        set3 = locals()[set3_name].select(**sel_contains(cyc_of_interest))

        max_1 = set1.max(axis=1)
        max_2 = set2.max(axis=1)
        max_3 = set3.max(axis=1)

        keep = median([max_1, max_2, max_3], axis=0) > expr_min

        set1 = set1.ix[keep] + eps
        set2 = set2.ix[keep] + eps
        set3 = set3.ix[keep] + eps


        dist_12 = pd.Series(index=set1.index, data=0)
        dist_13 = pd.Series(index=set1.index, data=0)
        dist_23 = pd.Series(index=set1.index, data=0)

        from multiprocessing import Pool, pool
        p = locals().get('p', Pool())
        if not isinstance(p, pool.Pool):
            p = Pool()
        dists = p.map(get_dists, dist_12.index)
        for gene, dists_for_gene in zip(dist_12.index, dists):
            dist_12.ix[gene] += dists_for_gene[0]
            dist_13.ix[gene] += dists_for_gene[1]
            dist_23.ix[gene] += dists_for_gene[2]

        xx = scoreatpercentile(dist_12, 80)
        yy = scoreatpercentile(dist_13, 80)
        nochgx = scoreatpercentile(dist_12, 50)
        nochgy = scoreatpercentile(dist_13, 50)
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
        set2_only = mpath.Path(
            [[0, yy],
             [xx * gene_specific_point, yy],
             [gene_specific_slope*(extent - yy) + gene_specific_point * xx, extent],
             [0, extent],
             [0, yy]
            ])
        set3_only = mpath.Path(
            [[xx, 0],
             [xx, yy * gene_specific_point],
             [extent, gene_specific_slope * (extent - xx) + gene_specific_point * yy],
             [extent, 0],
             [xx, 0]
            ])
        both_change = mpath.Path(
            [[xx * both_changed_point, yy],
             [both_changed_slope * (extent - yy) + both_changed_point * xx, extent],
             [extent,extent],
             [extent, both_changed_slope * (extent - xx) + both_changed_point * yy],
             [xx, both_changed_point * yy],
             [xx, yy],
             [xx * both_changed_point, yy]
            ])
        no_change = mpath.Path([[0, 0],
                                [0, nochgy],
                                [nochgx, nochgy],
                                [nochgx, 0],
                                [0, 0]])

        in_set2 = set2_only.contains_points(zip(dist_12, dist_13))
        in_set3 = set3_only.contains_points(zip(dist_12, dist_13))
        in_both = both_change.contains_points(zip(dist_12, dist_13))
        in_nochg =  no_change.contains_points(zip(dist_12, dist_13))

        obs = array([[sum(in_set2), sum(in_both)],
                             [sum(~(in_set2+in_set3+in_both)), sum(in_set3)]])
        chi2, p, dof, exp = chi2_contingency(obs)
        print(p)
        print(obs)
        print(exp)
        print(obs-exp)
        print()

        dist2 = dist_23 - abs(dist_12 - dist_13)
        dist2_sorted = dist2.copy()
        dist2_sorted.sort()

        dist2_95 = scoreatpercentile(dist2[in_both], 95)
        dist2_95 = .2
        dist2_hi = dist2.index[dist2 > dist2_95]
        dist2_95 = scoreatpercentile(dist2[in_both], 95)
        print("mean D2: {:.05%} \n95th pct: {:.05%}"
              .format(mean(dist2[in_both]), dist2_95))
        print("{} genes: {}".format(len(dist2_hi), ", ".join(dist2_hi)))

        ##################################
        ###         PLOTTING           ###
        ##################################

        mpl.close('all')

        bins = arange(0,0.2,.0033)
        mpl.hist(dist_12.dropna(), bins=bins,
                 #normed=True,
                 histtype='step',
                 label='{} vs {}'.format(set1_name, set2_name))
        mpl.hist(dist_13.dropna(), bins=bins,
                 #normed=True,
                 histtype='step',
                 label='{} vs {}'.format(set1_name, set3_name))
        mpl.hist(dist_23.dropna(), bins=bins,
                 #normed=True,
                 histtype='step',
                 label='{} vs {}'.format(set2_name, set3_name))
        mpl.legend()
        if screen:
            setcolor.set_foregroundcolor(mpl.gca(), 'w')
            setcolor.set_backgroundcolor(mpl.gca(), 'b')
        mpl.savefig('analysis/results/{}{}{}_Hist.png'.format(set1_name,
                                                              set2_name,
                                                              set3_name),
                    dpi=300, transparent=True)




        mpl.clf()
        mpl.scatter(x=dist_12.ix[dist2_sorted.index],
                    y=dist_13.ix[dist2_sorted.index],
                    c=dist2.ix[dist2_sorted.index],
                    edgecolors=(0., 0., 0., 0.))
        mpl.xlabel('{} vs {}'.format(set1_name, set2_name))
        mpl.ylabel('{} vs {}'.format(set1_name, set3_name))
        ax = mpl.gca()
        ax.add_patch(mpatches.PathPatch(set2_only,
                                        facecolor=(0., 0., 1., 0.1),
                                        edgecolor=(.5, .5, .5, 1))).set_zorder(5)
        ax.add_patch(mpatches.PathPatch(set3_only,
                                        facecolor=(0., 1., 0., 0.1),
                                        edgecolor=(.5, .5, .5, 1) )).set_zorder(5)
        ax.add_patch(mpatches.PathPatch(both_change,
                                        facecolor=(0., 1., 1., 0.1),
                                        edgecolor=(.5, .5, .5, 1))).set_zorder(5)
        ax.add_patch(mpatches.PathPatch(no_change,
                                        facecolor=(1., 0., 0., 0.1),
                                        edgecolor=(.5, .5, .5, 1))).set_zorder(5)
        mpl.xlim(0, min(1, max(dist_12)*1.1))
        mpl.ylim(0, min(1, max(dist_13)*1.1))
        ax.set_aspect(1)
        mpl.colorbar()
        if screen:
            setcolor.set_foregroundcolor(mpl.gca(), 'w')
            setcolor.set_backgroundcolor(mpl.gca(), 'b')
        mpl.savefig('analysis/results/{}{}{}_3way.png'.format(set1_name,
                                                              set2_name,
                                                              set3_name),
                    dpi=600,
                    transparent=True)

        mpl.clf()
        mpl.hist(dist2.ix[in_both], bins=bins, )
        if screen:
            setcolor.set_foregroundcolor(mpl.gca(), 'w')
            setcolor.set_backgroundcolor(mpl.gca(), 'b')
        mpl.savefig('analysis/results/{}{}{}_D2hist.png'.format(set1_name,
                                                                set2_name,
                                                                set3_name),
                    dpi=600,
                    transparent=True)

        ##################################
        ###        OUTPUT FILE         ###
        ##################################


        all_dists = pd.DataFrame(data={
            '{}_{}'.format(set1_name, set2_name): dist_12,
            '{}_{}'.format(set1_name, set3_name): dist_13,
            '{}_{}'.format(set2_name, set3_name): dist_23,
            'Dist2' : dist2
        },
            columns=[
                '{}_{}'.format(set1_name, set2_name),
                '{}_{}'.format(set1_name, set3_name),
                '{}_{}'.format(set2_name, set3_name),
                'Dist2'])
        all_dists.sort(columns='Dist2', inplace=True)
        all_dists.to_csv('analysis/results/{}{}{}_change.tsv'.format(set1_name,
                                                                     set2_name,
                                                                     set3_name),
                         sep='\t',
                        )
        all_d2s['_'.join([set1_name, set2_name, set3_name])] = dist2
    mpl.clf()
    fig1 = mpl.figure()
    fig2 = mpl.figure()
    for trio in all_d2s:
        if not trio.startswith('wt'): continue
        k = gaussian_kde(all_d2s[trio])
        fig1.gca().hist(all_d2s[trio],
                 bins=arange(0, .5, .0033),
                 normed=True,
                 label='{} 3-way distance'.format(trio),
                 histtype='step'
                )
        fig2.gca().plot(arange(0, .5, 1e-4),
                  k(arange(0, .5, 1e-4)),
                  label='{} 3-way distance'.format(trio),
                 )
    max_y = max(fig1.gca().get_ylim()[1], fig2.gca().get_ylim()[1])
    for f in (fig1, fig2):
        f.gca().set_xlabel('D2 score (difference of distances)')
        f.gca().set_ylabel('Density')
        f.gca().set_ylim(0, max_y)
        f.gca().legend()
    fig1.savefig('analysis/results/all_d2s.png', dpi=300)
    fig2.savefig('analysis/results/all_d2s_kde.png', dpi=300)
    mpl.close('all')
    assert False

    keep = 'TRACK DOWN THIS LINE'

    colors = pd.Series('y g b c'.split())
    cmaps = pd.Series(index=wt.index, data=0, dtype=int)
    for gene in cmaps.index:
        for i, has_bind in enumerate([has_zld, has_bcd]):
                cmaps.ix[gene] += 2**i * (gene in has_bind)

    cmaps.sort(ascending=False)
    cmaps2 = colors[cmaps]
    cmaps2.index = cmaps.index

    keep = keep.ix[cmaps.index]

    wt_zld = bcd_zld = wt_bcd = 0

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

