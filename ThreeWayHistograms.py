from __future__ import print_function, division
import pandas as pd
import matplotlib.pyplot as mpl
from numpy import arange, array, abs, median, mean
import DistributionDifference as DD
from scipy.stats import (scoreatpercentile, chi2_contingency, gaussian_kde,
                         spearmanr, pearsonr)
import setcolor
from Utils import load_to_locals, sel_contains, contains
from itertools import combinations
from sys import argv
from scipy.stats import linregress
import BindUtils as bu

cyc_of_interest = 'cyc14D'
eps = .1

read_table_args = dict(keep_default_na=False, na_values=['---', ''], index_col=0)

bind_dist = 1e4

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

#has_zld = bu.has_tfs['zld']
#has_bcd = bu.has_tfs['bcd']

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
        all_expr, (wt, bcd, zld, g20, hb), (wts, bcds, zlds, g20s, hbs), by_cycle\
        = load_to_locals(locals(), expr_min)
    print("Read expression in")

    combos = combinations( ['wt', 'bcd', 'zld', 'g20', 'hb'], 3)
    all_d2s = {}
    both_d2s = {}
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
        p.close()
        del p
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
            setcolor.set_screen(mpl.gcf())
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
            setcolor.set_screen(mpl.gcf())
        mpl.savefig('analysis/results/{}{}{}_3way.png'.format(set1_name,
                                                              set2_name,
                                                              set3_name),
                    dpi=600,
                    transparent=True)

        mpl.clf()
        mpl.hist(dist2.ix[in_both], bins=bins, )
        if screen:
            setcolor.set_screen(mpl.gcf())
        mpl.savefig('analysis/results/{}{}{}_D2hist.png'.format(set1_name,
                                                                set2_name,
                                                                set3_name),
                    dpi=600,
                    transparent=True)


        mpl.clf()
        x = dist2.ix[in_both]
        y = wts['cyc14D'].max(axis=1).ix[in_both]
        mpl.semilogy( x, y, '.')
        mpl.xlabel('$\\Delta D$')
        mpl.ylabel('Max FPKM')
        arr = pearsonr(x, y)
        rho = spearmanr(x, y)
        mpl.title('$r = {:0.4}{} ;\\ \\  \\rho={:0.4}{}$'
              .format(arr[0],
                      '^{*}' if arr[1]<0.05 else '',
                      rho[0],
                      '^{*}' if rho[1]<0.05 else '',
                     )
                 )
        if screen:
            setcolor.set_screen(mpl.gcf())
        mpl.savefig('analysis/results/{}{}{}_D2exprcorr.png'.format(set1_name,
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
        both_d2s['_'.join([set1_name, set2_name, set3_name])] = dist2.ix[in_both]


        binds = (
            bu.get_binding_matrix(dist2_sorted.index, dist=bind_dist)
            .ix[:, bu.ap_early_zld]
        )
        n = 200
        n_2 = 50
        n_mid = sum(in_both)/2
        hi_d2 = dist2_sorted.index[-50:]
        mid_d2 = (dist2_sorted
                  .ix[in_both]
                  .ix[int(n_mid - n/2):int(n_mid+n/2)]
                  .index)
        enriched_bind = pd.DataFrame(columns=['odds ratio', 'base freq',
                                              'p-value'])
        for tf in binds.columns:
            table = [
                [sum(binds.ix[hi_d2, tf]),
                 sum(1-binds.ix[hi_d2, tf])],
                [sum(binds.ix[mid_d2, tf]),
                 sum(1-binds.ix[mid_d2, tf])],
            ]
            try:
                chi2, p, dof, expected = chi2_contingency(table)
            except ValueError as e:
                print(e)
                print("Poor table for {}: {}".format(tf, table))
            if p * len(binds.columns) < 0.05:
                enriched_bind.ix[tf] = [
                    (table[0][0]/table[0][1])/(table[1][0]/table[1][1]),
                    sum(binds.ix[:, tf])/len(binds),
                    p
                ]
            else:
                print("Nonsig: {:>10}, {}, {}".format(tf, p, table))
        enriched_bind.sort('odds ratio', inplace=True, ascending=False)

        formatter = lambda x: lambda y: x.format(y)

        print(enriched_bind.to_latex(formatters=
                                     {
                                         'odds ratio':formatter('{:.3}'),
                                         'base freq':formatter('{:0.2%}'),
                                         'p-value':formatter('{:0.3g}'),
                                     }))
        window_size = 50
        print("Windows of {} genes".format(window_size))
        windows = [slice(i, i+window_size)
                   for i in range(0, len(dist2_sorted), window_size)]
        mpl.clf()
        for n in range(3, len(binds.columns)):
            xs = array([mean(dist2_sorted.ix[w]) for w in windows])
            ys = array([mean(binds.ix[w].sum(axis=1)>=n) for w in windows])
            slope, intercept, r_value, p_value, std_err = linregress(xs,ys)
            points, = mpl.plot(xs, ys, '.',
                               label='$\\geq$ {} TFs bound'.format(n))
            mpl.plot(xs * 15, slope * xs * 15 + intercept, color=points.get_color())
        mpl.legend(numpoints=1,
                   bbox_to_anchor=(1.1, 0),
                   loc='lower right',
                   borderaxespad=0.2)
        mpl.xlim(0, max(xs)*1.1)
        mpl.ylim(0, 1.1)
        mpl.xlabel('$\\Delta D$')
        mpl.ylabel('Fraction with level of binding')
        #mpl.tight_layout()
        if screen:
            setcolor.set_screen(mpl.gcf())
        mpl.savefig('analysis/results/{}{}{}-D2vsBind'.format(set1_name,
                                                              set2_name,
                                                              set3_name),
                   transparent=True, dpi=300)






    mpl.clf()
    for name, d2s in (('both', both_d2s), ('all', all_d2s)):
        fig1 = mpl.figure()
        fig2 = mpl.figure()
        for trio in d2s:
            if not trio.startswith('wt'): continue
            k = gaussian_kde(d2s[trio])
            fig1.gca().hist(d2s[trio],
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
        if screen:
            setcolor.set_screen(fig1)
            setcolor.set_screen(fig2)
        fig1.savefig('analysis/results/{}_d2s.png'.format(name),
                     dpi=300,transparent=True)
        fig2.savefig('analysis/results/{}_d2s_kde.png'.format(name),
                     dpi=300,transparent=True)
        #mpl.close('all')

