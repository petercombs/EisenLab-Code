from __future__ import division
import pandas as pd
import DistributionDifference as dd
import PlotUtils as pu
import BindUtils as bu
from numpy import log2
from Utils import sel_startswith, load_to_locals, center_of_mass
from multiprocessing import Pool
from matplotlib.pyplot import (violinplot, ylabel, savefig, xticks,
                               close)
from sys import argv
from scipy.stats import chi2_contingency
from progressbar import ProgressBar as pb
from Literature import all_hits

columns = (
    'WT_cyc13',
    'WT_cyc14A',
    'WT_cyc14C',
    'WT_cyc14D',
    'WT_cyc14E',
    'bcd_cyc13_rep1',
    'bcd_cyc13_rep2',
    'bcd_cyc14D_rep1',
    'bcd_cyc14D_rep2',
    'zld_cyc13_rep2',
    'zld_cyc13_rep3',
    'zld_cyc14D',
    'zld_cyc14B',
    'G20_cyc13_rep1',
    'G20_cyc13_rep2',
    'G20_cyc14D_rep1',
    'hb_cyc13_rep1',
    'hb_cyc14D_rep1',
    'hb_cyc14D_rep2',
)

def get_dists_mp(args):
    gene, row = args
    temp = pd.Series(data=1, index=row.index)
    return dd.earth_mover_multi_rep(row+eps, temp)

if __name__ == "__main__":
    expr_min = 15
    eps = 1
    max_diff = .04
    read_table_args = dict(index_col=0,
                           keep_default_na=False,
                           na_values=['---', ''])


    if 'all_expr' not in locals():
        all_expr, (wt, bcd, zld, g20, hb), (wts, bcds, zlds, g20s, hbs), by_cycle\
        = load_to_locals(locals(), expr_min)
    print("Read expression in")


    p = Pool()
    diff_from_ubiq = {}
    for genotype in 'wt bcd g20 zld hb'.split():
        diff_from_ubiq[genotype] = pd.Series(
            data=p.map(get_dists_mp,
                       by_cycle[genotype]['cyc14D'].iterrows()),
            index=all_expr.index
        )

    p.close()
    ubiq_all = ((diff_from_ubiq['wt'] < max_diff)
                & (diff_from_ubiq['bcd'] < max_diff)
                & (diff_from_ubiq['g20'] < max_diff))

    become_ubiq = {}
    set_become_ubiq = {}
    lose_ubiq = {}
    set_lose_ubiq = {}
    lose_pattern = {}
    set_lose_pattern = {}
    gain_pattern = {}
    set_gain_pattern={}
    for genotype in 'bcd g20 zld hb'.split():
        become_ubiq[genotype] = (
            (diff_from_ubiq['wt'] > 2*max_diff)
            & (diff_from_ubiq[genotype] < max_diff)
            & (by_cycle[genotype]['cyc14D'].max(axis=1) > expr_min)
            & (by_cycle['wt']['cyc14D'].max(axis=1) > expr_min)
        )
        lose_ubiq[genotype] = (
            (diff_from_ubiq['wt'] < max_diff)
            & (diff_from_ubiq[genotype] > 2*max_diff)
            & (by_cycle[genotype]['cyc14D'].max(axis=1) > expr_min)
            & (by_cycle['wt']['cyc14D'].max(axis=1) > expr_min)
        )
        lose_pattern[genotype] = (
            (diff_from_ubiq['wt'] > 2*max_diff)
            & (by_cycle[genotype]['cyc14D'].max(axis=1) < expr_min)
            & (by_cycle['wt'    ]['cyc14D'].max(axis=1) > expr_min)
        )
        gain_pattern[genotype] = (
            (diff_from_ubiq[genotype] > 2*max_diff)
            & (by_cycle[genotype]['cyc14D'].max(axis=1) > expr_min)
            & (by_cycle['wt'    ]['cyc14D'].max(axis=1) < expr_min)
        )
        for set_df, df in ((set_become_ubiq, become_ubiq),
                           (set_lose_ubiq, lose_ubiq),
                           (set_gain_pattern, gain_pattern),
                           (set_lose_pattern, lose_pattern)):
            set_df[genotype] = {gene for gene in df[genotype].index
                                if df[genotype][gene]}

    change_list = pd.DataFrame(
        {label:{k:the_dict[k] for k in the_dict}
         for label, the_dict in
         (('Patterned to Low Expr', set_lose_pattern ),
          ('Low Expr to  Patterned', set_gain_pattern ),
          ('Uniform to  Patterned', set_lose_ubiq ),
          ('Patterned to  Uniform', set_become_ubiq),
         )
        })
    change = change_list.applymap(len)
    print()
    print(change.to_string())
    print()
    print(change.to_latex())

    up_both = []
    down_both = []
    thresh = 1.5
    wt_bcd_avg = wts['cyc14D'].mean(axis=1)/(eps + bcds['cyc14D'].mean(axis=1))
    wt_g20_avg = wts['cyc14D'].mean(axis=1)/(eps + g20s['cyc14D'].mean(axis=1))
    bcd_g20_avg = ((eps+bcds['cyc14D'].mean(axis=1))
                   /(eps+g20s['cyc14D'].mean(axis=1)))
    for gene in ubiq_all.ix[ubiq_all].index:
        if (wt_bcd_avg.ix[gene] < 1/thresh) and (wt_g20_avg.ix[gene] < 1/thresh):
            up_both.append(gene)
        elif ((wt_bcd_avg.ix[gene] > thresh)
              and (wt_g20_avg.ix[gene] > thresh)):
            down_both.append(gene)


    print("Up in both: {:>12} \nDown in both: {:>10}"
          .format(len(up_both), len(down_both)))

    print '-'*40
    print('|{:^5}|{:>8}|{:>8}|'
          .format('type', 'become', 'lose'))
    print '-'*40
    for genotype in 'bcd g20 zld hb'.split():
        print('|{:^5}|{:>8}|{:>8}|'
              .format(genotype,
                      len(set_become_ubiq[genotype]),
                      len(set_lose_ubiq[genotype]),
                     ))

    print '-'*40
    close('all')
    violinplot([log2(bcd_g20_avg.ix[up_both]),
                log2(bcd_g20_avg.ix[down_both])],
               showmedians=True, showextrema=False)
    xticks([1, 2], ['Up in both', 'Down in both'])
    ylabel('$\\log_2 bcd^-/bcd^{+++}$')
    savefig('analysis/results/updownboth_log2')

    binds = bu.get_binding_matrix(all_expr.index)
    # Ignore the Hannon&Wieschaus Bicoid Dosage ChIPseq data
    binds.drop(labels=['steplike', 'bcd_insensitive', 'bcd_sensitive'], axis=1,
              inplace=True)
    binds = binds.ix[:, bu.ap_early_zld]
    diff_bind = pd.DataFrame(index=change_list.index,
                             columns=change_list.columns,
                             dtype=object)
    for row in pb()(diff_bind.index):
        for col in diff_bind.columns:
            diff_bind.ix[row, col] = list()
            in_set = change_list.ix[row, col]
            not_in_set = binds.index.difference(pd.Index(in_set))
            for tf in binds.columns:
                fbind = binds.ix[:, tf]
                chi2, p, dof, expected = chi2_contingency(
                    [[sum(fbind.ix[in_set]),
                      sum(1-fbind.ix[in_set])
                     ],
                     [sum(fbind.ix[not_in_set]),
                      sum(1-fbind.ix[not_in_set])
                     ]]
                )
                if (p
                    * len(diff_bind.index)
                    * len(diff_bind.columns)
                    * len(binds.columns)
                    < .05):
                    diff_bind.ix[row, col].append((tf, p))


    print (diff_bind
           .applymap(lambda x: ','.join([i[0] for i in x]) or '---')
           .to_latex())

    if len(argv) > 1:
        #cdt = pd.read_table(argv[1], index_col='NAME',
                            #keep_default_na=False,
                            #na_values=[''])
        all_expr_cdt = (all_expr
                        .ix[center_of_mass(by_cycle['wt']['cyc14D'])
                            .sort(inplace=False).index])
        for genotype in set_lose_ubiq:
            kwargs = dict(
                norm_rows_by='max',
                progress_bar=True,
                col_sep='_sl',
                total_width=40,
                box_height=10,
                split_columns=True,
                draw_box=True,
                draw_row_labels=True,
                color_row_labels=all_hits,
                draw_name=True,
                cmap_by_prefix=pu.cmap_by_prefix,
                make_hyperlinks=True,
                convert=True,
            )
            pu.svg_heatmap(
                data=(all_expr_cdt
                      .ix[center_of_mass(by_cycle[genotype]['cyc14D'])
                          .sort(inplace=False).index]
                      .select(**sel_startswith(columns))
                      .select(set_lose_ubiq[genotype].__contains__)),
                filename='analysis/results/{}_lose_ubiq.svg'.format(genotype),
                **kwargs)

            pu.svg_heatmap(
                data=(all_expr_cdt
                      .select(**sel_startswith(columns))
                      .select(set_become_ubiq[genotype].__contains__)),
                filename='analysis/results/{}_patterned_to_ubiq.svg'.format(genotype),
                **kwargs)

            kwargs['norm_rows_by'] = (by_cycle['wt']['cyc14D']
                                      .ix[all_expr_cdt.index]
                                      .select(set_become_ubiq[genotype].__contains__)
                                      .max(axis=1)+10)
            pu.svg_heatmap(
                data=(all_expr_cdt
                      .select(**sel_startswith(columns))
                      .select(set_become_ubiq[genotype].__contains__)),
                filename='analysis/results/{}_wtnorm_patterned_to_ubiq.svg'.format(genotype),
                **kwargs)




