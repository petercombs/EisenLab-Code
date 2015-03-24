from __future__ import print_function, division
import pandas as pd
from Utils import sel_startswith, sel_contains
from sys import argv
from scipy.stats import fisher_exact

if __name__ == "__main__":
    expr_min = 5
    eps = 1
    read_table_args = dict(index_col=0,
                           keep_default_na=False,
                           na_values=['---', ''])

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

    cdt = pd.read_table(argv[1], index_col='NAME',
                        keep_default_na=False,
                        na_values=[''])

    hannon = {line.strip() for line in open('prereqs/HannonGenes.txt')}

    ranges = {
        #('CG12398', 'Capa'),
        'ant1': ( 'gcm', 'CG10264'),
        'AntInZld_2': ( 'Cad88C', 'erm'),
        'both_poles': ( 'CG12986', 'tll'),
        #('CG17801', 'LysB'),
        'post1': ( 'CG2930', 'cas'),
        'RNAi_up_2': ( 'EfSec', 'modSP'),
        #('Esp', 'CG30286'),
        'RNAi_up_3': ( 'Esp', 'CG42808'),
        'tophalf': ( 'GstE2', 'Trissin'),
        #('CG12011', 'Obp49a'),
        'AntInZld_1': ( 'CG13712', 'bnb'),
    }

    in_range = all_expr.ix[cdt.index].index
    in_both = hannon.intersection(in_range)
    print(len(hannon))
    print("{:10s}: {}, {:.1%} of {}".format(
            'all',
            len(in_both), len(in_both)/len(in_range), len(in_range))
        )
    print("{:10s}: {:>5}, {:>5} of {:>5},\t({:^10.05})\t\t{}".format(
        'class', 'n', '%', 'N', 'pval', 'genes'))
    print('-'*60)
    for name, (gene1, gene2) in ranges.iteritems():
        in_range = all_expr.ix[cdt.index].ix[gene1:gene2].index
        in_both = hannon.intersection(in_range)
        n_in_both = len(in_both)
        in_either = len(hannon.union(in_range))
        print("{:10s}: {:5}, {:5.1%} of {:5},\t({:^10.05g})\t\t{}".format(
            name,
            (n_in_both), (n_in_both)/len(in_range), len(in_range),
            fisher_exact([[n_in_both, len(in_range)-n_in_both],
                          [100-n_in_both, len(all_expr) - in_either]])[1],
            ' '.join(in_both)
        )
        )
