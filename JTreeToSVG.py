import pandas as pd
import PlotUtils as pu
from Utils import sel_startswith, sel_contains
from sys import argv
from os import path

def cmap_by_prefix(prefix):

    cms = dict(
        WT = pu.ISH_CMS_5[0],
        bcd = pu.ISH_CMS_5[1],
        zld = pu.ISH_CMS_5[2],
        G20 = pu.ISH_CMS_5[3],
        hb = pu.ISH_CMS_5[4],
    )
    for p in cms:
        if prefix.startswith(p):
            return cms[p]
    return pu.ISH


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
            cyc_embs = {}
            for cyc in cycs:
                cyc_embs[cyc] = sub_df.select(**sel_contains(cyc))
            locals()[sub_df_name+'s'] = cyc_embs
    print("Read expression in")

    cdt = pd.read_table(argv[1], index_col='NAME',
                        keep_default_na=False,
                        na_values=[''])

    columns = (
        'WT_cyc13',
        'WT_cyc14D',
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


    ranges = [
        #('CG12398', 'Capa'),
        ('gcm', 'CG10264'),
        ('Cad88C', 'erm'),
        ('CG12986', 'tll'),
        #('CG17801', 'LysB'),
        ('CG2930', 'cas'),
        ('EfSec', 'modSP'),
        #('Esp', 'CG30286'),
        ('Esp', 'CG42808'),
        ('GstE2', 'Trissin'),
        #('CG12011', 'Obp49a'),
        ('CG13712', 'bnb'),
    ]

    for gene1, gene2 in ranges:
        outname = path.join(path.dirname(argv[1]),
                            'table_{}-{}.svg'.format(gene1, gene2))

        pu.svg_heatmap(
            data=all_expr.select(**sel_startswith(columns)).ix[cdt.index].ix[gene1:gene2],
            filename=outname,
            norm_rows_by='max',
            progress_bar=True,
            col_sep='_sl',
            total_width=150,
            box_height=30,
            split_columns=True,
            draw_box=True,
            draw_row_labels=True,
            draw_name=True,
            cmap_by_prefix=cmap_by_prefix,
        )

