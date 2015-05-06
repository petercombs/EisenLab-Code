from Utils import load_to_locals, sel_startswith
from numpy import linalg as la
import numpy as np
import pandas as pd
import PlotUtils as pu


if __name__ == "__main__":

    expr_min = 15
    if 'all_expr' not in locals():
        all_expr, (wt, bcd, zld, g20, hb), (wts, bcds, zlds, g20s, hbs), by_cycle\
        = load_to_locals(locals(), expr_min)

    time_points = ('WT_cyc13 '
                   'G20_cyc13_rep1 G20_cyc13_rep2 '
                   'bcd_cyc13_rep1 bcd_cyc13_rep2 '
                   'zld_cyc13_rep1 zld_cyc13_rep2 '
                   'hb_cyc13_rep1 '
                   'WT14D '
                   'G20_cyc14D_rep1 '
                   'bcd_cyc14D_rep1 bcd_cyc14D_rep2'
                   'zld_cyc14B zld_cyc14D'
                   'hb_cyc14D_rep1 hb_cyc14D_rep2').split()

    expr_points = all_expr.select(**sel_startswith(tuple(time_points)))
    svd = la.svd((expr_points
                             .subtract(expr_points.mean(axis=1), axis=0)
                             .dropna(how='all', axis=1)
                             .divide(expr_points.std(axis=1), axis=0)
                             .T))
    svd_df = pd.DataFrame(data=svd[0].T,
                          columns=expr_points.dropna(how='all', axis=0).columns)
    svd_signs = np.sign(svd_df.select(**sel_startswith('WT_cyc13')).mean(axis=1))

    kwargs = dict(
        norm_rows_by='center0all',
        progress_bar=True,
        col_sep='_sl',
        total_width=60,
        box_height=25,
        split_columns=True,
        draw_box=True,
        draw_row_labels=True,
        draw_name=True,
        make_hyperlinks=True,
        convert=True,
    )

    labels = ['{} -- {:.2%}'.format(i+1, pve)
              for i, pve in enumerate(svd[1]**2/sum(svd[1]**2))]
    pu.svg_heatmap(tuple(svd_df.select(**sel_startswith(t))
                         .multiply(svd_signs, axis=0)
                         .ix[:19]
                         for t in time_points),
                   'analysis/results/corr_svd.svg',
                   row_labels=labels[:20], cmap=pu.cm.RdBu,
                   **kwargs)


