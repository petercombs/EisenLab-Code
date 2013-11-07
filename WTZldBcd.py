#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
from matplotlib import cm
import PlotUtils
from PeakFinder import make_sort_num


startswith = lambda x: lambda y: y.startswith(x)

def splitns(n, s):
    return " ".join(s[i:i+n] for i in range(0, len(s), n))

new_cols = (list(coarse.columns[:12]) +
            list(coarse.columns[17:11:-1]) +
            list(coarse.columns[18:23]) +
            list(coarse.columns[29:22:-1]) +
            list(coarse.columns[30:]))

coarse_sort = coarse.reindex(columns=new_cols)


c1 = coarse_sort.select(startswith('cyc14_rep1'), axis=1)
c2 = coarse_sort.select(startswith('cyc14_rep2'), axis=1)
c3 = coarse_sort.select(startswith('cyc14_rep3'), axis=1)
b1 = coarse_sort.select(startswith('cycbcd14_rep1'), axis=1)
b2 = coarse_sort.select(startswith('cycbcd14_rep2'), axis=1)
b3 = coarse_sort.select(startswith('cycbcd14_rep3'), axis=1)

sort_num = make_sort_num((zld13, zld14A, zld14B), (b1, b2, b3), (c1, c2, c3),
                        fold=2).dropna()
sort_num.sort()
sort_num = sort_num[::-1]
sort_bin = pd.Series(index=sort_num.index, data=[splitns(3, np.base_repr(int(i),4)) for i in sort_num])

norm_rows = coarse_sort.max(axis=1)

n = 200

PlotUtils.svg_heatmap((coarse_sort.ix[sort_num.index[:200]]
                             .select(startswith('cyc14_'), axis=1),
                       coarse_sort.ix[sort_num.index[:200]]
                             .select(startswith('cycbcd14_'), axis=1),
                       zld13.join(zld14A).join(zld14B).ix[sort_num.index[:200]],
                      ),
                      'analysis/bcd/bcd_compare.svg',
                      cmap=(cm.Blues, cm.Greens, cm.Reds),
                      norm_rows_by=norm_rows.ix[sort_num.index[:200]],
                      col_sep='sl',
                      total_width=200, box_height=7,
                      draw_row_labels=True,
                     )
