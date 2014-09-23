from sklearn.decomposition import PCA
import PlotUtils
import pandas as pd
from matplotlib import cm
from BindUtils import get_binding_matrix

read_table_args = {'index_col': 0, 'keep_default_na': False,
                   'na_values': ['---', '']}

wt = pd.read_table('prereqs/WT6.01.summary.tsv', **read_table_args).sort_index()

zld = pd.read_table('prereqs/Zld6.01.summary_merged.tsv',
                    **read_table_args).sort_index()

bcd = pd.read_table('analysis/summary.tsv', **read_table_args).sort_index()

prepend = lambda y: lambda x: y+"_"+x
startswith = lambda y: lambda x: x.startswith(y)

change_list = pd.read_table('analysis/results/change_both.tsv',
                            **read_table_args)
both_bm = get_binding_matrix(change_list.index)


wt14 = wt.select(startswith('cyc14D'), axis=1).dropna(axis=1, how='all')
bcd14 = bcd.select(startswith('cyc14D'), axis=1).dropna(axis=1, how='all')
zld14 = zld.select(startswith('cyc14D'), axis=1).dropna(axis=1, how='all')

p = PCA()

bind_and_exp = (both_bm
                .rename(columns=prepend('bind_sl'))
                .join([wt14.rename(columns=prepend('wt'))
                       .divide(.01+wt14.max(axis=1), axis=0),
                       bcd14.rename(columns=prepend('bcd'))
                       .divide(.01+bcd14.max(axis=1), axis=0),
                       zld14.rename(columns=prepend('zld'))
                       .divide(.01+zld14.max(axis=1), axis=0)],
                      how='inner')
               )
bind_and_exp.dropna(how='all', axis=1, inplace=True)
f = p.fit(bind_and_exp)
comps = pd.DataFrame(data=f.components_, columns=bind_and_exp.columns)

PlotUtils.svg_heatmap(comps, col_sep='_sl',
                      filename='analysis/results/pca.svg',
                      box_size=10, box_height=20,
                      row_labels=p.explained_variance_ratio_,
                      draw_row_labels=True,
                      norm_rows_by='center0',
                      cmap=cm.RdBu)
comps.to_csv('analysis/results/pca.csv')

bind_pcas = (comps.select(startswith('bind_sl'), axis=1)
             .rename(columns=lambda x: x.split('_')[-1])
            )
exp_pcas = comps.select(startswith(('wt', 'bcd', 'zld')), axis=1)

