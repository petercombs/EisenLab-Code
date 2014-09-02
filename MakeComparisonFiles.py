import pandas as pd
import PlotUtils
from progressbar import ProgressBar
from matplotlib.cm import Greens, Reds
import os
from os.path import join
from numpy import nan

startswith = lambda x: lambda y: y.startswith(x)

def chunks(l, n):
    """ Return successive n-sized chunks from l.
    """
    return [l[i:i+n] for i in xrange(0, len(l), n)]

def make_comparison_file(gene):
    if max(wt.ix[gene]) < 3 and max(bcd.ix[gene] < 3):
        return
    if gene is nan:
        return
    data = (wt .select(startswith('cyc13'),       axis=1).ix[[gene]],
            wt .select(startswith('cyc14D'),      axis=1).ix[[gene]],
            bcd.select(startswith('cyc13_rep1'),  axis=1).ix[[gene]],
            bcd.select(startswith('cyc14D_rep1'), axis=1).ix[[gene]],
            bcd.select(startswith('cyc13_rep2'),  axis=1).ix[[gene]],
            bcd.select(startswith('cyc14D_rep2'), axis=1).ix[[gene]],
            zld.select(startswith('cyc13_rep1'),  axis=1).ix[[gene]],
            zld.select(startswith('cyc14D'),      axis=1).ix[[gene]],
            zld.select(startswith('cyc13_rep3'),  axis=1).ix[[gene]],
           )
    max_13s = max(float(dp.max(axis=1))*1.1 +10 for dp in data[0::2])
    max_14s = max(float(dp.max(axis=1))*1.1 +10 for dp in data[1::2])
    PlotUtils.svg_heatmap(data,
                          filename=join(outdir,
                                        gene + '.svg'),
                          cmap=(None,None,
                                Greens, Greens,
                                Greens, Greens,
                                Reds, Reds,
                                Reds
                               ),
                          norm_rows_by=((max_13s, max_14s)*5)[:-1],
                          draw_box=True,
                          draw_name=True,
                          data_names=('WT13 WT14D Bcd-13 Bcd-14D Bcd-13'
                                     ' Bcd-14D Zld-13 Zld-14D Zld-13').split(),
                          total_width=400,
                          box_height=130,
                          max_width=840,
                         )

read_table_args = dict(index_col=0,
                       keep_default_na=False,
                       na_values=['---', ''])
wt  = (pd.read_table('prereqs/WT6.01.summary.tsv', **read_table_args)
              #.dropna(how='any')
              .sort_index())
bcd = (pd.read_table('analysis/summary.tsv', **read_table_args)
              #.dropna(how='any', axis=1)
              .sort_index())
zld = (pd.read_table('prereqs/Zld6.01.summary_merged.tsv', **read_table_args)
              #.dropna(how='any', axis=1)
              .sort_index())
assert all(wt.index == bcd.index)
assert all(wt.index == zld.index)
assert all(zld.index == bcd.index)

outdir='analysis/results/svgs-withzld'

try:
        os.makedirs(outdir)
except OSError:
        pass

from multiprocessing import Pool
pbar = ProgressBar()
pool = Pool()
for genes_list in pbar(chunks(wt.index, 100)):
    pool.map(make_comparison_file, genes_list)

