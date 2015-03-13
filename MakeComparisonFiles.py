import pandas as pd
import PlotUtils
from progressbar import ProgressBar
import os
from os.path import join
from numpy import nan
from Utils import sel_startswith, contains, sel_contains
from JTreeToSVG import cmap_by_prefix

startswith = lambda x: lambda y: y.startswith(x)

def chunks(l, n):
    """ Return successive n-sized chunks from l.
    """
    return [l[i:i+n] for i in xrange(0, len(l), n)]

def make_comparison_file(gene, force=False):
    if not force and (max(wt.ix[gene].dropna()) < 3 and
                      max(bcd.ix[gene].dropna() < 3)):
        return
    if gene is nan:
        return
    data = (wt .select(contains('cyc13'),       axis=1).ix[[gene]],
            wt .select(contains('cyc14D'),      axis=1).ix[[gene]],
            bcd.select(contains('cyc13_rep1'),  axis=1).ix[[gene]],
            bcd.select(contains('cyc14D_rep1'), axis=1).ix[[gene]],
            bcd.select(contains('cyc13_rep2'),  axis=1).ix[[gene]],
            bcd.select(contains('cyc14D_rep2'), axis=1).ix[[gene]],
            g20.select(contains('cyc13_rep1'),  axis=1).ix[[gene]],
            g20.select(contains('cyc14D_rep1'), axis=1).ix[[gene]],
            g20.select(contains('cyc13_rep2'),  axis=1).ix[[gene]],
            None,
            zld.select(contains('cyc13_rep1'),  axis=1).ix[[gene]],
            zld.select(contains('cyc14D'),      axis=1).ix[[gene]],
            zld.select(contains('cyc13_rep3'),  axis=1).ix[[gene]],
            None,
            None,
            hb.select(**sel_contains('cyc14D_rep1')).ix[[gene]],
            None,
            hb.select(**sel_contains('cyc14D_rep2')).ix[[gene]],

           )
    max_13s = max(float(dp.max(axis=1))*1.1 +10 for dp in data[0::2] if dp is
                  not None)
    max_14s = max(float(dp.max(axis=1))*1.1 +10 for dp in data[1::2] if dp is
                  not None)
    normer = ((max_13s, max_14s)*5)[:-1]
    normer = (wt.ix[gene].select(contains('cyc13'  )).max() * 1.1 + 10,
              wt.ix[gene].select(contains('cyc14D' )).max() * 1.1 + 10,
              bcd.ix[gene].select(contains('cyc13' )).max() * 1.1 + 10,
              bcd.ix[gene].select(contains('cyc14D')).max() * 1.1 + 10,
              bcd.ix[gene].select(contains('cyc13' )).max() * 1.1 + 10,
              bcd.ix[gene].select(contains('cyc14D')).max() * 1.1 + 10,
              g20.ix[gene].select(contains('cyc13' )).max() * 1.1 + 10,
              g20.ix[gene].select(contains('cyc14D')).max() * 1.1 + 10,
              g20.ix[gene].select(contains('cyc13' )).max() * 1.1 + 10,
              g20.ix[gene].select(contains('cyc14D')).max() * 1.1 + 10,
              zld.ix[gene].select(contains('cyc13' )).max() * 1.1 + 10,
              zld.ix[gene].select(contains('cyc14D')).max() * 1.1 + 10,
              zld.ix[gene].select(contains('cyc13' )).max() * 1.1 + 10,
              zld.ix[gene].select(contains('cyc14D')).max() * 1.1 + 10,
              hb.ix[gene].select(contains('cyc13'  )).max() * 1.1 + 10,
              hb.ix[gene].select(contains('cyc14D' )).max() * 1.1 + 10,
              hb.ix[gene].select(contains('cyc13'  )).max() * 1.1 + 10,
              hb.ix[gene].select(contains('cyc14D' )).max() * 1.1 + 10,
             )
    names = [(stage + '- Max Expr {:.1f}'.format(float(dp.max(axis=1)))
              if dp is not None else '')
             for stage, dp in zip(('WT13 WT14D '
                                   'Bcd-13 Bcd-14D Bcd-13 Bcd-14D '
                                   'G20-13 G20-14D G20-13 G20-14D '
                                   'Zld-13 Zld-14D Zld-13 NONE '
                                   'NONE Hb-14D NONE Hb-14D').split(),
                                  data)
            ]
    # Use tuple multiplication for the correct number of each sample
    PlotUtils.svg_heatmap(data,
                          filename=join(outdir,
                                        gene + '.svg'),
                          cmap_by_prefix=cmap_by_prefix,
                          norm_rows_by=normer,
                          draw_box=True,
                          draw_name=True,
                          data_names=names,
                          spacers=10,
                          total_width=230,
                          box_height=130,
                          max_width=475,
                          col_sep='_sl',
                         )
    return join(outdir, gene+'.svg')

read_table_args = dict(index_col=0,
                       keep_default_na=False,
                       na_values=['---', ''])
all_expr = (pd.read_table('analysis/summary.tsv', **read_table_args)
            .sort_index())
wt  = all_expr.select(**sel_startswith('WT'))
bcd = all_expr.select(**sel_startswith('bcd'))
zld = all_expr.select(**sel_startswith('zld'))
g20 = all_expr.select(**sel_startswith('G20'))
hb  = all_expr.select(**sel_startswith('hb'))

assert all(wt.index == bcd.index)
assert all(wt.index == zld.index)
assert all(zld.index == bcd.index)
assert all(g20.index == wt.index)

outdir='analysis/results/svgs-withg20'

if __name__ == "__main__":
    try:
            os.makedirs(outdir)
    except OSError:
            pass

    from multiprocessing import Pool
    pbar = ProgressBar()
    pool = Pool()
    for genes_list in pbar(chunks(wt.index, 100)):
        pool.map(make_comparison_file, genes_list)

