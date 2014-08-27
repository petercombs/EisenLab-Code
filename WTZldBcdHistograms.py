from __future__ import print_function, division
import pandas as pd
import matplotlib.pyplot as mpl
from numpy import arange
from os.path import join
import DistributionDifference as DD

read_table_args = dict(keep_default_na=False, na_values='---', index_col=0)

try:
    wt = locals()['wt']
except KeyError:
    wt = pd.read_table('prereqs/WT6.01.summary.tsv', **read_table_args)
    zld = pd.read_table('prereqs/Zld6.01.summary_merged.tsv', **read_table_args)
    bcd = pd.read_table('analysis/summary.tsv', **read_table_args)

try:
    zld_bind = locals()['zld_bind']
except:
    zldbind = pd.read_table('prereqs/journal.pgen.1002266.s005.xls', skiprows=1)
    has_zld = set(item.strip() for item in zldbind.TSS_gene)

binding_directory = "prereqs/BDTNP_in_vivo_binding_Release.2.1/Supplemental_Tables"

try:
    mapping = locals()['mapping']

except KeyError:
    mapping = pd.read_table('prereqs/gene_map_table_fb_2014_04.tsv', skiprows=5)

    print("Reading in mapping table...")
    fbgn_to_name = {row['primary_FBid'] : row['##current_symbol']
                    for i, row in mapping.iterrows()}
    name_to_fbgn = {row['##current_symbol']: row['primary_FBid']
                    for i, row in mapping.iterrows()}
    print("...done")

try:
    has_bcd = locals()['has_bcd']
except KeyError:
    from glob import glob
    bcd_bind_files = glob(join(binding_directory, 'bcd*'))
    has_bcd = {fbgn_to_name[row.split()[-2]]
               for fname in bcd_bind_files
               for row in open(fname)
               if row.split()[-2] in fbgn_to_name
              }

startswith = lambda x: lambda y: y.startswith(x)

try:
    assert locals()['keep_old']
    wt_zld = locals()['wt_zld']
    wt_bcd = locals()['wt_bcd']
except (KeyError, AssertionError):
    wt_zld = pd.Series(index=wt.index)
    wt_bcd = pd.Series(index=wt.index)



    print("Calculating Distances")
    for gene in wt.index.intersection(zld.index):
        wt_zld.ix[gene] = (DD.earth_mover(wt.ix[gene].select(startswith('cyc14D'))+.01,
                                          zld.ix[gene].select(startswith('cyc14D'))+.01)
                          )

    for gene in wt.index.intersection(bcd.index):
        wt_bcd.ix[gene] = (DD.earth_mover(wt.ix[gene].select(startswith('cyc14D'))+.01,
                                          bcd.ix[gene].select(startswith('cyc14D_rep1'))+.01)
                           +DD.earth_mover(wt.ix[gene].select(startswith('cyc14D'))+.01,
                                           bcd.ix[gene].select(startswith('cyc14D_rep2'))+.01))/2



both = wt_zld.dropna().index.intersection(wt_bcd.dropna().index)
wt14D_hi = wt.select(startswith('cyc14D'), axis=1).max(axis=1) > 3

print("Plotting")

mpl.clf()
mpl.hist(wt_zld.dropna(), bins=arange(0, 1, .05), normed=True, label='WT vs Zld', histtype='step', color='r')
mpl.hist(wt_bcd.dropna(), bins=arange(0, 1, .05), normed=True, label='WT vs Bcd', histtype='step', color='g')
mpl.legend()
mpl.savefig('analysis/results/WTZldBcdHist.png')

colors = pd.Series('y g b c'.split())
cmaps = pd.Series(index=wt.index, data=0, dtype=int)
for gene in cmaps.index:
    for i, has_bind in enumerate([has_zld, has_bcd]):
            cmaps.ix[gene] += 2**i * (gene in has_bind)

cmaps2 = colors[cmaps]
cmaps2.index = cmaps.index

mpl.clf()
mpl.plot([0,1], [0,1], 'r:')
mpl.scatter(x=wt_zld.ix[both].ix[wt14D_hi],
            y=wt_bcd.ix[both].ix[wt14D_hi],
            c=cmaps2.ix[both].ix[wt14D_hi],
            marker = '.',
            edgecolors='none',
            alpha=0.2)
mpl.xlabel('WT vs Zld'); mpl.ylabel('WT vs Bcd')
mpl.savefig('analysis/results/WTBcdWTZldCorr.png', dpi=600)
