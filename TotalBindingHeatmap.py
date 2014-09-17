from __future__ import print_function, division
import pandas as pd
from numpy import isfinite, all
from numpy.linalg import LinAlgError
from glob import glob
from bisect import bisect
from PlotUtils import svg_heatmap
from matplotlib import pyplot as mpl
from os.path import basename
from collections import defaultdict
from statsmodels import api as sm

read_table_args = dict(keep_default_na=False, na_values=['---', ''], index_col=0)

eps = .1
cyc = 'cyc14D'
sw_cyc = lambda x: x.startswith(cyc)
bind_dist = 5000

try:
    wt = locals()['wt']
    assert min(wt.min()) == eps
except (KeyError, AssertionError):
    print("Importing data")
    keep_old = False
    wt = pd.read_table('prereqs/WT6.01.summary.tsv', **read_table_args)
    zld = pd.read_table('prereqs/Zld6.01.summary_merged.tsv', **read_table_args)
    bcd = pd.read_table('analysis/summary.tsv', **read_table_args)
    for dataset in (wt, zld, bcd):
        dataset.ix[:,:] += eps
        for i, column in enumerate(dataset.columns):
            if not all(isfinite(dataset.ix[:, column])):
                try:
                    dataset.ix[:, column] = (dataset.ix[:, i-1]
                                             + dataset.ix[:, i+1])/2
                except:
                    print("Whoops! {}".format(column))
all_fcs = (wt.select(sw_cyc, axis=1).mean(axis=1)/
           (bcd.select(sw_cyc, axis=1).mean(axis=1)+eps))
try:
    tss_dict = locals()['tss_dict']
except KeyError:
    tss_dict = defaultdict(list)
    tss_per_gene = defaultdict(list)
    for i, row in pd.read_table('Reference/tss').iterrows():
        tss_dict[row['chr'].replace('dmel_', '')].append((row['TSS_start'],
                                                          row['gene_name']))
        tss_per_gene[row['gene_name']].append((row['chr'].replace('dmel_', ''),
                                              row['TSS_start']))

def find_near(chrom, coord, dist):
    out = set()
    center = bisect([i[0] for i in chrom], coord)
    for (coord2, gene) in reversed(chrom[:center]):
        if coord - coord2 > dist: break
        out.add(gene)

    for (coord2, gene) in chrom[center:]:
        if coord2 - coord > dist: break
        out.add(gene)

    return out

def find_nearest_dist(chrom, coord):
    center = bisect([i[0] for i in chrom, coord])
    dist = abs(chrom[center][0] - coord)
    distp1 = abs(coord - chrom[center+1][0])
    distm1 = abs(coord - chrom[center-1][0])
    # I should only need one of these, but I'm too lazy to read the
    # documentation to see whether it's the +1 or -1
    return min(dist, distp1, distm1)

def get_binding_matrix(genes):
    binding_matrix = pd.DataFrame(index=genes,
                                  columns=tfs, data=0.0, dtype=float)
    for tf in tfs:
        for gene in binding_matrix.index:
            binding_matrix.ix[gene, tf] = float(gene in has_tfs[tf])
    binding_matrix.sort_index(axis=1, inplace=True)
    return binding_matrix


tf_files = glob('Reference/*_peaks')
tfs = ([i.split('/')[1].split('_')[0] for i in tf_files])
has_tfs = {}
for tf, tf_file in zip(tfs, tf_files):
    has_tf = set()
    for coord in pd.read_table(tf_file).NewPeak:
        chr, coord = coord.split(':')
        has_tf.update(find_near(tss_dict[chr], int(coord), bind_dist))
    has_tfs[tf] = has_tf


bms = []
fits = []

for fname in ['analysis/results/change_bcd.tsv',
              'analysis/results/change_both.tsv',]:
    genes = pd.read_table(fname, **read_table_args).index
    fcs = pd.Series(index=genes, data=0)
    for gene in fcs.index:
        fcs.ix[gene] = all_fcs.ix[gene]
    fcs.sort()
    #mpl.pcolormesh(binding_matrix)
    #mpl.savefig('analysis/results/{}_bindingmap.png'.format(basename(fname)),
                #dpi=300)
    binding_matrix = get_binding_matrix(fcs.index)
    svg_heatmap(binding_matrix,
                'analysis/results/{}_bindingmap.svg'.format(basename(fname)),
                draw_row_labels=True, box_size=10, box_height=10)
    bms.append(binding_matrix)


    binding_matrix['const'] = 1.0

    try:
        fit = sm.Logit(fcs>1, binding_matrix).fit()
        fits.append(fit)
        print(fit.summary())
    except LinAlgError:
        print("No luck on "+fname)

all_dist = pd.read_table('analysis/results/change_both.tsv', **read_table_args)
binding_matrix = get_binding_matrix(all_dist.index)
svg_heatmap(binding_matrix,
            'analysis/results/change_both_sorted.svg',
            draw_row_labels=True, box_size=10, box_height=10)
bms.append(binding_matrix)



cats = pd.read_table('analysis/results/change_cats.txt', header=None,
                     names=['genes', 'changes'],
                     **read_table_args).changes.dropna()
binding_matrix = get_binding_matrix(cats.index)

binding_matrix['const'] = 1.0

fitMN = sm.MNLogit(cats, binding_matrix).fit()
fits.append(fitMN)
print(fitMN.summary())

