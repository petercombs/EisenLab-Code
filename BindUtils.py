import pandas as pd
from bisect import bisect
from glob import glob
from collections import defaultdict

bind_dist = 5000

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

def find_near(chrom, coord, dist, nearest_only=False):
    out = set()
    center = bisect([i[0] for i in chrom], coord)
    for (coord2, gene) in reversed(chrom[:center]):
        if coord - coord2 > dist: break
        out.add(gene)
        if nearest_only:
            break

    for (coord2, gene) in chrom[center:]:
        if coord2 - coord > dist: break
        out.add(gene)
        if nearest_only:
            break

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
