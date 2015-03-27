import pandas as pd
from bisect import bisect
from glob import glob
from collections import defaultdict
import sys
from pickle import load, dump

bind_dist = 20000

ap_early_mat = {'bcd','cad'}
ap_early_gap = {'gt', 'hb', 'kni', 'kr', 'D'}
ap_early_term= {'hkb', 'tll'}

ap_early = ap_early_mat.union(ap_early_gap).union(ap_early_term)

ap_early_zld = ap_early.copy()
ap_early_zld.add('zld')


try:
    tss_dict = load(open('Reference/tss_dict.pkl'))
    tss_per_gene = load(open('Reference/tss_per_gene.pkl'))
except IOError:
    tss_dict = defaultdict(list)
    tss_per_gene = defaultdict(list)
    for i, row in pd.read_table('Reference/tss',
                                keep_default_na=False,
                                na_values=['-', '---', '']).iterrows():
        chr = row['chr'].replace('dmel_', '')
        gene_name = str(row['gene_name'])
        tss_start = row['TSS_start']
        tss_dict[chr].append((tss_start, gene_name))
        tss_per_gene[gene_name].append((chr, tss_start))
    dump(tss_dict, open('Reference/tss_dict.pkl', 'w'))
    dump(tss_per_gene, open('Reference/tss_per_gene.pkl', 'w'))


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
    if len(chrom) == 0:
        return sys.maxsize
    center = bisect(chrom, coord)
    if center < len(chrom):
        dist = abs(chrom[center] - coord)
    else:
        return abs(chrom[-1] - coord)
    if center == 0:
        return dist
    distm1 = abs(coord - chrom[center-1])
    # I should only need one of these, but I'm too lazy to read the
    # documentation to see whether it's the +1 or -1
    return min(dist, distm1)

def get_binding_matrix(genes, dtype=float, dist=bind_dist):
    binding_matrix = pd.DataFrame(index=genes,
                                  columns=tfs,
                                  data=distance_matrix < dist,
                                  dtype=dtype)
    return binding_matrix

tf_files = glob('Reference/*_peaks')
tfs = ([i.split('/')[1].rsplit('_',1)[0] for i in tf_files])

try:
    peak_data = open('Reference/peak_data.pkl')
    tf_sites = load(peak_data)
    distance_matrix = load(peak_data)
    peak_data.close()

except IOError:
    peak_data = open('Reference/peak_data.pkl', 'w')
    #has_tfs = {}
    tf_sites = {}
    for tf, tf_file in zip(tfs, tf_files):
        #has_tf = set()
        tf_chrs = defaultdict(list)
        for coord in pd.read_table(tf_file).NewPeak:
            chr, coord = coord.split(':')
            tf_chrs[chr].append(int(coord))
            #has_tf.update(find_near(tss_dict[chr], int(coord), bind_dist))
        #has_tfs[tf] = has_tf
        for chr in tf_chrs.itervalues():
            chr.sort()
        tf_sites[tf] = tf_chrs


    dump(tf_sites, peak_data)
    distance_matrix = pd.DataFrame(index=tss_per_gene.keys(),
                                   columns=tfs,
                                   dtype=int)

    for gene in distance_matrix.index:
        for tf in distance_matrix.columns:
            dists = []
            for chr, coord in tss_per_gene[gene]:
                dists.append(find_nearest_dist(tf_sites[tf][chr], coord))
            distance_matrix.ix[gene, tf] = min(dists)
    dump(distance_matrix, peak_data)
    peak_data.close()


