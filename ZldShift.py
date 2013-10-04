#!/usr/bin/env python

from __future__ import division
import pandas as pd
from os.path import join
from numpy import mean
from matplotlib import cm
from PlotUtils import svg_heatmap


mapping = pd.read_table('prereqs/gene_map_table_fb_2013_04.tsv', skiprows=5)

fbgn_to_name = {row['primary_FBid'] : row['##current_symbol']
                for i, row in mapping.iterrows()}
name_to_fbgn = {row['##current_symbol']: row['primary_FBid']
                for i, row in mapping.iterrows()}

binding_directory = "prereqs/BDTNP_in_vivo_binding_Release.2.1/Supplemental_Tables"
bcd1_fname = join(binding_directory, 'bcd_1_012505-sym-1_table.txt')
bcd2_fname = join(binding_directory, 'bcd_2_092005-sym-1_table.txt')

bcd1 = pd.read_table(bcd1_fname)
bcd2 = pd.read_table(bcd2_fname)

bcd_genes = {fbgn_to_name[gene]
             for gene in bcd1['Closest_transcribed_gene']
             if gene in fbgn_to_name}
bcd_genes.update({fbgn_to_name[gene]
                  for gene in bcd2['Closest_transcribed_gene']
                  if gene in fbgn_to_name})


def has_anterior_peak(data):
    n = len(data)
    first_third = data[:n//3]
    middle_third = data[n//3:2*(n//3)]
    if mean(first_third) > 3* mean(middle_third):
        return True
    return False

stage = 'cyc14B'
all_stages = ('cyc11', 'cyc14A', 'cyc14B')

zld_exp = pd.read_table('analysis/summary.tsv', index_col=0).sort_index().select(lambda x: x.startswith(all_stages), axis=1)
#wt_exp = pd.read_table('prereqs/journal.pone.0071820.s008.txt', index_col=0)
wt_exp = pd.read_table('prereqs/WT5.53_summary.tsv', index_col=0).sort_index().select(lambda x: x.startswith(all_stages), axis=1)
in_both = set(wt_exp.index).intersection(set(zld_exp.index))

antant = set()
n_antant = 0
antshift = set()
n_antshift = 0
antlost = set()
n_antlost = 0
noant = set()
n_noant = 0

for gene in in_both:
    wt_gene =  wt_exp.ix[gene].select(lambda x: x.startswith(stage))
    zld_gene =  zld_exp.ix[gene].select(lambda x: x.startswith(stage))
    if max(wt_gene) < 10 and max(zld_gene) < 10:
        continue
    if has_anterior_peak(zld_gene) and not has_anterior_peak(wt_gene):
        antshift.add(gene)
        n_antshift += gene in bcd_genes
    elif has_anterior_peak(zld_gene) and has_anterior_peak(wt_gene):
        antant.add(gene)
        n_antant += gene in bcd_genes
    elif not has_anterior_peak(zld_gene) and has_anterior_peak(wt_gene):
        antlost.add(gene)
        n_antlost += gene in bcd_genes
    elif not has_anterior_peak(zld_gene) and not has_anterior_peak(wt_gene):
        noant.add(gene)
        n_noant += gene in bcd_genes

print "The {} anterior shifted genes have {} sites ({:%})".format(len(antshift),
                                                                  n_antshift,
                                                                  n_antshift/len(antshift))

print "The {} anterior maintained genes have {} sites ({:%})".format(len(antant),
                                                                     n_antant,
                                                                     n_antant/len(antant))

print "The {} anterior lost genes have {} sites ({:%})".format(len(antlost),
                                                                  n_antlost,
                                                                  n_antlost/len(antlost))

print "The {} never anterior genes have {} sites ({:%})".format(len(noant),
                                                                  n_noant,
                                                                  n_noant/len(noant))

norm_col = wt_exp.max(axis=1)
norm_col[norm_col<10] = 10
svg_heatmap((wt_exp.select(antshift.__contains__),
             zld_exp.select(antshift.__contains__)),
            'analysis/results/antshift.svg',
            norm_rows_by=norm_col.select(antshift.__contains__),
            total_width=300,
            box_height=5,
            draw_row_labels=True,
            cmap=(cm.Blues, cm.Blues),
            col_sep='_sl',
           )

svg_heatmap((wt_exp.select(antlost.__contains__),
             zld_exp.select(antlost.__contains__)),
            'analysis/results/antlost.svg',
            norm_rows_by=norm_col.select(antlost.__contains__),
            total_width=300,
            box_height=5,
            draw_row_labels=True,
            cmap=(cm.Blues, cm.Blues),
            col_sep='_sl',
           )

svg_heatmap((wt_exp.select(antant.__contains__),
             zld_exp.select(antant.__contains__)),
            'analysis/results/antant.svg',
            norm_rows_by=norm_col.select(antant.__contains__),
            total_width=300,
            box_height=5,
            draw_row_labels=True,
            cmap=(cm.Blues, cm.Blues),
            col_sep='_sl',
           )


xlf = pd.ExcelFile('prereqs/ZLD_vs_expression-1203.xlsx')
zld_bind = xlf.parse(xlf.sheet_names[0])
zld_sorted = zld_bind.sort('ZLD Prom st3')
zld_sites = zld_sorted.Gene.dropna().unique()

svg_heatmap((wt_exp.ix[zld_sites][:500], zld_exp.ix[zld_sites][:500]),
            'analysis/results/sort_by_zld_st3.svg',
            norm_rows_by=norm_col.ix[zld_sites][:500],
            total_width=300,
            box_height=5,
            draw_row_labels=True,
            cmap=(cm.Blues, cm.Reds),
            col_sep='_sl',
           )

