#!/usr/bin/env python

from __future__ import division, print_function
import pandas as pd
from os.path import join
from numpy import mean
from matplotlib import cm
from PlotUtils import svg_heatmap
from PeakFinder import has_anterior_peak, has_posterior_peak, has_central_peak
from glob import glob


mapping = pd.read_table('prereqs/gene_map_table_fb_2013_04.tsv', skiprows=5)

print("Reading in mapping table...")
fbgn_to_name = {row['primary_FBid'] : row['##current_symbol']
                for i, row in mapping.iterrows()}
name_to_fbgn = {row['##current_symbol']: row['primary_FBid']
                for i, row in mapping.iterrows()}
print("...done")

binding_directory = "prereqs/BDTNP_in_vivo_binding_Release.2.1/Supplemental_Tables"
def get_TF_sites(TF):
    files = glob(join(binding_directory, TF+"_*.txt"))
    sites = set()
    for file in files:
        table = pd.read_table(file)
        sites.update({fbgn_to_name[fbgn]
                      for fbgn in table['Closest_transcribed_gene']
                      if fbgn in fbgn_to_name})
    return sites


#bcd1_fname = join(binding_directory, 'bcd_1_012505-sym-1_table.txt')
#bcd2_fname = join(binding_directory, 'bcd_2_092005-sym-1_table.txt')

#bcd1 = pd.read_table(bcd1_fname)
#bcd2 = pd.read_table(bcd2_fname)

#bcd_genes = {fbgn_to_name[gene]
             #for gene in bcd1['Closest_transcribed_gene']
             #if gene in fbgn_to_name}
##bcd_genes.update({fbgn_to_name[gene]
                  #for gene in bcd2['Closest_transcribed_gene']
                  #if gene in fbgn_to_name})


stagewt = 'cyc13'
stagezld = 'cyc13'
all_stages = ('cyc11', 'cyc13', 'cyc14A', 'cyc14B')

zld_exp = pd.read_table('analysis/summary.tsv', index_col=0).sort_index().select(lambda x: x.startswith(all_stages), axis=1)
#wt_exp = pd.read_table('prereqs/journal.pone.0071820.s008.txt', index_col=0)
wt_exp = pd.read_table('prereqs/WT5.53_summary.tsv', index_col=0).sort_index().select(lambda x: x.startswith(all_stages), axis=1)
wt_exp = wt_exp.drop(['cyc14B_sl06_FPKM'], axis=1)
in_both = set(wt_exp.index).intersection(set(zld_exp.index))

antant = set()
antshift = set()
antlost = set()
noant = set()

postpost = set()
postshift = set()
postlost = set()
nopost = set()

centcent = set()
centshift = set()
centlost = set()
nocent = set()

print("Assigning to categories")

for gene in in_both:
    wt_gene =  wt_exp.ix[gene].select(lambda x: x.startswith(stagewt))
    zld_gene =  zld_exp.ix[gene].select(lambda x: x.startswith(stagezld))
    if max(wt_gene) < 10 and max(zld_gene) < 10:
        continue
    zld_has_ant = has_anterior_peak(zld_gene)
    wt_has_ant = has_anterior_peak(wt_gene)
    ant_types = [[noant, antlost], [antshift, antant]]
    ant_types[zld_has_ant][wt_has_ant].add(gene)

    zld_has_post = has_posterior_peak(zld_gene)
    wt_has_post = has_posterior_peak(wt_gene)
    post_types = [[nopost, postlost], [postshift, postpost]]
    post_types[zld_has_post][wt_has_post].add(gene)

    zld_has_cent = has_central_peak(zld_gene)
    wt_has_cent = has_central_peak(wt_gene)
    cent_types = [[nocent, centlost], [centshift, centcent]]
    cent_types[zld_has_cent][wt_has_cent].add(gene)

print("...done")

list_of_TFs = ('D TFIIB bcd cad da dl ftz gt h hb hkb kni kr mad med polII '
               'prd run shn slp1 sna tll twi z').split()

for TF in list_of_TFs:
    n_antant = 0
    n_antshift = 0
    n_antlost = 0
    n_noant = 0

    TF_sites = get_TF_sites(TF)
    for gene in antant:
        if gene in TF_sites:
            n_antant+= 1
    for gene in antshift:
        if gene in TF_sites:
            n_antshift+= 1
    for gene in antlost:
        if gene in TF_sites:
            n_antlost+= 1
    for gene in noant:
        if gene in TF_sites:
            n_noant+= 1



    n_postpost = 0
    n_postshift = 0
    n_postlost = 0
    n_nopost = 0

    for gene in postpost:
        if gene in TF_sites:
            n_postpost+= 1
    for gene in postshift:
        if gene in TF_sites:
            n_postshift+= 1
    for gene in postlost:
        if gene in TF_sites:
            n_postlost+= 1
    for gene in nopost:
        if gene in TF_sites:
            n_nopost+= 1

    print('-'*32, '\n', TF,'\n', '-'*32)
    print("The {} anterior shifted genes have {} sites ({:%})".format(len(antshift),
                                                                      n_antshift,
                                                                      n_antshift/len(antshift))
         )

    print("The {} anterior maintained genes have {} sites ({:%})".format(len(antant),
                                                                         n_antant,
                                                                         n_antant/len(antant))
         )

    print("The {} anterior lost genes have {} sites ({:%})".format(len(antlost),
                                                                      n_antlost,
                                                                      n_antlost/len(antlost))
         )

    print("The {} never anterior genes have {} sites ({:%})".format(len(noant),
                                                                      n_noant,
                                                                      n_noant/len(noant))
         )

    print('-'*15)
    print("The {} posterior shifted genes have {} sites ({:%})".format(len(postshift),
                                                                      n_postshift,
                                                                      n_postshift/len(postshift))
         )

    print("The {} posterior maintained genes have {} sites ({:%})".format(len(postpost),
                                                                         n_postpost,
                                                                         n_postpost/(len(postpost)+1e-6))
         )

    print("The {} posterior lost genes have {} sites ({:%})".format(len(postlost),
                                                                      n_postlost,
                                                                      n_postlost/(len(postlost)+1e-6))
         )

    print("The {} never posterior genes have {} sites ({:%})".format(len(nopost),
                                                                      n_nopost,
                                                                      n_nopost/(len(nopost)
                                                                               +
                                                                               1e-6))
         )


norm_col = wt_exp.select(lambda x: not x.startswith('cyc11'), axis=1).max(axis=1)
norm_col[norm_col<10] = 10
svg_heatmap((wt_exp.select(antshift.__contains__),
             zld_exp.select(antshift.__contains__)),
            'analysis/results/antshift.svg',
            norm_rows_by=norm_col.select(antshift.__contains__),
            total_width=300,
            box_height=5,
            draw_row_labels=True,
            cmap=(cm.Blues, cm.Reds),
            col_sep='_sl',
           )

svg_heatmap((wt_exp.select(antlost.__contains__),
             zld_exp.select(antlost.__contains__)),
            'analysis/results/antlost.svg',
            norm_rows_by=norm_col.select(antlost.__contains__),
            total_width=300,
            box_height=5,
            draw_row_labels=True,
            cmap=(cm.Blues, cm.Reds),
            col_sep='_sl',
           )

svg_heatmap((wt_exp.select(antant.__contains__),
             zld_exp.select(antant.__contains__)),
            'analysis/results/antant.svg',
            norm_rows_by=norm_col.select(antant.__contains__),
            total_width=300,
            box_height=5,
            draw_row_labels=True,
            cmap=(cm.Blues, cm.Reds),
            col_sep='_sl',
           )


svg_heatmap((wt_exp.select(postshift.__contains__),
             zld_exp.select(postshift.__contains__)),
            'analysis/results/postshift.svg',
            norm_rows_by=norm_col.select(postshift.__contains__),
            total_width=300,
            box_height=5,
            draw_row_labels=True,
            cmap=(cm.Blues, cm.Reds),
            col_sep='_sl',
           )

svg_heatmap((wt_exp.select(postlost.__contains__),
             zld_exp.select(postlost.__contains__)),
            'analysis/results/postlost.svg',
            norm_rows_by=norm_col.select(postlost.__contains__),
            total_width=300,
            box_height=5,
            draw_row_labels=True,
            cmap=(cm.Blues, cm.Reds),
            col_sep='_sl',
           )

svg_heatmap((wt_exp.select(postpost.__contains__),
             zld_exp.select(postpost.__contains__)),
            'analysis/results/postpost.svg',
            norm_rows_by=norm_col.select(postpost.__contains__),
            total_width=300,
            box_height=5,
            draw_row_labels=True,
            cmap=(cm.Blues, cm.Reds),
            col_sep='_sl',
           )


svg_heatmap((wt_exp.select(centshift.__contains__),
             zld_exp.select(centshift.__contains__)),
            'analysis/results/centshift.svg',
            norm_rows_by=norm_col.select(centshift.__contains__),
            total_width=300,
            box_height=5,
            draw_row_labels=True,
            cmap=(cm.Blues, cm.Reds),
            col_sep='_sl',
           )

svg_heatmap((wt_exp.select(centlost.__contains__),
             zld_exp.select(centlost.__contains__)),
            'analysis/results/centlost.svg',
            norm_rows_by=norm_col.select(centlost.__contains__),
            total_width=300,
            box_height=5,
            draw_row_labels=True,
            cmap=(cm.Blues, cm.Reds),
            col_sep='_sl',
           )

svg_heatmap((wt_exp.select(centcent.__contains__),
             zld_exp.select(centcent.__contains__)),
            'analysis/results/centcent.svg',
            norm_rows_by=norm_col.select(centcent.__contains__),
            total_width=300,
            box_height=5,
            draw_row_labels=True,
            cmap=(cm.Blues, cm.Reds),
            col_sep='_sl',
           )

#xlf = pd.ExcelFile('prereqs/ZLD_vs_expression-1203.xlsx')
#zld_bind = xlf.parse(xlf.sheet_names[0])
#zld_sorted = zld_bind.sort('ZLD Prom st3')
#zld_sites = zld_sorted.Gene.dropna().unique()
#
#svg_heatmap((wt_exp.ix[zld_sites][:500], zld_exp.ix[zld_sites][:500]),
#            'analysis/results/sort_by_zld_st3.svg',
#            norm_rows_by=norm_col.ix[zld_sites][:500],
#            total_width=300,
#            box_height=5,
#            draw_row_labels=True,
#            cmap=(cm.Blues, cm.Reds),
#            col_sep='_sl',
#           )
#
