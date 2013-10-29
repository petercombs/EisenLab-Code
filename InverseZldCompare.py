from __future__ import division
import pandas as pd
from os.path import join
from numpy import mean
import numpy as np
from matplotlib import cm
from PlotUtils import svg_heatmap
from glob import glob
import sys

binding_directory = "prereqs/BDTNP_in_vivo_binding_Release.2.1/Supplemental_Tables"

mapping = pd.read_table('prereqs/gene_map_table_fb_2013_04.tsv', skiprows=5)

print "Reading in mapping table..."
fbgn_to_name = {row['primary_FBid'] : row['##current_symbol']
                for i, row in mapping.iterrows()}
name_to_fbgn = {row['##current_symbol']: row['primary_FBid']
                for i, row in mapping.iterrows()}
print "...done"



def get_TF_sites(TF):
    files = glob(join(binding_directory, TF+"_*.txt"))
    sites = set()
    for file in files:
        table = pd.read_table(file)
        sites.update({fbgn_to_name[fbgn]
                      for fbgn in table['Closest_transcribed_gene']
                      if fbgn in fbgn_to_name})
    return sites


def has_anterior_peak(data):
    n = np.shape(data)[-1]
    first_third = mean(data.ix[:,:n//3], axis=1)
    middle_third = mean(data.ix[:,n//3:-(n//3)], axis=1)
    return (first_third > 3* middle_third) * first_third > 3

def has_posterior_peak(data):
    n = np.shape(data)[-1]
    first_third = mean(data.ix[:,-n//3:], axis=1)
    middle_third = mean(data.ix[:,n//3:-(n//3)], axis=1)
    return (first_third > 3* middle_third) * (first_third > 3)

def has_central_peak(data):
    n = np.shape(data)[-1]
    first_third = mean(data.ix[:,:n//3], axis=1)
    middle_third = mean(data.ix[:,n//3:-(n//3)], axis=1)
    last_third = mean(data.ix[:,-(n//3):], axis=1)
    return (middle_third > 3*first_third)*(middle_third >  last_third) * (middle_third > 3)


def make_sort_num(wt_data, zld_data):
    return (1 * has_posterior_peak(zld_data)
            + 2 * has_posterior_peak(wt_data)
            + 4 * has_central_peak(zld_data)
            + 8 * has_central_peak(wt_data)
            + 16 * has_anterior_peak(zld_data)
            + 32 * has_anterior_peak(wt_data)
           )

stagewt = sys.argv[1] if (len(sys.argv) > 1) else 'cyc13'
stagezld = sys.argv[1] if (len(sys.argv) > 1) else 'cyc13'
all_stages = ('cyc11', 'cyc13', 'cyc14A', 'cyc14B')

zld_exp = pd.read_table('analysis/summary_merged.tsv', index_col=0).sort_index().select(lambda x: x.startswith(all_stages), axis=1)
#wt_exp = pd.read_table('prereqs/journal.pone.0071820.s008.txt', index_col=0)
wt_exp = pd.read_table('prereqs/WT5.53_summary.tsv', index_col=0).sort_index().select(lambda x: x.startswith(all_stages), axis=1)
wt_exp = wt_exp.drop(['cyc14B_sl06_FPKM'], axis=1)

wt_comp = wt_exp.select(lambda x: x.startswith(stagewt), axis=1)
zld_comp = zld_exp.select(lambda x: x.startswith(stagezld), axis=1)

assert all(wt_exp.index == zld_exp.index)

TFs = ('bcd cad D da dl ftz gt h hb hkb kni kr mad med polII prd run shn '
           'slp1 sna TFIIB tll twi z'.split())
zld_ant = 16
wt_ant = 32

print "Comparing WT {} to zld- {}".format(stagewt, stagezld)
print "TF\tpval     \tSig \t[[gain,nogain],[gain_TF, nogain_TF]]\t Odds Ratio"
for TF in TFs:
    TF_sites = get_TF_sites(TF)
    TF_sites.intersection_update(wt_exp.index)
    sort_num = make_sort_num(wt_comp, zld_comp)

    TF_genes = sort_num.ix[TF_sites]

    assert sort_num.dtype == TF_genes.dtype

    TF_genes.sort()
    TF_genes = TF_genes[::-1]

    antgain = sum(((sort_num & zld_ant) * ~(sort_num & wt_ant)) != 0)
    antgain_TF = sum(((TF_genes & zld_ant) * ~(TF_genes & wt_ant)) != 0)
    has_site = sum((wt_comp.min(axis=1)>3).ix[TF_sites])
    is_exp = sum(wt_comp.min(axis=1) > 3)

    from scipy import stats

    fisher_tab = [[antgain, is_exp-antgain], [antgain_TF, has_site-antgain_TF]]
    odds_ratio, pval =  stats.fisher_exact(fisher_tab,
                                           alternative='two-sided')
    print TF, "\t%9g"%pval, '\t',
    print '*' * int(np.ceil(max((np.log10(.05/len(TFs)) -
                                                    np.log10(pval)), 0))),
    if pval < .05/len(TFs):
        print '\t', fisher_tab, '\t', 1/odds_ratio
    else:
        print


