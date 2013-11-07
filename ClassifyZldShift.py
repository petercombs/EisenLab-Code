from __future__ import division
import pandas as pd
from os.path import join
from numpy import mean, all, array, random
from matplotlib import cm
from PlotUtils import svg_heatmap
from glob import glob
from sklearn import svm
from collections import defaultdict


print "Reading in mapping table..."
mapping = pd.read_table('prereqs/gene_map_table_fb_2013_04.tsv', skiprows=5)
fbgn_to_name = {row['primary_FBid'] : row['##current_symbol']
                for i, row in mapping.iterrows()}
name_to_fbgn = {row['##current_symbol']: row['primary_FBid']
                for i, row in mapping.iterrows()}
print "...done"

binding_directory = "prereqs/BDTNP_in_vivo_binding_Release.2.1/Supplemental_Tables"
def get_TF_sites(TF):
    files = glob(join(binding_directory, TF+"_*.txt"))
    sites = defaultdict(float)
    for file in files:
        table = pd.read_table(file)
        for i, row in table.iterrows():
            score = row['Peak_score']
            fbgn = row['Closest_transcribed_gene']
            name = fbgn_to_name.get(fbgn, 'NULL')
            if sites[name] < score:
                sites[name] = score
    return sites

def has_anterior_peak(data):
    n = len(data)
    first_third = data[:n//3]
    middle_third = data[n//3:2*(n//3)]
    if mean(first_third) > 3* mean(middle_third):
        return True
    return False

stage = 'cyc14B'
all_stages = ('cyc11', 'cyc13', 'cyc14A', 'cyc14B')

zld_exp = pd.read_table('analysis/summary.tsv', index_col=0).sort_index().select(lambda x: x.startswith(all_stages), axis=1)
#wt_exp = pd.read_table('prereqs/journal.pone.0071820.s008.txt', index_col=0)
wt_exp = pd.read_table('prereqs/WT5.53_summary.tsv', index_col=0).sort_index().select(lambda x: x.startswith(all_stages), axis=1)

assert all(zld_exp.index == wt_exp.index)


list_of_TFs = ('D TFIIB bcd cad da dl ftz gt h hb hkb kni kr mad med polII '
               'prd run shn slp1 sna tll twi z').split()


categories = pd.Series(0, index=wt_exp.index, dtype=int)
binding = pd.DataFrame(0, index=wt_exp.index, columns = list_of_TFs, dtype=float)

for TF in list_of_TFs:
    sites = get_TF_sites(TF)
    for gene, score in sites.iteritems():
        if gene in binding.index:
            binding.ix[gene][TF] = score

for gene in categories.index:
    wt_gene =  wt_exp.ix[gene].select(lambda x: x.startswith(stage))
    zld_gene =  zld_exp.ix[gene].select(lambda x: x.startswith(stage))
    zld_ant_peak = has_anterior_peak(zld_gene)
    wt_ant_peak = has_anterior_peak(wt_gene)
    if max(wt_gene) < 10 and max(zld_gene) < 10:
        continue
    if zld_ant_peak and not wt_ant_peak:
        categories[gene] = 1
    elif zld_ant_peak and wt_ant_peak:
        categories[gene] = 2
    elif not zld_ant_peak and wt_ant_peak:
        categories[gene] = 3
    elif not zld_ant_peak and not wt_ant_peak:
        categories[gene] = 4

print "Done classifying"

is_exp = binding[categories != 0]
is_exp_cat = categories[categories != 0]

clf = svm.SVC()
test = pd.Series(random.rand(len(is_exp_cat)) > .1, index=is_exp_cat.index)
clf.fit(array(is_exp.ix[test]), array(is_exp_cat[test]))

res = clf.predict(array(is_exp.ix[~test]))

print sum(res == is_exp_cat.ix[~test])
