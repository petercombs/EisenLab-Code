import pandas as pd
from glob import glob
from os import path

whole_tab = None


for fname in glob('*/genes.fpkm_tracking'):
        tab = pd.read_table(fname)
        tab = tab.rename(columns=lambda x: path.dirname(fname) + "_FPKM"
                         if x == "FPKM" else (x.replace("FPKM", 
                                                        path.dirname(fname)) 
                                              if x.startswith("FPKM")
                                              else x))
        if whole_tab is not None:
            whole_tab = pd.merge(whole_tab, tab)
        else:
            whole_tab = tab

whole_tab.to_csv('merged_genes.fpkm_tracking', sep='\t')
