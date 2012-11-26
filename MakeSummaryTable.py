"""MakeSummaryTable

Takes a collection of fpkm_tracking files and makes a summary table with all of
the data as a CSV file.  Arguments:
    1) superdirectory
    2) -c -- also include confidence intervals

"""
import pandas
import os
from os import path
from glob import glob
from sys import argv


fnames = glob(path.join(argv[1], '*', 'genes.fpkm_tracking'))

conf = '-c' in argv

df = None
for fname in sorted(fnames):
    table = pandas.read_table(fname)
    alldir, fname = path.split(fname)
    basedir, dirname = path.split(alldir)
    table = table.drop_duplicates('gene_short_name').dropna(how='any')
    table.set_index('gene_short_name', inplace=True, verify_integrity=True)
    if df is None:
        df = pandas.DataFrame({dirname+"_FPKM": table.FPKM})
        if conf:
            #df.insert(len(df.columns),
                      #dirname+"_conf_range", 
                      #(table.FPKM_conf_hi - table.FPKM_conf_lo))
            df.insert(len(df.columns),
                      dirname+"_conf_lo",
                      table.FPKM_conf_lo)
            df.insert(len(df.columns),
                      dirname+"_conf_hi",
                      table.FPKM_conf_hi)
    else:
        df.insert(len(df.columns), dirname+"_FPKM", table.FPKM)
        if conf:
            #df.insert(len(df.columns), dirname + "_conf_range",
                      #table.FPKM_conf_hi - table.FPKM_conf_lo)
            df.insert(len(df.columns),
                      dirname+"_conf_lo",
                      table.FPKM_conf_lo)
            df.insert(len(df.columns),
                      dirname+"_conf_hi",
                      table.FPKM_conf_hi)
    

df.to_csv(path.join(argv[1], 'summary' + ('_with_conf' * conf) + '.tsv'),
          sep='\t')

