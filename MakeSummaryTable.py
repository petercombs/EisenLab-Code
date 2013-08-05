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

def get_stagenum(name, series, dir):
    # Slice until the first digit
    name_base = name[:[i for i, c in enumerate(name) if c.isdigit()][0]]
    dir = {'+':1, '-':-1}[dir]
    return (sorted(ix for ix in series if ix.startswith(name_base))[::dir]
            .index(name)) + 1


fnames = glob(path.join(argv[1], '*', 'genes.fpkm_tracking'))
if len(argv) > 2 and argv[2] != '-c':
    has_params = argv[2]
    params = pandas.read_table(has_params, index_col='Label')
else:
    has_params = False

conf = '-c' in argv

df = None
for fname in sorted(fnames):
    table = pandas.read_table(fname)
    alldir, fname = path.split(fname)
    basedir, dirname = path.split(alldir)
    table = table.drop_duplicates('gene_short_name').dropna(how='any')
    table.set_index('gene_short_name', inplace=True, verify_integrity=True)

    if has_params:
        new_dirname = "cyc{stage}_{num:02}".format(
            stage=params.ix[dirname]['Stage'],
            num=get_stagenum(dirname, params.index,
                             params.ix[dirname]['Direction']))
        print dirname, '=', new_dirname
        dirname = new_dirname

    if df is None:
        df = pandas.DataFrame({dirname+"_FPKM": table.FPKM})
    else:
        df.insert(len(df.columns), dirname+"_FPKM", table.FPKM)

    if conf:
        df.insert(len(df.columns),
                  dirname+"_conf_lo",
                  table.FPKM_conf_lo)
        df.insert(len(df.columns),
                  dirname+"_conf_hi",
                  table.FPKM_conf_hi)


df.sort_index(axis=1).to_csv(path.join(argv[1], 'summary' + ('_with_conf' * conf) + '.tsv'),
                             sep='\t')

