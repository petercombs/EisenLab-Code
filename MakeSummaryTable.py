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

def parse_args():
    from argparse import ArgumentParser

    parser = ArgumentParser(description=
                           "Take a collection of fpkm_tracking files and makes "
                            "a summary table with all of the data as a TSV "
                            "file.")
    parser.add_argument('--confidence', '-c', default=False,
                        dest='conf', action='store_true',
                        help="Include confidence intervals") 
    parser.add_argument('--params', '-p', default=False, 
                        dest='has_params',
                        help="Parameters file including renaming conventions: "
                        'Directory in column "Label", stage in column "Stage"')
    parser.add_argument('--key', '-k', default='gene_short_name', 
                        help='The column to combine on (FBgn in tracking_id)')
    parser.add_argument('--strip-low-reads', '-s', default=0, type=int,
                        help='Remove slices with fewer than N counts (off by'
                        'default')
    parser.add_argument('basedir', 
                        help='The directory containing directories, which '
                        'contain genes.fpkm_tracking files')
    return parser.parse_args()


def get_stagenum(name, series, dir):
    # Slice until the first digit
    name_base = name[:[i for i, c in enumerate(name) if c.isdigit()][0]]
    dir = {'+':1, '-':-1}[dir]
    return (sorted(ix for ix in series if ix.startswith(name_base))[::dir]
            .index(name)) + 1


args = parse_args()
fnames = glob(path.join(args.basedir, '*', 'genes.fpkm_tracking'))
if args.has_params:
    has_params = argv[2]
    params = pandas.read_table(has_params, index_col='Label')
    params = params.dropna(how='any')


df = None
for fname in sorted(fnames):
    table = pandas.read_table(fname)
    alldir, fname = path.split(fname)
    basedir, dirname = path.split(alldir)
    table = table.drop_duplicates(args.key).dropna(how='any')
    table.set_index(args.key, inplace=True, verify_integrity=True)
    if args.has_params and dirname not in params.index:
        continue

    if args.has_params:
        new_dirname = "cyc{stage}_sl{num:02}".format(
            stage=params.ix[dirname]['Stage'],
            num=get_stagenum(dirname, params.index,
                             params.ix[dirname]['Direction']))
        print dirname, '=', new_dirname
        dirname = new_dirname

    if args.strip_low_reads:
        from pysam import Samfile
        sf = Samfile(path.join(alldir, 'assigned_dmelR.bam'))
        if sf.mapped < args.strip_low_reads:
            print "Skipping", dirname 
            continue
    if df is None:
        df = pandas.DataFrame({dirname+"_FPKM": table.FPKM})
    else:
        df.insert(len(df.columns), dirname+"_FPKM", table.FPKM)

    if args.conf:
        df.insert(len(df.columns),
                  dirname+"_conf_lo",
                  table.FPKM_conf_lo)
        df.insert(len(df.columns),
                  dirname+"_conf_hi",
                  table.FPKM_conf_hi)


df.sort_index(axis=1).to_csv(path.join(args.basedir, 
                                       'summary' + ('_with_conf' * args.conf) + '.tsv'),
                             sep='\t')

