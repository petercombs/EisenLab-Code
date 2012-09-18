""" Run all the processing code to do all the analyses

Ideally, this will start with the reference files (at the moment, the FlyBase
GFF and FASTA files) and the sequence files in some kind of intelligent
structure, and generates all figures, tables, etc.

Let's see how far we can get...


"""
import subprocess as sp
import os
from glob import glob
import argparse as ap

def make_multi():
    "Make a multi-genome fasta file"
    sp.Popen(['python', 'MakeMultigenome.py', 'Reference/AAA']).wait()
    btb = sp.Popen(['bowtie-build', 'multi.fa', 'multi'])
    bt2b = sp.Popen(['bowtie2-build', 'multi.fa', 'multi'])
    btb.wait()
    bt2b.wait()

def combine_gffs():
    "Combine all GFF files"
    root = os.getcwd()
    os.chdir('Reference/AAA')

    out_name = 'multi.gff'
    out_fh = open(out_name, 'w')
    for fname in glob('d*.gff'):
        print "Getting FlyBase records from: ", fname
        specname = fname.split('-')[0]
        sp.Popen(['awk', '/^[^#].*FlyBase/ {print "%s_" $0}' % specname, fname],
                 stdout=out_fh).wait()

    out_fh.close()
    print "Converting to GTF"
    sp.Popen(['gffread', out_name, '-g', 'multi.fa', '-C', '-F', '-E',
              '-T', '-o', 'multi.gtf']).wait()

    os.chdir(root)

def do_tuxedo_suite():
    "Call the do_tux program to run the tuxedo suite"
    sp.Popen(['python', 'do_tux.py']).wait()

def assign_multireads():
    """Assign reads to species, using appropriate cutoff"""

    commandstr = ['python', 'AssignReads2.py']
    files = glob('analysis/*/accepted_hits.bam')
    sp.Popen(commandstr + files).wait()

def main(args):
    """ Run all the sub-processing code"""

    if not('make_genome' in args.skiplist or '1' in args.skiplist):
        make_multi()

    if not('combine_gff' in args.skiplist or '2' in args.skiplist):
        combine_gffs()

    if not('do_tux' in args.skiplist or '3' in args.skiplist):
        do_tuxedo_suite()


def parse_args():
    """ Parse arguments from the commandline"""
    parser = ap.ArgumentParser()
    parser.add_argument('--skip', '-s', dest='skiplist', action='append')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
