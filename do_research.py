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
    sp.Popen(['python', 'MakeMultigenome.py', 'Reference/AAA', 'mel', 'per',
              'wil', 'moj', 'vir']).wait()
    #btb = sp.Popen(['bowtie-build', 'multi.fa', 'multi'])
    #bt2b = sp.Popen(['bowtie2-build', 'multi.fa', 'multi'])
    btb.wait()
    bt2b.wait()

def make_gtfs():
    "Make GTF files from all flybase GFF files"
    root = os.getcwd()
    os.chdir('Reference/AAA')

    for fname in glob('d*.gff'):
        print "Converting %s to GTF" % fname
        sp.Popen(['gffread', out_name, '-E', '-T', '-o', 
                  fname.replace('gff', 'gtf')]).wait()

    os.chdir(root)

def map_reads():
    "Call the do_tux program to run the tuxedo suite"
    sp.Popen(['python', 'MapReads.py']).wait()

def assign_multireads():
    """Assign reads to species, using appropriate cutoff"""

    commandstr = ['python', 'AssignReads2.py']
    files = glob('analysis/*/accepted_hits.bam')
    sp.Popen(commandstr + files).wait()

def main(args):
    """ Run all the sub-processing code"""

    if not('make_gtfs' in args.skiplist or '1' in args.skiplist):
        make_gtfs()

    if not('make_genome' in args.skiplist or '2' in args.skiplist):
        make_multi()

    if not('map_reads' in args.skiplist or '3' in args.skiplist):
        map_reads()


def parse_args():
    """ Parse arguments from the commandline"""
    parser = ap.ArgumentParser()
    parser.add_argument('--skip', '-s', dest='skiplist', action='append')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
