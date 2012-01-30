#from sys import stdout
import re
from os.path import splitext
from argparse import ArgumentParser
from Bio import SeqIO

from PointClouds import to_number

CHROM = 0
START = 1
STOP = 2
STRAND = 3

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('-u', '--upstream', type=int, default=200,
                        help="Number of upstream bases to include (default 200)")
    parser.add_argument('-d', '--downstream', type=int, default=0,
                        help="Number of downstream bases to include (default 0)")
    parser.add_argument('sequence', type=file,
                        help='The FASTA-formatted sequence file')
    parser.add_argument('targets', type=file, nargs="+",
                        help='List of TSSs of interest')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    chroms = {chrom.id:chrom for chrom in SeqIO.parse(args.sequence, 'fasta')}

    
    for target in args.targets:
        #print target.name
        #outfile = open(splitext(target.name)[0] + '.fa', 'w')
        for line in target:
            if not line.startswith('#'):
                continue
            comm, dir, chr, lo_pos, hi_pos, strand, samps = line.split()
            samps = samps[1:-1].replace(',', '+')
            outfile = open(dir+samps+'.fasta', 'a')
            if strand == '+':
                start = int(lo_pos)
                promoter = chroms[chr][start - args.upstream
                                              : start + args.downstream]
            else:
                start = int(hi_pos)
                promoter = chroms[chr][start - args.downstream
                                              : start + args.upstream]
            promoter.id = '_'.join([chr, lo_pos, hi_pos, strand,
                                    str(args.upstream), str(args.downstream)])
            SeqIO.write(promoter, outfile, 'fasta')


if __name__ == "__main__":
    main()


