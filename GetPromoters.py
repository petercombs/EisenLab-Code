#from sys import stdout
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
    parser.add_argument('tss', type=file,
                        help='Tab-delimited file with:\n'
                        'Chrom\tStart\tStop\tStrand\tid')
    parser.add_argument('targets', type=file, nargs="+",
                        help='List of TSSs of interest')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    chroms = {chrom.id:chrom for chrom in SeqIO.parse(args.sequence, 'fasta')}
    TSSs = {line.split()[-1] : [to_number(d) for d in line.split()]
            for line in args.tss}

    for target in args.targets:
        print target.name
        outfile = open(splitext(target.name)[0] + '.fa', 'w')
        for line in target:
            tss = TSSs[line.strip()]
            if tss[STRAND] == '+':
                start = tss[START]
                promoter = chroms[str(tss[CHROM])][start - args.upstream
                                              : start + args.downstream]
            else:
                start = tss[STOP]
                promoter = chroms[str(tss[CHROM])][start - args.downstream
                                              : start + args.upstream]
            promoter.id = line.strip()
            SeqIO.write(promoter, outfile, 'fasta')


if __name__ == "__main__":
    main()


