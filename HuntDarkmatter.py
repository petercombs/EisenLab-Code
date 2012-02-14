from __future__ import division
from os import path
import pysam
import numpy as np
import sys
import met
import inspect_shell




if __name__ == "__main__":
    GFF_file = sys.argv[1]
    print [path.basename(name) for name in sys.argv[2:]]
    bamfiles = [pysam.Samfile(name, 'rb') 
                for name in sys.argv[2:] if name.endswith('.bam')]

    counts = []

    gff_len = max(i for i, l in enumerate(open(GFF_file))) + 1

    for line in open(GFF_file):
        try:
            if line.startswith("#"): continue

            data = line.split('\t')
            chrom = 'dmel_'+data[0]
            start = int(data[3]) - 1
            stop = int(data[4])
            ID = data[-1].strip()

            counts = []
            for bamfile in bamfiles:
                counts.append(len(set(read.pos for read
                    in bamfile.fetch(chrom, start, stop))))

            N = sum(counts)
            if N > 1000:
                print "="*30+"\n"
                print "Large N on %s\n"% ID
                print "Counts: %s\n" % counts
                continue
            expect = [1/len(counts) for count in counts]
            m = met.Multinom(expect, counts)
            pval =  m.twosided_exact_test() 
            print '-'*30
            print ('SIG' if pval < (.05 / gff_len) else 'NON'), pval
            print ID
            print counts
            sys.stdout.flush()
        except Exception as err:
            print "E"*30
            print ID
            print counts
            print err

                           


            

