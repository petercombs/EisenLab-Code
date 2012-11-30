from __future__ import division
from os import path
import pysam
import numpy as np
import sys
import met
import inspect_shell
from bx.bbi import bigwig_file as bigwig
import signal
import re
from scipy import stats

readlen = 50

def term_catcher(*args):
    raise KeyboardInterrupt("Sigterm")
signal.signal(signal.SIGTERM, term_catcher)

if __name__ == "__main__":
    GFF_file = sys.argv[1]
    fpkm_file = sys.argv[2]
    mapped_files = sys.argv[3:]
    print [path.basename(name) for name in mapped_files]
    bamfiles = [pysam.Samfile(name, 'rb') if name.endswith('.bam')
                else bigwig.BigWigFile(open(name, 'rb'))
                for name in mapped_files ]

    counts = []

    fpkms = {}
    n_samples = len(sys.argv) - 3
    for line in open(sys.argv[2]):
        try:
            data = line.strip().split('\t')
            fbgn = data[3]
            fpkm = [float(data[i]) for i in range(13,(13+4*n_samples),4)]
            fpkms[fbgn] = np.array(fpkm)
        except ValueError:
            pass

    gff_len = max(i for i, l in enumerate(open(GFF_file))) + 1

    for line in open(GFF_file):
        try:
            if line.startswith("#"): continue

            data = line.split('\t')
            chrom = 'dmel_'+data[0]
            start = int(data[3]) - 1
            stop = int(data[4])
            ID = data[-1].strip()
            try:
                fbgn = re.findall('FBgn[0-9]*', ID)[0]
            except IndexError:
                continue

            counts = []
            for bamfile in bamfiles:
                try:
                    counts.append(len(set(read.pos for read
                        in bamfile.fetch(chrom, start, stop))))
                except AttributeError:
                    counts.append(int(sum(bamfile.get_as_array(chrom,
                                                               start,
                                                               stop))/readlen))

            assert len(counts)
            N = sum(counts)
            N += (N == 0)
            try:
                fpkm = fpkms[fbgn]
            except KeyError:
                fpkm = np.ones(len(counts))
            if max(fpkm) == 0:
                fpkm += 1
            fpkm += 0.01 * min(fpkm[fpkm > 0])
            fpkm /= np.mean(fpkm)
            expect = [float(x*N/len(counts)) for count, x in zip(counts, fpkm)]
            if N > 250:
                print "="*30+"\n"
                print "Large N on %s"% ID
                print "Counts: %s" % counts
                print "Expect: %s\n" % expect
                sys.stdout.flush()
                X2, p = stats.chisquare(np.array(counts), np.array(expect))
                print p
            else:
                print '-'*30
                print ID
                print "Counts:", counts
                print "Expect:", expect
                sys.stdout.flush()
                m = met.Multinom(expect, counts)
                pval =  m.twosided_exact_test()
                print ('SIG' if pval < (.05 / gff_len) else 'NON'), pval

        except (ValueError,OverflowError,KeyboardInterrupt) as err:
            print "E"*30
            print err
            print ID
            print counts
            print expect


