from __future__ import print_function, division
from Bio import SeqIO, Seq
from collections import Counter
import sys

def prob(kmer, bases):
    ret = 1.0
    for base in kmer:
        ret *= bases[base]
    return ret



if __name__ == "__main__":
    k = int(sys.argv[1])
    for fname in sys.argv[2:]:
        print('-'*72, fname, '-'*72, sep='\n')
        bases = Counter()
        kmers = Counter()
        for seqrec in SeqIO.parse(open(fname), 'fasta'):
            bases += Counter(seqrec.seq)
            
            for i in range(len(seqrec.seq) - k):
                seq = seqrec.seq[i:i+k]
                rc = str(seq.reverse_complement())
                if rc in kmers:
                    kmers[rc] += 1
                else:
                    kmers[str(seq)] += 1

        n_bases = sum(bases.itervalues())
        n_kmers = sum(kmers.itervalues())
        for base in bases:
            bases[base] /= float(n_bases)

        print(n_kmers, bases)
        for kmer, num in kmers.most_common()[:20]:
            e_val = prob(kmer, bases) * n_kmers
            if ((kmer.startswith('G') or kmer.startswith('T'))
                and (kmer.endswith('G') or kmer.endswith('T'))):
                kmer = Seq.reverse_complement(kmer)
            print(kmer, num, e_val, num / e_val)
        print('...')
        for kmer, num in kmers.most_common()[-20:]:
            e_val = prob(kmer, bases) * n_kmers
            print(kmer, num, e_val, num / e_val)
