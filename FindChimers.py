from pysam import Samfile
from collections import defaultdict
from glob import glob
from os import path

reads = defaultdict(lambda : [None, None])
seqs = defaultdict(lambda : [None, None])
chimlist = open('chimlist.txt', 'w')
for fname in glob('assigned_*.bam'):
    if fname.count('_') > 1:
        continue
    print fname
    samfile = Samfile(fname)
    for read in samfile:
        qn = read.qname
        if reads[qn][read.is_read2] != None:
            print qn
            assert False
        reads[qn][read.is_read2] = fname
        seqs[qn][read.is_read2] = read.query
        if reads[qn][0] != None and reads[qn][1] != None:
            if reads[qn][0] != reads[qn][1]:
                chimlist.write("%s\t%s\t%s\n" % (qn, seqs[qn][0], seqs[qn][1]))
            del reads[qn]
            del seqs[qn]


