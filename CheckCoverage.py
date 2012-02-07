from __future__ import division
import pysam
import re
from glob import glob
from os import path

gtf_fname = 'Reference/dmel-all-r5.42.gtf'
analysis_dir = 'analysis42'

starts = set()
curr_len = -1
coverage = 0

all_dirs = []
all_rpks = []
all_pct_uniques = []

for bam_fname in glob(path.join(analysis_dir, '*', 'accepted_hits.bam')):
    print bam_fname
    bam_file = pysam.Samfile(bam_fname, 'rb')

    coverages = {}

    for line in open(gtf_fname):
        if line.startswith('#'): continue
        data = line.split()
        chrom = data[0]
        kind = data[2]
        start = int(data[3]) - 1
        stop = int(data[4])
        fbtr_finder = re.compile('FBtr[0-9]*')
        if kind == 'mRNA':
            parent = fbtr_finder.findall(line)[0]
            curr_len = 0
            coverage = 0
            starts = set()

        elif kind == 'exon':
            fbtr = fbtr_finder.findall(line)[0]
            curr_len += (stop - start)
            if fbtr != parent: continue
            for read in bam_file.fetch(chrom, start, stop):
                coverage += 1
                starts.add(read.pos)

        elif kind == 'CDS':
            coverages[parent] = (coverage/curr_len, len(starts)/curr_len)

    rpks, pct_uniques = zip(*coverages.itervalues())
    dir, fname = path.split(bam_fname)
    all_dirs.append(dir)
    all_rpks.append(rpks)
    all_pct_uniques.append(pct_uniques)

                


