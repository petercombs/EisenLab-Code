from sys import stdin, stdout, argv
from collections import defaultdict, Counter
from os import path
import cPickle as pickle
import subprocess
import os
import re

chrs = defaultdict(lambda : defaultdict(list))
fbgn_finder = re.compile('FBgn[0-9]+')
fbtr_finder = re.compile('FBtr[0-9]+')

fbtr_to_fbgn = {}

for line in open(argv[1]):
    try:
        data = line.split('\t')
        chr = data[0]
        kind = data[2]
        if fbgn_finder.findall(line) and fbtr_finder.findall(line):
            fbtr_to_fbgn[fbtr_finder.findall(line)[0]] = \
                fbgn_finder.findall(line)[0]
        if kind != 'exon': continue
        start = int(data[3])
        stop = int(data[4])
        try:
            fbgn = fbtr_to_fbgn[fbtr_finder.findall(line)[0]]
        except KeyError:
            fbgn = 'ERR'
        if chr not in chrs:
            print chr
        for i in range(start, stop+1):
            if fbgn not in chrs[chr][i]:
                chrs[chr][i].append(fbgn)
    except Exception as err:
        print err

genes = Counter()
for dir in os.listdir('.'):
    if not path.isdir(dir): continue
    print dir
    print '-'*30
    samtools = subprocess.Popen(['samtools', 'view', 
                                 path.join(dir, 'accepted_hits.bam')],
                                stdout=subprocess.PIPE)
    
    for i, line in enumerate(samtools.stdout):
        data = line.split()
        chr = data[2]
        pos = int(data[3])
        try:
            for gene in chrs[chr][pos]:
                genes[gene] += 1
        except KeyError:
            genes['OTHER'] += 1

        if i % 1e6 == 0:
            print '.',
            stdout.flush()

    outfh = open(dir + '_coverage.pkl', 'w')

    pickle.dump(genes, outfh)


