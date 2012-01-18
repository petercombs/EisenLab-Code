from sys import stdin, stdout, argv
from collections import defaultdict, Counter
from os import path
import cPickle as pickle
import subprocess
import os

chrs = defaultdict(lambda : defaultdict(list))

for line in open(argv[1]):
    try:
        data = line.split('\t')
        chr = data[0]
        start = int(data[3])
        stop = int(data[4])
        fbgn = data[-1][-12:-2]
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

    outfh = open(dir + '.pkl', 'w')

    pickle.dump(genes, outfh)


