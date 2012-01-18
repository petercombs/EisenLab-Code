from sys import stdin, stdout, argv
from collections import defaultdict, Counter
import cPickle as pickle

chrs = defaultdict(lambda : defaultdict(list))

for line in open('../Reference/melpsevir-all.gtf'):
    try:
        data = line.split('\t')
        chr = data[0]
        start = int(data[3])
        stop = int(data[4])
        fbgn = data[-1][-12:-2]
        for i in range(start, stop+1):
            if fbgn not in chrs[chr][i]:
                chrs[chr][i].append(fbgn)
    except Exception as err:
        print err

genes = Counter()
for i, line in enumerate(stdin):
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

outfh = open(argv[1]+'.pkl', 'w')

pickle.dump(genes, outfh)


