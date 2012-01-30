import re
import copy
from collections import defaultdict

chrs = defaultdict(lambda : defaultdict(list))
fbgn_finder = re.compile('FBgn[0-9]+')
fbtr_finder = re.compile('FBtr[0-9]+')
for line in open('dmel-all-r5.42.gtf'):
    if line.startswith('#'): continue
    data = line.split()
    chr = data[0].replace('dmel_', 'chr')
    if fbgn_finder.findall(line) and fbtr_finder.findall(line):
        fbtr_to_fbgn[fbtr_finder.findall(line)[0]] = fbgn_finder.findall(line)[0]
    elif data[2] == 'exon':
        start = int(data[3])
        stop = int(data[4])
        fbtr = fbtr_finder.findall(line)[0]
        if fbtr not in fbtr_to_fbgn:
            print fbtr, 'has no fbgn designation'
            continue
        fbgn = fbtr_to_fbgn[fbtr]
        strand = data[6]
        for i in range(start, stop+1):
            if (fbgn, strand) not in chrs[chr][i]:
                chrs[chr][i].append((fbgn, strand))


trimmed = open('trimmed_exons.gtf', 'w')
for chrname, chr in chrs.iteritems():
    for pos in sorted(chr.keys()):
        for fbgn, strand in copy.copy(chr[pos]):
            trimmed.write('%s\texon\t%d\t' % (chrname, pos))
            while (fbgn, strand) in chr[pos]:
                chr[pos].remove((fbgn, strand))
                pos += 1
            trimmed.write('%d\t.\t%s\t.\tgene_id="%s"\n' % (pos, strand, fbgn))
            
