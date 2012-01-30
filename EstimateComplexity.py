from __future__ import division
from collections import Counter, defaultdict
from numpy import histogram
FBTR_LEN = 11
MAX_STOP = 5000000.

def read_gtf(fname):
    strands = {}
    exons = {}
    lens = {}

    for line in open(fname):
        if line.startswith('#'): continue
        data = line.split()
        chrom = data[0]
        kind = data[2]
        low = int(data[3])
        high = int(data[4])
        strand = data[6]
        annot = data[8]

        if annot.startswith('ID'):
            ID = annot[3:3+FBTR_LEN]
            exons[ID] = []
            strands[ID] = chrom+strand
        
        if kind == 'exon':
            ID = annot[7:8+FBTR_LEN]
            insert_pos = 0 if strand == '-' else len(exons[ID])
            exons[ID].insert(insert_pos, (low, high))

    return exons, strands

def map_to_strand(rname, pos, exons, strands):
    try:
        curr_exons = exons[rname][:] # Copy the list of exons
    except KeyError:
        #print rname, pos, 'Not Found in GTF'
        return None
    exon = curr_exons.pop(0)
    while pos > (exon[1] - exon[0] + 1):
        pos -= (exon[1] - exon[0] + 1)
        try:
            exon = curr_exons.pop(0)
        except IndexError:
            #print rname, pos, exons[rname]
            return None

    if strands[rname].endswith('+'):
        return strands[rname][:-1], exon[0] + pos
    else:
        return strands[rname][:-1], exon[1] - pos


if __name__ == "__main__":
    from sys import stdin, stdout
    import re

    match_string = re.compile('([0-9]+)M')
    ids = set()
    posns = set()
    best_matches = defaultdict(Counter)

    exons, strands = read_gtf('Reference/melpsevir-all.gtf')
    has_printed = []

    try:
        for line in stdin:
            data = line.split()
            posn = map_to_strand(data[2], int(data[3]), exons, strands)
            if not posn:
                continue

            species = posn[0][:4]
            best_matches[data[0]][species] = max(max([int(val) for val in
                                                  match_string.findall(data[5])]),
                                                 best_matches[data[0]][species])

            ids.add(data[0])
            posns.add(posn)

            if len(ids) % 2e5 == 0 and len(ids) not in has_printed:
                print len(posns) / len(ids), ',',
                has_printed.append(len(ids))
                stdout.flush()
            if len(ids) > MAX_STOP:
                raise KeyboardInterrupt
    except KeyboardInterrupt:
        pass

    
    best_differences = []
    species_count = Counter()

    for key in best_matches:
        if len(best_matches[key]) == 1:
            best_differences.append(100)
            species_count[best_matches[key].keys()[0]] += 1
        else:
            try:
                scores = best_matches[key].most_common()
                best_differences.append(scores[0][1] - scores[1][1])
                if scores[0][1] - scores[1][1] > 2:
                    species_count[scores[0][0]] += 1

            except IndexError:
                print key, scores

    print
    print species_count
    for spec in ['dmel', 'dpse', 'dvir']:
        print spec, species_count[spec]/MAX_STOP
    print (MAX_STOP - sum(species_count.values()))/MAX_STOP
    #print histogram(best_differences, 101)
    #print histogram(best_differences, 101, normed=True)
    from cPickle import dump
    print "dumping!"
    stdout.flush()
    dump(best_differences, open('difference_dump.pkl', 'w'))

