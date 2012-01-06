from __future__ import division
from collections import Counter, defaultdict
FBTR_LEN = 11

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
    ids = set()
    posns = set()
    perfect_matches = defaultdict(list)

    exons, strands = read_gtf('Reference/melpsevir-all.gtf')
    has_printed = []

    for line in stdin:
        data = line.split()
        if not ((data[5] == '100M') 
                or ('99M' in data[5]) 
                or ('98M' in data[5])):
            continue
        posn = map_to_strand(data[2], int(data[3]), exons, strands)

        ids.add(data[0])
        posns.add(posn)

        if len(ids) % 2e5 == 0 and len(ids) not in has_printed:
            print len(posns) / len(ids), ',',
            has_printed.append(len(ids))
            stdout.flush()

