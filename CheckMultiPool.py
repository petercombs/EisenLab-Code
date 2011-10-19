from __future__ import print_function, division

import sys

from Bio import AlignIO
from os import path
from glob import glob
from collections import defaultdict
from numpy import mean
from jellyfish import hamming_distance as distance

import multiprocessing as mp
import pickle


def is_different(aln, species1, species2, length=40):
    strings1 = [str(aln[r].seq).replace('-', 'X') for r in species1]
    strings2 = [str(aln[r].seq).replace('-', 'Y') for r in species2]
    for pos in range(aln.get_alignment_length() - length):
        for row1 in strings1:
            for row2 in strings2:
                if row1[pos:pos + length] == row2[pos:pos+length]:
                    return False
    return True

def count_ambiguous_stretches(aln, species1, species2, length=40):
    stretches = 0
    strings1 = [str(aln[r].seq) for r in species1]
    strings2 = [str(aln[r].seq) for r in species2]
    for pos in range(aln.get_alignment_length() - length):
        for row1 in strings1:
            for row2 in strings2:
                if ((distance(row1[pos:pos+length], row2[pos:pos+length]) < 4)
                    and row1[pos:pos+length].count('-') != length):
                    stretches += 1
                    break
            else:
                continue
            # If there was an ambiguous base at this position, break out of both
            # loops
            break
    return stretches

def find_ambiguous_stretches(aln, species1, species2, length=40):
    stretches = []
    strings1 = [str(aln[r].seq) for r in species1]
    strings2 = [str(aln[r].seq) for r in species2]
    for pos in range(aln.get_alignment_length() - length):
        for row1 in strings1:
            for row2 in strings2:
                if row1[pos:pos+length] == row2[pos:pos+length]:
                    stretches.append((pos, row1, row2))
                    break
            else:
                continue
            # If there was an ambiguous base at this position, break out of both
            # loops
            break
    return stretches

def map_to_args(args):
    return count_stretches_in_file(*args)

def count_stretches_in_file(args):
    fname, expr_dict = args
    melname = path.basename(fname.replace('-aligned.fasta', ''))
    expr = mean(expr_dict[melname]) if melname in expr_dict else 1
    expr = expr or 1
    print((melname, expr), file=sys.stderr)

    ambiguous = defaultdict(lambda: defaultdict(int))
    total = defaultdict(int)

    try:
        aln = AlignIO.read(fname, 'fasta')
    except ValueError:
        # At least one "alignment" is empty because it only had mel sequence
        return ambiguous, total

    species_rows = defaultdict(list)
    for i, transcript in enumerate(aln):
        species_rows[transcript.name[:4]].append(i)

    for spec1 in species_rows:
        for spec2 in species_rows:
            if spec1 is spec2:
                # Assuming species_rows is in the same order both times
                # through (I think this is safe), then this prevents both
                # duplicates and self-comparison.
                break


            # Use minimum of lengths to be conservative
            total[spec1] += expr * min(len(
                str(aln[s].seq).replace('-','')) for s in species_rows[spec1])
            total[spec2] += expr * min(len(
                str(aln[s].seq).replace('-','')) for s in species_rows[spec2])


            ambiguous_pair = count_ambiguous_stretches(aln, species_rows[spec1],
                                                 species_rows[spec2],
                                                 comp_length)

            ambiguous[spec1][spec2] += ambiguous_pair * expr
            ambiguous[spec2][spec1] += ambiguous_pair * expr

    return ambiguous, total

def print_summary(ambiguous, total):
    print("Percent Ambiguous")
    print('\t', end='')
    for col in total:
        print(col, end='\t')

    #print()
    #for r in total:
    #    print(r, end='\t')
    #    for c in total:
    #        print(ambiguous[r][c], end='\t')
    #    print()
    #
    for row in total:
        print(row, end='\t')
        for col in total:
            print('%.3f' % (100 * ambiguous[row][col]/total[row]), end='\t')
        print()




if __name__ == "__main__":
    data_dir = '/Users/pacombs/data/Orthologs/aligned/'
    expr_dict_file = '/Users/pacombs/data/susanexprdict.pkl'
    expr_dict = pickle.load(open(expr_dict_file))
    comp_length = 90

    total = defaultdict(int)
    ambiguous = defaultdict(lambda: defaultdict(int))

    pool = mp.Pool()
    files = glob(path.join(data_dir, '*.fasta'))
    res = pool.map(count_stretches_in_file, zip(files, [expr_dict]*len(files)),
                   chunksize = 4000)

    for ambig, lens in res:
        for key1 in ambig:
            total[key1] += lens[key1]

            for key2 in ambig[key1]:
                ambiguous[key1][key2] += ambig[key1][key2]

    print_summary(ambiguous, total)
