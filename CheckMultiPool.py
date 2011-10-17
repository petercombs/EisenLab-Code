from __future__ import print_function, division

import sys

from Bio import AlignIO
from os import path
from glob import glob
from collections import defaultdict
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
    strings1 = [str(aln[r].seq).replace('?', 'X') for r in species1]
    strings2 = [str(aln[r].seq).replace('?', 'Y') for r in species2]
    for pos in range(aln.get_alignment_length() - length):
        for row1 in strings1:
            for row2 in strings2:
                if row1[pos:pos+length] == row2[pos:pos+length] and \
                   row1[pos:pos+length].count('-') != length:
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
    strings1 = [str(aln[r].seq).replace('-', 'X') for r in species1]
    strings2 = [str(aln[r].seq).replace('-', 'Y') for r in species2]
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

if __name__ == "__main__":
    data_dir = '/Users/pacombs/data/Orthologs/aligned/'
    expr_dict_file = '/Users/pacombs/data/susanexprdict.pkl'
    expr_dict = pickle.load(expr_dict_file)
    comp_length = 40

    total_species = defaultdict(int)
    ambiguous_species = defaultdict(lambda: defaultdict(int))

    for alignment in glob(path.join(data_dir, '*.fasta')):
        melname = alignment.replace('-aligned.fasta', '')
        expr = expr_dict[expr] if expr in expr_dict else 1
        print((melname, expr), file=sys.stderr)
        try:
            aln = AlignIO.read(alignment, 'fasta')
        except ValueError:
            # At least one "alignment" is empty because it only had mel sequence
            continue

        species_rows = defaultdict(list)
        for i, r in enumerate(aln):
            species_rows[r.name[:4]].append(i)

        for s1 in species_rows:
            for s2 in species_rows:
                if s1 is s2:
                    # Assuming species_rows is in the same order both times
                    # through (I think this is safe), then this prevents both
                    # duplicates and self-comparison.
                    break


                # Use minimum of lengths to be conservative
                total_species[s1] += min(len(
                    str(aln[s].seq).replace('-','')) for s in species_rows[s1])
                total_species[s2] += min(len(
                    str(aln[s].seq).replace('-','')) for s in species_rows[s2])


                ambiguous = count_ambiguous_stretches(aln, species_rows[s1],
                                                     species_rows[s2],
                                                     comp_length) * expr

                ambiguous_species[s1][s2] += ambiguous * expr
                ambiguous_species[s2][s1] += ambiguous * expr

    print('\t', end='')
    for c in total_species:
        print(c, end='\t')

    print()
    for r in total_species:
        print(r, end='\t')
        for c in total_species:
            print(ambiguous_species[r][c], end='\t')
        print()

    print("Percent Ambiguous")
    for r in total_species:
        print(r, end='\t')
        for c in total_species:
            print('%.3f' % (100 * ambiguous_species[r][c]/total_species[r]), end='\t')
        print()



