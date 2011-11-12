from __future__ import print_function, division
from collections import defaultdict
from sys import argv
from numpy import array, shape, argmax, var
import os
from matplotlib.pyplot import show, legend, errorbar, gca, savefig, title, \
        close, figure
from os import path

def check_data(gene, transcripts, intervals):
    if len(transcripts) <= 1:
        # Only one transcript isn't interesting and/or will be picked up on the
        # gene scans
        return

    transcripts = array(transcripts)
    n_transcripts, n_samples = shape(transcripts)
    best = argmax(transcripts) // n_samples
    for col in range(n_samples):
        this = argmax(transcripts[:, col])
        if this != best and var(transcripts[:, col]) > 0:
            if ((transcripts[this, col] - intervals[this][0, col]) 
                > (transcripts[best, col] + intervals[best][1, col])):
                print("-"*72)
                print(gene)
                print(col)
                print(transcripts)
                print(intervals)
                return gene
            # Flag it and get out quickly.  Let the hoo-mans sort this gene out.

def main():
    gn_to_tr = defaultdict(list)
    tr_expr = {}

    gn_col = 4
    tr_col = 0
    first_fpkm_col = 9
    first_lo_col = 10
    first_hi_col = 11
    num_cols = 4

    for line in open(argv[1]):
        line = line.split()
        gene = line[gn_col]
        transcript = line[tr_col]
        if not gene.startswith('FBgn'):
            continue

        gn_to_tr[gene].append(transcript)
        expr = array([float(i) for i in  line[first_fpkm_col::num_cols]])
        los = expr - array([float(i) for i in line[first_lo_col::num_cols]])
        his = array([float(i) for i in line[first_hi_col::num_cols]]) - expr
        tr_expr[transcript] = array([expr, los, his])

    try:
        os.makedirs(argv[2])
    except OSError:
        pass

    for gene in gn_to_tr:
        exprs = [tr_expr[tr][0, :] for tr in gn_to_tr[gene]]
        intervals = [tr_expr[tr][1:,:] for tr in gn_to_tr[gene]]
        if check_data(gene, exprs, intervals):
            plot_gene(gene, gn_to_tr, tr_expr)
            savefig(path.join(argv[2], gene+'.png'))
            close()




def plot_gene(gene, gn_to_tr, tr_expr):
    figure()
    for transcript in gn_to_tr[gene]:
        expr = tr_expr[transcript][0, :]
        errs = tr_expr[transcript][1:, :]
        errorbar(range(len(expr)), expr, yerr = errs, label=transcript)
    axes = gca()
    axes.set_xlim((-0.5, len(expr) - .5))
    #show()
    title(gene)
    legend()



'''

def main():
    infile = open(argv[1])
    gene_column = 4
    first_fpkm = 9
    cols_per_sample = 4

    #print(infile.readline(), end="")
    infile.readline()


    current_gene = ""
    transcript_data = []
    for line in infile:
        line = line.split()
        if line[gene_column] != current_gene:
            check_data(current_gene, transcript_data)
            transcript_data = []
            current_gene = line[gene_column]
        transcript_data.append(map(float, line[first_fpkm::cols_per_sample]))
    check_data(current_gene, transcript_data)

'''

if __name__ == "__main__":
    main()
