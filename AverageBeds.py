import pdb
from numpy import argmax, argmin, mean, Inf
from glob import glob
from collections import defaultdict, Counter 
import sys

def any(my_list, test):
    for elem in my_list:
        if elem == test:
            return True
    return False

files = [open(fname) for fname in glob('*.bed')]


chrom = sys.argv[1]

bedlines = [file.readline().split() for file in files]

for i, line in enumerate(bedlines):
    while line[0] != chrom:
        line = files[i].readline().split()
    bedlines[i] = line

starts, stops, exprs = [list(i) for i 
                        in zip(*[[int(el) for el in line[1:]] for line in bedlines])]
chroms = [line[0] for line in bedlines]
while any(chroms, chrom):
    lowest_pos = min(stops)
    print '%s\t%d\t%d\t%d' % (chrom, starts[0], min(stops), mean(exprs))
    #pdb.set_trace()

    starts = [lowest_pos for elem in starts]
    for i, stop in enumerate(stops):
        if stop == lowest_pos:
            newline = files[i].readline().split()
            chroms[i] = newline[0]
            if newline[0] == chrom:
                stops[i] = int(newline[2])
                exprs[i] = int(newline[3])

            else:
                stops[i] = Inf
                exprs[i] = 0



