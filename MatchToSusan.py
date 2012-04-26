import os
import sys
from os import path
from collections import defaultdict
from scipy.stats import norm
from numpy import isinf, log, finfo, float64, mean, std


min_val = log(finfo(float64).tiny)
susan_exprs = {}

for fname in os.listdir('susan/by_cycle/'):
    fbase = fname[:-4]
    susan_exprs[fbase] = defaultdict(lambda : (0,100000))
    for line in open(path.join('susan/by_cycle/', fname)):
        if 'NAME' in line: continue
        data = line.split()
        gene_name = data[0]
        mean_expr = mean([float(i) for i in data[1:]])
        std_expr = std([float(i) for i in data[1:]]) or 1000
        susan_exprs[fbase][gene_name] = (mean_expr, std_expr)
        

fout = raw_input("Output file name? Leave blank for none. ")

if fout:
    fout = open(fout, 'w')

likelihoods = {time:0 for time in susan_exprs}

for line in open(sys.argv[1]):
    if "_name" in line: continue
    data = line.split()
    gene_name = data[0]
    mean_expr = mean([float(i) for i in data[1:]])
    for time in likelihoods:
        L = norm.logpdf(mean_expr, *susan_exprs[time][gene_name])
        if isinf(L):
            # Not ideal, but since we're subtracting the minimum value, having
            # something in range will always contribute.
            print gene_name, time, mean_expr, susan_exprs[time][gene_name]
            assert False
            continue
        else:
            likelihoods[time] += L
            if fout:
                fout.write('%f\t%s - %s\t%f\t%f -- %f\n' % (L, gene_name, time,
                                                       mean(susan_exprs[time][gene_name]),
                                                       std(susan_exprs[time][gene_name]),
                                                       mean_expr))

min_val = min(val for val in likelihoods.itervalues())
max_val = max(val for val in likelihoods.itervalues())

for time in likelihoods:
    likelihoods[time] -= min_val
    likelihoods[time] /= (max_val - min_val)

print likelihoods

if fout: fout.close()
