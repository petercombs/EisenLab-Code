import os
import sys
from os import path
from collections import defaultdict, Counter
from scipy.stats import norm
from numpy import isinf, isnan, log, finfo, float64, mean, std, median, abs
import progressbar as pb


pbar_widgets = [pb.Percentage(), ' ', pb.Bar(), ' ', pb.ETA()]
min_val = log(finfo(float64).tiny)
susan_exprs = {}

susan_files = os.listdir('susan/by_cycle')
pbar = pb.ProgressBar(widgets=["Loading Susan's Data: "] + pbar_widgets, maxval=len(susan_files))

for fname in pbar(susan_files):
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

likelihoods = Counter()
allLs = defaultdict(list)
allZs = defaultdict(list)
n = Counter()
bads = Counter()

f = open(sys.argv[1])
pbar = pb.ProgressBar(widgets=["Calculating Liklihoods: "]+pbar_widgets, maxval=len(f.read()))
pbar.start()

f.seek(0)

for line in f:
    pbar.update(f.tell())
    if "name" in line.lower(): continue
    data = line.split()
    gene_name = data[0]
    mean_expr = mean([float(i) for i in data[1:]])
    std_expr = std([float(i) for i in data[1:]])
    if mean_expr == 0 or std_expr == 0: continue
    if gene_name.startswith('CR'): continue
    for time in susan_exprs:
        susan_mean = mean(susan_exprs[time][gene_name])
        susan_std = std(susan_exprs[time][gene_name])
        if susan_mean == susan_std: continue
        L = log(norm.cdf(mean_expr + std_expr, susan_mean, susan_std+.01)
                - norm.cdf(mean_expr - std_expr, susan_mean, susan_std+.01)
                + 1e-200)
        Z = abs((mean_expr - susan_mean)/(susan_std+.001))
        #L = log(min(norm.cdf(mean_expr, susan_mean, susan_std+.01),
        #            norm.sf(mean_expr, susan_mean, susan_std+.01))
        #        + 1e-200)

        if L < -100 or isinf(L) or isnan(L):
            bads[time] += 1

        if isinf(L) or isnan(L):
            # Not ideal, but since we're subtracting the minimum value, having
            # something in range will always contribute.
            #print gene_name, time, mean_expr, susan_exprs[time][gene_name]
            print mean_expr, std_expr
            print susan_mean, susan_std
            print norm.cdf(mean_expr + std_expr, susan_mean, susan_std+.01)
            print norm.cdf(mean_expr - std_expr, susan_mean, susan_std+.01)
            assert False
            continue
        else:
            likelihoods[time] += L
            n[time] += 1
            allLs[time].append(L)
            allZs[time].append(Z)
            if fout:
                fout.write('%f\t%s - %s\t%f +/- %f -- %f\n' % (L, gene_name,
                                                               time,
                                                               susan_mean, susan_std,
                                                               mean_expr))


print

min_val = min(val for val in likelihoods.itervalues())
max_val = max(val for val in likelihoods.itervalues())

for key in sorted(likelihoods.keys()):
    print '\t', key, likelihoods[key]/n[key], n[key], bads[key],
    print median(allLs[key]), median(allZs[key])

#for time in likelihoods:
    #likelihoods[time] -= min_val
    #likelihoods[time] /= (max_val - min_val)
##
#for key in sorted(likelihoods.keys()):
    #print '\t', key, likelihoods[key]
#
if fout: fout.close()
