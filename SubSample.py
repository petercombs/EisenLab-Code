from glob import glob
from pysam import Samfile
from random import random, randint
from numpy import histogram, linspace
from heapq import heappush, heapreplace
from os import path, makedirs



n = int(1e7)

files = []
for pat in "*_V?? *_V???".split():
#for pat in "B??".split():
    files.extend(sorted(glob('analysis/' + pat + '/accepted_hits_sorted.bam')))
print files
for fn in files:
    n = min(n, Samfile(fn).mapped)

print "Keeping {} reads".format(n)
for fn in files:
    sample = [] 
    count = 0
    outdir = path.dirname(fn) + '_subset'
    try:
        makedirs(outdir)
    except OSError:
        pass
    print fn, '->', outdir
    sf = Samfile(fn)
    try:
        i_weight = float(sf.mapped)/n
    except ValueError:
        for i_weight, read in enumerate(sf): 
            pass
        print i_weight
        i_weight /= float(n)
        sf.seek(0)

    print fn, count, i_weight
    for i, read in enumerate(sf):
        key = random()**i_weight
        if len(sample) < n:
            heappush(sample, (key, read, i+count))
        else:
            heapreplace(sample, (key, read, i+count))

    count += i
    print "Kept {: >12,} of {: >12,} reads".format(len(sample), count)
    of = Samfile(path.join(outdir, 'accepted_hits.bam'), mode='wb', template=sf)
    sample.sort(key=lambda (key, read, pos): (read.tid, read.pos))
    for key, read, pos in sample:
        of.write(read)
    sf.close()
    of.close()
