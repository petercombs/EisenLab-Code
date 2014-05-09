from glob import glob
from pysam import Samfile
from numpy.random import random, randint
from numpy import histogram, linspace
from heapq import heappush, heapreplace
from os import path, makedirs
from sys import stdout



def main():
    from itertools import repeat
    from multiprocessing import Pool
    n = int(1e7)

    files = []
    for pat in "*_V?? *_V???".split():
        files.extend(sorted(glob('analysis/' 
                                 + pat 
                                 + '/accepted_hits_sorted.bam')))
    print files
    for fn in files:
        n = min(n, Samfile(fn).mapped)
    print "Keeping {} reads".format(n)
    stdout.flush()
    p = Pool()
    p.map(subsample, zip(files, repeat(n)))

def subsample(fn, n=None):
    if n is None:
        fn, n = fn
    sample = [] 
    count = 0
    outdir = path.join(path.dirname(fn), 'subset')
    try:
        makedirs(outdir)
    except OSError:
        pass
    print fn, '->', outdir
    sf = Samfile(fn)
    try:
        i_weight = float(sf.mapped)/n
        print "Read out ", i_weight
    except ValueError:
        i_weight=0.0
        for read in sf: 
            i_weight += 1
        print "Counted ", i_weight
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
    stdout.flush()
    of = Samfile(path.join(outdir, 'accepted_hits.bam'), mode='wb', template=sf)
    sample.sort(key=lambda (key, read, pos): (read.tid, read.pos))
    for key, read, pos in sample:
        of.write(read)
    sf.close()
    of.close()

if __name__ == "__main__":
    main()
