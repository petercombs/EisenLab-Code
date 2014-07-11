from glob import glob
from pysam import Samfile
from numpy.random import random, randint
from numpy import histogram, linspace
from heapq import heappush, heappushpop
from os import path, makedirs
from sys import argv, stdout



def main():
    from itertools import repeat
    from multiprocessing import Pool
    n = int(1e7)

    files = []
    for pat in "*_V?? *_V???".split():
        glob_pat = path.join(argv[1], pat, 'accepted_hits_sorted.bam')
        print glob_pat
                                 
        files.extend(sorted(glob(glob_pat)))
    print files
    for fn in files:
        try:
            n = min(n, Samfile(fn).mapped)
        except ValueError:
            i = 0
            for  read in Samfile(fn):
                i += 1
            n = min(n, i)
    print "Keeping {} reads".format(n)
    stdout.flush()
    p = Pool(20)
    return p.map(subsample, zip(files, repeat([n, 3e6, 5e6, 7.5e6, 10e6])))
    #return    map(subsample, zip(files, repeat([n, 3e6, 5e6, 7.5e6, 10e6])))
    # Note: Comment out line -3 and -2, and uncomment line -1 to de-parallelize

def subsample(fn, ns=None):
    if ns is None:
        fn, ns = fn
    sample = []
    count = 0
    outdir_base = path.join(path.dirname(fn), 'subset')
    sf = Samfile(fn)
    try:
        i_weight = float(sf.mapped)/max(ns)
        print "Read out ", i_weight
    except ValueError:
        i_weight=0.0
        for read in sf:
            i_weight += 1
        print "Counted ", i_weight
        i_weight /= float(max(ns))
        sf = Samfile(fn)

    print fn, count, i_weight
    for i, read in enumerate(sf):
        key = random()**i_weight
        if len(sample) < max(ns):
            heappush(sample, (key, read, i+count))
        else:
            heappushpop(sample, (key, read, i+count))

    count += i

    for n in ns:
        if n == min(ns):
            outdir = outdir_base + '_min'
        else:
            outdir = outdir_base + '{:04.1f}M'.format(n/1e6)
        try:
            makedirs(outdir)
        except OSError:
            pass
        sampN = sorted(sample, reverse=True)[:int(n)]
        print "Kept {: >12,} of {: >12,} reads".format(len(sampN), count)
        print fn, '->', outdir
        stdout.flush()
        of = Samfile(path.join(outdir, 'accepted_hits.bam'), mode='wb', template=sf)
        sample.sort(key=lambda (key, read, pos): (read.tid, read.pos))
        for key, read, pos in sampN:
            of.write(read)
        of.close()
    sf.close()
    return [count for key, read, count in sample]

if __name__ == "__main__":
    subset_results = main()
