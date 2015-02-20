#!/usr/bin/env python
from __future__ import division, print_function
from collections import defaultdict
import pysam
import re
import sys
from glob import glob
from os import path
from scipy import stats
from numpy import array, log10, exp
import progressbar as pbar

if sys.version_info.major < 3:
    res = raw_input('Code very memory heavy in less than Python3. Type "Y" to continue')
    assert res[0].upper() == 'Y'

if len(sys.argv) < 3:
    print("""Usage: python3 CheckCoverage.py <GTF-File> BAMfile [BAMfile ...]""")
    sys.exit(1)

gtf_fname = sys.argv[1]  # 'Reference/dmel-all-r5.42.gtf'
analysis_dir = 'analysis'

#starts = set()
#curr_len = -1
#coverage = 0

#cutoff = 0
fbtr_finder = re.compile('FBtr[0-9]*')


def analyze_bamfile(bam_fname):
    bam_file = pysam.Samfile(bam_fname, 'rb')

    coverages = defaultdict(lambda: [0, 0, set()])
    #parent = ''
    #curr_len = 0
    #coverage = 0
    #starts = set()

    f = open(gtf_fname)
    for i, line in enumerate(f):
        pass
    pb = pbar.ProgressBar(widgets=[bam_fname, pbar.Bar(), pbar.ETA()],
                          maxval=i).start()
    f.seek(0)

    current_fbtr = None
    last_fbtr = None
    #from guppy import hpy
    #hp = hpy()
    #before = hp.heap()
    for i, line in enumerate(f):
        pb.update(i)
        if line.startswith('#'):
            continue
        if line.startswith('>'):
            break
        data = line.split()
        chrom = data[0]
        kind = data[2]
        start = int(data[3]) - 1
        stop = int(data[4])
        # parent = fbtr_finder.findall(line)[0]

        if kind == 'exon':
            fbtrs = fbtr_finder.findall(line)
            if not fbtrs:
                continue
            fbtr = fbtrs[0]
            coverages[fbtr][1] += (stop - start)

            starts = set()
            coverage = 0
            for read in bam_file.fetch(chrom, start, stop):
                coverage += 1
                starts.add(read.pos)
            coverages[fbtr][0] += coverage
            try:
                coverages[fbtr][2].update(starts)
            except AttributeError:
                print(last_fbtr, current_fbtr, fbtr)
                continue
            del(starts)
            if current_fbtr and fbtr != current_fbtr:
                coverages[current_fbtr][2] = len(coverages[current_fbtr][2])
                last_fbtr = current_fbtr
                current_fbtr = fbtr
            elif current_fbtr is None:
                current_fbtr = fbtr

    pb.finish()
    coverages[fbtr][2] = len(coverages[fbtr][2])
    rpks, curr_lens, uniques = zip(*coverages.values())
    dir, fname = path.split(bam_fname)

    return bam_fname, rpks, uniques, curr_lens

if __name__ == "__main__":
    import multiprocessing as mp

    POOL = mp.Pool(20)
    res = POOL.map(analyze_bamfile,
                   [f for f in sys.argv[2:] if f.endswith('.bam')])
    #res = map(analyze_bamfile, [f for f in sys.argv[2:] if f.endswith('.bam')])
    all_dirs, all_rpks, all_pct_uniques, all_lens = zip(*res)

    import pickle
    #out_fh = open('checkcoverage.pkl', 'w')
    #pickle.dump({'dirs': all_dirs, 'rpks': all_rpks,
                 #'pct_uniques': all_pct_uniques, 'lens': all_lens},
                #out_fh)
    for fname, rpks, uniques, curr_lens in res:
        print(fname)
        try:
            xs = array([rpk/(curr_len+1)
                       for rpk, curr_len in zip(rpks, curr_lens)])
            ys = array([(u)/(curr_len + 1)
                       for u, curr_len in zip(uniques, curr_lens)])
            #cutoff = max(xs[ys < .1])
            #reg = stats.linregress(log10(xs[(xs < cutoff)*(xs > 0)*(ys > 0)]),
                                   #log10(ys[(xs < cutoff)*(xs > 0)*(ys > 0)]))
            reg = stats.linregress(log10(xs[ys < .1]),
                                   log10(ys[ys < .1]))
            print("exp(%f) * x ** %f" % (reg[1], reg[0]))
            print("Duplicate badness score: ", 10**(-reg[1]-.031))
        except Exception as exc:
            print(exc)
