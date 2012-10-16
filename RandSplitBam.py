from pysam import Samfile
from sys import argv
from os import path
from random import randint
import progressbar as pbar

num_splits = int(argv[1])
for in_fname in argv[2:]:
    in_file = Samfile(in_fname)
    out_fname = path.splitext(in_fname)[0] + '_randsplit%i.bam'
    out_files = [Samfile(out_fname % i, 'w', in_file) for i in range(num_splits)]

    print "Measuring file size"
    start = in_file.tell()
    maxval = path.getsize(in_fname) * 2**16 # I don't know why it's off by 2^16

    pb = pbar.ProgressBar(maxval=maxval - start + 2**16,
                          widgets=[in_fname, '-',pbar.Percentage(), ':',
                                   pbar.Bar(), pbar.ETA()]) 
    pb.start()
    for read in pb(in_file):
        out_files[randint(0,num_splits-1)].write(read)
        pb.update(in_file.tell() - start)

    pb.finish()

