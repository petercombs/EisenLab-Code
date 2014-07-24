from __future__ import print_function
from pysam import Samfile
from os import path
from glob import glob
from collections import Counter
import gzip
import shutil
from sys import argv

def fastq_write(outfile, read, force_r1=False):
    outfile.write('@{}/{}\n{}\n+\n{}\n'
                  .format(read.qname,
                          1 if force_r1 else [1,2][read.is_read2],
                          read.seq,
                          read.qual))


def convert_files(fn):
    """ Converts a set of files to fastq

    Makes a number of assumptions based on the structure of the files"""
    sf = Samfile(fn)
    call = sf.header['PG'][0]['CL']

    # Assumes that the directory will be the single most common thing that shows
    # up in the call line of the Bamfile
    # One could, in principle, try to parse this out manually in a less fragile
    # way, but I'm in a hurry.
    basename = Counter(call.split('/')).most_common(1)[0][0]

    orphans = [{}, {}]

    outname_r1 = path.join(
        #'/data3/fly/data/species_embryo/slice_carriers/',
        'carrier_fastqs2',
        basename+ '_R1')
    print("Converting {} to {}".format(fn, outname_r1))
    outname_r2 = outname_r1.replace('_R1', '_R2')
    outname_orphans = outname_r1.replace('_R1', '_orphans')
    fn2 = [f for f in glob(path.join(path.dirname(fn), 'ambig_d*.bam'))
           if 'mel' not in f][0]
    print(fn2)
    # Uncomment next line to just list files
    #continue
    outfile_r1 = gzip.GzipFile(outname_r1+'.fastq.gz', 'wb')
    outfile_r2 = gzip.GzipFile(outname_r2+'.fastq.gz', 'wb')
    outfile_orphans = gzip.GzipFile(outname_orphans+'.fastq.gz', 'wb')
    for read in sf:
        orphans[read.is_read2][read.qname] = read
        if read.qname in orphans[0] and read.qname in orphans[1]:
            fastq_write(outfile_r1, orphans[0].pop(read.qname))
            fastq_write(outfile_r2, orphans[1].pop(read.qname))

    outfile_r1.close()
    outfile_r2.close()

    while orphans[0]:
        fastq_write(outfile_orphans, orphans[0].popitem()[1])
    while orphans[1]:
        fastq_write(outfile_orphans, orphans[1].popitem()[1], force_r1=True)
    outfile_orphans.close()


    sf2 = Samfile(fn2)
    outfile_r1 = gzip.GzipFile(outname_r1+'_ambig.fastq.gz', 'wb')
    outfile_r2 = gzip.GzipFile(outname_r2+'_ambig.fastq.gz', 'wb')
    outfile_orphans = gzip.GzipFile(outname_orphans+'_ambig.fastq.gz', 'wb')
    for read in sf2:
        orphans[read.is_read2][read.qname] = read
        if read.qname in orphans[0] and read.qname in orphans[1]:
            fastq_write(outfile_r1, orphans[0].pop(read.qname))
            fastq_write(outfile_r2, orphans[1].pop(read.qname))

    outfile_r1.close()
    outfile_r2.close()

    outfile_r1.close()
    outfile_r2.close()

    while orphans[0]:
        fastq_write(outfile_orphans, orphans[0].popitem()[1])
    while orphans[1]:
        fastq_write(outfile_orphans, orphans[1].popitem()[1], force_r1=True)
    outfile_orphans.close()

    if path.exists(path.join(path.dirname(fn), 'unmapped.bam')):
        outfile_r1 = gzip.GzipFile(outname_r1+'_unmapped.fastq.gz', 'wb')
        outfile_r2 = gzip.GzipFile(outname_r2+'_unmapped.fastq.gz', 'wb')
        outfile_orphans = gzip.GzipFile(outname_orphans+'_unmapped.fastq.gz', 'wb')
        sf3 = Samfile(path.join(path.dirname(fn), 'unmapped.bam'))
        for read in sf3:
            orphans[read.is_read2][read.qname] = read
            if read.qname in orphans[0] and read.qname in orphans[1]:
                fastq_write(outfile_r1, orphans[0].pop(read.qname))
                fastq_write(outfile_r2, orphans[1].pop(read.qname))

        outfile_r1.close()
        outfile_r2.close()

        while orphans[0]:
            fastq_write(outfile_orphans, orphans[0].popitem()[1])
        while orphans[1]:
            fastq_write(outfile_orphans, orphans[1].popitem()[1], force_r1=True)
        outfile_orphans.close()

if __name__ == "__main__":
    fns = [f for f in glob('analysis-multi/*/assigned*.bam') if 'mel' not in f]

    if len(argv) > 1:
        fns = [f for f in fns if argv[1] in f]

    print(fns)


    import multiprocessing as mp

    p = mp.Pool(20)
    p.map(convert_files, fns)
    #for the_fn in fns:
        #convert_files(the_fn)
