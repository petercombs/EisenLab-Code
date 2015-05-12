import pysam
import gzip
from sys import argv
from os import path
import shutil
import tempfile

if __name__ == "__main__":
    sf = pysam.Samfile(argv[1])
    rfs = [entry for entry in
           sf.header['PG'][0]['CL'].split()
           if entry.endswith('.gz') or entry.endswith('.fastq')][0]
    rfs = sorted(rfs.split(','))


    with open(argv[2]+'_R1.fastq.gz', 'wb') as outf:
        for rf in rfs:
            if rf.endswith('.gz'):
                with open(rf) as readfile:
                    shutil.copyfileobj(readfile, outf)
            else:
                with gzip.open(tempfile.TemporaryFile(), 'w+b') as reads:
                    for line in open(rf):
                        reads.write(line)
                    reads.seek(0)
                    shutil.copyfileobj(reads, outf)

    if any([path.exists(rf.replace('R1', 'R2')) for rf in rfs]):
        outf2 = gzip.open(argv[2]+'R2.fastq.gz', 'w')
        for rf in rfs:
            rf = rf.replace('R1', 'R2')
            if not path.exists(rf): continue

            if rf.endswith('.gz'):
                infile = gzip.open(rf)
            else:
                infile = open(rf)
            for line in infile:
                outf2.write(line)
        outf2.close()




