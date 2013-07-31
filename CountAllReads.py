from os import path
from glob import glob
from subprocess import Popen, PIPE

for dirname in sorted([d for d in glob('sequence/*') if path.isdir(d)]):
    try:
        fnames = sorted(glob(path.join(dirname, '*_R1_*.fastq.gz')))
        n = 4e6 * (len(fnames) - 1)
        extra_gzip = Popen(['gzip', '-d', '-c', fnames[-1]], stdout=PIPE)
        extra_wc = Popen(['wc', '-l'], stdin=extra_gzip.stdout, stdout=PIPE)
        n += int(extra_wc.communicate()[0].strip()) / 4
        print dirname, '{:,}'.format(n)
    except:
        print dirname, "ERR!"

