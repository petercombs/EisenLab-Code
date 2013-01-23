"""
MapReads.py: Runs tophat to map RNA-seq reads, and assigns them to the
appropriate species using a sub-program.

Does its best to automatically calculate things like filenames, etc, based only
on the given indices.

"""
import sys

from glob import glob
from os.path import join
import os
from time import time
from datetime import timedelta
from subprocess import Popen, PIPE
from argparse import Namespace

ARGS = Namespace()
ARGS.analysis_dir = 'analysis-multi'
ARGS.base_GTF =  'Reference/AAA/dmel-all-r5.46.gtf'
ARGS.refbase = 'Reference/AAA/'
ARGS.transbase = 'Reference/AAA/transcriptome/'
ARGS.base_species = 'mel'
ARGS.notificationEmail = 'peter.combs@berkeley.edu'
ARGS.seq_dir = 'sequence'
ARGS.config_file = 'RunConfig.cfg'

########################################################################

BASE = Namespace()
BASE.tophat_base = ('tophat -p8 --no-novel-juncs --read-edit-dist 6 '
                    '--report-secondary-alignments '
                    '--no-sort-bam ')
BASE.cufflinks_base = 'cufflinks -p 8 -q -u '
BASE.assign_base = 'nice python AssignReads2.py {fname}'


########################################################################



def process_config_file(cfg_fname):
    """ Get data we want out of the configuration file"""
    cfg_fh = open(cfg_fname)
    cfg_data = dict(sample_to_lib = [], samples=[], libraries = [],
                    sample_to_carrier = {})
    for line in cfg_fh:
        if line.startswith('#'):
            continue
        try:
            line = line.strip().split('\t')
            lib, mbepc, slice, idx, carrier = line
            cfg_data['sample_to_lib'].append((lib + slice,
                                              mbepc + '*index' + idx))
            cfg_data['sample_to_carrier'][lib + slice] = carrier
            cfg_data['samples'].append(lib + slice)
            cfg_data['libraries'].append(mbepc + slice)


        except ValueError:
            print ("Line '%s' should have exactly "
                   "5 tab-separated elements." % line)
            print "Continuing..."
    return cfg_data

def get_readfiles(args, sample, libname):
    " Find the names of the read files, based on configuration"

    # Directory structure:
    # Project_Eisen_M
    #   Sample_MBE_PC_64A_index1
    #     MBE_PC_64A_index1_ATCACG_L001_R1_001.fastq.gz
    #     MBE_PC_64A_index1_ATCACG_L001_R1_002.fastq.gz
    #     MBE_PC_64A_index1_ATCACG_L001_R1_003.fastq.gz
    #     ...
    #     MBE_PC_64A_index1_ATCACG_L001_R2_001.fastq.gz
    #     MBE_PC_64A_index1_ATCACG_L001_R2_002.fastq.gz
    #     MBE_PC_64A_index1_ATCACG_L001_R2_003.fastq.gz
    #     ...
    #   Sample_MBE_PC_64B_index2
    #     ...
    readnames = {}
    print "Finding reads for", sample
    seq_dir = join(args.seq_dir, 'Sample*' + libname + '*')
    read_1s = glob(join(seq_dir, "*_R1_*.fastq*"))
    print read_1s
    read_2s = glob(join(seq_dir, "*_R2_*.fastq*"))
    print read_2s
    readnames[sample] = [sorted(read_1s), sorted(read_2s)]
    return readnames

def count_reads(read_files):
    ''' Count the reads in the files'''

    # Let's cheat and only look in the last file: the others should have 4M
    # per...
    wc_proc = Popen(['wc', '-l', read_files[-1]], stdout=PIPE)
    wcout, _ = wc_proc.communicate()
    first_reads = int(4e6) * (len(read_files) - 1)
    result = str(wcout).splitlines()[-1].strip()
    lines, _ = result.split()
    last_reads = int(lines) / 4
    return first_reads + last_reads

DATA = Namespace()
DATA.config_data = process_config_file(ARGS.config_file)
DATA.readnames = {} #get_readfiles(DATA.config_data)

DATA.samples = DATA.config_data['samples']

TIMES = Namespace()
TIMES.start = time()

TEMP = Namespace()
TEMP.assign_procs = []
TEMP.rezip_procs = []

#for libname, (rf1, rf2) in DATA.readnames.items():
for sample, libname in DATA.config_data['sample_to_lib']:


    TEMP.sample_reads = get_readfiles(ARGS, sample, libname)
    DATA.readnames.update(TEMP.sample_reads)
    rf1, rf2 = TEMP.sample_reads[sample]

    # Print the name of the files we're going through, as a progress bar
    print '-'*72, '\n', sample, '\n', '-'*72
    sys.stdout.flush()

    # Skip dealing with the gzipped files... this has been in Tophat
    # since version 1.3.0

    TEMP.od = join(ARGS.analysis_dir, sample)
    try:
        os.makedirs(TEMP.od)
    except OSError:
        print ("Directory '%s' already exists... shouldn't be a problem" %
               TEMP.od)

    TEMP.idxfile = join(ARGS.refbase, ARGS.base_species +
                   DATA.config_data['sample_to_carrier'][sample])

    # Do tophat
    print 'Tophatting...', '\n', '='*30
    TEMP.GTF = join(ARGS.refbase, ARGS.base_species +
               DATA.config_data['sample_to_carrier'][sample] + '.gtf')

    TEMP.transcriptome = join(ARGS.transbase, ARGS.base_species +
                              DATA.config_data['sample_to_carrier'][sample])
    TEMP.transarg = '--transcriptome-index %s '% TEMP.transcriptome

    TEMP.commandstr =  (BASE.tophat_base + '-G %(GTF)s -o %(od)s --rg-library '
                   '%(library)s'
                   ' --rg-center VCGSL --rg-sample %(library)s'
                   ' --rg-id %(library)s '
                   ' --rg-platform ILLUMINA '
                   ' %(transarg)s'
                ' %(idxfile)s %(rf1)s %(rf2)s'
           % {'GTF': TEMP.GTF,
              'od': TEMP.od,
              'idxfile': TEMP.idxfile,
              'rf1': ','.join(rf1),
              'rf2': ','.join(rf2),
              'library': sample,
              'transarg' : TEMP.transarg})

    print TEMP.commandstr
    sys.stdout.flush()
    TEMP.tophat_proc = Popen(str(TEMP.commandstr).split())
    TEMP.tophat_proc.wait()

    TEMP.commandstr = BASE.assign_base.format(fname = join(TEMP.od,
                                                           'accepted_hits.bam'))
    print TEMP.commandstr
    TEMP.assign_procs.append(Popen(TEMP.commandstr.split()))




# Stop the timing
TIMES.end = time()
print "Mapping elapsed time", timedelta(seconds = TIMES.end - TIMES.start)

for proc in TEMP.assign_procs:
    proc.wait()

print "Assignment extra time", timedelta(seconds = time() - TIMES.end)

print "Final time", timedelta(seconds=time() - TIMES.start)

import cPickle as pickle

pickle.dump(dict(data=DATA, args=ARGS),
            open('mapreads_dump.pkl', 'w'))

