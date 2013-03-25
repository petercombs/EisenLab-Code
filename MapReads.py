"""
MapReads.py: Runs tophat to map RNA-seq reads, and assigns them to the
appropriate species using a sub-program.

Does its best to automatically calculate things like filenames, etc, based only
on the given indices.

"""
import sys

from glob import glob
from os.path import join, abspath
import os
from time import time
from datetime import timedelta
from subprocess import Popen, PIPE, call
from argparse import Namespace
from pysam import view as samview
import pandas as pd

ARGS = Namespace()
ARGS.analysis_dir = abspath('analysis-multi')
ARGS.base_species = 'Dmel'
ARGS.seq_dir = abspath('sequence')
ARGS.config_file = abspath('Parameters/RunConfig.cfg')
ARGS.ref_base = abspath('Reference/AAA/')
ARGS.star_params = abspath('Parameters/STAR_params.in')

########################################################################

BASE = Namespace()
BASE.map_base = ('STAR --parametersFiles {params} '
                 '--genomeDir {genome} '
                 '--readFilesIn {rf1} {rf2}')
BASE.assign_base = 'nice python AssignReads2.py {fname}'


########################################################################



def process_config_file(cfg_fname):
    """ Get data we want out of the configuration file"""
    cfg_file = pd.read_table(cfg_fname, index_col='Label')
    return cfg_file


def get_readfiles(args, mbepc, index):
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
    print "Finding reads for", mbepc, 'index', index
    seq_dir = join(args.seq_dir, 
                   'Sample_MBEPC{}*_index{}'.format(mbepc, index))
    read_1s = glob(join(seq_dir, "*_R1_*.fastq*"))
    print read_1s
    read_2s = glob(join(seq_dir, "*_R2_*.fastq*"))
    print read_2s
    return  [[abspath(r) for r in sorted(read_1s)], 
             [abspath(r) for r in sorted(read_2s)]]

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
DATA.readnames = {} 

DATA.samples = DATA.config_data.index

TIMES = Namespace()
TIMES.start = time()

TEMP = Namespace()
TEMP.assign_procs = []
TEMP.rezip_procs = []

#for libname, (rf1, rf2) in DATA.readnames.items():
for sample in DATA.config_data.index:


    TEMP.sample_reads = get_readfiles(ARGS, 
                                      DATA.config_data['MBEPC'][sample], 
                                      DATA.config_data['Index'][sample])
    DATA.readnames[sample] = TEMP.sample_reads
    rf1, rf2 = TEMP.sample_reads

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

    # Do mapping
    print 'Tophatting...', '\n', '='*30

    TEMP.old_dir = os.getcwd()
    os.chdir(TEMP.od)
    print DATA.config_data['CarrierSpecies'][sample]
    TEMP.genome = join(ARGS.ref_base, 
                       ARGS.base_species +
                       DATA.config_data['CarrierSpecies'][sample])
    print TEMP.genome
    TEMP.commandstr =  BASE.map_base.format(
        rf1 = ','.join(rf1),
        rf2 = ','.join(rf2),
        genome = TEMP.genome,
        params = ARGS.star_params)
    print TEMP.commandstr
    sys.stdout.flush()

    call(str(TEMP.commandstr).split())
    print "Converting to bam"
    print abspath('../../ToBamAssign.sh')
    sys.stdout.flush()

    TEMP.assign_procs.append(Popen(['sh' , '../../ToBamAssign.sh']))

    os.chdir(TEMP.old_dir)



# Stop the timing
TIMES.end = time()
print "Mapping elapsed time", timedelta(seconds = TIMES.end - TIMES.start)

for proc in TEMP.assign_procs:
    proc.wait()

print "Assignment extra time", timedelta(seconds = time() - TIMES.end)

print "Final time", timedelta(seconds=time() - TIMES.start)

import cPickle as pickle


call('STAR --parametersFiles {params} '
     '--genomeDir {genome} '
     '--genomeLoad Remove'.format(params=ARGS.star_params,
                                 genome=join(ARGS.ref_base)
                                 ).split())

pickle.dump(dict(data=DATA, args=ARGS),
            open('mapreads_dump.pkl', 'w'))

