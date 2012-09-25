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
ARGS.GTF =  'Reference/dmel.gtf'
ARGS.refbase = 'Reference/AAA/'
ARGS.base_species = 'mel'
ARGS.notificationEmail = 'peter.combs@berkeley.edu'
ARGS.seq_dir = 'sequence'
ARGS.config_file = 'RunConfig.cfg'

########################################################################

BASE = Namespace()
BASE.tophat_base = ('tophat -p8 --no-novel-juncs --read-edit-dist 6 '
                '--report-secondary-alignments ')
BASE.cufflinks_base = 'cufflinks -p 8 -q -u '


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

def get_readfiles(cfg_data):
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
    for sample, libname in cfg_data['sample_to_lib']:
        print "Finding reads for", sample
        seq_dir = join(ARGS.seq_dir, 'Sample*' + libname + '*')
        read_1s = glob(join(seq_dir, "*_R1_*.fastq"))
        print read_1s
        read_2s = glob(join(seq_dir, "*_R2_*.fastq"))
        print read_2s
        readnames[sample] = [','.join(read_1s), ','.join(read_2s)]
    return readnames

def count_reads(read_files):
    ''' Count the reads in the files'''

    # Let's cheat and only look in the last file: the others should have 4M
    # per...
    wc_proc = Popen(['wc', '-l', read_files[-1]], stdout=PIPE)
    wcout, wcerr = wc_proc.communicate()
    first_reads = int(4e6) * (len(read_files) - 1)
    result = wcout.splitlines()[-1].strip()
    lines, fname = result.split()
    last_reads = int(lines) / 4
    return first_reads + last_reads

DATA = Namespace()
DATA.config_data = process_config_file(ARGS.config_file)
DATA.readnames = get_readfiles(DATA.config_data)

DATA.samples = DATA.config_data['samples']


DATA.assign_procs = []
DATA.num_reads = {}

TIMES = Namespace()
TIMES.start = time()

for libname, (rf1, rf2) in DATA.readnames.items():
    # Print the name of the files we're going through, as a progress bar
    print '-'*72, '\n', libname, '\n', '-'*72

    # Just grab the first file name (PE have the same number in both)
    rfs = rf1.split(',')
    DATA.num_reads[libname] = count_reads(rfs)

    od = join(ARGS.analysis_dir, libname)
    try:
        os.makedirs(od)
    except OSError:
        print "Directory '%s' already exists... shouldn't be a problem" % od

    # Figure out Read Group ID
    f = open(rf1.split(',')[0])
    l = f.readline()
    f.close()
    rgid = l.split(":")[0][1:]
    lane = l.split(":")[1]

    idxfile = join(ARGS.refbase, ARGS.base_species +
                   DATA.config_data['sample_to_carrier'][libname])

    # Do tophat
    print 'Tophatting...', '\n', '='*30
    GTF = join(ARGS.refbase, ARGS.base_species +
               DATA.config_data['sample_to_carrier'][libname] + '.gtf')
    commandstr =  (BASE.tophat_base + '-G %(GTF)s -o %(od)s --rg-library '
                   '%(library)s'
                   ' --rg-center VCGSL --rg-sample %(library)s'
                   ' --rg-platform'
                   ' ILLUMINA --rg-id %(rgid)s  --rg-platform-unit %(lane)s'
                ' %(idxfile)s %(rf1)s %(rf2)s'
           % {'GTF': GTF,
              'od': od,
              'idxfile': idxfile,
              'rf1': rf1,
              'rf2': rf2,
              'library': libname,
              'rgid': rgid,
              'lane': lane})
    print commandstr
    sys.stdout.flush()
    tophat_proc = Popen(commandstr.split())
    tophat_proc.wait()
    commandstr = ['nice' 'python', 'AssignReads2.py',
                  join(od, 'accepted_hits.bam')]
    DATA.assign_procs.append(Popen(commandstr))




# Stop the timing
TIMES.end = time()
print "Mapping elapsed time", timedelta(seconds = TIMES.end - TIMES.start)

for proc in DATA.assign_procs:
    proc.wait()

print "Assignment extra time", timedelta(seconds = time() - TIMES.end)

TIMES.sortstart = time()
fs = glob(join(ARGS.analysis_dir, '*', 'assigned_dmel.bam'))
sort = Popen(['parallel',
              'samtools sort {} -m %d {//}/dmel_sorted' %3e9, # 3GB of memory
              ':::'] + fs)

sort.wait()

TIMES.sortend = time()
print "Sorting time", timedelta(seconds = TIMES.sortend - TIMES.sortstart)

# Figure out how well everything mapped
DATA.mapped_reads = {}
DATA.all_bams = [join(ARGS.analysis_dir, sample_dir, 'dmel_sorted.bam')
            for sample_dir in ARGS.samples]

for sample, bam in zip(ARGS.samples, DATA.all_bams):
    print '='*30
    commandstr = ['samtools', 'flagstat', bam]
    print commandstr
    samtools_proc = Popen(commandstr, stdout=PIPE)
    samout, samerr = samtools_proc.communicate()

    for line in samout.splitlines():
        if "mapped" in line:
            DATA.mapped_reads[sample] = float(line.split()[0])
            print "% reads dmel in ", sample,
            print DATA.mapped_reads[sample] / DATA.num_reads[sample] * 100
            break



for sample in ARGS.samples:
    od = join(ARGS.analysis_dir, sample)
    GTF = join(ARGS.refbase, ARGS.base_species +
               DATA.config_data['sample_to_carrier'][sample] + '.gtf')
    commandstr = (BASE.cufflinks_base + '-G %(GTF)s -o %(od)s %(hits)s'
                  % dict(GTF=GTF, od=od, hits=join(od, 'dmel_sorted.bam')))

    print commandstr
    sys.stdout.flush()
    cufflinks_proc = Popen(commandstr.split())

    cufflinks_proc.wait()

print "Final time", timedelta(seconds=time() - TIMES.start)
print "Cufflinks time", timedelta(seconds=time() - TIMES.sortend)
