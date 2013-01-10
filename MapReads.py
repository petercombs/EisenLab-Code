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
#BASE.tophat_base = ('bowtie2 -p 8 --all --no-mixed --local ')
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


TEMP.assign_procs = []
#DATA.num_reads = {}

TIMES = Namespace()
TIMES.start = time()

TEMP = Namespace()
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
    ## Unzip anything that's zipped
    #TEMP.to_unzip = []
    #for i, fname in enumerate(rf1):
        #if fname.endswith('.gz'):
            #print "Unzipping", fname
            #TEMP.to_unzip.append(fname)
            #rf1[i] = fname.strip('.gz')
#
    #for i, fname in enumerate(rf2):
        #if fname.endswith('.gz'):
            #print "Unzipping", fname
            #TEMP.to_unzip.append(fname)
            #rf2[i] = fname.strip('.gz')
#
    #if TEMP.to_unzip:
        #Popen(['parallel', '-j', '2', 'gunzip {}', ':::'] +
              #TEMP.to_unzip).wait()



    ## Just grab the first file name (PE have the same number in both)
    #TEMP.rfs = rf1
    ##DATA.num_reads[sample] = count_reads(TEMP.rfs)

    TEMP.od = join(ARGS.analysis_dir, sample)
    try:
        os.makedirs(TEMP.od)
    except OSError:
        print ("Directory '%s' already exists... shouldn't be a problem" %
               TEMP.od)

    # Figure out Read Group ID
    #TEMP.f = open(rf1[0])
    #TEMP.l = TEMP.f.readline()
    #TEMP.f.close()
    #TEMP.rgid = TEMP.l.split(":")[0][1:]
    #TEMP.lane = TEMP.l.split(":")[1]

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
                   ' --rg-platform ILLUMINA '
                   #' --rg-id %(rgid)s --rg-platform-unit %(lane)s'
                   ' %(transarg)s'
                ' %(idxfile)s %(rf1)s %(rf2)s'
           % {'GTF': TEMP.GTF,
              'od': TEMP.od,
              'idxfile': TEMP.idxfile,
              'rf1': ','.join(rf1),
              'rf2': ','.join(rf2),
              'library': sample,
              #'rgid': TEMP.rgid,
              #'lane': TEMP.lane,
              'transarg' : TEMP.transarg})

    print TEMP.commandstr
    sys.stdout.flush()
    TEMP.tophat_proc = Popen(str(TEMP.commandstr).split())
    TEMP.tophat_proc.wait()

    #TEMP.rezip_procs.append(Popen(['parallel', '-j', '2', 'gzip {}', ':::']
                                  #+ rf1 + rf2))

    TEMP.commandstr = ['nice', 'python', 'AssignReads2.py',
                  join(TEMP.od, 'accepted_hits.bam')]
    print TEMP.commandstr
    TEMP.assign_procs.append(Popen(TEMP.commandstr))




# Stop the timing
TIMES.end = time()
print "Mapping elapsed time", timedelta(seconds = TIMES.end - TIMES.start)

for proc in TEMP.assign_procs:
    proc.wait()

print "Assignment extra time", timedelta(seconds = time() - TIMES.end)

#TIMES.sortstart = time()
#DATA.fs = glob(join(ARGS.analysis_dir, '*', 'assigned_dmel.bam'))
#TEMP.sort = Popen(['parallel',
              #'samtools sort {} -m %d {//}/dmel_sorted' %3e9, # 3GB of memory
              #':::'] + DATA.fs)
#
#TEMP.sort.wait()
#
#TIMES.sortend = time()
#print "Sorting time", timedelta(seconds = TIMES.sortend - TIMES.sortstart)
#
## Figure out how well everything mapped
#DATA.mapped_reads = {}
#DATA.all_bams = [join(ARGS.analysis_dir, sample_dir, 'dmel_sorted.bam')
            #for sample_dir in DATA.samples]
#
#for sample, bam in zip(DATA.samples, DATA.all_bams):
    #print '='*30
    #commandstr = ['samtools', 'flagstat', bam]
    #print commandstr
    #samtools_proc = Popen(commandstr, stdout=PIPE)
    #samout, samerr = samtools_proc.communicate()
#
    #for line in str(samout).splitlines():
        #if "mapped" in line:
            #DATA.mapped_reads[sample] = float(line.split()[0])
            #print "% reads dmel in ", sample,
            #print DATA.mapped_reads[sample] / DATA.num_reads[sample] * 100
            #break
#
#
#
#for sample in DATA.samples:
    #TEMP.od = join(ARGS.analysis_dir, sample)
    #TEMP.commandstr = (BASE.cufflinks_base +
                  #'-G %(GTF)s -o %(od)s %(hits)s'
                  #% dict(GTF=ARGS.base_GTF,
                         #od=TEMP.od,
                         #hits=join(TEMP.od, 'dmel_sorted.bam')))
#
    #print TEMP.commandstr
    #sys.stdout.flush()
    #TEMP.cufflinks_proc = Popen(str(TEMP.commandstr).split())
#
    #TEMP.cufflinks_proc.wait()
#
print "Final time", timedelta(seconds=time() - TIMES.start)
#print "Cufflinks time", timedelta(seconds=time() - TIMES.sortend)

import cPickle as pickle

pickle.dump(dict(data=DATA, args=ARGS),
            open('mapreads_dump.pkl', 'w'))

#for proc in TEMP.rezip_procs:
    #proc.wait()
