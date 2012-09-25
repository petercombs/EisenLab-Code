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

analysis_dir = 'analysis-multi'
GTF =  'Reference/dmel.gtf'
idxbase = 'Reference/AAA/'
base_species = 'mel'
notificationEmail = 'peter.combs@berkeley.edu'
seq_dir = 'sequence'
config_file = 'RunConfig.cfg'

########################################################################

tophat_base = ('tophat -p8 --no-novel-juncs --read-mismatches 4 '
                '--report-secondary-alignments ')
cufflinks_base = 'cufflinks -p 8 -q -u '


########################################################################



def process_config_file(cfg_fname):
    """ Get data we want out of the configuration file"""
    cfg_fh = open(cfg_fname)
    cfg_data = dict(sample_to_lib = [], samples=[], libraries = [],
                    sample_to_carrier = {})
    for line in cfg_fh:
        try:
            line = line.strip().split('\t')
            lib, mbepc, slice, idx, carrier = line
            cfg_data['sample_to_lib'].append((lib + slice, mbepc+slice))
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
        read_1s = glob(join(seq_dir, 'Sample*' + libname + '*', "*_R1_*.fastq"))
        read_2s = glob(join(seq_dir, 'Sample*' + libname + '*', "*_R2_*.fastq"))
        readnames[sample] = [','.join(read_1s), ','.join(read_2s)]
    return readnames

def count_reads(read_files):
    ''' Count the reads in the files'''

    # Let's cheat and only look in the last file: the others should have 4M
    # per...
    wc_proc = Popen(['wc', '-l', read_files[-1]], stdout=PIPE)
    wcout, wcerr = wc_proc.communicate()
    first_reads = int(4e6) * (len(read_files) - 1)
    last_reads = int(wcout.splitelines()[-1].strip()) / 4
    return first_reads + last_reads

config_data = process_config_file(config_file)
readnames = get_readfiles(config_data)

samples = config_data['samples']


assign_procs = []
num_reads = {}

start = time()

for libname, (rf1, rf2) in readnames.items():
    # Print the name of the files we're going through, as a progress bar
    print '-'*72, '\n', libname, '\n', '-'*72

    # Just grab the first file name (PE have the same number in both)
    rfs = rf1.split(',')
    num_reads[libname] = count_reads(rfs)

    od = join(analysis_dir, libname)
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

    idxfile = join(idxbase, base_species +
                   config_data['sample_to_carrier'][libname])

    # Do tophat
    print 'Tophatting...', '\n', '='*30
    commandstr =  (tophat_base + '-G %(GTF)s -o %(od)s --rg-library '
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
    assign_procs.append(Popen(commandstr))




# Stop the timing
end = time()
print "Mapping elapsed time", timedelta(seconds = end-start)

for proc in assign_procs:
    proc.wait()

print "Assignment extra time", timedelta(seconds = time() - end)

sortstart = time()
fs = glob(join(analysis_dir, '*', 'assigned_dmel.bam'))
sort = Popen(['parallel',
              'samtools sort {} -m %d {//}/dmel_sorted' %3e9, # 3GB of memory
              ':::'] + fs)

sort.wait()

sortend = time()
print "Sorting time", timedelta(seconds = sortend - sortstart)

# Figure out how well everything mapped
mapped_reads = {}
all_bams = [join(analysis_dir, sample_dir, 'dmel_sorted.bam')
            for sample_dir in samples]

for sample, bam in zip(samples, all_bams):
    print '='*30
    commandstr = ['samtools', 'flagstat',
                  join(analysis_dir, sample, 'dmel_sorted.bam')]
    samtools_proc = Popen(commandstr, stdout=PIPE)
    samout, samerr = samtools_proc.communicate()

    for line in samout.splitlines():
        if "mapped" in line:
            mapped_reads[sample] = int(line.split()[0])
            print "% reads dmel in ", sample,
            print mapped_reads[sample]/num_reads[sample]
            break



for sample in samples:
    od = join(analysis_dir, libname)
    commandstr = (cufflinks_base + '-G %(GTF)s -o %(od)s %(hits)s'
                  % dict(GTF=GTF, od=od, hits=join(od, 'dmel_sorted.bam')))

    print commandstr
    cufflinks_proc = Popen(commandstr.split())

    cufflinks_proc.wait()

