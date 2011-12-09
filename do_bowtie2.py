import sys
import cPickle as pickle

from glob import glob
from os.path import join
import os
from time import time
from subprocess import Popen, PIPE

indices_used = [2,4,5,6]
analysis_dir = 'analysis-express'
GTF =  'Reference/dmel-all-r5.32_transcripts_fixed.gtf'
bowtie_index = 'Reference/dmel-all-transcript-r5.42'
target_seqs = 'Reference/dmel-all-transcript-r5.42.fasta'
interest = 'GenesOfInterest.txt'
FBtoName = 'Reference/dmelfbgns.txt'
notificationEmail = 'peter.combs@berkeley.edu'
seq_dir = 'sequence'

########################################################################


bowtie_options = ['bowtie2', '-a', '-p', '8', '--very-sensitive-local', '-t', ]

########################################################################


readnames = {"index%d" % idx: 
                ",".join(sorted( glob(join(seq_dir, '*_index%d_*' % idx))))
             for idx in indices_used }

libraries = { "index%d" % idx : chr(ord('A') + i )
             for i, idx in enumerate(indices_used)}

# Dictionary with the number of reads in each file
numreads = {}

# Dictionary with the number of mapped reads.
mappedreads = {}


# Interconversion from FlyBase IDs to Gene Names
FBKey = {}
NameKey = {}

for line in file(FBtoName):
    line = line.split()
    FBKey[line[0]] = line[1]
    NameKey[line[1]] = line[0]


start = time()
if '-cdo' not in sys.argv:
    for readname, rf in sorted(readnames.items()):
        # Print the name of the files we're going through, as a rough progress bar
        print '-'*72
        print rf
        print '-'*72

        od = join(analysis_dir, readname)
        try:
            os.makedirs(od)
        except OSError:
            print "Directory '%s' already exists... shouldn't be a problem" % od

        # Figure out Read Group ID
        f = open(rf.split(',')[0])
        l = f.readline()
        f.close()
        rgid = l.split(":")[0][1:]
        lane = l.split(":")[1]


        # Do Bowtie
        commandstr = bowtie_options + [ '-x', bowtie_index, '-U', rf, 
                                       '-S', join(od, 'accepted_hits.sam')]
        print ' '.join(commandstr)
        sys.stdout.flush()
        bowtie_proc = Popen(commandstr)

        bowtie_proc.wait()

        print '='*30

# Stop the timing
end = time()
hours = int((end - start)/3600)
minutes = int((end - start) / 60 - hours * 60)
seconds = int((end - start) - minutes * 60 - hours * 3600)
print hours, "h", minutes, "m", seconds, "s"

email = open('to_email.tmp', 'w')
email.write("%dh:%dm:%ds\n"%(hours, minutes, seconds))

email.write("\n\n")


# Dump everything out to a file, so we can play with it later, maybe
try:
    pickle.dump(dict([(k,v) for k,v in locals().copy().iteritems()
                  if ((type(v) is not type(sys))
                     and (type(v) is not file))]),
                file('tuxedo_dump', 'w'))
except Exception as exc:
    print exc
    print "the pickling still doesn't work... skipping"




# Send an email to the user
email.close()
errormail_proc = Popen(['mail', '-s', "Done", notificationEmail],
                       stdin=open('to_email.tmp'))

