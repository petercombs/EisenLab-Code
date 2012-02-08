import sys
import cPickle as pickle

from glob import glob
from os.path import join
import os
from time import time
from subprocess import Popen, PIPE

analysis_dir = 'analysis42'
GTF =  'Reference/dmel-all-r5.42.gtf'
idxfile = 'Reference/dmel-all-chromosome-r5.42'
interest = 'GenesOfInterest.txt'
FBtoName = 'Reference/dmelfbgns.txt'
notificationEmail = 'peter.combs@berkeley.edu'
seq_dir = 'sequence'

########################################################################

tophat_base = ['tophat', 
               '--num-threads', '4', 
               '--no-novel-juncs',
               '--transcriptome-index', 'FlyBaseTrancriptome/known',
               '--transcriptome-only']
cufflinks_base = ['cufflinks', 
                  '--num-threads',  '4', 
                  '--quiet', '--multi-read-correct', 
                  '--frag-bias-correct', idxfile + '.fa ', 
                  '-M', mask_GTF] #Mask file
cuffdiff_base = ('cuffdiff -p 8 -v -o %(ad)s %(gtf)s '
                 % {'gtf':GTF, 'ad': analysis_dir})


########################################################################


indices_used = set([name.split('_')[1][5:] 
                    for name in glob(join(seq_dir, '*index*'))])
readnames = {"index%d" % idx: [",".join(sorted( glob(join(seq_dir,
                                                          '*_index%d_*_R1*'
                                                          % idx)))),
                               ",".join(sorted( glob(join(seq_dir,
                                                          '*_index%d_*_R2*' % idx))))
                              ]
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

print "Looking at indices:", indices_used

start = time()
if '-cdo' not in sys.argv:
    for readname, (rf1, rf2) in sorted(readnames.items()):
        # Print the name of the files we're going through, as a rough progress bar
        print '-'*72
        print readname
        print '-'*72
        print rf1, rf2
        print '-'*32
        sys.stdout.flush()

        # Just grab the first file name (paired ends have the same number in both)
        rfs = rf1.split(',')

        # This section will probably need to be fixed if/when I do paired-end reads.
        # Just splitting on commas will have a non-existant file with the last file
        # of the first end and the first file of the second end
        wc_proc = Popen(['wc', '-l']+ rfs, stdout=PIPE)
        wcout, wcerr = wc_proc.communicate()
        print wcout

        # Store the number of reads in the file
        numreads[readname] = int(wcout.splitlines()[-1].split()[0])
        assert numreads[readname] % 4 == 0
        numreads[readname] /= 4

        od = join(analysis_dir, readname)
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


        # Do tophat
        print 'Tophatting...', '\n', '='*30
        commandstr =  tophat_base + ['--GTF', GTF, 
                                     '--output-dir', od,
                                     '--rg-library', readname, 
                                     '--rg-center', 'VCGSL', 
                                     '--rg-sample', readname,
                                     '--rg-platform', 'ILLUMINA', 
                                     '--rg-id', rgid, 
                                     '--rg-platform-unit', lane,
                                     idxfile, rf1, rf2]
        print " ".join(commandstr)
        sys.stdout.flush()
        sys.stderr.flush()
        tophat_proc = Popen(commandstr)
        tophat_proc.wait()


        if tophat_proc.returncode:
            errormail_proc = Popen(['mail', '-s', "Failed on tophatting %s" %
                                    readname,
                                    notificationEmail], stdin=PIPE)
            errormail_proc.communicate('Oh no!')

        # Do cufflinks

        print 'Cufflinksing...', '\n', '='*30
        sys.stdout.flush()
        sys.stderr.flush()
        commandstr = cufflinks_base + ['--GTF', GTF, '--output-dir', od,
                                       join(od, 'accepted_hits.bam')]
        print " ".join(commandstr)
        cufflinks_proc = Popen(commandstr)

        cufflinks_proc.wait()
        if cufflinks_proc.returncode:
            errormail_proc = Popen(['mail', '-s', "Failed on cufflinksing %s" %
                                    readname,
                                    notificationEmail], stdin=PIPE)
            errormail_proc.communicate('Oh no!')


        # Figure out how well everything mapped
        print '='*30
        commandstr = ['samtools', 'flagstat', join(od, 'accepted_hits.bam')]
        samtools_proc = Popen(commandstr, stdout=PIPE)
        samout, samerr = samtools_proc.communicate()

        for line in samout.splitlines():
            if "mapped" in line:
                mappedreads[readname] = int(line.split()[0])
                break
        p2 = Popen(['samtools', 'rmdup', '-s', join(od, 'accepted_hits.bam'),
                    join(od, 'filtered_hits.bam'),],
                   stdout=file(join(od, 'hit_filtering.log'), 'w'))
        p2.wait()

all_bams = map(lambda s: join('analysis', s, 'accepted_hits.bam'),
               (s for s in readnames))


# Do Cuffdiff
#system(cuffdiff_base + " ".join(all_bams))
cuffdiff_call = (cuffdiff_base.split() 
                 + ['-L', ','.join(libraries[rf] for rf in sorted(readnames.keys()))]
                 + all_bams)

print ' '.join(cuffdiff_call)

cuffdiff_proc = Popen(cuffdiff_call)
cuffdiff_proc.wait()

# Stop the timing
end = time()
hours = int((end - start)/3600)
minutes = int((end - start) / 60 - hours * 60)
seconds = int((end - start) - minutes * 60 - hours * 3600)
print hours, "h", minutes, "m", seconds, "s"

# GOTerm Finder stuff
email = open('to_email.tmp', 'w')
email.write("%dh:%dm:%ds\n"%(hours, minutes, seconds))
golink = "http://go.princeton.edu/cgi-bin/GOTermFinder" \
        +"?geneAsscFile=gene_association.fb"\
        +"&email1=%(email)s&email2=%(email)s" % {"email": notificationEmail}
email.write(golink+"&aspect=F\n")
email.write(golink+"&aspect=P\n")
email.write(golink+"&aspect=C\n")

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


# Print out the actual read mapping percentages

for rf in numreads:
    print rf, numreads[rf], 100.0*mappedreads[rf]/numreads[rf]


# Write out to the email the list of significant genes

genediff = {}
for line in file(join(analysis_dir,'gene_exp.diff')):
    genediff[line.split()[0]] = line.strip()
    if "yes" in line:
        email.write(line.split()[0] + "\n")


# Send an email to the user
email.close()
errormail_proc = Popen(['mail', '-s', "Done", notificationEmail],
                       stdin=open('to_email.tmp'))


# If we're on a system that supports it, do the plotting

try:
    from numpy import array
    from matplotlib import pyplot as mpl

    s1, s2, gene, idx = zip(*[(float(line.split()[6]), float(line.split()[7]),
                              line.split()[0], lnum)
                              for lnum, line
                              in enumerate(file(join(analysis_dir, 'gene_exp.diff')))
                             if 'FBgn' in line])

    s1 = array(s1)
    s2 = array(s2)

    FBgnToIDX = dict(zip(gene, idx))

    mpl.loglog(s1,s2,'k.', label="All Genes")
    for fname in glob('Genes*.txt'):
        label = fname[5:-4]
        genes_of_interest = [l.strip() for l in file(fname)]
        s1i = s1[[FBgnToIDX[NameKey[gene]] for gene in genes_of_interest
                 if gene in NameKey and NameKey[gene] in FBgnToIDX]]
        s2i = s2[[FBgnToIDX[NameKey[gene]] for gene in genes_of_interest
                 if gene in NameKey and NameKey[gene] in FBgnToIDX]]
        mpl.loglog(s1i, s2i, '.', label=label)
        print "-"*72, "\n", fname, "\n", "-"*72
        for g in genes_of_interest:
            if g.strip() in NameKey and NameKey[g.strip()] in genediff:
                print g, genediff[NameKey[g.strip()]]

    mpl.legend(numpoints=1, loc='lower right')
    mpl.loglog([1e-2,2e4], [1e-2,2e4], 'r:') # Diagonal line to guide eye


    # Clean up and label the axes
    ax = mpl.gca()
    ax.set_xlim(1e-2,1e+4)
    ax.set_ylim(1e-2,1e+4)
    ax.set_xlabel('Ant. Expr (RPKM)')
    ax.set_ylabel('Pos. Expr (RPKM)')
    mpl.savefig('LogLog.pdf')

    # If we have mutt, attach file to an email and send it off
    graphmail_proc = Popen(['mutt', '-s', 'LogLog', notificationEmail, '-a',
                            'LogLog.pdf'], stdin=open('to_email.tmp'))

except ImportError:
    print "Could not import Matplotlib.  You are using Python version",
    print sys.version
