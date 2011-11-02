import sys
import cPickle as pickle

from glob import glob
from os.path import join
from time import time
from subprocess import Popen, PIPE

analysis_dir = 'analysis'
GTF =  'Reference/dmel-all-r5.32_transcripts_fixed.gtf'
idxfile = 'Reference/dmel-all-r5.23'
interest = 'GenesOfInterest.txt'
FBtoName = 'Reference/dmelfbgns.txt'
notificationEmail = 'peter.combs@berkeley.edu'
seq_dir = 'sequence'

########################################################################

tophat_base = 'tophat -p8 -r 200 --no-novel-juncs '
cufflinks_base = 'cufflinks -p 8 -q '
cuffdiff_base = ('cufflinks.cuffdiff -p 8 -v --FDR .001 -o %(ad)s %(gtf)s '
                 % {'gtf':GTF, 'ad': analysis_dir})


########################################################################

reads = [','.join(glob(join(seq_dir, '*index2*'))),
         ','.join(glob(join(seq_dir, '*index4*'))),
         ','.join(glob(join(seq_dir, '*index5*'))),
         ','.join(glob(join(seq_dir, '*index6*')))
        ]

readnames = {'index2' : ','.join(glob(join(seq_dir, '*index2*'))),
             'index4' : ','.join(glob(join(seq_dir, '*index4*'))),
             'index5' : ','.join(glob(join(seq_dir, '*index5*'))),
             'index6' : ','.join(glob(join(seq_dir, '*index6*')))
            }
        

libraries = {
            'index2' : 'A',
            'index4' : 'B',
            'index5' : 'C',
            'index6' : 'D',
            }

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

for readname, rf in readnames.iteritems():
    # Print the name of the files we're going through, as a rough progress bar
    print '-'*72
    print rf
    print '-'*72

    # Just grab the first file name (paired ends have the same number in both)
    rf2 = rf.split(',')

    # This section will probably need to be fixed if/when I do paired-end reads.
    # Just splitting on commas will have a non-existant file with the last file
    # of the first end and the first file of the second end
    wc_proc = Popen(['wc', '-l']+ rf2, stdout=PIPE)
    wcout, wcerr = wc_proc.communicate()
    print wcout

    # Store the number of reads in the file
    numreads[readname] = int(wcout.splitlines()[-1].split()[0])
    assert numreads[readname] % 4 == 0
    numreads[readname] /= 4

    od = join(analysis_dir, readname)

    # Figure out Read Group ID
    f = open(rf.split(',')[0])
    l = f.readline()
    f.close()
    rgid = l.split(":")[0][1:]
    lane = l.split(":")[1]


    # Do tophat
    print 'Tophatting...', '\n', '='*30
    commandstr =  (tophat_base + '-G %(GTF)s -o %(od)s --rg-library %(library)s'
                   ' --rg-center VCGSL --rg-sample %(library)s --rg-platform'
                   ' ILLUMINA --rg-id %(rgid)s  --rg-platform-unit %(lane)s %(idxfile)s %(rf)s'
           % {'GTF': GTF,
              'od': od,
              'idxfile': idxfile,
              'rf': rf,
              'library': readname,
              'rgid': rgid,
              'lane': lane})
    print commandstr
    sys.stdout.flush()
    tophat_proc = Popen(commandstr.split())
    tophat_proc.wait()


    if tophat_proc.returncode:
        errormail_proc = Popen(['mail', '-s', "Failed on tophatting %s" %
                                readname,
                                notificationEmail], stdin=PIPE)
        errormail_proc.communicate('Oh no!')

    # Do cufflinks

    print 'Cufflinksing...', '\n', '='*30
    sys.stdout.flush()
    commandstr = (cufflinks_base + '-G %(GTF)s -o %(od)s %(hits)s'
           % {'GTF': GTF, 'od': od,
              'hits': join(od, 'accepted_hits.bam')})
    print commandstr
    cufflinks_proc = Popen(commandstr.split())

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
    p2 = Popen(['samtools', 'rmdup', join(od, 'accepted_hits.bam'),
                join(od, 'filtered_hits.bam'),],
               stdout=file(join(od, 'hit_filtering.log'), 'w'))
    p2.wait()

all_bams = map(lambda s: join('analysis', s, 'accepted_hits.bam'),
               (s for s in readnames))


# Do Cuffdiff
#system(cuffdiff_base + " ".join(all_bams))
cuffdiff_proc = Popen(cuffdiff_base.split() +
                      ['-L', ','.join(libraries[rf] for rf in readnames)]
                      + all_bams)

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
