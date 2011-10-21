import sys
import string
import cPickle as pickle

from glob import glob
from os import system, popen
from os.path import join
from time import time

analysis_dir = 'analysis'
GTF =  'Reference/dmel-all-r5.32_transcripts.gtf'
idxfile = 'Reference/dmel-all-r5.23'
interest = 'GenesOfInterest.txt'
FBtoName = 'Reference/dmelfbgns.txt'
notificationEmail = 'peter.combs@berkeley.edu'

########################################################################

tophat_base = 'tophat -p8 -r 200 --no-novel-juncs '
cufflinks_base = 'cufflinks -p 8 -q '
cuffdiff_base = ('cufflinks.cuffdiff -p 8 -v --FDR .001 -o %(ad)s %(gtf)s '
                 % {'gtf':GTF, 'ad': analysis_dir})


########################################################################

reads = ['s_5_1_sequence.txt s_5_2_sequence.txt', 
         's_6_1_sequence.txt s_6_2_sequence.txt']

libraries = {'s_5_1_sequence.txt' : 'PAC03',
             's_6_1_sequence.txt' : 'PAC02'}

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

for rf in reads:
    # Print the name of the files we're going through, as a rough progress bar
    print '-'*72
    print rf
    print '-'*72

    # Just grab the first file name (paired ends have the same number in both)
    rf2 = rf.split()[0]

    wcout = popen('wc -l ' + rf2, 'r', 1024).readline()
    print wcout

    # Store the number of reads in the file
    numreads[rf2] = int(wcout.split()[0])
    assert numreads[rf2]%4 == 0
    numreads[rf2] /= 4
    
    od = join(analysis_dir, rf2.split('.')[0])

    # Figure out Read Group ID
    f = open(rf.split()[0])
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
              'library': libraries[rf.split()[0]],
              'rgid': rgid,
              'lane': lane})
    print commandstr
    sys.stdout.flush()
    res = system(commandstr)


    if res:
        system('echo "Oh no!" | mail -s "Failed on Tophatting %s" %s'
               % (rf, notificationEmail))

    # Do cufflinks

    print 'Cufflinksing...', '\n', '='*30
    sys.stdout.flush()
    res = system(cufflinks_base + '-G %(GTF)s -o %(od)s %(hits)s' 
           % {'GTF': GTF, 'od': od,
              'hits': join(od, 'accepted_hits.bam')})

    if res:
        system('echo "Oh no!" | mail -s "Failed on Cufflinksing %s" %s'
               % (rf, notificationEmail))


    # Figure out how well everything mapped
    print '='*30
    pipe = popen('samtools flagstat ' + join(od, 'accepted_hits.bam'), 'r',
                 1024)

    for line in pipe:
        if "mapped" in line:
            mappedreads[rf2] = int(line.split()[0])
            break

all_bams = map(lambda s: join('analysis', s, 'accepted_hits.bam'), 
               (s.split('.fq')[0] for s in reads))


# Do Cuffdiff
system(cuffdiff_base + " ".join(all_bams))

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
except:
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
system('cat to_email.tmp | mail -s "Done" ' + notificationEmail )


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
    system('mutt -s "LogLog" %s -a LogLog.pdf < to_email.tmp'% notificationEmail)

except ImportError:
    print "Could not import Matplotlib.  You are using Python version",
    print sys.version
