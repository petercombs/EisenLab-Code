from Bio import SeqIO
from os import listdir, chdir
from sys import argv

chdir(argv[1])
base = argv[2]
rest = argv[3:]

print "Base species: ", base
print "Other species: ", rest

for other in rest:
    print base+other
    out_fh = open(base + other + '.fa', 'w')

    for fname in listdir('.'):
        if fname.startswith('d') and '.fa' in fname and ((base in fname) or (other in fname)):
            specname = fname.split('-')[0]
            print specname

            for rec in SeqIO.parse(fname, 'fasta'):
                rec.id = specname + '_' + rec.id
                SeqIO.write(rec, out_fh, 'fasta')
    out_fh.close()


for other in rest:
    print "GTF", base+other
    out_fh = open(base + other + '.gtf', 'w')

    for fname in listdir('.'):
        if fname.startswith('d') and '.gtf' in fname and ((base in fname) or (other in fname)):
            specname = fname.split('-')[0]
            print specname

            for line in open(fname):
                out_fh.write(specname+"_"+line)
    out_fh.close()


