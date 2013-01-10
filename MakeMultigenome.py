from Bio import SeqIO
from os import listdir, chdir
from sys import argv
from subprocess import Popen

chdir(argv[1])
base = argv[2]
rest = argv[3:]

print "Base species: ", base
print "Other species: ", rest

## Make the BASE_only.fa file
only_out_fh = open(base+"_only.fa", "w")
for fname in listdir('.'):
    if fname.startswith('d') and ('.fa' in fname) and (base in fname):
        specname = fname.split('-')[0]
        print specname

        for rec in SeqIO.parse(fname, 'fasta'):
            rec.id = specname + '_' + rec.id
            SeqIO.write(rec, only_out_fh, 'fasta')

only_out_fh.close()

bowtie_build_procs = []

## Make the joint fasta files with each species
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
    bowtie_build_procs.append(Popen(['bowtie2-build', out_fh.name, base+other]))

## Make the BASE_only gtf file
only_out_fh = open(base + '_only.gtf', 'w')
for fname in listdir('.'):
    if fname.startswith('d') and ('.gtf' in fname) and (base in fname):
        specname = fname.split('-')[0]
        print specname

        for line in open(fname):
            only_out_fh.write(specname+"_"+line)
only_out_fh.close()

## Make the joint gtf files with each species
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


for proc in bowtie_build_procs:
    proc.wait()
