from Bio import SeqIO
from os import listdir, chdir
from sys import argv

chdir(argv[1])
out_fh = open('multi.fa', 'w')

for fname in listdir('.'):
    if not fname.startswith('d'): continue
    specname = fname.split('-')[0]
    print specname

    for rec in SeqIO.parse(fname, 'fasta'):
        rec.id = specname + '_' + rec.id
        SeqIO.write(rec, out_fh, 'fasta')

