import os
from Bio import SeqIO
from os import path
from glob import glob

def parse_description(record):
    """ Turn the description into a dictionary

    keys are the left hand side of the equals
    values are the right and side of the equals
    separate items are separated by spaces
    semicolons are trimmed from descriptions
    """
    description = {}
    for item in record.description.split(' '):
        if '=' not in item: continue
        item = item.strip('; ').split('=')
        description[item[0]] = item[1]
    return description



data_dir = '/Volumes/Plasmid/FlyBase/12genomes/'
ortholog_table = path.join(data_dir, 'gene_orthologs_fb_2011_09.tsv')
output_dir = path.expanduser('~/data/Orthologs/')
if not path.isdir(output_dir):
    os.makedirs(output_dir)

other2mel = {}

for line in open(ortholog_table):
    # Skip comments and empty lines
    if line.startswith("#"): continue
    if not line.strip(): continue

    line = line.split()

    # Magic numbers taken from
    # http://flybase.org/static_pages/docs/datafiles.html
    mel = line[0]
    other = line[5]
    other2mel[other] = mel
    # We also want the mel names to map to mel
    other2mel[mel] = mel

for fname in glob(path.join(data_dir, '*.fasta')):
    print fname
    # dmel => mel
    for i, record in enumerate(SeqIO.parse(fname, 'fasta')):
        description = parse_description(record)
        try:
            melequiv = other2mel[description['parent']]
        except:
            if description['species'] == 'Dmel':
                melequiv = description['parent']
            else:
                continue
        record.id = "%s_%d" % (description['species'], i)

        outfile = open(path.join(output_dir, melequiv + '.fasta'), 'a')
        SeqIO.write(record, outfile, 'fasta')
        outfile.close()

