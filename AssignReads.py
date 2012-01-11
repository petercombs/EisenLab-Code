from Bio import SeqIO
from sys import stdin
from collections import Counter, defaultdict
import inspect_shell

parser = SeqIO.parse('Reference/melpsevir-transcript.fasta', 'fasta')
specs = {}
for rec in parser:
    pos = rec.description.find('species') + len('species') + 1
    specs[rec.id] = rec.description[pos:pos + 4]

counts = Counter()

for line in stdin:
    data = line.split()
    if int(data[1]) & 0x100:
        continue
    counts[specs[data[2]]] += 1

print counts


                
