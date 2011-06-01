from __future__ import print_function

import collections as cs
from Bio.Blast import NCBIXML
import os

os.chdir('blaststuff')
blast_recs5 = [r for r in NCBIXML.parse(open('5.blastout.xml'))]
blast_recs6 = [r for r in NCBIXML.parse(open('6.blastout.xml'))]
c5 = cs.Counter([tuple(r.alignments[0].hit_def.split()[:2]) for r in blast_recs5])
c6 = cs.Counter([tuple(r.alignments[0].hit_def.split()[:2]) for r in blast_recs6])


print(c5)
print(c6)
