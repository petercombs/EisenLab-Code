blast_records
blast_recs = [r for r in blast_records]
len(blast_recs)
blast_recs[0]
dir(blast_recs[0])
dir(blast_recs[0].alignments[0]
)
blast_recs[0].alignments[0].hit_def
blast_recs[0].alignments[0].hit_id
blast_recs[0].alignments[0].hsps
dir(blast_recs[0].alignments[0].hsps[0])
blast_recs[0].alignments[0].hsps[0].sbjct
blast_recs[0].alignments[0].hsps[0].query
blast_recs[0].alignments[0].hit_def
[r.alignments[0].hit_def.split()[:2] for r in blast_recs]
import collections as cs
cs.Counter([r.alignments[0].hit_def.split()[:2] for r in blast_recs])
cs.Counter(set([r.alignments[0].hit_def.split()[:2] for r in blast_recs]))
#?cs.Counter
cs.Counter(iterable=[r.alignments[0].hit_def.split()[:2] for r in blast_recs])
cs.Counter(iterable=[set(r.alignments[0].hit_def.split()[:2]) for r in blast_recs])
fcs.Counter(iterable=[tuple(r.alignments[0].hit_def.split()[:2]) for r in blast_recs])
fcs.Counter([tuple(r.alignments[0].hit_def.split()[:2]) for r in blast_recs])
[tuple(r.alignments[0].hit_def.split()[:2]) for r in blast_recs]
[tuple(r.alignments[0].hit_def.split()) for r in blast_recs]
cs.Counter([tuple(r.alignments[0].hit_def.split()) for r in blast_recs])
cs.Counter([tuple(r.alignments[0].hit_def) for r in blast_recs])
cs.Counter([r.alignments[0].hit_def for r in blast_recs])
cs.Counter([tuple(r.alignments[0].hit_def.split()[:2]) for r in blast_recs])
blast_recs6 = [r for r in NCBIXML.parse(open('5.blastout.xml'))]
blast_recs6 = [r for r in NCBIXML.parse(open('6.blastout.xml'))]
c5 = cs.Counter([tuple(r.alignments[0].hit_def.split()[:2]) for r in blast_recs])
sum(c5)
len(c5)
len(blast_recs6)
c6 = cs.Counter([tuple(r.alignments[0].hit_def.split()[:2]) for r in blast_recs6])
c6
#?sum
help(history)
#?history
_ip.magic("history -n -f GetSpeciesFromBlast.py")
