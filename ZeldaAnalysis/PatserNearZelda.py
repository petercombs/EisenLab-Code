import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

from Bio import SeqIO
from os import makedirs
from subprocess import Popen, PIPE
from collections import Counter
from glob import glob
from progressbar import ProgressBar
import PeakFinder

N = 200
FA = 'fasta'

if 'dmel' not in locals():
    dmel = {r.name: r for r in
            SeqIO.parse('prereqs/dmel-all-chromosome-r5.52.fasta', FA)}

if 'zld_bind' not in locals():
    zld_bind = pd.read_table('prereqs/journal.pgen.1002266.s005.xls', skiprows=1)

try:
    makedirs('analysis/ZldSeqs')
except OSError:
    pass

base_rate = Counter()

seqs = []
outfh = open('analysis/ZldSeqs/forpatser2.txt', 'w')
for i, site in zld_bind.iterrows():
    chr = site["Chr"].replace('chr', '')
    pos = site["Peak"]
    pid = site["PID"]
    subseq =  dmel[chr][pos-N:pos+N]
    subseq.id = "zld{:05d}".format(pid)
    subseq.description = 'zld{pid:05d} {chr}:{start}-{stop}'.format(
        pid = pid,
        chr = chr,
        start = pos - N,
        stop = pos + N,
    )
    if 'N' not in subseq:
        seqs.append(subseq.id)
        outfh.write('{} \\{}\\\n'.format(subseq.id, subseq.seq))
    base_rate += Counter(subseq.seq)
    out_fname = 'analysis/ZldSeqs/patser{pid:05d}.txt'.format(pid=pid)
    #SeqIO.write(subseq, out_fname, FA)
    outfh2 = open(out_fname, 'w')
    outfh2.write('{} \\{}\\\n'.format(subseq.id, subseq.seq))
    outfh2.close()

outfh.close()
at_rate = base_rate['A'] + base_rate['T']
gc_rate = base_rate['G'] + base_rate['C']


TFs = ('bcd cad D dl ftz gt h hb hkb kni kr mad '
       'prd run shn slp1 sna tll twi').split()
TFs = 'bcd cad D ftz gt hb hkb kni kr prd run slp1 tll'.split()
TFs = 'bcd'.split()
all_binds = pd.DataFrame(index=seqs, data=np.zeros((len(seqs), 2*N), dtype=int))

prog = ProgressBar()

overlaps = []
for i, TF in enumerate(prog(TFs)):
    #print "\n".join(["-"*20,TF, "-"*20])
    TF_file = 'prereqs/BDTNP_pwms/{}.cm'.format(TF)
    patser = ['patser',
              '-m', TF_file,
              '-f', outfh.name,
              '-A', 'a:t {} c:g {}'.format(at_rate, gc_rate),
              '-li',
              '-c',
              #'-t', '5',
              #'-ds'
             ]
    proc = Popen(patser, stdout=PIPE)
    for line in proc.stdout:
        data = line.split()
        if not data or not data[0].startswith('zld'): continue
        name = data[0]
        if name not in all_binds.index: continue
        if data[2].endswith('C'):
            data[2] = data[2][:-1]
        pos = int(data[2])
        if all_binds.ix[name, pos]:
            overlaps.append((name, pos, TF, all_binds.ix[name, pos]))
        all_binds.ix[name,pos] = i + 1



zld_bind_TSS = pd.read_table('prereqs/journal.pgen.1002266.s005.xls', skiprows=1, index_col='TSS_gene')
zld_bind_TSS.index = map(str.strip, zld_bind_TSS.index)

wt_exp = pd.read_table('prereqs/WT5.53_summary.tsv', index_col=0)
zld_exp  = pd.read_table('analysis/summary.tsv', index_col=0)

startswith = lambda y: lambda x: x.startswith(y)

z13_1 = zld_exp.select(startswith('cyc13_rep1'), axis=1)
z13_2 = zld_exp.select(startswith('cyc13_rep2'), axis=1)
z13_3 = zld_exp.select(startswith('cyc13_rep3'), axis=1)
w13 = wt_exp.select(startswith('cyc13'), axis=1)

sortnum = PeakFinder.make_sort_num((z13_1, z13_2, z13_3), w13)
sortnum.sort()
sortnum = sortnum[::-1]


PID_list = []
for gene in sortnum.index:
    if sortnum[gene] == 0: break
    if gene not in zld_bind_TSS.index: continue
    gene_bind = zld_bind_TSS.ix[gene]
    if len(np.shape(gene_bind)) == 2:
        PID_list.extend(gene_bind.PID)
    else:
        PID_list.append(gene_bind["PID"])

zpid_list = ['zld{:05d}'.format(int(pid)) for pid in PID_list]

all_binds.ix[zpid_list]

all_binds_df = pd.DataFrame(all_binds)
plt.pcolormesh(np.array(all_binds_df.ix[np.array(PID_list[:1000]) + 1]))
