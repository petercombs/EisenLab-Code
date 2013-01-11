import pandas as pd
from collections import defaultdict
from os import path
from subprocess import Popen

fasta_ref = 'Reference/AAA/mel_only.fa',
gtf_ref = 'Reference/AAA/mel_only.gtf'

Popen(['cuffcompare',
       '-o', path.join(path.dirname(gtf_ref), 'cuffcmp.'),
       '-s', fasta_ref,
       '-CG', '-r', gtf_ref, gtf_ref]).wait()

gtf_ref = path.join(path.dirname(gtf_ref), 'cuffcmp.') + 'combined.gtf'

design_file = pd.read_table('analysis-multi/design.tab')
files = defaultdict(list)
for ix, row in design_file.iterrows():
    files[row['condition']].append(path.join('analysis-multi',
                                          row['Sample'],
                                          'assigned_dmel_rescued.bam')
                                )
conditions = sorted(files.keys())

cd_base = 'cuffdiff -p 8 -L {conditions} -o {outdir} -u -N -b {fasta} {gtf} {bams}'
cd = cd_base.format(conditions = ','.join(conditions),
                    outdir = 'analysis',
                    fasta = fasta_ref,
                    gtf = gtf_ref,
                    bams = ' '.join([','.join(files[key]) for key in conditions]))

print cd
Popen(cd).wait()
cl_base = 'cufflinks -p 8 -o {outdir} -u -N -b {fasta} -G {gtf} {bamfile}'
for condition in files:
    for file in files[condition]:
        cl = cl_base.format(outdir = path.dirname(file),
                           fasta = fasta_ref,
                           gtf = gtf_ref,
                           bamfile = file)
        Popen(cl).wait()
