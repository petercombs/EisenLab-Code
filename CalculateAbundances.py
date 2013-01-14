import pandas as pd
from collections import defaultdict
from os import path
from subprocess import Popen, call

fasta_ref = 'Reference/AAA/mel_only.fa'
gtf_ref = 'Reference/AAA/mel_only.gtf'
design_fname = 'analysis-multi/design.tab'
analysis_dir = 'analysis-multi'

bamfile_base = 'assigned_dmel_rescued.bam'
cd_base = 'cuffdiff -p 8 -L {conditions} -o {outdir} -u -N -b {fasta} {gtf} {bams}'
cl_base = 'cufflinks -p 8 -o {outdir} -u -N -b {fasta} -G {gtf} {bamfile}'
cuffcmp = 'cuffcompare -o {cuffname} -s {fasta} -CG -r {gtf} {gtf}'
cuffname = path.join(path.dirname(gtf_ref), 'cuffcmp')

runstr = cuffcmp.format(cuffname = cuffname,
                        fasta = fasta_ref,
                        gtf = gtf_ref)
print runstr
call(runstr.split())

gtf_ref =  cuffname + '.combined.gtf'

design_file = pd.read_table(design_fname)
files = defaultdict(list)
for ix, row in design_file.iterrows():
    files[row['condition']].append(path.join(analysis_dir,
                                             row['Sample'],
                                             bamfile_base)
                                  )
conditions = sorted(files.keys())

cd = cd_base.format(conditions = ','.join(conditions),
                    outdir = analysis_dir,
                    fasta = fasta_ref,
                    gtf = gtf_ref,
                    bams = ' '.join([','.join(files[key]) for key in conditions]))

print cd
call(cd.split())
for condition in files:
    for file in files[condition]:
        cl = cl_base.format(outdir = path.dirname(file),
                           fasta = fasta_ref,
                           gtf = gtf_ref,
                           bamfile = file)
        print '-'*30
        print cl
        call(cl.split())
