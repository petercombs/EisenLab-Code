import pandas as pd
from collections import defaultdict
from os import path
from subprocess import Popen, call
from sys import stdout

fasta_ref = 'Reference/AAA/mel_only.fa'
gtf_ref = 'Reference/AAA/mel_only.gtf'
design_fname = 'analysis-multi/design.tab'
analysis_dir = 'analysis-multi'

bamfile_base = 'assigned_dmel_rescued.bam'
cd_base = ('cuffdiff --num-threads 8 -L {conditions} --output-dir {outdir}'
           ' -u --frag-bias-correct {fasta} '
           ' {gtf} {bams}')
cl_base = ('cufflinks --num-threads 8 --output-dir {outdir} -u'
           ' --frag-bias-correct {fasta} '
           ' -G {gtf} {bamfile}')
cuffcmp = ('cuffcompare -o {cuffname} -s {fasta} -CG -r {gtf} {gtf}')
cuffname = path.join(path.dirname(gtf_ref), 'cuffcmp')

runstr = cuffcmp.format(cuffname=cuffname,
                        fasta=fasta_ref,
                        gtf=gtf_ref)
print runstr
stdout.flush()
call(runstr.split())

gtf_ref = cuffname + '.combined.gtf'

design_file = pd.read_table(design_fname)
files = defaultdict(list)
for ix, row in design_file.iterrows():
    files[row['condition']].append(path.join(analysis_dir,
                                             row['Sample'],
                                             bamfile_base))
conditions = sorted(files.keys())

cd = cd_base.format(conditions=','.join(conditions),
                    outdir=analysis_dir,
                    fasta=fasta_ref,
                    gtf=gtf_ref,
                    bams=' '.join([','.join(files[key])
                                   for key in conditions]))

conds_nobcd = [c for c in conditions if 'Bcd' not in c]
cd_nobcd = cd_base.format(conditions=','.join(conds_nobcd),
                          outdir=analysis_dir+'-nobcd',
                          fasta=fasta_ref,
                          gtf=gtf_ref,
                          bams=' '.join([','.join(files[key])
                                         for key in conds_nobcd]))
print cd_nobcd
stdout.flush()
call(cd_nobcd.split())
print cd
stdout.flush()
call(cd.split())
for condition in files:
    for file in files[condition]:
        cl = cl_base.format(outdir=path.dirname(file),
                            fasta=fasta_ref,
                            gtf=gtf_ref,
                            bamfile=file)
        print '-'*30
        print cl
        call(cl.split())
