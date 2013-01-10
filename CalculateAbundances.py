import pandas as pd
from collections import defaultdict
from os import path

design_file = pd.read_table('analysis-multi/design.tab')
files = defaultdict(list)
for ix, row in design_file.iterrows():
    files[row['condition']].append(path.join('analysis-multi',
                                          row['Sample'],
                                          'assigned_dmel_rescued.bam')
                                )
conditions = sorted(files.keys())

cd_base = 'cuffdiff -p 4 -L {conditions} -o {outdir} -u -N -b {fasta} {gtf} {bams}'
cd = cd_base.format(conditions = ','.join(conditions),
                    outdir = 'analysis',
                    fasta = 'Reference/AAA/mel_only.fa',
                    gtf = 'Reference/AAA/mel_only.gtf',
                    bams = ' '.join([','.join(files[key]) for key in conditions]))

print cd
