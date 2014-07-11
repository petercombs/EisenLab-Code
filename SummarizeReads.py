from __future__ import print_function, division
from glob import glob
from os import path
from pysam import Samfile
import gzip

import pandas as pd

def get_protocol(dname):
    return path.basename(dname).split('_')[0] 

dnames = glob('analysis/*_V*')
protocols = set(get_protocol(dname) for dname in dnames)

protocol_map = dict(
    Truseq = 'TruSeq',
    Clontech = 'Clontech',
    EpicentreOligodT = 'TotalScript',
    SMART2 = 'Smart-seq2',
    SMART2dil2 = 'Smart-seq2, 2.5 fold dilution',
    SMART2dil5 = 'Smart-seq2, 5 fold dilution',
)

cost_map = dict(
    Truseq = '\\$45',
    Clontech = '\\$105',
    EpicentreOligodT = '\\$115',
    SMART2 = '\\$55',
    SMART2dil2 = '\\$28',
    SMART2dil5 = '\\$20',
)



print(r'''
\begin{table}[htdp]

\caption{Summary statistics for samples. Cost is estimated per sample assuming a
     large number of libraries at US catalog prices as of May 2014. Total Reads
     and Mapped reads are for pooled libraries, run in a single HiSeq lane each
     for experiments 2 and 3}
\begin{center}
\begin{tabular}{|c|c|c|c|c|c|c|} \hline
      Expt & Protocol & \% {\em D. virilis} & Cost &  Yield & Total Reads & Mapped Reads \\\hline ''')

config_file = pd.read_table('Parameters/RunConfig.cfg',
                           comment='#').dropna(how='all')
for i, row in config_file.iterrows():
    label = row['Label']
    index = row['Index']
    mbepc = int(row['MBEPC'])
    carrier = row['CarrierID']
    species = row['CarrierSpecies']
    glob_spec = ('{seqdir}/*/*{label}*/*_{{read}}_*.fastq*'
               .format(seqdir='sequence', label=label, index=index, id=mbepc))
    rf1 = glob(glob_spec.format(read='R1'))
    if rf1 == []:
        glob_spec = ('{seqdir}/*{id}*{index}*/*_{{read}}*.fastq*'
                    .format(seqdir='sequence', label=label, index=index,
                            id=mbepc))
        rf1 = glob(glob_spec.format(read='R1'))
        if rf1 == []:
            print("Warning: no sequence for ", label, index)
            print(glob_spec.format(read='R1'))

    reads_last = 0
    for line in gzip.open(sorted(rf1)[-1]):
        reads_last += 1
    print('% {}: {}'.format(sorted(rf1)[-1], reads_last/4))
    total = int(4e6*(len(rf1) - 1) + reads_last/4)
    protocol = get_protocol(label)

    mapped = set()

    for read in Samfile(path.join('analysis', label,
                                  'accepted_hits.bam')):
        mapped.add(read.qname)

    print(r"{expt} & {prot} & {frac}\% & {cost} & {libyield:,} fmole & {total:,} &{mapped:,} \\"
          .format(prot = protocol_map[protocol],
                  expt = 2 if 'dil' not in label else 3,
                  cost = cost_map[protocol],
                  libyield = row['Yield'],
                  frac = int(label[label.rfind('V')+1:
                                   label.rfind('V')+3]),
                  total = total,
                  mapped = len(mapped)))
    cost_map[protocol] = '"'
    protocol_map[protocol] = '"'

print(r'''\hline

\end{tabular}
\label{tab:protocols}
\end{center}
\end{table}
      ''')
