from __future__ import print_function, division
from glob import glob
from os import path
from pysam import Samfile
import gzip

import pandas as pd


def get_protocol(dname):
    return path.basename(dname).split('_')[0]


print(r'''
\begin{table}[htdp]

\caption{Summary statistics for samples.  Mapped reads is the number of reads
      that map at least once.  Percent unique reads indicates, of the reads
      that mapped at least once, the fraction that mapped exactly once.}

\begin{center}
\begin{tabular}{|c|c|c|c|c|} \hline
      Protocol & Total Reads & Mapped Reads (N$\ge$1)
      & \% mapped ($N\ge 1$) & \%  unique\\\hline ''')

config_file = pd.read_table('Parameters/RunConfig.cfg',
                            comment='#').dropna(how='all')
for i, row in config_file.iterrows():
    label = row['Label']
    index = row['Index']
    mbepc = int(row['MBEPC'])
    glob_spec = ('{seqdir}/*/*{label}*/*_{{read}}_*.fastq*'
                 .format(seqdir='sequence', label=label, index=index,
                         id=mbepc))
    rf1 = glob(glob_spec.format(read='R1'))
    if rf1 == []:
        glob_spec = ('{seqdir}/*{index}*/*_{{read}}*.fastq*'
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

    samfile = Samfile(path.join('analysis', label, 'accepted_hits_sorted.bam'))
    n_reads = 0
    n_primary = 0
    n_unique = 0
    for read in samfile:
        n_reads += 1
        n_primary += not read.is_secondary
        n_unique += dict(read.tags)['NH'] == 1

    print(r" {prot} & {total:,} &{mapped:,} "
          r"&{frac:0.1f}\%  &{umapped:0.1f}\%\\"
          .format(prot=protocol,
                  total=total,
                  mapped=n_primary,
                  frac=n_primary/total*100,
                  umapped=n_unique/n_primary*100,))

print(r'''\hline

\end{tabular}
\label{tab:protocols}
\end{center}
\end{table}
      ''')
