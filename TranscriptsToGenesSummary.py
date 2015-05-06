import pandas as pd
from progressbar import ProgressBar
from sys import argv


if __name__ == "__main__":
    read_table_args = dict(keep_default_na=False, na_values=['---', ''], index_col=0)
    gtf = open('Reference/mel_good.gtf')
    transcripts = pd.read_table(argv[1], **read_table_args)
    genes = pd.DataFrame(columns=transcripts.columns)

    i = 0
    for i, row in enumerate(gtf):
        pass
    gtf.seek(0)

    fbtrs_seen = set()
    pbar = ProgressBar(maxval=i).start()
    for i, line in enumerate(gtf):
        pbar.update(i)
        annot = {
            entry.strip().split(' ')[0]: entry.strip().split(' ')[1].strip().strip('"')
            for entry in line.split('\t')[-1].split(';')
            if entry.strip()
                }
        if 'transcript_id' in annot and 'gene_name' in annot:
            fbtr = annot['transcript_id']
            gene_name = annot['gene_name']
            if fbtr in fbtrs_seen or fbtr not in transcripts.index:
                continue
            if gene_name not in genes.index:
                genes.ix[gene_name] = 0
            genes.ix[gene_name] += transcripts.ix[fbtr]
            fbtrs_seen.add(fbtr)
    genes.to_csv('analysis/summary_fbtrsum.tsv',
                 sep='\t',
                 na_rep='---',
                 float_format='%8.2f',)


