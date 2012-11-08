from __future__ import division
import pandas

def main():
    df = pandas.read_table('CountConfig.tab')
    df.dropna()

    for row in df:
        mel_reads = row.mel_reads
        carrier_reads = row.carrier_reads
        carrier_conc = row.carrier_conc

        mel_conc = mel_reads * carrier_reads / carrier_conc
        print row.sample,
        print mel_conc, 'ng total RNA'
