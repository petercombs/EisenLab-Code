from __future__ import division
import pandas


def main():
    df = pandas.read_table('CountConfig.tab', index_col=5)
    df.dropna(how='any')

    for rowname, row in df.iterrows():
        mel_reads = row['mel_reads']
        carrier_reads = row['carrier_reads']
        carrier_conc = row['carrier_conc']

        mel_conc = carrier_conc * mel_reads / carrier_reads
        print rowname,
        print mel_conc, 'ng total RNA'
    return df

if __name__ == "__main__":
    df = main()
