import numpy as np
import pandas as pd
from Bio import SeqIO
from DistributionDifference import earth_mover


is_in = lambda y: lambda x: y in x
startswith = lambda y: lambda x: x.startswith(y)

if __name__ == "__main__":
    genome = {rec.id: rec.seq for rec in
              SeqIO.parse('prereqs/dmel-all-chromosome-r5.53.fasta', 'fasta')}

    zld_exp = pd.read_table('analysis/summary.tsv', index_col=0).sort_index()
    wt_exp = pd.read_table('prereqs/WT5.53_summary.tsv', index_col=0).sort_index()
    assert np.all(zld_exp.index == wt_exp.index)

    zld_bind = pd.read_table('prereqs/journal.pgen.1002266.s005.xls', skiprows=1)
    zld_bind.TSS_gene = zld_bind.TSS_gene.apply(str.strip)
    by_gene = zld_bind.groupby('TSS_gene')


    types = {'Intergenic':'N', 'Intronic':'I', 'Promoter':'P', 'UTR5':'5',
             'CDS':'C', 'UTR3':'3'}

    emd_array = pd.Series(data=0.0, index=zld_exp.index)
    n_emd =     pd.Series(data=0,   index=zld_exp.index)
    for zcyc, wcyc in zip('cyc11 cyc13_rep1 cyc14A cyc14B cyc14D'.split(),
                          'cyc11 cyc13      cyc14A cyc14B cyc14E'.split()):
        print "Processing cycle", zcyc, wcyc
        zld_comp = zld_exp.select(is_in(zcyc), axis=1)
        wt_comp  =  wt_exp.select(is_in(wcyc), axis=1)

        for gene in emd_array.index:
            z = zld_comp.ix[gene]
            w =  wt_comp.ix[gene]
            if z.max() > 3 and w.max() > 3:
                emd_array.ix[gene] += earth_mover(w, z)
                n_emd.ix[gene] += 1
        #for gene in wt_exp.index:
            #assert gene in zld_exp.index
            #diff_col[gene] += emd(zld_comp.ix[gene], wt_comp.ix[gene])


    emd_array /= n_emd
    emd_array = emd_array.ix[np.isfinite(emd_array)]
    emd_array.sort()

    dist = 250
    outfh = open('analysis/meme/zld_binds_{}.fa'.format(dist), 'w')
    posfh = open('analysis/meme/pos_{}.fa'.format(dist), 'w')
    negfh = open('analysis/meme/neg_{}.fa'.format(dist), 'w')
    for gene in emd_array.index:
        if emd_array.ix[gene] < 0:
            continue
        if gene not in by_gene.indices:
            print "No zld binding for ", gene
            continue
        for i, peak in enumerate(by_gene.indices[gene]):
            site = zld_bind.ix[peak]
            chr = site['Chr'][3:]
            peak = site['Peak']
            start = peak - dist
            stop = peak + dist
            seq = genome[chr][start:stop]
            seq.id = "{}_{}".format(gene, i)
            r = SeqIO.SeqRecord(seq, id = seq.id, description="{};{}:{}-{}".format(emd_array.ix[gene], chr, start, stop)) 
            SeqIO.write(r, outfh, 'fasta')
            if emd_array.ix[gene] < 0.05:
                SeqIO.write(r, negfh, 'fasta')
            elif emd_array.ix[gene] > 0.2:
                SeqIO.write(r, posfh, 'fasta')

    outfh.close()
    posfh.close()
    negfh.close()

    #zld_fig_genes = zld_exp.select(lambda x: '14A' in x or '11' in x, axis=1)
    #wt_fig_genes = wt_exp.select(lambda x: '14A' in x or '11' in x, axis=1)
#
    #zld_fig_genes = zld_fig_genes[wt_fig_genes.max(axis=1) > 10][:120]
    #wt_fig_genes = wt_fig_genes[wt_fig_genes.max(axis=1) > 10][:120]
#
    #assert (zld_fig_genes.index == wt_fig_genes.index).all()
