#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python
import cgi
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from StringIO import StringIO
from sys import argv, exit
from os import chmod
import progressbar as pb


startswith = lambda x: lambda y: y.startswith(x)
bad_slices = {'emb4_sl10_FPKM', 'emb6_sl01_FPKM', 'emb7_sl01_FPKM',
              'emb7_sl07_FPKM', 'emb7_sl08_FPKM'}


data = pd.read_table('genes.cuff', index_col=0).select(lambda x: x not in
                                                       bad_slices, axis=1)

gene_index = {}
for gene in data.index:
    gene_index[gene] = gene
    if ',' in gene:
        for gn in gene.split(','):
            gene_index[gn] = gene

if len(argv) == 1:
    genes = data.index
else:
    genes = argv[1:]

widgets = ['Something: ', pb.Percentage(), ' ', pb.Bar(marker=pb.RotatingMarker()),
                      ' ', pb.ETA(), ]
pbar = pb.ProgressBar(widgets=widgets)
for gene in pbar(genes):
    if gene not in gene_index:
        continue
    all_max = min(max(data.ix[gene_index[gene]]), 1000)
    fig = plt.figure(figsize=(8,1))
    for emb in range(1,8):
        ax = plt.subplot(1,7,emb)
        dat = data.ix[gene_index[gene]].select(startswith('emb%d' % emb))
        plt.pcolormesh(np.reshape(dat, (1, -1)),
                      cmap=plt.cm.Blues,
                      vmin=0, vmax=all_max)
        ax.set_xlim([0, len(dat)])
        ax.set_xticks([0, len(dat)])
        ax.set_xticklabels('AP')
        ax.set_yticks([])
        plt.tight_layout()

    for gene in gene.split(','):
        gene = ''.join((l+'_' if l.isupper() else l) for l in gene)
        plt.savefig('imgs/'+gene+'.png', format='png')
        chmod('imgs/'+gene+'.png', 0644)
    plt.close()
