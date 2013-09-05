#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from StringIO import StringIO
from sys import argv
from os import chmod, path
import progressbar as pb
import PlotUtils


startswith = lambda x: lambda y: y.startswith(x)
bad_slices = {'emb4_sl10_FPKM', 'emb6_sl01_FPKM', 'emb7_sl01_FPKM',
              'emb7_sl07_FPKM', 'emb7_sl08_FPKM'}


data = pd.read_table(argv[1], index_col=0).select(lambda x: x not in
                                                       bad_slices, axis=1)

headers = sorted({column[:column.find('_sl')] for column in data.columns})

gene_index = {}
for gene in data.index:
    gene_index[gene] = gene
    if ',' in gene:
        for gn in gene.split(','):
            gene_index[gn] = gene


if len(argv) == 2:
    genes = data.index
else:
    genes = argv[1:]

widgets = ['Something: ', pb.Percentage(), ' ', pb.Bar(marker=pb.RotatingMarker()),
                      ' ', pb.ETA(), ]
pbar = pb.ProgressBar(widgets=widgets)
fig = plt.figure(figsize=(len(headers)+1,1))
for gene in pbar(genes):
    if gene not in gene_index:
        continue
    all_max = min(max(data.ix[gene_index[gene]]), 10000)
    dats = []
    for emb, cyc in enumerate(headers):
        ax = plt.subplot(1,len(headers),emb+1)
        dat = data.ix[(gene_index[gene],),].select(startswith(cyc),axis=1)
        # indexing with a tuple gives us a dataframe
        dats.append(dat)
        #plt.pcolormesh(np.reshape(dat, (1, -1)),
        #              cmap=plt.cm.Blues,
        #              vmin=0, vmax=all_max)
        #ax.set_xlim([0, len(dat)])
        #ax.set_xticks([0, len(dat)])
        #ax.set_xticklabels('AP')
        #ax.set_yticks([])
        #ax.set_title(cyc)
        #plt.tight_layout()

    for gene in gene.split(','):
        gene = ''.join((l+'_' if l.isupper() else l) for l in gene)
        out_path = path.join(path.dirname(argv[1]),
                              'imgs',
                              gene+'.png')
        #plt.savefig(out_path, format='png')
        #plt.clf()
        #chmod(out_path, 0644)
        PlotUtils.svg_heatmap(tuple(dats), out_path+'.svg',
                              norm_rows_by=all_max,
                              box_height=40, total_width=100,
                              draw_box=True,
                              first_col='', last_col='')
        chmod(out_path+'.svg', 0644)
