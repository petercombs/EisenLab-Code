from __future__ import print_function
#python3 compatibility

import numpy as np
import csv
from PointClouds import PointCloudReader




pcr = PointCloudReader(open('../Data/D_mel_wt_atlas_r2.vpc'))

protidxs, prots= zip(*[(i, c) for i, c in enumerate(pcr.column) if 'P' in c])
rnas = [p.replace('P', '') for p in prots]
rnas = [r for r in rnas if r in pcr.column]
rnaidxs = [pcr.column.index(r) for r in rnas]

geneset = set([r.split('_')[0] for r in rnas])

protdata = [[] for p in prots]
rnadata = [[] for r in rnas]

for line in pcr:
    for i, p in enumerate(protidxs):
        protdata[i].append(line[p])
    for i, r in enumerate(rnaidxs):
        rnadata[i].append(line[r])


for gene in geneset:
    print('', end='\t')
    for prot in prots:
        if gene in prot:
            print(prot, end='\t')
    print('')
    for r, rna in enumerate(rnas):
        if gene in rna:
            print(rna, end='\t')
            for p, prot in enumerate(prots):
                if gene in prot:
                    print('%0.3f' %np.corrcoef(protdata[p], rnadata[r])[0,1], 
                          end='\t')
            print('')
