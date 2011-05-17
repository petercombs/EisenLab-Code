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
    rn = 0
    for r, rna in enumerate(rnas):
        if gene in rna:
            print(rna, end='\t')
            bestr = 0
            pn = 0
            for p, prot in enumerate(prots):
                if gene in prot:
                    r2 = np.corrcoef(protdata[p], rnadata[r])[0,1]
                    pn += 1
                    print('%0.3f' % r2,
                          end='\t')
                    if r2 > bestr:
                        bestdiff = pn - rn - 1
                        bestr = r2
            print(bestdiff)
            rn += 1
