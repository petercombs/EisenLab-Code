from __future__ import print_function
#python3 compatibility

import numpy as np
from PointClouds import PointCloudReader
import time

pcr = PointCloudReader(open('../Data/D_mel_wt_atlas_r2.vpc'))

protidxs, prots= zip(*[(i, c) for i, c in enumerate(pcr.column) if 'P' in c])
rnas = [p.replace('P', '') for p in prots]
rnas = [r for r in rnas if r in pcr.column]
rnaidxs = [pcr.column.index(r) for r in rnas]

geneset = set([r.split('_')[0] for r in rnas])

protdata = [[] for p in prots]
rnadata = [[] for r in rnas]

expr, pos = pcr.data_to_arrays()

nuclei, genes, times = np.shape(expr)

start_dist = time.time()
print("Ready")
dists = np.empty((nuclei, nuclei, times), dtype=np.float16)
pos.resize((nuclei, 1, 3, times))
posT = pos.reshape((1, nuclei, 3, times))
print("Steady")
for i in range(times):
    print(i)
    dists[:,:,i] = np.sum((pos[:,:,:,i] - posT[:,:,:,i])**2, axis=2)
stop_dist = time.time()
print("Took %fs to calculate distances" % (stop_dist - start_dist))

for line in pcr:
    for i, p in enumerate(protidxs):
        protdata[i].append(line[p])
    for i, r in enumerate(rnaidxs):
        rnadata[i].append(line[r])

for time in range(times-1):
    rna_diffused = []

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
