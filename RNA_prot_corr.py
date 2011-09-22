from __future__ import print_function
#python3 compatibility

import numpy as np
from PointClouds import PointCloudReader
import time
import pickle
import multiprocessing as mp

def calc_diffusion(S):
    print("Caculating diffusion")
    print("-"*70)
    print(S)
    print("-"*70)


    try:
        f = open('diffused%03d.pkl'%S)
        diffused = pickle.load(f)
        print("Sweet, grabbed it from a file")
    except IOError:
        diffused = {}
        f = open('diffused%03d.pkl'%S, 'w')
        for gene in gene_names:
            if gene+"P" not in gene_names:
                continue
            coord = gene_name_list.index(gene)
            print("Working on gene: " + gene)
            rna_diffused = np.zeros((nuclei, times))
            for t in range(times):
                print("Timepoint %d" % t)
                for n in range(nuclei):
                    rna_diffused[n, t] = sum(1 / np.sqrt(np.pi * S)
                                             * np.exp(-dists[n,:,t] / S)
                                             * expr[:,coord, t])
            diffused[gene] = rna_diffused
        pickle.dump(diffused, f)

    print("Sweet, diffusion is done")
    gene_corrs = {}
    for gene in diffused:
        print("-"*50)
        print(gene)
        print("-"*50)
        prot = gene_name_list.index(gene+"P")
        corrs = []
        for rna_time in range(times):
            corrs.append([])
            for prot_time in range(times):
                if sum(expr[:, prot, prot_time]) and sum(diffused[gene][:,rna_time]):
                    corr = np.corrcoef(expr[:,prot, prot_time], diffused[gene][:,rna_time])[0,1]
                    corrs[-1].append(corr)
                    print(corr, end='\t')
            print()
        gene_corrs[gene] = corrs
    return gene_corrs



pcr = PointCloudReader(open('../Data/D_mel_wt_atlas_r2.vpc'))

# use proteins to get the list of mRNAs.
protidxs, prots= zip(*[(i, c) for i, c in enumerate(pcr.column) if 'P' in c])

rnas = [p.replace('P', '') for p in prots]
# Some RNAs aren't measured at certain time points
rnas = [r for r in rnas if r in pcr.column]
rnaidxs = [pcr.column.index(r) for r in rnas]

geneset = set([r.split('_')[0] for r in rnas])

protdata = [[] for p in prots]
rnadata = [[] for r in rnas]

expr, pos = pcr.data_to_arrays()

nuclei, genes, times = np.shape(expr)

start_dist = time.time()
print("Ready")
dists = np.empty((nuclei, nuclei, times), dtype=np.float32)
pos.resize((nuclei, 1, 3, times))
posT = pos.reshape((1, nuclei, 3, times))
print("Steady")
for i in range(times):
    print(i)
    dists[:,:,i] = np.sum((pos[:,:,:,i] - posT[:,:,:,i])**2, axis=2)
stop_dist = time.time()
print("Took %fs to calculate distances" % (stop_dist - start_dist))

print("Loading expression data old-school")
for line in pcr:
    for i, p in enumerate(protidxs):
        protdata[i].append(line[p])
    for i, r in enumerate(rnaidxs):
        rnadata[i].append(line[r])
print("Done with that loading as well")

gene_names = pcr.get_gene_names()
gene_name_list = list(gene_names)
#S = 10 # 4 D t
Ss = [1, 3, 5, 10, 15, 20, 40, 80]

p = mp.Pool(3)
p.map(calc_diffusion, Ss)

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
