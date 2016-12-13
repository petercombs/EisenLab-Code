import DistributionDifference as dd
import PlotUtils as pu
import numpy as np
from matplotlib.pyplot import subplot, pcolormesh, xlim, savefig, title

gene = 'run'
flows = dd.earth_mover(bcds['cyc14D_rep1'].ix[gene], wts['cyc14D'].ix[gene], True)

from_array = np.zeros(max(f[0] for f in flows)+1)
to_array = np.zeros(max(f[1] for f in flows)+1)
for f in flows:
    from_array[f[0]] += f[2]*abs(float(f[0])/(len(from_array)-1) - float(f[1])/(len(to_array)-1))
    to_array[f[1]] += f[2]*abs(float(f[0])/(len(from_array)-1) - float(f[1])/(len(to_array)-1))

subplot(2, 1, 1)
title(gene)
pcolormesh(to_array.reshape((1,-1)), cmap=pu.ISH)
xlim(0, len(to_array))
subplot(2, 1, 2)
pcolormesh(from_array.reshape((1,-1)), cmap=pu.ISH_CMS_6[1])
xlim(0, len(from_array))
savefig('analysis/results/tmp-flow.png')
