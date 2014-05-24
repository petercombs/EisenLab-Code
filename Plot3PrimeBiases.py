from __future__ import print_function
from os import path
from collections import defaultdict
from matplotlib.pyplot import plot, legend, figure, xlabel, ylabel
from glob import glob
import numpy as np

xs = defaultdict(list)
ys = defaultdict(list)
for fname in glob('analysis/results/picard/*.txt'):
    basename = path.basename(fname)
    protocol = basename[:basename.find('_V')]
    print(protocol)
    results_file = open(fname)
    xs[protocol].append([])
    ys[protocol].append([])
    for line in results_file:
        if 'HISTOGRAM' in line:
            break
    results_file.readline()
    for line in results_file:
        if not line.strip(): continue
        x, y = line.strip().split()
        xs[protocol][-1].append(float(x))
        ys[protocol][-1].append(float(y))


cs = dict(zip(xs.keys(), 'rgbckm'))
seen = set()
figure()
for protocol in xs:
    ys_prot = np.array(ys[protocol]).T
    plot(xs[protocol][0], ys_prot[:,0],
         cs[protocol], label=protocol_map[protocol])
    plot(xs[protocol][0], ys_prot[:,1:],
         cs[protocol])

legend(loc='upper left', ncol=2)
xlabel('Normalized transcript position')
ylabel('Normalized coverage')
