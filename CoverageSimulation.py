from __future__ import division
from scipy import random
from numpy import (zeros, zeros_like, arange,
                   unique, mean, log10)
from scipy.stats import linregress
from progressbar import ProgressBar

avg_size = 1500.
n_reps = 100

regs = []
all_sizes = [100., 500., 1000., 1500., 3000., 5000., 10000.]
all_sizes = [1500.]

for avg_size in all_sizes:
    cs = arange(1, 5*avg_size)
    fs = zeros_like(cs)

    pb = ProgressBar()

    all_fracs = []
    for c in pb(cs):
        fracs = zeros(n_reps)
        for i in range(n_reps):
            fracs[i] = len(unique(
                random.randint(0, avg_size, c))
                          )/avg_size
            # Use randint to generate random read positions,
            # then count the number of unique start sites
            # and normalize by size of the transcript
            all_fracs.append((c/avg_size, fracs[i]))
        fs[c-1] = mean(fracs)

    print('-'*30)
    print(avg_size)
    print('-'*30)
    regs.append(linregress(log10(cs[fs<.1]/avg_size),
                           log10(fs[fs<.1])))
    print(regs[-1])
