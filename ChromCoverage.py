from __future__ import print_function
import svgfig
import sys
import numpy as np

step = 10000
widthscale = 10000

for fname in sys.argv[1:]:
    print('Working on ' + fname)
    chromcov = []
    for line in open(fname):
        cov = int(line.split()[0])
        pos = int(line.split()[1])

        if pos > len(chromcov):
            chromcov.extend([0]*(pos - len(chromcov)))
            chromcov[pos-1] = cov

    polvals = [(start/widthscale, 50*np.mean(chromcov[start:start+step])) 
               for start in range(0, len(chromcov), step)]

    polvals.insert(0,(0,0))
    polvals.append((start/widthscale,0))

    pgon = svgfig.Poly(polvals, mode='lines', loop=True, fill='gray')

    plot = svgfig.Plot(0, start/widthscale, 0, 5, pgon)

    svg = plot.SVG()
    svg.save(fname+'.svg')
        


