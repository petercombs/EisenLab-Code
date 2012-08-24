import sys
import numpy as np
from os import path
from matplotlib import pyplot as plt


line_num = -1

for line in open(sys.argv[1]):
    if line_num != -1:
        if line_num == 0:
            ID = line.strip().split('"')[1]
        elif line_num == 1:
            pos1 = line.find('[')+1
            pos2 = line.find(']')
            counts = [int(n) for n in line[pos1:pos2].split(', ')]
        elif line_num == 2:
            pos1 = line.find('[')+1
            pos2 = line.find(']')
            expects = [float(n) for n in line[pos1:pos2].split(', ')]
        elif line_num == 3:
            is_sig = line.startswith("SIG")
            if is_sig:
                print "Significant!", ID

        line_num += 1

        if line_num == 4:
            line_num = -1
            if is_sig:
                plt.figure()
                left = np.arange(len(counts))
                plt.bar(left, counts, .3, color='b')
                plt.bar(left+.3, expects, .3, color='g')
                plt.legend(['Observed', 'Expected'])
                plt.title(ID)
                plt.savefig(path.join('figs', ID+'.png'))

    elif line_num == -1 and line.startswith('------'):
        line_num = 0

