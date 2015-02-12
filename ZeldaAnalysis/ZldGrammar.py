import pandas as pd
from glob import glob
from collections import defaultdict
from os import path

TFs = defaultdict(list)
zld_bind = pd.read_table('prereqs/journal.pgen.1002266.s005.xls', skiprows=1)

for file in glob('Reference/peaks-25-dm3/*.bed'):
    tf = path.basename(file).split('_')[0]
    print '"{}"'.format(tf)
    try:
        for i, row in pd.read_table(file, header=None).iterrows():
            TFs[row[0], row[1]].append(tf)
    except:
        print "Failed on file", file, tf

    for value in TFs.values():
        assert value

nearest = []

for chr, pos in zip(zld_bind.Chr, zld_bind.Peak):
    i = 0
    while i < 200:
        if (chr, pos + i) in TFs:
            nearest.append((chr, pos+i, TFs[chr, pos+i]))
            break
        elif (chr, pos - i) in TFs:
            nearest.append((chr, pos+i, TFs[chr, pos-i]))
            break
        else:
            i+= 1
    else:
        nearest.append(None)

