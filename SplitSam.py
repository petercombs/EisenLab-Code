from __future__ import print_function
from collections import defaultdict

import sys

class OpenDefaultDict(dict):
    def __missing__(self, key):
        self[key] = value = open(str(key), 'w')
        return value


files = OpenDefaultDict()

for line in open(sys.argv[1]):
    chrom = line.split()[2].split('|')[-2]
    files[chrom].write(line)


