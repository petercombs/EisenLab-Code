#!/usr/bin/env python

from sys import stdin

chroms = {'scer_ref|NC_001133|': 'chrI',
          'scer_ref|NC_001134|': 'chrII',
          'scer_ref|NC_001135|': 'chrIII',
          'scer_ref|NC_001136|': 'chrIV',
          'scer_ref|NC_001137|': 'chrV',
          'scer_ref|NC_001138|': 'chrVI',
          'scer_ref|NC_001139|': 'chrVII',
          'scer_ref|NC_001140|': 'chrVIII',
          'scer_ref|NC_001141|': 'chrIX',
          'scer_ref|NC_001142|': 'chrX',
          'scer_ref|NC_001143|': 'chrXI',
          'scer_ref|NC_001144|': 'chrXII',
          'scer_ref|NC_001145|': 'chrXIII',
          'scer_ref|NC_001146|': 'chrXIV',
          'scer_ref|NC_001147|': 'chrXV',
          'scer_ref|NC_001148|': 'chrXVI',
          'scer_ref|NC_001224|': 'chrMito',
         }

for line in stdin:
    for old, new in chroms.iteritems():
        line = line.replace(old, new)
    print line.strip()

