from sys import stdin, stdout, argv, stderr
import pandas as pd
import re
from collections import defaultdict

ToFilter = re.compile('(snoRNA|CR[0-9]{4}|Rp[ILS]|mir-|tRNA|unsRNA|snRNA|'
            'snmRNA|scaRNA|rRNA|RNA:|mt:)' )
FBgn_finder = re.compile('FBgn[0-9]{7}')

orthologs = pd.read_table(argv[1], skiprows=4, converters={'GeneSymbol':str})

orthos = defaultdict(list)
for i, row in orthologs.iterrows():
    orthos[row['Ortholog_FBgn_ID']].append(row['GeneSymbol'])

for line in stdin:
    if ToFilter.findall(line):
        continue
    fbgn = FBgn_finder.findall(line)[0]
    if fbgn not in orthos:
        stdout.write(line)
        continue
    try:
        for ortholog in orthos[fbgn]:
            if not ToFilter.findall(ortholog):
                stdout.write(line)
                continue
    except Exception as exc:
        stderr.write(fbgn)
        stderr.write(str(orthos[fbgn]))
        stderr.write(str(exc))
        raise exc

