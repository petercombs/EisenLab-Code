from __future__ import print_function
import sys
import time

start_lib = time.time()
# Sequence is in the 9th column, ID in the 0th
ids = {line.split()[9] for i, line in enumerate(file(sys.argv[1]))
      if ((i % 1000000 == 0 and print('.', file=sys.stderr, end=' ')) or True)}

print(time.time() - start_lib, file=sys.stderr)

printlines = 0

grab_next = False
for i, line in enumerate(file(sys.argv[2])):
    # Sequence line
    if i % 4 == 1 and line.strip() not in ids:
        seq = line.strip()
        grab_next = True
    if grab_next:
        id = line.strip().split()[0].replace('+', '>')
        print id
        print seq
        grab_next = False

