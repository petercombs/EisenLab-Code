from __future__ import print_function
import sys
import time

start_lib = time.time()
ids = {line.split()[0] for i, line in enumerate(file(sys.argv[1]))
      if ((i % 1000000 == 0 and print('.', file=sys.stderr, end=' ')) or True)}

print(time.time() - start_lib, file=sys.stderr)

printlines = 0

for line in file(sys.argv[2]):
    if printlines == 0:
        id = line.split()[0][1:-2]
        if id not in ids:
            printlines = 4
    else:
        print(line.strip())
        printlines -= 1
