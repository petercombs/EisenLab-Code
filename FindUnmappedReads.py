from __future__ import print_function
import sys


ids = {line.split()[0] for i, line in enumerate(file(sys.argv[1]))
      if ((i % 1000000 == 0 and print('.', file=sys.stderr, end=' ')) or True)}

printlines = 0

for line in file(sys.argv[2]):
    id = line.split()[0][1:]
    if id not in ids:
        printlines = 4

    if printlines:
        print(line.strip())
        printlines -= 1
