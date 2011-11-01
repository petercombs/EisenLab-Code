import sys


ids = {line.split()[0] for line in file(sys.argv[1])}

printlines = 0

for line in file(sys.argv[2]):
    id = line.split()[0][1:]
    if id not in ids:
        printlines = 4

    if printlines:
        print line.strip()
        printlines -= 1
