from __future__ import print_function
from sys import argv

infile = open(argv[1])

print(infile.readline(), end="")

last_sig_locus = ""
last_sig_line = ""
last_sig_fold = 0
for line in infile:
    data = line.split()
    fold = float(data[8])
    if data[-1] == "yes":
        if last_sig_locus == data[2] and (last_sig_fold * fold) < 0:
            print(last_sig_line, end="")
            print(line, end="")
        last_sig_locus = data[2]
        last_sig_line = line
        last_sig_fold = fold

