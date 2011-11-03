from random import randint
import sys

def read_n(fobj, n):
    ret = [True]
    while ret[-1]:
        ret = [fobj.readline() for i in range(n)]
        yield ret

if __name__ == "__main__":
    for (h1, seq, h2, q) in read_n(file('unmapped_s_5_1.fq'), 4):
        print "H1:", h1, "Seq: ", seq, "H2: ", h2, "Q: ", q
        sys.stdout.flush()
        sys.stderr.flush()
        if not h1.startswith('@'):
            print h1
            sys.stdout.flush()
        if randint(0,99) == 0:
            print '>'+h1.replace('+', '>').strip()
            print seq.strip()
