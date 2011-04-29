from argparse import ArgumentParser, FileType
import sys

parser = ArgumentParser()
parser.add_argument('-r', dest="reference", 
                  help="File with conversions between names",
                  default="Reference/dmelfbgns.txt")
parser.add_argument('-k', dest="keys", action='append',
                  help="columns to output", type=int, default = [])
parser.add_argument('-i', dest='idkey', type=int, default=0,
                  help="The key with the FBgn number (printed by default)")
parser.add_argument('filelist', type=FileType('r'), default=[sys.stdin],
                    nargs='*', action='store')

options = parser.parse_args()
NameFile = options.reference or 'Reference/dmelfbgns.txt' 
idkey = options.idkey
fs = '\t'

FB2name = {}

for line in file(NameFile):
    s = line.split()
    FB2name[s[0]] = s[1]

sys.stderr.write(str(options.filelist) + '\n')

for infile in options.filelist:
    sys.stderr.write(str(infile) + "\n")
    for line in infile:
        s = line.split()
        print (s[idkey] in FB2name and FB2name[s[idkey]]) or s[idkey],
        print fs,
        print fs.join(s[key] for key in options.keys)
