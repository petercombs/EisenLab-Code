from argparse import ArgumentParser
from glob import glob
from fnmatch import fnmatch, filter
from subprocess import Popen, PIPE
from os import path 

def remote_glob(server, dirpath, indexstr):
    proc = Popen(['ssh', server, 'ls %s' % dirpath], stdout=PIPE)
    stdout, stderr = proc.communicate()
    results = [path.join(dirpath,fname.strip()) 
               for fname in stdout.splitlines() 
               if fnmatch(fname, indexstr)]

    return results

def remote_bowtie(server, bowtie='bowtie', args):
    pass

def parse_args():
    parser = ArgumentParser(description='Run bowtie on given indices, then '
                            'pipes results through eXpress')
    parser.add_argument('--remote-server', '-r')
    parser.add_argument('--bowtie-args', '-a', 
                        default='--very-sensitive-local -p 12 -o 1 -t')
    parser.add_argument('--bowtie-index', '-x')
print remote_glob('deme.qb3.berkeley.edu', 'data/PAC03-4Slices/sequence/', '*index2*')


