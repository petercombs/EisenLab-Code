from __future__ import print_function, division
import PointClouds as pc
from argparse import ArgumentParser

def parse_args()
    argparser = ArgumentParser('Match gene expression data from a single slice'
                               'to the corresponding slice from the BDTNP data')
    argparser.add_argument('-w', '--slice_width', default=50.0, type=float, 
                           help='Thickness of the slice in RNAseq data '
                           '(default: 50um)')
    argparser.add_argument('-x', dest='axis', action='store_const', const='x',
                           default='x',
                          help='Slices along the x axis in BDTNP data (DEFAULT)')
    argparser.add_argument('-y', dest='axis', action='store_const', const='y',
                          help='Slices along the y axis in BDTNP data')
    argparser.add_argument('-z', dest='axis', action='store_const', const='z'
                          help='slices along the z axis in BDTNP data')

    args = argparser.parse_args()
    print(args)
    return args


if __name__ == "__main__":
    args = parse_args()
