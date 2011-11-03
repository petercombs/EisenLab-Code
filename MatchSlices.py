"""
Takes RNAseq data from fly slices and matches those slices to virtual slices
taken from a specified axis of data from a Virtual Embryo file (e.g. from the
BDTN project, http://bdtnp.lbl.gov/).  Will report the best match for each time
point.
"""
from __future__ import print_function, division
import PointClouds as pc
import sys
import numpy as np
from argparse import ArgumentParser
from scipy import stats

def parse_args():
    """Parse command line arguments"""
    argparser = ArgumentParser(description='Match gene expression data from a'
                               'single slice to the corresponding slice from '
                               'the BDTNP data')
    argparser.add_argument('-p', '--pre_calc_data', default=None, type=open,
                           help="A pickle file containing the starts and"
                           " virtual slices")
    argparser.add_argument('-w', '--slice_width', default=50.0, type=float, 
                           help='Thickness of the slice in RNAseq data '
                           '(default: 50um)')
    argparser.add_argument('-x', dest='axis', action='store_const', const='x',
                           default='x',
                          help='Slices along the x axis in BDTNP data '
                           '(DEFAULT)')
    argparser.add_argument('-y', dest='axis', action='store_const', const='y',
                          help='Slices along the y axis in BDTNP data')
    argparser.add_argument('-z', dest='axis', action='store_const', const='z',
                          help='slices along the z axis in BDTNP data')
    argparser.add_argument(dest='pointcloud', type=open,
                          help='Point Cloud (vpc) file from BDTNP)')
    argparser.add_argument(dest='expr_file', nargs='+', type=open,
                          help='The genes.expr files output by cufflinks (or '
                           'other similar data with gene names in column 1 and '
                           'FPKM in column 6')
    argparser.add_argument('-t', '--translation-table', type=open,
                           default=False, help="A table of FBgn to real name "
                           'correspondences')
    argparser.add_argument('--sum', dest='reduce_fcn', default=np.sum,
                           action='store_const', const=np.sum,
                           help="Use Sum for virtual slicing (default)")
    argparser.add_argument('--mean', dest='reduce_fcn', action='store_const',
                           const=np.mean, help="Use Mean for virtual slicing")

    argparser.add_argument('--spearman', dest='comp_fcn', action='store_const',
                           default=stats.spearmanr, const=stats.spearmanr,
                           help="Use Spearman correlation for comparison "
                           "(default)")
    argparser.add_argument('--pearson', dest='comp_fcn', action='store_const',
                           const=stats.pearsonr,
                           help='Use Pearson correlation for comparison')
    argparser.add_argument('--name-col', '-n', default=0, type=int,
                           help="0-indexed Column with the gene name")
    argparser.add_argument('--expr-col', '-e', default=5, type=int,
                           help="0-indexed Column with the expresion")

    return argparser.parse_args()


def choose_axis(axis):
    """Converts x,y,z (as strings) to appropriate axis number"""
    axischooser = {'x':0, 0:0, 'y':1, 1:1, 'z':2, 2:2}
    return axischooser[axis]

def get_slicestarts(posarray, axis, resolution, width):
    "Returns the desired starts of virtual slices"
    # Start and stop of slices
    start = posarray[:, axis, :].min() - width
    stop = posarray[:, axis, :].max() + width

    return np.arange(start, stop, resolution)

def get_datasize(exparray):
    """Returns the number of genes and timepoints.

    Like shape(exparray), but with exactly 2 outputs"""

    datasize = np.shape(exparray)
    ngenes = datasize[1]

    if len(datasize) < 3:
        ntimes = 1
    else:
        ntimes = datasize[2]

    return ngenes, ntimes



def virtual_slice(exparray, posarray, axis='x', width=50.0, resolution=1.0,
                  reduce_fcn = np.sum):
    """ Calculate possible slices along an axis.

    Returns an array of slice start positions and an nslices * ngenes * ntimes
    array containing the slice data.

    width is the width of each slice, resolution is the spacing between virtual
    slices.

    If exparray is only a 2-dimensional array (i.e. only from a single time
    point), this will automatically reshape it to be a 3 dimensional array, with
    a singleton dimension.

    """
    # Convert axis specification to a usable value.
    axis = choose_axis(axis)


    slicestarts = get_slicestarts(posarray, axis, resolution, width)

    nslices = len(slicestarts)
    ngenes, ntimes = get_datasize(exparray)

    # Expand into 3D if only 1 time point
    if ntimes == 1:
        exparray = np.reshape(exparray, (-1, ngenes, 1))

    allslices = np.zeros((nslices, ngenes, ntimes))
    for k in range(ntimes):
        # Poor man's progress bar
        print("\nTime %d: " % k, end='')
        next_print = 0
        sys.stdout.flush()

        # Actually do the slicing:
        for i, slicepos in enumerate(slicestarts):
            allslices[i, :] = reduce_fcn(exparray[
                                  (slicepos <= posarray[:, axis, k])
                                * (posarray[:, axis, k] < slicepos + width)],
                              axis=0)

            # More progress bar
            if i / nslices > next_print:
                print('.', end='')
                sys.stdout.flush()
                next_print += .05

    print()

    return slicestarts, allslices


def parse_translation_table(translation_table):
    """ Turns a translation table into a pair of dictionaries

    Will automatically determine which column has the FBgn id numbers.
    """
    fbgn2name = {}
    name2fbgn = {}
    line = "Failed before reading!"
    try:
        for line in translation_table:
            fbgn, name = line.split()
            if 'FBgn' in name:
                name, fbgn = fbgn, name
            fbgn2name[fbgn] = name
            name2fbgn[name] = fbgn

    except ValueError:
        print('Translation table must be whitespace separated lines with '
              'one element as the FBgn designation, and the other the '
              'canonical name.  Failed on line:', file=sys.stderr)
        print(line, file=sys.stderr)
        print('Continuing without full table!', file=sys.stderr)

    return fbgn2name, name2fbgn


if __name__ == "__main__":
    args = parse_args()
    bdtnp_parser = pc.PointCloudReader(args.pointcloud)

    gene_names = bdtnp_parser.get_gene_names()
    gn_list = list(gene_names)

    fbgn2name = {}
    name2fbgn = {}
    if args.translation_table:
        fbgn2name, name2fbgn = parse_translation_table(args.translation_table)

    print("Loading Data...")

    if args.pre_calc_data:
        import cPickle as pickle
        starts = pickle.load(args.pre_calc_data)
        slices = pickle.load(args.pre_calc_data)
    else:
        sys.stdout.flush()
        all_data = [row for row in bdtnp_parser]
        exparray, posarray = bdtnp_parser.data_to_arrays()

        print("Doing virtual slicing")
        starts, slices = virtual_slice(exparray, posarray, axis=args.axis,
                                       width=args.slice_width,
                                       reduce_fcn=args.reduce_fcn)
    nslices, ngenes, ntimes = np.shape(slices)

    for expr_file in args.expr_file:
        print('*'*len(expr_file.name))
        print(expr_file.name)
        print('*'*len(expr_file.name))
        expr = {}
        for line in expr_file:
            linedat = line.split()
            name = linedat[args.name_col]
            if ((name not in gene_names) and
                (name not in fbgn2name or (fbgn2name[name] not in gene_names))):
                continue
            fpkm = float(linedat[args.expr_col])


            expr[fbgn2name[name] if name in fbgn2name else name] = fpkm

        # Turn expression matrix into a list
        expr_l = [(expr[gene] if gene in expr else 0) for gene in gn_list]
        corrmat = np.empty((ntimes, nslices))

        for time in range(ntimes):
            bestcorr = 0
            bestslice = -1
            for slice in range(nslices):
                #corr = np.abs(np.corrcoef(expr_l, slices[slice, :, time]))[1,0]
                #corr = stats.pearsonr(expr_l, slices[slice, :, time])[0]
                corr = stats.spearmanr(expr_l, slices[slice, :, time])[0]
                corrmat[time, slice] = corr
                if corr > bestcorr:
                    bestcorr = corr
                    bestslice = slice
            print("At time %2d, best at #%3d (%+8f), r=%f"
                  % (time, bestslice, starts[bestslice], bestcorr))

