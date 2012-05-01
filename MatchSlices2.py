"""
Takes RNAseq data from fly slices and matches those slices to virtual slices
taken from a specified axis of data from a Virtual Embryo file (e.g. from the
BDTN project, http://bdtnp.lbl.gov/).  Will report the best match for each time
point.

In version 2, uses estimates of absolute mRNA abundance from staged RNA-seq
data.  Ideally, something like from Lott et al, 2011
"""
from __future__ import print_function, division
import cPickle as pickle
import PointClouds as pc
import sys
import numpy as np
import progressbar as pbar
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
    argparser.add_argument('-W', '--std-width', default=25.0, type=float,
                           help='Thickness of standard deviation window '
                           '(default: 25um)')
    argparser.add_argument('-R', '--rna-seq', dest='rnaseq', type=open,
                           required=True,
                           help='RNA-seq expression levels for the expected '
                           'time point')
    argparser.add_argument('-T', '--time', type=int, required=True,
                           help='Time point in the Virtual Embryo to compare')


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
    argparser.add_argument(dest='expr_file', type=open,
                          help='The genes.expr files output by cufflinks')

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
    argparser.add_argument('--expr-col', '-e', default=[1], type=int, nargs='+',
                           help="0-indexed Columns with the expresion")

    return argparser.parse_args()


def choose_axis(axis):
    """Converts x,y,z (as strings) to appropriate axis number"""
    axischooser = {'x':0, 0:0, 'y':1, 1:1, 'z':2, 2:2}
    return axischooser[axis]

def get_slicestarts(posarray, axis, resolution, width, std_width):
    "Returns the desired starts of virtual slices"
    # Start and stop of slices
    start = posarray[:, axis, :].min() - width + std_width
    stop = posarray[:, axis, :].max() + width - std_width

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
                  std_width=25.0, reduce_fcn = np.sum):
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


    slicestarts = get_slicestarts(posarray, axis, resolution, width, std_width)

    nslices = len(slicestarts)
    ngenes, ntimes = get_datasize(exparray)

    # Expand into 3D if only 1 time point
    if ntimes == 1:
        exparray = np.reshape(exparray, (-1, ngenes, 1))

    allslices = np.zeros((nslices, ngenes, ntimes))
    print("Calculating means")
    for k in range(ntimes):
        progress = pbar.ProgressBar(widgets=['Time %d: ' % k, pbar.Percentage(),
                                             ' ', pbar.Bar(), pbar.ETA()],
                                    maxval=nslices)

        # Actually do the slicing:
        for i, slicepos in progress(enumerate(slicestarts)):
            allslices[i, :] = reduce_fcn(exparray[
                                  (slicepos <= posarray[:, axis, k])
                                * (posarray[:, axis, k] < slicepos + width)],
                              axis=0)


    allstds = np.zeros((nslices, ngenes, ntimes))
    print("Calculating variation")
    for k in range(ntimes):
        progress = pbar.ProgressBar(widgets=['Time %d' % k, pbar.Percentage(),
                                             ' ', pbar.Bar(), pbar.ETA()],
                                    maxval=nslices)

        for i, slicepos in progress(enumerate(slicestarts)):
            allstds[i,:] = np.std(exparray[
                                            (slicepos
                                             < (posarray[:,axis,k] - std_width))
                                          * (posarray[:,axis,k]
                                             < (slicepos + width +  std_width))
                                          ],
                                  axis=0)


    return slicestarts, allslices, allstds


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

def nan_safe_spearman(a, b):
    """ Calculate the Spearman rank correlation, ignoring any points where one
    or both entries are NaN

    This ought to generally improve comparison across time points, since before
    the differing number of genes that are measured at each time point could
    affect the weighting.
    """
    assert np.shape(a) == np.shape(b)
    nums = ~(np.isnan(a) + np.isnan(b))
    return stats.spearmanr(np.array(a)[nums], np.array(b)[nums])

def load_rnaseq_expr(rnaseq, translate):
    rnaseq_expr = {}
    for line in rnaseq:
        if line.startswith('NAME'): continue
        data = line.split()
        gn_name = data[0]
        if gn_name in translate:
            gn_name = translate[gn_name]
        expr = np.mean([float(i) for i in data[1:]])
        rnaseq_expr[gn_name] = expr

    return rnaseq_expr


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
    rnaseq_expr = load_rnaseq_expr(args.rnaseq, fbgn2name)

    if args.pre_calc_data:
        starts = pickle.load(args.pre_calc_data)
        slices = pickle.load(args.pre_calc_data)
        stds = pickle.load(args.pre_calc_data)
    else:
        exparray, posarray = bdtnp_parser.data_to_arrays()

        print("Doing virtual slicing")
        starts, slices, stds = virtual_slice(exparray, posarray, axis=args.axis,
                                             width=args.slice_width,
                                             std_width=args.std_width,
                                             reduce_fcn=args.reduce_fcn)
        tmp = open('lastslice.pkl', 'w')
        pickle.dump(starts, tmp)
        pickle.dump(slices, tmp)
        pickle.dump(stds, tmp)
        tmp.close()
    nslices, ngenes, ntimes = np.shape(slices)
    assert ntimes > args.time
    nsamples = len(args.expr_col)

    all_likelihoods = np.zeros((nsamples, nslices, ngenes))
    likelihoods = np.zeros((nsamples, nslices))
    column_headers = args.expr_file.readline().split()
    column_headers = [column_headers[i] for i in args.expr_col]

    for line in args.expr_file:
        data = line.split()
        gn_name = data[args.name_col]
        try:
            gn_index = gn_list.index(gn_name)
            gn_exprs = [float(data[i]) for i in args.expr_col]
            if np.NaN in slices[:, gn_index, args.time]: continue
            #print(gn_name)
            for i, expr in enumerate(gn_exprs):
                gn_slices = slices[:, gn_index, args.time]
                gn_stds = stds[:, gn_index, args.time]
                gn_stds = (gn_stds / np.sqrt(sum(gn_slices)) * args.slice_width
                           * nsamples * rnaseq_expr[gn_name])
                gn_slices = (gn_slices / sum(gn_slices) * args.slice_width
                             * nsamples * rnaseq_expr[gn_name])

                delta_L = stats.norm.logpdf(expr, gn_slices, gn_stds)
                delta_L = delta_L.clip(-1e9, 100)
                all_likelihoods[i,:,gn_index] = delta_L
                likelihoods[i,:] += delta_L
        except ValueError:
            pass

    for i, column_header in enumerate(column_headers):
        print(column_header)
        print(starts[np.argmax(likelihoods[i,:])])

