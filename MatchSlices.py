from __future__ import print_function, division
import PointClouds as pc
import sys
import numpy as np
from argparse import ArgumentParser
from scipy import stats

def parse_args():
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
                          help='Slices along the x axis in BDTNP data (DEFAULT)')
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

    args = argparser.parse_args()
    print(args)
    return args


def get_gene_names(vpc_reader):
    defined_names = ('id', 'x', 'y', 'z', 'Nx', 'Ny', 'Nz')
    names = set(name.split('_')[0] for name in vpc_reader.column
                if not name.startswith(defined_names))
    return names

def data_to_arrays(all_data, columns, genes):
    times = set(name.split('_')[-1] for name in columns if name != 'id')
    exparray = np.zeros((len(all_data), len(genes), len(times)))
    for j, gene in enumerate(genes):
        for k, time in enumerate(times):
            try:
                colnum = columns.index(gene + "__" + time)
                for i, row in enumerate(all_data):
                    exparray[i, j, k] = row[colnum]
            except ValueError:
                # No data for this gene at this time!
                pass

    posarray = np.zeros([len(all_data), 3, len(times)])
    for k, time in enumerate(times):
        for j, dim in enumerate(['x', 'y', 'z']):
            colnum = columns.index(dim + '__' + time)
            for i, row in enumerate(all_data):
                posarray[i,j,k] = row[colnum]

    return exparray, posarray

def virtual_slice(exparray, posarray, axis='x', width=50.0, resolution=1.0,
                  reduce_fcn = np.sum):
    # Convert axis specification to a usable value.
    axischooser = {'x':0, 0:0, 'y':1, 1:1, 'z':2, 2:2}
    axis = axischooser[axis]

    # Start and stop of slices
    start = posarray[:,axis,:].min() - width
    stop = posarray[:,axis,:].max() + width

    slicestarts = np.arange(start, stop, resolution)
    nslices = len(slicestarts)
    datasize = np.shape(exparray)
    ngenes = datasize[1]

    if len(datasize) < 3:
        exparray = np.reshape(exparray, (datasize[0], ngenes, 1))
        ntimes = 1
    else:
        ntimes = datasize[2]

    allslices = np.zeros((nslices, ngenes, ntimes))
    for k in range(ntimes):
        print("\nTime %d: " % k, end='')
        next_print = 0
        sys.stdout.flush()
        for i, pos in enumerate(slicestarts):
            expr = reduce_fcn(exparray[(pos <= posarray[:,axis,k])
                                   * (posarray[:,axis,k] < pos + width)],
                          axis=0)
            allslices[i,:] = expr
            if i / nslices > next_print:
                print('.', end='')
                sys.stdout.flush()
                next_print += .05
    print ()
    return slicestarts, allslices


def parse_translation_table(translation_table):
    fbgn2name = {}
    name2fbgn = {}
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

    gene_names = get_gene_names(bdtnp_parser)
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
        exparray, posarray = data_to_arrays(all_data, bdtnp_parser.column, gn_list)

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
            name = linedat[0]
            if ((name not in gene_names) and
                (name not in fbgn2name or (fbgn2name[name] not in gene_names))):
                continue
            fpkm = float(linedat[7])


            expr[fbgn2name[name] if name in fbgn2name else name] = fpkm

        expr_l = [(expr[gene] if gene in expr else 0) for gene in gn_list]
        for time in range(ntimes):
            bestcorr = 0
            bestslice = -1
            for slice in range(nslices):
                #corr = np.abs(np.corrcoef(expr_l, slices[slice, :, time]))[1,0]
                #corr = stats.pearsonr(expr_l, slices[slice, :, time])[0]
                corr = stats.spearmanr(expr_l, slices[slice, :, time])[0]
                if corr > bestcorr:
                    bestcorr = corr
                    bestslice = slice
            print("At time %f, best at #%d (%f), r=%f"
                  % (time, bestslice, starts[bestslice], bestcorr))

