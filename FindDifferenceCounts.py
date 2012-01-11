import sys
from cPickle import load
from scipy import stats
from numpy import array, shape

if __name__ == "__main__":
    slices = []
    for fname in sys.argv[1:]:
        slices.append(load(open(fname)))

    for gene in slices[0]:
        data = [slice[gene] for slice in slices]
        n = max(len(i) for i in data)
        for i in data:
            while len(i) < n:
                i.append(0)

        data = array(data, dtype=float)
        rows, cols = shape(data)
        # Pseudocounts for fixing low-count data
        data += 2.0 / cols

        try:
            chi2, p, dof, expected = stats.chi2_contingency(data)
            if p < (.05 / (len(slices[0])* len(sys.argv[1:]))):
                print gene, data
        except ValueError:
            print "Failed on: ", gene, data


