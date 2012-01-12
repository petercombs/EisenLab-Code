import sys
from cPickle import load
from scipy import stats
from numpy import array, shape, sum, sqrt

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
                cutoff = sqrt(stats.chi2.ppf(.95, dof))
                rt_chi2_signed = (data - expected) / sqrt(expected)
                print '='*30
                print gene, cols, "TSSs", sum(data, axis=1), 'dof=', dof,
                print 'cutoff=', cutoff
                print data
                print (data - expected)/sqrt(expected)
                for sample in range(rows):
                    for tss in range(cols):
                        if rt_chi2_signed[sample, tss] > cutoff:
                            print "Hi: Sample %d TSS %d" % (sample, tss)
                        elif rt_chi2_signed[sample, tss] < -cutoff:
                            print "Lo: Sample %d TSS %d" % (sample, tss)
        except ValueError:
            print "Failed on: ", gene, data


