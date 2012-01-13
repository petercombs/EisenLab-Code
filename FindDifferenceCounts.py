import sys
from cPickle import load
from scipy import stats
from numpy import array, shape, sum, sqrt

if __name__ == "__main__":
    slices = []
    tss_data = False
    try:
        tss_data = load(open('inv_tss.pkl'))
    except IOError:
        pass
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
                for tss in range(cols):
                    if tss_data:
                        this_data = tss_data[gene, tss]
                    else:
                        this_data = False
                    high_in = []
                    low_in = []
                    for sample in range(rows):
                        if rt_chi2_signed[sample, tss] > cutoff:
                            high_in.append(sample)
                        elif rt_chi2_signed[sample, tss] < -cutoff:
                            low_in.append(sample)

                    print '# LO',
                    for item in this_data:
                        print item,
                    print str(low_in).replace(' ', '')
                    print '# HI',
                    for item in this_data:
                        print item,
                    print str(high_in).replace(' ', '')
        except ValueError:
            print "Failed on: ", gene, data


