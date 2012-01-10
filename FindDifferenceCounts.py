import sys
from cPickle import load
from scipy import stats

if __name__ == "__main__":
    slices = []
    for fname in sys.argv[1:]:
        genes.append(load(open(fname)))

    for gene in slices[0]:
        data = [slice[gene] for slice in slices]
        chi2, p, dof, expected = stats.chi2_contingency(data)
        if p < (.05 / len(slices[0])):
            print gene, data


