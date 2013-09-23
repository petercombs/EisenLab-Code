import numpy as np
from scipy import interpolate

import collections
import functools
import itertools

class memoized(object):
   '''Decorator. Caches a function's return value each time it is called.
   If called later with the same arguments, the cached value is returned
   (not reevaluated).
   '''
   def __init__(self, func):
      self.func = func
      self.cache = {}
   def __call__(self, *args):
      if not isinstance(args, collections.Hashable):
         # uncacheable. a list, for instance.
         # better to not cache than blow up.
         return self.func(*args)
      if tuple(*args) in self.cache:
         return self.cache[tuple(*args)]
      else:
         value = self.func(*args)
         self.cache[tuple(*args)] = value
         return value
   def __repr__(self):
      '''Return the function's docstring.'''
      return self.func.__doc__
   def __get__(self, obj, objtype):
      '''Support instance methods.'''
      return functools.partial(self.__call__, obj)

@memoized
def convert_to_distribution(points):
    """Converts the given data points to a smoothed distribution from 0-100%

    """

    x = np.linspace(0,1, len(points), endpoint=True)
    f = interpolate.interp1d(x, points, kind='cubic')
    retval =  np.cumsum(f(np.linspace(0, 1, 30, endpoint=True)).clip(0,1e5))
    return retval / (sum(retval)+1e-10)

def diff_stat(points1, points2):

    dist1 = convert_to_distribution(points1)
    dist2 = convert_to_distribution(points2)


    normfac = np.log(max(max(points1), max(points2)) + 1)

    return np.max(np.abs(dist1 - dist2)) * normfac

divmat = np.zeros([0,0])

def tang_stat(points1, points2):
    assert len(points1) == len(points2)
    points1 = np.array(points1 / np.mean(points1))
    points2 = np.array(points2 / np.mean(points2))

    va = np.reshape(np.repeat(points1, len(points2)), (len(points2), -1),
                    order='C')
    vb = np.reshape(np.repeat(points2, len(points1)), (-1, len(points2)),
                    order='F')

    global divmat
    if np.shape(divmat) != (len(points1), len(points2)):
        x, y = np.mgrid[0:len(points1), 0:len(points2)]
        divmat = 1/(np.abs(x - y) + 1)
    return np.sqrt(np.sum(np.triu((va - vb)**2 * divmat)))

#    stat = 0
#    for i in range(len(points1)):
#        for j in range(len(points2)):
#            stat += (points1[i] - points2[j])**2 / (np.abs(i - j)+1)
#            if i == j:
#                break
#
#    return np.sqrt(stat)



import progressbar as pb
def pdist(X, metric, p=2, w=None, V=None, VI=None):
    X = np.asarray(X, order='c')


    s = X.shape
    if len(s) != 2:
        raise ValueError('A 2-dimensional array must be passed.')

    m, n = s
    dm = np.zeros((m * (m - 1) / 2,), dtype=np.double)



    k = 0
    prog = pb.ProgressBar(widgets=['calculating distances', pb.Bar(),
                                   pb.Percentage(), pb.ETA()])
    for i in prog(range(0, m - 1)):
        for j in xrange(i + 1, m):
            dm[k] = metric(X[i], X[j])
            k = k + 1
    return dm

if __name__ == "__main__":
    import pandas as pd
    import matplotlib.pyplot as mpl

    zld_exp = pd.read_table('zeldaSummary.tsv', index_col=0)
    wt_exp = pd.read_table('wt_summary.tsv', index_col=0)

    zld_bind = pd.read_table('journal.pgen.1002266.s005.xls', skiprows=1)
    zld_bind.TSS_gene = zld_bind.TSS_gene.apply(str.strip)
    by_gene = zld_bind.groupby('TSS_gene')


    types = {'Intergenic':'N', 'Intronic':'I', 'Promoter':'P', 'UTR5':'5',
             'CDS':'C', 'UTR3':'3'}

    zld_comp = zld_exp.select(lambda x: 'cyc14A' in x, axis=1)
    wt_comp = wt_exp.select(lambda x: 'cyc14A' in x, axis=1)

    diff_col = pd.Series(index=zld_comp.index)
    for gene in wt_exp.index:
        assert gene in zld_exp.index
        diff_col[gene] = diff_stat(zld_comp.ix[gene], wt_comp.ix[gene])

    zld_exp['diff_col'] = diff_col
    wt_exp['diff_col'] = diff_col

    zld_exp.sort(column='diff_col', ascending=False, inplace=True)
    wt_exp.sort(column='diff_col', ascending=False, inplace=True)

    zld_fig_genes = zld_exp.select(lambda x: '14A' in x or '11' in x, axis=1)
    wt_fig_genes = wt_exp.select(lambda x: '14A' in x or '11' in x, axis=1)

    zld_fig_genes = zld_fig_genes[wt_fig_genes.max(axis=1) > 10][:120]
    wt_fig_genes = wt_fig_genes[wt_fig_genes.max(axis=1) > 10][:120]

    assert (zld_fig_genes.index == wt_fig_genes.index).all()

    mpl.figure(figsize=(2,8))
    mpl.subplot(1,2,1)
    mpl.pcolormesh(wt_fig_genes.divide(wt_fig_genes.max(axis=1)+1, axis=0).as_matrix(), cmap = mpl.cm.Blues)
    mpl.yticks(np.arange(len(zld_fig_genes.index)), zld_fig_genes.index)
    mpl.xlim(0, len(wt_fig_genes.columns))
    mpl.xticks([])


    mpl.subplot(1,2,2)
    mpl.pcolormesh(zld_fig_genes.divide(zld_fig_genes.max(axis=1)+1, axis=0).as_matrix(), cmap=mpl.cm.Reds)
    mpl.xlim(0, len(zld_fig_genes.columns))
    mpl.xticks([])
    ticktypes = [(i+0.5, ''.join(types[t]
                             for t in by_gene.get_group(gene).Label.unique()))
                 for i, gene in enumerate(zld_fig_genes.index)
                 if gene in by_gene.indices]
    mpl.yticks(*zip(*ticktypes))
    mpl.tight_layout()

