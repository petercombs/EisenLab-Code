from matplotlib.pyplot import (figure, plot, title, yticks, xticks, savefig,
                               subplot)
from numpy import linspace, shape, all
import pandas as pd

startswith = lambda y: lambda x: x.startswith(y)

wt_exp = pd.read_table('prereqs/WT5.53_summary.tsv',
                       na_values='-', keep_default_na=False)
zld_exp = pd.read_table('analysis/summary_merged.tsv',
                       na_values='-', keep_default_na=False)

wt_exp.sort_index(inplace=True)
zld_exp.sort_index(inplace=True)
assert(all(wt_exp.index == zld_exp.index))

gl = {line.strip() for line in open('prereqs/mikezyg.txt')}
figure()
wt_embs = set(column.split('_sl')[0] for column in wt_exp.columns)
for i, emb in enumerate(sorted(wt_embs)):
    subplot(2, len(wt_embs), i + 1)
    cd = wt_exp.select(startswith(emb), axis=1).ix[gl]
    old_shape = shape(cd)
    cd = cd.divide(cd.max(axis=1), axis=0)
    assert shape(cd) == old_shape
    plot(cd.median(axis=0))
    title('CaS\n' + emb)
    yticks([])
    xt = xticks()[0]
    xticks(linspace(0, xt[-1], 5, endpoint=True),
           ["%2d" % (t) for t in linspace(0, 100, 5, endpoint=True)[::-1]])

zld_embs = set(column.split('_sl')[0] for column in zld_exp.columns)

for i, emb in enumerate(sorted(zld_embs)):
    subplot(2, len(zld_embs), len(zld_embs) + i + 1)
    cd = zld_exp.select(startswith(emb), axis=1).ix[gl]
    old_shape = shape(cd)
    cd = cd.divide(cd.max(axis=1), axis=0)
    assert shape(cd) == old_shape
    plot(cd.median(axis=0))
    title("MZld-\n" + emb)
    yticks([])
    xt = xticks()[0]
    xticks(linspace(0, xt[-1], 5, endpoint=True),
           ["%2d" % (t) for t in linspace(0, 100, 5, endpoint=True)[::-1]])


savefig('analysis/results/zygotic_anteriorized.png')
