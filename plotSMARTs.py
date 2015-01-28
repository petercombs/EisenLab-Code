from matplotlib.pyplot import (figure, plot, xlabel, ylabel, savefig, xlim,
                               ylim,gca, clf, legend, semilogx, hist, xticks)
import setcolor
from numpy import arange, std
from scipy.stats import scoreatpercentile, bartlett

try:
    outdir = locals()['outdir']
    all_samples = locals()['all_samples']

except KeyError:
    print("Need to run AnalyzeFoldChange.py first")
    assert False

print(outdir)

x_truseq = [0, 5, 10, 20]
x_values_few = [0, 50, 100, 200]
x_values_many = [0, 10, 50, 100, 200]

print("Done with setup")

figure()
p = plot(x_truseq, all_samples['Truseq'].T, 'k.-', alpha=0.01)
xticks(x_truseq)
xlabel('Concentration of D. virilis (pg)')
ylabel('Predicted amount of D. virilis (pg)')
#setcolor.set_foregroundcolor(gca(), 'w')
#setcolor.set_backgroundcolor(gca(), 'k')
ylim(0, scoreatpercentile(all_samples['Truseq'], 99))
gca().set_aspect(1)
savefig('{outdir}/TruSeqLines.png'.format(outdir=outdir),
        transparent=True,
        dpi=300)

print("Done with TruSeq Lines")


figure()
p = plot(x_values_few, all_samples['SMART2'].T, 'r', alpha=0.01)
xticks(x_values_few)
xlabel('Concentration of D. virilis (pg)')
ylabel('Predicted amount of D. virilis (pg)')
ylim(0, scoreatpercentile(all_samples['SMART2'], 99))
#setcolor.set_foregroundcolor(gca(), 'w')
#setcolor.set_backgroundcolor(gca(), 'k')
gca().set_aspect(1)
savefig('{outdir}/SMART2_Lines.png'.format(outdir=outdir), transparent=True, dpi=300)
print("Done with SMART2 Lines")


clf()
p = plot(x_values_many, all_samples['SMART2dil2'].T, 'g', alpha=0.01)
xticks(x_values_many)
ylim(0, scoreatpercentile(all_samples['SMART2dil2'], 99))
xlabel('Concentration of D. virilis (pg)')
ylabel('Predicted amount of D. virilis (pg)')
#setcolor.set_foregroundcolor(gca(), 'w')
#setcolor.set_backgroundcolor(gca(), 'k')
gca().set_aspect(1)
savefig('{outdir}/SMART2-dil2Lines.png'.format(outdir=outdir), transparent=True, dpi=300)
print("Done with 2.5x dilution Lines")


clf()
p = plot(x_values_many, all_samples['SMART2dil5'].T, 'c', alpha=0.01)
xticks(x_values_many)
ylim(0, scoreatpercentile(all_samples['SMART2dil5'], 99))
xlabel('Concentration of D. virilis (pg)')
ylabel('Predicted amount of D. virilis (pg)')
#setcolor.set_backgroundcolor(gca(), 'k')
#setcolor.set_foregroundcolor(gca(), 'w')
gca().set_aspect(1)
savefig('{outdir}/SMART2-dil5Lines.png'.format(outdir=outdir), transparent=True, dpi=300)
print("Done with 5x dilution Lines")

clf()
semilogx([1], [0]);
hist([samples.ix[r_values > .9, -1], samples.ix[r_values<=.9, -1]],
     bins=10**(arange(np.floor(np.log10(expr_min)), 4, .1)),
     normed=False, histtype='step', stacked=False, label=['r > .9', 'r < .9'])
legend()
xlabel('expression')
#setcolor.set_foregroundcolor(gca(), 'w')
#setcolor.set_backgroundcolor(gca(), 'k')
savefig('{outdir}/{protocol}-r90.png'.format(protocol=protocol,outdir=outdir), transparent=True, dpi=300)
print("Done with r=.90 histogram")

clf()
for c, p in zip('k r g c'.split(), 'Truseq SMART2 SMART2dil2 SMART2dil5'.split()):

    fpkms = arange(expr_min, 100, .1)
    pct5 = [scoreatpercentile(all_rs[p][all_samples_nonorm[p].ix[:,-1] > n], 5) for n in fpkms]
    semilogx(fpkms, pct5, c, label=p.replace('dil', '--'))
legend(loc='lower right')
xlabel('Raw Expression Cutoff')
ylabel('5th Percentile Correlation Coefficient')
#setcolor.set_foregroundcolor(gca(), 'w')
#setcolor.set_backgroundcolor(gca(), 'k')
savefig('{outdir}/r5th-SMARTs-{expmin}.png'.format(outdir=outdir, expmin=expr_min), transparent=True, dpi=300)
print("Done with cutoff estimation")

base = 'SMART2dil2'
base_slopes = all_slopes[base]

figure()
betweens = []
withins = []
for other in 'SMART2 SMART2dil5'.split():
    other_slopes = all_slopes[other]
    in_both = set(base_slopes.index).intersection(other_slopes.index)
    labelstr = ('Same pre-amplification'
                if 'dil' in other else
                'Different pre-amplification')
    plot(base_slopes.ix[in_both], other_slopes.ix[in_both], '.',
         label=labelstr)
    print("SMART2dil2 vs {}".format(other))
    withins.append(base_slopes.ix[in_both] + other_slopes.ix[in_both])
    betweens.append(base_slopes.ix[in_both] - other_slopes.ix[in_both])
    print("std within  ", std(withins[-1]))
    print("std between ", std(betweens[-1]))
    print(bartlett(withins[-1], betweens[-1]))
xlabel('Slope 1')
ylabel('Slope 2')

legend(numpoints=1)

xlim(.4, 1.5)
ylim(.4, 1.5)
#setcolor.set_foregroundcolor(gca(), 'w')
#setcolor.set_backgroundcolor(gca(), 'k')
savefig('{outdir}/noisesource_{expmin}.png'.format(outdir=outdir,
                                                   expmin=expr_min), transparent=True, dpi=300)
print("Done with noise sourcing")


