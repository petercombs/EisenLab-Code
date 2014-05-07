from matplotlib.pyplot import figure, plot, xlabel, ylabel, savefig, \
        xlim, ylim,gca, clf, legend, semilogx, hist
import setcolor
from numpy import arange
from scipy.stats import scoreatpercentile

print(outdir)

x_truseq = [0, 5, 10, 20]
x_values_few = [0, 50, 100, 200]
x_values_many = [0, 10, 50, 100, 200]

print("Done with setup")

figure()
p = plot(x_truseq, all_samples['Truseq'].T, 'w', alpha=0.01)
xlabel('ng D. virilis')
ylabel('normalized FPKM')
setcolor.set_foregroundcolor(gca(), 'w')
setcolor.set_backgroundcolor(gca(), 'k')
ylim(0, 40)
savefig('{outdir}/TruSeqLines.png'.format(outdir=outdir), transparent=True, dpi=300)

print("Done with TruSeq Lines")


figure()
p = plot(x_values_few, all_samples['SMART2'].T, 'r', alpha=0.01)
xlabel('pg D. virilis')
ylabel('normalized FPKM')
ylim(0, 400)
setcolor.set_foregroundcolor(gca(), 'w')
setcolor.set_backgroundcolor(gca(), 'k')
savefig('{outdir}/SMART2Lines.png'.format(outdir=outdir), transparent=True, dpi=300)
print("Done with SMART2 Lines")


clf()
p = plot(x_values_many, all_samples['SMART2dil2'].T, 'g', alpha=0.01)
ylim(0, 400)
xlabel('pg D. virilis')
ylabel('normalized FPKM')
setcolor.set_foregroundcolor(gca(), 'w')
setcolor.set_backgroundcolor(gca(), 'k')
savefig('{outdir}/SMART2-dil2Lines.png'.format(outdir=outdir), transparent=True, dpi=300)
print("Done with 2.5x dilution Lines")


clf()
p = plot(x_values_many, all_samples['SMART2dil5'].T, 'c', alpha=0.01)
ylim(0, 400)
xlabel('pg D. virilis')
ylabel('normalized FPKM')
setcolor.set_backgroundcolor(gca(), 'k')
setcolor.set_foregroundcolor(gca(), 'w')
savefig('{outdir}/SMART2-dil5Lines.png'.format(outdir=outdir), transparent=True, dpi=300)
print("Done with 5x dilution Lines")

clf()
semilogx([1], [0]);
hist([samples.ix[r_values > .9, -1], samples.ix[r_values<=.9, -1]], bins=10**(arange(-1, 4, .1)), normed=False, histtype='step', stacked=False, label=['r > .9', 'r < .9'])
legend()
xlabel('FPKM')
setcolor.set_foregroundcolor(gca(), 'w')
setcolor.set_backgroundcolor(gca(), 'k')
savefig('{outdir}/{protocol}-r90.png'.format(protocol=protocol,outdir=outdir), transparent=True, dpi=300)
print("Done with r=.90 histogram")

clf()
for c, p in zip('w r g c'.split(), 'Truseq SMART2 SMART2dil2 SMART2dil5'.split()):

    fpkms = arange(0, 100, .1)
    pct5 = [scoreatpercentile(all_rs[p][all_samples_nonorm[p].ix[:,-1] > n], 5) for n in fpkms]
    semilogx(fpkms, pct5, c, label=p.replace('dil', '--'))
legend(loc='lower right')
xlabel('FPKM Cutoff')
ylabel('5th Percentile Correlation Coefficient')
setcolor.set_foregroundcolor(gca(), 'w')
setcolor.set_backgroundcolor(gca(), 'k')
savefig('{outdir}/r5th-SMARTs.png'.format(outdir=outdir), transparent=True, dpi=300)
print("Done with cutoff estimation")

base = 'SMART2dil2'
base_slopes = all_slopes[base]

figure()
for other in 'SMART2 SMART2dil5'.split():
    other_slopes = all_slopes[other]
    in_both = set(base_slopes.index).intersection(other_slopes.index)
    plot(base_slopes.ix[in_both], other_slopes.ix[in_both], '.',
         label='{} vs {}'.format(base, other))

legend(numpoints=1)

xlim(.4, 1.5)
ylim(.4, 1.5)
setcolor.set_foregroundcolor(gca(), 'w')
setcolor.set_backgroundcolor(gca(), 'k')
savefig('{outdir}/noisesource.png'.format(outdir=outdir), transparent=True, dpi=300)
print("Done with cutoff estimation")


