from __future__ import print_function
from matplotlib.pyplot import (figure, plot, xlabel, ylabel, savefig, xlim,
                               ylim,gca, clf, legend, semilogx, hist, xticks,
                               close, subplots)
import setcolor
import numpy as np
from numpy import arange, std
from scipy.stats import scoreatpercentile, bartlett

close('all')

try:
    outdir = locals()['outdir']
    all_samples = locals()['all_samples']
    all_samples_nonorm = locals()['all_samples_nonorm']
    all_slopes = locals()['all_slopes']
    samples = locals()['samples']
    expr_min = locals()['expr_min']

except KeyError:
    print("Need to run AnalyzeFoldChange.py first")
    assert False

print(outdir)

x_truseq = [0, 5, 10, 20]
x_values_few = [0, 50, 100, 200]
x_values_many = [0, 10, 50, 100, 200]

print("Done with setup")

for protocol, x_vals in [['Truseq', x_truseq],
                         ['EpicentreOligodT', x_values_few],
                         ['Clontech', x_values_few],
                         ['SMART2', x_values_few],
                         ['SMART2dil2', x_values_many],
                         ['SMART2dil5', x_values_many],]:
    print('-'*40, '\n', protocol, '\n','-'*40)
    clf()
    p = plot(x_vals, all_samples[protocol].T, 'k.-', alpha=0.01)
    xticks(x_vals)
    units = 'ng' if protocol == 'Truseq' else 'pg'
    xlabel('Concentration of D. virilis ({})'.format(units))
    ylabel('Predicted amount of D. virilis ({})'.format(units))
    ylim(0, scoreatpercentile(all_samples[protocol].ix[:,-1], 99))
    gca().set_aspect(1)
    savefig('{outdir}/{protocol}_Lines.png'.format(outdir=outdir,
                                                   protocol=protocol),
            transparent=True,
            dpi=300,
           )

    clf()
    f, (ax1, ax2) = subplots(2,1,sharex=True)
    max_n = 0
    for i, x in enumerate(x_vals):
        for ax in (ax1, ax2):
            n, bins, patches = ax.hist(all_samples[protocol].ix[:,i],
                    bins=arange(0,
                                scoreatpercentile(all_samples[protocol].ix[:,-1], 99),
                                max(x_vals)/100.),
                    histtype='step',
                    label='{} ng'.format(x),
                   )
            max_n = max(max_n, max(n))
    ax2.set_ylim(0, 250)
    ax1.set_ylim(300, 1.1*max_n)
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax1.xaxis.tick_top()
    ax1.tick_params(labeltop='off') # don't put tick labels at the top
    ax2.xaxis.tick_bottom()
    d = .015 # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d,+d),(-d,+d), **kwargs)      # top-left diagonal
    ax1.plot((1-d,1+d),(-d,+d), **kwargs)    # top-right diagonal

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d,+d),(1-d,1+d), **kwargs)   # bottom-left diagonal
    ax2.plot((1-d,1+d),(1-d,1+d), **kwargs) # bottom-right diagonal


    ax1.legend(loc='upper right')
    xlabel('Normalized expression')
    ylabel('Number of genes')
    ax2.yaxis.set_label_coords(-.075, 1.05)
    savefig('{outdir}/{protocol}_Hists.png'.format(outdir=outdir,
                                                   protocol=protocol),
            transparent=True,
            dpi=300)


comment='''

figure()
p = plot(x_truseq, all_samples['Truseq'].T, 'k.-', alpha=0.01)
xticks(x_truseq)
xlabel('Concentration of D. virilis (pg)')
ylabel('Predicted amount of D. virilis (pg)')
#setcolor.set_foregroundcolor(gca(), 'w')
#setcolor.set_backgroundcolor(gca(), 'k')
ylim(0, scoreatpercentile(all_samples['Truseq'].ix[:,-1], 99))
gca().set_aspect(1)
savefig('{outdir}/TruSeqLines.png'.format(outdir=outdir),
        transparent=True,
        dpi=300)

clf()
for i,x in enumerate(x_truseq):
    hist(all_samples['Truseq'].ix[:,i],
         bins=arange(0, scoreatpercentile(all_samples['Truseq'].ix[:,-1], 99), .25),
         histtype='step',
         label='{} ng'.format(x),
        )
legend(loc='upper right')
xlabel('Normalized expression')
ylabel('Number of genes')
savefig('{outdir}/TruSeqHists.png'.format(outdir=outdir),
        transparent=True,
        dpi=300)


print("Done with TruSeq Lines")

clf()
p = plot(x_values_few, all_samples['EpicentreOligodT'].T, 'm', alpha=0.01)
xticks(x_values_few)
xlabel('Concentration of D. virilis (pg)')
ylabel('Predicted amount of D. virilis (pg)')
ylim(0, scoreatpercentile(all_samples['EpicentreOligodT'].ix[:,-1], 99))
#setcolor.set_foregroundcolor(gca(), 'w')
#setcolor.set_backgroundcolor(gca(), 'k')
gca().set_aspect(1)
savefig('{outdir}/EpicentreOligodT_Lines.png'.format(outdir=outdir), transparent=True, dpi=300)
print("Done with Epicentre Lines")

clf()
p = plot(x_values_few, all_samples['SMART2'].T, 'r', alpha=0.01)
xticks(x_values_few)
xlabel('Concentration of D. virilis (pg)')
ylabel('Predicted amount of D. virilis (pg)')
ylim(0, scoreatpercentile(all_samples['SMART2'].ix[:,-1], 99))
#setcolor.set_foregroundcolor(gca(), 'w')
#setcolor.set_backgroundcolor(gca(), 'k')
gca().set_aspect(1)
savefig('{outdir}/SMART2_Lines.png'.format(outdir=outdir), transparent=True, dpi=300)
print("Done with SMART2 Lines")

clf()
for i,x in enumerate(x_values_few):
    hist(all_samples['SMART2'].ix[:,i],
         bins=arange(0, scoreatpercentile(all_samples['SMART2'].ix[:,-1], 99), 2.5),
         histtype='step',
         label='{} ng'.format(x),
        )
legend(loc='upper right')
xlabel('Normalized expression')
ylabel('Number of genes')
savefig('{outdir}/SMART2_Hists.png'.format(outdir=outdir),
        transparent=True,
        dpi=300)

clf()
p = plot(x_values_many, all_samples['SMART2dil2'].T, 'g', alpha=0.01)
xticks(x_values_many)
ylim(0, scoreatpercentile(all_samples['SMART2dil2'].ix[:,-1], 99))
xlabel('Concentration of D. virilis (pg)')
ylabel('Predicted amount of D. virilis (pg)')
#setcolor.set_foregroundcolor(gca(), 'w')
#setcolor.set_backgroundcolor(gca(), 'k')
gca().set_aspect(1)
savefig('{outdir}/SMART2-dil2Lines.png'.format(outdir=outdir), transparent=True, dpi=300)
print("Done with 2.5x dilution Lines")

clf()
for i,x in enumerate(x_values_many):
    hist(all_samples['SMART2dil2'].ix[:,i],
         bins=arange(0, scoreatpercentile(all_samples['SMART2dil2'].ix[:,-1], 99), 2.5),
         histtype='step',
         label='{} ng'.format(x),
        )
legend(loc='upper right')
xlabel('Normalized expression')
ylabel('Number of genes')
savefig('{outdir}/SMART2-dil2Hists.png'.format(outdir=outdir),
        transparent=True,
        dpi=300)

clf()
p = plot(x_values_many, all_samples['SMART2dil5'].T, 'c', alpha=0.01)
xticks(x_values_many)
ylim(0, scoreatpercentile(all_samples['SMART2dil5'].ix[:,-1], 99))
xlabel('Concentration of D. virilis (pg)')
ylabel('Predicted amount of D. virilis (pg)')
#setcolor.set_backgroundcolor(gca(), 'k')
#setcolor.set_foregroundcolor(gca(), 'w')
gca().set_aspect(1)
savefig('{outdir}/SMART2-dil5Lines.png'.format(outdir=outdir), transparent=True, dpi=300)
print("Done with 5x dilution Lines")

clf()
for i,x in enumerate(x_values_many):
    hist(all_samples['SMART2dil5'].ix[:,i],
         bins=arange(0, scoreatpercentile(all_samples['SMART2dil5'].ix[:,-1], 99), 2.5),
         histtype='step',
         label='{} ng'.format(x),
        )
legend(loc='upper right')
xlabel('Normalized expression')
ylabel('Number of genes')
savefig('{outdir}/SMART2-dil5Hists.png'.format(outdir=outdir),
        transparent=True,
        dpi=300)

'''

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
xlabel('Slope of gene in SS2-2.5x')
ylabel('Slope of gene in SS2 or SS2-5x')

legend(numpoints=1, loc='lower left')

xlim(.4, 1.5)
ylim(.4, 1.5)
#setcolor.set_foregroundcolor(gca(), 'w')
#setcolor.set_backgroundcolor(gca(), 'k')
savefig('{outdir}/noisesource_{expmin}.png'.format(outdir=outdir,
                                                   expmin=expr_min), transparent=True, dpi=300)
print("Done with noise sourcing")


