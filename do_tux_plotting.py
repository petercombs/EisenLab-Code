""" Possibly obsolete plotting system"""

# If we're on a system that supports it, do the plotting

try:
    from numpy import array
    from matplotlib import pyplot as mpl

    s1, s2, gene, idx = zip(*[(float(line.split()[6]), float(line.split()[7]),
                              line.split()[0], lnum)
                              for lnum, line
                              in enumerate(file(join(analysis_dir, 'gene_exp.diff')))
                             if 'FBgn' in line])

    s1 = array(s1)
    s2 = array(s2)

    FBgnToIDX = dict(zip(gene, idx))

    mpl.loglog(s1,s2,'k.', label="All Genes")
    for fname in glob('Genes*.txt'):
        label = fname[5:-4]
        genes_of_interest = [l.strip() for l in file(fname)]
        s1i = s1[[FBgnToIDX[NameKey[gene]] for gene in genes_of_interest
                 if gene in NameKey and NameKey[gene] in FBgnToIDX]]
        s2i = s2[[FBgnToIDX[NameKey[gene]] for gene in genes_of_interest
                 if gene in NameKey and NameKey[gene] in FBgnToIDX]]
        mpl.loglog(s1i, s2i, '.', label=label)
        print "-"*72, "\n", fname, "\n", "-"*72
        for g in genes_of_interest:
            if g.strip() in NameKey and NameKey[g.strip()] in genediff:
                print g, genediff[NameKey[g.strip()]]

    mpl.legend(numpoints=1, loc='lower right')
    mpl.loglog([1e-2,2e4], [1e-2,2e4], 'r:') # Diagonal line to guide eye


    # Clean up and label the axes
    ax = mpl.gca()
    ax.set_xlim(1e-2,1e+4)
    ax.set_ylim(1e-2,1e+4)
    ax.set_xlabel('Ant. Expr (RPKM)')
    ax.set_ylabel('Pos. Expr (RPKM)')
    mpl.savefig('LogLog.pdf')

    # If we have mutt, attach file to an email and send it off
    graphmail_proc = Popen(['mutt', '-s', 'LogLog', notificationEmail, '-a',
                            'LogLog.pdf'], stdin=open('to_email.tmp'))

except ImportError:
    print "Could not import Matplotlib.  You are using Python version",
    print sys.version
