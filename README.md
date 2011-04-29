Eisen Lab Code
==============

This is a collection of code that I've written to do various tasks in lab.  I
make no claims to suitability for any purpose, and all code is released under
the CRAPL v0 license.  Please contact me direclty (peter.combs@berkeley.edu)
for any data I've generated, as it's most likely too large to fit on github
anyways.

do_tux.py
---------

This is an automated script I first wrote to run the Tuxedo suite (tophat,
bowtie, and cufflinks; http://www.cbcb.umd.edu/software/) of bioinformatics
tools on a couple of RNA-seq runs that I did almost literally right before the
2011 Drosophila Conference.

All of the modifications to get it to work on a separate data set should be
below the first line of #'s, and all of the installation-specific software
calls below the second line of #'s. Some day I'll use optparse or argparse to
make these configurable options, but that day is not today.

Other than those lines, it expects the following:

 * Two FASTQ files corresponding to the reads (any name should work)

 * A number of files with names matching  Genes*.txt, containing lists of genes
   to plot on the loglog graph.

FB2name.py
----------

This is a tool to convert FlyBase gene identifiers (e.g. FBgn0000166 to bcd).
It generally assumes that data will have a columnar format, so you give it the
column (0 indexed)  that the FlyBase identifier is in with the -i flag, then
the other columns you want to output as well with individual -k flags.  For
example, to convert everything in gene_exp.diff, then take columns 6 and 7:

	python FB2name.py -r Reference/dmelfbgns.txt -i 0 -k 6 -k 7 gene_exp.diff

Use `python FB2name.py -h` for full options