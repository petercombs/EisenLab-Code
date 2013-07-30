SliceSeq code
==============
__Zelda mutant analysis version__

This is a collection of code that has evolved, poorly, to do various generic
and specific analyses for SliceSeq, which at the time of writing, consists of
pooling a carrier RNA with the total RNA from each of the slices, then
computationally un-pooling the resulting RNA-seq reads. 

I make no claims to suitability for any purpose, and all code (unless otherwise
noted) is released under the CRAPL v0 license.  Please contact me direclty
(peter.combs@berkeley.edu) for any data I've generated, as it's most likely too
large to fit on github anyways.

I'm attempting to use a makefile as much as possible (See, for example, [this
blog
post](http://www.bioinformaticszen.com/post/decomplected-workflows-makefiles/)).
Thus, to do "everything", for some minimal definition of "everything", simply
run: $ make


Dependencies
------------

Software dependencies, at time of writing, include:
  * [STAR RNA-Seq](http://code.google.com/p/rna-star/) for mapping. Selected
    primarily for speed over Bowtie/Tophat, with no obvious major differences
in mapping results. Version 2.3.0.1
  * [Cufflinks](http://cufflinks.cbcb.umd.edu/) for abundance estimation (i.e.
    FPKM values for each gene).  Version 2.1.1

As much as possible, these data files will be publicly available, standardized
sets.  Known data dependencies are:

 * FASTA and GFF files for all species.  
 * [`journal.pbio.1000590.s002.txt`](http://www.plosbiology.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pbio.1000590.s002),
   the supplementary data file from Lott, et al 2011
 * `RunConfig.cfg` A tab-delimited file indicating, for each sample, the
   carrier species and sequencing index, among other statistics.  Please
contact me for my copy if there's any trouble

AssignReads2.py
---------------

Utility to demultiplex reads from a pooled RNA-seq sample.  Given a list of
`accepted_hits.bam` files from Tophat, this will assign the reads into bam
files for uniquely mapping to either species, or to an `ambiguous.bam` file,
all in the same directory as the original bam file.  There is a variable
`ambig_threshold` that determines what is counted as uniquely mapping, and is
currently set to 3, meaning that reads require 4 or more mismatches to prefer
one species over another.

CheckCoverage.py
----------------

Utility to estimate the relative level of PCR Duplication found in a sample.
By looking at the number of unique read positions in low-to-moderately
expressed genes, and comparing to a simulated model, a Badness score can be
calculated, which roughly corresponds to the number of reads per fragment.  A
perfect Badness score would be 1, with higher scores indicating a higher level
of PCR duplication in the sample. 

Usage:

    $ python CheckCoverage.py GTF-file bamfile.bam [bamfile.bam ...]

The GTF file works best when using a FlyBase derived file, and assumes the
following order of annotation types:

 * mRNA: Should have both the FBtr ID and the FBgn ID in the annotation field

 * exon: One or more exons per transcript, containingi the FBtr ID in the
  annotation field

 * CDS: Used by this program as a signal that there are no more exons for this
  transcript.  If there are, it will confuse the program.

PointClouds.py
--------------

This is a utility class for reading Point Cloud files from the Berkeley
Drosophila Transcription Network
(http://bdtnp.lbl.gov/Fly-Net/bioimaging.jsp?w=vpcFormat). Thus far, it's only
been tested on VirtualEmbryo files, but it should also work on single embryo
point clouds.  Metadata is loaded into the PointCloudReader object as
variables.  Example usage:

    from PointClouds import PointCloudReader
    pcr = PointCloudReader(open('D_mel.vpc'))
    bcddata = []
    for line in pcr:
        bcddata.append(line[pcr.column.index('bcd__1')])

which loads all of the data from the first timepoint with Bicoid into the
bcddata list.

