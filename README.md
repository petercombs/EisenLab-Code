Eisen Lab Code
==============

This is a collection of code that I've written to do various tasks in lab.  I
make no claims to suitability for any purpose, and all code (unless otherwise
noted) is released under the CRAPL v0 license.  Please contact me direclty
(peter.combs@berkeley.edu) for any data I've generated, as it's most likely too
large to fit on github anyways.

The Fall2012 and Spring2013 branches are used for RNA-seq analysis of sliced
single _Drosophila melanogaster_ embryos, in particular for the submission of
the paper "Sequencing mRNA from cryo-siced _Drosophila_ embryos to determine
genome-wide spatial patterns of gene expression".  With the right configuration
and data files[\*], everything from raw reads to final figures should be able
to be accomplished with a run of

    $ ./configure
	$ make

The data for that paper is available at the Gene Expression Omnibus, under
accession GSE43506

[\*] As much as possible, these data files will be publicly available, standardized sets.  Known dependencies are:

 * FlyBase FASTA and GFF files for all species.  I believe they have to be unzipped. 
 * `journal.pbio.1000590.s002.txt`, the supplementary data file from Lott, et al 2011
 * `RunConfig.cfg` A tab-delimited file indicating, for each sample, the carrier
   species and sequencing index, among other statistics.  Please contact me for
   my copy if there's any trouble
 * `fig2_list.txt` A list of genes for making the table comparing FlyExpress
   thumbnails to the sliced expression patterns.

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

BayesMatch.py
--------------
This takes RNA-seq data from slices of drosophila embryos and attempts to find
the location in a BDTNP virtual embryo that has the best (Speaman) correlation
to that slice.


FB2name.py
----------

This is a tool to convert FlyBase gene identifiers (e.g. FBgn0000166 to bcd).
It generally assumes that data will have a columnar format, so you give it the
column (0 indexed)  that the FlyBase identifier is in with the -i flag, then
the other columns you want to output as well with individual -k flags.  For
example, to convert everything in gene_exp.diff, then take columns 6 and 7:

	python FB2name.py -r Reference/dmelfbgns.txt -i 0 -k 6 -k 7 gene_exp.diff

Use `python FB2name.py -h` for full options
