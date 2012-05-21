import pysam
import sys
from os import path
from collections import defaultdict, Counter
from progressbar import ProgressBar, ETA, Bar, Percentage
from numpy import histogram

def get_nh(read):
    return [t[1] for t in read.tags if t[0] == 'NH'][0]

def get_species(read):
    return references[read.rname].split('_')[0]

def process_read(read):
    if not read.tags:
        print read
        print read.tags
        assert False
    nh = get_nh(read)
    species = get_species(read)
    if nh == 1:
        assigned.write(read)
        species_counts[species] += 1
        return
    else:
        resolve_multiread(read, nh, species)

def resolve_multiread(read, nh, species):
    nm = [t[1] for t in read.tags if t[0] == 'NM'][0]
    to_be_resolved_counts[read.qname] += 1
    if nm < to_be_resolved_vals[read.qname][species]:
        to_be_resolved_vals[read.qname][species] = nm
        to_be_resolved_reads[read.qname][species] = read

    if nh == to_be_resolved_counts[read.qname]:
        on_last_multiread(read)

def on_last_multiread(read):
    # Sort out the reads
    if len(to_be_resolved_reads[read.qname]) == 1:
        # Looks like there were multiple hits from the same species.
        # Report the best, or if equal quality, the first (which
        # tophat would've given anyways)
        assigned.write(read)
        species = references[read.rname].split('_')[0]
        species_counts[species] += 1
    else:
        # Hits from multiple species
        vals = sorted([(val, spec) for spec, val in
                        to_be_resolved_vals[read.qname].iteritems()])
        diff_val = vals[1][0] - vals[0][0]
        ambig_counts[diff_val] += 1
        if diff_val > ambig_threshold:
            assigned.write(to_be_resolved_reads[read.qname][vals[0][1]])
            species_counts[vals[0][1]] += 1
        else:
            species_counts['ambig'] += 1
            for amb_read in to_be_resolved_reads[read.qname].itervalues():
                ambig.write(amb_read)

    # Clean up the dictionaries
    to_be_resolved_vals.pop(read.qname)
    to_be_resolved_counts.pop(read.qname)
    to_be_resolved_reads.pop(read.qname)

ambig_threshold = 3

for fname in sys.argv[1:]:
    print fname
    samfile = pysam.Samfile(fname, 'rb')
    references = samfile.references
    dir = path.dirname(fname)
    assigned = pysam.Samfile(path.join(dir, 'assigned.bam'), 'wb',
                             template=samfile)
    ambig = pysam.Samfile(path.join(dir, 'ambiguous.bam'), 'wb',
                             template=samfile)

    to_be_resolved_reads = defaultdict(dict)
    to_be_resolved_vals = defaultdict(lambda : defaultdict(lambda : 1000))
    to_be_resolved_counts = Counter()
    species_counts = Counter()
    ambig_counts = Counter()

    start = samfile.tell()
    for read in samfile: pass
    maxval = samfile.tell()
    pbar = ProgressBar(maxval=maxval, widgets = [fname, ': ', Percentage(), ' ',
                                                Bar(), ' ', ETA(), ' '])
    samfile.seek(start)
    pbar.start()

    for read in samfile:
        pbar.update(samfile.tell() - start)
        process_read(read)
    print
    print "Species assignments: ", species_counts
    print "Ambiguity distribution: ", ambig_counts

            

