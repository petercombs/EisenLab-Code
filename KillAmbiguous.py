#!/usr/bin/env python

import pysam
from collections import defaultdict
from Utils import get_bam_length
from progressbar import ProgressBar, Percentage, Bar, ETA
import sys
from glob import glob
from os import path

# SAM TAG Table
is_multiple = 0x1           #template having multiple fragments in sequencing
is_aligned = 0x2            #each fragment properly aligned according to the aligner
is_unmapped = 0x4           #fragment unmapped
is_next_unmapped = 0x8      #next fragment in the template unmapped
is_revcomp = 0x10           #SEQ being reverse complemented
is_next_reversed = 0x20     #SEQ of the next fragment in the template being reversed
is_read1 = 0x40             #the first fragment in the template
is_read2 = 0x80             #the last fragment in the template
is_secondary = 0x100        #secondary alignment
is_failqc = 0x200           #not passing quality controls
is_dupe = 0x400             #PCR or optical duplicate

def reheader(in_sam, keepstr = 'dmel'):
    new_header = in_sam.header.copy()
    new_SQ = [sub_dict for sub_dict in new_header['SQ']
              if keepstr in sub_dict['SN']]
    new_header['SQ'] = new_SQ
    return new_header



#Operate on each file independently
for fname in sys.argv[1:]:
    all_files = glob(path.join(path.dirname(fname), 'assigned_d???.bam'))
    all_files.append(path.join(path.dirname(fname), 'ambiguous.bam'))
    print all_files
    print "Getting read names..."
    sys.stdout.flush()
    qnames_in_file = [{read.qname for read in pysam.Samfile(fn)}
                      for fn in all_files]
    print "All read names complete"
    sys.stdout.flush()

    failed_strict = pysam.Samfile(path.join(path.dirname(fname),
                                            'ambiguous_strict.bam'),
                                  'wb', header=pysam.Samfile(all_files[-1]))
    for i, fname in enumerate(all_files):
        data = defaultdict(lambda : [None, None])
        #if 'ambiguous' in fname:
        if 'dmel' not in fname:
            continue

        # Build a single set of qnames found in other files
        other_qnames = set().union(*qnames_in_file[:i])
        other_qnames.update(*qnames_in_file[i+1:])

        infile = pysam.Samfile(fname)
        outfile = pysam.Samfile(fname[:-4] + '_strict_unsorted.bam', 'wb',
                                header=reheader(infile, fname[-8:-4]))
        irefs = infile.references
        orefs = outfile.references
        maxval, start = get_bam_length(infile) # For progress bar goodness


        pbar = ProgressBar(maxval=maxval - start,
                           widgets = [fname, ': ', Percentage(), ' ', Bar(), ' ',
                                      ETA(), ' '])
        pbar.start()

        # Load up all the read pairs that mapped, but don't load in any read
        # pair that failed to map
        for read in infile:
            pbar.update(infile.tell() - start)
            qname = read.qname
            if qname in other_qnames:
                failed_strict.write(read)
                continue
            read.rname = orefs.index(irefs[read.rname])
            data[qname][read.is_read2] = read

            # When we have both ends of the read,
            # fix all the SAM tags
            if None not in data[qname]:
                read1, read2 = data.pop(qname)
                is_same = read1.rname == read2.rname
                read1.rnext = read2.rname
                read2.rnext = read1.rname
                read1.pnext = read2.pos
                read2.pnext = read1.pos
                dist = abs(read2.pos - read1.pos) if is_same else 0
                read1.tlen = dist
                read2.tlen = -dist
                read1.flag = (is_multiple + is_aligned +
                              read1.is_reverse * is_revcomp +
                              read2.is_reverse * is_next_reversed +
                              is_read1)
                read2.flag = (is_multiple + is_aligned +
                              read2.is_reverse * is_revcomp +
                              read1.is_reverse * is_next_reversed +
                              is_read2)
                assert read1.flag & 0x1
                outfile.write(read1)
                outfile.write(read2)

        pbar.finish()



        # Finish up singleton reads
        pbar = ProgressBar(widgets=['Finishing :', Percentage(), ' ', Bar(), ' ',
                                    ETA(), ' '])
        for key in pbar(data):
            read1, read2 = data[key]
            if read1 and read1.flag & is_unmapped:
                print read1
            if read2 and read2.flag & is_unmapped:
                print read

            if read1 is None and read2 is None:
                print "Somehow we ended up with a double empty"
                assert False
            elif read1 is not None and read2 is not None:
                print "Somehow we didn't clear out the double-map"
                assert False
            elif read2 is None:
                read1.rnext = -1
                read1.pnext = 0
                read1.tlen = 0
                read1.flag = read1.flag | is_next_unmapped
                outfile.write(read1)
            elif read1 is None:
                read2.rnext = -1
                read2.pnext = 0
                read2.tlen = 0
                read2.flag = read2.flag | is_next_unmapped
                outfile.write(read2)
            else:
                print "How did we get here?"
                assert False
        filename = outfile.filename
        outfile.close()

        sorted_filename = filename[:filename.index('_unsorted')]
        print "Sorting into", sorted_filename
        pysam.sort('-m', "%d" % 3e9, filename, sorted_filename)