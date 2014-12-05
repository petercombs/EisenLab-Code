from os import path


def get_bam_length(samfile):
    start = samfile.tell()
    maxval = path.getsize(samfile.filename) * 2**16
    # I don't know why that 2**16 factor is there!
    return maxval + 2**16, start

def get_rfs_from_bam(bamfile):
    rfs = [entry for entry in
           bamfile.header['PG'][0]['CL'].split()
           if entry.endswith('.gz') or entry.endswith('.fastq')][0]
    return rfs.split(',')



def strip_to_number(dataval, chars='\'" \t #'):
    return to_number(dataval.strip(chars))


def to_number(dataval):
    """ A forgiving number converter.

    Will convert to int if possible, float otherwise, and if neither, will
    return the input.
    """
    try:
        datavalf = float(dataval)
        # If we could convert it to a float, it might have been an
        # int
        try:
            return int(dataval)
        except ValueError:
            # not an int, but since we got to the inner try, it is a
            # float
            return datavalf

    except ValueError:
        return dataval


def contains(string_or_iterable):
    if isinstance(string_or_iterable, str):
        return lambda x: string_or_iterable in x
    else:
        return lambda x: any(i in x for i in string_or_iterable)

def startswith(string_or_iterable):
    if not isinstance(string_or_iterable, str):
        string_or_iterable = tuple(string_or_iterable)
    return lambda x: x.startswith(string_or_iterable)

def sel_contains(string_or_iterable):
    return dict(crit=contains(string_or_iterable), axis=1)

def sel_startswith(string_or_iterable):
    return dict(crit=startswith(string_or_iterable), axis=1)
