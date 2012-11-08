from os import path

def get_bam_length(samfile):
    start = samfile.tell()
    maxval = path.getsize(samfile.filename) * 2**16 
    # I don't know why that 2**16 factor is there!
    return maxval +  2**16, start

def strip_to_number(dataval, chars = '\'" \t #'):
    return to_number(dataval.strip(chars))


def to_number(dataval):
    """ A forgiving number converter.

    Will convert to int if possible, float otherwise, and if neither, will return
    the input.
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

