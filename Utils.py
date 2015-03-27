from os import path
from numpy import shape, linspace, sum, isfinite, copy


def get_bam_length(samfile):
    start = samfile.tell()
    maxval = path.getsize(samfile.filename) * 2**16
    # I don't know why that 2**16 factor is there!
    return maxval + 2**16, start


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

def center_of_mass(data):
    if 'columns' in dir(data) and 'rep' in data.columns[0]:
        reps = {c.split('_sl')[0] for c in data.columns}
        retval = 0
        for rep in reps:
            retval += center_of_mass_onerep(data.select(**sel_startswith(rep)))
        return retval / len(reps)
    elif 'index' in dir(data) and 'rep' in data.index[0]:
        reps = {c.split('_sl')[0] for c in data.index}
        retval = 0
        for rep in reps:
            retval += center_of_mass_onerep(data.select(startswith(rep)))
        return retval / len(reps)

    else:
        return center_of_mass_onerep(data)

def center_of_mass_onerep(data):
    dims = shape(data)
    cols = dims[-1]
    xs = linspace(0, 1, cols, endpoint=True)
    data_clean = data.copy()
    data_clean[~isfinite(data_clean)] = 0
    data_clean += 0.01
    return sum(data_clean * xs, axis=len(dims)-1)/sum(data_clean, axis=len(dims)-1)
