""" A reader for the VirtualEmbryo File Format

This module provides classes that assist in the reading and writing of
VirtualEmbryo files from the Berkeley Drosophila Transcription Network Project
(http://bdtnp.lbl.gov/Fly-Net/bioimaging.jsp?w=vpcFormat).  It attempts to
preserve as much of the header information as possible in an easily interactible
way.
"""


class PointCloudReader:
    def __init__(self, fh):
        self.__filehandle__ = fh
        pos = fh.tell()
        line = fh.readline()
        while line[0] == '#':
            # Read the header info
            if line[1] == '#':
                #Double # == comment
                pos = fh.tell()
                line = fh.readline()
                continue
            
            metadata = line.split('=')
            dataname = metadata[0].strip(' \t#')
            dataval = metadata[1].strip()

            if ';' in dataval:
                # 2D array, ';' delimits rows, ',' delimits columns
                assert dataval[0] == '['
                assert dataval[-1] == ']'

                datarows = dataval[1:-1].split(';')
                datarowscols = list(map(lambda s: s.split(','), datarows))
                dataval = map(lambda r: list(map(strip_to_number, r)), datarowscols)
                dataval = list(dataval)

            elif dataval[0] == '[' and dataval[-1] == ']':
                # 1D array
                dataval = map(strip_to_number,
                              dataval[1:-1].split(','))
                dataval = list(dataval)

            else:
                # It's a scalar.
                dataval = to_number(dataval)

            self.__setattr__(dataname, dataval)
            pos = fh.tell()
            line = fh.readline()

        # Clean up afterwards
        self.__filehandle__.seek(pos)

    def __next__(self):
        line = self.__filehandle__.readline()
        if line == '':
            raise StopIteration
        return list(map(strip_to_number,
                        line.split(',')))
    
    def __iter__(self):
        return self


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

