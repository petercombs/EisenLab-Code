from numpy import mean, shape, ceil, log2, sign
import pandas as pd

def has_anterior_peak(data, fold=3, thresh=3):
    n = shape(data)[-1]
    if isinstance(data, pd.DataFrame):
        first_third = data.ix[:,:n//3].mean(axis=1)
        middle_third = data.ix[:,n//3:-(n//3)].mean(axis=1)
    else:
        first_third = mean(data[:n//3])
        middle_third = mean(data[n//3:-(n//3)])
    return ((sign(fold) * first_third > fold * middle_third) *
            first_third > thresh)

def has_posterior_peak(data, fold=3, thresh=3):
    n = shape(data)[-1]
    if isinstance(data, pd.DataFrame):
        first_third = data.ix[:,-n//3:].mean(axis=1)
        middle_third = data.ix[:,n//3:-(n//3)].mean(axis=1)
    else:
        first_third = mean(data[-n//3:])
        middle_third = mean(data[n//3:-(n//3)])

    return ((first_third > fold * middle_third) *
            (first_third > thresh))

def has_central_peak(data, fold=3, thresh=3):
    n = shape(data)[-1]
    if isinstance(data, pd.DataFrame):
        first_third = data.ix[:,:n//3].mean(axis=1)
        middle_third = data.ix[:,n//3:-(n//3)].mean(axis=1)
        last_third = data.ix[:,-(n//3):].mean(axis=1)
    else:
        first_third = mean(data[:n//3])
        middle_third = mean(data[n//3:-(n//3)])
        last_third = mean(data[-(n//3):])
    return ((middle_third > fold * first_third)*
            (middle_third > fold * last_third) *
            (middle_third > thresh))



def make_sort_num(*args, **kwargs):
    retval = 0
    pow = 0
    for func in (has_posterior_peak, has_central_peak, has_anterior_peak):
        for exp_set in args:
            if isinstance(exp_set, tuple):
                retval += 2**pow * sum(func(s, **kwargs) for s in exp_set)
                pow += int(ceil(log2(len(exp_set) + 1)))
            else:
                retval += 2**pow * func(exp_set, **kwargs)
                pow += 1
    return retval
