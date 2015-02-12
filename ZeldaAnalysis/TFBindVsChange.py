from __future__ import division,print_function
from maptlotlib import pyplot as mpl
import numpy as np
import pandas as pd


try:
    np.shape(gene_tab)
except NameError:
    gene_tab = pd.read_table('prereqs/gene_map_table_fb_2013_04.tsv',
                             skiprows=5,
                             skipfooter=3,
                             index_col='##current_symbol',
                            ).dropna(how='all')

def main():
    for current_tf in 'bcd hb'.split():




if __name__ == "__main__":
    main()
