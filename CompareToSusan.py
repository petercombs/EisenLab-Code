from collections import defaultdict
from numpy import array, log2, size, shape, mean
from scipy import stats

susan_data_fname = '../journal.pbio.1000590.s002.txt'
susan_data = open(susan_data_fname)
susan_name_col = 0
susan_col_start = 6
susan_col_end = 30
susan_colnames = susan_data.readline().split('\t')[susan_col_start: susan_col_end]

my_data_fname = 'analysis/genes.fpkm_tracking'
my_data = open(my_data_fname)
my_name_col = 4
my_col_headers = my_data.readline().split()
my_cols, my_colnames = zip(*[(idx, colname)
                             for idx, colname in enumerate(my_col_headers)
                             if colname.endswith('_FPKM')])
my_colnames = list(my_colnames)
my_colnames.append('averaged')

susan_expr = {}
my_expr = {}

for line in susan_data:
    cols = line.split('\t')
    if len(cols) < susan_col_end:
        print line
        continue
    susan_expr[cols[susan_name_col]] = map(float, cols[susan_col_start:
                                                       susan_col_end])

for line in my_data:
    cols = line.split()
    my_expr[cols[my_name_col]] = map(float,
                                     [cols[i] for i in my_cols])
    my_expr[cols[my_name_col]].append(mean(my_expr[cols[my_name_col]]))

#print sorted([(max(susan_expr[key]), key) for key in susan_expr if key not in
              #my_expr])

keys = sorted([key for key in susan_expr if key in my_expr])
my_expr_arr = array([my_expr[key] for key in keys])
susan_expr_arr = array([susan_expr[key] for key in keys])

for my_column, my_colname in enumerate(my_colnames):
    print '-'*len(my_colname)
    print my_colname
    print '-'*len(my_colname)
    print '\tSprR\tPrsnR'
    my_sample = my_expr_arr[:, my_column]
    best = -1
    bestr = -1
    for susan_column, susan_colname in enumerate(susan_colnames):
        susan_sample = susan_expr_arr[:, susan_column]
        pearson_r, p = stats.pearsonr(my_sample, susan_sample)
        spearman_r, p = stats.spearmanr(my_sample, susan_sample)
        if spearman_r > bestr:
            best = susan_column
            bestr = spearman_r

        print '%s\t%4f\t%4f' % (susan_colname, spearman_r, pearson_r)

    print "Best hit: ", susan_colnames[best], bestr
#### Male vs Female identification

males = ['M14C', 'M14C_r2']
male_cols = [num for num, col in enumerate(susan_colnames) if col in males]
fems = ['F14C', 'F14C_r2']
fem_cols = [num for num, col in enumerate(susan_colnames) if col in fems]

num_rows, num_cols = shape(susan_expr_arr)
for i in range(num_rows):
    male_vals = sorted(susan_expr_arr[i,male_cols])
    fem_vals = sorted(susan_expr_arr[i,fem_cols])
    if abs(log2(mean(male_vals)/mean(fem_vals))) > 1:
        my_val = my_expr_arr[i,-1]
        if ((fem_vals[0] < male_vals[1] < fem_vals[1])
            or (male_vals[0] < fem_vals[1] < fem_vals[1])):
            #print "Ambiguous: ", keys[i]
            continue
        if male_vals[0] < my_val < male_vals[1] and male_vals[0] != 0:
            print keys[i], 'MALE', male_vals, my_val, fem_vals
        if fem_vals[0] < my_val < fem_vals[1] and fem_vals[0] != 0:
            print keys[i], 'FEM', fem_vals, my_val, male_vals
