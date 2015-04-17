import sqlite3
import Utils

Utils.load_to_locals(locals(), expr_min=0)
con = sqlite3.connect('analysis/summary.db')
all_expr.to_sql('all_expr', con, if_exists='replace')


