#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python
import cgi
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from StringIO import StringIO

startswith = lambda x: lambda y: y.startswith(x)
bad_slices = {'emb4_sl10_FPKM', 'emb6_sl01_FPKM', 'emb7_sl01_FPKM',
              'emb7_sl07_FPKM', 'emb7_sl08_FPKM'}

#### HTTP HEADER
print "Content-Type: image/png\n"

data = pd.read_table('genes.cuff', index_col=0).select(lambda x: x not in
                                                       bad_slices, axis=1)
gene = cgi.FieldStorage()['gene'].value.strip()

all_max = min(max(data.ix[gene]), 1000)

fig = plt.figure(figsize=(8,1))
for emb in range(1,8):
    ax = plt.subplot(1,7,emb)
    dat = data.ix[gene].select(startswith('emb%d' % emb))
    plt.pcolormesh(np.reshape(dat, (1, -1)),
                  cmap=plt.cm.Blues,
                  vmin=0, vmax=all_max)
    ax.set_xlim([0, len(dat)])
    ax.set_xticks([0, len(dat)])
    ax.set_xticklabels('AP')
    ax.set_yticks([])
    plt.tight_layout()

out = StringIO()
plt.savefig(out, format='png')
print out.getvalue()

