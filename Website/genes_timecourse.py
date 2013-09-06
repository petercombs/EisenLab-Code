#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import cgi
import sys
import pandas as pd
from os import path, environ
from subprocess import Popen

data = pd.read_table('genes.cuff', index_col=0)
procs = []

print "Content-Type: text/html"     # HTML is following
print                               # blank line, end of headers
print "<html><head>"
print "<TITLE>Gene Expression Across Time</TITLE>"
print "</head><body>"

print "<table>"
print "<thead align='center'><tr><td>Gene</td><td>Expression Timecourse</td>"
print "<td>Max FPKM</td></tr></thead>"


genes = [gene.strip() for gene in cgi.FieldStorage().getfirst('genes').split()
         if gene.strip()]
no_imgs = [gene for gene in genes
           if not path.exists(path.join('imgs', gene+'.png.svg'))
           and gene in data.index]
outfh = open('searches.log', 'a')
outfh.write(str(environ['REMOTE_ADDR']))
outfh.write(str(genes))
outfh.write('\n')
outfh.close()
if no_imgs:
    procs.append(Popen(['./draw_to_gene.py']
                       + no_imgs))

print """<tr><td></td>
<td><embed src="header.svg" width="750" height="50" /></td>
<td></td></tr>
"""

gene_index = {}
for gene in data.index:
    max_val = max(data.ix[gene])
    if ',' in gene:
        for gn in gene.split(','):
            gene_index[gn] = max_val
    gene_index[gene] = max_val

for gene in  genes:
    print "<tr><td>"
    print gene, 
    print '</td><td align="center">'
    if gene in gene_index:
        img = ''.join((l + '_' if l.isupper() else l) for l in gene)

        print '<embed src="'+path.join('imgs', img + '.png.svg')+'"'
        print 'width="750" height="50" />' 
        print '</td>'
        print '<td>%s</td>' % gene_index[gene]
    else:
        print 'No gene by that name</td><td>N/A</td>'
    print '</tr>'

print "</table>"
print """
Data is shown with Anterior to the left and Posterior to the right. Heatmaps are
normalized to maximum expression in any timepoint for each gene. Hover to
display Cufflinks' FPKM value in each slice datapoint.
"""
print "</body>"

print "Data calculated based on <a href="http://flybase.org/">FlyBase</a> files: "
print open("versions.txt").read()

print "</html>"
for proc in procs:
    proc.wait()
