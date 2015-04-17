#!/usr/bin/python
# -*- coding: utf-8 -*-

import cgi
import sys
import pandas as pd
from os import path, environ, tempnam
import sqlite3
import PlotUtils as pu
import time

tick = time.time()

startswith = lambda y: lambda x: x.startswith(y)
sel_startswith = lambda x: dict(crit=startswith(x), axis=1)

con = sqlite3.connect('summary.db')
index = pd.read_sql_query('SELECT "index" FROM all_expr',
                          con,
                          index_col='index').index

fbgns = {line.split()[1]: line.split()[0]
         for line in open('gene_table.tsv')
         if (not line.startswith('#')) and (line.strip())}
procs = []

print "Content-Type: text/html"     # HTML is following
print                               # blank line, end of headers
print "<html><head>"
print "<TITLE>Gene Expression Across Time</TITLE>"
print "</head><body>"

#print "<table>"
#print "<thead align='center'><tr><td>Gene</td><td>Expression Timecourse</td>"
#print "<td>Max FPKM</td></tr></thead>"


form = cgi.FieldStorage()
genes = [gene.strip() for gene in form.getfirst('genes').split()
         if gene.strip()]
good_genes = []
no_imgs = []
for i, gene in enumerate(genes):
    if gene.startswith('FBgn') and gene in fbgns:
        genes[i] = fbgns[gene]
    if genes[i] in index:
        good_genes.append(genes[i])
    else:
        no_imgs.append(genes[i])


if no_imgs:
    print "<p />Could not find: {}".format(','.join(no_imgs))

columns = tuple(form.getlist("samples"))
colnames = []
for column in columns:
    column = (column
              .replace('G20', '2.4x bicoid')
              .replace('bcd', '0x bicoid')
              .replace('hb', '0x hb')
              .replace('zld', '0x zld')
             )
    colnames.append(column)

data = pd.read_sql_query('SELECT * FROM all_expr WHERE ("index" in ("{}"))'
                         .format('", "'.join(genes)),
                         con, index_col='index')

data = (data
        .ix[good_genes]
        .fillna(value=pd.np.nan)
       )

data2 = tuple(data.select(**sel_startswith(column))
              for column in columns)


#print columns, data.columns
#print form.getfirst("format")
#print data.to_html()

is_svg = form.getfirst("format").lower() != 'png'
by_row = form.getfirst("orientation").lower() != 'genecol'

kwargs = dict(
    cmap_by_prefix=pu.cmap_by_prefix,
    col_sep='sl',
    norm_rows_by='max',
    total_width=min(1000/len(columns), 150),
    box_height=20,
    make_hyperlinks=True,
    convert = not is_svg,
    )

if by_row:
    outname = tempnam('imgs')
    kwargs['draw_box'] = True
    kwargs['draw_name'] = True
    kwargs['data_names'] = colnames
    kwargs['draw_row_labels'] = True
    pu.svg_heatmap(data2, filename=outname+'.svg', **kwargs)
    if is_svg:
        print '<embed src="{}" width="{}" height="{}">'.format(
            outname+'.svg',
            kwargs['total_width']*1.1*len(columns)+100,
            kwargs['box_height']*len(good_genes)+100,)
    else:
        print '<img src="{}">'.format(outname+'.png')


else:
    kwargs['total_width'] = 150
    kwargs['max_width'] = kwargs['total_width']
    kwargs['vspacer'] = 0
    for gene in good_genes:
        outname = tempnam('imgs', gene+'-')
        pu.svg_heatmap(tuple(d.ix[[gene]] for d in data2),
                       filename=outname+'.svg', **kwargs
                      )
        print "<table border=1><thead align=center><td>{}</td></thead>".format(gene)
        print "<tr><td>"
        if is_svg:
            print '<embed src="{}" width="150" height={}'.format(
                outname+'.svg',
                kwargs['box_height']*len(columns),
            )
        else:
            print '<img src="{}">'.format(outname+'.png')
        print "</td></tr>"
        print "</table><br>"

print """

<p />Data is shown with Anterior to the left and Posterior to the right. Heatmaps are
normalized to maximum expression in the given embryo. Slices with
poor quality data are masked with hashes, and colors are estimated as the
average of neighboring slices.

"""
if is_svg:
    print"""
<p />Hover to display Cufflinks' FPKM value in each slice datapoint.

"""
print "</body>"
'''
gene_index = {}
for gene in data.index:
    max_val = max(data.ix[gene].dropna())
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

        print '<embed src="'+path.join('imgs', img + '.svg')+'"'
        print 'width="770" height="50" />'
        print '</td>'
        print '<td>%s</td>' % gene_index[gene]
    else:
        print 'No gene by that name</td><td>N/A</td>'
    print '</tr>'

print "</table><br>"
'''
print '<p />Data calculated based on <a href="http://flybase.org/">FlyBase</a> files: '
print '<br><pre>'
print open("versions.txt").read()
print '<br></pre>'
'''
print """Underlying sequence expression data from: <br>
			<h3>
                <a href="http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0071820">
				Sequencing mRNA from cryo-sliced <i>Drosophila</i> embryos to determine genome-wide spatial patterns of gene expression
                </a>
			</h3>
			<p>
			Peter A. Combs, Michael B. Eisen
			</p>
		</center>
		<br>

		<p>
		<h3>Abstract</h3></p><p>
		Complex spatial and temporal patterns of gene expression underlie embryo differentiation, yet methods do not yet exist for the efficient genome-wide determination of spatial expression patterns during development. <em>In situ</em> imaging of transcripts and proteins is the gold-standard, but it is difficult and time consuming to apply to an entire genome, even when highly automated. Sequencing, in contrast, is fast and genome-wide, but is generally applied to homogenized tissues, thereby discarding spatial information. It is likely that these methods will ultimately converge, and we will be able to sequence RNAs in situ, simultaneously determining their identity and location. As a step along this path, we developed methods to cryosection individual blastoderm stage <em>Drosophila melanogaster</em> embryos along the anterior-posterior axis and sequence the mRNA isolated from each 25μm slice. The spatial patterns of gene expression we infer closely match patterns previously determined by<em> in situ</em> hybridization and microscopy. We applied this method to generate a genome-wide timecourse of spatial gene expression from shortly after fertilization through gastrulation. We identify numerous genes with spatial patterns that have not yet been described in the several ongoing systematic in situ based projects. This simple experiment demonstrates the potential for combining careful anatomical dissection with high-throughput sequencing to obtain spatially resolved gene expression on a genome-wide scale.
		</p>

"""

for proc in procs:
    proc.wait()


outfh = open('searches.log', 'a')
outfh.write(str(environ['REMOTE_ADDR']))
outfh.write(str(genes))
outfh.write('\n')
outfh.close()
'''
tock = time.time()
print "<hr>"
print "<p />Used {:.1f} ms".format((tock-tick)*1000)
print "</body>"
print "</html>"
