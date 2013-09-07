#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python
# -*- coding: utf-8 -*-

import cgi
import sys
import pandas as pd
from os import path, environ
from subprocess import Popen

data = pd.read_table('genes.cuff', index_col=0)
fbgns = {line.split()[1]: line.split()[0] 
         for line in open('gene_table.tsv') 
         if (not line.startswith('#')) and (line.strip())}
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
for i, gene in enumerate(genes):
    if gene.startswith('FBgn') and gene in fbgns:
        genes[i] = fbgns[gene]

no_imgs = [gene for gene in genes
           if not path.exists(path.join('imgs', gene+'.png.svg'))
           and gene in data.index]
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

print 'Data calculated based on <a href="http://flybase.org/">FlyBase</a> files: '
print '<br><pre>'
print open("versions.txt").read()
print '<br></pre>'

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

print "</html>"
for proc in procs:
    proc.wait()


outfh = open('searches.log', 'a')
outfh.write(str(environ['REMOTE_ADDR']))
outfh.write(str(genes))
outfh.write('\n')
outfh.close()
