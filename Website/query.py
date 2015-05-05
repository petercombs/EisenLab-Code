#!/usr/bin/python
# -*- coding: utf-8 -*-

import cgi
import pandas as pd
from os import tempnam
import sqlite3
import PlotUtils as pu
from time import time
import warnings
import json
import sys

warnings.filterwarnings("ignore", message=".*tempnam.*")

print """Content-Type: text/html

<html><head>
<meta charset=utf-8>
<TITLE>Gene expression in Drosophila Mutants</TITLE>
</head><body>
"""
sys.stdout.flush()

tick = time()

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



form = cgi.FieldStorage()

input = dict(genes=form.getfirst('genes', "eve\nhb"),
             norm_by=form.getfirst('norm_by', "max"),
             columns=form.getlist('samples'),
             format=form.getfirst("format", "svg"),
             orientation=form.getfirst('orientation', "GeneRow"),
            )

if not input['columns']:
    input['columns'] = ['WT_cyc13']

try:
    input['avg_norm'] = float(form.getfirst('avg_norm'))
except TypeError:
    input['avg_norm'] = 1.0


genes = input['genes']
if genes:
    genes = [gene.strip() for gene in genes.split()
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
        print "<p />Could not find: {}".format(', '.join(no_imgs))

    if good_genes:
        columns = tuple(input['columns'])
        norm_by = input['norm_by'] or 'max'
        colnames = []
        for column in columns:
            column = (column
                      .replace('bcd', '0x bcd')
                      .replace('G20', '2.4x bcd')
                      .replace('hb', '0x hb')
                      .replace('zld', '0x zld')
                     )
            colnames.append(column)

        #print "<hr> Used {:0.1f} ms<hr>".format(1000*(time() - tick))

        data = pd.read_sql_query('SELECT * FROM all_expr WHERE ("index" in ("{}"))'
                                 .format('", "'.join(good_genes)),
                                 con, index_col='index')
        data = (data
                .ix[good_genes]
                .fillna(value=pd.np.nan)
               )

        if norm_by == 'max_all':
            norm_by = columns

        if norm_by != 'max':
            norm_by = data.select(**sel_startswith(norm_by)).max(axis=1) + 10
        else:
            # Something weird is going on with the type?
            norm_by = 'max'



        data2 = tuple(data.select(**sel_startswith(column))
                      for column in columns)



        is_svg = input["format"].lower() != 'png'
        by_row = input["orientation"].lower() != 'genecol'

        kwargs = dict(
            cmap_by_prefix=pu.cmap_by_prefix,
            col_sep='sl',
            norm_rows_by=norm_by,
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
            kwargs['draw_average'] = True
            kwargs['average_scale'] = input['avg_norm']
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
                if isinstance(norm_by, str) and norm_by == 'max':
                    kwargs['norm_rows_by'] = 'max'
                else:
                    kwargs['norm_rows_by'] = norm_by.ix[gene]
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

print '<p />Data calculated based on <a href="http://flybase.org/">FlyBase</a> files: '
print '<br><pre>'
print open("versions.txt").read()
print '<br></pre>'

'''
outfh = open('searches.log', 'a')
outfh.write(str(environ['REMOTE_ADDR']))
outfh.write(str(genes))
outfh.write('\n')
outfh.close()
'''

print open('form.html').read()


tock = time()
print "<hr>"
print """<script type="text/javascript">
var inputs={};
document.getElementById("genes").setRangeText(inputs.genes, 0, 1000000);
document.getElementById(inputs.format).checked=true;
document.getElementById(inputs.orientation).checked=true;
document.getElementById("norm_by").value=inputs.norm_by;
document.getElementById("avg_norm").value=inputs.avg_norm;
var samps = document.getElementsByName("samples")
//samps is not iterable???
for (var i=0; i < samps.length; i++){{
    var el = samps[i];
    if (inputs.columns.indexOf(el.value) >= 0){{
        el.checked=true;
    }}else{{
        el.checked=false;
    }}
}}

</script>
""".format(json.dumps(input))

print "<p />Used {:.1f} ms".format((tock-tick)*1000)
print "</body>"
print "</html>"
