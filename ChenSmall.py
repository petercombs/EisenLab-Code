import pandas as pd
import BindUtils as bu
import PlotUtils as pu
from Utils import sel_startswith, sel_contains, load_to_locals
from numpy import mean, linspace, nan
from progressbar import ProgressBar
from Literature import staller_2015

def get_chrom(coord_spec):
    return coord_spec.split(':')[0]

def com_calculator(the_series):
    new_series = the_series.replace(nan, 0)
    return sum(linspace(0, 1, len(the_series), endpoint=True)
            * new_series)/sum(new_series)


def get_com(expr_data):
    coms = pd.Series(index=expr_data.index, data=0)
    for gene in coms.index:
        coms.ix[gene] = com_calculator(expr_data.ix[gene])
    return coms

def get_average_pos(coord_spec):
    coords = coord_spec.split(':')[1]
    return mean([float(i) for i in coords.split('..')])

def get_max(data, set, stage):
    return data.ix[set].dropna(axis=1).select(**sel_contains(stage)).max(axis=1)


if __name__=="__main__":
    gff_cols = ['chr', 'src', 'type', 'start', 'stop', 'score', 'strand',
                'frame', 'attribute']
    chen_crms = pd.read_csv('prereqs/Chen_Bcd6.csv')
    chen_crms = pd.read_table('prereqs/chensmall_r6.gff',
                              names=gff_cols, header=None)

    read_table_args = dict(index_col=0,
                           keep_default_na=False,
                           na_values=['---', ''])
    expr_min = 5
    eps = 1
    read_table_args = dict(index_col=0,
                           keep_default_na=False,
                           na_values=['---', ''])

    if 'all_expr' not in locals():
        all_expr, (wt, bcd, zld, g20, hb), (wts, bcds, zlds, g20s, hbs), by_cycle\
        = load_to_locals(locals(), expr_min)
    print("Read expression in")

    dist = 10000
    min_expr = 10

    genes = {gene
             for c in chen_crms.index
             for gene in bu.find_near(bu.tss_dict[chen_crms.chr[c]],
                                      chen_crms.start[c],
                                      dist,
                                      nearest_only=True,
                                     )
            }

    hi_genes = {gene
                for gene in genes
                if (
                    (gene in wt.index)
                    and (max(wt.dropna(axis=1).ix[gene]) > min_expr)
                   )
               }
    print({gene:max(wt.dropna(axis=1).ix[gene]) for gene in hi_genes})

    print(genes)
    print(genes.difference(hi_genes))

    coms = get_com(wt.ix[hi_genes].select(**sel_contains('cyc14D')))
    coms.sort()
    hi_genes = coms.index

    data = (
        bcd.select(**sel_contains('cyc14D_rep1')).ix[hi_genes],
        bcd.select(**sel_contains('cyc14D_rep2')).ix[hi_genes],
        wt .select(**sel_contains('cyc14D')).ix[hi_genes],
        g20.select(**sel_contains('cyc14D_rep1')).ix[hi_genes],
        #g20.select(**sel_contains('cyc14D_rep2')).ix[hi_genes],
        #zld.select(**sel_contains('cyc14D')).ix[hi_genes],
        #hb.select(**sel_contains('cyc14D_rep1')).ix[hi_genes],
        #hb.select(**sel_contains('cyc14D_rep2')).ix[hi_genes],
    )
    #max_13s = max(float(dp.max(axis=1))*1.1 +10 for dp in data[0::2])
    #max_14s = max(float(dp.max(axis=1))*1.1 +10 for dp in data[1::2])
    #normer = ((max_13s, max_14s)*5)[:-1]
    normer = (
        get_max(bcd, hi_genes, 'cyc14D') * 1.1 + 10,
        get_max(bcd, hi_genes, 'cyc14D') * 1.1 + 10,
        get_max(wt, hi_genes, 'cyc14D' ) * 1.1 + 10,
        get_max(g20, hi_genes, 'cyc14D') * 1.1 + 10,
        #get_max(g20, hi_genes, 'cyc14D') * 1.1 + 10,
        #get_max(zld, hi_genes, 'cyc14D') * 1.1 + 10,
        #get_max(hb, hi_genes, 'cyc14D') * 1.1 + 10,
        #get_max(hb, hi_genes, 'cyc14D') * 1.1 + 10,
    )
    wt_pos = 2
    names = (
        'Bcd-14D Bcd-14D '
        'WT14D '
        #'G20-14D G20-14D '
        'G20-14D '
        #'Zld-14D '
        #'Hb-14D Hb-14D '
    )
    #names = [stage + '- Max Expr {:.1f}'.format(float(dp.max(axis=1)))
             #for stage, dp in zip(names.split(),
                                  #data)
            #]
    kwargs = dict(
        cmap_by_prefix=pu.cmap_by_prefix,
        col_sep='sl',
        draw_name=True,
        data_names=names.split(),
        draw_box=True,
        total_width=150,
        box_height=10,
        make_hyperlinks=True,
        draw_row_labels=True,
        color_row_labels=staller_2015,
        convert=True,
    )
    pu.svg_heatmap(data,
                   filename='analysis/results/chen.svg',
                   norm_rows_by=normer,
                   **kwargs
                  )

    pu.svg_heatmap(data,
                   filename='analysis/results/chen_wtnorm.svg',
                   norm_rows_by=normer[wt_pos],
                   **kwargs
                  )
    pb = ProgressBar()
    for gene in pb(hi_genes):
        kwargs['draw_box'] = False
        kwargs['draw_name'] = False
        kwargs['max_width'] = 150
        kwargs['box_height'] = 20
        kwargs['vspacer'] = 0

        pu.svg_heatmap(tuple(d.ix[[gene]] for d in data),
                       filename='analysis/results/chen/{}.svg'.format(gene),
                       norm_rows_by=tuple(n.ix[gene] for n in normer),
                       **kwargs
                      )
        pu.svg_heatmap(tuple(d.ix[[gene]] for d in data),
                       filename='analysis/results/chen_wtnorm/{}.svg'.format(gene),
                       norm_rows_by=normer[wt_pos].ix[[gene]],
                       **kwargs
                      )

'''
convert -background none -gravity center \( chen.png \) \( -mattecolor none
-size 10x10 xc:none \) \
             -bordercolor none -background none -mattecolor black -frame 0 \
             \( -gravity SouthEast -frame 2 -splice 10x10
chen/pngs/{noc,bowl,slp2,slp1}.png +append \) \
             \( -gravity SouthEast -frame 2 -splice 10x10
chen/pngs/{CG16813,prd,h,kni}.png +append \) \
             \( -gravity SouthEast -frame 2 -splice 10x10
chen/pngs/{CG16815,CG13713,Bsg25D}.png +append \) \
          -gravity Center -append chenfig.png
'''
