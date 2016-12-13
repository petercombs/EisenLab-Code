from __future__ import division
from Literature import all_hits
from Utils import sel_startswith
import PlotUtils as pu
import DistributionDifference as dd
import pandas as pd
from matplotlib.pyplot import subplot, pcolormesh, xlim, savefig, xticks, yticks

same_scale = dict(bcd=.7, zld=0.7, hb=0.7, G20=0.7)
diff_scale = dict(bcd=2, zld=1, hb=2, G20=2)

kwargs = dict(
    norm_rows_by=wts['cyc14D'].max(axis=1) + 10,
    progress_bar=False,
    col_sep='_sl',
    total_width=200,
    box_height=60,
    split_columns=True,
    draw_box=True,
    draw_row_labels=True,
    color_row_labels=all_hits,
    draw_name=False,
    cmap_by_prefix=pu.cmap_by_prefix,
    make_hyperlinks=True,
    convert=True,
    vspacer=0,
)

wtexp = wts['cyc14D']
wtcols = wtexp.dropna(how='all', axis=1).columns

for genotype in ['bcd', 'zld', 'hb', 'G20']:
    chg_ix = all_dists['WT_cyc14D-{}_cyc14D'.format(genotype)].sort(inplace=False).index[-100:]
    pu.svg_heatmap((wts['cyc14C'],
                    wts['cyc14D'],
                    bcds['cyc14D_rep1'],
                    bcds['cyc14D_rep2'],
                    g20s['cyc14D_rep1'],
                    g20s['cyc14D_rep2'],
                    hbs['cyc14D_rep1'],
                    hbs['cyc14D_rep2'],
                    zlds['cyc14A'],
                    zlds['cyc14B'],
                    zlds['cyc14D']),
                   'analysis/results/average_{}_diff.svg'.format(genotype),
                   index=chg_ix,
                   max_width=200,
                   draw_average_only=True,
                   average_scale=diff_scale[genotype],
                   **kwargs)

    from_array = pd.Series(data=pd.np.nan, index=wtcols)
    from_array.ix[wtcols] = 0
    genotype_expr = locals()[genotype.lower()+'s']['cyc14D']
    reps = {col.split('_sl')[0] for col in genotype_expr.columns}
    to_arrays = []
    for rep in reps:
        print(rep)
        rep_exp = genotype_expr.select(**sel_startswith(rep))
        rep_cols = rep_exp.dropna(how='all', axis=1).columns

        to_array = pd.Series(data=pd.np.nan, index=rep_cols)
        to_array.ix[rep_cols] = 0
        for gene in chg_ix:
            flows = dd.earth_mover(wtexp.ix[gene]+5,
                                   rep_exp.ix[gene]+5,
                                   True)
            for from_ix, to_ix, amt in flows:
                weight = abs(from_ix/(len(wtcols)-1) - to_ix/(len(rep_cols)-1))
                to_array[rep_cols[to_ix]] += amt * weight
                from_array[wtcols[from_ix]] += amt * weight / len(reps)
        to_arrays.append(to_array)

    subplot(1+len(reps), 1, 1)
    pcolormesh(from_array.as_matrix().reshape((1,-1)), cmap=pu.ISH)
    xlim(0, len(from_array))
    xticks([])
    yticks([])
    for i, rep_array in enumerate(to_arrays):
        subplot(1+len(reps), 1, 2+i)
        print(2+i)
        pcolormesh(rep_array.as_matrix().reshape((1,-1)),
                   cmap=pu.cmap_by_prefix(genotype))
        xlim(0, len(rep_array))
        xticks([])
        yticks([])
    savefig('analysis/results/average_flow_{}.png'.format(genotype))
    savefig('analysis/results/average_flow_{}.eps'.format(genotype))



for genotype in ['bcd', 'zld', 'hb', 'G20']:
    chg_ix = all_dists['WT_cyc14D-{}_cyc14D'.format(genotype)].sort(inplace=False).index[:500]
    pu.svg_heatmap((wts['cyc14C'],
                    wts['cyc14D'],
                    bcds['cyc14D_rep1'],
                    bcds['cyc14D_rep2'],
                    g20s['cyc14D_rep1'],
                    hbs['cyc14D_rep1'],
                    hbs['cyc14D_rep2'],
                    zlds['cyc14A'],
                    zlds['cyc14B'],
                    zlds['cyc14D']),
                   'analysis/results/average_{}_same.svg'.format(genotype),
                   index=chg_ix,
                   max_width=200,
                   draw_average_only=True,
                   average_scale=same_scale[genotype],
                   **kwargs)
