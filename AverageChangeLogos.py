from Literature import all_hits
import PlotUtils as pu

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
