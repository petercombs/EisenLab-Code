startswith = lambda x: lambda y: y.startswith(x)

for set_name, gene_set in (("antant", antant),
                           ("antlost", antlost),
                           ("antshift", antshift),
                           ("noant",noant)):
    axs = []
    max_y = min_y = 0
    figure()
    for i, emb in enumerate(embs):
        ax = subplot(1,4,i+1)
        axs.append(ax)
        n_wt = sum([n.startswith(emb) for n in wt_exp.columns])
        n_zld = sum([n.startswith(emb) for n in zld_exp.columns])
        y_wt = ( wt_exp.ix[gene_set]
                   .select(startswith(emb), axis=1)
                   .divide(norm_col.ix[gene_set], axis=0)
                   .mean(axis=0))
        y_zld = ( zld_exp.ix[gene_set]
                   .select(startswith(emb), axis=1)
                   .divide(norm_col.ix[gene_set], axis=0)
                   .mean(axis=0))
        title(emb)
        plot(linspace(0, 100, n_wt, True)[::-1], y_wt, 'b-')
        plot(linspace(0, 100, n_zld, True)[::-1], y_zld, 'r-')
        max_y = max(max_y, ax.get_ylim()[1])
        min_y = min(min_y, ax.get_ylim()[0])

    for ax in axs:
        ax.set_ylim(min_y, max_y)
        yticks = ax.get_yticks()
        ax.set_yticks([])
        ax.invert_xaxis()

    axs[0].set_yticks(yticks)

    draw_if_interactive()
    savefig('analysis/results/zldshift/'+set_name + ".png")


