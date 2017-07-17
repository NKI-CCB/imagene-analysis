def plot_ds(ds, fdr, le_prop=0.0, abs=True):
    ds = ds.copy()
    ds['significance_mask'] = (ds['fdr'] > fdr) | (ds['le_prop'] < le_prop)
    ds = ds.sel(
        gene_set=np.logical_not(ds['significance_mask'])
                 .sum('mri_feature') > 0,
    )
    if abs:
        zlim = [0, np.max(ds['nes'])]
        cmap='viridis'
    else:
        zlim = np.max(np.abs(ds['nes']))
        zlim = [-zlim, zlim]
        cmap = 'coolwarm'

    with plot.subplots(1, 1) as (fig, ax):
        plot.heatmap(
            ds['nes'], mask=ds['significance_mask'],
            zlim=zlim, cmap=cmap,
            row_dendrogram=True, col_dendrogram=True,
            ax=ax,
        )
        if len(ds['gene_set']) < 50:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        else:
            ax.set_xticklabels("")
