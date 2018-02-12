import click
import matplotlib
import numpy as np
import pandas as pd
import xarray as xr

from lib import click_utils
import plot


def wf_plot(vals, highlight, ax, ylabel="", yscale="linear", ylim=None,
            xbaseline=None, reverse=False):
    if reverse:
        vals_order = np.argsort(-vals.values)
    else:
        vals_order = np.argsort(vals.values)

    vals = vals[vals_order]

    x = np.arange(len(vals)) / len(vals) * 100

    hl_mask = np.isin(vals['gene_set'], highlight)
    x_hl = x[hl_mask]
    vals_hl = vals[hl_mask]

    y_invert = False
    if yscale == 'mlog10':
        yscale = 'log'
        y_invert = True

    if xbaseline is None:
        if yscale == 'log':
            xbaseline = np.max(x)
        else:
            xbaseline = 0

    ax.set_yscale(yscale)
    if y_invert:
        ax.invert_yaxis()

    ax.fill_between(x, xbaseline, vals, facecolors='#777777',
                    step='mid', edgecolors='none')
    ax.vlines(x_hl, xbaseline, vals_hl, colors='#44ff44', linewidth=.8)

    if ylim is not None:
        ax.set_ylim(ylim)

    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(ylabel)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(True)
    ax.spines['right'].set_visible(False)
    ax.tick_params(bottom='off', top='off', left='on', right='off')
    ax.set_xticklabels("")


class SFDRNormalize(matplotlib.colors.Normalize):
    def __init__(self, vmin=None, vmax=None, clip=False,
                 sig_threshold=0.05):
        self.sig_threshold = sig_threshold
        matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        data_points = [self.vmin, -self.sig_threshold, 0, self.sig_threshold,
                       self.vmax]
        norm_points = [0, 0.499, 0.5, 0.501, 1]
        return np.ma.masked_array(np.interp(value, data_points, norm_points))
#####


def plot_gsea_heatmap(gsea, genesets_annot, fig):
    genesets = genesets_annot['gene_set'].values
    assert all(np.isin(genesets, gsea['gene_set']))
    wf_prop = 0.3
    table_prop = 0.3
    hm_prop = 1-wf_prop-table_prop
    cbar_vmargin = 0.25
    wf_hmargin = 0.05
    wf_vmargin = 0.01

    sel_gsea = gsea.reindex_like(genesets_annot)

    # Heatmap
    sel_gsea['slogfdr'] = np.sign(sel_gsea['nes']) * -np.log10(sel_gsea['fdr'])
    ax = fig.add_axes([wf_hmargin, table_prop,
                       1-wf_hmargin, hm_prop])
    hm = plot.heatmap(
        sel_gsea['slogfdr'],
        zlim=[-2, 2],
        norm=SFDRNormalize(sig_threshold=-np.log10(0.25)),
        method='pcolormesh',
        cbar=False,
        ax=ax,
    )
    ax_cbar = fig.add_axes(
        [wf_hmargin+cbar_vmargin, 0.02,
         1-wf_hmargin-2*cbar_vmargin, 0.03],
    )
    cbar = fig.colorbar(hm, ax_cbar, orientation='horizontal')
    fdr_ticks_at = np.array([0.25, 0.05, 0.01])
    lfdr_ticks_at = -np.log10(fdr_ticks_at)
    cbar_tick_lv = np.append(np.append(-lfdr_ticks_at[::-1], [0.0]),
                             lfdr_ticks_at)
    cbar_tick_v = np.append(-np.append(fdr_ticks_at[::-1], [0.0]),
                            fdr_ticks_at)
    cbar.set_ticks(cbar_tick_lv)
    ax_cbar.set_xlabel("signed FDR")
    cbar.ax.set_xticklabels([f"{v}" for v in cbar_tick_v])
    ax.set_xticklabels("")
    ax.tick_params(bottom='off')
    ax.set_xlabel("Gene Set")
    ax.set_ylabel("CAD Factor")
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Top waterfall plots
    ax_nes = fig.add_axes([(0/3)+wf_hmargin, 1-wf_prop+wf_vmargin,
                           (1/3)-wf_hmargin, wf_prop-wf_vmargin])
    wf_plot(gsea['nes'][0, :], genesets, ax_nes, 'NES')

    ax_mesa = fig.add_axes([(1/3)+wf_hmargin, 1-wf_prop+wf_vmargin,
                           (1/3)-wf_hmargin, wf_prop-wf_vmargin])
    mesa_mid = int(gsea['max_es_at'].max() / 2)
    wf_plot(gsea['max_es_at'][0, :], genesets, ax_mesa, 'Max. ES at',
            xbaseline=mesa_mid, reverse=True)
    ax_mesa.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

    ax_le = fig.add_axes([(2/3)+wf_hmargin, 1-wf_prop+wf_vmargin,
                          (1/3)-wf_hmargin, wf_prop-wf_vmargin])
    wf_plot(gsea['le_prop'][0, :], genesets, ax_le, 'Leading Edge')

    # Bottom table
    ga = genesets_annot.copy()
    sy = zip(ga['source'].values, ga['source_year'].values)
    ga['source'] = ('gene_set', [f"{s} ({y})" for s, y in sy])
    del ga['source_year']
    gaa = ga.to_array()
    table = ax.table(
        cellText=gaa.values,
        cellLoc='center',
        rowLabels=gaa['variable'].values,
        loc='bottom',
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    for (col, row), cell in table.get_celld().items():
        cell.set_linewidth(0.1)
        cell.set_edgecolor('w')
        if col % 2 == 1:
            cell.set_facecolor('#eeeeee')
        if col == 2:
            cell.set_height(2*cell.get_height())
        if row == -1:
            cell.set_text_props(weight='bold')


@click.command()
@click.argument('gsea_results', type=click_utils.in_path)
@click.argument('sel_genesets', type=click_utils.in_path)
@click.argument('out', type=click_utils.out_path)
def plot_gsea_heatmap_(gsea_results, sel_genesets, out):
    gsea = xr.open_dataset(gsea_results).load()
    gsea['gene_set'] = (xr.apply_ufunc(np.char.decode, gsea['gene_set'])
                        .astype('object'))
    gsea['mri_feature'] = gsea['mri_feature'].astype(int)

    geneset_annot = (pd.read_table(sel_genesets, sep='\t', quotechar='"',
                                   comment='#').
                     set_index('gene_set').to_xarray())

    with plot.figure(figsize=(7.0, 3.5)) as fig:
        plot_gsea_heatmap(gsea, geneset_annot, fig)
        fig.savefig(out, format="svg")


if __name__ == '__main__':
    plot_gsea_heatmap_()