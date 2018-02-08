import click
import xarray as xr

import plot

from visualization.labels import feature_order, feature_display_names
from visualization.style import set_style


click_in_path = click.Path(exists=True, dir_okay=False, resolve_path=True)
click_out_path = click.Path(exists=False, dir_okay=False, resolve_path=True)


@click.command()
@click.argument('cad_factors', type=click_in_path)
@click.argument('out', type=click_out_path)
def plot_mri_cad_factors(cad_factors, out):
    fa_dataset = xr.open_dataset(cad_factors).load()

    assert all(f in feature_order for f in fa_dataset['cad_feature'].values)
    fa_dataset = fa_dataset.reindex(cad_feature=feature_order)

    with plot.subplots(figsize=(3.5, 3.5)) as (fig, ax):
        plot.heatmap(
            fa_dataset['loadings'].T,
            aspect='equal',
            xlabel="CAD Factor",
            ylabel="CAD Feature",
            yticklabels=[feature_display_names[f]
                         for f in fa_dataset.coords['cad_feature'].values],
            zlabel="Loading",
            ax=ax,
        )

        fig.savefig(out, format="svg")


if __name__ == '__main__':
    set_style()
    plot_mri_cad_factors()
