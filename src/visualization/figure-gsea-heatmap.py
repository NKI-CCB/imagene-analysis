import click
import matplotlib
import xarray as xr

import plot


click_in_path = click.Path(exists=True, dir_okay=False, resolve_path=True)
click_out_path = click.Path(exists=False, dir_okay=False, resolve_path=True)


@click.command()
@click.argument('gsea_results', type=click_in_path)
@click.argument('out', type=click_out_path)
def plot_gsea_heatmap(cad_factors, out):

    with plot.subplots(figsize=(3.5, 3.5)) as (fig, ax):
        fig.savefig(out, format="svg")


if __name__ == '__main__':
    plot_gsea_heatmap()
