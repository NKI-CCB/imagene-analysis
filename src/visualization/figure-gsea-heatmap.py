import click

from lib import click_utils
import plot


@click.command()
@click.argument('gsea_results', type=click_utils.in_path)
@click.argument('out', type=click_utils.out_path)
def plot_gsea_heatmap(cad_factors, out):

    with plot.subplots(figsize=(3.5, 3.5)) as (fig, ax):
        fig.savefig(out, format="svg")


if __name__ == '__main__':
    plot_gsea_heatmap()
