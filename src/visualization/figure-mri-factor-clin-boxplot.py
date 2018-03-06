import click
import numpy as np
import scipy.stats
import yaml
import xarray as xr

from lib import click_utils
import plot
from visualization.style import set_style
from visualization.labels import factor_display_names, clin_display_names


def only(lst):
    assert len(lst) == 1
    return lst[0]


def split_by(x, y):
    x_split = dict()
    for cat in set(y):
        if cat is np.nan:
            continue
        x_split[cat] = x[y == cat]
    return x_split


@click.command()
@click.argument('cad_factors', type=click_utils.in_path)
@click.argument('factor_id', type=str)
@click.argument('factor_annotation', type=click_utils.in_path)
@click.argument('clinical_annotation', type=click_utils.in_path)
@click.argument('clinical_var', type=str)
@click.argument('out', type=click_utils.out_path)
@click.argument('stats_out', type=click_utils.out_path)
def plot_factor_in_subtype(cad_factors, factor_id, factor_annotation,
                           clinical_annotation, clinical_var, out, stats_out):
    factor_da = xr.open_dataset(cad_factors)['factors']
    with open(factor_annotation) as f:
        factor_index = only([i for i, v in yaml.load(f).items()
                             if v['id'] == factor_id])
    factor = factor_da.sel(factor=factor_index).load()
    factor.name = factor_id

    clin = xr.open_dataset(clinical_annotation)[clinical_var].load()
    factor, clin = xr.align(factor[factor.notnull()], clin[clin.notnull()])

    with plot.subplots(figsize=(3.5, 3.5)) as (fig, ax):
        plot.boxplot(
            clin, factor,
            title="",
            xlabel=clin_display_names[clinical_var],
            ylabel=factor_display_names[factor_id],
            ax=ax,
        )
        fig.savefig(out, format='svg')

    factor_by_clin = split_by(factor.values, clin.values)
    h, p = scipy.stats.kruskal(*factor_by_clin.values())
    with open(stats_out, 'w') as f:
        f.write(f"h: {h:.6e}\n")
        f.write(f"p: {p:.6e}\n")


if __name__ == '__main__':
    set_style()
    plot_factor_in_subtype()
