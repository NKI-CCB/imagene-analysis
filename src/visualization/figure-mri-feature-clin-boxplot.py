import click
import numpy as np
import scipy.stats
import yaml
import xarray as xr

from lib import click_utils
import plot
from visualization.style import set_style
from visualization.labels import feature_display_names, clin_display_names


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
@click.argument('mri_features', type=click_utils.in_path)
@click.argument('feature_id', type=str)
@click.argument('clinical_annotation', type=click_utils.in_path)
@click.argument('clinical_var', type=str)
@click.argument('out', type=click_utils.out_path)
@click.argument('stats_out', type=click_utils.out_path)
def plot_factor_in_subtype(mri_features, feature_id, clinical_annotation,
                           clinical_var, out, stats_out):
    mri_ds = xr.open_dataset(mri_features)
    feature = mri_ds[feature_id]

    clin = xr.open_dataset(clinical_annotation)[clinical_var].load()
    feature, clin = xr.align(feature[feature.notnull()], clin[clin.notnull()])

    with plot.subplots(figsize=(3.5, 3.5)) as (fig, ax):
        plot.boxplot(
            clin, feature,
            title="",
            xlabel=clin_display_names[clinical_var],
            ylabel=feature_display_names[feature_id],
            ax=ax,
        )
        fig.savefig(out, format='svg')

    feature_by_clin = split_by(feature.values, clin.values)
    h, p = scipy.stats.kruskal(*feature_by_clin.values())
    with open(stats_out, 'w') as f:
        f.write(f"h: {h:.6e}\n")
        f.write(f"p: {p:.6e}\n")


if __name__ == '__main__':
    set_style()
    plot_factor_in_subtype()
