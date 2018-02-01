import click
import matplotlib
import xarray as xr

import plot


click_in_path = click.Path(exists=True, dir_okay=False, resolve_path=True)
click_out_path = click.Path(exists=False, dir_okay=False, resolve_path=True)


@click.command()
@click.argument('cad_factors', type=click_in_path)
@click.argument('tumor_size_factor', type=int)
@click.argument('clinical_annotation', type=click_in_path)
@click.argument('out', type=click_out_path)
def plot_tumor_size_in_subtype(cad_factors, tumor_size_factor,
                               clinical_annotation, out):
    factor_da = xr.open_dataset(cad_factors)['factors']
    tumor_size = factor_da.sel(factor=tumor_size_factor).load()

    subtype = xr.open_dataset(clinical_annotation)['ihc_subtype']
    cases = tumor_size.coords['case'][tumor_size.notnull()]
    tumor_size = tumor_size.reindex(case=cases)
    subtype = subtype.reindex(case=cases)
    subtype = xr.DataArray([s.item() for s in subtype],
                           dims=['case'],
                           coords={'case': cases})

    print(subtype)
    print(tumor_size)

    with plot.subplots(figsize=(3.5, 3.5)) as (fig, ax):
        plot.boxplot(
            subtype, tumor_size,
            title="",
            xlabel="IHC Subtype",
            ylabel="MRI CAD Factor Size",
            ax=ax,
        )
        fig.savefig(out, format="svg")


if __name__ == '__main__':
    matplotlib.rcParams['font.family'] = 'sans-serif'
    matplotlib.rcParams['font.sans-serif'] = ['Arial']
    matplotlib.rcParams['font.weight'] = 'normal'
    matplotlib.rcParams['figure.dpi'] = 300
    matplotlib.rcParams['figure.facecolor'] = 'none'
    matplotlib.rcParams['axes.labelsize'] = 11
    matplotlib.rcParams['xtick.labelsize'] = 9
    matplotlib.rcParams['ytick.labelsize'] = 9

    plot_tumor_size_in_subtype()
