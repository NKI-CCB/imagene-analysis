import click
import xarray as xr

from lib import click_utils
import plot
from visualization.style import set_style


@click.command()
@click.argument('cad_factors', type=click_utils.in_path)
@click.argument('tumor_size_factor', type=int)
@click.argument('clinical_annotation', type=click_utils.in_path)
@click.argument('out', type=click_utils.out_path)
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
    set_style()
    plot_tumor_size_in_subtype()
