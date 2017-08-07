import click
import xarray as xr


@click.command()
@click.option('--er-positive', is_flag=True)
@click.argument('filename', type=click.Path(exists=True))
@click.argument('clinical', type=click.Path(exists=True))
@click.argument('out_filename', type=click.Path())
def select_samples(filename, out_filename, clinical, er_positive):
    """Select a subset of samples."""
    mri = xr.open_dataset(filename).load()
    clin = xr.open_dataset(clinical).load()

    sel_cases = set(mri['case'].values)
    if er_positive:
        sel_cases &=\
            set(clin['case'].loc[clin['ihc_subtype'] == 'ER+/HER2-'].values)

    mri = mri.sel(case=list(sel_cases))

    mri.to_netcdf(out_filename)


if __name__ == '__main__':
    select_samples()
