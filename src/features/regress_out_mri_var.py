from datetime import datetime, timezone

import click
import numpy as np
import xarray as xr


def regress_out(var, expl_var):
    X = np.vstack([expl_var, np.ones(len(expl_var))]).T
    X[np.isnan(expl_var), 0] = np.nanmean(expl_var)
    y = var
    y[np.isnan(y)] = np.nanmean(y)
    coeff, _, _, _ = np.linalg.lstsq(X, var)
    y_hat = X @ coeff
    residuals = y - y_hat
    residuals[np.isnan(var)] = np.nan
    residuals[np.isnan(expl_var)] = np.nan
    return residuals


@click.command()
@click.argument('filename', type=click.Path(exists=True))
@click.argument('out_filename', type=click.Path())
def regress_out_mri_var(filename, out_filename):
    """Regress volume out of MRI features."""
    data_set = xr.open_dataset(filename).load()

    expl_var = np.cbrt(data_set['volume'])
    for var in data_set.data_vars.values():
        # Skip multidimensional variables
        if len(var.dims) != 1:
            continue
        # Skip variables with other dimensions
        if var.dims[0] != expl_var.dims[0]:
            continue
        # Skip non numeric variables
        if not np.issubdtype(var.dtype, np.number):
            continue

        if var.name[0:3] == 'vol':
            vv = np.cbrt(var.values)
        elif var.name[0:3] == 'var':
            vv = np.sqrt(var.values)
        else:
            vv = var.values
        var.values = regress_out(vv, expl_var.values)

    time_str = (datetime.utcnow()
                .replace(microsecond=0, tzinfo=timezone.utc)
                .isoformat())
    data_set.attrs['history'] = (
        "{time} regress_out_mri_var.py Regressed volume out\n"
        .format(time=time_str) +
        data_set.attrs['history']
    )

    data_set.to_netcdf(out_filename)


if __name__ == '__main__':
    regress_out_mri_var()
