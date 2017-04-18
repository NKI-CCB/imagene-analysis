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

transformation_functions = {
    'linear': lambda x: x,
    'square root': np.sqrt,
    'cubic root': np.cbrt,
}


@click.command()
@click.argument('filename', type=click.Path(exists=True))
@click.argument('variable_name')
@click.argument('out_filename', type=click.Path())
@click.option('--trans-lin', 'transformation', flag_value='linear',
              default=True,
              help="Apply no transformation to explanatory variable. "
                   "[default]")
@click.option('--trans-sqrt', 'transformation', flag_value='square root',
              help="Apply square root transformation to explanatory variable.")
@click.option('--trans-cbrt', 'transformation', flag_value='cubic root',
              help="Apply cubic root transformation to explanatory variable.")
@click.pass_context
def regress_out_mri_var(ctx, filename, variable_name, out_filename,
                        transformation):
    """Regress a variable out of the MRI features to make other features
       independent of it.
    """
    data_set = xr.open_dataset(filename).load()

    if variable_name not in data_set:
        ctx.fail(
            '"{var}" not in data set "{fn}".'
            .format(var=variable_name, fn=filename)
        )

    trans_f = transformation_functions[transformation]
    expl_var = trans_f(data_set[variable_name])
    if len(expl_var.dims) != 1:
        ctx.fail(
            '"{var}" is not a one-dimensional variable.'
            .format(var=variable_name, err=True)
        )
        ctx.abort()

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
        var.values = regress_out(var.values, expl_var.values)

    time_str = (datetime.utcnow()
                .replace(microsecond=0, tzinfo=timezone.utc)
                .isoformat())
    data_set.attrs['history'] = (
        "{time} regress_out_mri_var.py "
        "Regressed {variable_name} with {transformation} transformation out.\n"
        .format(time=time_str, variable_name=variable_name,
                transformation=transformation) +
        data_set.attrs['history']
    )

    data_set.to_netcdf(out_filename)

if __name__ == '__main__':
    regress_out_mri_var()
