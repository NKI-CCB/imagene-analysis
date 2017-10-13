from math import log

import click
import click_log
import numpy as np
import netCDF4
import xarray as xr

import funcsfa


def calc_bic(deviance, degrees_of_freedom, n_samples):
    return deviance + log(n_samples) * degrees_of_freedom


def dof_elastic_net(X, l2, B, eps=1e-6):
    df = 0.0
    for i in range(B.shape[0]):
        active_set = (B[i, :] > eps) | (B[i, :] < -eps)
        if sum(active_set) == 0:
            continue
        X_a = X[active_set, :].T
        I = np.identity(X_a.shape[1])
        df += np.trace(X_a @ np.linalg.inv(X_a.T @ X_a + l2*I) @ X_a.T)
    return(df)


def sparsity(coeff, eps=1e-6):
    return np.mean(np.abs(coeff) < eps)


click_in_path = click.Path(exists=True, dir_okay=False, resolve_path=True)
click_out_path = click.Path(exists=False, dir_okay=False, resolve_path=True)


@click.command()
@click.argument('parameter_sweep', type=click_in_path)
@click.argument('data', type=click_in_path)
@click.argument('out', type=click_out_path)
@click_log.simple_verbosity_option()
@click_log.init(__name__)
def eval_sfa_bic(parameter_sweep, data, out):
    data = funcsfa.StackedDataMatrix.from_netcdf(data)
    gexp = data.dt('gexp')
    mri = data.dt('mri_cad')
    n_samples = gexp.shape[0]
    assert mri.shape[0] == n_samples

    with netCDF4.Dataset(parameter_sweep, 'r') as ds:
        models = list(ds['models'].groups)

        nan_a = np.full(len(models), np.nan, np.float64)
        out_ds = xr.Dataset({
             'deviance_gexp': xr.DataArray(nan_a.copy(), dims=['model']),
             'dof_gexp': xr.DataArray(nan_a.copy(), dims=['model']),
             'deviance_mri': xr.DataArray(nan_a.copy(), dims=['model']),
             'dof_mri': xr.DataArray(nan_a.copy(), dims=['model']),
             'bic': xr.DataArray(nan_a.copy(), dims=['model']),
             'alpha': xr.DataArray(nan_a.copy(), dims=['model']),
             'l_gexp': xr.DataArray(nan_a.copy(), dims=['model']),
             'l_mri': xr.DataArray(nan_a.copy(), dims=['model']),
             'k': xr.DataArray(np.full(len(models), -1, np.int8),
                               dims=['model']),
             'sparsity_gexp': xr.DataArray(nan_a.copy(), dims=['model']),
             'sparsity_mri': xr.DataArray(nan_a.copy(), dims=['model']),
             'max_diff_coefficients': xr.DataArray(nan_a.copy(),
                                                   dims=['model']),
             'max_diff_factors': xr.DataArray(nan_a.copy(), dims=['model']),
             'n_iter': xr.DataArray(np.full(len(models), -1, np.int32),
                                    dims=['model']),
            },
            coords={
            'model': np.array(models, object),
            },
        )

        with click.progressbar(models) as model_bar:
            for model_id in model_bar:
                model_g = ds['models'][model_id]
                if 'factor_value' not in list(model_g.variables):
                    continue
                B_gexp = np.array(model_g['coefficient_gexp'])
                B_mri = np.array(model_g['coefficient_mri'])
                Z = np.array(model_g['factor_value'])
                dev_gexp = np.sum(np.square((gexp - (Z @ B_gexp)))) / n_samples
                out_ds['deviance_gexp'].loc[model_id] = dev_gexp
                dev_mri = np.sum(np.square((mri - (Z @ B_mri)))) / n_samples
                out_ds['deviance_mri'].loc[model_id] = dev_mri

                k = model_g.getncattr('k')
                out_ds['k'].loc[model_id] = k
                alpha = model_g.getncattr('alpha')
                out_ds['alpha'].loc[model_id] = alpha
                l_gexp = model_g.getncattr('l_gexp')
                out_ds['l_gexp'].loc[model_id] = l_gexp
                l2_gexp = 1e-6 + l_gexp * (1 - alpha)
                dof_gexp = dof_elastic_net(Z.T, l2_gexp / Z.shape[0], B_gexp.T)
                out_ds['dof_gexp'].loc[model_id] = dof_gexp
                l_mri = model_g.getncattr('l_mri')
                out_ds['l_mri'].loc[model_id] = l_mri
                l2_mri = 1e-6 + l_mri * (1 - alpha)
                dof_mri = dof_elastic_net(Z.T, l2_mri / Z.shape[0], B_mri.T)
                out_ds['dof_mri'].loc[model_id] = dof_mri

                out_ds['bic'].loc[model_id] = calc_bic(
                    dev_gexp + dev_mri,
                    dof_gexp + dof_mri + (Z.shape[0] * Z.shape[1]),
                    Z.shape[0])

                out_ds['sparsity_gexp'].loc[model_id] = sparsity(B_gexp)
                out_ds['sparsity_mri'].loc[model_id] = sparsity(B_mri)

                out_ds['max_diff_coefficients'].loc[model_id] =\
                    model_g['monitor']['max_diff_coefficients'][-1]
                out_ds['max_diff_factors'].loc[model_id] =\
                    model_g['monitor']['max_diff_factors'][-1]
                out_ds['n_iter'].loc[model_id] =\
                    model_g['monitor']['iteration'][-1]

        out_ds.attrs['empty_model_bic'] = calc_bic(
            np.sum(np.square(gexp)) + np.sum(np.square(mri)),
            0,
            Z.shape[0],
        )

    out_ds.to_netcdf(out)


if __name__ == '__main__':
    eval_sfa_bic()
