from datetime import datetime, timezone

import click
import factor_rotation
import numpy as np
import sklearn.decomposition
import xarray as xr


def read_mri(data_set):
    data_set = data_set.copy()
    del data_set['Comment']
    del data_set['MultiFocal']
    mri = data_set.to_array('cad_feature', 'mri_cad_features')
    data_set.close()

    mri = mri.isel(case=np.where(mri.isnull().sum('cad_feature') == 0)[0])
    mri = mri.transpose('case', 'cad_feature')

    return mri


def adjust_scale(mri):
    mri_adj = xr.DataArray(np.full(mri.shape, np.nan), mri.coords, mri.dims)
    for feature in mri['cad_feature'].values:
        if feature[0:3] == 'vol':
            mri_adj.loc[:, feature] = np.cbrt(mri.loc[:, feature])
        elif feature[0:3] == 'var':
            mri_adj.loc[:, feature] = np.sqrt(mri.loc[:, feature])
        else:
            mri_adj.loc[:, feature] = mri.loc[:, feature]
    return mri_adj


def compute_factors(mri, n_components):
    pca = sklearn.decomposition.PCA()
    mri_a = (mri / mri.std('case')).values
    pca.fit(mri_a)
    rotated_components, rotation = factor_rotation.rotate_factors(
        pca.components_[:n_components, :].T,
        'varimax',
    )
    rotated_components = rotated_components.T
    factors_a = pca.transform(mri_a)[:, :n_components] @ rotation

    # Check results
    X = (mri / mri.std('case'))
    X = X - X.mean('case')
    X = X.values
    X_rec = factors_a @ rotated_components
    err = np.mean((X - X_rec)**2) / np.mean(X**2)
    assert err < 0.05, "High reconstruction error {}".format(err)

    # Add metadata back
    factors = xr.DataArray(
        factors_a,
        dims=['case', 'factor'],
        coords={
            'case': mri.coords['case'],
            'factor': np.array(range(1, n_components+1), 'int8'),
        },
    )
    loadings = xr.DataArray(
        rotated_components,
        dims=['factor', 'cad_feature'],
        coords={
            'factor': np.array(range(1, n_components+1), 'int8'),
            'cad_feature': mri.coords['cad_feature'],
        }
    )

    return factors, loadings


@click.command()
@click.argument('n_components', type=int)
@click.argument('filename', type=click.Path(exists=True))
@click.argument('out_filename', type=click.Path())
def fa_mri_features(filename, out_filename, n_components):
    """Regress volume out of MRI features."""
    mri_data_set = xr.open_dataset(filename).load()

    mri = read_mri(mri_data_set)
    mri = adjust_scale(mri)

    factors, loadings = compute_factors(mri, n_components)

    fa_data_set = xr.Dataset({'factors': factors, 'loadings': loadings})

    time_str = (datetime.utcnow()
                .replace(microsecond=0, tzinfo=timezone.utc)
                .isoformat())
    fa_data_set.attrs['history'] = (
        "{time} fa_mri_features.py Factor analysis\n"
        .format(time=time_str) +
        mri_data_set.attrs['history']
    )

    fa_data_set.to_netcdf(out_filename)


if __name__ == '__main__':
    fa_mri_features()
