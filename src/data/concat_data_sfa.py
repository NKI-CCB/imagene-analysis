import click
import numpy as np
import xarray as xr

import funcsfa

click_in_path = click.Path(exists=True, dir_okay=False, resolve_path=True)
click_out_path = click.Path(exists=False, dir_okay=False, resolve_path=True)


@click.command()
@click.argument('gexp', type=click_in_path)
@click.argument('mri', type=click_in_path)
@click.argument('out', type=click_out_path)
def prepare_data(gexp, mri, out):
    gexp = xr.open_dataset(gexp).load()
    mri = xr.open_dataset(mri).load()

    mri_features = list(set(mri.data_vars) - {'Comment', 'MultiFocal'})
    mri_da = mri[mri_features].to_array('cad_feature', 'mri_cad')
    mri_da = mri_da.reindex(cad_feature=mri_features)

    mri_da = mri_da.transpose('case', 'cad_feature')

    pce_null_v = mri_da.loc[:, 'PCE_top10percent'].isnull()
    pce_mean = mri_da.loc[:, 'PCE_top10percent'].mean()
    mri_da.loc[pce_null_v, 'PCE_top10percent'] = pce_mean

    sel_cases = list(
        set(gexp['case'].values) &
        set(mri_da['case'][mri_da.isnull().sum('cad_feature') == 0].values)
    )
    mri_da = mri_da.reindex(case=sel_cases)
    gexp = gexp.reindex(case=sel_cases)

    gexp['gexpw'] = gexp['log2_cpm'] * gexp['weight']
    gexp['gexpw'] = gexp['gexpw'] - gexp['gexpw'].mean('case')
    gene_mad = np.abs(gexp['gexpw']).median('case')
    sel_genes = gexp['gene'][np.argsort(-gene_mad)[0:1000]]
    gexp = gexp.reindex(gene=sel_genes)

    gexp_std = np.std(gexp['gexpw']).values.item()
    gexp_dm = funcsfa.DataMatrix(
        gexp['log2_cpm'].values,
        samples=sel_cases,
        features=gexp['gene'].values,
        weights=gexp['weight'].values / gexp_std,
    )
    assert(np.allclose(np.var(gexp_dm.dataW), 1.0))
    mri_dm = funcsfa.DataMatrix(
        mri_da.values,
        samples=sel_cases,
        features=mri_features,
        weights=np.full(mri_da.shape, 1 / (mri_da - mri_da.mean('case')).std())
    )
    assert(np.allclose(np.var(mri_dm.dataW), 1.0))

    res = funcsfa.StackedDataMatrix([gexp_dm, mri_dm], ['gexp', 'mri_cad'])

    res.to_netcdf(out, ['gene', 'cad_feature'])


if __name__ == '__main__':
    prepare_data()
