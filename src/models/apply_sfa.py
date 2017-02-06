import argparse
from pathlib import Path

import h5py
import numpy as np
import xarray as xr


def parse_args():
    parser = argparse.ArgumentParser(
        description="Apply sparse-factor analysis (SFA) model"
    )

    parser.add_argument('gene_expression', help="Gene expression data in"
                                                "netCDF format")
    parser.add_argument('sfa_model', help="SFA model in hdf5 file.")
    parser.add_argument('out', help="Applied SFA factors NetCDF output.")

    args = parser.parse_args()
    args.gene_expression = Path(args.gene_expression).resolve()
    args.sfa_model = Path(args.sfa_model).resolve()
    args.out = Path(args.out)
    args.out = args.out.parent.resolve() / args.out.name

    return args


def read_coefficients(sfa_model_path):
    from itertools import accumulate, chain

    with h5py.File(str(sfa_model_path)) as f:
        dt_names = [s.decode() for s in f.attrs['datatype names']]
        dt_len = list(f.attrs['number of features per data type'])
        dt_slice = [slice(s, s+l)
                    for s, l in zip(accumulate(chain([0], dt_len)), dt_len)]
        rnaseq_slice = dt_slice[dt_names.index('rnaseq')]

        entrezgene_id = [int(s.decode().split('|')[1])
                         for s in f['feature names'][rnaseq_slice]]

        coeff_array = np.array(f['coefficients'])

    factor = ['Factor {}'.format(i+1)
              for i in range(coeff_array.shape[1])]
    coefficients = xr.DataArray(
        data=coeff_array,
        dims=['entrez_gene_id', 'factor'],
        coords={
            'entrez_gene_id': entrezgene_id,
            'factor': np.array(factor, dtype='object')
        })

    return coefficients


def sel_most_var_gene(ds):
    return ds['gene'][np.argmax(ds['log2_cpm'].var('patient').values)]


def read_gene_expression(gene_expression_path):
    gexp_ds = xr.open_dataset(str(gene_expression_path), decode_cf=False)
    gexp_ds = gexp_ds.load()
    gexp_ds = gexp_ds.sel(gene=gexp_ds['entrez_gene_id'] !=
                          gexp_ds['entrez_gene_id'].attrs['_FillValue'])
    # There are Entrez gene ids with multiple Ensembl genes. Most of these
    # enseml genes are not expressed.
    gexp_ds = gexp_ds.sel(
        gene=gexp_ds['read_count'].mean('patient').values > 1.0
    )
    # If there are multiple expressed, select the most variable
    sel_genes = gexp_ds.groupby('entrez_gene_id').apply(sel_most_var_gene)
    gexp_ds = gexp_ds.sel(gene=sel_genes)

    gexp_ds = gexp_ds.swap_dims({'gene': 'entrez_gene_id'})

    return gexp_ds


def apply_sfa(gene_expression, coefficients):
    shared_genes = (gene_expression.indexes['entrez_gene_id'] &
                    coefficients.indexes['entrez_gene_id'])
    gene_expression = gene_expression.reindex(entrez_gene_id=shared_genes)
    gexp = gene_expression['log2_cpm']
    coefficients = coefficients.reindex(entrez_gene_id=shared_genes)

    X = np.array(gexp.transpose('entrez_gene_id', 'patient').values)
    X = X - np.mean(X, 1, keepdims=True)
    B = np.array(coefficients.transpose('entrez_gene_id', 'factor').values)

    btb = (B.T @ B + np.identity(B.shape[1]))
    o = B @ np.linalg.inv(btb)
    Z = o.T @ X

    return xr.DataArray(
        data=Z,
        dims=('factor', 'patient'),
        coords={'factor': coefficients.coords['factor'],
                'patient': gene_expression.coords['patient']},
    )


if __name__ == "__main__":

    args = parse_args()

    coefficients = read_coefficients(args.sfa_model)
    gene_expression = read_gene_expression(args.gene_expression)
    factors = apply_sfa(gene_expression, coefficients)

    ds = xr.Dataset({'factors': factors})
    ds.to_netcdf(str(args.out))
