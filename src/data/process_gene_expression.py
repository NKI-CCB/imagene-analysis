import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


def parse_args():
    parser = argparse.ArgumentParser(
        description="Clean up and annotate gene expression data set."
    )
    parser.add_argument(
        'gene_expression_data',
        help="NetCDF file gene expression data",
    )
    parser.add_argument(
        'sample_tracking',
        help="Sample tracking sheet to map identifiers.",
    )
    parser.add_argument(
        'gene_annotation',
        help="Annotation of genes",
    )
    parser.add_argument('out', help='Output NetCDF file.')
    args = parser.parse_args()

    # Check for existence for paths
    args.gene_expression_data = Path(args.gene_expression_data).resolve()
    args.sample_tracking = Path(args.sample_tracking).resolve()
    args.gene_annotation = Path(args.gene_annotation).resolve()
    args.out = Path(args.out)
    args.out = args.out.parent.resolve() / args.out.name

    return args


def counts_to_log2_cpm(counts):
    cpm = xr.DataArray(
        data=np.full(counts.shape, np.nan, dtype=np.float64),
        coords=counts.coords,
        dims=counts.dims,
    )
    cpm.attrs['unit'] = "lb(re 1)"
    cpm.attrs['long_name'] = "log(counts per million)"
    library_size = counts.sum('gene')
    cpm.values = np.log2((0.5 + counts) / (library_size+1)*1e6)

    return cpm


def map_sample_to_case(samples, sample_tracking_path):
    with sample_tracking_path.open() as f:
        reader = csv.DictReader(f, delimiter='\t')
        sample_case_map = {l['rna_sample']: l['margins_patient']
                           for l in reader}

    cases = [sample_case_map[s] for s in samples.values]
    cases = xr.DataArray(
        data=np.array(cases, dtype='int64'),
        coords=samples.coords,
        dims=samples.dims
    )
    cases.attrs['long_name'] = 'MARGINS Study Number'
    return cases


def merge_vals(v):
    v = set(v)
    if len(v) > 1:
        return ''
    else:
        return v.pop()


def annotate_genes(data_set, annot_path):
    data_set['gene'].values = [s.split('.')[0]
                               for s in data_set['gene'].values]
    annot = pd.read_table(annot_path, dtype=str, na_filter=False)
    annot = annot.groupby('Ensembl Gene ID').aggregate(merge_vals)
    annot = annot.loc[data_set['gene']]

    data_set['entrez_gene_id'] = xr.DataArray(
        data=np.array(annot['EntrezGene ID'].fillna(-1).replace('', -1),
                      dtype='int64'),
        dims=('gene',),
    )
    data_set['entrez_gene_id'].encoding['_FillValue'] = -1
    data_set['entrez_gene_id'].attrs['long_name'] = 'EntrezGene ID'

    data_set['hgnc_symbol'] = xr.DataArray(
        data=np.array(annot['HGNC symbol'].fillna(''), dtype='object'),
        dims=('gene',),
    )
    data_set['hgnc_symbol'].encoding['_FillValue'] = ''
    data_set['hgnc_symbol'].attrs['long_name'] = 'HGNC Symbol'

    return data_set


if __name__ == "__main__":
    args = parse_args()

    data_set = xr.open_dataset(str(args.gene_expression_data))

    data_set['log2_cpm'] = counts_to_log2_cpm(data_set['read_count'])
    data_set['case'] = map_sample_to_case(data_set['sample'],
                                          args.sample_tracking)
    data_set = data_set.swap_dims({'sample': 'case'})
    data_set = data_set.reset_coords(['sample'])

    data_set = annotate_genes(data_set, args.gene_annotation)

    time_str = (datetime.utcnow()
                .replace(microsecond=0, tzinfo=timezone.utc)
                .isoformat())
    data_set.attrs['history'] = (
        "{date} process_gene_expression.py Provide extra sample and gene "
        "annotation\n"
        .format(date=time_str) +
        data_set.attrs['history']
    )
    data_set.attrs['date_metadata_modified'] = time_str

    data_set.to_netcdf(str(args.out))
