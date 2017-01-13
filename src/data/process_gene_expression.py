import argparse
import csv
from pathlib import Path

import numpy as np
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
    parser.add_argument('out', help='Output NetCDF file.')
    args = parser.parse_args()

    # Check for existence for paths
    args.gene_expression_data = Path(args.gene_expression_data).resolve()
    args.sample_tracking = Path(args.sample_tracking).resolve()
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


def map_sample_to_patient(samples, sample_tracking_path):
    with sample_tracking_path.open() as f:
        reader = csv.DictReader(f, delimiter='\t')
        sample_patient_map = {l['rna_sample']: l['margins_patient']
                              for l in reader}

    patients = [sample_patient_map[s] for s in samples.values]
    patients = xr.DataArray(
        data=np.array(patients, dtype='int64'),
        coords=samples.coords,
        dims=samples.dims
    )
    patients.attrs['long_name'] = 'MARGINS Study Number'
    return patients


if __name__ == "__main__":
    args = parse_args()

    data_set = xr.open_dataset(str(args.gene_expression_data))

    data_set = data_set.rename({'gene_expression': 'read_count'})
    data_set['log2_cpm'] = counts_to_log2_cpm(data_set['read_count'])
    data_set['patient'] = map_sample_to_patient(data_set['sample'],
                                                args.sample_tracking)
    data_set = data_set.swap_dims({'sample': 'patient'})
    data_set = data_set.reset_coords(['batch', 'sample'])

    data_set.to_netcdf(str(args.out))
