import argparse
from pathlib import Path

import pandas as pd
import xarray as xr


def parse_args():
    parser = argparse.ArgumentParser(
        description="Clean up MRI dataset and convert to NetCDF"
    )
    parser.add_argument('mri_data', help="xlsx file with MRI features.")
    parser.add_argument('out', help='Output NetCDF file.')
    parser.add_argument('--study-nr-col',
                        default='MARGINSstudyNr',
                        help='Column name with MARGINS Study Number')
    args = parser.parse_args()

    args.mri_data = Path(args.mri_data).resolve()
    args.out = Path(args.out)
    args.out = args.out.parent.resolve() / args.out.name

    return args


def read_mri_xlsx(mri_path, study_nr_col):
    with mri_path.open('rb') as f:
        mri_df = pd.read_excel(f)

    if study_nr_col not in mri_df.columns:
        raise Exception('Could not find margins study number column')

    mri_df = mri_df.rename(columns={study_nr_col: 'patient'})
    mri_df = mri_df.set_index('patient')

    mri_ds = mri_df.to_xarray()
    mri_ds['patient'].attrs['long_name'] = "MARGINS Study Number"

    return mri_ds


if __name__ == "__main__":
    args = parse_args()

    data_set = read_mri_xlsx(args.mri_data, args.study_nr_col)
    data_set.to_netcdf(str(args.out))
