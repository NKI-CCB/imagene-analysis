import argparse
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


def parse_args():
    parser = argparse.ArgumentParser(
        description="Clean up MRI eigenbreasts dataset and convert to NetCDF"
    )
    parser.add_argument('mri_data', nargs="+",
                        help="xlsx files with MRI features.")
    parser.add_argument('out', help='Output NetCDF file.')
    parser.add_argument('--study-nr-col',
                        default='StudyID',
                        help='Column name with MARGINS Study Number')
    args = parser.parse_args()

    args.mri_data = [Path(s).resolve() for s in args.mri_data]
    args.out = Path(args.out)
    args.out = args.out.parent.resolve() / args.out.name

    return args


def read_mri_xlsx(mri_path, study_nr_col):
    with mri_path.open('rb') as f:
        mri_df = pd.read_excel(f)

    if study_nr_col not in mri_df.columns:
        raise Exception("Could not find margins study number column")

    mri_df = mri_df.rename(columns={study_nr_col: 'case'})
    mri_df = mri_df.set_index('case')

    mri_da = xr.DataArray(
        mri_df,
        dims=['case', 'PC'],
        coords={
            'case': mri_df.index,
            'PC': np.array([int(i) for i in mri_df.columns], 'int16'),
        }
    ).T

    return mri_da


if __name__ == "__main__":
    args = parse_args()

    data_set = dict()
    for pth in args.mri_data:
        var_name = pth.stem.split('_', 1)[1]
        data_set[var_name] = read_mri_xlsx(pth, args.study_nr_col)
    data_set = xr.Dataset(data_set)
    data_set['case'].attrs['long_name'] = "MARGINS Study Number"
    data_set.attrs['title'] = ("MRI eigenbreast features from Margins of "
                               "samples with gene expression data from "
                               "Imagene")

    time_str = (datetime.utcnow()
                .replace(microsecond=0, tzinfo=timezone.utc)
                .isoformat())
    data_set.attrs['history'] = (
        "{time} process_mri_eigenbreasts.py Converted from {fns}."
        .format(time=time_str, fns=", ".join([str(p) for p in args.mri_data]))
    )

    data_set.to_netcdf(args.out)
