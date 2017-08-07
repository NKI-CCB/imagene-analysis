import click
import numpy as np
import pandas as pd
import xarray as xr


def convert_to_boolean(da, values=('F', 'T')):
    a = da.values
    bin_a = np.full(a.shape, -1, dtype='i1')
    bin_a[a == values[0]] = 0
    bin_a[a == values[1]] = 1
    bin_da = xr.DataArray(bin_a, da.coords, da.dims, da.name, da.attrs)
    bin_da.encoding['_FillValue'] = -1

    return bin_da


@click.command()
@click.argument('filename', type=click.Path(exists=True))
@click.argument('out_filename', type=click.Path())
def process_clinical(filename, out_filename):
    """Convert clinical tsv to netcdf."""
    df = pd.read_table(filename)
    df = df.rename(columns={
        'margins_patient': 'case',
        'AdjRT': 'adjuvant_radiotherapy',
        'AdjChemo': 'adjuvant_chemotherapy',
        'AdjHormo': 'adjuvant_hormonal_therapy',
        'AdjAntiHER2': 'adjuvant_anti_her2_therapy',
        'AdjSystemic': 'adjuvant_systemic_therapy',
        'pos_LN': 'positive_lymph_nodes',
        'largest_diameter_MRI': 'largest_diameter_mri',
        'histograde': 'grade',
        'age_at_diag': 'age_at_diagnosis',
    })
    del df['rna_sample']  # Shouldn't be in this data set.
    df = df.set_index('case')

    ds = df.to_xarray()
    ds['age_at_diagnosis'].attrs['units'] = 'year'
    ds['largest_diameter_mri'].attrs['units'] = 'cm'

    ds['adjuvant_radiotherapy'] =\
        convert_to_boolean(ds['adjuvant_radiotherapy'])
    ds['adjuvant_chemotherapy'] =\
        convert_to_boolean(ds['adjuvant_chemotherapy'])
    ds['adjuvant_hormonal_therapy'] =\
        convert_to_boolean(ds['adjuvant_hormonal_therapy'])
    ds['adjuvant_anti_her2_therapy'] =\
        convert_to_boolean(ds['adjuvant_anti_her2_therapy'])
    ds['adjuvant_systemic_therapy'] =\
        convert_to_boolean(ds['adjuvant_systemic_therapy'])

    ds.to_netcdf(out_filename)


if __name__ == '__main__':
    process_clinical()
