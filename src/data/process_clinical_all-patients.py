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
    columns = {
        'StudyNumber': 'case',
        'AdjRT': 'adjuvant_radiotherapy',
        'AdjChemo': 'adjuvant_chemotherapy',
        'AdjHormo': 'adjuvant_hormonal_therapy',
        'AdjHER2': 'adjuvant_anti_her2_therapy',
        'LymphNodePos_BV': 'positive_lymph_nodes',
        'Diameter_BV': 'largest_diameter_mri',
        'Histograde_BV': 'grade',
        'Age': 'age_at_diagnosis',
        'IHC_1erpos_2her2pos_3tripneg': 'ihc_subtype',
    }
    df = df.loc[:, columns.keys()]
    df = df.rename(columns=columns)
    df = df.set_index('case')

    ds = df.to_xarray()
    ds['age_at_diagnosis'].attrs['units'] = 'year'
    ds['largest_diameter_mri'].attrs['units'] = 'mm'

    ds['adjuvant_radiotherapy'] =\
        convert_to_boolean(ds['adjuvant_radiotherapy'], ('N', 'J'))
    ds['adjuvant_chemotherapy'] =\
        convert_to_boolean(ds['adjuvant_chemotherapy'], ('N', 'J'))
    ds['adjuvant_hormonal_therapy'] =\
        convert_to_boolean(ds['adjuvant_hormonal_therapy'], ('N', 'J'))
    ds['adjuvant_anti_her2_therapy'] =\
        convert_to_boolean(ds['adjuvant_anti_her2_therapy'], ('N', 'J'))

    assert np.all(ds['grade'] > 0)
    assert np.all((ds['grade'] < 4) |
                  (ds['grade'] == 777) |
                  (ds['grade'] == 999))
    ds['grade'][ds['grade'] == 777] = -1
    ds['grade'][ds['grade'] == 999] = -1
    ds['grade'] = ds['grade'].astype('i1')
    ds['grade'].encoding['_FillValue'] = -1

    ds['largest_diameter_mri'] = ds['largest_diameter_mri'].astype('f8')
    assert np.all(ds['largest_diameter_mri'] > 0)
    assert np.all((ds['largest_diameter_mri'] < 500) |
                  (ds['largest_diameter_mri'] == 999))
    ds['largest_diameter_mri'][ds['largest_diameter_mri'] == 999] = np.nan
    ds['largest_diameter_mri'].encoding['_FillValue'] = np.nan

    assert np.all(ds['positive_lymph_nodes'].values >= 0)
    assert np.all((ds['positive_lymph_nodes'] < 50) |
                  (ds['positive_lymph_nodes'] == 999))
    ds['positive_lymph_nodes'][ds['positive_lymph_nodes'] == 999] = -1
    ds['positive_lymph_nodes'] = ds['positive_lymph_nodes'].astype('i1')
    ds['positive_lymph_nodes'].encoding['_FillValue'] = -1

    assert np.all(ds['age_at_diagnosis'].values >= 0)
    assert np.all((ds['age_at_diagnosis'] < 200) |
                  (ds['age_at_diagnosis'] == 999))
    ds['age_at_diagnosis'][ds['age_at_diagnosis'] == 999] = -1
    ds['age_at_diagnosis'] = ds['age_at_diagnosis'].astype('i1')
    ds['age_at_diagnosis'].encoding['_FillValue'] = -1

    ihc_subtype = np.full(ds['ihc_subtype'].shape, '', object)
    assert np.all(np.isin(ds['ihc_subtype'].values, [0, 1, 2, 3, 555, 999]))
    ihc_subtype[ds['ihc_subtype'] == 1] = 'ER+'
    ihc_subtype[ds['ihc_subtype'] == 2] = 'HER2+/ER-'
    ihc_subtype[ds['ihc_subtype'] == 3] = 'TN'
    ihc_subtype[ds['ihc_subtype'] == 555] = ''
    ihc_subtype[ds['ihc_subtype'] == 999] = ''
    ds['ihc_subtype'] = ('case', ihc_subtype)

    ds.to_netcdf(out_filename)


if __name__ == '__main__':
    process_clinical()
