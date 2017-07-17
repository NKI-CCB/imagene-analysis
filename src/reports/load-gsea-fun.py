def load_gsea_ds(fn):
    ds = xr.open_dataset(fn).load()
    ds['mri_feature'] = [s.decode() for s in ds['mri_feature'].values]
    ds['mri_feature'] = [" ".join(s.split('_')).title()
                         for s in ds['mri_feature'].values]
    ds['gene_set_code'] = ('gene_set',
                           [s.decode() for s in ds['gene_set'].values])
    ds['gene_set'] = [" ".join(s.split('_')).title()
                      for s in ds['gene_set_code'].values]
    return ds
