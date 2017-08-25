from collections import namedtuple
from datetime import datetime, timezone
import logging

import click
import click_log
import numpy as np
import xarray as xr


logger = logging.getLogger(__name__)


VoomResult = namedtuple('VoomResult', ['expression', 'weights'])


def voom(counts, library_size):
    from rpy2.robjects.packages import importr
    from rpy2.robjects.numpy2ri import numpy2ri

    logger.info("Running limma voom in R")

    limma = importr('limma')
    edgeR = importr('edgeR')
    base_r = importr('base')
    r_dollar = getattr(base_r, '$')

    library_size_r = base_r.c(numpy2ri(library_size.values))
    counts_r = edgeR.DGEList(counts=numpy2ri(counts.values.T),
                             lib_size=library_size_r)

    counts_r = edgeR.calcNormFactors(counts_r)
    v = limma.voom(counts_r, plot=False)

    gexp = xr.DataArray(
        np.array(r_dollar(v, 'E')).T,
        coords=counts.coords,
        attrs={
            'units': 'lb(re 1)',
            'long_name': "Gene expression in log2 range"
        })
    weights = xr.DataArray(
        np.array(r_dollar(v, 'weights')).T,
        coords=counts.coords,
        attrs={
            'units': 'lb(re 1)',
            'long_name': "Limma voom weights"
        })

    return VoomResult(gexp, weights)


click_in_path = click.Path(exists=True, dir_okay=False, resolve_path=True)
click_out_path = click.Path(exists=False, dir_okay=False, resolve_path=True)


@click.command()
@click.argument('gexp', type=click_in_path)
@click.argument('out', type=click_out_path)
@click_log.simple_verbosity_option()
@click_log.init(__name__)
def run_sfa(gexp, out):
    logging.info("Running limma voom")
    ds = xr.open_dataset(gexp).load()
    if 'log2_cpm' in ds:
        del ds['log2_cpm']
    library_size = (ds['read_count'].sum('gene') + ds['N_unmapped'] +
                    ds['N_multimapping'] + ds['N_noFeature'] +
                    ds['N_ambiguous'])
    log2_cpm, weights = voom(ds['read_count'], library_size)

    logger.info("Preparing output")
    ds['log2_cpm'] = log2_cpm
    ds['weight'] = weights
    del ds['read_count']
    del ds['N_unmapped']
    del ds['N_multimapping']
    del ds['N_noFeature']
    del ds['N_ambiguous']

    time_str = (datetime.utcnow()
                .replace(microsecond=0, tzinfo=timezone.utc)
                .isoformat())
    ds.attrs['history'] = (
        "{date} process_gene_expression_voom.py Apply Limma-Voom\n"
        .format(date=time_str) +
        ds.attrs['history']
    )

    logger.info("Writing result to {}".format(out))
    ds.to_netcdf(out)


if __name__ == '__main__':
    run_sfa()
