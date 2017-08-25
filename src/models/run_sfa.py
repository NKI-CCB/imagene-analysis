from collections import namedtuple
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import product
import logging

import click
import click_log
import numpy as np
import netCDF4

import funcsfa


logger = logging.getLogger(__name__)


_Params = namedtuple('Params', ['k', 'alpha', 'l_gexp', 'l_mri'])


class Params(_Params):

    def __str__(self):
        return (
            "Params(k={k:d}, alpha={alpha:.2e}, l_gexp={l_gexp:.2e}, "
            "l_mri={l_mri:.2e})"
            .format(k=self.k, alpha=self.alpha, l_gexp=self.l_gexp,
                    l_mri=self.l_mri)
        )


class SFAResult():

    def __init__(self, params):
        self.params = params
        self.monitor = None
        self.factors = None
        self.coefficients = None
        self.error = None

    def set_result(self, monitor, factors, coefficients):
        self.monitor = monitor
        self.factors = factors
        self.coefficients = coefficients

    def set_error(self, error):
        self.error = error


def sfa(data, params, max_iter, eps=1e-6, l2_eps=1e-6):
    l1_gexp = params.l_gexp * params.alpha
    l2_gexp = l2_eps + params.l_gexp * (1-params.alpha)
    l1_mri = params.l_mri * params.alpha
    l2_mri = l2_eps + params.l_mri * (1-params.alpha)

    sfa = funcsfa.SFA()
    result = SFAResult(params)
    try:
        logger.info("Starting with {}".format(params))
        mon = sfa.monitored_fit(data, params.k,
                                (l1_gexp, l1_mri), (l2_gexp, l2_mri),
                                max_iter=max_iter, eps=eps)
        logger.info("Finished with {}".format(params))

        result.set_result(mon, sfa.transform(data), sfa.coefficients)
    except Exception as e:
        logger.exception("Error with {}".format(params))
        result.set_error(e)

    return result


def parse_int_range(ctx, param, range_str):
    """Parse a string that can be interpreted as a range.

    For example:
        "1:10" = range(0, 10)
        "1:100:10" = range(0, 100, 10)
    """
    import argparse
    range_a = [int(n) for n in range_str.split(':')]
    if len(range_a) == 1:
        range_a = [range_a[0], range_a[0], 1]
    if len(range_a) == 2:
        range_a += [1]  # default step size 1
    if len(range_a) != 3:
        raise argparse.ArgumentTypeError('Range improperly formatted')

    return list(range(range_a[0], range_a[1]+1, range_a[2]))


def parse_zero_one_range(ctx, param, range_str):
    """Parse a string that can be interpreted as a range between zero and one.

    For example:
        "0.5" = [0.5]
        "0:1:.25" = [0.0, .25, .5, .75, 1.0]
    """
    import argparse
    range_split = range_str.split(':')
    if len(range_split) == 1:
        return [float(range_split[0])]
    if len(range_split) == 2:
        range_split += [1]  # default step size 1
    if len(range_split) != 3:
        raise argparse.ArgumentError('Range improperly formatted')
    start = float(range_split[0])
    stop = float(range_split[1])
    if start < 0:
        raise argparse.ArgumentTypeError('Range is under 0')
    if stop > 1:
        raise argparse.ArgumentTypeError('Range is above 1')
    step_size = float(range_split[2])
    n_steps = round((stop - start) / step_size, 6)
    if abs(int(n_steps) - n_steps) > 0:
        raise argparse.ArgumentError('Range does not fit')
    return [(start + (i*step_size)) for i in range(int(n_steps)+1)]


def parse_log_range(ctx, param, range_str):
    """Parse a string that can be interpreted as a log range.
    0 in linear space is always included. Includes endpoint. Step size needs to
    fit.

    For example:
        "-2:2:1" = [0, 2**-2, 2**-1, 1, 2, 2**2, 2**3]
        "2:3:.25" = [0, 2**2, 2**2.25, 2**2.5, 2**2.75, 2**3]
    """
    import argparse
    range_split = range_str.split(':')
    if len(range_split) == 1:
        if range_split[0] == 'z':
            return [0]
        else:
            return [2**float(range_split[0])]
    if len(range_split) == 2:
        range_split += [1]  # default step size 1
    if len(range_split) != 3:
        raise argparse.ArgumentError('Range improperly formatted')
    start = float(range_split[0])
    stop = float(range_split[1])
    step_size = float(range_split[2])
    n_steps = (stop - start) / step_size
    if abs(int(n_steps) - n_steps) > 0:
        raise argparse.ArgumentError('Range does not fit')
    return [0] + [2**(start + (i*step_size)) for i in range(int(n_steps)+1)]


def write_result(result, group, zlib=True):
    group.createDimension('factor', result.params.k)

    for param_name, param_val in zip(Params._fields, result.params):
        group.setncattr(param_name, param_val)

    if result.coefficients is not None:
        gexp_coef = group.createVariable('coefficient_gexp', 'f8',
                                         ('factor', 'gene'), zlib=zlib)
        gexp_coef[:, :] = result.coefficients[0]

        mri_coef = group.createVariable('coefficient_mri', 'f8',
                                        ('factor', 'cad_feature'), zlib=zlib)
        mri_coef[:] = result.coefficients[1].T

    if result.factors is not None:
        factors = group.createVariable('factor_value', 'f8',
                                       ('case', 'factor'), zlib=zlib)
        factors[:] = result.factors

    if result.monitor is not None:
        m_group = group.createGroup('monitor')
        m_group.createDimension('iteration', len(result.monitor['iteration']))
        iteration = m_group.createVariable('iteration', 'i4', ('iteration', ),
                                           zlib=zlib)
        iteration[:] = result.monitor['iteration']

        for m_name, m_value in result.monitor.items():
            if m_name == 'iteration':
                continue
            m_array = np.array(m_value)
            m_var = m_group.createVariable(m_name, m_array.dtype,
                                           ('iteration', ), zlib=zlib)
            m_var[:] = m_array
            m_array = np.array(m_value)


click_in_path = click.Path(exists=True, dir_okay=False, resolve_path=True)
click_out_path = click.Path(exists=False, dir_okay=False, resolve_path=True)


@click.command()
@click.argument('data', type=click_in_path)
@click.argument('out', type=click_out_path)
@click.option('--k', default="2", callback=parse_int_range)
@click.option('--alpha', default="0.5", callback=parse_zero_one_range)
@click.option('--l-gexp', default="1.0", callback=parse_log_range)
@click.option('--l-mri', default="1.0", callback=parse_log_range)
@click.option('--max-iter', default=1000)
@click.option('--eps', default=1e-6)
@click.option('--threads', default=1)
@click_log.simple_verbosity_option()
@click_log.init(__name__)
def run_sfa(data, out, k, alpha, l_gexp, l_mri, max_iter, eps, threads):
    data = funcsfa.StackedDataMatrix.from_netcdf(data)

    logger.info("Sweeping k: {}".format(", ".join([str(i) for i in k])))
    logger.info("Sweeping alpha: {}".format(", ".join(str(a) for a in alpha)))
    logger.info("Sweeping l_gexp: {}".format(
        ", ".join(str(l) for l in l_gexp)))
    logger.info("Sweeping l_mri: {}".format(
        ", ".join(str(l) for l in l_mri)))

    with ThreadPoolExecutor(threads) as executor:
        futures = list()
        for param in (Params(*p) for p in product(k, alpha, l_gexp, l_mri)):
            fut = executor.submit(sfa, data, param, max_iter, eps)
            futures.append(fut)

        out_ds = netCDF4.Dataset(out, 'w', format='NETCDF4')
        models = out_ds.createGroup("models")
        out_ds.createDimension('case', len(data.samples))
        case = out_ds.createVariable('case', 'i8', ('case', ))
        case[:] = data.samples
        out_ds.createDimension('gene', data.dt_n_features[0])
        gene = out_ds.createVariable('gene', str, ('gene', ))
        for i, g in enumerate(data.features[data.slices[0]]):
            gene[i] = g
        out_ds.createDimension('cad_feature', data.dt_n_features[1])
        cad_feature = out_ds.createVariable('cad_feature', str,
                                            ('cad_feature', ))
        for i, f in enumerate(data.features[data.slices[1]]):
            cad_feature[i] = f

        n_results = len(futures)
        for i, fut in enumerate(as_completed(futures)):
            res = fut.result()
            logger.info("({:.2%}) Writing results of {}".format(
                i/n_results,
                res.params,
            ))
            group_name = "M_" + "_".join([str(p) for p in res.params])
            group = models.createGroup(group_name)
            write_result(res, group, zlib=True)

        out_ds.close()


if __name__ == '__main__':
    run_sfa()
