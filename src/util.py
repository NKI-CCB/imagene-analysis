from itertools import product

import numpy as np
import scipy.stats
import xarray as xr


def swivel_dim(x, dim):
    dims_order = [dim] + [i for i in range(len(x.shape)) if i != dim]
    return x.transpose(*dims_order)


_cor_funs = {
    'pearson': scipy.stats.pearsonr,
    'spearman': scipy.stats.spearmanr,
}


def cor(x, y, dim=0, *, nan_policy='omit', method='pearson'):
    cor_fun = _cor_funs.get(method, None)
    if cor_fun is None:
        raise ValueError("cor_fun must be one of {" +
                         ",".join(_cor_funs.keys()) + "}")

    if isinstance(dim, str):
        x_dim = x.dims.index(dim)
        y_dim = y.dims.index(dim)
    else:
        x_dim = dim
        y_dim = dim
        if hasattr(x, 'dims'):
            assert x.dims[dim] == y.dims[dim]
            dim = x.dims[x_dim]
        else:
            dim = 'sample'

    if hasattr(x, 'dims'):
        fdim_x = [d for d in x.dims if d != dim]
    else:
        if len(x.dim) == 2:
            fdim_x = ['x']
        else:
            fdim_x = ['x{}'.format(i) for i in range(len(x.dim)-1)]
    if hasattr(y, 'dims'):
        fdim_y = [d for d in y.dims if d != dim]
    else:
        if len(y.dim) == 2:
            fdim_y = ['y']
        else:
            fdim_y = ['y{}'.format(i) for i in range(len(y.dim)-1)]
    assert all([d not in fdim_x for d in fdim_y])

    coords = dict()
    if hasattr(x, 'coords'):
        for d in fdim_x:
            assert d not in coords
            coords[d] = x.coords[d]
    if hasattr(y, 'coords'):
        for d in fdim_y:
            assert d not in coords
            coords[d] = y.coords[d]

    x_a = np.asarray(x)
    x_a = swivel_dim(x_a, x_dim)
    y_a = np.asarray(y)
    y_a = swivel_dim(y_a, y_dim)

    fshape_x = x_a.shape[1:]
    fshape_y = y_a.shape[1:]
    cor_a = np.full(fshape_x + fshape_y, np.nan)
    p_a = np.full(fshape_x + fshape_y, np.nan)
    fiter_x = product(*[range(n) for n in fshape_x])
    fiter_y = product(*[range(n) for n in fshape_y])

    for x_i, y_i in product(fiter_x, fiter_y):
        x_slice = x_a[(slice(None),) + x_i]
        y_slice = y_a[(slice(None),) + y_i]
        if nan_policy == 'omit':
            not_na_idx = np.isfinite(x_slice) & np.isfinite(y_slice)
            x_slice = x_slice[not_na_idx]
            y_slice = y_slice[not_na_idx]
        elif nan_policy == 'propagate':
            pass
        else:
            raise ValueError("nan_policy must be one of {'propagate', 'omit'}")
        r, p = cor_fun(x_slice, y_slice)
        cor_a[x_i + y_i] = r
        p_a[x_i + y_i] = p

    p_da = xr.DataArray(p_a, dims=fdim_x+fdim_y)
    cor_da = xr.DataArray(cor_a, dims=fdim_x+fdim_y)
    return xr.Dataset({'correlation': cor_da, 'nominal_p': p_da},
                      coords=coords)
