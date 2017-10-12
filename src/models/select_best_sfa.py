import click
import netCDF4
import xarray as xr


def nc_copy(g, h, recursive=False):
    for dim in g.dimensions.values():
        h.createDimension(dim.name, dim.size)

    for var in g.variables.values():
        h_var = h.createVariable(var.name, var.datatype, var.dimensions)
        h_var[:] = var[:]

    if recursive:
        for group in g.groups.values():
            nc_copy(group, h.createGroup(group.name))


@click.command()
@click.argument('models', type=click.Path(exists=True, dir_okay=False,
                                          resolve_path=True))
@click.argument('bics', type=click.Path(exists=True, dir_okay=False,
                                        resolve_path=True))
@click.argument('out', type=click.Path(exists=False, dir_okay=False,
                                       writable=True, resolve_path=True))
def select_best(models, bics, out):
    bics_ds = xr.open_dataset(bics).load()

    max_iter = bics_ds['n_iter'].max()
    bics_ds_converged = bics_ds.sel(model=bics_ds['n_iter'] < max_iter)
    sel_model_idx = int(bics_ds_converged['bic'].argmin())
    sel_model = str(bics_ds_converged['model'][sel_model_idx].values)

    with netCDF4.Dataset(models) as models, netCDF4.Dataset(out, 'w') as out:
        nc_copy(models, out, recursive=False)
        model = models['models'][sel_model]
        nc_copy(model, out, recursive=True)


if __name__ == '__main__':
    select_best()
