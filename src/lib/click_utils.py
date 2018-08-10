import click

in_path = click.Path(exists=True, dir_okay=False, resolve_path=True)
out_path = click.Path(exists=False, dir_okay=False, resolve_path=True)
