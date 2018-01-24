import click
import matplotlib
import numpy as np
import sklearn
import xarray as xr

import plot

from features.fa_mri_features import read_mri, adjust_scale


click_in_path = click.Path(exists=True, dir_okay=False, resolve_path=True)
click_out_path = click.Path(exists=False, dir_okay=False, resolve_path=True)


@click.command()
@click.argument('mri_features', type=click_in_path)
@click.argument('out', type=click_out_path)
def plot_fa_variance_explained(mri_features, out):
    mri_data_set = xr.open_dataset(mri_features).load()
    mri = read_mri(mri_data_set)
    mri = adjust_scale(mri)

    pca = sklearn.decomposition.PCA()
    mri_a = (mri / mri.std('case')).values
    pca.fit(mri_a)

    total_var = np.sum(pca.explained_variance_)

    with plot.subplots(3, 1, sharex=True, figsize=(3, 3)) as (fig, axs):

        axs[0].plot(
            range(1, len(pca.explained_variance_)+1),
            np.cumsum(pca.explained_variance_) / total_var,
            clip_on=False, zorder=100,
        )
        axs[0].set_ylim(top=1.0)
        axs[0].set_ylabel("$\Sigma$EV")
        axs[1].plot(
            range(1, len(pca.explained_variance_)+1),
            pca.explained_variance_ / total_var,
            clip_on=False, zorder=100,
        )
        axs[1].set_ylim(bottom=0.0)
        axs[1].set_ylabel("EV")
        axs[2].plot(
            range(2, len(pca.explained_variance_)+1),
            abs(pca.explained_variance_[1:] -
                pca.explained_variance_[:-1]) /
            total_var,
            clip_on=False, zorder=100,
        )
        axs[2].set_ylim(bottom=0.0)
        axs[2].set_ylabel("$\Delta$EV")

        for ax in axs:
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['left'].set_position(('outward', 10))

        axs[2].set_xlabel("Principal Component")
        axs[2].set_xlim(1, len(pca.explained_variance_))
        axs[2].set_xticks([1, 5, 10, 15, 20, len(pca.explained_variance_)])

        fig.savefig(out, format="svg")


if __name__ == '__main__':
    matplotlib.rcParams['font.family'] = 'sans-serif'
    matplotlib.rcParams['font.sans-serif'] = ['Arial']
    matplotlib.rcParams['font.weight'] = 'normal'
    matplotlib.rcParams['figure.dpi'] = 300
    matplotlib.rcParams['figure.facecolor'] = 'none'
    matplotlib.rcParams['axes.labelsize'] = 11
    matplotlib.rcParams['xtick.labelsize'] = 9
    matplotlib.rcParams['ytick.labelsize'] = 9

    plot_fa_variance_explained()
