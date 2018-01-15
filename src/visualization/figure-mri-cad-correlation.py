import click
import matplotlib
import numpy as np
import xarray as xr

import plot


click_in_path = click.Path(exists=True, dir_okay=False, resolve_path=True)
click_out_path = click.Path(exists=False, dir_okay=False, resolve_path=True)

feature_order = [
    'volume', 'largest_diameter',
    'vol_init_enhancement_GT100', 'ld_init_enhancement_GT100',
    'vol_late_LT0', 'ld_late_LT0',
    'mean_smoothness_uptake', 'variation_smoothness_uptake',
    'mean_smoothness_all_timeframes', 'variation_smoothness_all_timeframes',
    'mean_sharpness_uptake', 'var_sharpness_uptake',
    'mean_sharpness_all_timeframes', 'var_sharpness_all_timeframes',
    'uptake_speed',
    'top_init_enhancement', 'top_late_enhancement',
    'ser', 'washout',
    'mean_vox_val', 'variance_vox_val',
    'circularity',
    'irregularity',
    'PCE_top10percent',
]

feature_display_names = {
    'volume': "Volume",
    'largest_diameter': "Diameter",
    'vol_init_enhancement_GT100': "Volume Initial Enhancement > 100",
    'ld_init_enhancement_GT100': "Diameter Initial Enhancement > 100",
    'vol_late_LT0': "Volume Late Enhancement < 0",
    'ld_late_LT0': "Diameter Late Enhancement < 0",
    'mean_smoothness_uptake': "Smoothness Uptake (mean)",
    'variation_smoothness_uptake': "Smoothness Uptake (variation)",
    'mean_smoothness_all_timeframes': "Smoothness Maximum (mean)",
    'variation_smoothness_all_timeframes': "Smoothness Maximum (variation)",
    'mean_sharpness_uptake': "Sharpness Uptake (mean)",
    'var_sharpness_uptake': "Sharpness Uptake (variation)",
    'mean_sharpness_all_timeframes': "Sharpness Maximum (mean)",
    'var_sharpness_all_timeframes': "Sharpness Maximum (variation)",
    'uptake_speed': "Uptake Speed",
    'top_init_enhancement': "Top Initial Enhancement",
    'top_late_enhancement': "Top Late Enhancement",
    'ser': "SER",
    'washout': "Washout",
    'mean_vox_val': "Average Pre-Contrast Voxel Value",
    'variance_vox_val': "Variance Pre-Contrast Voxel Value",
    'circularity': "Circularity",
    'irregularity': "Irregularity",
    'PCE_top10percent': "PCE",
}


@click.command()
@click.argument('mri_features', type=click_in_path)
@click.argument('out', type=click_out_path)
def plot_mri_cad_factor_correlation(mri_features, out):
    mri_ds = xr.open_dataset(mri_features)
    del mri_ds['Comment']
    del mri_ds['MultiFocal']
    mri = mri_ds.to_array('cad_feature', 'mri_cad_features')
    mri_ds.close()
    mri = mri.isel(case=np.where(mri.isnull().sum('cad_feature') == 0)[0])
    mri = mri.transpose('case', 'cad_feature')

    assert all(f.item() in feature_order for f in mri['cad_feature'].values)
    mri = mri.reindex(cad_feature=feature_order)

    cor = xr.DataArray(
        data=np.corrcoef(mri.values, rowvar=False),
        dims=('cad_feature', 'cad_feature'),
        coords={'cad_feature': mri.coords['cad_feature']},
    )
    cor.name = 'correlation'

    with plot.subplots(figsize=(3.5, 3.5)) as (fig, ax):
        c = plot.heatmap(
            cor,
            mask=np.tri(cor.shape[0]) < 0.5,
            aspect='equal',
            xlabel="CAD Feature", ylabel="CAD Feature",
            xticklabels=['']*cor.shape[0],
            yticklabels=[feature_display_names[f]
                         for f in cor.coords['cad_feature'].values],
            cbar=False,
            ax=ax,
        )
        cax = fig.add_axes([.7, .5, .05, .25])
        cbar = fig.colorbar(c, cax=cax)
        cbar.set_ticks([-1.0, -0.5, 0.0, 0.5, 1.0])
        cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=9)
        cbar.ax.set_title("Pearson Correlation", fontsize=10)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
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

    plot_mri_cad_factor_correlation()
