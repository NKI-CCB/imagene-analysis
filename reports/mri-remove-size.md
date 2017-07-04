---
title: Exploration to remove size from MRI CAD features
author: Tycho Bismeijer
date: 2017-07-03
---

## Setup ## {.collapsed}

Load libraries.


```python
from math import sqrt, ceil
import sys

from IPython.display import display, Markdown
import numpy as np
import scipy.stats
import sklearn.decomposition
import xarray as xr

sys.path.append('../src/lib/')
import plot
```



Setup style of plots.


```python
import matplotlib
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Alegreya Sans']
matplotlib.rcParams['font.weight'] = 'regular'
matplotlib.rcParams['figure.dpi'] = 300
```




```python
mri_ds = xr.open_dataset('../data/processed/mri-features.nc')
del mri_ds['Comment']
del mri_ds['MultiFocal']
mri = mri_ds.to_array('cad_feature', 'mri_cad_features')
mri_ds.close()
mri = mri.isel(case=np.where(mri.isnull().sum('cad_feature') == 0)[0])
mri = mri.transpose('case', 'cad_feature')
print(mri)
```

```
<xarray.DataArray 'mri_cad_features' (case: 282, cad_feature: 24)>
array([[  8.343770e-01,   4.046130e-01,   7.845880e+02, ...,
3.197430e+00,
          5.604590e-01,   3.728000e-01],
       [  7.080830e-01,   5.653050e-01,   1.351080e+04, ...,
1.607050e+00,
          5.389980e-01,   2.770000e-01],
       [  8.127700e-01,   5.227460e-01,   5.224870e+03, ...,
3.299920e+00,
          5.780610e-01,   4.661000e-01],
       ...,
       [  8.156050e-01,   4.480070e-01,   4.401050e+03, ...,
3.501010e+00,
          5.770520e-01,   2.941000e-01],
       [  7.208790e-01,   6.401100e-01,   1.725120e+04, ...,
3.899760e+00,
          5.887040e-01,   5.646000e-01],
       [  6.442930e-01,   4.969370e-01,   1.826620e+03, ...,
1.412870e+00,
          5.200650e-01,   5.387000e-01]])
Coordinates:
  * case         (case) int64 192 196 199 207 208 217 219 269 273 274
288 ...
  * cad_feature  (cad_feature) <U35 'circularity' 'irregularity'
'volume' ...
Attributes:
    title: MRI features from Margins of samples with gene expression
data from Imagene
    history: 2017-04-12T14:56:13+00:00 process_mri.py Converted from
/home/tycho/Projects/Imagene/mri-rnaseq-integration/data/raw/mri-
features.xlsx.
```



## Only linear scale adjustment ##

### PCA ###


```python
pca = sklearn.decomposition.PCA(10)
pca.fit((mri / mri.std('case')).values)
```

```
PCA(copy=True, iterated_power='auto', n_components=10,
random_state=None,
  svd_solver='auto', tol=0.0, whiten=False)
```




```python
plot.heatmap(
    pca.components_.T,
    xticklabels=[str(i+1) for i in range(pca.components_.shape[0])],
    yticklabels=mri.cad_feature.values,
)
```

![](figures/mri-remove-size_heatmap-pca_1.png){#heatmap-pca }\



```python
with plot.subplots(1, 1, figsize=(7, 7)) as (fig, ax):
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    for i in range(pca.components_.shape[1]):
        ax.arrow(0, 0, pca.components_[0, i], pca.components_[1, i])
        alength = sqrt(pca.components_[0, i]**2 + pca.components_[1, i]**2)
        if alength > 0.25:
            ax.text(pca.components_[0, i], pca.components_[1, i],
                    mri['cad_feature'].values[i])
```

![](figures/mri-remove-size_plot-pca-a_1.png){#plot-pca-a }\


### Normal QQ-plots of MRI CAD features ###


```python
nf = len(mri['cad_feature'])
nc = len(mri['case'])
n_row = ceil(nf / 3)
n_col = 3
fs = (n_col*3, n_row*3)
with plot.subplots(n_row, n_col, figsize=fs) as (fig, axs):
    norm = scipy.stats.norm(0, 1)
    q_norm = norm.ppf(np.linspace(1 / nc, 1 - (1 / nc), nc))
    for i in range(nf):
        row, col = divmod(i, n_col)
        q_feature = mri.isel(cad_feature=i).values
        plot.qqplot(
            q_norm, q_feature,
            ax=axs[row, col],
            xlabel='Normal Q',
            ylabel=str(mri['cad_feature'].values[i]) + " Q",
            title=""
        )
    fig.tight_layout()
```

![](figures/mri-remove-size_qqplot-mri-features_1.png){#qqplot-mri-features }\


## Adjusted Scale ##


```python
mri_adj = xr.DataArray(np.full(mri.shape, np.nan), mri.coords, mri.dims)
for feature in mri['cad_feature'].values:
    if feature[0:3] == 'vol':
        mri_adj.loc[:, feature] = np.cbrt(mri.loc[:, feature])
    elif feature[0:3] == 'var':
        mri_adj.loc[:, feature] = np.sqrt(mri.loc[:, feature])
    else:
        mri_adj.loc[:, feature] = mri.loc[:, feature]
mri_adj
```

```
<xarray.DataArray (case: 282, cad_feature: 24)>
array([[  0.834377,   0.404613,   9.223177, ...,   1.788136,
0.560459,
          0.3728  ],
       [  0.708083,   0.565305,  23.817364, ...,   1.267695,
0.538998,   0.277   ],
       [  0.81277 ,   0.522746,  17.352358, ...,   1.816568,
0.578061,
          0.4661  ],
       ...,
       [  0.815605,   0.448007,  16.387729, ...,   1.871099,
0.577052,
          0.2941  ],
       [  0.720879,   0.64011 ,  25.838845, ...,   1.974781,
0.588704,
          0.5646  ],
       [  0.644293,   0.496937,  12.224077, ...,   1.188642,
0.520065,
          0.5387  ]])
Coordinates:
  * case         (case) int64 192 196 199 207 208 217 219 269 273 274
288 ...
  * cad_feature  (cad_feature) <U35 'circularity' 'irregularity'
'volume' ...
```



### Normal QQ-plots of MRI CAD features ###


```python
nf = len(mri_adj['cad_feature'])
nc = len(mri_adj['case'])
n_row = ceil(nf / 3)
n_col = 3
fs = (n_col*3, n_row*3)
with plot.subplots(n_row, n_col, figsize=fs) as (fig, axs):
    norm = scipy.stats.norm(0, 1)
    q_norm = norm.ppf(np.linspace(1 / nc, 1 - (1 / nc), nc))
    for i in range(nf):
        row, col = divmod(i, n_col)
        q_feature = mri_adj.isel(cad_feature=i).values
        plot.qqplot(
            q_norm, q_feature,
            ax=axs[row, col],
            xlabel='Normal Q',
            ylabel=str(mri_adj['cad_feature'].values[i]) + " Q",
            title=""
        )
    fig.tight_layout()
```

![](figures/mri-remove-size_qqplot-mri-features-adj_1.png){#qqplot-mri-features-adj }\


### PCA ###


```python
pca_adj = sklearn.decomposition.PCA(10)
pca_adj.fit((mri_adj / mri_adj.std('case')).values)
```

```
PCA(copy=True, iterated_power='auto', n_components=10,
random_state=None,
  svd_solver='auto', tol=0.0, whiten=False)
```




```python
plot.heatmap(
    pca_adj.components_.T,
    xticklabels=[str(i+1) for i in range(pca_adj.components_.shape[0])],
    yticklabels=mri_adj.cad_feature.values,
)
```

![](figures/mri-remove-size_heatmap-pca-adj_1.png){#heatmap-pca-adj }\



```python
with plot.subplots(1, 1, figsize=(7, 7)) as (fig, ax):
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    for i in range(pca_adj.components_.shape[1]):
        ax.arrow(0, 0, pca_adj.components_[0, i], pca_adj.components_[1, i])
        alength = sqrt(pca_adj.components_[0, i]**2 +
                       pca_adj.components_[1, i]**2)
        if alength > 0.25:
            ax.text(pca_adj.components_[0, i], pca_adj.components_[1, i],
                    mri_adj['cad_feature'].values[i])
```

![](figures/mri-remove-size_plot-pca-adj-a_1.png){#plot-pca-adj-a }\

