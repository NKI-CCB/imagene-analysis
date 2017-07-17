---
title: Gene Set Enrichment Analysis of CAD features
author: Tycho Bismeijer
date: 2017-06-06
---

## Setup ## {.collapsed}

Load libraries.


```python
from IPython.display import display, Markdown
import numpy as np
import pandas as pd
from tabulate import tabulate
import xarray as xr

import plot
```



Setup style of plots.


```python
import matplotlib
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Alegreya Sans']
matplotlib.rcParams['font.weight'] = 'regular'
matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['figure.figsize'] = (5.2, 5.2)

```



Load gene set enrichment analysis (GSEA) results on the MsigDB Hallmarks set.


```python
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

```




```python
df_h = load_gsea_ds("../analyses/gsea/mri-features_h.all_T.nc")
```



A function to plot a heatmap with normalized enrichment statistics (NES) that
are significant.

A function to plot a heatmap with normalized enrichment statistics (NES) that
are significant.


```python
def plot_ds(ds, fdr, le_prop=0.0, abs=True):
    ds = ds.copy()
    ds['significance_mask'] = (ds['fdr'] > fdr) | (ds['le_prop'] < le_prop)
    ds = ds.sel(
        gene_set=np.logical_not(ds['significance_mask'])
                 .sum('mri_feature') > 0,
    )
    if abs:
        zlim = [0, np.max(ds['nes'])]
        cmap='viridis'
    else:
        zlim = np.max(np.abs(ds['nes']))
        zlim = [-zlim, zlim]
        cmap = 'coolwarm'

    with plot.subplots(1, 1) as (fig, ax):
        plot.heatmap(
            ds['nes'], mask=ds['significance_mask'],
            zlim=zlim, cmap=cmap,
            row_dendrogram=True, col_dendrogram=True,
            ax=ax,
        )
        if len(ds['gene_set']) < 50:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        else:
            ax.set_xticklabels("")

```



A utility function to display tables with [DataTables](https://datatables.net).
And a function to show significant pathway enrichments in a table.


```python
def display_table(t):
    display(Markdown(
        '<div class="datatable">' +
        tabulate(t, headers='keys') +
        "\n\n</div>"
    ))
def table_ds(ds, fdr, le_prop=0.0):
    df = ds.to_dataframe()
    df.reset_index(level=0, inplace=True)
    display_table(df.loc[(df['fdr'] < fdr) & (df['le_prop'] > le_prop)])

```




## Results ##

Hallmarks from MSigDB.


```python
plot_ds(df_h, fdr=0.25)
```

![](figures/gsea_hallmarks-es-plot_1){#hallmarks-es-plot }\



```python
table_ds(df_h, fdr=0.25)
```


<div class="datatable">mri_feature                 gene_set                                          es            p      nes        fdr       fwer    max_es_at    le_prop  gene_set_code
--------------------------  ------------------------------------------  --------  -----------  -------  ---------  ---------  -----------  ---------  ------------------------------------------
Largest Diameter            Hallmark Tnfa Signaling Via Nfkb            0.423123  0.0728927    1.30912  0.170349   1                 2662   0.380952  HALLMARK_TNFA_SIGNALING_VIA_NFKB
Variance Vox Val            Hallmark Tnfa Signaling Via Nfkb            0.424242  0.070493     1.31151  0.215989   1                 3380   0.440476  HALLMARK_TNFA_SIGNALING_VIA_NFKB
Ld Init Enhancement Gt100   Hallmark Tnfa Signaling Via Nfkb            0.431836  0.0614939    1.33536  0.15937    1                 3595   0.470238  HALLMARK_TNFA_SIGNALING_VIA_NFKB
Ld Late Lt0                 Hallmark Tnfa Signaling Via Nfkb            0.401436  0.114489     1.24356  0.229609   1                 4354   0.52381   HALLMARK_TNFA_SIGNALING_VIA_NFKB
Irregularity                Hallmark Hypoxia                            0.436426  0.010499     1.33794  0.195698   0.524948          3576   0.490909  HALLMARK_HYPOXIA
Largest Diameter            Hallmark Hypoxia                            0.47027   0.00139986   1.44223  0.106688   0.069993          3847   0.551515  HALLMARK_HYPOXIA
Ld Init Enhancement Gt100   Hallmark Hypoxia                            0.456109  0.00319968   1.3984   0.141125   0.159984          4145   0.569697  HALLMARK_HYPOXIA
Ld Late Lt0                 Hallmark Hypoxia                            0.441323  0.00919908   1.35359  0.210174   0.459954          4529   0.593939  HALLMARK_HYPOXIA
Largest Diameter            Hallmark Cholesterol Homeostasis            0.42097   0.0821918    1.23775  0.196002   1                 4959   0.615385  HALLMARK_CHOLESTEROL_HOMEOSTASIS
Ld Init Enhancement Gt100   Hallmark Cholesterol Homeostasis            0.426082  0.070493     1.25155  0.200161   1                 4229   0.538462  HALLMARK_CHOLESTEROL_HOMEOSTASIS
Ld Late Lt0                 Hallmark Cholesterol Homeostasis            0.432656  0.0526       1.27417  0.229609   1                 4319   0.538462  HALLMARK_CHOLESTEROL_HOMEOSTASIS
Irregularity                Hallmark Mitotic Spindle                    0.468216  0.00519948   1.45198  0.182654   0.259974          3192   0.416667  HALLMARK_MITOTIC_SPINDLE
Volume                      Hallmark Mitotic Spindle                    0.458632  0.00959904   1.42425  0.169069   0.479952          1315   0.255208  HALLMARK_MITOTIC_SPINDLE
Largest Diameter            Hallmark Mitotic Spindle                    0.479571  0.0029997    1.49148  0.100434   0.149985          2210   0.317708  HALLMARK_MITOTIC_SPINDLE
Variance Vox Val            Hallmark Mitotic Spindle                    0.449459  0.0132987    1.3971   0.175313   0.664934          2502   0.338542  HALLMARK_MITOTIC_SPINDLE
Vol Init Enhancement Gt100  Hallmark Mitotic Spindle                    0.454314  0.010199     1.41175  0.162906   0.509949          1311   0.255208  HALLMARK_MITOTIC_SPINDLE
Ld Init Enhancement Gt100   Hallmark Mitotic Spindle                    0.482646  0.00249975   1.5014   0.11648    0.124988          2068   0.307292  HALLMARK_MITOTIC_SPINDLE
Ld Late Lt0                 Hallmark Mitotic Spindle                    0.488948  0.00149985   1.51969  0.104714   0.0749925         1274   0.255208  HALLMARK_MITOTIC_SPINDLE
Irregularity                Hallmark Tgf Beta Signaling                 0.481729  0.0166066    1.40965  0.187118   0.830332          3252   0.458333  HALLMARK_TGF_BETA_SIGNALING
Largest Diameter            Hallmark Tgf Beta Signaling                 0.467622  0.026508     1.37199  0.155047   1                 2087   0.375     HALLMARK_TGF_BETA_SIGNALING
Variance Vox Val            Hallmark Tgf Beta Signaling                 0.488881  0.0113011    1.43582  0.152756   0.565057          3126   0.479167  HALLMARK_TGF_BETA_SIGNALING
Ld Init Enhancement Gt100   Hallmark Tgf Beta Signaling                 0.457773  0.0372074    1.3414   0.15937    1                 4020   0.5625    HALLMARK_TGF_BETA_SIGNALING
Ld Late Lt0                 Hallmark Tgf Beta Signaling                 0.438461  0.0635318    1.28793  0.229609   1                 3444   0.479167  HALLMARK_TGF_BETA_SIGNALING
Largest Diameter            Hallmark Il6 Jak Stat3 Signaling            0.426439  0.1327       1.24472  0.196002   1                 3686   0.442623  HALLMARK_IL6_JAK_STAT3_SIGNALING
Ld Init Enhancement Gt100   Hallmark Il6 Jak Stat3 Signaling            0.413739  0.1664       1.20931  0.221814   1                 3612   0.42623   HALLMARK_IL6_JAK_STAT3_SIGNALING
Largest Diameter            Hallmark Dna Repair                         0.413265  0.0241976    1.2874   0.187021   1                 4381   0.514925  HALLMARK_DNA_REPAIR
Ld Init Enhancement Gt100   Hallmark Dna Repair                         0.406111  0.0340966    1.26466  0.199533   1                 3960   0.470149  HALLMARK_DNA_REPAIR
Ld Late Lt0                 Hallmark Dna Repair                         0.398047  0.0455954    1.2387   0.229609   1                 3969   0.470149  HALLMARK_DNA_REPAIR
Irregularity                Hallmark G2M Checkpoint                     0.65688   0.0009999    2.01081  0.0174147  0.049995          2008   0.6       HALLMARK_G2M_CHECKPOINT
Volume                      Hallmark G2M Checkpoint                     0.659236  0.00059994   2.02399  0.0164795  0.029997          1991   0.594444  HALLMARK_G2M_CHECKPOINT
Largest Diameter            Hallmark G2M Checkpoint                     0.661803  0.00039996   2.0377   0.0136601  0.019998          1274   0.522222  HALLMARK_G2M_CHECKPOINT
Mean Vox Val                Hallmark G2M Checkpoint                     0.63622   0.00189981   1.95806  0.0258701  0.0949905         1209   0.483333  HALLMARK_G2M_CHECKPOINT
Variance Vox Val            Hallmark G2M Checkpoint                     0.63118   0.0019998    1.94054  0.0273211  0.09999           1454   0.5       HALLMARK_G2M_CHECKPOINT
Top Late Enhancement        Hallmark G2M Checkpoint                     0.54457   0.030497     1.6739   0.181895   1                 2006   0.444444  HALLMARK_G2M_CHECKPOINT
Vol Init Enhancement Gt100  Hallmark G2M Checkpoint                     0.657133  0.00049995   2.01328  0.0173137  0.0249975         2029   0.6       HALLMARK_G2M_CHECKPOINT
Ld Init Enhancement Gt100   Hallmark G2M Checkpoint                     0.659652  0.0009999    2.03168  0.0152105  0.049995          1056   0.5       HALLMARK_G2M_CHECKPOINT
Vol Late Lt0                Hallmark G2M Checkpoint                     0.642938  0.00179982   1.97414  0.0244195  0.089991          1726   0.55      HALLMARK_G2M_CHECKPOINT
Ld Late Lt0                 Hallmark G2M Checkpoint                     0.676242  0.00019998   2.08124  0.0111076  0.009999          1589   0.566667  HALLMARK_G2M_CHECKPOINT
Largest Diameter            Hallmark Apoptosis                          0.427161  0.0119988    1.31617  0.170349   0.59994           4880   0.627737  HALLMARK_APOPTOSIS
Vol Init Enhancement Gt100  Hallmark Apoptosis                          0.443099  0.00539946   1.36192  0.203981   0.269973          3048   0.437956  HALLMARK_APOPTOSIS
Ld Init Enhancement Gt100   Hallmark Apoptosis                          0.438886  0.00579942   1.35351  0.15937    0.289971          3303   0.481752  HALLMARK_APOPTOSIS
Ld Late Lt0                 Hallmark Apoptosis                          0.420238  0.0163984    1.29461  0.229609   0.819918          4255   0.569343  HALLMARK_APOPTOSIS
Largest Diameter            Hallmark Adipogenesis                       0.393601  0.0345965    1.23139  0.196002   1                 3604   0.432432  HALLMARK_ADIPOGENESIS
Ld Init Enhancement Gt100   Hallmark Adipogenesis                       0.401419  0.0205979    1.25697  0.200161   1                 3635   0.443243  HALLMARK_ADIPOGENESIS
Ld Late Lt0                 Hallmark Adipogenesis                       0.399633  0.0242976    1.25009  0.229609   1                 4249   0.502703  HALLMARK_ADIPOGENESIS
Largest Diameter            Hallmark Myogenesis                         0.432089  0.0414959    1.31583  0.170349   1                 3888   0.496063  HALLMARK_MYOGENESIS
Variance Vox Val            Hallmark Myogenesis                         0.426804  0.0478952    1.30248  0.215989   1                 3970   0.519685  HALLMARK_MYOGENESIS
Ld Init Enhancement Gt100   Hallmark Myogenesis                         0.439666  0.0322968    1.33822  0.15937    1                 4634   0.582677  HALLMARK_MYOGENESIS
Variance Vox Val            Hallmark Interferon Alpha Response          0.51592   0.0926707    1.6048   0.0876278  1                 2360   0.5       HALLMARK_INTERFERON_ALPHA_RESPONSE
Variance Vox Val            Hallmark Interferon Gamma Response          0.487252  0.0721928    1.53632  0.109285   1                 3332   0.547619  HALLMARK_INTERFERON_GAMMA_RESPONSE
Ld Init Enhancement Gt100   Hallmark Unfolded Protein Response          0.39327   0.070293     1.21684  0.219236   1                 3165   0.376147  HALLMARK_UNFOLDED_PROTEIN_RESPONSE
Ld Late Lt0                 Hallmark Unfolded Protein Response          0.40635   0.0458954    1.2573   0.229609   1                 2967   0.376147  HALLMARK_UNFOLDED_PROTEIN_RESPONSE
Irregularity                Hallmark Mtorc1 Signaling                   0.438523  0.00989901   1.37591  0.187118   0.494951          2936   0.417989  HALLMARK_MTORC1_SIGNALING
Volume                      Hallmark Mtorc1 Signaling                   0.480353  0.0019998    1.50236  0.139425   0.09999           1996   0.375661  HALLMARK_MTORC1_SIGNALING
Largest Diameter            Hallmark Mtorc1 Signaling                   0.467799  0.00259974   1.46739  0.100434   0.129987          3089   0.460317  HALLMARK_MTORC1_SIGNALING
Vol Init Enhancement Gt100  Hallmark Mtorc1 Signaling                   0.480969  0.00169983   1.50148  0.13906    0.0849915         2118   0.386243  HALLMARK_MTORC1_SIGNALING
Ld Init Enhancement Gt100   Hallmark Mtorc1 Signaling                   0.469849  0.00209979   1.47343  0.120396   0.10499           3484   0.502646  HALLMARK_MTORC1_SIGNALING
Vol Late Lt0                Hallmark Mtorc1 Signaling                   0.461441  0.00489951   1.44115  0.200927   0.244976          2205   0.375661  HALLMARK_MTORC1_SIGNALING
Ld Late Lt0                 Hallmark Mtorc1 Signaling                   0.494309  0.00079992   1.55067  0.102153   0.039996          2967   0.470899  HALLMARK_MTORC1_SIGNALING
Irregularity                Hallmark E2F Targets                        0.661312  0.0025       2.03679  0.0174147  0.125             2165   0.65      HALLMARK_E2F_TARGETS
Volume                      Hallmark E2F Targets                        0.683821  0.0006       2.11311  0.0164795  0.03              1597   0.611111  HALLMARK_E2F_TARGETS
Largest Diameter            Hallmark E2F Targets                        0.695368  0.0004       2.15754  0.0107079  0.02              1717   0.638889  HALLMARK_E2F_TARGETS
Mean Vox Val                Hallmark E2F Targets                        0.62874   0.00539946   1.94596  0.0258701  0.269973          2376   0.638889  HALLMARK_E2F_TARGETS
Variance Vox Val            Hallmark E2F Targets                        0.645477  0.00279972   1.99707  0.0273211  0.139986          2846   0.694444  HALLMARK_E2F_TARGETS
Top Late Enhancement        Hallmark E2F Targets                        0.531999  0.0480952    1.64419  0.181895   1                 2788   0.516667  HALLMARK_E2F_TARGETS
Vol Init Enhancement Gt100  Hallmark E2F Targets                        0.666194  0.00160016   2.05337  0.0173137  0.080008          1586   0.594444  HALLMARK_E2F_TARGETS
Ld Init Enhancement Gt100   Hallmark E2F Targets                        0.685042  0.00059994   2.12568  0.0142098  0.029997          1634   0.611111  HALLMARK_E2F_TARGETS
Vol Late Lt0                Hallmark E2F Targets                        0.670707  0.0014       2.07074  0.0236188  0.07              1934   0.622222  HALLMARK_E2F_TARGETS
Ld Late Lt0                 Hallmark E2F Targets                        0.697474  0.00019998   2.1605   0.0111076  0.009999          1589   0.616667  HALLMARK_E2F_TARGETS
Irregularity                Hallmark Myc Targets V1                     0.49056   0.013        1.58026  0.102426   0.65              3092   0.507772  HALLMARK_MYC_TARGETS_V1
Volume                      Hallmark Myc Targets V1                     0.536135  0.0034       1.71899  0.0594213  0.17              2563   0.544041  HALLMARK_MYC_TARGETS_V1
Largest Diameter            Hallmark Myc Targets V1                     0.548356  0.00229977   1.76789  0.0345054  0.114989          3059   0.595855  HALLMARK_MYC_TARGETS_V1
Mean Vox Val                Hallmark Myc Targets V1                     0.44926   0.0346965    1.44659  0.249834   1                 2571   0.409326  HALLMARK_MYC_TARGETS_V1
Variance Vox Val            Hallmark Myc Targets V1                     0.410154  0.0860914    1.31991  0.215989   1                 1979   0.310881  HALLMARK_MYC_TARGETS_V1
Vol Init Enhancement Gt100  Hallmark Myc Targets V1                     0.526336  0.00449955   1.68317  0.0609682  0.224978          3178   0.601036  HALLMARK_MYC_TARGETS_V1
Ld Init Enhancement Gt100   Hallmark Myc Targets V1                     0.532298  0.00269973   1.71756  0.046412   0.134987          3045   0.569948  HALLMARK_MYC_TARGETS_V1
Vol Late Lt0                Hallmark Myc Targets V1                     0.563096  0.001        1.80364  0.0368544  0.05              2812   0.601036  HALLMARK_MYC_TARGETS_V1
Ld Late Lt0                 Hallmark Myc Targets V1                     0.558002  0.00179982   1.79931  0.0365     0.089991          2767   0.57513   HALLMARK_MYC_TARGETS_V1
Irregularity                Hallmark Myc Targets V2                     0.471714  0.0917284    1.38255  0.187118   1                 3615   0.581818  HALLMARK_MYC_TARGETS_V2
Volume                      Hallmark Myc Targets V2                     0.672466  0.000600661  1.96549  0.0164795  0.030033          2960   0.818182  HALLMARK_MYC_TARGETS_V2
Largest Diameter            Hallmark Myc Targets V2                     0.616533  0.00460599   1.81457  0.032599   0.230299          2487   0.654545  HALLMARK_MYC_TARGETS_V2
Vol Init Enhancement Gt100  Hallmark Myc Targets V2                     0.628815  0.00350631   1.83813  0.0394312  0.175316          2741   0.727273  HALLMARK_MYC_TARGETS_V2
Ld Init Enhancement Gt100   Hallmark Myc Targets V2                     0.609964  0.00550661   1.79503  0.0375259  0.27533           3107   0.709091  HALLMARK_MYC_TARGETS_V2
Vol Late Lt0                Hallmark Myc Targets V2                     0.65343   0.00140154   1.90539  0.0255537  0.0700771         2093   0.672727  HALLMARK_MYC_TARGETS_V2
Ld Late Lt0                 Hallmark Myc Targets V2                     0.637843  0.00240216   1.87667  0.0314548  0.120108          2533   0.690909  HALLMARK_MYC_TARGETS_V2
Irregularity                Hallmark Epithelial Mesenchymal Transition  0.510312  0.06         1.61505  0.102426   1                 3495   0.573034  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
Volume                      Hallmark Epithelial Mesenchymal Transition  0.532664  0.0468953    1.68298  0.0595464  1                 2802   0.55618   HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
Largest Diameter            Hallmark Epithelial Mesenchymal Transition  0.600656  0.0121988    1.90182  0.0254187  0.609939          2256   0.589888  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
Variance Vox Val            Hallmark Epithelial Mesenchymal Transition  0.543811  0.0364       1.73218  0.0609472  1                 3019   0.589888  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
Vol Init Enhancement Gt100  Hallmark Epithelial Mesenchymal Transition  0.530256  0.0474953    1.67756  0.0609682  1                 3154   0.58427   HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
Ld Init Enhancement Gt100   Hallmark Epithelial Mesenchymal Transition  0.592565  0.0158984    1.87685  0.0298206  0.794921          3104   0.674157  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
Vol Late Lt0                Hallmark Epithelial Mesenchymal Transition  0.484227  0.088        1.53229  0.139371   1                 3670   0.595506  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
Ld Late Lt0                 Hallmark Epithelial Mesenchymal Transition  0.512036  0.0588941    1.62101  0.0817359  1                 2671   0.52809   HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
Irregularity                Hallmark Glycolysis                         0.421999  0.0069993    1.29873  0.238278   0.349965          3228   0.427711  HALLMARK_GLYCOLYSIS
Largest Diameter            Hallmark Glycolysis                         0.439438  0.00209979   1.35222  0.163797   0.10499           3643   0.457831  HALLMARK_GLYCOLYSIS
Ld Init Enhancement Gt100   Hallmark Glycolysis                         0.423937  0.00609939   1.30417  0.175133   0.30497           3798   0.457831  HALLMARK_GLYCOLYSIS
Ld Late Lt0                 Hallmark Glycolysis                         0.417464  0.010499     1.28296  0.229609   0.524948          3959   0.46988   HALLMARK_GLYCOLYSIS
Largest Diameter            Hallmark Reactive Oxigen Species Pathway    0.433971  0.0905362    1.24934  0.196002   1                 4505   0.547619  HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY
Ld Init Enhancement Gt100   Hallmark Reactive Oxigen Species Pathway    0.431177  0.10013      1.2398   0.206088   1                 4924   0.595238  HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY
Largest Diameter            Hallmark P53 Pathway                        0.39948   0.0258974    1.23533  0.196002   1                 3804   0.446927  HALLMARK_P53_PATHWAY
Ld Init Enhancement Gt100   Hallmark P53 Pathway                        0.423824  0.00679932   1.30988  0.175133   0.339966          4000   0.486034  HALLMARK_P53_PATHWAY
Ld Late Lt0                 Hallmark P53 Pathway                        0.395802  0.0312969    1.22227  0.23661    1                 3000   0.368715  HALLMARK_P53_PATHWAY
Irregularity                Hallmark Uv Response Up                     0.441632  0.00369963   1.34892  0.195698   0.184982          2623   0.380952  HALLMARK_UV_RESPONSE_UP
Largest Diameter            Hallmark Uv Response Up                     0.428395  0.0110989    1.30855  0.170349   0.554945          4824   0.579365  HALLMARK_UV_RESPONSE_UP
Ld Init Enhancement Gt100   Hallmark Uv Response Up                     0.415952  0.0209979    1.26967  0.199533   1                 3141   0.404762  HALLMARK_UV_RESPONSE_UP
Ld Late Lt0                 Hallmark Uv Response Up                     0.418325  0.0177982    1.27562  0.229609   0.889911          4350   0.52381   HALLMARK_UV_RESPONSE_UP
Largest Diameter            Hallmark Uv Response Dn                     0.385683  0.0882912    1.19942  0.23067    1                 4822   0.530303  HALLMARK_UV_RESPONSE_DN
Ld Init Enhancement Gt100   Hallmark Uv Response Dn                     0.391873  0.069593     1.21871  0.219236   1                 5427   0.613636  HALLMARK_UV_RESPONSE_DN
Irregularity                Hallmark Angiogenesis                       0.50339   0.0879676    1.37015  0.187118   1                 2725   0.517241  HALLMARK_ANGIOGENESIS
Volume                      Hallmark Angiogenesis                       0.497823  0.0972433    1.35683  0.23306    1                 3694   0.62069   HALLMARK_ANGIOGENESIS
Largest Diameter            Hallmark Angiogenesis                       0.550272  0.0356855    1.49721  0.100434   1                 2029   0.482759  HALLMARK_ANGIOGENESIS
Mean Vox Val                Hallmark Angiogenesis                       0.527581  0.0578788    1.43756  0.249834   1                 4450   0.724138  HALLMARK_ANGIOGENESIS
Variance Vox Val            Hallmark Angiogenesis                       0.549025  0.036108     1.50207  0.115432   1                 3712   0.655172  HALLMARK_ANGIOGENESIS
Vol Init Enhancement Gt100  Hallmark Angiogenesis                       0.52425   0.0639952    1.43247  0.161177   1                 3896   0.689655  HALLMARK_ANGIOGENESIS
Ld Init Enhancement Gt100   Hallmark Angiogenesis                       0.53125   0.0506099    1.44813  0.126276   1                 3119   0.586207  HALLMARK_ANGIOGENESIS
Ld Late Lt0                 Hallmark Angiogenesis                       0.508244  0.0804192    1.3827   0.191475   1                 2920   0.551724  HALLMARK_ANGIOGENESIS
Largest Diameter            Hallmark Coagulation                        0.496778  0.0088       1.47505  0.100434   0.44              3178   0.494505  HALLMARK_COAGULATION
Variance Vox Val            Hallmark Coagulation                        0.450299  0.0469953    1.33796  0.215094   1                 3467   0.483516  HALLMARK_COAGULATION
Ld Init Enhancement Gt100   Hallmark Coagulation                        0.473078  0.0221       1.40602  0.141125   1                 3525   0.505495  HALLMARK_COAGULATION
Ld Late Lt0                 Hallmark Coagulation                        0.424813  0.0930907    1.2619   0.229609   1                 3418   0.461538  HALLMARK_COAGULATION
Largest Diameter            Hallmark Il2 Stat5 Signaling                0.400957  0.0425957    1.23923  0.196002   1                 4159   0.486667  HALLMARK_IL2_STAT5_SIGNALING
Ld Init Enhancement Gt100   Hallmark Il2 Stat5 Signaling                0.389525  0.0669933    1.20487  0.221814   1                 4100   0.473333  HALLMARK_IL2_STAT5_SIGNALING
Largest Diameter            Hallmark Peroxisome                         0.42527   0.039896     1.2411   0.196002   1                 3707   0.406977  HALLMARK_PEROXISOME
Ld Init Enhancement Gt100   Hallmark Peroxisome                         0.435172  0.0245975    1.26918  0.199533   1                 3701   0.418605  HALLMARK_PEROXISOME
Ld Late Lt0                 Hallmark Peroxisome                         0.426401  0.0375962    1.2432   0.229609   1                 2313   0.313953  HALLMARK_PEROXISOME
Irregularity                Hallmark Spermatogenesis                    0.571019  0.00319968   1.60095  0.102426   0.159984          1332   0.385965  HALLMARK_SPERMATOGENESIS
Volume                      Hallmark Spermatogenesis                    0.524219  0.0185981    1.47199  0.143555   0.929907           604   0.298246  HALLMARK_SPERMATOGENESIS
Largest Diameter            Hallmark Spermatogenesis                    0.549927  0.00819918   1.54873  0.100434   0.409959          2848   0.508772  HALLMARK_SPERMATOGENESIS
Mean Vox Val                Hallmark Spermatogenesis                    0.510354  0.030097     1.4326   0.249834   1                 2464   0.438596  HALLMARK_SPERMATOGENESIS
Variance Vox Val            Hallmark Spermatogenesis                    0.608164  0.0002       1.71178  0.0609472  0.01              1706   0.438596  HALLMARK_SPERMATOGENESIS
Top Late Enhancement        Hallmark Spermatogenesis                    0.54021   0.0112989    1.51811  0.247764   0.564944          2275   0.438596  HALLMARK_SPERMATOGENESIS
Vol Init Enhancement Gt100  Hallmark Spermatogenesis                    0.522971  0.019        1.46824  0.146559   0.95               630   0.298246  HALLMARK_SPERMATOGENESIS
Ld Init Enhancement Gt100   Hallmark Spermatogenesis                    0.55011   0.0071       1.54793  0.104422   0.355             2571   0.473684  HALLMARK_SPERMATOGENESIS
Ld Late Lt0                 Hallmark Spermatogenesis                    0.529765  0.0144986    1.48961  0.109763   0.724928          1790   0.385965  HALLMARK_SPERMATOGENESIS
Largest Diameter            Hallmark Kras Signaling Up                  0.407948  0.0857914    1.24662  0.196002   1                 3780   0.454545  HALLMARK_KRAS_SIGNALING_UP
Variance Vox Val            Hallmark Kras Signaling Up                  0.444082  0.030197     1.35367  0.212404   1                 4025   0.552448  HALLMARK_KRAS_SIGNALING_UP
Ld Init Enhancement Gt100   Hallmark Kras Signaling Up                  0.404424  0.0967903    1.23585  0.206088   1                 3643   0.440559  HALLMARK_KRAS_SIGNALING_UP
Largest Diameter            Hallmark Pancreas Beta Cells                0.500599  0.191638     1.21086  0.21978    1                 4434   0.647059  HALLMARK_PANCREAS_BETA_CELLS
Ld Init Enhancement Gt100   Hallmark Pancreas Beta Cells                0.489316  0.218753     1.18453  0.248889   1                 4982   0.705882  HALLMARK_PANCREAS_BETA_CELLS
Ld Late Lt0                 Hallmark Pancreas Beta Cells                0.509492  0.16483      1.22965  0.234378   1                 3968   0.588235  HALLMARK_PANCREAS_BETA_CELLS

</div>
