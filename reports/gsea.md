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



A utility function to display tables with [DataTables](https://datatables.net).


```python
def display_table(t):
    display(Markdown(
        '<div class="datatable">' +
        tabulate(t, headers='keys') +
        "\n\n</div>"
    ))
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
df_h = load_gsea_ds("../analyses/gsea/mri-features_h.all_T.nc")
df_h_rv = load_gsea_ds("../analyses/gsea/mri-features-reg-volume_h.all_T.nc")
df_cgp_rv = load_gsea_ds("../analyses/gsea/"
    "mri-features-reg-volume_c2.cgp_F.nc")
df_cp_rv = load_gsea_ds("../analyses/gsea/"
    "mri-features-reg-volume_c2.cp_T.nc")
```



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
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
```



A function to show significant pathway enrichments in a table.


```python
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

## Results---Volume regressed out ##

Hallmarks from MSigDB.


```python
plot_ds(df_h_rv, fdr=0.25)
```

![](figures/gsea_hallmarks-es-plot-rv_1){#hallmarks-es-plot-rv }\



```python
table_ds(df_h_rv, fdr=0.25)
```


<div class="datatable">mri_feature                     gene_set                                          es           p      nes       fdr       fwer    max_es_at    le_prop  gene_set_code
------------------------------  ------------------------------------------  --------  ----------  -------  --------  ---------  -----------  ---------  ------------------------------------------
Mean Smoothness All Timeframes  Hallmark G2M Checkpoint                     0.423469  0.155684    1.3269   0.241149  1                 4137   0.543956  HALLMARK_G2M_CHECKPOINT
Ld Late Lt0                     Hallmark Notch Signaling                    0.510357  0.0393495   1.39116  0.237305  1                 2548   0.392857  HALLMARK_NOTCH_SIGNALING
Mean Smoothness All Timeframes  Hallmark Adipogenesis                       0.431181  0.00349965  1.35268  0.237582  0.174983          1548   0.264865  HALLMARK_ADIPOGENESIS
Ld Late Lt0                     Hallmark Myogenesis                         0.509159  0.00259974  1.54987  0.207988  0.129987          2418   0.409449  HALLMARK_MYOGENESIS
Mean Smoothness All Timeframes  Hallmark Myogenesis                         0.444774  0.0234977   1.35872  0.237582  1                 2855   0.417323  HALLMARK_MYOGENESIS
Mean Smoothness All Timeframes  Hallmark E2F Targets                        0.439658  0.1456      1.38613  0.237582  1                 4405   0.595628  HALLMARK_E2F_TARGETS
Mean Smoothness All Timeframes  Hallmark Myc Targets V1                     0.42645   0.0661934   1.3709   0.237582  1                 4765   0.634021  HALLMARK_MYC_TARGETS_V1
Mean Smoothness All Timeframes  Hallmark Myc Targets V2                     0.539792  0.0274384   1.58631  0.237582  1                 3861   0.727273  HALLMARK_MYC_TARGETS_V2
Ld Late Lt0                     Hallmark Epithelial Mesenchymal Transition  0.492551  0.0820918   1.55582  0.207988  1                 3842   0.564246  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
Mean Smoothness All Timeframes  Hallmark Epithelial Mesenchymal Transition  0.531395  0.0424958   1.68018  0.237582  1                 2811   0.558659  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
Mean Smoothness All Timeframes  Hallmark Oxidative Phosphorylation          0.42228   0.0505949   1.35322  0.237582  1                 3487   0.473958  HALLMARK_OXIDATIVE_PHOSPHORYLATION
Ld Late Lt0                     Hallmark Uv Response Dn                     0.443049  0.0105989   1.3685   0.24173   0.529947          3452   0.462687  HALLMARK_UV_RESPONSE_DN
Mean Smoothness All Timeframes  Hallmark Uv Response Dn                     0.489173  0.00069993  1.51279  0.237582  0.0349965         2890   0.485075  HALLMARK_UV_RESPONSE_DN
Ld Late Lt0                     Hallmark Angiogenesis                       0.645808  0.00332527  1.76558  0.184437  0.166264          1836   0.448276  HALLMARK_ANGIOGENESIS
Mean Smoothness All Timeframes  Hallmark Angiogenesis                       0.532335  0.0484262   1.4483   0.237582  1                 2804   0.517241  HALLMARK_ANGIOGENESIS
Ld Late Lt0                     Hallmark Coagulation                        0.474589  0.0222978   1.40746  0.237305  1                 3324   0.494505  HALLMARK_COAGULATION
Ld Late Lt0                     Hallmark Bile Acid Metabolism               0.497515  0.00249975  1.43503  0.237305  0.124988          2947   0.457143  HALLMARK_BILE_ACID_METABOLISM
Ld Late Lt0                     Hallmark Spermatogenesis                    0.490717  0.0306      1.41177  0.237305  1                 1759   0.263158  HALLMARK_SPERMATOGENESIS
Mean Smoothness All Timeframes  Hallmark Spermatogenesis                    0.459858  0.0772923   1.32161  0.241149  1                 4434   0.596491  HALLMARK_SPERMATOGENESIS

</div>

Gene signatures (c2.cgp) from MSigDB.


```python
plot_ds(df_cgp_rv, fdr=0.25, le_prop=0.0, abs=False)
```

![](figures/gsea_cgp-es-plot-rv-a_1){#cgp-es-plot-rv-a }\



```python
plot_ds(df_cgp_rv, fdr=0.25, le_prop=0.8, abs=False)
```

![](figures/gsea_cgp-es-plot-rv-b_1){#cgp-es-plot-rv-b }\



```python
table_ds(df_cgp_rv, fdr=0.25)
```


<div class="datatable">mri_feature                          gene_set                                                           es            p       nes        fdr      fwer    max_es_at    le_prop  gene_set_code
-----------------------------------  ----------------------------------------------------------  ---------  -----------  --------  ---------  --------  -----------  ---------  ----------------------------------------------------------
Top Late Enhancement                 Nakamura Cancer Microenvironment Dn                         -0.701482  0.00603015   -1.78569  0.213855   1               10021   0.560976  NAKAMURA_CANCER_MICROENVIRONMENT_DN
Mean Smoothness All Timeframes       West Adrenocortical Tumor Markers Dn                        -0.723733  0.0513644    -1.58421  0.216167   1               10536   0.533333  WEST_ADRENOCORTICAL_TUMOR_MARKERS_DN
Top Late Enhancement                 Winter Hypoxia Up                                           -0.518715  0.0312065    -1.67149  0.23811    1                8803   0.512195  WINTER_HYPOXIA_UP
Washout                              Nakamura Tumor Zone Peripheral Vs Central Up                 0.508932  0.00204165    1.89197  0.193047   1                2538   0.53913   NAKAMURA_TUMOR_ZONE_PERIPHERAL_VS_CENTRAL_UP
Top Late Enhancement                 Nakamura Tumor Zone Peripheral Vs Central Up                -0.504224  0.00338174   -1.89259  0.213855   1                8877   0.508696  NAKAMURA_TUMOR_ZONE_PERIPHERAL_VS_CENTRAL_UP
Mean Smoothness All Timeframes       Piccaluga Angioimmunoblastic Lymphoma Up                    -0.554146  0.0759646    -1.61015  0.202604   1                9289   0.52356   PICCALUGA_ANGIOIMMUNOBLASTIC_LYMPHOMA_UP
Variation Smoothness All Timeframes  Liu Sox4 Targets Up                                         -0.429264  0.0203633    -1.63473  0.236712   1                8255   0.408333  LIU_SOX4_TARGETS_UP
Mean Smoothness All Timeframes       Liu Sox4 Targets Up                                         -0.421462  0.0252708    -1.6019   0.204175   1                8344   0.425     LIU_SOX4_TARGETS_UP
Variation Smoothness All Timeframes  Bertucci Medullary Vs Ductal Breast Cancer Dn               -0.654727  0.00837989   -1.92066  0.168147   1                9313   0.577465  BERTUCCI_MEDULLARY_VS_DUCTAL_BREAST_CANCER_DN
Mean Smoothness All Timeframes       Bertucci Medullary Vs Ductal Breast Cancer Dn               -0.673118  0.00383219   -1.96362  0.109401   1                8871   0.676056  BERTUCCI_MEDULLARY_VS_DUCTAL_BREAST_CANCER_DN
Irregularity                         Davicioni Pax Foxo1 Signature In Arms Up                     0.562401  0.00280505    1.85001  0.218099   1                1203   0.355556  DAVICIONI_PAX_FOXO1_SIGNATURE_IN_ARMS_UP
Mean Smoothness All Timeframes       Davicioni Pax Foxo1 Signature In Arms Up                    -0.476906  0.0382382    -1.56644  0.21886    1                7896   0.533333  DAVICIONI_PAX_FOXO1_SIGNATURE_IN_ARMS_UP
Variation Smoothness All Timeframes  Fournier Acinar Development Early Up                        -0.591486  0.0166633    -1.6829   0.231284   1                9808   0.368421  FOURNIER_ACINAR_DEVELOPMENT_EARLY_UP
Mean Smoothness All Timeframes       Fournier Acinar Development Early Up                        -0.584787  0.0200723    -1.66158  0.181553   1                9710   0.368421  FOURNIER_ACINAR_DEVELOPMENT_EARLY_UP
Mean Smoothness All Timeframes       Watanabe Colon Cancer Msi Vs Mss Dn                         -0.459035  0.0327416    -1.54531  0.227613   1                8731   0.414634  WATANABE_COLON_CANCER_MSI_VS_MSS_DN
Washout                              Sotiriou Breast Cancer Grade 1 Vs 3 Up                       0.822161  0.0105242     1.74502  0.24889    1                1469   0.85      SOTIRIOU_BREAST_CANCER_GRADE_1_VS_3_UP
Top Late Enhancement                 Sotiriou Breast Cancer Grade 1 Vs 3 Up                      -0.858786  0.001417     -1.81813  0.213855   1                9993   0.883333  SOTIRIOU_BREAST_CANCER_GRADE_1_VS_3_UP
Mean Smoothness All Timeframes       Sotiriou Breast Cancer Grade 1 Vs 3 Dn                      -0.627779  0.0675248    -1.54765  0.227442   1                7559   0.886364  SOTIRIOU_BREAST_CANCER_GRADE_1_VS_3_DN
Top Late Enhancement                 Chemnitz Response To Prostaglandin E2 Up                    -0.597293  0.0308567    -1.75535  0.213855   1                9173   0.561983  CHEMNITZ_RESPONSE_TO_PROSTAGLANDIN_E2_UP
Variation Smoothness All Timeframes  Chemnitz Response To Prostaglandin E2 Dn                    -0.408299  0.00798563   -1.66671  0.236189   1                7250   0.521429  CHEMNITZ_RESPONSE_TO_PROSTAGLANDIN_E2_DN
Mean Smoothness All Timeframes       Chemnitz Response To Prostaglandin E2 Dn                    -0.408617  0.00697767   -1.66476  0.181553   1                7374   0.510714  CHEMNITZ_RESPONSE_TO_PROSTAGLANDIN_E2_DN
Mean Smoothness All Timeframes       Igarashi Atf4 Targets Dn                                    -0.389166  0.0219386    -1.54457  0.228069   1                8914   0.358025  IGARASHI_ATF4_TARGETS_DN
Variation Smoothness All Timeframes  Zhong Response To Azacitidine And Tsa Up                    -0.402931  0.0153815    -1.62113  0.244212   1                8108   0.44186   ZHONG_RESPONSE_TO_AZACITIDINE_AND_TSA_UP
Mean Smoothness All Timeframes       Zhong Response To Azacitidine And Tsa Up                    -0.375708  0.0380024    -1.5095   0.249559   1                7580   0.48062   ZHONG_RESPONSE_TO_AZACITIDINE_AND_TSA_UP
Mean Smoothness All Timeframes       Davicioni Molecular Arms Vs Erms Dn                         -0.441355  0.0647295    -1.55237  0.226989   1                9670   0.302013  DAVICIONI_MOLECULAR_ARMS_VS_ERMS_DN
Mean Smoothness All Timeframes       Davicioni Targets Of Pax Foxo1 Fusions Up                   -0.448461  0.0531382    -1.58055  0.216167   1                8252   0.479452  DAVICIONI_TARGETS_OF_PAX_FOXO1_FUSIONS_UP
Mean Smoothness All Timeframes       Davicioni Rhabdomyosarcoma Pax Foxo1 Fusion Up              -0.440917  0.0321792    -1.53842  0.231366   1                8175   0.54      DAVICIONI_RHABDOMYOSARCOMA_PAX_FOXO1_FUSION_UP
Irregularity                         Sengupta Nasopharyngeal Carcinoma With Lmp1 Dn               0.595662  0.000615637   2.01824  0.12903    1                1655   0.423529  SENGUPTA_NASOPHARYNGEAL_CARCINOMA_WITH_LMP1_DN
Variation Smoothness All Timeframes  Sengupta Nasopharyngeal Carcinoma With Lmp1 Dn              -0.489487  0.0116653    -1.66048  0.236572   1                6882   0.623529  SENGUPTA_NASOPHARYNGEAL_CARCINOMA_WITH_LMP1_DN
Mean Smoothness All Timeframes       Sengupta Nasopharyngeal Carcinoma With Lmp1 Dn              -0.483036  0.0149676    -1.64156  0.194025   1                7836   0.505882  SENGUPTA_NASOPHARYNGEAL_CARCINOMA_WITH_LMP1_DN
Mean Smoothness All Timeframes       Turashvili Breast Normal Ductal Vs Lobular Up               -0.420202  0.0277108    -1.54962  0.226989   1                9010   0.360656  TURASHVILI_BREAST_NORMAL_DUCTAL_VS_LOBULAR_UP
Irregularity                         Turashvili Breast Ductal Carcinoma Vs Ductal Normal Dn       0.597016  0.0113246     1.83014  0.218099   1                2498   0.522013  TURASHVILI_BREAST_DUCTAL_CARCINOMA_VS_DUCTAL_NORMAL_DN
Variation Smoothness All Timeframes  Turashvili Breast Ductal Carcinoma Vs Ductal Normal Dn      -0.542701  0.0364949    -1.67509  0.233816   1                8029   0.566038  TURASHVILI_BREAST_DUCTAL_CARCINOMA_VS_DUCTAL_NORMAL_DN
Mean Smoothness All Timeframes       Turashvili Breast Ductal Carcinoma Vs Ductal Normal Dn      -0.578174  0.0181525    -1.78116  0.149917   1                7981   0.622642  TURASHVILI_BREAST_DUCTAL_CARCINOMA_VS_DUCTAL_NORMAL_DN
Irregularity                         Turashvili Breast Ductal Carcinoma Vs Lobular Normal Dn      0.671952  0.00462405    1.88478  0.191431   1                2242   0.607143  TURASHVILI_BREAST_DUCTAL_CARCINOMA_VS_LOBULAR_NORMAL_DN
Variation Smoothness All Timeframes  Turashvili Breast Ductal Carcinoma Vs Lobular Normal Dn     -0.598976  0.0276824    -1.68987  0.2224     1                8792   0.535714  TURASHVILI_BREAST_DUCTAL_CARCINOMA_VS_LOBULAR_NORMAL_DN
Mean Smoothness All Timeframes       Turashvili Breast Ductal Carcinoma Vs Lobular Normal Dn     -0.647827  0.00929105   -1.82827  0.133542   1                8298   0.678571  TURASHVILI_BREAST_DUCTAL_CARCINOMA_VS_LOBULAR_NORMAL_DN
Irregularity                         Turashvili Breast Lobular Carcinoma Vs Ductal Normal Dn      0.641105  0.00040032    2.02687  0.12903    1                2265   0.54878   TURASHVILI_BREAST_LOBULAR_CARCINOMA_VS_DUCTAL_NORMAL_DN
Variation Smoothness All Timeframes  Turashvili Breast Lobular Carcinoma Vs Ductal Normal Dn     -0.558535  0.0138554    -1.76026  0.193127   1                8313   0.573171  TURASHVILI_BREAST_LOBULAR_CARCINOMA_VS_DUCTAL_NORMAL_DN
Mean Smoothness All Timeframes       Turashvili Breast Lobular Carcinoma Vs Ductal Normal Dn     -0.58164   0.00602894   -1.83475  0.133542   1                8305   0.609756  TURASHVILI_BREAST_LOBULAR_CARCINOMA_VS_DUCTAL_NORMAL_DN
Irregularity                         Turashvili Breast Lobular Carcinoma Vs Lobular Normal Up     0.562542  0.00645422    1.83601  0.218099   1                1458   0.383721  TURASHVILI_BREAST_LOBULAR_CARCINOMA_VS_LOBULAR_NORMAL_UP
Variation Smoothness All Timeframes  Turashvili Breast Lobular Carcinoma Vs Lobular Normal Up    -0.512762  0.0256462    -1.66621  0.236189   1                8313   0.511628  TURASHVILI_BREAST_LOBULAR_CARCINOMA_VS_LOBULAR_NORMAL_UP
Mean Smoothness All Timeframes       Turashvili Breast Lobular Carcinoma Vs Lobular Normal Up    -0.551955  0.00770166   -1.79474  0.141179   1                8305   0.569767  TURASHVILI_BREAST_LOBULAR_CARCINOMA_VS_LOBULAR_NORMAL_UP
Mean Smoothness All Timeframes       Wilcox Response To Progesterone Dn                          -0.601481  0.0410246    -1.69561  0.169834   1                9507   0.477273  WILCOX_RESPONSE_TO_PROGESTERONE_DN
Mean Smoothness All Timeframes       Zhou Inflammatory Response Live Up                          -0.383455  0.0322838    -1.53202  0.233789   1                8133   0.407018  ZHOU_INFLAMMATORY_RESPONSE_LIVE_UP
Mean Smoothness All Timeframes       Hooi St7 Targets Up                                         -0.424361  0.0332736    -1.54699  0.227613   1                9098   0.384615  HOOI_ST7_TARGETS_UP
Top Late Enhancement                 Pramoonjago Sox4 Targets Dn                                 -0.561779  0.0132318    -1.74244  0.213855   1                8881   0.547619  PRAMOONJAGO_SOX4_TARGETS_DN
Mean Smoothness All Timeframes       Corre Multiple Myeloma Up                                   -0.546334  0.0264211    -1.66422  0.181553   1                9020   0.488372  CORRE_MULTIPLE_MYELOMA_UP
Mean Smoothness All Timeframes       Corre Multiple Myeloma Dn                                   -0.494     0.0408786    -1.60678  0.204175   1                9568   0.4375    CORRE_MULTIPLE_MYELOMA_DN
Variation Smoothness All Timeframes  Charafe Breast Cancer Basal Vs Mesenchymal Dn               -0.66496   0.0135377    -1.75461  0.193127   1                9110   0.605263  CHARAFE_BREAST_CANCER_BASAL_VS_MESENCHYMAL_DN
Mean Smoothness All Timeframes       Charafe Breast Cancer Basal Vs Mesenchymal Dn               -0.713759  0.00256765   -1.89045  0.117961   1                8602   0.763158  CHARAFE_BREAST_CANCER_BASAL_VS_MESENCHYMAL_DN
Washout                              Borczuk Malignant Mesothelioma Up                            0.509844  0.0146469     1.78793  0.221302   1                1794   0.418919  BORCZUK_MALIGNANT_MESOTHELIOMA_UP
Mean Smoothness All Timeframes       Roy Wound Blood Vessel Up                                   -0.57503   0.0902722    -1.54945  0.226989   1                9244   0.422222  ROY_WOUND_BLOOD_VESSEL_UP
Mean Smoothness All Timeframes       Newman Ercc6 Targets Dn                                     -0.642854  0.0643751    -1.57572  0.216178   1                9734   0.48      NEWMAN_ERCC6_TARGETS_DN
Mean Smoothness All Timeframes       Horiuchi Wtap Targets Up                                    -0.397874  0.0378083    -1.55008  0.226989   1                8629   0.374502  HORIUCHI_WTAP_TARGETS_UP
Top Late Enhancement                 Basaki Ybx1 Targets Up                                      -0.539472  0.0282678    -1.75236  0.213855   1                9932   0.428571  BASAKI_YBX1_TARGETS_UP
Variation Smoothness All Timeframes  Wikman Asbestos Lung Cancer Up                              -0.624841  0.0230738    -1.6314   0.238002   1                9330   0.533333  WIKMAN_ASBESTOS_LUNG_CANCER_UP
Mean Smoothness All Timeframes       Wikman Asbestos Lung Cancer Up                              -0.624235  0.0216303    -1.63275  0.194532   1                8563   0.666667  WIKMAN_ASBESTOS_LUNG_CANCER_UP
Mean Smoothness All Timeframes       Rodrigues Dcc Targets Dn                                    -0.413149  0.035801     -1.53793  0.231366   1                7358   0.571429  RODRIGUES_DCC_TARGETS_DN
Mean Smoothness All Timeframes       Rodrigues Ntn1 Targets Up                                   -0.66959   0.0252525    -1.61557  0.200977   1                7852   0.9       RODRIGUES_NTN1_TARGETS_UP
Circularity                          Vecchi Gastric Cancer Early Dn                              -0.595015  0.000605449  -1.97153  0.223353   1                8649   0.548023  VECCHI_GASTRIC_CANCER_EARLY_DN
Irregularity                         Vecchi Gastric Cancer Early Dn                               0.585274  0.00181855    1.94135  0.12903    1                2190   0.502825  VECCHI_GASTRIC_CANCER_EARLY_DN
Variation Smoothness All Timeframes  Vecchi Gastric Cancer Early Dn                              -0.613843  0.000402091  -2.03233  0.147564   1                9555   0.474576  VECCHI_GASTRIC_CANCER_EARLY_DN
Mean Smoothness All Timeframes       Vecchi Gastric Cancer Early Dn                              -0.605596  0.00040347   -2.00578  0.109401   1                8657   0.587571  VECCHI_GASTRIC_CANCER_EARLY_DN
Washout                              Jaeger Metastasis Up                                         0.613213  0.0115689     1.77537  0.235363   1                1777   0.542857  JAEGER_METASTASIS_UP
Top Late Enhancement                 Jaeger Metastasis Up                                        -0.617926  0.00958466   -1.79234  0.213855   1                9884   0.457143  JAEGER_METASTASIS_UP
Mean Smoothness All Timeframes       Jaeger Metastasis Dn                                        -0.448239  0.0208629    -1.63846  0.194025   1                7528   0.534884  JAEGER_METASTASIS_DN
Variation Smoothness All Timeframes  Ginestier Breast Cancer Znf217 Amplified Up                 -0.565826  0.0102968    -1.76669  0.193127   1                8329   0.56      GINESTIER_BREAST_CANCER_ZNF217_AMPLIFIED_UP
Mean Smoothness All Timeframes       Ginestier Breast Cancer Znf217 Amplified Up                 -0.54161   0.0205285    -1.68035  0.175441   1                8716   0.493333  GINESTIER_BREAST_CANCER_ZNF217_AMPLIFIED_UP
Mean Smoothness All Timeframes       Gargalovic Response To Oxidized Phospholipids Turquoise Up  -0.478866  0.0420052    -1.59106  0.211209   1                7922   0.5       GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_TURQUOISE_UP
Mean Smoothness All Timeframes       Gargalovic Response To Oxidized Phospholipids Red Dn        -0.683499  0.00959808   -1.74564  0.155059   1                9043   0.75      GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_RED_DN
Variation Smoothness All Timeframes  Gargalovic Response To Oxidized Phospholipids Blue Dn       -0.512593  0.00702388   -1.7408   0.198874   1                9786   0.326087  GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_BLUE_DN
Mean Smoothness All Timeframes       Gargalovic Response To Oxidized Phospholipids Blue Dn       -0.519783  0.00560448   -1.76735  0.153387   1                9318   0.369565  GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_BLUE_DN
Variation Smoothness All Timeframes  Gargalovic Response To Oxidized Phospholipids Grey Up       -0.763767  0.000404122  -1.92597  0.168147   1                9444   0.866667  GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_GREY_UP
Mean Smoothness All Timeframes       Gargalovic Response To Oxidized Phospholipids Grey Up       -0.703707  0.004671     -1.77428  0.149917   1                8846   0.866667  GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_GREY_UP
Variation Smoothness All Timeframes  Gargalovic Response To Oxidized Phospholipids Grey Dn       -0.490983  0.0103958    -1.7246   0.211541   1                8767   0.4       GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_GREY_DN
Mean Smoothness All Timeframes       Gargalovic Response To Oxidized Phospholipids Grey Dn       -0.469411  0.0169082    -1.64667  0.19228    1                9459   0.345455  GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_GREY_DN
Variation Smoothness All Timeframes  Bilban B Cll Lpl Up                                         -0.482322  0.00921659   -1.71135  0.213798   1                8058   0.433962  BILBAN_B_CLL_LPL_UP
Mean Smoothness All Timeframes       Bilban B Cll Lpl Up                                         -0.434938  0.0346       -1.5417   0.230896   1                7689   0.490566  BILBAN_B_CLL_LPL_UP
Variation Smoothness All Timeframes  Mullighan Mll Signature 1 Dn                                -0.391742  0.00528885   -1.67306  0.233816   1                8551   0.394444  MULLIGHAN_MLL_SIGNATURE_1_DN
Mean Smoothness All Timeframes       Mullighan Mll Signature 1 Dn                                -0.380515  0.00705219   -1.62921  0.194532   1                8632   0.372222  MULLIGHAN_MLL_SIGNATURE_1_DN
Variation Smoothness All Timeframes  Tonks Targets Of Runx1 Runx1T1 Fusion Hsc Up                -0.431306  0.0234798    -1.61961  0.244212   1                8695   0.411765  TONKS_TARGETS_OF_RUNX1_RUNX1T1_FUSION_HSC_UP
Mean Smoothness All Timeframes       Tonks Targets Of Runx1 Runx1T1 Fusion Hsc Up                -0.47578   0.00442834   -1.77882  0.149917   1                9148   0.418301  TONKS_TARGETS_OF_RUNX1_RUNX1T1_FUSION_HSC_UP
Irregularity                         Tonks Targets Of Runx1 Runx1T1 Fusion Granulocyte Up         0.53489   0.00497315    1.83165  0.218099   1                3161   0.590909  TONKS_TARGETS_OF_RUNX1_RUNX1T1_FUSION_GRANULOCYTE_UP
Variation Smoothness All Timeframes  Tonks Targets Of Runx1 Runx1T1 Fusion Granulocyte Up        -0.51674   0.00841852   -1.75704  0.193127   1               10281   0.318182  TONKS_TARGETS_OF_RUNX1_RUNX1T1_FUSION_GRANULOCYTE_UP
Mean Smoothness All Timeframes       Tonks Targets Of Runx1 Runx1T1 Fusion Granulocyte Up        -0.525415  0.00665054   -1.78377  0.149917   1                9822   0.363636  TONKS_TARGETS_OF_RUNX1_RUNX1T1_FUSION_GRANULOCYTE_UP
Top Late Enhancement                 Papaspyridonos Unstable Aterosclerotic Plaque Up            -0.594269  0.0478745    -1.67911  0.230757   1                9596   0.434783  PAPASPYRIDONOS_UNSTABLE_ATEROSCLEROTIC_PLAQUE_UP
Irregularity                         Papaspyridonos Unstable Aterosclerotic Plaque Dn             0.721079  0.00197316    1.94626  0.12903    1                1158   0.5       PAPASPYRIDONOS_UNSTABLE_ATEROSCLEROTIC_PLAQUE_DN
Variation Smoothness All Timeframes  Papaspyridonos Unstable Aterosclerotic Plaque Dn            -0.735615  0.000595002  -1.99076  0.168147   1               10205   0.578947  PAPASPYRIDONOS_UNSTABLE_ATEROSCLEROTIC_PLAQUE_DN
Mean Smoothness All Timeframes       Papaspyridonos Unstable Aterosclerotic Plaque Dn            -0.760839  0.00040016   -2.04648  0.109401   1               10011   0.631579  PAPASPYRIDONOS_UNSTABLE_ATEROSCLEROTIC_PLAQUE_DN
Variation Smoothness All Timeframes  Senese Hdac1 Targets Dn                                     -0.480555  0.00663984   -1.79181  0.193127   1                7979   0.533679  SENESE_HDAC1_TARGETS_DN
Mean Smoothness All Timeframes       Senese Hdac1 Targets Dn                                     -0.518965  0.00180868   -1.94044  0.109401   1                8310   0.53886   SENESE_HDAC1_TARGETS_DN
Variation Smoothness All Timeframes  Senese Hdac1 And Hdac2 Targets Dn                           -0.50644   0.0115538    -1.80763  0.192967   1                9059   0.447059  SENESE_HDAC1_AND_HDAC2_TARGETS_DN
Mean Smoothness All Timeframes       Senese Hdac1 And Hdac2 Targets Dn                           -0.532365  0.00521879   -1.89159  0.117961   1                9187   0.470588  SENESE_HDAC1_AND_HDAC2_TARGETS_DN
Variation Smoothness All Timeframes  Senese Hdac2 Targets Dn                                     -0.494201  0.0603258    -1.60286  0.247184   1                8542   0.52381   SENESE_HDAC2_TARGETS_DN
Mean Smoothness All Timeframes       Senese Hdac2 Targets Dn                                     -0.53913   0.022334     -1.74837  0.153387   1                9234   0.5       SENESE_HDAC2_TARGETS_DN
Variation Smoothness All Timeframes  Lee Neural Crest Stem Cell Up                               -0.503758  0.0582192    -1.59974  0.249895   1                8328   0.505263  LEE_NEURAL_CREST_STEM_CELL_UP
Mean Smoothness All Timeframes       Lee Neural Crest Stem Cell Up                               -0.491289  0.0715434    -1.56228  0.221717   1                9102   0.4       LEE_NEURAL_CREST_STEM_CELL_UP
Variation Smoothness All Timeframes  Lee Neural Crest Stem Cell Dn                               -0.475429  0.00862415   -1.70634  0.216651   1                9356   0.320513  LEE_NEURAL_CREST_STEM_CELL_DN
Mean Smoothness All Timeframes       Lee Neural Crest Stem Cell Dn                               -0.457878  0.0142342    -1.64123  0.194025   1                9828   0.282051  LEE_NEURAL_CREST_STEM_CELL_DN
Variation Smoothness Uptake          Tien Intestine Probiotics 6Hr Up                             0.656584  0.0149078     1.85034  0.173515   1                1995   0.622642  TIEN_INTESTINE_PROBIOTICS_6HR_UP
Variation Smoothness All Timeframes  Kinsey Targets Of Ewsr1 Flii Fusion Dn                      -0.51008   0.00841178   -1.82079  0.191445   1                8254   0.521739  KINSEY_TARGETS_OF_EWSR1_FLII_FUSION_DN
Mean Smoothness All Timeframes       Kinsey Targets Of Ewsr1 Flii Fusion Dn                      -0.546606  0.00278773   -1.95276  0.109401   1                8449   0.539855  KINSEY_TARGETS_OF_EWSR1_FLII_FUSION_DN
Irregularity                         Sabates Colorectal Adenoma Dn                                0.62918   0.00079602    1.96427  0.12903    1                2388   0.561983  SABATES_COLORECTAL_ADENOMA_DN
Variation Smoothness All Timeframes  Sabates Colorectal Adenoma Dn                               -0.704799  0.000203128  -2.19026  0.0298014  0.54093          9839   0.53719   SABATES_COLORECTAL_ADENOMA_DN
Mean Smoothness All Timeframes       Sabates Colorectal Adenoma Dn                               -0.679654  0.000203004  -2.11175  0.109401   0.540601         9901   0.528926  SABATES_COLORECTAL_ADENOMA_DN
Mean Smoothness All Timeframes       Nagashima Egf Signaling Up                                  -0.590907  0.115756     -1.51518  0.246789   1                8762   0.538462  NAGASHIMA_EGF_SIGNALING_UP
Mean Smoothness All Timeframes       Kim Wt1 Targets Up                                          -0.421335  0.0538045    -1.53759  0.231366   1                8762   0.384211  KIM_WT1_TARGETS_UP
Variation Smoothness All Timeframes  Kim Wt1 Targets 12Hr Up                                     -0.497131  0.0020024    -1.88586  0.168147   1                8690   0.446043  KIM_WT1_TARGETS_12HR_UP
Mean Smoothness All Timeframes       Kim Wt1 Targets 12Hr Up                                     -0.517535  0.000802729  -1.95859  0.109401   1                8871   0.453237  KIM_WT1_TARGETS_12HR_UP
Mean Smoothness All Timeframes       Elvidge Hypoxia Up                                          -0.44579   0.0622382    -1.54571  0.227613   1                8673   0.394558  ELVIDGE_HYPOXIA_UP
Top Late Enhancement                 Elvidge Hypoxia Dn                                          -0.47402   0.0191958    -1.68567  0.230757   1                8567   0.536232  ELVIDGE_HYPOXIA_DN
Variation Smoothness All Timeframes  Jaatinen Hematopoietic Stem Cell Up                         -0.383096  0.0146134    -1.6087   0.244618   1                9004   0.34322   JAATINEN_HEMATOPOIETIC_STEM_CELL_UP
Mean Smoothness All Timeframes       Jaatinen Hematopoietic Stem Cell Up                         -0.410935  0.00405515   -1.72434  0.16069    1                8533   0.423729  JAATINEN_HEMATOPOIETIC_STEM_CELL_UP
Washout                              Graham Cml Quiescent Vs Normal Quiescent Up                  0.580678  0.0141815     1.81583  0.206704   1                1439   0.476923  GRAHAM_CML_QUIESCENT_VS_NORMAL_QUIESCENT_UP
Top Late Enhancement                 Graham Cml Quiescent Vs Normal Quiescent Up                 -0.559536  0.0246266    -1.75323  0.213855   1                9993   0.446154  GRAHAM_CML_QUIESCENT_VS_NORMAL_QUIESCENT_UP
Top Late Enhancement                 Graham Cml Dividing Vs Normal Quiescent Up                  -0.623486  0.042946     -1.74807  0.213855   1                9993   0.554745  GRAHAM_CML_DIVIDING_VS_NORMAL_QUIESCENT_UP
Top Late Enhancement                 Graham Normal Quiescent Vs Normal Dividing Dn               -0.748067  0.0374721    -1.68507  0.230757   1                9993   0.716418  GRAHAM_NORMAL_QUIESCENT_VS_NORMAL_DIVIDING_DN
Variation Smoothness All Timeframes  Wamunyokoli Ovarian Cancer Lmp Dn                           -0.560738  0.00280786   -1.91634  0.168147   1                9094   0.458101  WAMUNYOKOLI_OVARIAN_CANCER_LMP_DN
Mean Smoothness All Timeframes       Wamunyokoli Ovarian Cancer Lmp Dn                           -0.576903  0.00141729   -1.96034  0.109401   1                8934   0.502793  WAMUNYOKOLI_OVARIAN_CANCER_LMP_DN
Variation Smoothness All Timeframes  Wamunyokoli Ovarian Cancer Grades 1 2 Dn                    -0.641488  0.00159681   -1.93656  0.168147   1                9302   0.568966  WAMUNYOKOLI_OVARIAN_CANCER_GRADES_1_2_DN
Mean Smoothness All Timeframes       Wamunyokoli Ovarian Cancer Grades 1 2 Dn                    -0.615175  0.0045954    -1.85576  0.129188   1                9148   0.534483  WAMUNYOKOLI_OVARIAN_CANCER_GRADES_1_2_DN
Mean Smoothness All Timeframes       Watanabe Ulcerative Colitis With Cancer Up                  -0.647826  0.0294294    -1.64022  0.194025   1                8686   0.533333  WATANABE_ULCERATIVE_COLITIS_WITH_CANCER_UP
Variation Smoothness All Timeframes  Rodrigues Thyroid Carcinoma Dn                              -0.473239  0.0256102    -1.61705  0.244212   1                8586   0.422535  RODRIGUES_THYROID_CARCINOMA_DN
Irregularity                         Gozgit Esr1 Targets Up                                       0.523227  0.00522193    1.84751  0.218099   1                3532   0.588785  GOZGIT_ESR1_TARGETS_UP
Variation Smoothness All Timeframes  Gozgit Esr1 Targets Up                                      -0.497379  0.00881587   -1.75464  0.193127   1                8132   0.53271   GOZGIT_ESR1_TARGETS_UP
Mean Smoothness All Timeframes       Gozgit Esr1 Targets Up                                      -0.498674  0.008        -1.76071  0.153387   1                7967   0.588785  GOZGIT_ESR1_TARGETS_UP
Variation Smoothness All Timeframes  Landis Breast Cancer Progression Dn                         -0.5625    0.0283759    -1.71274  0.213798   1                8841   0.508197  LANDIS_BREAST_CANCER_PROGRESSION_DN
Mean Smoothness All Timeframes       Landis Breast Cancer Progression Dn                         -0.581947  0.0177858    -1.7749   0.149917   1                8870   0.540984  LANDIS_BREAST_CANCER_PROGRESSION_DN
Irregularity                         Delys Thyroid Cancer Dn                                      0.586114  0.00318916    1.93755  0.12903    1                2787   0.566038  DELYS_THYROID_CANCER_DN
Variation Smoothness All Timeframes  Delys Thyroid Cancer Dn                                     -0.617946  0.000803536  -2.04259  0.147564   1                9605   0.509434  DELYS_THYROID_CANCER_DN
Mean Smoothness All Timeframes       Delys Thyroid Cancer Dn                                     -0.627896  0.000601443  -2.0784   0.109401   1                9711   0.522013  DELYS_THYROID_CANCER_DN
Mean Smoothness All Timeframes       Chiaradonna Neoplastic Transformation Kras Cdc25 Dn         -0.624207  0.0372168    -1.67838  0.175463   1                9321   0.536585  CHIARADONNA_NEOPLASTIC_TRANSFORMATION_KRAS_CDC25_DN
Washout                              Chiaradonna Neoplastic Transformation Kras Up                0.516626  0.00487112    1.80502  0.217568   1                2312   0.490909  CHIARADONNA_NEOPLASTIC_TRANSFORMATION_KRAS_UP
Top Late Enhancement                 Chiaradonna Neoplastic Transformation Kras Up               -0.518501  0.00443191   -1.81188  0.213855   1                8973   0.463636  CHIARADONNA_NEOPLASTIC_TRANSFORMATION_KRAS_UP
Variation Smoothness All Timeframes  Chiaradonna Neoplastic Transformation Kras Dn               -0.471288  0.0284284    -1.67444  0.233816   1                7847   0.528455  CHIARADONNA_NEOPLASTIC_TRANSFORMATION_KRAS_DN
Mean Smoothness All Timeframes       Chiaradonna Neoplastic Transformation Kras Dn               -0.48235   0.0204778    -1.70761  0.164501   1                8998   0.406504  CHIARADONNA_NEOPLASTIC_TRANSFORMATION_KRAS_DN
Variation Smoothness All Timeframes  Chebotaev Gr Targets Up                                     -0.516587  0.00564744   -1.77121  0.193127   1                9487   0.375     CHEBOTAEV_GR_TARGETS_UP
Mean Smoothness All Timeframes       Chebotaev Gr Targets Up                                     -0.551516  0.00120773   -1.8799   0.117961   1                8947   0.464286  CHEBOTAEV_GR_TARGETS_UP
Mean Smoothness All Timeframes       Chebotaev Gr Targets Dn                                     -0.586667  0.00863974   -1.8269   0.133542   1                9080   0.549296  CHEBOTAEV_GR_TARGETS_DN
Mean Smoothness All Timeframes       Berenjeno Rock Signaling Not Via Rhoa Dn                    -0.533592  0.0236568    -1.67409  0.17801    1                9040   0.489362  BERENJENO_ROCK_SIGNALING_NOT_VIA_RHOA_DN
Variation Smoothness All Timeframes  Berenjeno Transformed By Rhoa Reversibly Dn                 -0.59301   0.0475519    -1.60892  0.244618   1                8980   0.407407  BERENJENO_TRANSFORMED_BY_RHOA_REVERSIBLY_DN
Mean Smoothness All Timeframes       Berenjeno Transformed By Rhoa Reversibly Dn                 -0.647794  0.0136656    -1.75602  0.153387   1                8852   0.481481  BERENJENO_TRANSFORMED_BY_RHOA_REVERSIBLY_DN
Washout                              Creighton Akt1 Signaling Via Mtor Dn                         0.739675  0.00221864    1.89579  0.193047   1                2101   0.727273  CREIGHTON_AKT1_SIGNALING_VIA_MTOR_DN
Ser                                  Creighton Akt1 Signaling Via Mtor Dn                         0.770095  0.000608026   1.98214  0.21078    1                1949   0.727273  CREIGHTON_AKT1_SIGNALING_VIA_MTOR_DN
Top Late Enhancement                 Creighton Akt1 Signaling Via Mtor Dn                        -0.675719  0.0143376    -1.72785  0.217384   1                8903   0.636364  CREIGHTON_AKT1_SIGNALING_VIA_MTOR_DN
Vol Late Lt0                         Creighton Akt1 Signaling Via Mtor Dn                         0.773874  0.000598563   1.99325  0.206783   1                2060   0.818182  CREIGHTON_AKT1_SIGNALING_VIA_MTOR_DN
Variation Smoothness All Timeframes  Naderi Breast Cancer Prognosis Dn                           -0.659962  0.0316951    -1.64249  0.236712   1               10323   0.5       NADERI_BREAST_CANCER_PROGNOSIS_DN
Mean Smoothness All Timeframes       Naderi Breast Cancer Prognosis Dn                           -0.668469  0.0295537    -1.65999  0.181553   1                9875   0.5625    NADERI_BREAST_CANCER_PROGNOSIS_DN
Variation Smoothness All Timeframes  Vanharanta Uterine Fibroid Dn                               -0.520665  0.0407258    -1.63449  0.236712   1                7046   0.7       VANHARANTA_UTERINE_FIBROID_DN
Mean Smoothness All Timeframes       Vanharanta Uterine Fibroid Dn                               -0.526127  0.0380567    -1.65257  0.187178   1               10099   0.3       VANHARANTA_UTERINE_FIBROID_DN
Irregularity                         Landis Erbb2 Breast Tumors 65 Dn                             0.589284  0.00588928    1.83105  0.218099   1                1541   0.424242  LANDIS_ERBB2_BREAST_TUMORS_65_DN
Variation Smoothness All Timeframes  Landis Erbb2 Breast Tumors 65 Dn                            -0.577157  0.00884244   -1.79099  0.193127   1               10179   0.333333  LANDIS_ERBB2_BREAST_TUMORS_65_DN
Mean Smoothness All Timeframes       Landis Erbb2 Breast Tumors 65 Dn                            -0.587656  0.00702952   -1.82416  0.133542   1                9800   0.393939  LANDIS_ERBB2_BREAST_TUMORS_65_DN
Mean Smoothness All Timeframes       Gaussmann Mll Af4 Fusion Targets A Up                       -0.382886  0.0181167    -1.55658  0.223946   1                8910   0.319728  GAUSSMANN_MLL_AF4_FUSION_TARGETS_A_UP
Irregularity                         Gaussmann Mll Af4 Fusion Targets D Up                        0.663348  0.000986777   1.97836  0.12903    1                1374   0.4       GAUSSMANN_MLL_AF4_FUSION_TARGETS_D_UP
Variation Smoothness All Timeframes  Gaussmann Mll Af4 Fusion Targets D Up                       -0.56559   0.00972994   -1.69491  0.219274   1                8578   0.5       GAUSSMANN_MLL_AF4_FUSION_TARGETS_D_UP
Mean Smoothness All Timeframes       Gaussmann Mll Af4 Fusion Targets D Up                       -0.550916  0.0134122    -1.6533   0.187023   1                8299   0.55      GAUSSMANN_MLL_AF4_FUSION_TARGETS_D_UP
Variation Smoothness All Timeframes  Gaussmann Mll Af4 Fusion Targets E Up                       -0.571708  0.00527276   -1.84656  0.171536   1                9349   0.445946  GAUSSMANN_MLL_AF4_FUSION_TARGETS_E_UP
Mean Smoothness All Timeframes       Gaussmann Mll Af4 Fusion Targets E Up                       -0.595476  0.00344688   -1.91994  0.117961   1                9391   0.459459  GAUSSMANN_MLL_AF4_FUSION_TARGETS_E_UP
Variation Smoothness All Timeframes  Gaussmann Mll Af4 Fusion Targets E Dn                       -0.622474  0.0208249    -1.64112  0.236712   1                9420   0.466667  GAUSSMANN_MLL_AF4_FUSION_TARGETS_E_DN
Variation Smoothness All Timeframes  Gaussmann Mll Af4 Fusion Targets F Up                       -0.468256  0.0361955    -1.63309  0.237952   1                8220   0.496689  GAUSSMANN_MLL_AF4_FUSION_TARGETS_F_UP
Mean Smoothness All Timeframes       Gaussmann Mll Af4 Fusion Targets F Up                       -0.459331  0.0408369    -1.60171  0.204175   1                8513   0.456954  GAUSSMANN_MLL_AF4_FUSION_TARGETS_F_UP
Irregularity                         Gaussmann Mll Af4 Fusion Targets G Up                        0.46005   0.000402495   1.8619   0.218099   1                1392   0.239437  GAUSSMANN_MLL_AF4_FUSION_TARGETS_G_UP
Variation Smoothness All Timeframes  Gaussmann Mll Af4 Fusion Targets G Up                       -0.452323  0.000201491  -1.82722  0.188929   0.536571         8008   0.464789  GAUSSMANN_MLL_AF4_FUSION_TARGETS_G_UP
Mean Smoothness All Timeframes       Gaussmann Mll Af4 Fusion Targets G Up                       -0.473665  0.000202224  -1.90971  0.117961   0.538524         8842   0.394366  GAUSSMANN_MLL_AF4_FUSION_TARGETS_G_UP
Mean Smoothness All Timeframes       Berenjeno Transformed By Rhoa Forever Dn                    -0.515535  0.0492326    -1.53726  0.231366   1               10113   0.36      BERENJENO_TRANSFORMED_BY_RHOA_FOREVER_DN
Washout                              Ouellet Ovarian Cancer Invasive Vs Lmp Up                    0.6321    0.00240964    1.88968  0.193047   1                2485   0.589286  OUELLET_OVARIAN_CANCER_INVASIVE_VS_LMP_UP
Top Late Enhancement                 Ouellet Ovarian Cancer Invasive Vs Lmp Up                   -0.644014  0.00221106   -1.91707  0.213855   1                9347   0.526786  OUELLET_OVARIAN_CANCER_INVASIVE_VS_LMP_UP
Vol Late Lt0                         Ouellet Ovarian Cancer Invasive Vs Lmp Up                    0.628519  0.00298686    1.88753  0.214618   1                2271   0.580357  OUELLET_OVARIAN_CANCER_INVASIVE_VS_LMP_UP
Variation Smoothness Uptake          Ouellet Ovarian Cancer Invasive Vs Lmp Up                    0.589075  0.0161831     1.74825  0.201345   1                2981   0.714286  OUELLET_OVARIAN_CANCER_INVASIVE_VS_LMP_UP
Ld Late Lt0                          Landis Erbb2 Breast Tumors 324 Dn                            0.531165  0.00240288    1.9709   0.180214   1                1908   0.419847  LANDIS_ERBB2_BREAST_TUMORS_324_DN
Variation Smoothness All Timeframes  Landis Erbb2 Breast Tumors 324 Dn                           -0.470432  0.0134269    -1.75055  0.193127   1                9098   0.374046  LANDIS_ERBB2_BREAST_TUMORS_324_DN
Mean Smoothness All Timeframes       Landis Erbb2 Breast Tumors 324 Dn                           -0.502118  0.00444444   -1.86282  0.129188   1                8870   0.435115  LANDIS_ERBB2_BREAST_TUMORS_324_DN
Variation Smoothness Uptake          Barrier Cancer Relapse Tumor Sample Up                       0.752534  0.00339118    1.82981  0.173515   1                2315   0.666667  BARRIER_CANCER_RELAPSE_TUMOR_SAMPLE_UP
Ld Late Lt0                          Wang Esophagus Cancer Vs Normal Dn                           0.482032  0.0001998     2.04614  0.0973724  0.532068         2701   0.43956   WANG_ESOPHAGUS_CANCER_VS_NORMAL_DN
Variation Smoothness All Timeframes  Roylance Breast Cancer 16Q Copy Number Dn                   -0.785279  0.00915192   -1.7711   0.193127   1                9839   0.545455  ROYLANCE_BREAST_CANCER_16Q_COPY_NUMBER_DN
Mean Smoothness All Timeframes       Roylance Breast Cancer 16Q Copy Number Dn                   -0.753682  0.0219423    -1.7051   0.164501   1                9534   0.545455  ROYLANCE_BREAST_CANCER_16Q_COPY_NUMBER_DN
Top Late Enhancement                 Johansson Gliomagenesis By Pdgfb Up                         -0.538674  0.0116771    -1.76405  0.213855   1                9271   0.490909  JOHANSSON_GLIOMAGENESIS_BY_PDGFB_UP
Mean Smoothness All Timeframes       Mcbryan Pubertal Breast 3 4Wk Up                            -0.398575  0.0169966    -1.61156  0.202062   1                7955   0.455556  MCBRYAN_PUBERTAL_BREAST_3_4WK_UP
Irregularity                         Mcbryan Pubertal Breast 4 5Wk Dn                             0.429118  0.000590435   1.86267  0.218099   1                2832   0.377358  MCBRYAN_PUBERTAL_BREAST_4_5WK_DN
Ld Late Lt0                          Mcbryan Pubertal Breast 4 5Wk Dn                             0.50197   0.0002        2.17002  0.0483711  0.5326           1477   0.345912  MCBRYAN_PUBERTAL_BREAST_4_5WK_DN
Variation Smoothness All Timeframes  Mcbryan Pubertal Breast 4 5Wk Dn                            -0.41412   0.00263051   -1.78449  0.193127   1                9459   0.27044   MCBRYAN_PUBERTAL_BREAST_4_5WK_DN
Mean Smoothness All Timeframes       Mcbryan Pubertal Breast 4 5Wk Dn                            -0.403382  0.00303521   -1.73389  0.16069    1                9775   0.251572  MCBRYAN_PUBERTAL_BREAST_4_5WK_DN
Mean Smoothness All Timeframes       Mcbryan Pubertal Breast 5 6Wk Up                            -0.366424  0.012495     -1.55874  0.223275   1                7506   0.44      MCBRYAN_PUBERTAL_BREAST_5_6WK_UP
Variation Smoothness Uptake          Lui Thyroid Cancer Pax8 Pparg Dn                             0.574305  0.0163001     1.75935  0.192501   1                1864   0.454545  LUI_THYROID_CANCER_PAX8_PPARG_DN
Mean Smoothness All Timeframes       Ouellet Cultured Ovarian Cancer Invasive Vs Lmp Dn          -0.603016  0.00967615   -1.74245  0.155356   1                8865   0.47619   OUELLET_CULTURED_OVARIAN_CANCER_INVASIVE_VS_LMP_DN
Mean Smoothness All Timeframes       Mcbryan Pubertal Tgfb1 Targets Up                           -0.458537  0.0299919    -1.64893  0.190729   1                8689   0.429448  MCBRYAN_PUBERTAL_TGFB1_TARGETS_UP
Mean Smoothness All Timeframes       Mcbryan Pubertal Tgfb1 Targets Dn                           -0.415233  0.0334543    -1.51643  0.246789   1                9822   0.277778  MCBRYAN_PUBERTAL_TGFB1_TARGETS_DN
Mean Smoothness All Timeframes       Zirn Tretinoin Response Up                                  -0.604125  0.038392     -1.61229  0.202062   1                9053   0.444444  ZIRN_TRETINOIN_RESPONSE_UP
Variation Smoothness All Timeframes  Mcbryan Terminal End Bud Up                                 -0.647077  0.0324074    -1.61277  0.244212   1                9808   0.545455  MCBRYAN_TERMINAL_END_BUD_UP
Mean Smoothness All Timeframes       Mcbryan Terminal End Bud Up                                 -0.641251  0.0377662    -1.59883  0.205744   1                9893   0.545455  MCBRYAN_TERMINAL_END_BUD_UP
Mean Smoothness All Timeframes       Mahadevan Imatinib Resistance Up                            -0.728666  0.0158348    -1.72186  0.16069    1                9604   0.6       MAHADEVAN_IMATINIB_RESISTANCE_UP
Variation Smoothness All Timeframes  Mahadevan Imatinib Resistance Dn                            -0.685082  0.0224742    -1.66037  0.236572   1                9605   0.5       MAHADEVAN_IMATINIB_RESISTANCE_DN
Mean Smoothness All Timeframes       Mahadevan Imatinib Resistance Dn                            -0.678966  0.0238386    -1.64548  0.19228    1                9989   0.5       MAHADEVAN_IMATINIB_RESISTANCE_DN
Top Late Enhancement                 Tomida Metastasis Up                                        -0.656148  0.0189583    -1.70251  0.226338   1                8773   0.65      TOMIDA_METASTASIS_UP
Vol Late Lt0                         Tomida Metastasis Up                                         0.707242  0.00218862    1.85119  0.214618   1                2434   0.75      TOMIDA_METASTASIS_UP
Washout                              Hu Angiogenesis Dn                                           0.624237  0.00606061    1.84471  0.193047   1                1308   0.428571  HU_ANGIOGENESIS_DN
Top Late Enhancement                 Hu Angiogenesis Dn                                          -0.585199  0.0176849    -1.72124  0.219645   1                9253   0.457143  HU_ANGIOGENESIS_DN
Variation Smoothness All Timeframes  Begum Targets Of Pax3 Foxo1 Fusion Up                       -0.56802   0.0240696    -1.69279  0.220105   1                7682   0.702128  BEGUM_TARGETS_OF_PAX3_FOXO1_FUSION_UP
Mean Smoothness All Timeframes       Begum Targets Of Pax3 Foxo1 Fusion Up                       -0.588121  0.0127042    -1.74923  0.153387   1                9718   0.446809  BEGUM_TARGETS_OF_PAX3_FOXO1_FUSION_UP
Mean Smoothness All Timeframes       Begum Targets Of Pax3 Foxo1 Fusion Dn                       -0.549712  0.0633039    -1.57586  0.216178   1                7755   0.631579  BEGUM_TARGETS_OF_PAX3_FOXO1_FUSION_DN
Variation Smoothness All Timeframes  Wong Endmetrium Cancer Dn                                   -0.759225  0.0116349    -1.7969   0.193127   1                9174   0.803922  WONG_ENDMETRIUM_CANCER_DN
Mean Smoothness All Timeframes       Wong Endmetrium Cancer Dn                                   -0.82483   0.00059988   -1.95799  0.109401   1               10342   0.745098  WONG_ENDMETRIUM_CANCER_DN
Variation Smoothness All Timeframes  Chassot Skin Wound                                          -0.857296  0.0194098    -1.65016  0.236572   1               10412   0.7       CHASSOT_SKIN_WOUND
Mean Smoothness All Timeframes       Chassot Skin Wound                                          -0.808934  0.0481518    -1.5573   0.223946   1               10153   0.7       CHASSOT_SKIN_WOUND
Variation Smoothness All Timeframes  Riz Erythroid Differentiation Apobec2                       -0.63833   0.0143059    -1.67174  0.234704   1               10382   0.357143  RIZ_ERYTHROID_DIFFERENTIATION_APOBEC2
Mean Smoothness All Timeframes       Riz Erythroid Differentiation Apobec2                       -0.596393  0.0380644    -1.56889  0.21886    1                9871   0.357143  RIZ_ERYTHROID_DIFFERENTIATION_APOBEC2
Irregularity                         Riz Erythroid Differentiation 6Hr                            0.670411  0.000606428   1.94605  0.12903    1                 713   0.35      RIZ_ERYTHROID_DIFFERENTIATION_6HR
Variation Smoothness All Timeframes  Riz Erythroid Differentiation 12Hr                          -0.59381   0.00913423   -1.75834  0.193127   1                9631   0.419355  RIZ_ERYTHROID_DIFFERENTIATION_12HR
Mean Smoothness All Timeframes       Riz Erythroid Differentiation 12Hr                          -0.617443  0.00457985   -1.82463  0.133542   1                9173   0.483871  RIZ_ERYTHROID_DIFFERENTIATION_12HR
Mean Smoothness All Timeframes       Ebauer Targets Of Pax3 Foxo1 Fusion Up                      -0.416656  0.0197015    -1.61595  0.200977   1                9555   0.333333  EBAUER_TARGETS_OF_PAX3_FOXO1_FUSION_UP
Variation Smoothness All Timeframes  Perez Tp63 Targets                                          -0.452975  0.00279888   -1.81258  0.191445   1                8063   0.46087   PEREZ_TP63_TARGETS
Mean Smoothness All Timeframes       Perez Tp63 Targets                                          -0.427855  0.00778443   -1.71386  0.160858   1                8238   0.430435  PEREZ_TP63_TARGETS
Variation Smoothness All Timeframes  Perez Tp53 And Tp63 Targets                                 -0.506264  0.00160224   -1.88603  0.168147   1                8642   0.451852  PEREZ_TP53_AND_TP63_TARGETS
Mean Smoothness All Timeframes       Perez Tp53 And Tp63 Targets                                 -0.462397  0.0122147    -1.72138  0.16069    1                8252   0.466667  PEREZ_TP53_AND_TP63_TARGETS
Variation Smoothness All Timeframes  Dacosta Uv Response Via Ercc3 Ttd Dn                        -0.548942  0.0125174    -1.73277  0.210109   1                8281   0.576923  DACOSTA_UV_RESPONSE_VIA_ERCC3_TTD_DN
Mean Smoothness All Timeframes       Dacosta Uv Response Via Ercc3 Ttd Dn                        -0.568545  0.00638595   -1.79514  0.141179   1                9249   0.461538  DACOSTA_UV_RESPONSE_VIA_ERCC3_TTD_DN
Mean Smoothness All Timeframes       Breuhahn Growth Factor Signaling In Liver Cancer            -0.550456  0.0685635    -1.51518  0.246789   1                8839   0.619048  BREUHAHN_GROWTH_FACTOR_SIGNALING_IN_LIVER_CANCER
Mean Smoothness All Timeframes       Roversi Glioma Loh Regions                                  -0.549462  0.00937627   -1.74244  0.155356   1                9980   0.346154  ROVERSI_GLIOMA_LOH_REGIONS
Washout                              Li Amplified In Lung Cancer                                  0.585873  0.00668964    1.88778  0.193047   1                2724   0.632258  LI_AMPLIFIED_IN_LUNG_CANCER
Top Late Enhancement                 Li Amplified In Lung Cancer                                 -0.59583   0.00323363   -1.90424  0.213855   1                8873   0.548387  LI_AMPLIFIED_IN_LUNG_CANCER
Vol Late Lt0                         Li Amplified In Lung Cancer                                  0.565137  0.010326      1.82833  0.240607   1                2849   0.664516  LI_AMPLIFIED_IN_LUNG_CANCER
Variation Smoothness Uptake          Li Amplified In Lung Cancer                                  0.579428  0.00558882    1.85931  0.173515   1                1822   0.470968  LI_AMPLIFIED_IN_LUNG_CANCER
Mean Smoothness All Timeframes       Ebauer Myogenic Targets Of Pax3 Foxo1 Fusion                -0.544105  0.0580581    -1.555    0.22495    1                8942   0.592593  EBAUER_MYOGENIC_TARGETS_OF_PAX3_FOXO1_FUSION
Top Late Enhancement                 Mattioli Mgus Vs Pcl                                        -0.53107   0.0231519    -1.68411  0.230757   1                9241   0.483871  MATTIOLI_MGUS_VS_PCL
Variation Smoothness Uptake          Mattioli Mgus Vs Pcl                                         0.555115  0.00807087    1.77253  0.192501   1                1560   0.483871  MATTIOLI_MGUS_VS_PCL
Mean Smoothness All Timeframes       Lui Thyroid Cancer Cluster 2                                -0.493627  0.0500603    -1.53522  0.232143   1                9326   0.361111  LUI_THYROID_CANCER_CLUSTER_2
Variation Smoothness All Timeframes  Lui Thyroid Cancer Cluster 3                                 0.732357  0.00238854    1.92411  0.206872   1                1752   0.62963   LUI_THYROID_CANCER_CLUSTER_3
Variation Smoothness Uptake          Lui Thyroid Cancer Cluster 3                                 0.735783  0.00277943    1.94006  0.173515   1                1788   0.62963   LUI_THYROID_CANCER_CLUSTER_3
Variation Smoothness All Timeframes  Lui Targets Of Pax8 Pparg Fusion                             0.645671  0.0130667     1.81509  0.237345   1                1752   0.545455  LUI_TARGETS_OF_PAX8_PPARG_FUSION
Variation Smoothness Uptake          Lui Targets Of Pax8 Pparg Fusion                             0.653388  0.00840504    1.83098  0.173515   1                1510   0.515152  LUI_TARGETS_OF_PAX8_PPARG_FUSION
Vol Late Lt0                         Schlosser Myc And Serum Response Synergy                     0.659812  0.00503221    1.84617  0.214618   1                1934   0.678571  SCHLOSSER_MYC_AND_SERUM_RESPONSE_SYNERGY
Variation Smoothness All Timeframes  Schlosser Myc And Serum Response Synergy                     0.674736  0.00519896    1.89124  0.212371   1                1166   0.5       SCHLOSSER_MYC_AND_SERUM_RESPONSE_SYNERGY
Variation Smoothness Uptake          Schlosser Myc And Serum Response Synergy                     0.676905  0.00400561    1.89198  0.173515   1                1867   0.571429  SCHLOSSER_MYC_AND_SERUM_RESPONSE_SYNERGY
Top Late Enhancement                 Rosty Cervical Cancer Proliferation Cluster                 -0.815014  0.0129686    -1.76673  0.213855   1                9973   0.771186  ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER
Irregularity                         Kobayashi Response To Romidepsin                             0.749369  0.0002002     2.01275  0.12903    0.533133         1584   0.466667  KOBAYASHI_RESPONSE_TO_ROMIDEPSIN
Mean Smoothness All Timeframes       Kobayashi Response To Romidepsin                            -0.563047  0.055589     -1.51092  0.248827   1                8215   0.533333  KOBAYASHI_RESPONSE_TO_ROMIDEPSIN
Variation Smoothness All Timeframes  Li Cisplatin Resistance Dn                                  -0.640735  0.0103917    -1.74782  0.194054   1                9871   0.478261  LI_CISPLATIN_RESISTANCE_DN
Mean Smoothness All Timeframes       Li Cisplatin Resistance Dn                                  -0.627217  0.0156752    -1.71242  0.160858   1                9772   0.478261  LI_CISPLATIN_RESISTANCE_DN
Mean Smoothness All Timeframes       Li Cisplatin Resistance Up                                  -0.590179  0.0389402    -1.59771  0.206656   1                8731   0.529412  LI_CISPLATIN_RESISTANCE_UP
Mean Smoothness All Timeframes       Hernandez Aberrant Mitosis By Docetacel 4Nm Up              -0.589095  0.0376485    -1.57021  0.218767   1                9149   0.5       HERNANDEZ_ABERRANT_MITOSIS_BY_DOCETACEL_4NM_UP
Mean Smoothness All Timeframes       Hernandez Mitotic Arrest By Docetaxel 1 Up                  -0.584139  0.0130864    -1.71704  0.16069    1                9634   0.535714  HERNANDEZ_MITOTIC_ARREST_BY_DOCETAXEL_1_UP
Mean Smoothness All Timeframes       Hernandez Aberrant Mitosis By Docetacel 2Nm Dn              -0.566933  0.0424671    -1.56709  0.21886    1                8985   0.5       HERNANDEZ_ABERRANT_MITOSIS_BY_DOCETACEL_2NM_DN
Irregularity                         Hernandez Mitotic Arrest By Docetaxel 2 Dn                   0.70032   0.0046748     1.81546  0.232085   1                1708   0.466667  HERNANDEZ_MITOTIC_ARREST_BY_DOCETAXEL_2_DN
Variation Smoothness All Timeframes  Hernandez Mitotic Arrest By Docetaxel 2 Dn                  -0.647206  0.0192659    -1.6744   0.233816   1               10174   0.4       HERNANDEZ_MITOTIC_ARREST_BY_DOCETAXEL_2_DN
Mean Smoothness All Timeframes       Hernandez Mitotic Arrest By Docetaxel 2 Dn                  -0.63092   0.027722     -1.63229  0.194532   1                9657   0.466667  HERNANDEZ_MITOTIC_ARREST_BY_DOCETAXEL_2_DN
Top Late Enhancement                 Xu Hgf Signaling Not Via Akt1 48Hr Dn                       -0.687273  0.00673401   -1.8359   0.213855   1               10442   0.388889  XU_HGF_SIGNALING_NOT_VIA_AKT1_48HR_DN
Mean Smoothness All Timeframes       Shi Sparc Targets Up                                        -0.549921  0.0255216    -1.6009   0.204175   1                8093   0.666667  SHI_SPARC_TARGETS_UP
Mean Smoothness All Timeframes       Snijders Amplified In Head And Neck Tumors                  -0.501978  0.0540705    -1.54536  0.227613   1                9014   0.40625   SNIJDERS_AMPLIFIED_IN_HEAD_AND_NECK_TUMORS
Mean Smoothness All Timeframes       Korkola Teratoma                                            -0.477496  0.0609562    -1.50879  0.249559   1               10443   0.354839  KORKOLA_TERATOMA
Variation Smoothness All Timeframes  Korkola Yolk Sac Tumor                                      -0.581067  0.0115713    -1.7179   0.213798   1                9265   0.46875   KORKOLA_YOLK_SAC_TUMOR
Mean Smoothness All Timeframes       Korkola Yolk Sac Tumor                                      -0.572636  0.0133065    -1.70592  0.164501   1                7997   0.625     KORKOLA_YOLK_SAC_TUMOR
Variation Smoothness All Timeframes  Teramoto Opn Targets Cluster 6                              -0.591506  0.00549339   -1.81604  0.191445   1                9705   0.375     TERAMOTO_OPN_TARGETS_CLUSTER_6
Mean Smoothness All Timeframes       Teramoto Opn Targets Cluster 6                              -0.541392  0.0211013    -1.66939  0.180021   1                9819   0.375     TERAMOTO_OPN_TARGETS_CLUSTER_6
Top Late Enhancement                 Meinhold Ovarian Cancer Low Grade Dn                        -0.685231  0.01642      -1.70186  0.226338   1               10209   0.55      MEINHOLD_OVARIAN_CANCER_LOW_GRADE_DN
Variation Smoothness Uptake          Meinhold Ovarian Cancer Low Grade Dn                         0.68816   0.013766      1.72725  0.234559   1                2829   0.75      MEINHOLD_OVARIAN_CANCER_LOW_GRADE_DN
Variation Smoothness All Timeframes  Dacosta Ercc3 Allele Xpcs Vs Ttd Dn                         -0.586691  0.0191764    -1.69642  0.219274   1                9200   0.551724  DACOSTA_ERCC3_ALLELE_XPCS_VS_TTD_DN
Mean Smoothness All Timeframes       Dacosta Ercc3 Allele Xpcs Vs Ttd Dn                         -0.552903  0.0417252    -1.60705  0.204175   1                8947   0.586207  DACOSTA_ERCC3_ALLELE_XPCS_VS_TTD_DN
Variation Smoothness All Timeframes  Yanagihara Esx1 Targets                                     -0.565227  0.0335451    -1.61257  0.244212   1                8751   0.44      YANAGIHARA_ESX1_TARGETS
Mean Smoothness All Timeframes       Yanagihara Esx1 Targets                                     -0.586344  0.0188268    -1.67522  0.177783   1                8134   0.56      YANAGIHARA_ESX1_TARGETS
Mean Smoothness All Timeframes       Tsunoda Cisplatin Resistance Up                             -0.634293  0.063332     -1.53671  0.231371   1                8991   0.642857  TSUNODA_CISPLATIN_RESISTANCE_UP
Variation Smoothness All Timeframes  Mohankumar Tlx1 Targets Dn                                  -0.476945  0.0294235    -1.6411   0.236712   1                8977   0.420168  MOHANKUMAR_TLX1_TARGETS_DN
Mean Smoothness All Timeframes       Mohankumar Tlx1 Targets Dn                                  -0.502363  0.013756     -1.72651  0.16069    1                9057   0.453782  MOHANKUMAR_TLX1_TARGETS_DN
Variation Smoothness All Timeframes  Dairkee Tert Targets Dn                                     -0.473069  0.00438334   -1.734    0.209758   1                6544   0.719512  DAIRKEE_TERT_TARGETS_DN
Mean Smoothness All Timeframes       Dairkee Tert Targets Dn                                     -0.470605  0.0047553    -1.7248   0.16069    1                6910   0.670732  DAIRKEE_TERT_TARGETS_DN
Variation Smoothness All Timeframes  Dawson Methylated In Lymphoma Tcl1                          -0.520177  0.0310263    -1.61454  0.244212   1                7557   0.677419  DAWSON_METHYLATED_IN_LYMPHOMA_TCL1
Mean Smoothness All Timeframes       Dawson Methylated In Lymphoma Tcl1                          -0.527485  0.0270862    -1.63844  0.194025   1                7832   0.645161  DAWSON_METHYLATED_IN_LYMPHOMA_TCL1
Mean Smoothness All Timeframes       Simbulan Uv Response Normal Dn                              -0.50217   0.0478691    -1.53727  0.231366   1                9556   0.387097  SIMBULAN_UV_RESPONSE_NORMAL_DN
Variation Smoothness All Timeframes  Simbulan Uv Response Immortalized Dn                        -0.545962  0.0241276    -1.64512  0.236712   1                9105   0.433333  SIMBULAN_UV_RESPONSE_IMMORTALIZED_DN
Mean Smoothness All Timeframes       Simbulan Uv Response Immortalized Dn                        -0.55703   0.0195453    -1.68162  0.175441   1                9556   0.433333  SIMBULAN_UV_RESPONSE_IMMORTALIZED_DN
Variation Smoothness All Timeframes  Wattel Autonomous Thyroid Adenoma Dn                        -0.606391  0.00865714   -1.84137  0.171536   1                8073   0.636364  WATTEL_AUTONOMOUS_THYROID_ADENOMA_DN
Mean Smoothness All Timeframes       Wattel Autonomous Thyroid Adenoma Dn                        -0.598608  0.00948537   -1.817    0.136277   1                8010   0.681818  WATTEL_AUTONOMOUS_THYROID_ADENOMA_DN
Variation Smoothness All Timeframes  Sasai Resistance To Neoplastic Transfromation               -0.603498  0.0706357    -1.61146  0.244212   1                8595   0.595745  SASAI_RESISTANCE_TO_NEOPLASTIC_TRANSFROMATION
Mean Smoothness All Timeframes       Sasai Resistance To Neoplastic Transfromation               -0.643636  0.0346844    -1.71809  0.16069    1                9098   0.617021  SASAI_RESISTANCE_TO_NEOPLASTIC_TRANSFROMATION
Mean Smoothness All Timeframes       Kim Mycn Amplification Targets Dn                           -0.420807  0.0183246    -1.6061   0.204175   1                8647   0.355556  KIM_MYCN_AMPLIFICATION_TARGETS_DN
Variation Smoothness All Timeframes  Kim Mycl1 Amplification Targets Dn                          -0.610606  0.0129817    -1.69389  0.219598   1                9439   0.4       KIM_MYCL1_AMPLIFICATION_TARGETS_DN
Mean Smoothness All Timeframes       Kim Mycl1 Amplification Targets Dn                          -0.561473  0.0385243    -1.55996  0.223217   1                9157   0.4       KIM_MYCL1_AMPLIFICATION_TARGETS_DN
Variation Smoothness All Timeframes  Hasina Nol7 Targets Up                                      -0.701889  0.0315369    -1.61545  0.244212   1               10306   0.454545  HASINA_NOL7_TARGETS_UP
Mean Smoothness All Timeframes       Hasina Nol7 Targets Up                                      -0.709244  0.0273344    -1.62451  0.196797   1                9947   0.454545  HASINA_NOL7_TARGETS_UP
Mean Smoothness All Timeframes       Hatada Methylated In Lung Cancer Up                         -0.379294  0.0275614    -1.55922  0.223217   1                8834   0.336683  HATADA_METHYLATED_IN_LUNG_CANCER_UP
Top Late Enhancement                 Furukawa Dusp6 Targets Pci35 Dn                             -0.65112   0.0364225    -1.73483  0.213855   1                9440   0.647059  FURUKAWA_DUSP6_TARGETS_PCI35_DN
Mean Smoothness All Timeframes       Lien Breast Carcinoma Metaplastic                           -0.685154  0.0783373    -1.5776   0.216178   1                8272   0.866667  LIEN_BREAST_CARCINOMA_METAPLASTIC
Irregularity                         Tomlins Prostate Cancer Dn                                   0.69298   0.0141491     1.80742  0.243479   1                2240   0.684211  TOMLINS_PROSTATE_CANCER_DN
Variation Smoothness All Timeframes  Tomlins Prostate Cancer Dn                                  -0.726592  0.00397456   -1.90567  0.168147   1                9560   0.657895  TOMLINS_PROSTATE_CANCER_DN
Mean Smoothness All Timeframes       Tomlins Prostate Cancer Dn                                  -0.715064  0.0056304    -1.85811  0.129188   1                9433   0.657895  TOMLINS_PROSTATE_CANCER_DN
Variation Smoothness All Timeframes  Theodorou Mammary Tumorigenesis                             -0.741619  0.00120048   -1.87487  0.168147   1                9268   0.666667  THEODOROU_MAMMARY_TUMORIGENESIS
Mean Smoothness All Timeframes       Theodorou Mammary Tumorigenesis                             -0.748384  0.000801603  -1.89219  0.117961   1                8774   0.722222  THEODOROU_MAMMARY_TUMORIGENESIS
Variation Smoothness All Timeframes  Schlesinger Methylated De Novo In Cancer                    -0.553377  0.0137807    -1.72477  0.211541   1                7891   0.613636  SCHLESINGER_METHYLATED_DE_NOVO_IN_CANCER
Mean Smoothness All Timeframes       Schlesinger Methylated De Novo In Cancer                    -0.555819  0.0134376    -1.72782  0.16069    1                7235   0.75      SCHLESINGER_METHYLATED_DE_NOVO_IN_CANCER
Variation Smoothness All Timeframes  Cowling Mycn Targets                                        -0.633754  0.0507686    -1.62658  0.242889   1                8812   0.653846  COWLING_MYCN_TARGETS
Mean Smoothness All Timeframes       Cowling Mycn Targets                                        -0.685238  0.0176292    -1.74891  0.153387   1                9457   0.615385  COWLING_MYCN_TARGETS
Variation Smoothness All Timeframes  Gross Elk3 Targets Up                                       -0.540289  0.0286918    -1.64265  0.236712   1                8874   0.576923  GROSS_ELK3_TARGETS_UP
Mean Smoothness All Timeframes       Gross Elk3 Targets Up                                       -0.521067  0.0405487    -1.57709  0.216178   1                9403   0.538462  GROSS_ELK3_TARGETS_UP
Irregularity                         Gross Hypoxia Via Elk3 Only Up                               0.659718  0.00569221    1.82493  0.218099   1                2139   0.52      GROSS_HYPOXIA_VIA_ELK3_ONLY_UP
Variation Smoothness All Timeframes  Ingram Shh Targets Up                                       -0.547609  0.00182927   -1.90281  0.168147   1                8331   0.543478  INGRAM_SHH_TARGETS_UP
Mean Smoothness All Timeframes       Ingram Shh Targets Up                                       -0.528734  0.00500701   -1.84008  0.133542   1                9232   0.434783  INGRAM_SHH_TARGETS_UP
Mean Smoothness All Timeframes       Ingram Shh Targets Dn                                       -0.467126  0.0403079    -1.56128  0.222247   1               10216   0.285714  INGRAM_SHH_TARGETS_DN
Mean Smoothness All Timeframes       Mantovani Viral Gpcr Signaling Up                           -0.435503  0.0432259    -1.51283  0.248519   1               10011   0.288136  MANTOVANI_VIRAL_GPCR_SIGNALING_UP
Washout                              Rahman Tp53 Targets Phosphorylated                           0.723344  0.00121433    1.83961  0.193047   1                1796   0.619048  RAHMAN_TP53_TARGETS_PHOSPHORYLATED
Top Late Enhancement                 Rahman Tp53 Targets Phosphorylated                          -0.710034  0.00502614   -1.80114  0.213855   1                8127   0.904762  RAHMAN_TP53_TARGETS_PHOSPHORYLATED
Vol Late Lt0                         Rahman Tp53 Targets Phosphorylated                           0.76624   0.000601564   1.95016  0.214618   1                1817   0.619048  RAHMAN_TP53_TARGETS_PHOSPHORYLATED
Variation Smoothness Uptake          Rahman Tp53 Targets Phosphorylated                           0.72861   0.0030036     1.85496  0.173515   1                1794   0.52381   RAHMAN_TP53_TARGETS_PHOSPHORYLATED
Mean Smoothness All Timeframes       Guenther Growth Spherical Vs Adherent Dn                    -0.618458  0.0195382    -1.71665  0.16069    1                9106   0.56      GUENTHER_GROWTH_SPHERICAL_VS_ADHERENT_DN
Mean Smoothness All Timeframes       Schwab Targets Of Bmyb Polymorphic Variants Dn              -0.594362  0.0520431    -1.51816  0.245417   1                8947   0.416667  SCHWAB_TARGETS_OF_BMYB_POLYMORPHIC_VARIANTS_DN
Top Late Enhancement                 Caffarel Response To Thc 24Hr 5 Up                          -0.571861  0.0159141    -1.69289  0.227222   1               10021   0.4375    CAFFAREL_RESPONSE_TO_THC_24HR_5_UP
Variation Smoothness Uptake          Caffarel Response To Thc 24Hr 5 Up                           0.614023  0.0059125     1.81637  0.174539   1                1767   0.53125   CAFFAREL_RESPONSE_TO_THC_24HR_5_UP
Variation Smoothness All Timeframes  Bertucci Invasive Carcinoma Ductal Vs Lobular Dn            -0.821446  0.000396197  -1.9747   0.168147   1               10328   0.771429  BERTUCCI_INVASIVE_CARCINOMA_DUCTAL_VS_LOBULAR_DN
Mean Smoothness All Timeframes       Bertucci Invasive Carcinoma Ductal Vs Lobular Dn            -0.774546  0.00620869   -1.85352  0.129779   1               10170   0.742857  BERTUCCI_INVASIVE_CARCINOMA_DUCTAL_VS_LOBULAR_DN
Variation Smoothness All Timeframes  Garcia Targets Of Fli1 And Dax1 Dn                           0.552469  0.00874404    1.79487  0.247899   1                2673   0.616     GARCIA_TARGETS_OF_FLI1_AND_DAX1_DN
Variation Smoothness Uptake          Garcia Targets Of Fli1 And Dax1 Dn                           0.540936  0.0109366     1.7612   0.192501   1                2232   0.544     GARCIA_TARGETS_OF_FLI1_AND_DAX1_DN
Mean Smoothness All Timeframes       Fridman Senescence Dn                                       -0.659118  0.0306306    -1.61944  0.200738   1                9822   0.545455  FRIDMAN_SENESCENCE_DN
Variation Smoothness All Timeframes  Schaeffer Prostate Development 6Hr Up                       -0.383021  0.00731856   -1.63479  0.236712   1                9057   0.302158  SCHAEFFER_PROSTATE_DEVELOPMENT_6HR_UP
Mean Smoothness All Timeframes       Schaeffer Prostate Development 6Hr Up                       -0.395666  0.00429097   -1.67952  0.175441   1                9043   0.309353  SCHAEFFER_PROSTATE_DEVELOPMENT_6HR_UP
Mean Smoothness All Timeframes       Schaeffer Prostate Development 12Hr Up                      -0.453677  0.0167564    -1.64562  0.19228    1                9043   0.362637  SCHAEFFER_PROSTATE_DEVELOPMENT_12HR_UP
Variation Smoothness All Timeframes  Schaeffer Prostate Development 48Hr Dn                      -0.50045   0.00642055   -1.81633  0.191445   1                8542   0.478417  SCHAEFFER_PROSTATE_DEVELOPMENT_48HR_DN
Mean Smoothness All Timeframes       Schaeffer Prostate Development 48Hr Dn                      -0.520587  0.00320256   -1.89076  0.117961   1                8805   0.478417  SCHAEFFER_PROSTATE_DEVELOPMENT_48HR_DN
Variation Smoothness All Timeframes  Schaeffer Prostate Development And Cancer Box5 Up           -0.838679  0.00179677   -1.7894   0.193127   1                9841   0.9       SCHAEFFER_PROSTATE_DEVELOPMENT_AND_CANCER_BOX5_UP
Mean Smoothness All Timeframes       Schaeffer Prostate Development And Cancer Box5 Up           -0.839687  0.00140112   -1.79759  0.141179   1                9822   0.9       SCHAEFFER_PROSTATE_DEVELOPMENT_AND_CANCER_BOX5_UP
Irregularity                         Rickman Head And Neck Cancer A                               0.641435  0.00121212    1.93846  0.12903    1                2005   0.547619  RICKMAN_HEAD_AND_NECK_CANCER_A
Variation Smoothness All Timeframes  Rickman Head And Neck Cancer A                              -0.538052  0.0232935    -1.63682  0.236712   1                8330   0.547619  RICKMAN_HEAD_AND_NECK_CANCER_A
Mean Smoothness All Timeframes       Rickman Head And Neck Cancer A                              -0.562191  0.0119192    -1.71231  0.160858   1                9669   0.452381  RICKMAN_HEAD_AND_NECK_CANCER_A
Mean Smoothness All Timeframes       Rickman Head And Neck Cancer B                              -0.631041  0.0297352    -1.62845  0.194532   1                8493   0.681818  RICKMAN_HEAD_AND_NECK_CANCER_B
Mean Smoothness All Timeframes       Rickman Head And Neck Cancer C                              -0.586586  0.0249196    -1.59642  0.20775    1               10228   0.346154  RICKMAN_HEAD_AND_NECK_CANCER_C
Variation Smoothness All Timeframes  Rickman Head And Neck Cancer F                              -0.772911  0.000401849  -1.91161  0.168147   1               10457   0.6       RICKMAN_HEAD_AND_NECK_CANCER_F
Mean Smoothness All Timeframes       Rickman Head And Neck Cancer F                              -0.759935  0.00119379   -1.89013  0.117961   1               10303   0.64      RICKMAN_HEAD_AND_NECK_CANCER_F
Variation Smoothness All Timeframes  Colin Pilocytic Astrocytoma Vs Glioblastoma Up              -0.569696  0.0291701    -1.60295  0.247184   1                9430   0.478261  COLIN_PILOCYTIC_ASTROCYTOMA_VS_GLIOBLASTOMA_UP
Mean Smoothness All Timeframes       Colin Pilocytic Astrocytoma Vs Glioblastoma Up              -0.640193  0.00390304   -1.80257  0.141179   1                9983   0.478261  COLIN_PILOCYTIC_ASTROCYTOMA_VS_GLIOBLASTOMA_UP
Variation Smoothness All Timeframes  Baldwin Prkci Targets Up                                    -0.52715   0.0221729    -1.63113  0.238002   1                9478   0.419355  BALDWIN_PRKCI_TARGETS_UP
Mean Smoothness All Timeframes       Baldwin Prkci Targets Up                                    -0.552187  0.0115082    -1.70806  0.164501   1                8052   0.645161  BALDWIN_PRKCI_TARGETS_UP
Washout                              Den Interact With Lca5                                       0.693699  0.00423985    1.84149  0.193047   1                2423   0.8       DEN_INTERACT_WITH_LCA5
Top Late Enhancement                 Den Interact With Lca5                                      -0.671486  0.00929105   -1.77448  0.213855   1                8766   0.8       DEN_INTERACT_WITH_LCA5
Vol Late Lt0                         Den Interact With Lca5                                       0.755798  0.000201776   2.01053  0.206783   0.537328         1906   0.8       DEN_INTERACT_WITH_LCA5
Top Late Enhancement                 Benporath Es 2                                              -0.71786   0.013241     -1.68239  0.230757   1                9182   0.625     BENPORATH_ES_2
Variation Smoothness All Timeframes  Benporath Prc2 Targets                                      -0.492871  0.00853852   -1.74238  0.198874   1                6967   0.634409  BENPORATH_PRC2_TARGETS
Mean Smoothness All Timeframes       Benporath Prc2 Targets                                      -0.51027   0.00376014   -1.8068   0.141179   1                8417   0.467742  BENPORATH_PRC2_TARGETS
Washout                              Benporath Proliferation                                      0.759029  0.00262361    1.8653   0.193047   1                1544   0.707317  BENPORATH_PROLIFERATION
Top Late Enhancement                 Benporath Proliferation                                     -0.764524  0.0020141    -1.8759   0.213855   1                9878   0.707317  BENPORATH_PROLIFERATION
Mean Smoothness All Timeframes       Stark Prefrontal Cortex 22Q11 Deletion Up                   -0.432699  0.037337     -1.56738  0.21886    1                7624   0.490566  STARK_PREFRONTAL_CORTEX_22Q11_DELETION_UP
Variation Smoothness All Timeframes  Stark Hyppocampus 22Q11 Deletion Up                         -0.593288  0.000783239  -1.8756   0.168147   1                8657   0.605263  STARK_HYPPOCAMPUS_22Q11_DELETION_UP
Mean Smoothness All Timeframes       Stark Hyppocampus 22Q11 Deletion Up                         -0.527386  0.0195618    -1.66739  0.180662   1                8495   0.552632  STARK_HYPPOCAMPUS_22Q11_DELETION_UP
Mean Smoothness All Timeframes       Amit Egf Response 40 Hela                                   -0.558989  0.0895342    -1.54759  0.227442   1                9339   0.459459  AMIT_EGF_RESPONSE_40_HELA
Washout                              Amit Egf Response 20 Mcf10A                                  0.672635  0.00490697    1.84377  0.193047   1                2558   0.642857  AMIT_EGF_RESPONSE_20_MCF10A
Variation Smoothness All Timeframes  Amit Egf Response 40 Mcf10A                                 -0.58732   0.0436467    -1.62641  0.242889   1               10412   0.421053  AMIT_EGF_RESPONSE_40_MCF10A
Mean Smoothness All Timeframes       Amit Egf Response 40 Mcf10A                                 -0.563094  0.0636107    -1.56209  0.221717   1                9638   0.473684  AMIT_EGF_RESPONSE_40_MCF10A
Mean Smoothness All Timeframes       Amit Serum Response 60 Mcf10A                               -0.544461  0.0382987    -1.62643  0.196012   1                9575   0.37037   AMIT_SERUM_RESPONSE_60_MCF10A
Variation Smoothness All Timeframes  Lin Silenced By Tumor Microenvironment                      -0.518141  0.00928355   -1.77388  0.193127   1                7908   0.561644  LIN_SILENCED_BY_TUMOR_MICROENVIRONMENT
Mean Smoothness All Timeframes       Lin Silenced By Tumor Microenvironment                      -0.500162  0.0127711    -1.71775  0.16069    1                8226   0.506849  LIN_SILENCED_BY_TUMOR_MICROENVIRONMENT
Mean Smoothness All Timeframes       Sung Metastasis Stroma Up                                   -0.502699  0.0164362    -1.75326  0.153387   1                7386   0.613861  SUNG_METASTASIS_STROMA_UP
Washout                              Mori Pre Bi Lymphocyte Up                                    0.676604  0.0169146     1.7922   0.221302   1                1702   0.621212  MORI_PRE_BI_LYMPHOCYTE_UP
Top Late Enhancement                 Mori Pre Bi Lymphocyte Up                                   -0.676885  0.0188031    -1.79891  0.213855   1               10215   0.515152  MORI_PRE_BI_LYMPHOCYTE_UP
Washout                              Mori Large Pre Bii Lymphocyte Up                             0.706812  0.0247336     1.75105  0.246582   1                1585   0.641026  MORI_LARGE_PRE_BII_LYMPHOCYTE_UP
Top Late Enhancement                 Mori Large Pre Bii Lymphocyte Up                            -0.717438  0.021936     -1.76985  0.213855   1               10215   0.564103  MORI_LARGE_PRE_BII_LYMPHOCYTE_UP
Washout                              Mori Mature B Lymphocyte Dn                                  0.635664  0.00101399    1.96101  0.193047   1                1321   0.548387  MORI_MATURE_B_LYMPHOCYTE_DN
Top Late Enhancement                 Mori Mature B Lymphocyte Dn                                 -0.61135   0.00403796   -1.8855   0.213855   1                9964   0.532258  MORI_MATURE_B_LYMPHOCYTE_DN
Vol Late Lt0                         Mori Mature B Lymphocyte Dn                                  0.598302  0.00525784    1.84595  0.214618   1                1210   0.483871  MORI_MATURE_B_LYMPHOCYTE_DN
Variation Smoothness All Timeframes  Lee Targets Of Ptch1 And Sufu Dn                            -0.476013  0.00925367   -1.63108  0.238002   1                6768   0.644444  LEE_TARGETS_OF_PTCH1_AND_SUFU_DN
Mean Smoothness All Timeframes       Lee Targets Of Ptch1 And Sufu Dn                            -0.492122  0.00621741   -1.69223  0.169834   1                7144   0.644444  LEE_TARGETS_OF_PTCH1_AND_SUFU_DN
Uptake Speed                         Nikolsky Breast Cancer 17Q21 Q25 Amplicon                   -0.683145  0.00414365   -2.02137  0.187252   1                9456   0.61674   NIKOLSKY_BREAST_CANCER_17Q21_Q25_AMPLICON
Mean Smoothness All Timeframes       Whitehurst Paclitaxel Sensitivity                           -0.500016  0.0124373    -1.66148  0.181553   1                9661   0.4       WHITEHURST_PACLITAXEL_SENSITIVITY
Mean Smoothness All Timeframes       Ji Metastasis Repressed By Stk11                            -0.595733  0.0315832    -1.637    0.194408   1                8749   0.545455  JI_METASTASIS_REPRESSED_BY_STK11
Variation Smoothness All Timeframes  Onder Cdh1 Targets 2 Up                                     -0.572897  0.0335717    -1.75608  0.193127   1                8368   0.611374  ONDER_CDH1_TARGETS_2_UP
Mean Smoothness All Timeframes       Onder Cdh1 Targets 2 Up                                     -0.611419  0.0119593    -1.87067  0.12188    1                9244   0.554502  ONDER_CDH1_TARGETS_2_UP
Mean Smoothness All Timeframes       Onder Cdh1 Signaling Via Ctnnb1                             -0.50683   0.0902361    -1.5321   0.233789   1                8945   0.426667  ONDER_CDH1_SIGNALING_VIA_CTNNB1
Variation Smoothness All Timeframes  Cervera Sdhb Targets 2                                      -0.48855   0.0126608    -1.70413  0.217973   1                8349   0.45      CERVERA_SDHB_TARGETS_2
Mean Smoothness All Timeframes       Cervera Sdhb Targets 2                                      -0.530038  0.00221908   -1.84268  0.133542   1                8876   0.475     CERVERA_SDHB_TARGETS_2
Variation Smoothness All Timeframes  Rozanov Mmp14 Targets Dn                                    -0.57175   0.00843373   -1.76814  0.193127   1                8248   0.565217  ROZANOV_MMP14_TARGETS_DN
Mean Smoothness All Timeframes       Rozanov Mmp14 Targets Dn                                    -0.55729   0.0133038    -1.71413  0.160858   1                8287   0.565217  ROZANOV_MMP14_TARGETS_DN
Variation Smoothness All Timeframes  Cervera Sdhb Targets 1 Up                                   -0.534058  0.0126455    -1.76321  0.193127   1                8315   0.5       CERVERA_SDHB_TARGETS_1_UP
Mean Smoothness All Timeframes       Cervera Sdhb Targets 1 Up                                   -0.545754  0.00759089   -1.80088  0.141179   1                8277   0.527027  CERVERA_SDHB_TARGETS_1_UP
Variation Smoothness All Timeframes  Nelson Response To Androgen Dn                              -0.601317  0.00633413   -1.76624  0.193127   1               10772   0.25      NELSON_RESPONSE_TO_ANDROGEN_DN
Mean Smoothness All Timeframes       Nelson Response To Androgen Dn                              -0.600713  0.00589391   -1.7656   0.153387   1               10665   0.25      NELSON_RESPONSE_TO_ANDROGEN_DN
Mean Smoothness All Timeframes       Ross Aml With Pml Rara Fusion                               -0.433727  0.0301359    -1.57518  0.216178   1                8559   0.410714  ROSS_AML_WITH_PML_RARA_FUSION
Top Late Enhancement                 Houstis Ros                                                 -0.58879   0.0137152    -1.71899  0.221263   1                9813   0.566667  HOUSTIS_ROS
Top Late Enhancement                 Basso B Lymphocyte Network                                  -0.496671  0.014377     -1.71073  0.226338   1                9762   0.48062   BASSO_B_LYMPHOCYTE_NETWORK
Variation Smoothness Uptake          Basso B Lymphocyte Network                                   0.513056  0.00670876    1.77134  0.192501   1                1453   0.449612  BASSO_B_LYMPHOCYTE_NETWORK
Washout                              Menssen Myc Targets                                          0.645461  0.0128023     1.7648   0.237919   1                2465   0.705882  MENSSEN_MYC_TARGETS
Top Late Enhancement                 Menssen Myc Targets                                         -0.635664  0.0192737    -1.73317  0.213855   1                9204   0.607843  MENSSEN_MYC_TARGETS
Vol Late Lt0                         Menssen Myc Targets                                          0.677846  0.00366077    1.85259  0.214618   1                1577   0.627451  MENSSEN_MYC_TARGETS
Variation Smoothness All Timeframes  Menssen Myc Targets                                          0.653462  0.0130156     1.79757  0.247899   1                2570   0.705882  MENSSEN_MYC_TARGETS
Variation Smoothness Uptake          Menssen Myc Targets                                          0.659561  0.0109113     1.8086   0.174539   1                2360   0.705882  MENSSEN_MYC_TARGETS
Ld Late Lt0                          Sana Tnf Signaling Dn                                        0.591725  0.00181014    1.9758   0.180214   1                1735   0.481013  SANA_TNF_SIGNALING_DN
Variation Smoothness All Timeframes  Sana Tnf Signaling Dn                                       -0.523576  0.017734     -1.76911  0.193127   1                8931   0.493671  SANA_TNF_SIGNALING_DN
Mean Smoothness All Timeframes       Sana Tnf Signaling Dn                                       -0.556941  0.00576083   -1.87229  0.12188    1                9487   0.468354  SANA_TNF_SIGNALING_DN
Washout                              Peng Leucine Deprivation Dn                                  0.592511  0.00686176    1.84249  0.193047   1                2382   0.596591  PENG_LEUCINE_DEPRIVATION_DN
Top Late Enhancement                 Peng Leucine Deprivation Dn                                 -0.625474  0.000994233  -1.94591  0.213855   1                8861   0.625     PENG_LEUCINE_DEPRIVATION_DN
Vol Late Lt0                         Peng Leucine Deprivation Dn                                  0.614825  0.00219868    1.92236  0.214618   1                2278   0.630682  PENG_LEUCINE_DEPRIVATION_DN
Variation Smoothness Uptake          Peng Leucine Deprivation Dn                                  0.579207  0.0113253     1.80202  0.174539   1                2521   0.596591  PENG_LEUCINE_DEPRIVATION_DN
Top Late Enhancement                 Yagi Aml Relapse Prognosis                                  -0.564861  0.00645161   -1.73962  0.213855   1               10633   0.321429  YAGI_AML_RELAPSE_PROGNOSIS
Top Late Enhancement                 Manalo Hypoxia Dn                                           -0.574353  0.0504643    -1.66149  0.248326   1                9019   0.54918   MANALO_HYPOXIA_DN
Variation Smoothness All Timeframes  Lee Liver Cancer Myc E2F1 Dn                                -0.531172  0.00632911   -1.77135  0.193127   1                8713   0.433333  LEE_LIVER_CANCER_MYC_E2F1_DN
Washout                              Zucchi Metastasis Up                                         0.557074  0.00159458    1.87734  0.193047   1                1192   0.388889  ZUCCHI_METASTASIS_UP
Top Late Enhancement                 Zucchi Metastasis Up                                        -0.508689  0.00856915   -1.70673  0.226338   1               10049   0.388889  ZUCCHI_METASTASIS_UP
Top Late Enhancement                 Vantveer Breast Cancer Metastasis Dn                        -0.651586  0.0324018    -1.69825  0.226338   1                9457   0.546392  VANTVEER_BREAST_CANCER_METASTASIS_DN
Washout                              Tarte Plasma Cell Vs Plasmablast Dn                          0.597051  0.0042262     1.89799  0.193047   1                1760   0.515789  TARTE_PLASMA_CELL_VS_PLASMABLAST_DN
Top Late Enhancement                 Tarte Plasma Cell Vs Plasmablast Dn                         -0.610673  0.00180469   -1.93681  0.213855   1                9241   0.550877  TARTE_PLASMA_CELL_VS_PLASMABLAST_DN
Vol Late Lt0                         Tarte Plasma Cell Vs Plasmablast Dn                          0.590736  0.004995      1.88117  0.214618   1                2006   0.547368  TARTE_PLASMA_CELL_VS_PLASMABLAST_DN
Variation Smoothness Uptake          Tarte Plasma Cell Vs Plasmablast Dn                          0.566338  0.0141844     1.79765  0.17746    1                2337   0.596491  TARTE_PLASMA_CELL_VS_PLASMABLAST_DN
Mean Smoothness All Timeframes       Jechlinger Epithelial To Mesenchymal Transition Dn          -0.442371  0.0577769    -1.52762  0.236945   1                8178   0.483333  JECHLINGER_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_DN
Washout                              Le Egr2 Targets Up                                           0.631376  0.0297352     1.77464  0.235363   1                2088   0.646465  LE_EGR2_TARGETS_UP
Top Late Enhancement                 Le Egr2 Targets Up                                          -0.64487   0.0218094    -1.81985  0.213855   1                9270   0.656566  LE_EGR2_TARGETS_UP
Top Late Enhancement                 Mootha Voxphos                                              -0.657098  0.0444712    -1.66737  0.24362    1                9044   0.691358  MOOTHA_VOXPHOS
Variation Smoothness All Timeframes  Mootha Voxphos                                               0.698768  0.0156839     1.7909   0.247899   1                2234   0.740741  MOOTHA_VOXPHOS
Variation Smoothness Uptake          Mootha Voxphos                                               0.717322  0.00817221    1.83702  0.173515   1                1857   0.728395  MOOTHA_VOXPHOS
Mean Smoothness All Timeframes       Pomeroy Medulloblastoma Prognosis Up                        -0.496358  0.0146651    -1.67243  0.178057   1                9283   0.375     POMEROY_MEDULLOBLASTOMA_PROGNOSIS_UP
Variation Smoothness Uptake          Schuhmacher Myc Targets Up                                   0.61958   0.0184258     1.76184  0.192501   1                2491   0.671233  SCHUHMACHER_MYC_TARGETS_UP
Mean Smoothness All Timeframes       Sweet Kras Targets Up                                       -0.57212   0.0666667    -1.62556  0.196481   1                8738   0.519481  SWEET_KRAS_TARGETS_UP
Variation Smoothness All Timeframes  Shipp Dlbcl Vs Follicular Lymphoma Dn                       -0.560209  0.0281833    -1.65124  0.236572   1               10095   0.354839  SHIPP_DLBCL_VS_FOLLICULAR_LYMPHOMA_DN
Mean Smoothness All Timeframes       Shipp Dlbcl Vs Follicular Lymphoma Dn                       -0.513842  0.065867     -1.51436  0.24733    1                8575   0.483871  SHIPP_DLBCL_VS_FOLLICULAR_LYMPHOMA_DN
Washout                              Sana Response To Ifng Dn                                     0.594394  0.00101276    2.00261  0.193047   1                1820   0.5       SANA_RESPONSE_TO_IFNG_DN
Top Late Enhancement                 Sana Response To Ifng Dn                                    -0.592853  0.00120048   -1.99895  0.213855   1                9275   0.513158  SANA_RESPONSE_TO_IFNG_DN
Vol Late Lt0                         Sana Response To Ifng Dn                                     0.560316  0.00201248    1.89544  0.214618   1                1311   0.394737  SANA_RESPONSE_TO_IFNG_DN
Mean Smoothness All Timeframes       Lee Liver Cancer Myc Dn                                     -0.452614  0.0374554    -1.51591  0.246789   1                9849   0.277778  LEE_LIVER_CANCER_MYC_DN
Irregularity                         Lee Liver Cancer Ciprofibrate Up                             0.54878   0.00474308    1.8446   0.218099   1                 536   0.232558  LEE_LIVER_CANCER_CIPROFIBRATE_UP
Ld Late Lt0                          Lee Liver Cancer Ciprofibrate Up                             0.628932  0.00020012    2.10507  0.0853978  0.53292          1449   0.488372  LEE_LIVER_CANCER_CIPROFIBRATE_UP
Variation Smoothness All Timeframes  Lee Liver Cancer Ciprofibrate Up                            -0.5026    0.0136793    -1.68175  0.231284   1                8578   0.395349  LEE_LIVER_CANCER_CIPROFIBRATE_UP
Mean Smoothness All Timeframes       Lee Liver Cancer Ciprofibrate Up                            -0.471162  0.0270054    -1.57987  0.216167   1                7519   0.55814   LEE_LIVER_CANCER_CIPROFIBRATE_UP
Variation Smoothness All Timeframes  Iizuka Liver Cancer Progression L1 G1 Up                    -0.630292  0.0317009    -1.63703  0.236712   1                9839   0.5       IIZUKA_LIVER_CANCER_PROGRESSION_L1_G1_UP
Top Late Enhancement                 Pomeroy Medulloblastoma Prognosis Dn                        -0.611605  0.0258065    -1.67607  0.232024   1                9929   0.375     POMEROY_MEDULLOBLASTOMA_PROGNOSIS_DN
Variation Smoothness All Timeframes  Pomeroy Medulloblastoma Prognosis Dn                         0.692389  0.000990688   1.91988  0.206872   1                1860   0.65      POMEROY_MEDULLOBLASTOMA_PROGNOSIS_DN
Variation Smoothness Uptake          Pomeroy Medulloblastoma Prognosis Dn                         0.703249  0.000594295   1.94834  0.173515   1                1822   0.675     POMEROY_MEDULLOBLASTOMA_PROGNOSIS_DN
Mean Smoothness All Timeframes       Manalo Hypoxia Up                                           -0.485274  0.0462626    -1.60957  0.202694   1                8749   0.435028  MANALO_HYPOXIA_UP
Washout                              Shipp Dlbcl Vs Follicular Lymphoma Up                        0.720942  0.00547556    1.81547  0.206704   1                1796   0.682927  SHIPP_DLBCL_VS_FOLLICULAR_LYMPHOMA_UP
Top Late Enhancement                 Shipp Dlbcl Vs Follicular Lymphoma Up                       -0.69825   0.0137374    -1.75719  0.213855   1                8892   0.731707  SHIPP_DLBCL_VS_FOLLICULAR_LYMPHOMA_UP
Vol Late Lt0                         Shipp Dlbcl Vs Follicular Lymphoma Up                        0.727374  0.00519481    1.84649  0.214618   1                2000   0.707317  SHIPP_DLBCL_VS_FOLLICULAR_LYMPHOMA_UP
Variation Smoothness Uptake          Shipp Dlbcl Vs Follicular Lymphoma Up                        0.695177  0.0149931     1.7646   0.192501   1                1618   0.634146  SHIPP_DLBCL_VS_FOLLICULAR_LYMPHOMA_UP
Washout                              Peng Rapamycin Response Dn                                   0.567568  0.00543588    1.84721  0.193047   1                2506   0.596491  PENG_RAPAMYCIN_RESPONSE_DN
Ser                                  Peng Rapamycin Response Dn                                   0.59957   0.00100888    1.95901  0.21078    1                2233   0.592105  PENG_RAPAMYCIN_RESPONSE_DN
Top Late Enhancement                 Peng Rapamycin Response Dn                                  -0.587056  0.00339051   -1.90874  0.213855   1                8929   0.592105  PENG_RAPAMYCIN_RESPONSE_DN
Vol Late Lt0                         Peng Rapamycin Response Dn                                   0.565359  0.00664385    1.84227  0.214618   1                2144   0.561404  PENG_RAPAMYCIN_RESPONSE_DN
Variation Smoothness All Timeframes  Peng Rapamycin Response Dn                                   0.568078  0.00878594    1.8432   0.212371   1                2785   0.605263  PENG_RAPAMYCIN_RESPONSE_DN
Variation Smoothness Uptake          Peng Rapamycin Response Dn                                   0.581299  0.00401768    1.88481  0.173515   1                2521   0.605263  PENG_RAPAMYCIN_RESPONSE_DN
Variation Smoothness All Timeframes  Goldrath Immune Memory                                      -0.509564  0.0035693    -1.80629  0.192967   1                9419   0.365385  GOLDRATH_IMMUNE_MEMORY
Mean Smoothness All Timeframes       Goldrath Immune Memory                                      -0.470873  0.0163152    -1.66592  0.181553   1                9349   0.384615  GOLDRATH_IMMUNE_MEMORY
Top Late Enhancement                 Bhattacharya Embryonic Stem Cell                            -0.57297   0.0187689    -1.73597  0.213855   1                9204   0.549296  BHATTACHARYA_EMBRYONIC_STEM_CELL
Variation Smoothness All Timeframes  Yao Hoxa10 Targets Via Progesterone Up                      -0.545461  0.012249     -1.77183  0.193127   1               10017   0.344828  YAO_HOXA10_TARGETS_VIA_PROGESTERONE_UP
Mean Smoothness All Timeframes       Yao Hoxa10 Targets Via Progesterone Up                      -0.546635  0.0106234    -1.77851  0.149917   1                9838   0.37931   YAO_HOXA10_TARGETS_VIA_PROGESTERONE_UP
Mean Smoothness All Timeframes       Zhang Targets Of Ewsr1 Fli1 Fusion                          -0.420519  0.0369318    -1.5496   0.226989   1                8068   0.459459  ZHANG_TARGETS_OF_EWSR1_FLI1_FUSION
Top Late Enhancement                 Zhan Multiple Myeloma Subgroups                             -0.740196  0.00460184   -1.79229  0.213855   1                9230   0.724138  ZHAN_MULTIPLE_MYELOMA_SUBGROUPS
Variation Smoothness All Timeframes  Haddad T Lymphocyte And Nk Progenitor Up                    -0.527206  0.00100301   -1.83978  0.171536   1                7741   0.640625  HADDAD_T_LYMPHOCYTE_AND_NK_PROGENITOR_UP
Mean Smoothness All Timeframes       Haddad T Lymphocyte And Nk Progenitor Up                    -0.569458  0.000201532  -1.98764  0.109401   0.536679         7899   0.6875    HADDAD_T_LYMPHOCYTE_AND_NK_PROGENITOR_UP
Mean Smoothness All Timeframes       Vernell Retinoblastoma Pathway Dn                           -0.590006  0.0136766    -1.70298  0.165416   1               10045   0.4       VERNELL_RETINOBLASTOMA_PATHWAY_DN
Variation Smoothness All Timeframes  Wang Immortalized By Hoxa9 And Meis1 Dn                     -0.765497  0.00239904   -1.9149   0.168147   1               10560   0.470588  WANG_IMMORTALIZED_BY_HOXA9_AND_MEIS1_DN
Mean Smoothness All Timeframes       Wang Immortalized By Hoxa9 And Meis1 Dn                     -0.707276  0.0120895    -1.76925  0.15255    1                9102   0.588235  WANG_IMMORTALIZED_BY_HOXA9_AND_MEIS1_DN
Variation Smoothness All Timeframes  Alcalay Aml By Npm1 Localization Up                         -0.412683  0.00548446   -1.67377  0.233816   1                7682   0.45283   ALCALAY_AML_BY_NPM1_LOCALIZATION_UP
Mean Smoothness All Timeframes       Alcalay Aml By Npm1 Localization Up                         -0.412079  0.00543916   -1.67655  0.176752   1                8498   0.396226  ALCALAY_AML_BY_NPM1_LOCALIZATION_UP
Top Late Enhancement                 Yu Myc Targets Up                                           -0.789517  0.00720865   -1.79237  0.213855   1               10215   0.684211  YU_MYC_TARGETS_UP
Ld Late Lt0                          Guo Hex Targets Dn                                           0.567994  0.000607656   2.00687  0.156851   1                2499   0.46      GUO_HEX_TARGETS_DN
Top Late Enhancement                 Ferrando T All With Mll Enl Fusion Dn                       -0.573463  0.0139085    -1.7471   0.213855   1                8963   0.578947  FERRANDO_T_ALL_WITH_MLL_ENL_FUSION_DN
Top Late Enhancement                 Zhan Multiple Myeloma Pr Up                                 -0.87614   0.0139281    -1.68216  0.230757   1               10197   0.911765  ZHAN_MULTIPLE_MYELOMA_PR_UP
Mean Smoothness All Timeframes       Kang Immortalized By Tert Up                                -0.440785  0.0437258    -1.5439   0.228407   1                9057   0.412698  KANG_IMMORTALIZED_BY_TERT_UP
Ld Late Lt0                          Li Wilms Tumor Vs Fetal Kidney 1 Up                          0.524333  0.00140562    1.93342  0.24158    1                2062   0.410714  LI_WILMS_TUMOR_VS_FETAL_KIDNEY_1_UP
Variation Smoothness All Timeframes  Li Wilms Tumor Vs Fetal Kidney 1 Up                         -0.44932   0.0272287    -1.65287  0.236572   1                9295   0.363095  LI_WILMS_TUMOR_VS_FETAL_KIDNEY_1_UP
Mean Smoothness All Timeframes       Li Wilms Tumor Vs Fetal Kidney 1 Up                         -0.456414  0.0223063    -1.68483  0.172418   1                8271   0.488095  LI_WILMS_TUMOR_VS_FETAL_KIDNEY_1_UP
Mean Smoothness All Timeframes       Affar Yy1 Targets Up                                        -0.388998  0.0198965    -1.58477  0.216099   1                8493   0.367188  AFFAR_YY1_TARGETS_UP
Top Late Enhancement                 Affar Yy1 Targets Dn                                        -0.471053  0.0491607    -1.65973  0.248326   1               10180   0.341615  AFFAR_YY1_TARGETS_DN
Variation Smoothness All Timeframes  Kang Immortalized By Tert Dn                                -0.549152  0.00139553   -1.86524  0.168147   1                8252   0.596774  KANG_IMMORTALIZED_BY_TERT_DN
Mean Smoothness All Timeframes       Kang Immortalized By Tert Dn                                -0.524652  0.00579768   -1.77806  0.149917   1                9061   0.5       KANG_IMMORTALIZED_BY_TERT_DN
Variation Smoothness All Timeframes  Yao Hoxa10 Targets Via Progesterone Dn                      -0.723857  0.00637196   -1.78229  0.193127   1                8984   0.583333  YAO_HOXA10_TARGETS_VIA_PROGESTERONE_DN
Mean Smoothness All Timeframes       Yao Hoxa10 Targets Via Progesterone Dn                      -0.697532  0.0124901    -1.71659  0.16069    1                9539   0.583333  YAO_HOXA10_TARGETS_VIA_PROGESTERONE_DN
Variation Smoothness All Timeframes  Dorsey Gab2 Targets                                         -0.690169  0.0274955    -1.67336  0.233816   1                9365   0.611111  DORSEY_GAB2_TARGETS
Mean Smoothness All Timeframes       Dorsey Gab2 Targets                                         -0.715605  0.0145418    -1.74076  0.156443   1               10199   0.555556  DORSEY_GAB2_TARGETS
Variation Smoothness All Timeframes  Zhan Multiple Myeloma Ms Up                                 -0.539867  0.0425237    -1.61011  0.244212   1                8830   0.5       ZHAN_MULTIPLE_MYELOMA_MS_UP
Mean Smoothness All Timeframes       Zhan Multiple Myeloma Ms Up                                 -0.588215  0.0157622    -1.75013  0.153387   1                9818   0.444444  ZHAN_MULTIPLE_MYELOMA_MS_UP
Mean Smoothness All Timeframes       Nouzova Methylated In Apl                                   -0.468721  0.00974817   -1.62909  0.194532   1                8969   0.416667  NOUZOVA_METHYLATED_IN_APL
Variation Smoothness All Timeframes  Verhaak Aml With Npm1 Mutated Dn                            -0.498444  0.010649     -1.77049  0.193127   1                8928   0.421622  VERHAAK_AML_WITH_NPM1_MUTATED_DN
Mean Smoothness All Timeframes       Verhaak Aml With Npm1 Mutated Dn                            -0.493377  0.0152795    -1.75351  0.153387   1                8657   0.454054  VERHAAK_AML_WITH_NPM1_MUTATED_DN
Mean Smoothness All Timeframes       Rorie Targets Of Ewsr1 Fli1 Fusion Dn                       -0.774682  0.00373057   -1.6352   0.194532   1                9897   0.5       RORIE_TARGETS_OF_EWSR1_FLI1_FUSION_DN
Top Late Enhancement                 Greenbaum E2A Targets Up                                    -0.76681   0.0212422    -1.75287  0.213855   1                9885   0.71875   GREENBAUM_E2A_TARGETS_UP
Top Late Enhancement                 Petrova Prox1 Targets Up                                    -0.661881  0.043039     -1.66477  0.24649    1                9497   0.625     PETROVA_PROX1_TARGETS_UP
Mean Smoothness All Timeframes       Hoegerkorp Cd44 Targets Direct Up                           -0.577908  0.0466219    -1.56427  0.219929   1                9643   0.375     HOEGERKORP_CD44_TARGETS_DIRECT_UP
Variation Smoothness All Timeframes  Martinez Response To Trabectedin                            -0.532986  0.0116255    -1.71141  0.213798   1                8413   0.477273  MARTINEZ_RESPONSE_TO_TRABECTEDIN
Mean Smoothness All Timeframes       Martinez Response To Trabectedin                            -0.515455  0.0195604    -1.66305  0.181553   1                7440   0.613636  MARTINEZ_RESPONSE_TO_TRABECTEDIN
Variation Smoothness All Timeframes  Takao Response To Uvb Radiation Up                           0.627485  0.0123821     1.81173  0.237345   1                2404   0.630137  TAKAO_RESPONSE_TO_UVB_RADIATION_UP
Variation Smoothness Uptake          Takao Response To Uvb Radiation Up                           0.657174  0.00413793    1.89703  0.173515   1                1711   0.575342  TAKAO_RESPONSE_TO_UVB_RADIATION_UP
Irregularity                         Urs Adipocyte Differentiation Up                             0.560664  0.000588697   1.9491   0.12903    1                1345   0.313725  URS_ADIPOCYTE_DIFFERENTIATION_UP
Variation Smoothness All Timeframes  Urs Adipocyte Differentiation Up                            -0.544667  0.00277447   -1.89447  0.168147   1                9359   0.45098   URS_ADIPOCYTE_DIFFERENTIATION_UP
Mean Smoothness All Timeframes       Urs Adipocyte Differentiation Up                            -0.525463  0.0042245    -1.827    0.133542   1                9142   0.470588  URS_ADIPOCYTE_DIFFERENTIATION_UP
Variation Smoothness All Timeframes  Traynor Rett Syndrom Up                                     -0.617688  0.0154472    -1.75166  0.193127   1                9003   0.666667  TRAYNOR_RETT_SYNDROM_UP
Mean Smoothness All Timeframes       Traynor Rett Syndrom Up                                     -0.667002  0.00362611   -1.8963   0.117961   1                9976   0.606061  TRAYNOR_RETT_SYNDROM_UP
Variation Smoothness All Timeframes  Suzuki Response To Tsa And Decitabine 1A                    -0.745328  0.00259688   -1.86158  0.168147   1               10114   0.529412  SUZUKI_RESPONSE_TO_TSA_AND_DECITABINE_1A
Mean Smoothness All Timeframes       Suzuki Response To Tsa And Decitabine 1A                    -0.753007  0.00340955   -1.87774  0.117961   1               10001   0.588235  SUZUKI_RESPONSE_TO_TSA_AND_DECITABINE_1A
Variation Smoothness All Timeframes  Su Pancreas                                                 -0.572191  0.00059512   -1.80487  0.192967   1                7574   0.625     SU_PANCREAS
Mean Smoothness All Timeframes       Su Pancreas                                                 -0.522756  0.00955604   -1.64599  0.19228    1                7529   0.583333  SU_PANCREAS
Washout                              Zamora Nos2 Targets Up                                       0.605636  0.00488003    1.82216  0.206704   1                2088   0.571429  ZAMORA_NOS2_TARGETS_UP
Top Late Enhancement                 Zamora Nos2 Targets Up                                      -0.600846  0.00785498   -1.80641  0.213855   1                9285   0.539683  ZAMORA_NOS2_TARGETS_UP
Variation Smoothness All Timeframes  Zamora Nos2 Targets Up                                       0.624426  0.00298924    1.8708   0.212371   1                1646   0.571429  ZAMORA_NOS2_TARGETS_UP
Variation Smoothness Uptake          Zamora Nos2 Targets Up                                       0.62123   0.00455446    1.86958  0.173515   1                2066   0.619048  ZAMORA_NOS2_TARGETS_UP
Top Late Enhancement                 Weigel Oxidative Stress Response                            -0.502654  0.0104691    -1.6789   0.230757   1                9586   0.4375    WEIGEL_OXIDATIVE_STRESS_RESPONSE
Variation Smoothness All Timeframes  Joseph Response To Sodium Butyrate Dn                       -0.465831  0.0205778    -1.62524  0.243827   1                7741   0.545455  JOSEPH_RESPONSE_TO_SODIUM_BUTYRATE_DN
Mean Smoothness All Timeframes       Joseph Response To Sodium Butyrate Dn                       -0.474965  0.0144354    -1.66163  0.181553   1                8163   0.5       JOSEPH_RESPONSE_TO_SODIUM_BUTYRATE_DN
Variation Smoothness All Timeframes  Ivanova Hematopoiesis Stem Cell Long Term                   -0.410447  0.00259326   -1.71258  0.213798   1                7463   0.556075  IVANOVA_HEMATOPOIESIS_STEM_CELL_LONG_TERM
Mean Smoothness All Timeframes       Ivanova Hematopoiesis Stem Cell Long Term                   -0.425536  0.000602894  -1.76966  0.15255    1                8220   0.481308  IVANOVA_HEMATOPOIESIS_STEM_CELL_LONG_TERM
Mean Smoothness All Timeframes       Weigel Oxidative Stress By Hne And H2O2                     -0.424241  0.0294821    -1.52054  0.242693   1               10443   0.264706  WEIGEL_OXIDATIVE_STRESS_BY_HNE_AND_H2O2
Mean Smoothness All Timeframes       Browne Hcmv Infection 16Hr Dn                               -0.464353  0.0320048    -1.5935   0.210308   1                8847   0.444444  BROWNE_HCMV_INFECTION_16HR_DN
Mean Smoothness All Timeframes       Weston Vegfa Targets 6Hr                                    -0.625835  0.033386     -1.72819  0.16069    1                8928   0.627907  WESTON_VEGFA_TARGETS_6HR
Washout                              Rhodes Cancer Meta Signature                                 0.704675  0.00121704    1.93197  0.193047   1                1341   0.59322   RHODES_CANCER_META_SIGNATURE
Ser                                  Rhodes Cancer Meta Signature                                 0.699697  0.00183225    1.93409  0.249199   1                2417   0.728814  RHODES_CANCER_META_SIGNATURE
Top Late Enhancement                 Rhodes Cancer Meta Signature                                -0.714751  0.000604839  -1.95742  0.213855   1                9375   0.677966  RHODES_CANCER_META_SIGNATURE
Vol Late Lt0                         Rhodes Cancer Meta Signature                                 0.702675  0.00080402    1.93936  0.214618   1                2759   0.864407  RHODES_CANCER_META_SIGNATURE
Variation Smoothness All Timeframes  Rhodes Cancer Meta Signature                                 0.666107  0.00769231    1.84082  0.212371   1                2123   0.677966  RHODES_CANCER_META_SIGNATURE
Variation Smoothness Uptake          Rhodes Cancer Meta Signature                                 0.681996  0.00415677    1.87918  0.173515   1                1966   0.661017  RHODES_CANCER_META_SIGNATURE
Mean Smoothness All Timeframes       Nielsen Malignat Fibrous Histiocytoma Up                    -0.701285  0.0883526    -1.51264  0.248519   1                9420   0.666667  NIELSEN_MALIGNAT_FIBROUS_HISTIOCYTOMA_UP
Ld Late Lt0                          Ruan Response To Tnf Troglitazone Dn                         0.620938  0.00181635    1.9614   0.186615   1                2258   0.533333  RUAN_RESPONSE_TO_TNF_TROGLITAZONE_DN
Mean Smoothness All Timeframes       Ruan Response To Tnf Troglitazone Dn                        -0.50989   0.0285829    -1.62107  0.199788   1               11196   0.133333  RUAN_RESPONSE_TO_TNF_TROGLITAZONE_DN
Variation Smoothness All Timeframes  Browne Hcmv Infection 10Hr Dn                               -0.482982  0.0215032    -1.61875  0.244212   1                8711   0.560976  BROWNE_HCMV_INFECTION_10HR_DN
Mean Smoothness All Timeframes       Mcclung Creb1 Targets Dn                                    -0.487474  0.00516899   -1.71578  0.160771   1                9862   0.3       MCCLUNG_CREB1_TARGETS_DN
Mean Smoothness All Timeframes       Kaab Heart Atrium Vs Ventricle Up                           -0.448749  0.0147941    -1.69945  0.168479   1                8787   0.402062  KAAB_HEART_ATRIUM_VS_VENTRICLE_UP
Variation Smoothness Uptake          Hu Genotoxic Damage 24Hr                                     0.608607  0.00633162    1.79381  0.17972    1                1711   0.535714  HU_GENOTOXIC_DAMAGE_24HR
Variation Smoothness All Timeframes  Nielsen Gist Vs Synovial Sarcoma Up                         -0.722625  0.0247414    -1.67463  0.233816   1                9482   0.642857  NIELSEN_GIST_VS_SYNOVIAL_SARCOMA_UP
Mean Smoothness All Timeframes       Nielsen Gist Vs Synovial Sarcoma Up                         -0.792014  0.00587282   -1.83657  0.133542   1                9925   0.642857  NIELSEN_GIST_VS_SYNOVIAL_SARCOMA_UP
Mean Smoothness All Timeframes       Cui Tcf21 Targets Up                                        -0.6675    0.021761     -1.74952  0.153387   1                9265   0.642857  CUI_TCF21_TARGETS_UP
Washout                              Vietor Ifrd1 Targets                                         0.63853   0.00793493    1.76891  0.237566   1                1069   0.444444  VIETOR_IFRD1_TARGETS
Top Late Enhancement                 Vietor Ifrd1 Targets                                        -0.629171  0.00856574   -1.74176  0.213855   1               10549   0.388889  VIETOR_IFRD1_TARGETS
Mean Smoothness All Timeframes       Browne Hcmv Infection 24Hr Dn                               -0.442027  0.0786787    -1.52733  0.236945   1                8519   0.467742  BROWNE_HCMV_INFECTION_24HR_DN
Mean Smoothness All Timeframes       Chiba Response To Tsa Up                                    -0.520181  0.0463297    -1.60018  0.204483   1                9634   0.369565  CHIBA_RESPONSE_TO_TSA_UP
Irregularity                         Ramaswamy Metastasis Dn                                      0.595768  0.000393701   2.02179  0.12903    1                1849   0.444444  RAMASWAMY_METASTASIS_DN
Variation Smoothness All Timeframes  Ramaswamy Metastasis Dn                                     -0.500718  0.00949555   -1.70932  0.213798   1                8132   0.527778  RAMASWAMY_METASTASIS_DN
Washout                              Mody Hippocampus Neonatal                                    0.649477  0.00708502    1.78723  0.221302   1                1777   0.607143  MODY_HIPPOCAMPUS_NEONATAL
Top Late Enhancement                 Mody Hippocampus Neonatal                                   -0.606772  0.0264166    -1.66289  0.247946   1                9454   0.535714  MODY_HIPPOCAMPUS_NEONATAL
Variation Smoothness All Timeframes  Browne Hcmv Infection 14Hr Dn                               -0.374065  0.00316957   -1.66633  0.236189   1                8458   0.3875    BROWNE_HCMV_INFECTION_14HR_DN
Mean Smoothness All Timeframes       Browne Hcmv Infection 14Hr Dn                               -0.380307  0.00218644   -1.68913  0.172226   1                8351   0.420833  BROWNE_HCMV_INFECTION_14HR_DN
Mean Smoothness All Timeframes       Delaserna Myod Targets Dn                                   -0.470829  0.0289591    -1.60324  0.204175   1                9122   0.36      DELASERNA_MYOD_TARGETS_DN
Mean Smoothness All Timeframes       Yamazaki Tceb3 Targets Up                                   -0.431678  0.0603258    -1.53807  0.231366   1                8266   0.449664  YAMAZAKI_TCEB3_TARGETS_UP
Mean Smoothness All Timeframes       Cui Tcf21 Targets Dn                                        -0.621367  0.00434869   -1.84113  0.133542   1                8400   0.545455  CUI_TCF21_TARGETS_DN
Mean Smoothness All Timeframes       Weston Vegfa Targets 12Hr                                   -0.638019  0.0358646    -1.67134  0.178706   1                8698   0.692308  WESTON_VEGFA_TARGETS_12HR
Top Late Enhancement                 Natsume Response To Interferon Beta Dn                      -0.573793  0.0117385    -1.76632  0.213855   1                8687   0.681818  NATSUME_RESPONSE_TO_INTERFERON_BETA_DN
Variation Smoothness Uptake          Natsume Response To Interferon Beta Dn                       0.589059  0.00970874    1.81355  0.174539   1                1114   0.431818  NATSUME_RESPONSE_TO_INTERFERON_BETA_DN
Variation Smoothness All Timeframes  Lindvall Immortalized By Tert Up                            -0.468626  0.0155317    -1.65791  0.236572   1                8805   0.424242  LINDVALL_IMMORTALIZED_BY_TERT_UP
Mean Smoothness All Timeframes       Lindvall Immortalized By Tert Up                            -0.486302  0.00956747   -1.71785  0.16069    1                9657   0.363636  LINDVALL_IMMORTALIZED_BY_TERT_UP
Variation Smoothness All Timeframes  Browne Hcmv Infection 20Hr Dn                               -0.52918   0.00518031   -1.82258  0.191445   1                7914   0.604651  BROWNE_HCMV_INFECTION_20HR_DN
Mean Smoothness All Timeframes       Browne Hcmv Infection 20Hr Dn                               -0.553027  0.00159521   -1.90079  0.117961   1                8266   0.616279  BROWNE_HCMV_INFECTION_20HR_DN
Mean Smoothness All Timeframes       Simbulan Parp1 Targets Up                                   -0.678234  0.0179104    -1.74443  0.155059   1                9775   0.571429  SIMBULAN_PARP1_TARGETS_UP
Mean Smoothness All Timeframes       Ule Splicing Via Nova2                                      -0.479561  0.0361277    -1.5543   0.225193   1                7355   0.657143  ULE_SPLICING_VIA_NOVA2
Variation Smoothness All Timeframes  Ivanova Hematopoiesis Stem Cell                             -0.409587  0.000802407  -1.7431   0.198874   1                8274   0.416216  IVANOVA_HEMATOPOIESIS_STEM_CELL
Mean Smoothness All Timeframes       Ivanova Hematopoiesis Stem Cell                             -0.383044  0.00658288   -1.62861  0.194532   1                7930   0.437838  IVANOVA_HEMATOPOIESIS_STEM_CELL
Variation Smoothness All Timeframes  Kuninger Igf1 Vs Pdgfb Targets Dn                           -0.52489   0.0210186    -1.66507  0.236572   1                9421   0.371429  KUNINGER_IGF1_VS_PDGFB_TARGETS_DN
Mean Smoothness All Timeframes       Kuninger Igf1 Vs Pdgfb Targets Dn                           -0.496974  0.0368496    -1.57663  0.216178   1                9347   0.4       KUNINGER_IGF1_VS_PDGFB_TARGETS_DN
Mean Smoothness All Timeframes       Mahajan Response To Il1A Dn                                 -0.456449  0.0146117    -1.67374  0.17801    1                8762   0.41791   MAHAJAN_RESPONSE_TO_IL1A_DN
Circularity                          Ruan Response To Tnf Dn                                     -0.566071  0.00160128   -2.00092  0.205876   1                8045   0.560606  RUAN_RESPONSE_TO_TNF_DN
Irregularity                         Ruan Response To Tnf Dn                                      0.516724  0.0088141     1.82785  0.218099   1                3172   0.5       RUAN_RESPONSE_TO_TNF_DN
Ld Late Lt0                          Ruan Response To Tnf Dn                                      0.58558   0.000605938   2.05504  0.0973724  1                2860   0.590909  RUAN_RESPONSE_TO_TNF_DN
Variation Smoothness All Timeframes  Ruan Response To Tnf Dn                                     -0.465036  0.0267412    -1.65174  0.236572   1                9087   0.348485  RUAN_RESPONSE_TO_TNF_DN
Mean Smoothness All Timeframes       Ruan Response To Tnf Dn                                     -0.455695  0.0326522    -1.61179  0.202062   1                8878   0.363636  RUAN_RESPONSE_TO_TNF_DN
Mean Smoothness All Timeframes       Mariadason Regulated By Histone Acetylation Up              -0.391573  0.0329039    -1.5111   0.248827   1                8025   0.426471  MARIADASON_REGULATED_BY_HISTONE_ACETYLATION_UP
Mean Smoothness All Timeframes       Lindvall Immortalized By Tert Dn                            -0.543412  0.0670572    -1.59444  0.209864   1                9445   0.478261  LINDVALL_IMMORTALIZED_BY_TERT_DN
Top Late Enhancement                 Burton Adipogenesis 5                                       -0.502048  0.0151909    -1.73186  0.213855   1                8994   0.530435  BURTON_ADIPOGENESIS_5
Mean Smoothness All Timeframes       Weston Vegfa Targets                                        -0.490549  0.0613274    -1.58522  0.216099   1                8928   0.469136  WESTON_VEGFA_TARGETS
Mean Smoothness All Timeframes       Chen Lvad Support Of Failing Heart Up                       -0.539407  0.0254867    -1.69839  0.168989   1                8036   0.630952  CHEN_LVAD_SUPPORT_OF_FAILING_HEART_UP
Irregularity                         Browne Hcmv Infection 8Hr Dn                                 0.545958  0.00100583    1.89293  0.185049   1                1893   0.351351  BROWNE_HCMV_INFECTION_8HR_DN
Washout                              Rhodes Undifferentiated Cancer                               0.7986    0.00427003    1.80315  0.217568   1                1321   0.806452  RHODES_UNDIFFERENTIATED_CANCER
Top Late Enhancement                 Rhodes Undifferentiated Cancer                              -0.801373  0.00284611   -1.80084  0.213855   1                9866   0.806452  RHODES_UNDIFFERENTIATED_CANCER
Irregularity                         Lu Aging Brain Up                                            0.460082  0.00518341    1.80275  0.246546   1                2944   0.448276  LU_AGING_BRAIN_UP
Variation Smoothness All Timeframes  Lu Aging Brain Up                                           -0.416436  0.0217042    -1.63516  0.236712   1                8124   0.431034  LU_AGING_BRAIN_UP
Mean Smoothness All Timeframes       Lu Aging Brain Up                                           -0.447886  0.00727273   -1.75543  0.153387   1                7896   0.50431   LU_AGING_BRAIN_UP
Top Late Enhancement                 Burton Adipogenesis 3                                       -0.670308  0.0431698    -1.71106  0.226338   1               10157   0.549451  BURTON_ADIPOGENESIS_3
Variation Smoothness All Timeframes  Ma Pituitary Fetal Vs Adult Up                              -0.575203  0.003003     -1.81446  0.191445   1                9539   0.375     MA_PITUITARY_FETAL_VS_ADULT_UP
Mean Smoothness All Timeframes       Ma Pituitary Fetal Vs Adult Up                              -0.482756  0.0435304    -1.52045  0.242693   1                9204   0.375     MA_PITUITARY_FETAL_VS_ADULT_UP
Mean Smoothness All Timeframes       Verrecchia Early Response To Tgfb1                          -0.540816  0.078806     -1.58875  0.212416   1                9696   0.436364  VERRECCHIA_EARLY_RESPONSE_TO_TGFB1
Irregularity                         Xu Gh1 Exogenous Targets Up                                  0.516393  0.00140168    1.83012  0.218099   1                2623   0.425     XU_GH1_EXOGENOUS_TARGETS_UP
Variation Smoothness All Timeframes  Xu Gh1 Exogenous Targets Up                                 -0.52982   0.0006       -1.87687  0.168147   1                8151   0.575     XU_GH1_EXOGENOUS_TARGETS_UP
Mean Smoothness All Timeframes       Xu Gh1 Exogenous Targets Up                                 -0.448672  0.0157653    -1.58819  0.212532   1                7796   0.6       XU_GH1_EXOGENOUS_TARGETS_UP
Variation Smoothness All Timeframes  Zhu Cmv 8 Hr Dn                                             -0.504589  0.0117273    -1.70115  0.219274   1                8912   0.470588  ZHU_CMV_8_HR_DN
Mean Smoothness All Timeframes       Zhu Cmv 8 Hr Dn                                             -0.493176  0.0199315    -1.65569  0.184559   1                9098   0.45098   ZHU_CMV_8_HR_DN
Mean Smoothness All Timeframes       Zhu Cmv 24 Hr Dn                                            -0.547011  0.0946353    -1.57177  0.218767   1                9321   0.525     ZHU_CMV_24_HR_DN
Variation Smoothness All Timeframes  Pal Prmt5 Targets Dn                                        -0.572595  0.013206     -1.70881  0.213798   1                8513   0.5       PAL_PRMT5_TARGETS_DN
Mean Smoothness All Timeframes       Pal Prmt5 Targets Dn                                        -0.548349  0.0222089    -1.64158  0.194025   1               10579   0.277778  PAL_PRMT5_TARGETS_DN
Irregularity                         Sato Silenced Epigenetically In Pancreatic Cancer            0.649073  0.00121433    1.90347  0.176552   1                 960   0.416667  SATO_SILENCED_EPIGENETICALLY_IN_PANCREATIC_CANCER
Variation Smoothness All Timeframes  Browne Hcmv Infection 18Hr Dn                               -0.469657  0.0332206    -1.65657  0.236572   1                8288   0.507692  BROWNE_HCMV_INFECTION_18HR_DN
Mean Smoothness All Timeframes       Browne Hcmv Infection 18Hr Dn                               -0.486358  0.0225055    -1.712    0.160858   1                8253   0.530769  BROWNE_HCMV_INFECTION_18HR_DN
Washout                              Moreira Response To Tsa Up                                   0.679061  0.00265794    1.88917  0.193047   1                1932   0.62963   MOREIRA_RESPONSE_TO_TSA_UP
Top Late Enhancement                 Moreira Response To Tsa Up                                  -0.643154  0.00545124   -1.80571  0.213855   1               10139   0.481481  MOREIRA_RESPONSE_TO_TSA_UP
Variation Smoothness All Timeframes  Gentile Uv Response Cluster D6                              -0.481851  0.0236253    -1.61662  0.244212   1                9939   0.323529  GENTILE_UV_RESPONSE_CLUSTER_D6
Mean Smoothness All Timeframes       Gentile Uv Response Cluster D6                              -0.494377  0.0166093    -1.66077  0.181553   1                9535   0.352941  GENTILE_UV_RESPONSE_CLUSTER_D6
Variation Smoothness Uptake          Zhang Response To Cantharidin Dn                             0.61156   0.00698603    1.80441  0.174539   1                1390   0.5625    ZHANG_RESPONSE_TO_CANTHARIDIN_DN
Top Late Enhancement                 Burton Adipogenesis Peak At 24Hr                            -0.778578  0.00885847   -1.80807  0.213855   1                9879   0.789474  BURTON_ADIPOGENESIS_PEAK_AT_24HR
Variation Smoothness All Timeframes  Kang Doxorubicin Resistance Dn                              -0.69009   0.00299341   -1.8686   0.168147   1                9232   0.555556  KANG_DOXORUBICIN_RESISTANCE_DN
Mean Smoothness All Timeframes       Kang Doxorubicin Resistance Dn                              -0.630532  0.014568     -1.7051   0.164501   1                8285   0.611111  KANG_DOXORUBICIN_RESISTANCE_DN
Variation Smoothness All Timeframes  Mclachlan Dental Caries Dn                                  -0.573108  0.00261517   -1.85715  0.168147   1                7460   0.689655  MCLACHLAN_DENTAL_CARIES_DN
Mean Smoothness All Timeframes       Mclachlan Dental Caries Dn                                  -0.601893  0.00120797   -1.95092  0.109401   1                8493   0.603448  MCLACHLAN_DENTAL_CARIES_DN
Mean Smoothness All Timeframes       Tseng Irs1 Targets Dn                                       -0.469643  0.00621367   -1.77165  0.152227   1                8749   0.405941  TSENG_IRS1_TARGETS_DN
Mean Smoothness All Timeframes       Rodwell Aging Kidney No Blood Up                            -0.461268  0.0734188    -1.55491  0.22495    1                8839   0.447368  RODWELL_AGING_KIDNEY_NO_BLOOD_UP
Mean Smoothness All Timeframes       Chang Immortalized By Hpv31 Up                              -0.42639   0.0378389    -1.5254   0.239156   1                7509   0.559322  CHANG_IMMORTALIZED_BY_HPV31_UP
Variation Smoothness Uptake          Mody Hippocampus Prenatal                                    0.630316  0.0086785     1.80944  0.174539   1                1879   0.631579  MODY_HIPPOCAMPUS_PRENATAL
Top Late Enhancement                 Burton Adipogenesis 4                                       -0.554694  0.0058752    -1.79088  0.213855   1                8684   0.613636  BURTON_ADIPOGENESIS_4
Mean Smoothness All Timeframes       Zhu Cmv All Dn                                              -0.529878  0.0518563    -1.66032  0.181553   1                9321   0.482759  ZHU_CMV_ALL_DN
Mean Smoothness All Timeframes       Chiba Response To Tsa Dn                                    -0.598316  0.0313131    -1.62427  0.196797   1                9301   0.391304  CHIBA_RESPONSE_TO_TSA_DN
Irregularity                         Abe Vegfa Targets                                            0.73235   0.00478851    1.83306  0.218099   1                2034   0.625     ABE_VEGFA_TARGETS
Variation Smoothness All Timeframes  Liu Smarca4 Targets                                         -0.516823  0.0380857    -1.60188  0.247184   1                8674   0.456522  LIU_SMARCA4_TARGETS
Mean Smoothness All Timeframes       Liu Smarca4 Targets                                         -0.49642   0.0600278    -1.53883  0.231366   1                9256   0.391304  LIU_SMARCA4_TARGETS
Variation Smoothness All Timeframes  Zhang Antiviral Response To Ribavirin Dn                    -0.51887   0.0183523    -1.67821  0.233816   1                8440   0.571429  ZHANG_ANTIVIRAL_RESPONSE_TO_RIBAVIRIN_DN
Mean Smoothness All Timeframes       Burton Adipogenesis Peak At 0Hr                             -0.522114  0.0485903    -1.62834  0.194532   1                9246   0.446429  BURTON_ADIPOGENESIS_PEAK_AT_0HR
Mean Smoothness All Timeframes       Nielsen Gist And Synovial Sarcoma Up                        -0.600957  0.0614919    -1.50916  0.249559   1                9989   0.454545  NIELSEN_GIST_AND_SYNOVIAL_SARCOMA_UP
Mean Smoothness All Timeframes       Burton Adipogenesis 11                                      -0.483183  0.019798     -1.64052  0.194025   1                8593   0.436364  BURTON_ADIPOGENESIS_11
Variation Smoothness Uptake          Dazard Uv Response Cluster G1                                0.565973  0.0188125     1.75731  0.192501   1                1510   0.528302  DAZARD_UV_RESPONSE_CLUSTER_G1
Top Late Enhancement                 Ly Aging Old Dn                                             -0.681173  0.0196794    -1.74769  0.213855   1                9974   0.574074  LY_AGING_OLD_DN
Mean Smoothness All Timeframes       Tseng Adipogenic Potential Dn                               -0.603128  0.0100705    -1.79124  0.144516   1                8945   0.515152  TSENG_ADIPOGENIC_POTENTIAL_DN
Variation Smoothness All Timeframes  Nielsen Leiomyosarcoma Cnn1 Up                              -0.811856  0.0148344    -1.71031  0.213798   1                9313   0.857143  NIELSEN_LEIOMYOSARCOMA_CNN1_UP
Mean Smoothness All Timeframes       Nielsen Leiomyosarcoma Cnn1 Up                              -0.843385  0.00519584   -1.79564  0.141179   1                9809   0.857143  NIELSEN_LEIOMYOSARCOMA_CNN1_UP
Washout                              Simbulan Parp1 Targets Dn                                    0.812785  0.00792361    1.74442  0.24889    1                1259   0.785714  SIMBULAN_PARP1_TARGETS_DN
Circularity                          Li Adipogenesis By Activated Pparg                          -0.807968  0.0007857    -2.00329  0.205876   1               10368   0.571429  LI_ADIPOGENESIS_BY_ACTIVATED_PPARG
Top Late Enhancement                 Ly Aging Middle Dn                                          -0.909251  0.00543478   -1.70312  0.226338   1               10360   0.866667  LY_AGING_MIDDLE_DN
Variation Smoothness All Timeframes  Bacolod Resistance To Alkylating Agents Up                  -0.688566  0.00242326   -1.86425  0.168147   1                9709   0.555556  BACOLOD_RESISTANCE_TO_ALKYLATING_AGENTS_UP
Mean Smoothness All Timeframes       Bacolod Resistance To Alkylating Agents Up                  -0.706297  0.00122499   -1.90437  0.117961   1                8952   0.666667  BACOLOD_RESISTANCE_TO_ALKYLATING_AGENTS_UP
Top Late Enhancement                 Hu Genotoxic Damage 4Hr                                     -0.659771  0.0247038    -1.70794  0.226338   1                9946   0.580645  HU_GENOTOXIC_DAMAGE_4HR
Mean Smoothness All Timeframes       Mcdowell Acute Lung Injury Dn                               -0.490262  0.0389066    -1.57674  0.216178   1                8982   0.529412  MCDOWELL_ACUTE_LUNG_INJURY_DN
Variation Smoothness All Timeframes  Gentile Uv Low Dose Dn                                      -0.525083  0.016559     -1.70371  0.217973   1                8637   0.491228  GENTILE_UV_LOW_DOSE_DN
Mean Smoothness All Timeframes       Gentile Uv Low Dose Dn                                      -0.561837  0.00523349   -1.82213  0.133542   1                8369   0.596491  GENTILE_UV_LOW_DOSE_DN
Variation Smoothness All Timeframes  Wang Smarce1 Targets Up                                     -0.553012  0.0429933    -1.69923  0.219274   1                8254   0.572687  WANG_SMARCE1_TARGETS_UP
Mean Smoothness All Timeframes       Wang Smarce1 Targets Up                                     -0.587553  0.017803     -1.80752  0.141179   1                8910   0.54185   WANG_SMARCE1_TARGETS_UP
Mean Smoothness All Timeframes       Lee Aging Muscle Up                                         -0.475862  0.0192503    -1.64464  0.192647   1                8803   0.485714  LEE_AGING_MUSCLE_UP
Washout                              Visala Response To Heat Shock And Aging Dn                   0.729751  0.00525996    1.81308  0.206704   1                1952   0.769231  VISALA_RESPONSE_TO_HEAT_SHOCK_AND_AGING_DN
Top Late Enhancement                 Visala Response To Heat Shock And Aging Dn                  -0.702552  0.0107244    -1.73664  0.213855   1               10151   0.538462  VISALA_RESPONSE_TO_HEAT_SHOCK_AND_AGING_DN
Vol Late Lt0                         Visala Response To Heat Shock And Aging Dn                   0.748903  0.00304507    1.86069  0.214618   1                2110   0.846154  VISALA_RESPONSE_TO_HEAT_SHOCK_AND_AGING_DN
Mean Smoothness All Timeframes       Browne Hcmv Infection 12Hr Dn                               -0.453601  0.0180254    -1.63773  0.194228   1                8408   0.453488  BROWNE_HCMV_INFECTION_12HR_DN
Washout                              Jiang Aging Hypothalamus Up                                  0.597785  0.00503423    1.84632  0.193047   1                1471   0.547619  JIANG_AGING_HYPOTHALAMUS_UP
Vol Late Lt0                         Jiang Aging Hypothalamus Up                                  0.600201  0.00518031    1.85973  0.214618   1                1406   0.547619  JIANG_AGING_HYPOTHALAMUS_UP
Variation Smoothness Uptake          Jiang Aging Hypothalamus Up                                  0.594955  0.00569409    1.83259  0.173515   1                1666   0.571429  JIANG_AGING_HYPOTHALAMUS_UP
Mean Smoothness All Timeframes       Su Placenta                                                 -0.538487  0.0547836    -1.52184  0.242693   1                8292   0.5       SU_PLACENTA
Washout                              Zhong Secretome Of Lung Cancer And Macrophage                0.546595  0.00790113    1.76702  0.237566   1                1777   0.492537  ZHONG_SECRETOME_OF_LUNG_CANCER_AND_MACROPHAGE
Top Late Enhancement                 Zhong Secretome Of Lung Cancer And Macrophage               -0.565722  0.00383065   -1.83074  0.213855   1                9567   0.447761  ZHONG_SECRETOME_OF_LUNG_CANCER_AND_MACROPHAGE
Washout                              Zhong Secretome Of Lung Cancer And Endothelium               0.550598  0.00285249    1.85119  0.193047   1                1777   0.5       ZHONG_SECRETOME_OF_LUNG_CANCER_AND_ENDOTHELIUM
Top Late Enhancement                 Zhong Secretome Of Lung Cancer And Endothelium              -0.532798  0.00727126   -1.79672  0.213855   1                9416   0.45      ZHONG_SECRETOME_OF_LUNG_CANCER_AND_ENDOTHELIUM
Washout                              Zhong Secretome Of Lung Cancer And Fibroblast                0.525551  0.00342259    1.84274  0.193047   1                2227   0.5       ZHONG_SECRETOME_OF_LUNG_CANCER_AND_FIBROBLAST
Top Late Enhancement                 Zhong Secretome Of Lung Cancer And Fibroblast               -0.535108  0.00220928   -1.88186  0.213855   1                9393   0.432203  ZHONG_SECRETOME_OF_LUNG_CANCER_AND_FIBROBLAST
Mean Smoothness All Timeframes       Clasper Lymphatic Vessels During Metastasis Dn              -0.708924  0.101408     -1.5477   0.227442   1                9466   0.7       CLASPER_LYMPHATIC_VESSELS_DURING_METASTASIS_DN
Top Late Enhancement                 Singh Nfe2L2 Targets                                        -0.710761  0.0138015    -1.72128  0.219645   1                9520   0.642857  SINGH_NFE2L2_TARGETS
Ser                                  Pellicciotta Hdac In Antigen Presentation Up                 0.629409  0.00244499    1.96323  0.21078    1                1292   0.442623  PELLICCIOTTA_HDAC_IN_ANTIGEN_PRESENTATION_UP
Top Late Enhancement                 Pellicciotta Hdac In Antigen Presentation Up                -0.555965  0.0231715    -1.72696  0.217384   1                9514   0.459016  PELLICCIOTTA_HDAC_IN_ANTIGEN_PRESENTATION_UP
Variation Smoothness Uptake          Pellicciotta Hdac In Antigen Presentation Dn                 0.631989  0.023569      1.71853  0.239548   1                2115   0.510638  PELLICCIOTTA_HDAC_IN_ANTIGEN_PRESENTATION_DN
Mean Smoothness All Timeframes       Harris Brain Cancer Progenitors                             -0.620688  0.0578926    -1.58012  0.216167   1                9198   0.55      HARRIS_BRAIN_CANCER_PROGENITORS
Variation Smoothness All Timeframes  Lein Neuron Markers                                         -0.540431  0.00809061   -1.66148  0.236572   1                8733   0.4       LEIN_NEURON_MARKERS
Irregularity                         Lein Oligodendrocyte Markers                                 0.550856  0.000605694   1.98095  0.12903    1                1185   0.27451   LEIN_OLIGODENDROCYTE_MARKERS
Variation Smoothness All Timeframes  Lein Choroid Plexus Markers                                 -0.523168  0.00601443   -1.81374  0.191445   1                8284   0.537313  LEIN_CHOROID_PLEXUS_MARKERS
Mean Smoothness All Timeframes       Lein Choroid Plexus Markers                                 -0.526036  0.0038291    -1.8203   0.133542   1                8953   0.462687  LEIN_CHOROID_PLEXUS_MARKERS
Irregularity                         Lein Cerebellum Markers                                      0.538654  0.00178962    1.89786  0.182442   1                2029   0.408163  LEIN_CEREBELLUM_MARKERS
Variation Smoothness All Timeframes  Lein Cerebellum Markers                                     -0.476855  0.0138166    -1.68198  0.231284   1                8558   0.44898   LEIN_CEREBELLUM_MARKERS
Mean Smoothness All Timeframes       Lein Cerebellum Markers                                     -0.478188  0.0124849    -1.68606  0.172248   1                9804   0.306122  LEIN_CEREBELLUM_MARKERS
Variation Smoothness All Timeframes  Lein Localized To Distal And Proximal Dendrites             -0.597721  0.0259171    -1.60394  0.247184   1                8618   0.5       LEIN_LOCALIZED_TO_DISTAL_AND_PROXIMAL_DENDRITES
Mean Smoothness All Timeframes       Lein Localized To Distal And Proximal Dendrites             -0.595598  0.0246129    -1.60281  0.204175   1                8975   0.5       LEIN_LOCALIZED_TO_DISTAL_AND_PROXIMAL_DENDRITES
Variation Smoothness All Timeframes  Gavin Foxp3 Targets Cluster T7                               0.518335  0.00275536    1.84672  0.212371   1                1900   0.422222  GAVIN_FOXP3_TARGETS_CLUSTER_T7
Variation Smoothness Uptake          Gavin Foxp3 Targets Cluster T7                               0.542903  0.00058812    1.93823  0.173515   1                2241   0.488889  GAVIN_FOXP3_TARGETS_CLUSTER_T7
Mean Smoothness All Timeframes       Gavin Foxp3 Targets Cluster P7                              -0.467492  0.0498709    -1.53713  0.231366   1                8148   0.5       GAVIN_FOXP3_TARGETS_CLUSTER_P7
Mean Smoothness All Timeframes       Riggins Tamoxifen Resistance Dn                             -0.390529  0.0301265    -1.55095  0.226989   1                7048   0.524272  RIGGINS_TAMOXIFEN_RESISTANCE_DN
Mean Smoothness All Timeframes       Wang Lsd1 Targets Up                                        -0.581235  0.0737206    -1.51088  0.248827   1                8870   0.470588  WANG_LSD1_TARGETS_UP
Mean Smoothness All Timeframes       Ji Carcinogenesis By Kras And Stk11 Dn                      -0.732257  0.0639806    -1.57391  0.216793   1                9634   0.625     JI_CARCINOGENESIS_BY_KRAS_AND_STK11_DN
Mean Smoothness All Timeframes       Kyng Dna Damage By 4Nqo                                     -0.467562  0.0263631    -1.58932  0.212416   1                9690   0.4       KYNG_DNA_DAMAGE_BY_4NQO
Mean Smoothness All Timeframes       Kondo Prostate Cancer With H3K27Me3                         -0.518181  0.061154     -1.51232  0.248519   1                9080   0.384615  KONDO_PROSTATE_CANCER_WITH_H3K27ME3
Variation Smoothness All Timeframes  Kondo Ezh2 Targets                                          -0.511966  0.00139581   -1.89863  0.168147   1                8098   0.532164  KONDO_EZH2_TARGETS
Mean Smoothness All Timeframes       Kondo Ezh2 Targets                                          -0.528297  0.00100604   -1.94813  0.109401   1                8639   0.48538   KONDO_EZH2_TARGETS
Irregularity                         Lee Liver Cancer                                             0.697049  0.00119928    1.90446  0.176552   1                3006   0.8       LEE_LIVER_CANCER
Variation Smoothness All Timeframes  Wong Mitochondria Gene Module                                0.60301   0.0134787     1.8148   0.237345   1                2257   0.592965  WONG_MITOCHONDRIA_GENE_MODULE
Variation Smoothness Uptake          Wong Mitochondria Gene Module                                0.598795  0.0138861     1.80411  0.174539   1                1857   0.547739  WONG_MITOCHONDRIA_GENE_MODULE
Washout                              Wong Proteasome Gene Module                                  0.598725  0.0166566     1.74867  0.247737   1                1310   0.478261  WONG_PROTEASOME_GENE_MODULE
Top Late Enhancement                 Wong Proteasome Gene Module                                 -0.596154  0.0157874    -1.74145  0.213855   1               10072   0.456522  WONG_PROTEASOME_GENE_MODULE
Variation Smoothness Uptake          Wong Proteasome Gene Module                                  0.588698  0.0225074     1.71973  0.239548   1                2112   0.543478  WONG_PROTEASOME_GENE_MODULE
Top Late Enhancement                 Amundson Gamma Radiation Response                           -0.849801  0.00586688   -1.78636  0.213855   1               10154   0.823529  AMUNDSON_GAMMA_RADIATION_RESPONSE
Mean Smoothness All Timeframes       Massarweh Response To Estradiol                             -0.480321  0.0588351    -1.55651  0.223946   1                7972   0.5       MASSARWEH_RESPONSE_TO_ESTRADIOL
Washout                              Sarrio Epithelial Mesenchymal Transition Up                  0.630999  0.0327435     1.75217  0.246582   1                2102   0.625     SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_UP
Top Late Enhancement                 Sarrio Epithelial Mesenchymal Transition Up                 -0.670725  0.0122274    -1.85428  0.213855   1                9992   0.544118  SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_UP
Mean Smoothness All Timeframes       Ouillette Cll 13Q14 Deletion Dn                             -0.435254  0.032509     -1.55087  0.226989   1                9347   0.311111  OUILLETTE_CLL_13Q14_DELETION_DN
Variation Smoothness All Timeframes  Iwanaga Carcinogenesis By Kras Pten Dn                      -0.41766   0.00217994   -1.76504  0.193127   1                9418   0.305556  IWANAGA_CARCINOGENESIS_BY_KRAS_PTEN_DN
Mean Smoothness All Timeframes       Iwanaga Carcinogenesis By Kras Pten Dn                      -0.416413  0.00321156   -1.75363  0.153387   1                9146   0.333333  IWANAGA_CARCINOGENESIS_BY_KRAS_PTEN_DN
Variation Smoothness All Timeframes  Iwanaga Carcinogenesis By Kras Dn                           -0.407832  0.00639872   -1.63801  0.236712   1                7957   0.448276  IWANAGA_CARCINOGENESIS_BY_KRAS_DN
Mean Smoothness All Timeframes       Camps Colon Cancer Copy Number Dn                           -0.533914  0.0532319    -1.53269  0.233789   1                8894   0.388889  CAMPS_COLON_CANCER_COPY_NUMBER_DN
Variation Smoothness All Timeframes  Dairkee Cancer Prone Response Bpa E2                         0.545861  0.00118437    1.90878  0.206872   1                2242   0.483516  DAIRKEE_CANCER_PRONE_RESPONSE_BPA_E2
Variation Smoothness Uptake          Dairkee Cancer Prone Response Bpa E2                         0.528004  0.00256309    1.85154  0.173515   1                2064   0.450549  DAIRKEE_CANCER_PRONE_RESPONSE_BPA_E2
Mean Smoothness All Timeframes       Riggi Ewing Sarcoma Progenitor Up                           -0.404833  0.0252016    -1.58004  0.216167   1                7384   0.566787  RIGGI_EWING_SARCOMA_PROGENITOR_UP
Variation Smoothness All Timeframes  Riggi Ewing Sarcoma Progenitor Dn                           -0.538464  0.0347041    -1.72658  0.211541   1                9003   0.5       RIGGI_EWING_SARCOMA_PROGENITOR_DN
Mean Smoothness All Timeframes       Riggi Ewing Sarcoma Progenitor Dn                           -0.577683  0.0136848    -1.84921  0.130786   1                8652   0.608696  RIGGI_EWING_SARCOMA_PROGENITOR_DN
Variation Smoothness All Timeframes  Ambrosini Flavopiridol Treatment Tp53                       -0.396755  0.00674068   -1.62879  0.24086    1                7198   0.477273  AMBROSINI_FLAVOPIRIDOL_TREATMENT_TP53
Irregularity                         Huang Foxa2 Targets Up                                       0.607199  0.00060241    1.93877  0.12903    1                1486   0.405405  HUANG_FOXA2_TARGETS_UP
Variation Smoothness All Timeframes  Huang Foxa2 Targets Up                                      -0.520626  0.0148565    -1.65292  0.236572   1                9184   0.459459  HUANG_FOXA2_TARGETS_UP
Mean Smoothness All Timeframes       Masri Resistance To Tamoxifen And Aromatase Inhibitors Up   -0.619693  0.0333001    -1.57107  0.218767   1                9909   0.444444  MASRI_RESISTANCE_TO_TAMOXIFEN_AND_AROMATASE_INHIBITORS_UP
Mean Smoothness All Timeframes       Mccabe Hoxc6 Targets Dn                                     -0.594756  0.0528356    -1.5514   0.226989   1                9735   0.466667  MCCABE_HOXC6_TARGETS_DN
Variation Smoothness All Timeframes  Ting Silenced By Dicer                                      -0.578022  0.0295115    -1.63882  0.236712   1                8493   0.56      TING_SILENCED_BY_DICER
Mean Smoothness All Timeframes       Ting Silenced By Dicer                                      -0.567694  0.0330264    -1.60469  0.204175   1                9577   0.44      TING_SILENCED_BY_DICER
Variation Smoothness All Timeframes  Molenaar Targets Of Ccnd1 And Cdk4 Up                       -0.463136  0.0212892    -1.60204  0.247184   1                8431   0.4375    MOLENAAR_TARGETS_OF_CCND1_AND_CDK4_UP
Mean Smoothness All Timeframes       Molenaar Targets Of Ccnd1 And Cdk4 Up                       -0.451465  0.0333399    -1.55951  0.223217   1                8195   0.458333  MOLENAAR_TARGETS_OF_CCND1_AND_CDK4_UP
Variation Smoothness All Timeframes  Acevedo Liver Cancer With H3K27Me3 Up                       -0.448498  0.00604473   -1.72317  0.212304   1                7891   0.511628  ACEVEDO_LIVER_CANCER_WITH_H3K27ME3_UP
Mean Smoothness All Timeframes       Acevedo Liver Cancer With H3K27Me3 Up                       -0.469574  0.00160901   -1.80315  0.141179   1                8417   0.48062   ACEVEDO_LIVER_CANCER_WITH_H3K27ME3_UP
Mean Smoothness All Timeframes       Acevedo Liver Cancer With H3K27Me3 Dn                       -0.394511  0.0140647    -1.56772  0.21886    1                8499   0.382716  ACEVEDO_LIVER_CANCER_WITH_H3K27ME3_DN
Mean Smoothness All Timeframes       Hoque Methylated In Cancer                                  -0.492815  0.0214558    -1.61771  0.200977   1                9094   0.416667  HOQUE_METHYLATED_IN_CANCER
Variation Smoothness All Timeframes  Smid Breast Cancer Relapse In Lung Dn                       -0.75076   0.00457438   -1.7909   0.193127   1               10130   0.533333  SMID_BREAST_CANCER_RELAPSE_IN_LUNG_DN
Mean Smoothness All Timeframes       Smid Breast Cancer Relapse In Lung Dn                       -0.735273  0.0101231    -1.76055  0.153387   1               10456   0.5       SMID_BREAST_CANCER_RELAPSE_IN_LUNG_DN
Variation Smoothness All Timeframes  Smid Breast Cancer Luminal A Up                             -0.715306  0.00654892   -1.84692  0.171536   1                9467   0.643836  SMID_BREAST_CANCER_LUMINAL_A_UP
Mean Smoothness All Timeframes       Smid Breast Cancer Luminal A Up                             -0.734388  0.00357214   -1.89821  0.117961   1               10178   0.575342  SMID_BREAST_CANCER_LUMINAL_A_UP
Mean Smoothness All Timeframes       Le Ski Targets Up                                           -0.63522   0.0434695    -1.57003  0.218767   1                9039   0.571429  LE_SKI_TARGETS_UP
Variation Smoothness All Timeframes  Bonome Ovarian Cancer Poor Survival Up                      -0.592096  0.0219669    -1.6353   0.236712   1                7854   0.64      BONOME_OVARIAN_CANCER_POOR_SURVIVAL_UP
Mean Smoothness All Timeframes       Bonome Ovarian Cancer Poor Survival Up                      -0.657676  0.00360072   -1.82132  0.133542   1                8528   0.64      BONOME_OVARIAN_CANCER_POOR_SURVIVAL_UP
Variation Smoothness All Timeframes  Izadpanah Stem Cell Adipose Vs Bone Dn                      -0.626926  0.00199641   -1.97168  0.168147   1                8595   0.597826  IZADPANAH_STEM_CELL_ADIPOSE_VS_BONE_DN
Mean Smoothness All Timeframes       Izadpanah Stem Cell Adipose Vs Bone Dn                      -0.62474   0.00281237   -1.96379  0.109401   1                8348   0.641304  IZADPANAH_STEM_CELL_ADIPOSE_VS_BONE_DN
Variation Smoothness All Timeframes  Zheng Glioblastoma Plasticity Dn                            -0.55696   0.0339932    -1.64981  0.236572   1                8728   0.477273  ZHENG_GLIOBLASTOMA_PLASTICITY_DN
Mean Smoothness All Timeframes       Zheng Glioblastoma Plasticity Dn                            -0.60528   0.00980981   -1.78626  0.148577   1                8380   0.568182  ZHENG_GLIOBLASTOMA_PLASTICITY_DN
Irregularity                         Buckanovich T Lymphocyte Homing On Tumor Dn                  0.604277  0.00180072    1.90366  0.176552   1                1596   0.45      BUCKANOVICH_T_LYMPHOCYTE_HOMING_ON_TUMOR_DN
Variation Smoothness All Timeframes  Buckanovich T Lymphocyte Homing On Tumor Dn                 -0.634447  0.000601805  -2.01196  0.147564   1                9966   0.4       BUCKANOVICH_T_LYMPHOCYTE_HOMING_ON_TUMOR_DN
Mean Smoothness All Timeframes       Buckanovich T Lymphocyte Homing On Tumor Dn                 -0.583371  0.00279776   -1.85122  0.130455   1                9512   0.4       BUCKANOVICH_T_LYMPHOCYTE_HOMING_ON_TUMOR_DN
Mean Smoothness All Timeframes       Zheng Il22 Signaling Dn                                     -0.491648  0.0455903    -1.516    0.246789   1                8975   0.454545  ZHENG_IL22_SIGNALING_DN
Washout                              Honma Docetaxel Resistance                                   0.702655  0.00343365    1.84224  0.193047   1                2667   0.764706  HONMA_DOCETAXEL_RESISTANCE
Top Late Enhancement                 Honma Docetaxel Resistance                                  -0.701898  0.00525359   -1.82968  0.213855   1                8645   0.764706  HONMA_DOCETAXEL_RESISTANCE
Variation Smoothness All Timeframes  Matzuk Central For Female Fertility                         -0.645     0.0328843    -1.63848  0.236712   1               10568   0.5       MATZUK_CENTRAL_FOR_FEMALE_FERTILITY
Mean Smoothness All Timeframes       Matzuk Central For Female Fertility                         -0.663554  0.0230984    -1.67807  0.175463   1               10678   0.5       MATZUK_CENTRAL_FOR_FEMALE_FERTILITY
Top Late Enhancement                 Whiteford Pediatric Cancer Markers                          -0.779092  0.0159151    -1.75364  0.213855   1                9973   0.733333  WHITEFORD_PEDIATRIC_CANCER_MARKERS
Washout                              Grade Metastasis Dn                                          0.667737  0.00442567    1.84398  0.193047   1                2029   0.609756  GRADE_METASTASIS_DN
Top Late Enhancement                 Grade Metastasis Dn                                         -0.677631  0.00380685   -1.87298  0.213855   1                9544   0.560976  GRADE_METASTASIS_DN
Vol Late Lt0                         Grade Metastasis Dn                                          0.677569  0.00278884    1.86974  0.214618   1                2072   0.658537  GRADE_METASTASIS_DN
Variation Smoothness All Timeframes  Grade Metastasis Dn                                          0.666657  0.00357072    1.84068  0.212371   1                1830   0.560976  GRADE_METASTASIS_DN
Variation Smoothness Uptake          Grade Metastasis Dn                                          0.726656  0.000198689   2.00696  0.173515   0.529108         1447   0.585366  GRADE_METASTASIS_DN
Variation Smoothness All Timeframes  Grade Colon Cancer Dn                                       -0.492356  0.0197723    -1.61347  0.244212   1               10095   0.296296  GRADE_COLON_CANCER_DN
Washout                              Grade Colon And Rectal Cancer Up                             0.560133  0.00425273    1.86707  0.193047   1                2204   0.509728  GRADE_COLON_AND_RECTAL_CANCER_UP
Top Late Enhancement                 Grade Colon And Rectal Cancer Up                            -0.585738  0.00120724   -1.95014  0.213855   1                9147   0.525292  GRADE_COLON_AND_RECTAL_CANCER_UP
Variation Smoothness All Timeframes  Grade Colon And Rectal Cancer Up                             0.551072  0.00514037    1.84359  0.212371   1                2250   0.544747  GRADE_COLON_AND_RECTAL_CANCER_UP
Variation Smoothness Uptake          Grade Colon And Rectal Cancer Up                             0.551792  0.00614713    1.84258  0.173515   1                2231   0.536965  GRADE_COLON_AND_RECTAL_CANCER_UP
Variation Smoothness All Timeframes  Grade Colon And Rectal Cancer Dn                            -0.495867  0.00416749   -1.74086  0.198874   1                8780   0.421875  GRADE_COLON_AND_RECTAL_CANCER_DN
Mean Smoothness All Timeframes       Grade Colon And Rectal Cancer Dn                            -0.516093  0.00181232   -1.80067  0.141179   1                9314   0.375     GRADE_COLON_AND_RECTAL_CANCER_DN
Variation Smoothness All Timeframes  Boquest Stem Cell Up                                        -0.643486  0.0360306    -1.75055  0.193127   1                8237   0.727273  BOQUEST_STEM_CELL_UP
Mean Smoothness All Timeframes       Boquest Stem Cell Up                                        -0.687919  0.00960384   -1.88045  0.117961   1                9659   0.61244   BOQUEST_STEM_CELL_UP
Ld Late Lt0                          Boquest Stem Cell Dn                                         0.597005  0.00277118    1.98654  0.180214   1                1486   0.374269  BOQUEST_STEM_CELL_DN
Variation Smoothness All Timeframes  Boquest Stem Cell Cultured Vs Fresh Dn                      -0.65461   0.0547314    -1.61394  0.244212   1               10096   0.5       BOQUEST_STEM_CELL_CULTURED_VS_FRESH_DN
Mean Smoothness All Timeframes       Boquest Stem Cell Cultured Vs Fresh Dn                      -0.6807    0.0315895    -1.6794   0.175441   1               10336   0.5       BOQUEST_STEM_CELL_CULTURED_VS_FRESH_DN
Irregularity                         Labbe Wnt3A Targets Dn                                       0.495686  0.00321543    1.80805  0.243479   1                2474   0.377358  LABBE_WNT3A_TARGETS_DN
Variation Smoothness All Timeframes  Labbe Wnt3A Targets Dn                                      -0.445312  0.016805     -1.62343  0.244212   1                9018   0.320755  LABBE_WNT3A_TARGETS_DN
Mean Smoothness All Timeframes       Labbe Tgfb1 Targets Up                                      -0.47156   0.0663602    -1.53016  0.235273   1                9012   0.402985  LABBE_TGFB1_TARGETS_UP
Variation Smoothness All Timeframes  Labbe Targets Of Tgfb1 And Wnt3A Up                         -0.484578  0.0357001    -1.60555  0.247184   1                9245   0.351351  LABBE_TARGETS_OF_TGFB1_AND_WNT3A_UP
Mean Smoothness All Timeframes       Labbe Targets Of Tgfb1 And Wnt3A Up                         -0.514224  0.0183579    -1.69423  0.169834   1                8756   0.432432  LABBE_TARGETS_OF_TGFB1_AND_WNT3A_UP
Washout                              West Adrenocortical Tumor Up                                 0.522665  0.00824618    1.82885  0.206704   1                1840   0.450758  WEST_ADRENOCORTICAL_TUMOR_UP
Top Late Enhancement                 West Adrenocortical Tumor Up                                -0.524518  0.00625      -1.82269  0.213855   1                9216   0.465909  WEST_ADRENOCORTICAL_TUMOR_UP
Top Late Enhancement                 West Adrenocortical Carcinoma Vs Adenoma Up                 -0.69259   0.00705219   -1.75291  0.213855   1                9640   0.583333  WEST_ADRENOCORTICAL_CARCINOMA_VS_ADENOMA_UP
Variation Smoothness All Timeframes  West Adrenocortical Carcinoma Vs Adenoma Up                  0.757095  0.00221729    1.92121  0.206872   1                1245   0.666667  WEST_ADRENOCORTICAL_CARCINOMA_VS_ADENOMA_UP
Variation Smoothness Uptake          West Adrenocortical Carcinoma Vs Adenoma Up                  0.681498  0.0108086     1.72806  0.234559   1                1181   0.583333  WEST_ADRENOCORTICAL_CARCINOMA_VS_ADENOMA_UP
Irregularity                         Ray Tumorigenesis By Erbb2 Cdc25A Dn                         0.474707  0.0012053     1.85995  0.218099   1                2088   0.357143  RAY_TUMORIGENESIS_BY_ERBB2_CDC25A_DN
Variation Smoothness All Timeframes  Ray Tumorigenesis By Erbb2 Cdc25A Dn                        -0.455811  0.00262361   -1.79319  0.193127   1                7570   0.5       RAY_TUMORIGENESIS_BY_ERBB2_CDC25A_DN
Mean Smoothness All Timeframes       Ray Tumorigenesis By Erbb2 Cdc25A Dn                        -0.437377  0.00663183   -1.71802  0.16069    1                8068   0.428571  RAY_TUMORIGENESIS_BY_ERBB2_CDC25A_DN
Mean Smoothness All Timeframes       Yague Pretumor Drug Resistance Dn                           -0.818372  0.0346       -1.56517  0.219911   1                9596   0.833333  YAGUE_PRETUMOR_DRUG_RESISTANCE_DN
Washout                              Podar Response To Adaphostin Dn                              0.713945  0.00964049    1.75682  0.242201   1                1919   0.666667  PODAR_RESPONSE_TO_ADAPHOSTIN_DN
Vol Late Lt0                         Podar Response To Adaphostin Dn                              0.769208  0.00160901    1.88707  0.214618   1                1282   0.666667  PODAR_RESPONSE_TO_ADAPHOSTIN_DN
Top Late Enhancement                 Chen Hoxa5 Targets 9Hr Dn                                   -0.502245  0.0143319    -1.68042  0.230757   1                9595   0.4375    CHEN_HOXA5_TARGETS_9HR_DN
Circularity                          Cadwell Atg16L1 Targets Up                                  -0.604823  0.00080032   -1.99662  0.205876   1               10555   0.319149  CADWELL_ATG16L1_TARGETS_UP
Irregularity                         Cadwell Atg16L1 Targets Up                                   0.592599  0.00241449    1.94794  0.12903    1                 840   0.319149  CADWELL_ATG16L1_TARGETS_UP
Vol Init Enhancement Gt100           Cadwell Atg16L1 Targets Up                                  -0.648389  0.00020016   -2.11838  0.139336   0.533026         9360   0.489362  CADWELL_ATG16L1_TARGETS_UP
Variation Smoothness All Timeframes  Cadwell Atg16L1 Targets Up                                  -0.560463  0.00403796   -1.84295  0.171536   1                8869   0.446809  CADWELL_ATG16L1_TARGETS_UP
Mean Smoothness All Timeframes       Cadwell Atg16L1 Targets Up                                  -0.514656  0.0160243    -1.69501  0.169834   1                8884   0.425532  CADWELL_ATG16L1_TARGETS_UP
Variation Smoothness All Timeframes  Bredemeyer Rag Signaling Not Via Atm Up                     -0.441367  0.00928355   -1.61696  0.244212   1                7010   0.581395  BREDEMEYER_RAG_SIGNALING_NOT_VIA_ATM_UP
Mean Smoothness All Timeframes       Bredemeyer Rag Signaling Not Via Atm Up                     -0.462003  0.00421433   -1.69221  0.169834   1                9243   0.348837  BREDEMEYER_RAG_SIGNALING_NOT_VIA_ATM_UP
Mean Smoothness All Timeframes       Thum Mir21 Targets Heart Disease Up                         -0.795919  0.0662772    -1.52196  0.242693   1                9659   0.866667  THUM_MIR21_TARGETS_HEART_DISEASE_UP
Variation Smoothness All Timeframes  Chng Multiple Myeloma Hyperploid Up                          0.792     0.0110041     1.85665  0.212371   1                1651   0.86      CHNG_MULTIPLE_MYELOMA_HYPERPLOID_UP
Variation Smoothness Uptake          Chng Multiple Myeloma Hyperploid Up                          0.790179  0.0122288     1.85475  0.173515   1                1452   0.84      CHNG_MULTIPLE_MYELOMA_HYPERPLOID_UP
Washout                              Maloney Response To 17Aag Dn                                 0.585472  0.0140732     1.76216  0.23934    1                1932   0.554054  MALONEY_RESPONSE_TO_17AAG_DN
Variation Smoothness Uptake          Maloney Response To 17Aag Dn                                 0.605482  0.00734564    1.81222  0.174539   1                1991   0.581081  MALONEY_RESPONSE_TO_17AAG_DN
Variation Smoothness All Timeframes  Mueller Plurinet                                             0.603169  0.0031552     1.90367  0.206872   1                1959   0.547325  MUELLER_PLURINET
Variation Smoothness Uptake          Mueller Plurinet                                             0.5655    0.0142772     1.7818   0.187272   1                1774   0.489712  MUELLER_PLURINET
Mean Smoothness All Timeframes       Tsai Response To Radiation Therapy                          -0.573945  0.0733144    -1.53235  0.233789   1                9301   0.482759  TSAI_RESPONSE_TO_RADIATION_THERAPY
Mean Smoothness All Timeframes       Vart Kshv Infection Angiogenic Markers Dn                   -0.493573  0.0429056    -1.60895  0.202854   1                8013   0.542169  VART_KSHV_INFECTION_ANGIOGENIC_MARKERS_DN
Washout                              Lee Liver Cancer Survival Dn                                 0.600188  0.00541625    1.87731  0.193047   1                2125   0.551282  LEE_LIVER_CANCER_SURVIVAL_DN
Top Late Enhancement                 Lee Liver Cancer Survival Dn                                -0.577587  0.0110024    -1.78611  0.213855   1                8869   0.583333  LEE_LIVER_CANCER_SURVIVAL_DN
Variation Smoothness Uptake          Lee Liver Cancer Survival Dn                                 0.557675  0.0215458     1.72295  0.2391     1                2257   0.557692  LEE_LIVER_CANCER_SURVIVAL_DN
Variation Smoothness All Timeframes  Wu Silenced By Methylation In Bladder Cancer                -0.600221  0.0259274    -1.70035  0.219274   1                8964   0.541667  WU_SILENCED_BY_METHYLATION_IN_BLADDER_CANCER
Mean Smoothness All Timeframes       Wu Silenced By Methylation In Bladder Cancer                -0.611754  0.0203187    -1.72915  0.16069    1                8842   0.604167  WU_SILENCED_BY_METHYLATION_IN_BLADDER_CANCER
Mean Smoothness All Timeframes       Gu Pdef Targets Up                                          -0.564875  0.0593389    -1.63101  0.194532   1                9535   0.484375  GU_PDEF_TARGETS_UP
Variation Smoothness All Timeframes  Gu Pdef Targets Dn                                          -0.528768  0.030895     -1.61111  0.244212   1                7377   0.518519  GU_PDEF_TARGETS_DN
Mean Smoothness All Timeframes       Gu Pdef Targets Dn                                          -0.514483  0.0400802    -1.56634  0.21886    1                6350   0.777778  GU_PDEF_TARGETS_DN
Variation Smoothness All Timeframes  Wagner Apo2 Sensitivity                                     -0.71477   0.0234657    -1.64815  0.236712   1                8509   0.818182  WAGNER_APO2_SENSITIVITY
Mean Smoothness All Timeframes       Wagner Apo2 Sensitivity                                     -0.723559  0.0167098    -1.67245  0.178057   1                8709   0.818182  WAGNER_APO2_SENSITIVITY
Variation Smoothness All Timeframes  Seki Inflammatory Response Lps Dn                           -0.651599  0.0222133    -1.66102  0.236572   1                7903   0.785714  SEKI_INFLAMMATORY_RESPONSE_LPS_DN
Mean Smoothness All Timeframes       Seki Inflammatory Response Lps Dn                           -0.630656  0.0298745    -1.61554  0.200977   1                7622   0.857143  SEKI_INFLAMMATORY_RESPONSE_LPS_DN
Variation Smoothness All Timeframes  Frasor Response To Serm Or Fulvestrant Up                   -0.573097  0.0277833    -1.61626  0.244212   1                8200   0.473684  FRASOR_RESPONSE_TO_SERM_OR_FULVESTRANT_UP
Mean Smoothness All Timeframes       Frasor Response To Serm Or Fulvestrant Up                   -0.572927  0.0298358    -1.61261  0.202062   1                7870   0.526316  FRASOR_RESPONSE_TO_SERM_OR_FULVESTRANT_UP
Variation Smoothness All Timeframes  Ichiba Graft Versus Host Disease 35D Dn                     -0.526657  0.0152152    -1.65036  0.236572   1                9338   0.37037   ICHIBA_GRAFT_VERSUS_HOST_DISEASE_35D_DN
Mean Smoothness All Timeframes       Ichiba Graft Versus Host Disease 35D Dn                     -0.517289  0.0184634    -1.62825  0.194532   1                8857   0.444444  ICHIBA_GRAFT_VERSUS_HOST_DISEASE_35D_DN
Variation Smoothness All Timeframes  Ray Targets Of P210 Bcr Abl Fusion Up                       -0.59466   0.0273147    -1.61031  0.244212   1                8230   0.615385  RAY_TARGETS_OF_P210_BCR_ABL_FUSION_UP
Mean Smoothness All Timeframes       Ray Targets Of P210 Bcr Abl Fusion Up                       -0.60155   0.0222668    -1.63239  0.194532   1                9429   0.461538  RAY_TARGETS_OF_P210_BCR_ABL_FUSION_UP
Variation Smoothness All Timeframes  Ruiz Tnc Targets Up                                         -0.507127  0.00806127   -1.78442  0.193127   1                8049   0.546099  RUIZ_TNC_TARGETS_UP
Mean Smoothness All Timeframes       Ruiz Tnc Targets Up                                         -0.498873  0.0109157    -1.75057  0.153387   1                8400   0.496454  RUIZ_TNC_TARGETS_UP
Mean Smoothness All Timeframes       Schraets Mll Targets Up                                     -0.570038  0.0252597    -1.66289  0.181553   1                9965   0.413793  SCHRAETS_MLL_TARGETS_UP
Washout                              Chang Core Serum Response Up                                 0.572482  0.00947772    1.82721  0.206704   1                2382   0.594595  CHANG_CORE_SERUM_RESPONSE_UP
Top Late Enhancement                 Chang Core Serum Response Up                                -0.595631  0.00363123   -1.88785  0.213855   1                9168   0.562162  CHANG_CORE_SERUM_RESPONSE_UP
Variation Smoothness All Timeframes  Chang Core Serum Response Dn                                -0.451707  0.00980981   -1.70916  0.213798   1                8812   0.413265  CHANG_CORE_SERUM_RESPONSE_DN
Mean Smoothness All Timeframes       Chang Core Serum Response Dn                                -0.458168  0.00758937   -1.73277  0.16069    1                9019   0.387755  CHANG_CORE_SERUM_RESPONSE_DN
Top Late Enhancement                 Chang Cycling Genes                                         -0.68276   0.0357578    -1.7349   0.213855   1                9982   0.606838  CHANG_CYCLING_GENES
Uptake Speed                         Aguirre Pancreatic Cancer Copy Number Up                    -0.547524  0.00060423   -2.11769  0.068206   1                8784   0.48855   AGUIRRE_PANCREATIC_CANCER_COPY_NUMBER_UP
Washout                              Hoffmann Immature To Mature B Lymphocyte Dn                  0.523397  0.00527063    1.76704  0.237566   1                1668   0.454545  HOFFMANN_IMMATURE_TO_MATURE_B_LYMPHOCYTE_DN
Top Late Enhancement                 Hoffmann Immature To Mature B Lymphocyte Dn                 -0.501326  0.0112721    -1.69833  0.226338   1                9195   0.515152  HOFFMANN_IMMATURE_TO_MATURE_B_LYMPHOCYTE_DN
Top Late Enhancement                 Lee Early T Lymphocyte Up                                   -0.698908  0.0581138    -1.69505  0.226338   1               10103   0.642857  LEE_EARLY_T_LYMPHOCYTE_UP
Variation Smoothness All Timeframes  Lee Intrathymic T Progenitor                                -0.599872  0.027688     -1.61294  0.244212   1                8128   0.625     LEE_INTRATHYMIC_T_PROGENITOR
Mean Smoothness All Timeframes       Lee Intrathymic T Progenitor                                -0.584748  0.0407422    -1.57614  0.216178   1                8176   0.625     LEE_INTRATHYMIC_T_PROGENITOR
Variation Smoothness All Timeframes  Shedden Lung Cancer Good Survival A4                        -0.450663  0.00574941   -1.72629  0.211541   1                7958   0.492958  SHEDDEN_LUNG_CANCER_GOOD_SURVIVAL_A4
Mean Smoothness All Timeframes       Shedden Lung Cancer Good Survival A4                        -0.420924  0.016097     -1.60539  0.204175   1                7445   0.542254  SHEDDEN_LUNG_CANCER_GOOD_SURVIVAL_A4
Mean Smoothness All Timeframes       Meissner Npc Hcp With H3K4Me3 And H3K27Me3                  -0.475314  0.0176587    -1.64993  0.190145   1                8362   0.467742  MEISSNER_NPC_HCP_WITH_H3K4ME3_AND_H3K27ME3
Mean Smoothness All Timeframes       Meissner Npc Hcp With H3 Unmethylated                       -0.380652  0.0227813    -1.56946  0.21886    1                8723   0.347826  MEISSNER_NPC_HCP_WITH_H3_UNMETHYLATED
Variation Smoothness All Timeframes  Tavazoie Metastasis                                         -0.5638    0.000802568  -1.89888  0.168147   1                9127   0.442308  TAVAZOIE_METASTASIS
Mean Smoothness All Timeframes       Tavazoie Metastasis                                         -0.541672  0.00322386   -1.82128  0.133542   1                8996   0.461538  TAVAZOIE_METASTASIS
Mean Smoothness All Timeframes       Mili Pseudopodia                                            -0.436891  0.0397032    -1.53485  0.232143   1                8288   0.5       MILI_PSEUDOPODIA
Top Late Enhancement                 Croonquist Il6 Deprivation Dn                               -0.76839   0.0198582    -1.74754  0.213855   1                9993   0.658824  CROONQUIST_IL6_DEPRIVATION_DN
Top Late Enhancement                 Croonquist Nras Signaling Dn                                -0.750563  0.0344199    -1.69461  0.226338   1                9864   0.677419  CROONQUIST_NRAS_SIGNALING_DN
Variation Smoothness All Timeframes  Iritani Mad1 Targets Dn                                      0.625435  0.0124481     1.79871  0.247899   1                2135   0.681818  IRITANI_MAD1_TARGETS_DN
Variation Smoothness Uptake          Iritani Mad1 Targets Dn                                      0.621501  0.0132463     1.78768  0.186182   1                2314   0.704545  IRITANI_MAD1_TARGETS_DN
Mean Smoothness All Timeframes       Zhan Late Differentiation Genes Up                          -0.52583   0.0292654    -1.63356  0.194532   1                9368   0.392857  ZHAN_LATE_DIFFERENTIATION_GENES_UP
Top Late Enhancement                 Ishida E2F Targets                                          -0.797135  0.0266531    -1.65932  0.248326   1               10417   0.690476  ISHIDA_E2F_TARGETS
Variation Smoothness All Timeframes  Valk Aml Cluster 1                                          -0.630519  0.0199362    -1.70003  0.219274   1                9706   0.45      VALK_AML_CLUSTER_1
Mean Smoothness All Timeframes       Valk Aml Cluster 1                                          -0.619223  0.0234516    -1.66795  0.180662   1                9434   0.5       VALK_AML_CLUSTER_1
Irregularity                         Valk Aml Cluster 2                                           0.663559  0.00336035    1.89035  0.185049   1                2233   0.5625    VALK_AML_CLUSTER_2
Variation Smoothness All Timeframes  Valk Aml Cluster 2                                          -0.578164  0.0217391    -1.6496   0.236572   1                8230   0.5625    VALK_AML_CLUSTER_2
Mean Smoothness All Timeframes       Valk Aml Cluster 2                                          -0.615027  0.00819345   -1.76137  0.153387   1                8673   0.5625    VALK_AML_CLUSTER_2
Variation Smoothness All Timeframes  Valk Aml Cluster 4                                          -0.570024  0.0239856    -1.63662  0.236712   1                9362   0.4       VALK_AML_CLUSTER_4
Mean Smoothness All Timeframes       Valk Aml Cluster 4                                          -0.566997  0.0276221    -1.62808  0.194532   1                9058   0.45      VALK_AML_CLUSTER_4
Variation Smoothness All Timeframes  Valk Aml Cluster 8                                          -0.556253  0.0248608    -1.61612  0.244212   1                8345   0.6875    VALK_AML_CLUSTER_8
Mean Smoothness All Timeframes       Valk Aml Cluster 8                                          -0.560173  0.0211588    -1.62271  0.198215   1                8066   0.6875    VALK_AML_CLUSTER_8
Variation Smoothness All Timeframes  Valk Aml Cluster 10                                         -0.651043  0.000994827  -1.91021  0.168147   1                8924   0.535714  VALK_AML_CLUSTER_10
Mean Smoothness All Timeframes       Valk Aml Cluster 10                                         -0.617478  0.00277888   -1.81332  0.139792   1                9278   0.5       VALK_AML_CLUSTER_10
Variation Smoothness All Timeframes  Valk Aml Cluster 12                                         -0.577537  0.0226595    -1.6667   0.236189   1                8509   0.6       VALK_AML_CLUSTER_12
Mean Smoothness All Timeframes       Valk Aml Cluster 12                                         -0.625578  0.0085777    -1.79778  0.141179   1                8839   0.6       VALK_AML_CLUSTER_12
Irregularity                         Valk Aml Cluster 16                                          0.620011  0.00376089    1.82053  0.224248   1                1677   0.352941  VALK_AML_CLUSTER_16
Mean Smoothness All Timeframes       Valk Aml With Evi1                                          -0.635929  0.0051731    -1.7967   0.141179   1                9278   0.55      VALK_AML_WITH_EVI1
Variation Smoothness All Timeframes  Valk Aml With Flt3 Itd                                      -0.547937  0.0123133    -1.69573  0.219274   1                8214   0.56      VALK_AML_WITH_FLT3_ITD
Mean Smoothness All Timeframes       Valk Aml With Flt3 Itd                                      -0.549122  0.0119904    -1.70545  0.164501   1                8673   0.52      VALK_AML_WITH_FLT3_ITD
Variation Smoothness All Timeframes  Poola Invasive Breast Cancer Dn                             -0.522191  0.0198079    -1.68026  0.232491   1                6598   0.788462  POOLA_INVASIVE_BREAST_CANCER_DN
Mean Smoothness All Timeframes       Poola Invasive Breast Cancer Dn                             -0.527199  0.0152183    -1.69229  0.169834   1                7542   0.625     POOLA_INVASIVE_BREAST_CANCER_DN
Variation Smoothness All Timeframes  Ding Lung Cancer By Mutation Rate                           -0.595717  0.0189356    -1.64839  0.236712   1                7813   0.555556  DING_LUNG_CANCER_BY_MUTATION_RATE
Ld Late Lt0                          Boyault Liver Cancer Subclass G2                             0.601429  0.000811194   1.93028  0.24158    1                2691   0.541667  BOYAULT_LIVER_CANCER_SUBCLASS_G2
Variation Smoothness All Timeframes  Boyault Liver Cancer Subclass G2                            -0.529233  0.00885312   -1.69494  0.219274   1               10376   0.291667  BOYAULT_LIVER_CANCER_SUBCLASS_G2
Mean Smoothness All Timeframes       Boyault Liver Cancer Subclass G2                            -0.500696  0.0201613    -1.61033  0.202604   1                9429   0.375     BOYAULT_LIVER_CANCER_SUBCLASS_G2
Top Late Enhancement                 Boyault Liver Cancer Subclass G3 Up                         -0.563211  0.0301558    -1.71008  0.226338   1                8713   0.527778  BOYAULT_LIVER_CANCER_SUBCLASS_G3_UP
Mean Smoothness All Timeframes       Boyault Liver Cancer Subclass G56 Up                        -0.613694  0.0373604    -1.58016  0.216167   1               10172   0.4       BOYAULT_LIVER_CANCER_SUBCLASS_G56_UP
Variation Smoothness All Timeframes  Chiang Liver Cancer Subclass Ctnnb1 Up                      -0.499155  0.000199402  -1.90276  0.168147   0.531007         8328   0.5       CHIANG_LIVER_CANCER_SUBCLASS_CTNNB1_UP
Mean Smoothness All Timeframes       Chiang Liver Cancer Subclass Ctnnb1 Up                      -0.487129  0.000598205  -1.85756  0.129188   1                8054   0.548077  CHIANG_LIVER_CANCER_SUBCLASS_CTNNB1_UP
Mean Smoothness All Timeframes       Chiang Liver Cancer Subclass Ctnnb1 Dn                      -0.514531  0.0450596    -1.63944  0.194025   1                9533   0.381679  CHIANG_LIVER_CANCER_SUBCLASS_CTNNB1_DN
Washout                              Chiang Liver Cancer Subclass Proliferation Up                0.612532  0.0281231     1.77439  0.235363   1                2124   0.630872  CHIANG_LIVER_CANCER_SUBCLASS_PROLIFERATION_UP
Top Late Enhancement                 Chiang Liver Cancer Subclass Proliferation Up               -0.639958  0.0126837    -1.8689   0.213855   1                9433   0.604027  CHIANG_LIVER_CANCER_SUBCLASS_PROLIFERATION_UP
Variation Smoothness All Timeframes  Chiang Liver Cancer Subclass Proliferation Dn               -0.503173  0.0143814    -1.70959  0.213798   1                7932   0.514286  CHIANG_LIVER_CANCER_SUBCLASS_PROLIFERATION_DN
Mean Smoothness All Timeframes       Chiang Liver Cancer Subclass Proliferation Dn               -0.469785  0.0369574    -1.59268  0.210308   1                9693   0.3       CHIANG_LIVER_CANCER_SUBCLASS_PROLIFERATION_DN
Variation Smoothness All Timeframes  Chiang Liver Cancer Subclass Interferon Dn                  -0.475763  0.021135     -1.60292  0.247184   1               10064   0.363636  CHIANG_LIVER_CANCER_SUBCLASS_INTERFERON_DN
Mean Smoothness All Timeframes       Chiang Liver Cancer Subclass Interferon Dn                  -0.513268  0.00770751   -1.71939  0.16069    1               10045   0.363636  CHIANG_LIVER_CANCER_SUBCLASS_INTERFERON_DN
Washout                              Chiang Liver Cancer Subclass Unannotated Dn                  0.55053   0.0178143     1.75706  0.242201   1                2331   0.538462  CHIANG_LIVER_CANCER_SUBCLASS_UNANNOTATED_DN
Top Late Enhancement                 Chiang Liver Cancer Subclass Unannotated Dn                 -0.54131   0.0255379    -1.71625  0.223749   1                8948   0.538462  CHIANG_LIVER_CANCER_SUBCLASS_UNANNOTATED_DN
Variation Smoothness Uptake          Chiang Liver Cancer Subclass Unannotated Dn                  0.5518    0.0178923     1.75461  0.193482   1                2234   0.565934  CHIANG_LIVER_CANCER_SUBCLASS_UNANNOTATED_DN
Washout                              Kaposi Liver Cancer Met Up                                   0.667923  0.00322191    1.85566  0.193047   1                3167   0.833333  KAPOSI_LIVER_CANCER_MET_UP
Top Late Enhancement                 Kaposi Liver Cancer Met Up                                  -0.656011  0.00579768   -1.8241   0.213855   1                8411   0.777778  KAPOSI_LIVER_CANCER_MET_UP
Vol Late Lt0                         Kaposi Liver Cancer Met Up                                   0.691636  0.00221283    1.92145  0.214618   1                2663   0.722222  KAPOSI_LIVER_CANCER_MET_UP
Variation Smoothness All Timeframes  Meissner Npc Hcp With H3K4Me2 And H3K27Me3                  -0.471988  0.00679864   -1.71661  0.213798   1                8995   0.350515  MEISSNER_NPC_HCP_WITH_H3K4ME2_AND_H3K27ME3
Mean Smoothness All Timeframes       Meissner Npc Hcp With H3K4Me2 And H3K27Me3                  -0.441961  0.0199919    -1.60111  0.204175   1                7443   0.494845  MEISSNER_NPC_HCP_WITH_H3K4ME2_AND_H3K27ME3
Mean Smoothness All Timeframes       Meissner Npc Hcp With H3K4Me2                               -0.382354  0.0167866    -1.5815   0.216167   1                8252   0.398374  MEISSNER_NPC_HCP_WITH_H3K4ME2
Variation Smoothness All Timeframes  Mikkelsen Mcv6 Hcp With H3K27Me3                            -0.476163  0.015625     -1.66863  0.236189   1                8241   0.495726  MIKKELSEN_MCV6_HCP_WITH_H3K27ME3
Mean Smoothness All Timeframes       Mikkelsen Mcv6 Hcp With H3K27Me3                            -0.491328  0.00896414   -1.71741  0.16069    1                8755   0.452991  MIKKELSEN_MCV6_HCP_WITH_H3K27ME3
Mean Smoothness All Timeframes       Mikkelsen Mcv6 Icp With H3K4Me3 And H3K27Me3                -0.642604  0.0499401    -1.54917  0.226989   1                7517   0.833333  MIKKELSEN_MCV6_ICP_WITH_H3K4ME3_AND_H3K27ME3
Variation Smoothness All Timeframes  Mikkelsen Ips Icp With H3K4Me3 And H327Me3                  -0.55883   0.0153264    -1.72672  0.211541   1                9124   0.408163  MIKKELSEN_IPS_ICP_WITH_H3K4ME3_AND_H327ME3
Mean Smoothness All Timeframes       Mikkelsen Ips Icp With H3K4Me3 And H327Me3                  -0.610407  0.00256815   -1.88517  0.117961   1                9106   0.44898   MIKKELSEN_IPS_ICP_WITH_H3K4ME3_AND_H327ME3
Mean Smoothness All Timeframes       Nakayama Fgf2 Targets                                       -0.560665  0.0460633    -1.58923  0.212416   1                9819   0.454545  NAKAYAMA_FGF2_TARGETS
Washout                              Ngo Malignant Glioma 1P Loh                                  0.761614  0.00262203    1.89474  0.193047   1                2038   0.8125    NGO_MALIGNANT_GLIOMA_1P_LOH
Top Late Enhancement                 Ngo Malignant Glioma 1P Loh                                 -0.687269  0.0210294    -1.70076  0.226338   1                8884   0.8125    NGO_MALIGNANT_GLIOMA_1P_LOH
Vol Late Lt0                         Ngo Malignant Glioma 1P Loh                                  0.767738  0.00101502    1.90335  0.214618   1                1828   0.75      NGO_MALIGNANT_GLIOMA_1P_LOH
Mean Smoothness All Timeframes       Li Prostate Cancer Epigenetic                               -0.556001  0.0231913    -1.63411  0.194532   1                8680   0.56      LI_PROSTATE_CANCER_EPIGENETIC
Mean Smoothness All Timeframes       Gyorffy Doxorubicin Resistance                              -0.517883  0.0122416    -1.68794  0.172248   1                9758   0.354839  GYORFFY_DOXORUBICIN_RESISTANCE
Mean Smoothness All Timeframes       Gyorffy Mitoxantrone Resistance                             -0.482475  0.0265754    -1.59158  0.211187   1                9690   0.351351  GYORFFY_MITOXANTRONE_RESISTANCE
Mean Smoothness All Timeframes       Korkola Embryonic Carcinoma Vs Seminoma Dn                  -0.644398  0.0323612    -1.58318  0.216167   1                7954   0.75      KORKOLA_EMBRYONIC_CARCINOMA_VS_SEMINOMA_DN
Top Late Enhancement                 Kobayashi Egfr Signaling 24Hr Dn                            -0.678308  0.0264389    -1.77002  0.213855   1                9993   0.566327  KOBAYASHI_EGFR_SIGNALING_24HR_DN
Top Late Enhancement                 Fournier Acinar Development Late 2                          -0.535986  0.0270871    -1.72249  0.219645   1                8754   0.552     FOURNIER_ACINAR_DEVELOPMENT_LATE_2
Irregularity                         Cairo Hepatoblastoma Dn                                      0.504685  0.00040858    1.96704  0.12903    1                3006   0.527027  CAIRO_HEPATOBLASTOMA_DN
Variation Smoothness All Timeframes  Cairo Hepatoblastoma Classes Dn                             -0.541362  0.00677831   -1.84584  0.171536   1                9104   0.429487  CAIRO_HEPATOBLASTOMA_CLASSES_DN
Mean Smoothness All Timeframes       Cairo Hepatoblastoma Classes Dn                             -0.559347  0.00357498   -1.90944  0.117961   1                9098   0.461538  CAIRO_HEPATOBLASTOMA_CLASSES_DN
Mean Smoothness All Timeframes       Cairo Liver Development Up                                  -0.417355  0.0486476    -1.55656  0.223946   1                7533   0.55303   CAIRO_LIVER_DEVELOPMENT_UP
Top Late Enhancement                 Winnepenninckx Melanoma Metastasis Up                       -0.663627  0.019806     -1.79389  0.213855   1                9089   0.664234  WINNEPENNINCKX_MELANOMA_METASTASIS_UP
Variation Smoothness All Timeframes  Uzonyi Response To Leukotriene And Thrombin                 -0.68885   0.0620093    -1.61502  0.244212   1               10039   0.5       UZONYI_RESPONSE_TO_LEUKOTRIENE_AND_THROMBIN
Mean Smoothness All Timeframes       Uzonyi Response To Leukotriene And Thrombin                 -0.674636  0.0805342    -1.5751   0.216178   1                9434   0.5625    UZONYI_RESPONSE_TO_LEUKOTRIENE_AND_THROMBIN
Variation Smoothness All Timeframes  Zhan Multiple Myeloma Mf Dn                                 -0.574744  0.00579073   -1.78015  0.193127   1                8939   0.5       ZHAN_MULTIPLE_MYELOMA_MF_DN
Mean Smoothness All Timeframes       Zhan Multiple Myeloma Mf Dn                                 -0.590979  0.00200642   -1.83189  0.133542   1                8854   0.533333  ZHAN_MULTIPLE_MYELOMA_MF_DN
Variation Smoothness Uptake          Dang Regulated By Myc Up                                     0.570503  0.00756369    1.78387  0.187272   1                2008   0.515152  DANG_REGULATED_BY_MYC_UP
Variation Smoothness All Timeframes  Dang Myc Targets Up                                          0.618832  0.00534442    1.86558  0.212371   1                2570   0.688     DANG_MYC_TARGETS_UP
Variation Smoothness Uptake          Dang Myc Targets Up                                          0.629018  0.00257222    1.89948  0.173515   1                2491   0.688     DANG_MYC_TARGETS_UP
Washout                              Wong Embryonic Stem Cell Core                                0.618509  0.0146015     1.79648  0.221302   1                2388   0.598662  WONG_EMBRYONIC_STEM_CELL_CORE
Top Late Enhancement                 Wong Embryonic Stem Cell Core                               -0.629037  0.011953     -1.81547  0.213855   1                9083   0.568562  WONG_EMBRYONIC_STEM_CELL_CORE
Variation Smoothness All Timeframes  Wong Embryonic Stem Cell Core                                0.678034  0.00119474    1.95717  0.206872   1                2407   0.715719  WONG_EMBRYONIC_STEM_CELL_CORE
Variation Smoothness Uptake          Wong Embryonic Stem Cell Core                                0.67924   0.000791609   1.97038  0.173515   1                2257   0.698997  WONG_EMBRYONIC_STEM_CELL_CORE
Irregularity                         Nakayama Soft Tissue Tumors Pca1 Dn                          0.607529  0.00119       1.9543   0.12903    1                2471   0.574468  NAKAYAMA_SOFT_TISSUE_TUMORS_PCA1_DN
Variation Smoothness All Timeframes  Nakayama Soft Tissue Tumors Pca1 Dn                         -0.526288  0.0160178    -1.69129  0.221342   1                8189   0.489362  NAKAYAMA_SOFT_TISSUE_TUMORS_PCA1_DN
Mean Smoothness All Timeframes       Nakayama Soft Tissue Tumors Pca1 Dn                         -0.543042  0.00915937   -1.74821  0.153387   1                9761   0.319149  NAKAYAMA_SOFT_TISSUE_TUMORS_PCA1_DN
Circularity                          Nakayama Soft Tissue Tumors Pca2 Dn                         -0.74361   0.00298626   -1.9641   0.223353   1               10349   0.6       NAKAYAMA_SOFT_TISSUE_TUMORS_PCA2_DN
Irregularity                         Nakayama Soft Tissue Tumors Pca2 Dn                          0.738409  0.00360288    1.94933  0.12903    1                1840   0.709091  NAKAYAMA_SOFT_TISSUE_TUMORS_PCA2_DN
Ld Late Lt0                          Nakayama Soft Tissue Tumors Pca2 Dn                          0.77709   0.000995619   2.05126  0.0973724  1                 605   0.563636  NAKAYAMA_SOFT_TISSUE_TUMORS_PCA2_DN
Variation Smoothness All Timeframes  Nakayama Soft Tissue Tumors Pca2 Dn                         -0.764187  0.00119261   -2.01298  0.147564   1               10183   0.636364  NAKAYAMA_SOFT_TISSUE_TUMORS_PCA2_DN
Mean Smoothness All Timeframes       Nakayama Soft Tissue Tumors Pca2 Dn                         -0.730452  0.00533913   -1.94401  0.109401   1                9730   0.654545  NAKAYAMA_SOFT_TISSUE_TUMORS_PCA2_DN
Variation Smoothness All Timeframes  Chandran Metastasis Dn                                      -0.491466  0.00218993   -1.87961  0.168147   1                9314   0.385246  CHANDRAN_METASTASIS_DN
Mean Smoothness All Timeframes       Chandran Metastasis Dn                                      -0.513439  0.000601082  -1.95725  0.109401   1                9083   0.430328  CHANDRAN_METASTASIS_DN
Mean Smoothness All Timeframes       Mikkelsen Es Hcp With H3 Unmethylated                       -0.618225  0.0266508    -1.61934  0.200738   1               10658   0.307692  MIKKELSEN_ES_HCP_WITH_H3_UNMETHYLATED
Variation Smoothness All Timeframes  Mikkelsen Es Icp With H3K4Me3 And H3K27Me3                  -0.533115  0.00869917   -1.75344  0.193127   1                8233   0.517241  MIKKELSEN_ES_ICP_WITH_H3K4ME3_AND_H3K27ME3
Mean Smoothness All Timeframes       Mikkelsen Es Icp With H3K4Me3 And H3K27Me3                  -0.502362  0.0174209    -1.65558  0.184559   1                9106   0.413793  MIKKELSEN_ES_ICP_WITH_H3K4ME3_AND_H3K27ME3
Variation Smoothness All Timeframes  Mikkelsen Npc Hcp With H3K4Me3 And H3K27Me3                 -0.455384  0.018898     -1.6428   0.236712   1                9038   0.358025  MIKKELSEN_NPC_HCP_WITH_H3K4ME3_AND_H3K27ME3
Mean Smoothness All Timeframes       Mikkelsen Npc Hcp With H3K4Me3 And H3K27Me3                 -0.476257  0.0090569    -1.72096  0.16069    1                8362   0.45679   MIKKELSEN_NPC_HCP_WITH_H3K4ME3_AND_H3K27ME3
Variation Smoothness All Timeframes  Yao Temporal Response To Progesterone Cluster 0             -0.507038  0.023781     -1.67048  0.235528   1                9358   0.380952  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_0
Mean Smoothness All Timeframes       Yao Temporal Response To Progesterone Cluster 0             -0.514595  0.0166366    -1.69329  0.169834   1                8077   0.52381   YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_0
Mean Smoothness All Timeframes       Yao Temporal Response To Progesterone Cluster 1             -0.454527  0.0369261    -1.54066  0.231366   1                8059   0.516667  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_1
Washout                              Yao Temporal Response To Progesterone Cluster 13             0.597194  0.0160869     1.79046  0.221302   1                2083   0.550633  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_13
Top Late Enhancement                 Yao Temporal Response To Progesterone Cluster 13            -0.566594  0.0372       -1.69828  0.226338   1                8795   0.582278  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_13
Variation Smoothness All Timeframes  Yao Temporal Response To Progesterone Cluster 13             0.63983   0.00336901    1.9266   0.206872   1                2159   0.613924  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_13
Variation Smoothness Uptake          Yao Temporal Response To Progesterone Cluster 13             0.643336  0.00219036    1.93579  0.173515   1                1887   0.56962   YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_13
Washout                              Yao Temporal Response To Progesterone Cluster 14             0.622434  0.000605938   2.02875  0.193047   1                2614   0.664122  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_14
Ser                                  Yao Temporal Response To Progesterone Cluster 14             0.624478  0.000405927   2.02902  0.21078    1                2005   0.541985  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_14
Top Late Enhancement                 Yao Temporal Response To Progesterone Cluster 14            -0.644173  0.000202388  -2.09202  0.131489   0.53896          9131   0.59542   YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_14
Vol Late Lt0                         Yao Temporal Response To Progesterone Cluster 14             0.623858  0.000201329   2.04096  0.206783   0.536139         1885   0.541985  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_14
Variation Smoothness Uptake          Yao Temporal Response To Progesterone Cluster 14             0.542264  0.0123236     1.76289  0.192501   1                2632   0.580153  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_14
Washout                              Yao Temporal Response To Progesterone Cluster 17             0.589444  0.00161878    1.94694  0.193047   1                2077   0.526627  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_17
Ser                                  Yao Temporal Response To Progesterone Cluster 17             0.599108  0.00121803    1.98627  0.21078    1                2010   0.502959  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_17
Top Late Enhancement                 Yao Temporal Response To Progesterone Cluster 17            -0.563061  0.00544465   -1.85206  0.213855   1                8764   0.56213   YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_17
Vol Late Lt0                         Yao Temporal Response To Progesterone Cluster 17             0.583858  0.00182149    1.92731  0.214618   1                1605   0.455621  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_17
Variation Smoothness Uptake          Yao Temporal Response To Progesterone Cluster 17             0.556284  0.00438247    1.83764  0.173515   1                2278   0.52071   YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_17
Variation Smoothness All Timeframes  Nakamura Adipogenesis Early Up                              -0.496836  0.0399521    -1.60435  0.247184   1                8104   0.583333  NAKAMURA_ADIPOGENESIS_EARLY_UP
Mean Smoothness All Timeframes       Nakamura Adipogenesis Early Up                              -0.545189  0.0104292    -1.75587  0.153387   1                8067   0.666667  NAKAMURA_ADIPOGENESIS_EARLY_UP
Mean Smoothness All Timeframes       Nakamura Adipogenesis Early Dn                              -0.622083  0.0520104    -1.61182  0.202062   1                9174   0.588235  NAKAMURA_ADIPOGENESIS_EARLY_DN
Variation Smoothness All Timeframes  Nakamura Adipogenesis Late Up                               -0.454684  0.00120919   -1.80031  0.193127   1                8955   0.380435  NAKAMURA_ADIPOGENESIS_LATE_UP
Mean Smoothness All Timeframes       Nakamura Adipogenesis Late Up                               -0.460697  0.00138999   -1.83302  0.133542   1                8443   0.467391  NAKAMURA_ADIPOGENESIS_LATE_UP
Mean Smoothness All Timeframes       Nakamura Adipogenesis Late Dn                               -0.586739  0.0706608    -1.57022  0.218767   1                9174   0.617647  NAKAMURA_ADIPOGENESIS_LATE_DN
Mean Smoothness All Timeframes       Zembutsu Sensitivity To Cyclophosphamide                    -0.588453  0.0283133    -1.60409  0.204175   1                9759   0.428571  ZEMBUTSU_SENSITIVITY_TO_CYCLOPHOSPHAMIDE
Top Late Enhancement                 Dazard Uv Response Cluster G5                               -0.779713  0.0126506    -1.68754  0.230757   1                9330   0.6       DAZARD_UV_RESPONSE_CLUSTER_G5
Variation Smoothness All Timeframes  Ruan Response To Troglitazone Dn                            -0.631263  0.00944154   -1.77396  0.193127   1               11068   0.1875    RUAN_RESPONSE_TO_TROGLITAZONE_DN
Mean Smoothness All Timeframes       Ruan Response To Troglitazone Dn                            -0.612764  0.0123703    -1.72201  0.16069    1               11003   0.1875    RUAN_RESPONSE_TO_TROGLITAZONE_DN
Variation Smoothness All Timeframes  Nielsen Gist                                                -0.470688  0.00462591   -1.76195  0.193127   1                8947   0.39759   NIELSEN_GIST
Mean Smoothness All Timeframes       Nielsen Gist                                                -0.4504    0.0127608    -1.68563  0.172248   1                8441   0.421687  NIELSEN_GIST
Mean Smoothness All Timeframes       Nielsen Leiomyosarcoma Cnn1 Dn                              -0.691655  0.042293     -1.60422  0.204175   1                9626   0.588235  NIELSEN_LEIOMYOSARCOMA_CNN1_DN
Irregularity                         Nielsen Leiomyosarcoma Dn                                    0.659617  0.00520104    1.80487  0.245668   1                2281   0.625     NIELSEN_LEIOMYOSARCOMA_DN
Variation Smoothness All Timeframes  Nielsen Leiomyosarcoma Dn                                   -0.623925  0.0136218    -1.71787  0.213798   1                8830   0.6875    NIELSEN_LEIOMYOSARCOMA_DN
Mean Smoothness All Timeframes       Nielsen Leiomyosarcoma Dn                                   -0.644455  0.00861033   -1.77435  0.149917   1                8066   0.8125    NIELSEN_LEIOMYOSARCOMA_DN
Mean Smoothness All Timeframes       Nielsen Malignat Fibrous Histiocytoma Dn                    -0.53281   0.0521967    -1.52971  0.23533    1                8744   0.4375    NIELSEN_MALIGNAT_FIBROUS_HISTIOCYTOMA_DN
Mean Smoothness All Timeframes       Nielsen Schwannoma Dn                                       -0.725112  0.011131     -1.76337  0.153387   1                9167   0.705882  NIELSEN_SCHWANNOMA_DN
Mean Smoothness All Timeframes       Browne Hcmv Infection 30Min Up                              -0.489345  0.0599637    -1.52123  0.242693   1                9476   0.40625   BROWNE_HCMV_INFECTION_30MIN_UP
Top Late Enhancement                 Karlsson Tgfb1 Targets Up                                   -0.499641  0.0260469    -1.67587  0.232024   1                8948   0.5       KARLSSON_TGFB1_TARGETS_UP
Mean Smoothness All Timeframes       Noushmehr Gbm Silenced By Methylation                       -0.536663  0.0129712    -1.66913  0.180021   1                8258   0.571429  NOUSHMEHR_GBM_SILENCED_BY_METHYLATION
Variation Smoothness All Timeframes  Stambolsky Response To Vitamin D3 Up                        -0.457104  0.0120797    -1.64     0.236712   1                7332   0.571429  STAMBOLSKY_RESPONSE_TO_VITAMIN_D3_UP
Mean Smoothness All Timeframes       Stambolsky Response To Vitamin D3 Up                        -0.472761  0.0076       -1.69129  0.170151   1                7621   0.555556  STAMBOLSKY_RESPONSE_TO_VITAMIN_D3_UP
Mean Smoothness All Timeframes       Kim All Disorders Oligodendrocyte Number Corr Dn            -0.651603  0.0128257    -1.70351  0.165416   1                8762   0.625     KIM_ALL_DISORDERS_OLIGODENDROCYTE_NUMBER_CORR_DN
Top Late Enhancement                 Chicas Rb1 Targets Low Serum                                -0.550224  0.0237902    -1.68862  0.230757   1                9809   0.486486  CHICAS_RB1_TARGETS_LOW_SERUM
Variation Smoothness All Timeframes  Hoelzel Nf1 Targets Up                                      -0.491113  0.022262     -1.65248  0.236572   1                9386   0.356322  HOELZEL_NF1_TARGETS_UP
Mean Smoothness All Timeframes       Hoelzel Nf1 Targets Up                                      -0.560329  0.00420336   -1.87863  0.117961   1                8773   0.471264  HOELZEL_NF1_TARGETS_UP
Variation Smoothness All Timeframes  Hoelzel Nf1 Targets Dn                                      -0.577877  0.0142114    -1.77894  0.193127   1                9107   0.506667  HOELZEL_NF1_TARGETS_DN
Mean Smoothness All Timeframes       Hoelzel Nf1 Targets Dn                                      -0.643423  0.00160514   -1.98289  0.109401   1                9058   0.573333  HOELZEL_NF1_TARGETS_DN
Variation Smoothness All Timeframes  Dutertre Estradiol Response 6Hr Dn                          -0.541298  0.00261833   -1.86325  0.168147   1                9436   0.376471  DUTERTRE_ESTRADIOL_RESPONSE_6HR_DN
Mean Smoothness All Timeframes       Dutertre Estradiol Response 6Hr Dn                          -0.558993  0.00181232   -1.9251   0.117961   1                9173   0.411765  DUTERTRE_ESTRADIOL_RESPONSE_6HR_DN
Irregularity                         Figueroa Aml Methylation Cluster 1 Up                        0.491063  0.000586969   1.96283  0.12903    1                1762   0.328571  FIGUEROA_AML_METHYLATION_CLUSTER_1_UP
Variation Smoothness All Timeframes  Figueroa Aml Methylation Cluster 1 Up                       -0.436939  0.00140168   -1.72767  0.211541   1                9369   0.342857  FIGUEROA_AML_METHYLATION_CLUSTER_1_UP
Mean Smoothness All Timeframes       Figueroa Aml Methylation Cluster 1 Up                       -0.399472  0.0119856    -1.5812   0.216167   1                7987   0.457143  FIGUEROA_AML_METHYLATION_CLUSTER_1_UP
Variation Smoothness All Timeframes  Figueroa Aml Methylation Cluster 3 Dn                       -0.561532  0.00569633   -1.75456  0.193127   1                8821   0.5       FIGUEROA_AML_METHYLATION_CLUSTER_3_DN
Mean Smoothness All Timeframes       Figueroa Aml Methylation Cluster 3 Dn                       -0.50661   0.0299306    -1.58293  0.216167   1                7769   0.615385  FIGUEROA_AML_METHYLATION_CLUSTER_3_DN
Circularity                          Wang Classic Adipogenic Targets Of Pparg                    -0.777567  0.000400962  -2.03782  0.205876   1                9868   0.5625    WANG_CLASSIC_ADIPOGENIC_TARGETS_OF_PPARG
Irregularity                         Wang Classic Adipogenic Targets Of Pparg                     0.696219  0.00784235    1.82533  0.218099   1                 840   0.375     WANG_CLASSIC_ADIPOGENIC_TARGETS_OF_PPARG
Variation Smoothness All Timeframes  Pangas Tumor Suppression By Smad1 And Smad5 Up              -0.495507  0.0456182    -1.62164  0.244212   1                8857   0.471698  PANGAS_TUMOR_SUPPRESSION_BY_SMAD1_AND_SMAD5_UP
Mean Smoothness All Timeframes       Pangas Tumor Suppression By Smad1 And Smad5 Up              -0.493619  0.0511182    -1.61525  0.200977   1                8690   0.509434  PANGAS_TUMOR_SUPPRESSION_BY_SMAD1_AND_SMAD5_UP
Top Late Enhancement                 Li Dcp2 Bound Mrna                                          -0.574524  0.0263528    -1.6966   0.226338   1                8931   0.641026  LI_DCP2_BOUND_MRNA
Variation Smoothness All Timeframes  Li Dcp2 Bound Mrna                                           0.661316  0.000598444   1.96225  0.206872   1                2199   0.641026  LI_DCP2_BOUND_MRNA
Variation Smoothness Uptake          Li Dcp2 Bound Mrna                                           0.62774   0.00520938    1.86054  0.173515   1                2407   0.641026  LI_DCP2_BOUND_MRNA
Irregularity                         Ohguchi Liver Hnf4A Targets Up                               0.652377  0.000800641   1.96921  0.12903    1                1230   0.458333  OHGUCHI_LIVER_HNF4A_TARGETS_UP
Variation Smoothness All Timeframes  Ohguchi Liver Hnf4A Targets Up                              -0.569946  0.0140562    -1.71487  0.213798   1               10034   0.333333  OHGUCHI_LIVER_HNF4A_TARGETS_UP
Mean Smoothness All Timeframes       Ohguchi Liver Hnf4A Targets Up                              -0.577782  0.0120943    -1.73467  0.16069    1                7532   0.666667  OHGUCHI_LIVER_HNF4A_TARGETS_UP
Mean Smoothness All Timeframes       Gabriely Mir21 Targets                                      -0.494011  0.0452906    -1.63033  0.194532   1                8433   0.456929  GABRIELY_MIR21_TARGETS
Mean Smoothness All Timeframes       Chyla Cbfa2T3 Targets Dn                                    -0.349144  0.0114353    -1.53595  0.2318     1                7880   0.440789  CHYLA_CBFA2T3_TARGETS_DN
Variation Smoothness All Timeframes  Bhat Esr1 Targets Not Via Akt1 Dn                           -0.532525  0.0160643    -1.71555  0.213798   1                8379   0.512821  BHAT_ESR1_TARGETS_NOT_VIA_AKT1_DN
Mean Smoothness All Timeframes       Bhat Esr1 Targets Not Via Akt1 Dn                           -0.553247  0.00684242   -1.78021  0.149917   1                8268   0.551282  BHAT_ESR1_TARGETS_NOT_VIA_AKT1_DN
Variation Smoothness All Timeframes  Bhat Esr1 Targets Via Akt1 Dn                               -0.540191  0.00538492   -1.78662  0.193127   1                7387   0.675676  BHAT_ESR1_TARGETS_VIA_AKT1_DN
Mean Smoothness All Timeframes       Bhat Esr1 Targets Via Akt1 Dn                               -0.537447  0.00421856   -1.77455  0.149917   1                7506   0.662162  BHAT_ESR1_TARGETS_VIA_AKT1_DN
Mean Smoothness All Timeframes       Johnstone Parvb Targets 1 Dn                                -0.493337  0.0350524    -1.60214  0.204175   1                8567   0.508772  JOHNSTONE_PARVB_TARGETS_1_DN
Mean Smoothness All Timeframes       Miyagawa Targets Of Ewsr1 Ets Fusions Up                    -0.390965  0.0401598    -1.52813  0.236945   1                7616   0.510101  MIYAGAWA_TARGETS_OF_EWSR1_ETS_FUSIONS_UP
Variation Smoothness All Timeframes  Miyagawa Targets Of Ewsr1 Ets Fusions Dn                    -0.501564  0.0605938    -1.61095  0.244212   1                8262   0.520468  MIYAGAWA_TARGETS_OF_EWSR1_ETS_FUSIONS_DN
Mean Smoothness All Timeframes       Miyagawa Targets Of Ewsr1 Ets Fusions Dn                    -0.557323  0.0171544    -1.78647  0.148577   1                8567   0.54386   MIYAGAWA_TARGETS_OF_EWSR1_ETS_FUSIONS_DN
Irregularity                         Steger Adipogenesis Up                                       0.767694  0.0070922     1.82794  0.218099   1                 978   0.538462  STEGER_ADIPOGENESIS_UP
Variation Smoothness All Timeframes  Steger Adipogenesis Up                                      -0.710437  0.0263799    -1.69679  0.219274   1                9731   0.538462  STEGER_ADIPOGENESIS_UP
Mean Smoothness All Timeframes       Steger Adipogenesis Dn                                      -0.698886  0.0902469    -1.53951  0.231366   1                8351   0.857143  STEGER_ADIPOGENESIS_DN
Mean Smoothness All Timeframes       Zhu Skil Targets Up                                         -0.67064   0.0367647    -1.63939  0.194025   1               10319   0.421053  ZHU_SKIL_TARGETS_UP
Mean Smoothness All Timeframes       Issaeva Mll2 Targets                                        -0.566326  0.0458476    -1.61491  0.200977   1                8019   0.615385  ISSAEVA_MLL2_TARGETS
Mean Smoothness All Timeframes       Pasini Suz12 Targets Dn                                     -0.507335  0.0256462    -1.73003  0.16069    1                7890   0.564912  PASINI_SUZ12_TARGETS_DN
Variation Smoothness All Timeframes  Sumi Hnf4A Targets                                          -0.640269  0.0187438    -1.65602  0.236572   1                7407   0.75      SUMI_HNF4A_TARGETS
Mean Smoothness All Timeframes       Sumi Hnf4A Targets                                          -0.631734  0.0172621    -1.63208  0.194532   1                9360   0.5       SUMI_HNF4A_TARGETS
Mean Smoothness All Timeframes       Azare Neoplastic Transformation By Stat3 Dn                 -0.61767   0.0574599    -1.52112  0.242693   1                8945   0.545455  AZARE_NEOPLASTIC_TRANSFORMATION_BY_STAT3_DN
Mean Smoothness All Timeframes       Wiederschain Targets Of Bmi1 And Pcgf2                      -0.544871  0.0586117    -1.57379  0.216793   1                8358   0.543478  WIEDERSCHAIN_TARGETS_OF_BMI1_AND_PCGF2
Variation Smoothness Uptake          Bilanges Rapamycin Sensitive Via Tsc1 And Tsc2               0.593023  0.0189469     1.76315  0.192501   1                2344   0.609375  BILANGES_RAPAMYCIN_SENSITIVE_VIA_TSC1_AND_TSC2
Washout                              Kim Tial1 Targets                                            0.616972  0.0064843     1.81562  0.206704   1                1781   0.5       KIM_TIAL1_TARGETS
Top Late Enhancement                 Kim Tial1 Targets                                           -0.617424  0.00549115   -1.80987  0.213855   1                9090   0.571429  KIM_TIAL1_TARGETS
Mean Smoothness All Timeframes       Torchia Targets Of Ewsr1 Fli1 Fusion Up                     -0.352464  0.0217218    -1.50898  0.249559   1                8049   0.40625   TORCHIA_TARGETS_OF_EWSR1_FLI1_FUSION_UP
Variation Smoothness All Timeframes  Winzen Degraded Via Khsrp                                   -0.513511  0.0434521    -1.60589  0.247184   1                8987   0.461538  WINZEN_DEGRADED_VIA_KHSRP
Mean Smoothness All Timeframes       Winzen Degraded Via Khsrp                                   -0.552624  0.017593     -1.72356  0.16069    1                9208   0.476923  WINZEN_DEGRADED_VIA_KHSRP
Variation Smoothness All Timeframes  Ikeda Mir30 Targets Up                                      -0.478915  0.0242866    -1.65395  0.236572   1                7753   0.495327  IKEDA_MIR30_TARGETS_UP
Mean Smoothness All Timeframes       Ikeda Mir30 Targets Up                                      -0.489284  0.0166968    -1.69576  0.169834   1                7315   0.598131  IKEDA_MIR30_TARGETS_UP
Mean Smoothness All Timeframes       Servitja Islet Hnf1A Targets Up                             -0.520147  0.0793491    -1.56752  0.21886    1                9167   0.483333  SERVITJA_ISLET_HNF1A_TARGETS_UP
Variation Smoothness All Timeframes  Servitja Liver Hnf1A Targets Dn                             -0.48628   0.00321027   -1.76857  0.193127   1                9484   0.352941  SERVITJA_LIVER_HNF1A_TARGETS_DN
Mean Smoothness All Timeframes       Servitja Liver Hnf1A Targets Dn                             -0.479359  0.0045835    -1.74895  0.153387   1                9630   0.352941  SERVITJA_LIVER_HNF1A_TARGETS_DN
Irregularity                         Madan Dppa4 Targets                                          0.68407   0.00216535    1.85134  0.218099   1                 386   0.25      MADAN_DPPA4_TARGETS
Variation Smoothness All Timeframes  Kohoutek Ccnt2 Targets                                      -0.481074  0.0082495    -1.66007  0.236572   1                9060   0.419355  KOHOUTEK_CCNT2_TARGETS
Mean Smoothness All Timeframes       Kohoutek Ccnt2 Targets                                      -0.478895  0.00999201   -1.65778  0.18313    1                9807   0.387097  KOHOUTEK_CCNT2_TARGETS
Variation Smoothness All Timeframes  Pedersen Metastasis By Erbb2 Isoform 1                      -0.578502  0.0535997    -1.61064  0.244212   1                8766   0.472222  PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_1
Mean Smoothness All Timeframes       Pedersen Metastasis By Erbb2 Isoform 1                      -0.583423  0.0504487    -1.61495  0.200977   1                9109   0.444444  PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_1
Variation Smoothness All Timeframes  Pedersen Metastasis By Erbb2 Isoform 6                      -0.586328  0.026629     -1.63479  0.236712   1                8245   0.681818  PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_6
Mean Smoothness All Timeframes       Pedersen Metastasis By Erbb2 Isoform 6                      -0.576449  0.0312437    -1.61551  0.200977   1                8764   0.636364  PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_6
Mean Smoothness All Timeframes       Plasari Tgfb1 Targets 1Hr Up                                -0.599558  0.0822408    -1.54545  0.227613   1                9709   0.482759  PLASARI_TGFB1_TARGETS_1HR_UP
Mean Smoothness All Timeframes       Plasari Tgfb1 Targets 10Hr Up                               -0.471051  0.0561122    -1.56449  0.219929   1                8349   0.479167  PLASARI_TGFB1_TARGETS_10HR_UP
Irregularity                         Plasari Tgfb1 Targets 10Hr Dn                                0.487449  0.00341365    1.84946  0.218099   1                3153   0.463542  PLASARI_TGFB1_TARGETS_10HR_DN
Variation Smoothness All Timeframes  Plasari Tgfb1 Targets 10Hr Dn                               -0.467346  0.00588355   -1.78237  0.193127   1                7929   0.557292  PLASARI_TGFB1_TARGETS_10HR_DN
Mean Smoothness All Timeframes       Plasari Tgfb1 Targets 10Hr Dn                               -0.450058  0.0105434    -1.71925  0.16069    1                8470   0.494792  PLASARI_TGFB1_TARGETS_10HR_DN
Variation Smoothness All Timeframes  Plasari Nfic Targets Basal Up                               -0.634016  0.017762     -1.65768  0.236572   1                7347   0.857143  PLASARI_NFIC_TARGETS_BASAL_UP
Mean Smoothness All Timeframes       Plasari Nfic Targets Basal Up                               -0.651352  0.0126432    -1.70154  0.166496   1                7444   0.857143  PLASARI_NFIC_TARGETS_BASAL_UP
Irregularity                         Plasari Tgfb1 Signaling Via Nfic 1Hr Dn                      0.526802  0.0013785     1.85338  0.218099   1                2531   0.460674  PLASARI_TGFB1_SIGNALING_VIA_NFIC_1HR_DN
Variation Smoothness All Timeframes  Plasari Tgfb1 Signaling Via Nfic 1Hr Dn                     -0.528314  0.000998004  -1.85831  0.168147   1                7601   0.629213  PLASARI_TGFB1_SIGNALING_VIA_NFIC_1HR_DN
Mean Smoothness All Timeframes       Plasari Tgfb1 Signaling Via Nfic 1Hr Dn                     -0.539035  0.000401284  -1.89666  0.117961   1                8250   0.561798  PLASARI_TGFB1_SIGNALING_VIA_NFIC_1HR_DN
Variation Smoothness All Timeframes  Plasari Tgfb1 Signaling Via Nfic 10Hr Up                    -0.600427  0.00239952   -1.86291  0.168147   1                9841   0.44186   PLASARI_TGFB1_SIGNALING_VIA_NFIC_10HR_UP
Mean Smoothness All Timeframes       Plasari Tgfb1 Signaling Via Nfic 10Hr Up                    -0.599816  0.00160772   -1.857    0.129188   1                9772   0.44186   PLASARI_TGFB1_SIGNALING_VIA_NFIC_10HR_UP
Variation Smoothness All Timeframes  Wang Mll Targets                                            -0.530967  0.00722456   -1.833    0.181077   1                9137   0.420455  WANG_MLL_TARGETS
Mean Smoothness All Timeframes       Wang Mll Targets                                            -0.543979  0.00561235   -1.8783   0.117961   1                9202   0.443182  WANG_MLL_TARGETS
Variation Smoothness All Timeframes  Delacroix Rar Targets Dn                                    -0.59342   0.0181855    -1.65188  0.236572   1               10676   0.315789  DELACROIX_RAR_TARGETS_DN
Mean Smoothness All Timeframes       Delacroix Rar Targets Dn                                    -0.607033  0.0130559    -1.68763  0.172248   1               10729   0.315789  DELACROIX_RAR_TARGETS_DN
Mean Smoothness All Timeframes       Liu Il13 Priming Model                                      -0.662298  0.0217563    -1.68115  0.175441   1                9109   0.538462  LIU_IL13_PRIMING_MODEL
Washout                              Fu Interact With Alkbh8                                      0.803507  0.00281407    1.80071  0.217568   1                 948   0.769231  FU_INTERACT_WITH_ALKBH8
Top Late Enhancement                 Fu Interact With Alkbh8                                     -0.785817  0.00579884   -1.75474  0.213855   1                9896   0.769231  FU_INTERACT_WITH_ALKBH8
Irregularity                         Guillaumond Klf10 Targets Up                                 0.536966  0.00141015    1.87679  0.203797   1                2780   0.457143  GUILLAUMOND_KLF10_TARGETS_UP
Variation Smoothness All Timeframes  Pedrioli Mir31 Targets Dn                                   -0.434895  0.00377734   -1.74854  0.194054   1                8049   0.445283  PEDRIOLI_MIR31_TARGETS_DN
Mean Smoothness All Timeframes       Pedrioli Mir31 Targets Dn                                   -0.420265  0.00597967   -1.68658  0.172248   1                8689   0.377358  PEDRIOLI_MIR31_TARGETS_DN
Variation Smoothness All Timeframes  Zhang Adipogenesis By Bmp7                                  -0.662847  0.0022115    -1.89221  0.168147   1                8987   0.461538  ZHANG_ADIPOGENESIS_BY_BMP7
Mean Smoothness All Timeframes       Zhang Adipogenesis By Bmp7                                  -0.633061  0.0056191    -1.80189  0.141179   1                9217   0.461538  ZHANG_ADIPOGENESIS_BY_BMP7
Variation Smoothness All Timeframes  Foster Kdm1A Targets Up                                     -0.429426  0.0152843    -1.63634  0.236712   1                8390   0.392857  FOSTER_KDM1A_TARGETS_UP
Mean Smoothness All Timeframes       Foster Kdm1A Targets Up                                     -0.413134  0.0223124    -1.57835  0.216178   1                9029   0.348214  FOSTER_KDM1A_TARGETS_UP
Variation Smoothness All Timeframes  Pece Mammary Stem Cell Up                                    0.623318  0.0181138     1.81838  0.237345   1                1324   0.54918   PECE_MAMMARY_STEM_CELL_UP
Variation Smoothness Uptake          Pece Mammary Stem Cell Up                                    0.622788  0.0166145     1.82804  0.173515   1                1399   0.565574  PECE_MAMMARY_STEM_CELL_UP
Mean Smoothness All Timeframes       Pece Mammary Stem Cell Dn                                   -0.442313  0.0135655    -1.65769  0.18313    1                8103   0.465116  PECE_MAMMARY_STEM_CELL_DN
Irregularity                         Kumar Autophagy Network                                      0.533104  0.000199481   2.0212   0.12903    0.531219         2030   0.348837  KUMAR_AUTOPHAGY_NETWORK
Irregularity                         Acevedo Fgfr1 Targets In Prostate Cancer Model Dn            0.519284  0.00798722    1.83528  0.218099   1                2218   0.410811  ACEVEDO_FGFR1_TARGETS_IN_PROSTATE_CANCER_MODEL_DN
Variation Smoothness All Timeframes  Acevedo Fgfr1 Targets In Prostate Cancer Model Dn           -0.52418   0.00518341   -1.85467  0.168479   1                8237   0.524324  ACEVEDO_FGFR1_TARGETS_IN_PROSTATE_CANCER_MODEL_DN
Mean Smoothness All Timeframes       Acevedo Fgfr1 Targets In Prostate Cancer Model Dn           -0.551567  0.00238095   -1.96034  0.109401   1                8942   0.475676  ACEVEDO_FGFR1_TARGETS_IN_PROSTATE_CANCER_MODEL_DN
Mean Smoothness All Timeframes       Lim Mammary Luminal Mature Dn                               -0.658452  0.0249497    -1.72103  0.16069    1                9012   0.682927  LIM_MAMMARY_LUMINAL_MATURE_DN
Variation Smoothness All Timeframes  Durand Stroma S Up                                          -0.499703  0.0116536    -1.80544  0.192967   1                7974   0.542453  DURAND_STROMA_S_UP
Mean Smoothness All Timeframes       Durand Stroma S Up                                          -0.497033  0.0118236    -1.79666  0.141179   1                8043   0.54717   DURAND_STROMA_S_UP
Ld Late Lt0                          Smirnov Response To Ir 2Hr Up                               -0.643554  0.00121507   -2.01599  0.210341   1               10549   0.333333  SMIRNOV_RESPONSE_TO_IR_2HR_UP
Mean Smoothness All Timeframes       Ghandhi Direct Irradiation Dn                               -0.536981  0.0362089    -1.5928   0.210308   1                8095   0.521739  GHANDHI_DIRECT_IRRADIATION_DN
Top Late Enhancement                 Zhou Cell Cycle Genes In Ir Response 24Hr                   -0.689335  0.0372093    -1.7394   0.213855   1               10197   0.561905  ZHOU_CELL_CYCLE_GENES_IN_IR_RESPONSE_24HR
Variation Smoothness All Timeframes  Zwang Class 2 Transiently Induced By Egf                    -0.601626  0.037336     -1.64981  0.236572   1                7930   0.625     ZWANG_CLASS_2_TRANSIENTLY_INDUCED_BY_EGF
Mean Smoothness All Timeframes       Zwang Class 2 Transiently Induced By Egf                    -0.646394  0.0116184    -1.76074  0.153387   1                8852   0.5625    ZWANG_CLASS_2_TRANSIENTLY_INDUCED_BY_EGF
Mean Smoothness All Timeframes       Zwang Class 3 Transiently Induced By Egf                    -0.451661  0.0776119    -1.53073  0.235036   1                8550   0.373626  ZWANG_CLASS_3_TRANSIENTLY_INDUCED_BY_EGF
Variation Smoothness All Timeframes  Eppert Hsc R                                                -0.485149  0.010198     -1.72571  0.211541   1                7849   0.568421  EPPERT_HSC_R
Mean Smoothness All Timeframes       Eppert Hsc R                                                -0.490127  0.00901623   -1.74439  0.155059   1                9347   0.389474  EPPERT_HSC_R
Variation Smoothness All Timeframes  Eppert Ce Hsc Lsc                                           -0.589425  0.0215225    -1.6831   0.231284   1                7942   0.689655  EPPERT_CE_HSC_LSC
Mean Smoothness All Timeframes       Eppert Ce Hsc Lsc                                           -0.635845  0.00579421   -1.82069  0.133542   1                9363   0.517241  EPPERT_CE_HSC_LSC

</div>

Pathways (c2.cp) from MSigDB.


```python
plot_ds(df_cp_rv, fdr=0.25)
```

![](figures/gsea_cp-es-plot-rv_1){#cp-es-plot-rv }\



```python
table_ds(df_cp_rv, fdr=0.25)
```


<div class="datatable">mri_feature           gene_set                                                                                                                 es           p      nes       fdr    fwer    max_es_at    le_prop  gene_set_code
--------------------  -----------------------------------------------------------------------------------------------------------------  --------  ----------  -------  --------  ------  -----------  ---------  -----------------------------------------------------------------------------------------------------------------
Top Late Enhancement  Kegg Proteasome                                                                                                    0.604507  0.00738941  1.76684  0.228139       1         3437   0.736842  KEGG_PROTEASOME
Top Late Enhancement  Biocarta Salmonella Pathway                                                                                        0.7041    0.00670712  1.71326  0.228139       1          966   0.583333  BIOCARTA_SALMONELLA_PATHWAY
Top Late Enhancement  Biocarta Proteasome Pathway                                                                                        0.684692  0.00246837  1.93668  0.228139       1         3076   0.851852  BIOCARTA_PROTEASOME_PATHWAY
Top Late Enhancement  Biocarta Actiny Pathway                                                                                            0.697663  0.00133142  1.79405  0.228139       1          966   0.529412  BIOCARTA_ACTINY_PATHWAY
Top Late Enhancement  Reactome Cross Presentation Of Soluble Exogenous Antigens Endosomes                                                0.581219  0.0108641   1.70627  0.228139       1         3437   0.7       REACTOME_CROSS_PRESENTATION_OF_SOLUBLE_EXOGENOUS_ANTIGENS_ENDOSOMES
Top Late Enhancement  Reactome Spry Regulation Of Fgf Signaling                                                                          0.680575  0.00321477  1.68243  0.228139       1         2507   0.461538  REACTOME_SPRY_REGULATION_OF_FGF_SIGNALING
Top Late Enhancement  Reactome Er Phagosome Pathway                                                                                      0.56234   0.0121486   1.6877   0.228139       1         3692   0.648148  REACTOME_ER_PHAGOSOME_PATHWAY
Top Late Enhancement  Reactome P53 Independent G1 S Dna Damage Checkpoint                                                                0.604889  0.00522088  1.77814  0.228139       1         3437   0.744186  REACTOME_P53_INDEPENDENT_G1_S_DNA_DAMAGE_CHECKPOINT
Top Late Enhancement  Reactome Cdk Mediated Phosphorylation And Removal Of Cdc6                                                          0.594129  0.00703871  1.74521  0.228139       1         3437   0.714286  REACTOME_CDK_MEDIATED_PHOSPHORYLATION_AND_REMOVAL_OF_CDC6
Top Late Enhancement  Reactome Regulation Of Ornithine Decarboxylase Odc                                                                 0.562853  0.0121535   1.66086  0.238874       1         3437   0.674419  REACTOME_REGULATION_OF_ORNITHINE_DECARBOXYLASE_ODC
Top Late Enhancement  Reactome Regulation Of Apoptosis                                                                                   0.567539  0.0065052   1.69173  0.228139       1         3437   0.66      REACTOME_REGULATION_OF_APOPTOSIS
Top Late Enhancement  Reactome P53 Dependent G1 Dna Damage Response                                                                      0.566718  0.00772085  1.68287  0.228139       1         3437   0.6875    REACTOME_P53_DEPENDENT_G1_DNA_DAMAGE_RESPONSE
Top Late Enhancement  Reactome Prefoldin Mediated Transfer Of Substrate To Cct Tric                                                      0.642802  0.00624488  1.71903  0.228139       1         1817   0.55      REACTOME_PREFOLDIN_MEDIATED_TRANSFER_OF_SUBSTRATE_TO_CCT_TRIC
Top Late Enhancement  Reactome Formation Of Tubulin Folding Intermediates By Cct Tric                                                    0.709025  0.00391162  1.80898  0.228139       1         1507   0.6       REACTOME_FORMATION_OF_TUBULIN_FOLDING_INTERMEDIATES_BY_CCT_TRIC
Top Late Enhancement  Reactome Cdt1 Association With The Cdc6 Orc Origin Complex                                                         0.567759  0.00781563  1.67823  0.228139       1         3437   0.680851  REACTOME_CDT1_ASSOCIATION_WITH_THE_CDC6_ORC_ORIGIN_COMPLEX
Top Late Enhancement  Reactome Autodegradation Of The E3 Ubiquitin Ligase Cop1                                                           0.575319  0.00984925  1.69919  0.228139       1         3437   0.697674  REACTOME_AUTODEGRADATION_OF_THE_E3_UBIQUITIN_LIGASE_COP1
Top Late Enhancement  Reactome Regulation Of Mitotic Cell Cycle                                                                          0.557018  0.00620062  1.67558  0.228139       1         3262   0.623188  REACTOME_REGULATION_OF_MITOTIC_CELL_CYCLE
Top Late Enhancement  Reactome Assembly Of The Pre Replicative Complex                                                                   0.564952  0.00680681  1.6893   0.228139       1         3437   0.690909  REACTOME_ASSEMBLY_OF_THE_PRE_REPLICATIVE_COMPLEX
Top Late Enhancement  Reactome Destabilization Of Mrna By Auf1 Hnrnp D0                                                                  0.572738  0.00742972  1.69995  0.228139       1         3437   0.695652  REACTOME_DESTABILIZATION_OF_MRNA_BY_AUF1_HNRNP_D0
Top Late Enhancement  Reactome Apc C Cdh1 Mediated Degradation Of Cdc20 And Other Apc C Cdh1 Targeted Proteins In Late Mitosis Early G1  0.564667  0.0060042   1.68625  0.228139       1         3437   0.649123  REACTOME_APC_C_CDH1_MEDIATED_DEGRADATION_OF_CDC20_AND_OTHER_APC_C_CDH1_TARGETED_PROTEINS_IN_LATE_MITOSIS_EARLY_G1
Top Late Enhancement  Reactome Apc C Cdc20 Mediated Degradation Of Mitotic Proteins                                                      0.566883  0.00560504  1.70335  0.228139       1         3437   0.644068  REACTOME_APC_C_CDC20_MEDIATED_DEGRADATION_OF_MITOTIC_PROTEINS
Top Late Enhancement  Reactome Scf Beta Trcp Mediated Degradation Of Emi1                                                                0.561254  0.0123333   1.65976  0.238874       1         3437   0.688889  REACTOME_SCF_BETA_TRCP_MEDIATED_DEGRADATION_OF_EMI1
Top Late Enhancement  Reactome Scfskp2 Mediated Degradation Of P27 P21                                                                   0.568042  0.00821561  1.67845  0.228139       1         3437   0.6875    REACTOME_SCFSKP2_MEDIATED_DEGRADATION_OF_P27_P21
Top Late Enhancement  Reactome Vif Mediated Degradation Of Apobec3G                                                                      0.5713    0.00914757  1.6937   0.228139       1         3947   0.755556  REACTOME_VIF_MEDIATED_DEGRADATION_OF_APOBEC3G

</div>
