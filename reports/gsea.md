---
title: Gene Set Enrichment Analysis of CAD features
author: Tycho Bismeijer
date: 2017-04-12
---

## Setup ## {.collapsed}

Load libraries.


```python
import sys

from IPython.display import display, Markdown
import numpy as np
import pandas as pd
from tabulate import tabulate
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
    ds = xr.open_dataset(fn)
    ds['mri_feature'] = [s.decode() for s in ds['mri_feature'].values]
    ds['mri_feature'] = [" ".join(s.split('_')).title()
                         for s in ds['mri_feature'].values]
    ds['gene_set_code'] = ('gene_set',
                           [s.decode() for s in ds['gene_set'].values])
    ds['gene_set'] = [" ".join(s.split('_')).title()
                      for s in ds['gene_set_code'].values]
    return ds
df_h = load_gsea_ds("../analyses/gsea/mri-features_h.all_T.nc")
```



A function to plot a heatmap with normalized enrichment statistics (NES) that
are significant.


```python
def plot_ds(ds, fdr, le_prop=0.0):
    ds = ds.copy()
    ds['significance_mask'] = (ds['fdr'] > fdr) | (ds['le_prop'] < le_prop)
    ds = ds.sel(
        gene_set=np.logical_not(ds['significance_mask'])
                 .sum('mri_feature') > 0,
    )
    with plot.subplots(1, 1, figsize=(8, 5.5), dpi=192) as (fig, ax):
        plot.heatmap(ds['nes'], mask=ds['significance_mask'], ax=ax)
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

![](figures/gsea_hallmarks-es-plot_1.png){#hallmarks-es-plot }\



```python
table_ds(df_h, fdr=0.25)
```


<div class="datatable">mri_feature                 gene_set                                          es           p      nes         fdr       fwer    max_es_at    le_prop  gene_set_code
--------------------------  ------------------------------------------  --------  ----------  -------  ----------  ---------  -----------  ---------  ------------------------------------------
Largest Diameter            Hallmark Tnfa Signaling Via Nfkb            0.423123  0.069593    1.31066  0.17295     1                 2662   0.380952  HALLMARK_TNFA_SIGNALING_VIA_NFKB
Variance Vox Val            Hallmark Tnfa Signaling Via Nfkb            0.424242  0.0659934   1.31555  0.210306    1                 3380   0.440476  HALLMARK_TNFA_SIGNALING_VIA_NFKB
Ld Init Enhancement Gt100   Hallmark Tnfa Signaling Via Nfkb            0.431836  0.0577942   1.33905  0.155851    1                 3595   0.470238  HALLMARK_TNFA_SIGNALING_VIA_NFKB
Ld Late Lt0                 Hallmark Tnfa Signaling Via Nfkb            0.401436  0.112289    1.24523  0.221385    1                 4354   0.52381   HALLMARK_TNFA_SIGNALING_VIA_NFKB
Irregularity                Hallmark Hypoxia                            0.436426  0.00959904  1.34224  0.190522    0.479952          3576   0.490909  HALLMARK_HYPOXIA
Largest Diameter            Hallmark Hypoxia                            0.47027   0.00209979  1.44506  0.103834    0.10499           3847   0.551515  HALLMARK_HYPOXIA
Ld Init Enhancement Gt100   Hallmark Hypoxia                            0.456109  0.00339966  1.40084  0.139897    0.169983          4145   0.569697  HALLMARK_HYPOXIA
Ld Late Lt0                 Hallmark Hypoxia                            0.441323  0.00769923  1.35728  0.204204    0.384962          4529   0.593939  HALLMARK_HYPOXIA
Largest Diameter            Hallmark Cholesterol Homeostasis            0.42097   0.0808919   1.23951  0.19224     1                 4959   0.615385  HALLMARK_CHOLESTEROL_HOMEOSTASIS
Ld Init Enhancement Gt100   Hallmark Cholesterol Homeostasis            0.426082  0.0682932   1.25588  0.193519    1                 4229   0.538462  HALLMARK_CHOLESTEROL_HOMEOSTASIS
Ld Late Lt0                 Hallmark Cholesterol Homeostasis            0.432656  0.0579942   1.27565  0.221385    1                 4319   0.538462  HALLMARK_CHOLESTEROL_HOMEOSTASIS
Irregularity                Hallmark Mitotic Spindle                    0.468216  0.00439956  1.45792  0.176628    0.219978          3192   0.416667  HALLMARK_MITOTIC_SPINDLE
Volume                      Hallmark Mitotic Spindle                    0.458632  0.00739926  1.42589  0.164462    0.369963          1315   0.255208  HALLMARK_MITOTIC_SPINDLE
Largest Diameter            Hallmark Mitotic Spindle                    0.479571  0.00209979  1.49469  0.0986456   0.10499           2210   0.317708  HALLMARK_MITOTIC_SPINDLE
Variance Vox Val            Hallmark Mitotic Spindle                    0.449459  0.0115988   1.39672  0.174099    0.579942          2502   0.338542  HALLMARK_MITOTIC_SPINDLE
Vol Init Enhancement Gt100  Hallmark Mitotic Spindle                    0.454314  0.00949905  1.41064  0.162916    0.474953          1311   0.255208  HALLMARK_MITOTIC_SPINDLE
Ld Init Enhancement Gt100   Hallmark Mitotic Spindle                    0.482646  0.00179982  1.50363  0.114258    0.089991          2068   0.307292  HALLMARK_MITOTIC_SPINDLE
Ld Late Lt0                 Hallmark Mitotic Spindle                    0.488948  0.00139986  1.51967  0.103221    0.069993          1274   0.255208  HALLMARK_MITOTIC_SPINDLE
Irregularity                Hallmark Tgf Beta Signaling                 0.481729  0.0148015   1.41513  0.180171    0.740074          3252   0.458333  HALLMARK_TGF_BETA_SIGNALING
Largest Diameter            Hallmark Tgf Beta Signaling                 0.467622  0.0244098   1.37393  0.152142    1                 2087   0.375     HALLMARK_TGF_BETA_SIGNALING
Variance Vox Val            Hallmark Tgf Beta Signaling                 0.488881  0.0136014   1.43143  0.15503     0.680068          3126   0.479167  HALLMARK_TGF_BETA_SIGNALING
Ld Init Enhancement Gt100   Hallmark Tgf Beta Signaling                 0.457773  0.0387039   1.34395  0.155851    1                 4020   0.5625    HALLMARK_TGF_BETA_SIGNALING
Ld Late Lt0                 Hallmark Tgf Beta Signaling                 0.438461  0.0659      1.28801  0.221385    1                 3444   0.479167  HALLMARK_TGF_BETA_SIGNALING
Largest Diameter            Hallmark Il6 Jak Stat3 Signaling            0.426439  0.13384     1.24668  0.19224     1                 3686   0.442623  HALLMARK_IL6_JAK_STAT3_SIGNALING
Ld Init Enhancement Gt100   Hallmark Il6 Jak Stat3 Signaling            0.413739  0.170134    1.20869  0.22304     1                 3612   0.42623   HALLMARK_IL6_JAK_STAT3_SIGNALING
Largest Diameter            Hallmark Dna Repair                         0.413265  0.0231977   1.28972  0.183985    1                 4381   0.514925  HALLMARK_DNA_REPAIR
Ld Init Enhancement Gt100   Hallmark Dna Repair                         0.406111  0.0327967   1.26606  0.193519    1                 3960   0.470149  HALLMARK_DNA_REPAIR
Ld Late Lt0                 Hallmark Dna Repair                         0.398047  0.0452955   1.2423   0.221385    1                 3969   0.470149  HALLMARK_DNA_REPAIR
Irregularity                Hallmark G2M Checkpoint                     0.65688   0.00069993  2.02283  0.0144105   0.0349965         2008   0.6       HALLMARK_G2M_CHECKPOINT
Volume                      Hallmark G2M Checkpoint                     0.659236  0.00049995  2.02675  0.0154117   0.0249975         1991   0.594444  HALLMARK_G2M_CHECKPOINT
Largest Diameter            Hallmark G2M Checkpoint                     0.661803  0.00029997  2.03293  0.0135604   0.0149985         1274   0.522222  HALLMARK_G2M_CHECKPOINT
Mean Vox Val                Hallmark G2M Checkpoint                     0.63622   0.00189981  1.95362  0.0278705   0.0949905         1209   0.483333  HALLMARK_G2M_CHECKPOINT
Variance Vox Val            Hallmark G2M Checkpoint                     0.63118   0.0019998   1.94187  0.0296225   0.09999           1454   0.5       HALLMARK_G2M_CHECKPOINT
Top Late Enhancement        Hallmark G2M Checkpoint                     0.54457   0.0310969   1.66838  0.181943    1                 2006   0.444444  HALLMARK_G2M_CHECKPOINT
Vol Init Enhancement Gt100  Hallmark G2M Checkpoint                     0.657133  0.00069993  2.0199   0.0169633   0.0349965         2029   0.6       HALLMARK_G2M_CHECKPOINT
Ld Init Enhancement Gt100   Hallmark G2M Checkpoint                     0.659652  0.00069993  2.02288  0.0153616   0.0349965         1056   0.5       HALLMARK_G2M_CHECKPOINT
Vol Late Lt0                Hallmark G2M Checkpoint                     0.642938  0.00219978  1.97364  0.0219185   0.109989          1726   0.55      HALLMARK_G2M_CHECKPOINT
Ld Late Lt0                 Hallmark G2M Checkpoint                     0.676242  0.00039996  2.08     0.00930702  0.019998          1589   0.566667  HALLMARK_G2M_CHECKPOINT
Largest Diameter            Hallmark Apoptosis                          0.427161  0.00919908  1.32072  0.17295     0.459954          4880   0.627737  HALLMARK_APOPTOSIS
Vol Init Enhancement Gt100  Hallmark Apoptosis                          0.443099  0.00619938  1.36873  0.195113    0.309969          3048   0.437956  HALLMARK_APOPTOSIS
Ld Init Enhancement Gt100   Hallmark Apoptosis                          0.438886  0.00389961  1.35813  0.155851    0.194981          3303   0.481752  HALLMARK_APOPTOSIS
Ld Late Lt0                 Hallmark Apoptosis                          0.420238  0.0128987   1.2983   0.221385    0.644936          4255   0.569343  HALLMARK_APOPTOSIS
Largest Diameter            Hallmark Adipogenesis                       0.393601  0.030497    1.23419  0.19224     1                 3604   0.432432  HALLMARK_ADIPOGENESIS
Ld Init Enhancement Gt100   Hallmark Adipogenesis                       0.401419  0.019798    1.25854  0.193519    0.989901          3635   0.443243  HALLMARK_ADIPOGENESIS
Ld Late Lt0                 Hallmark Adipogenesis                       0.399633  0.0231977   1.25303  0.221385    1                 4249   0.502703  HALLMARK_ADIPOGENESIS
Largest Diameter            Hallmark Myogenesis                         0.432089  0.0362964   1.3232   0.17295     1                 3888   0.496063  HALLMARK_MYOGENESIS
Variance Vox Val            Hallmark Myogenesis                         0.426804  0.049795    1.30018  0.21883     1                 3970   0.519685  HALLMARK_MYOGENESIS
Ld Init Enhancement Gt100   Hallmark Myogenesis                         0.439666  0.0276972   1.34547  0.155851    1                 4634   0.582677  HALLMARK_MYOGENESIS
Variance Vox Val            Hallmark Interferon Alpha Response          0.51592   0.091475    1.60654  0.0884871   1                 2360   0.5       HALLMARK_INTERFERON_ALPHA_RESPONSE
Variance Vox Val            Hallmark Interferon Gamma Response          0.487252  0.0683932   1.54007  0.107014    1                 3332   0.547619  HALLMARK_INTERFERON_GAMMA_RESPONSE
Ld Init Enhancement Gt100   Hallmark Unfolded Protein Response          0.39327   0.0716928   1.21331  0.22304     1                 3165   0.376147  HALLMARK_UNFOLDED_PROTEIN_RESPONSE
Ld Late Lt0                 Hallmark Unfolded Protein Response          0.40635   0.0448955   1.25401  0.221385    1                 2967   0.376147  HALLMARK_UNFOLDED_PROTEIN_RESPONSE
Irregularity                Hallmark Mtorc1 Signaling                   0.438523  0.0113989   1.37975  0.180171    0.569943          2936   0.417989  HALLMARK_MTORC1_SIGNALING
Volume                      Hallmark Mtorc1 Signaling                   0.480353  0.00159984  1.50355  0.135653    0.079992          1996   0.375661  HALLMARK_MTORC1_SIGNALING
Largest Diameter            Hallmark Mtorc1 Signaling                   0.467799  0.00339966  1.4672   0.0986456   0.169983          3089   0.460317  HALLMARK_MTORC1_SIGNALING
Vol Init Enhancement Gt100  Hallmark Mtorc1 Signaling                   0.480969  0.00169983  1.50483  0.13664     0.0849915         2118   0.386243  HALLMARK_MTORC1_SIGNALING
Ld Init Enhancement Gt100   Hallmark Mtorc1 Signaling                   0.469849  0.00319968  1.47288  0.115854    0.159984          3484   0.502646  HALLMARK_MTORC1_SIGNALING
Vol Late Lt0                Hallmark Mtorc1 Signaling                   0.461441  0.00359964  1.44465  0.197984    0.179982          2205   0.375661  HALLMARK_MTORC1_SIGNALING
Ld Late Lt0                 Hallmark Mtorc1 Signaling                   0.494309  0.00119988  1.55195  0.0993749   0.059994          2967   0.470899  HALLMARK_MTORC1_SIGNALING
Irregularity                Hallmark E2F Targets                        0.661312  0.00169983  2.05076  0.0144105   0.0849915         2165   0.65      HALLMARK_E2F_TARGETS
Volume                      Hallmark E2F Targets                        0.683821  0.00039996  2.11392  0.0154117   0.019998          1597   0.611111  HALLMARK_E2F_TARGETS
Largest Diameter            Hallmark E2F Targets                        0.695368  0.00029997  2.15005  0.0107082   0.0149985         1717   0.638889  HALLMARK_E2F_TARGETS
Mean Vox Val                Hallmark E2F Targets                        0.62874   0.0047      1.94435  0.0278705   0.235             2376   0.638889  HALLMARK_E2F_TARGETS
Variance Vox Val            Hallmark E2F Targets                        0.645477  0.00369963  2.00222  0.0296225   0.184982          2846   0.694444  HALLMARK_E2F_TARGETS
Top Late Enhancement        Hallmark E2F Targets                        0.531999  0.0503      1.63907  0.181943    1                 2788   0.516667  HALLMARK_E2F_TARGETS
Vol Init Enhancement Gt100  Hallmark E2F Targets                        0.666194  0.00129987  2.05892  0.0169633   0.0649935         1586   0.594444  HALLMARK_E2F_TARGETS
Ld Init Enhancement Gt100   Hallmark E2F Targets                        0.685042  0.0006      2.11588  0.0138104   0.03              1634   0.611111  HALLMARK_E2F_TARGETS
Vol Late Lt0                Hallmark E2F Targets                        0.670707  0.00129987  2.0708   0.021318    0.0649935         1934   0.622222  HALLMARK_E2F_TARGETS
Ld Late Lt0                 Hallmark E2F Targets                        0.697474  0.0003      2.15791  0.00930702  0.015             1589   0.616667  HALLMARK_E2F_TARGETS
Irregularity                Hallmark Myc Targets V1                     0.49056   0.0129987   1.58217  0.100693    0.649935          3092   0.507772  HALLMARK_MYC_TARGETS_V1
Volume                      Hallmark Myc Targets V1                     0.536135  0.00249975  1.72294  0.050318    0.124988          2563   0.544041  HALLMARK_MYC_TARGETS_V1
Largest Diameter            Hallmark Myc Targets V1                     0.548356  0.00149985  1.76419  0.0349868   0.0749925         3059   0.595855  HALLMARK_MYC_TARGETS_V1
Mean Vox Val                Hallmark Myc Targets V1                     0.44926   0.0349965   1.44783  0.247161    1                 2571   0.409326  HALLMARK_MYC_TARGETS_V1
Variance Vox Val            Hallmark Myc Targets V1                     0.410154  0.0805919   1.32425  0.210306    1                 1979   0.310881  HALLMARK_MYC_TARGETS_V1
Vol Init Enhancement Gt100  Hallmark Myc Targets V1                     0.526336  0.00339966  1.68904  0.056264    0.169983          3178   0.601036  HALLMARK_MYC_TARGETS_V1
Ld Init Enhancement Gt100   Hallmark Myc Targets V1                     0.532298  0.00309969  1.71254  0.0467353   0.154985          3045   0.569948  HALLMARK_MYC_TARGETS_V1
Vol Late Lt0                Hallmark Myc Targets V1                     0.563096  0.00089991  1.80526  0.0339537   0.0449955         2812   0.601036  HALLMARK_MYC_TARGETS_V1
Ld Late Lt0                 Hallmark Myc Targets V1                     0.558002  0.00069993  1.79927  0.0343509   0.0349965         2767   0.57513   HALLMARK_MYC_TARGETS_V1
Irregularity                Hallmark Myc Targets V2                     0.471714  0.0885975   1.3883   0.180171    1                 3615   0.581818  HALLMARK_MYC_TARGETS_V2
Volume                      Hallmark Myc Targets V2                     0.672466  0.00050045  1.97029  0.0154117   0.0250225         2960   0.818182  HALLMARK_MYC_TARGETS_V2
Largest Diameter            Hallmark Myc Targets V2                     0.616533  0.00480577  1.81347  0.0330003   0.240288          2487   0.654545  HALLMARK_MYC_TARGETS_V2
Vol Init Enhancement Gt100  Hallmark Myc Targets V2                     0.628815  0.00310279  1.84455  0.0370623   0.15514           2741   0.727273  HALLMARK_MYC_TARGETS_V2
Ld Init Enhancement Gt100   Hallmark Myc Targets V2                     0.609964  0.00580697  1.79427  0.0355519   0.290348          3107   0.709091  HALLMARK_MYC_TARGETS_V2
Vol Late Lt0                Hallmark Myc Targets V2                     0.65343   0.00170153  1.91731  0.0219185   0.0850766         2093   0.672727  HALLMARK_MYC_TARGETS_V2
Ld Late Lt0                 Hallmark Myc Targets V2                     0.637843  0.00280477  1.87528  0.0289218   0.140238          2533   0.690909  HALLMARK_MYC_TARGETS_V2
Irregularity                Hallmark Epithelial Mesenchymal Transition  0.510312  0.0561944   1.62164  0.100693    1                 3495   0.573034  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
Volume                      Hallmark Epithelial Mesenchymal Transition  0.532664  0.039996    1.70594  0.050318    1                 2802   0.55618   HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
Largest Diameter            Hallmark Epithelial Mesenchymal Transition  0.600656  0.0115988   1.92012  0.0218501   0.579942          2256   0.589888  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
Variance Vox Val            Hallmark Epithelial Mesenchymal Transition  0.543811  0.0388      1.71844  0.0637984   1                 3019   0.589888  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
Vol Init Enhancement Gt100  Hallmark Epithelial Mesenchymal Transition  0.530256  0.0396      1.69675  0.056264    1                 3154   0.58427   HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
Ld Init Enhancement Gt100   Hallmark Epithelial Mesenchymal Transition  0.592565  0.0126987   1.89556  0.0254192   0.634937          3104   0.674157  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
Vol Late Lt0                Hallmark Epithelial Mesenchymal Transition  0.484227  0.0814      1.54722  0.126887    1                 3670   0.595506  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
Ld Late Lt0                 Hallmark Epithelial Mesenchymal Transition  0.512036  0.0538946   1.63627  0.0723946   1                 2671   0.52809   HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
Irregularity                Hallmark Glycolysis                         0.421999  0.00789921  1.29951  0.238173    0.394961          3228   0.427711  HALLMARK_GLYCOLYSIS
Largest Diameter            Hallmark Glycolysis                         0.439438  0.00239976  1.35163  0.163548    0.119988          3643   0.457831  HALLMARK_GLYCOLYSIS
Ld Init Enhancement Gt100   Hallmark Glycolysis                         0.423937  0.00689931  1.30292  0.177469    0.344966          3798   0.457831  HALLMARK_GLYCOLYSIS
Ld Late Lt0                 Hallmark Glycolysis                         0.417464  0.009999    1.28371  0.221385    0.49995           3959   0.46988   HALLMARK_GLYCOLYSIS
Largest Diameter            Hallmark Reactive Oxigen Species Pathway    0.433971  0.0875525   1.25086  0.19224     1                 4505   0.547619  HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY
Ld Init Enhancement Gt100   Hallmark Reactive Oxigen Species Pathway    0.431177  0.0942471   1.24448  0.2012      1                 4924   0.595238  HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY
Largest Diameter            Hallmark P53 Pathway                        0.39948   0.0236976   1.23654  0.19224     1                 3804   0.446927  HALLMARK_P53_PATHWAY
Ld Init Enhancement Gt100   Hallmark P53 Pathway                        0.423824  0.00619938  1.31225  0.176333    0.309969          4000   0.486034  HALLMARK_P53_PATHWAY
Ld Late Lt0                 Hallmark P53 Pathway                        0.395802  0.0284972   1.22533  0.231592    1                 3000   0.368715  HALLMARK_P53_PATHWAY
Irregularity                Hallmark Uv Response Up                     0.441632  0.00329967  1.34867  0.190522    0.164984          2623   0.380952  HALLMARK_UV_RESPONSE_UP
Largest Diameter            Hallmark Uv Response Up                     0.428395  0.010099    1.30593  0.17295     0.50495           4824   0.579365  HALLMARK_UV_RESPONSE_UP
Ld Init Enhancement Gt100   Hallmark Uv Response Up                     0.415952  0.0214979   1.26868  0.193519    1                 3141   0.404762  HALLMARK_UV_RESPONSE_UP
Ld Late Lt0                 Hallmark Uv Response Up                     0.418325  0.0170983   1.27762  0.221385    0.854915          4350   0.52381   HALLMARK_UV_RESPONSE_UP
Largest Diameter            Hallmark Uv Response Dn                     0.385683  0.0852915   1.20008  0.230266    1                 4822   0.530303  HALLMARK_UV_RESPONSE_DN
Ld Init Enhancement Gt100   Hallmark Uv Response Dn                     0.391873  0.0681932   1.21921  0.22304     1                 5427   0.613636  HALLMARK_UV_RESPONSE_DN
Irregularity                Hallmark Angiogenesis                       0.50339   0.0879742   1.37691  0.180171    1                 2725   0.517241  HALLMARK_ANGIOGENESIS
Volume                      Hallmark Angiogenesis                       0.497823  0.0934947   1.36907  0.215029    1                 3694   0.62069   HALLMARK_ANGIOGENESIS
Largest Diameter            Hallmark Angiogenesis                       0.550272  0.0324308   1.51094  0.0986456   1                 2029   0.482759  HALLMARK_ANGIOGENESIS
Mean Vox Val                Hallmark Angiogenesis                       0.527581  0.0582299   1.44208  0.247161    1                 4450   0.724138  HALLMARK_ANGIOGENESIS
Variance Vox Val            Hallmark Angiogenesis                       0.549025  0.037732    1.49637  0.11939     1                 3712   0.655172  HALLMARK_ANGIOGENESIS
Vol Init Enhancement Gt100  Hallmark Angiogenesis                       0.52425   0.0551348   1.44256  0.149967    1                 3896   0.689655  HALLMARK_ANGIOGENESIS
Ld Init Enhancement Gt100   Hallmark Angiogenesis                       0.53125   0.0514245   1.4604   0.115854    1                 3119   0.586207  HALLMARK_ANGIOGENESIS
Ld Late Lt0                 Hallmark Angiogenesis                       0.508244  0.0788252   1.39559  0.175199    1                 2920   0.551724  HALLMARK_ANGIOGENESIS
Largest Diameter            Hallmark Coagulation                        0.496778  0.00989901  1.47949  0.0986456   0.494951          3178   0.494505  HALLMARK_COAGULATION
Variance Vox Val            Hallmark Coagulation                        0.450299  0.0478952   1.33315  0.210306    1                 3467   0.483516  HALLMARK_COAGULATION
Ld Init Enhancement Gt100   Hallmark Coagulation                        0.473078  0.0213979   1.40714  0.139897    1                 3525   0.505495  HALLMARK_COAGULATION
Ld Late Lt0                 Hallmark Coagulation                        0.424813  0.080192    1.26463  0.221385    1                 3418   0.461538  HALLMARK_COAGULATION
Largest Diameter            Hallmark Il2 Stat5 Signaling                0.400957  0.0460954   1.23883  0.19224     1                 4159   0.486667  HALLMARK_IL2_STAT5_SIGNALING
Ld Init Enhancement Gt100   Hallmark Il2 Stat5 Signaling                0.389525  0.0692931   1.20431  0.22304     1                 4100   0.473333  HALLMARK_IL2_STAT5_SIGNALING
Largest Diameter            Hallmark Peroxisome                         0.42527   0.0421958   1.23684  0.19224     1                 3707   0.406977  HALLMARK_PEROXISOME
Ld Init Enhancement Gt100   Hallmark Peroxisome                         0.435172  0.0259974   1.26548  0.193519    1                 3701   0.418605  HALLMARK_PEROXISOME
Ld Late Lt0                 Hallmark Peroxisome                         0.426401  0.0385961   1.23984  0.221385    1                 2313   0.313953  HALLMARK_PEROXISOME
Irregularity                Hallmark Spermatogenesis                    0.571019  0.00359964  1.6107   0.100693    0.179982          1332   0.385965  HALLMARK_SPERMATOGENESIS
Volume                      Hallmark Spermatogenesis                    0.524219  0.0182982   1.47324  0.139505    0.914909           604   0.298246  HALLMARK_SPERMATOGENESIS
Largest Diameter            Hallmark Spermatogenesis                    0.549927  0.0066      1.54664  0.0986456   0.33              2848   0.508772  HALLMARK_SPERMATOGENESIS
Mean Vox Val                Hallmark Spermatogenesis                    0.510354  0.030397    1.43376  0.247161    1                 2464   0.438596  HALLMARK_SPERMATOGENESIS
Variance Vox Val            Hallmark Spermatogenesis                    0.608164  0.00029997  1.70805  0.0637984   0.0149985         1706   0.438596  HALLMARK_SPERMATOGENESIS
Top Late Enhancement        Hallmark Spermatogenesis                    0.54021   0.0105      1.51705  0.243525    0.525             2275   0.438596  HALLMARK_SPERMATOGENESIS
Vol Init Enhancement Gt100  Hallmark Spermatogenesis                    0.522971  0.0163      1.46912  0.144184    0.815              630   0.298246  HALLMARK_SPERMATOGENESIS
Ld Init Enhancement Gt100   Hallmark Spermatogenesis                    0.55011   0.00529947  1.54548  0.102477    0.264974          2571   0.473684  HALLMARK_SPERMATOGENESIS
Ld Late Lt0                 Hallmark Spermatogenesis                    0.529765  0.0141986   1.48903  0.109357    0.709929          1790   0.385965  HALLMARK_SPERMATOGENESIS
Largest Diameter            Hallmark Kras Signaling Up                  0.407948  0.0906909   1.24532  0.19224     1                 3780   0.454545  HALLMARK_KRAS_SIGNALING_UP
Variance Vox Val            Hallmark Kras Signaling Up                  0.444082  0.0320968   1.35069  0.210306    1                 4025   0.552448  HALLMARK_KRAS_SIGNALING_UP
Ld Init Enhancement Gt100   Hallmark Kras Signaling Up                  0.404424  0.0988901   1.23521  0.206815    1                 3643   0.440559  HALLMARK_KRAS_SIGNALING_UP
Largest Diameter            Hallmark Pancreas Beta Cells                0.500599  0.184024    1.21764  0.209768    1                 4434   0.647059  HALLMARK_PANCREAS_BETA_CELLS
Ld Init Enhancement Gt100   Hallmark Pancreas Beta Cells                0.489316  0.215651    1.18849  0.242055    1                 4982   0.705882  HALLMARK_PANCREAS_BETA_CELLS
Ld Late Lt0                 Hallmark Pancreas Beta Cells                0.509492  0.15896     1.23723  0.221385    1                 3968   0.588235  HALLMARK_PANCREAS_BETA_CELLS

</div>
