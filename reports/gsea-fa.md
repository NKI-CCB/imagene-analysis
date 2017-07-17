---
title: Gene Set Enrichment Analysis of CAD factors
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
df_h = load_gsea_ds("../analyses/gsea/mri-features-fa_h.all_T.nc")
df_cgp = load_gsea_ds("../analyses/gsea/mri-features-fa_c2.cgp_F.nc")
df_cp = load_gsea_ds("../analyses/gsea/mri-features-fa_c2.cp_T.nc")
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

### Factor Analysis ###


```python
fa_dataset = xr.open_dataset("../data/processed/mri-features-fa.nc").load()
plot.heatmap(fa_dataset['loadings'].T)
fa_dataset.close()
```

![](figures/gsea-fa_loadings-heatmap_1){#loadings-heatmap }\


### GSEA ###

Hallmarks from MSigDB.


```python
plot_ds(df_h, fdr=0.25)
```

![](figures/gsea-fa_hallmarks-es-plot_1){#hallmarks-es-plot }\



```python
table_ds(df_h, fdr=0.25)
```


<div class="datatable">  mri_feature  gene_set                                          es           p      nes        fdr       fwer    max_es_at    le_prop  gene_set_code
-------------  ------------------------------------------  --------  ----------  -------  ---------  ---------  -----------  ---------  ------------------------------------------
            1  Hallmark Tnfa Signaling Via Nfkb            0.404068  0.09989     1.26032  0.243126   1                 2943   0.375     HALLMARK_TNFA_SIGNALING_VIA_NFKB
            1  Hallmark Hypoxia                            0.442754  0.00839916  1.36392  0.18121    0.419958          4087   0.539394  HALLMARK_HYPOXIA
           10  Hallmark Hypoxia                            0.446571  0.00669933  1.37539  0.221197   0.334967          3062   0.460606  HALLMARK_HYPOXIA
            1  Hallmark Cholesterol Homeostasis            0.433461  0.0556      1.2731   0.243126   1                 4896   0.615385  HALLMARK_CHOLESTEROL_HOMEOSTASIS
            1  Hallmark Mitotic Spindle                    0.474633  0.00359964  1.46669  0.120276   0.179982          1498   0.260417  HALLMARK_MITOTIC_SPINDLE
            5  Hallmark Mitotic Spindle                    0.468536  0.00509949  1.45224  0.187469   0.254975          2483   0.364583  HALLMARK_MITOTIC_SPINDLE
           10  Hallmark Mitotic Spindle                    0.457638  0.00919908  1.41975  0.200093   0.459954          2391   0.34375   HALLMARK_MITOTIC_SPINDLE
            1  Hallmark Tgf Beta Signaling                 0.415936  0.129265    1.21175  0.245631   1                 1353   0.291667  HALLMARK_TGF_BETA_SIGNALING
           10  Hallmark Tgf Beta Signaling                 0.449681  0.0513      1.30943  0.244034   1                 3514   0.458333  HALLMARK_TGF_BETA_SIGNALING
            1  Hallmark Il6 Jak Stat3 Signaling            0.41441   0.153746    1.22026  0.245631   1                 2939   0.377049  HALLMARK_IL6_JAK_STAT3_SIGNALING
            1  Hallmark G2M Checkpoint                     0.668156  0.0009999   2.04023  0.0150636  0.049995          1661   0.561111  HALLMARK_G2M_CHECKPOINT
            5  Hallmark G2M Checkpoint                     0.666704  0.00039996  2.04919  0.0153624  0.019998          1236   0.511111  HALLMARK_G2M_CHECKPOINT
           10  Hallmark G2M Checkpoint                     0.652118  0.00149985  1.9985   0.0205183  0.0749925         2352   0.633333  HALLMARK_G2M_CHECKPOINT
            1  Hallmark Apoptosis                          0.434595  0.00749925  1.34578  0.188287   0.374963          4597   0.613139  HALLMARK_APOPTOSIS
            1  Hallmark Adipogenesis                       0.385017  0.0526947   1.2012   0.24749    1                 3663   0.437838  HALLMARK_ADIPOGENESIS
            1  Hallmark Myogenesis                         0.414695  0.0724928   1.26357  0.243126   1                 4748   0.566929  HALLMARK_MYOGENESIS
            1  Hallmark Unfolded Protein Response          0.39315   0.0811919   1.20818  0.245631   1                 3132   0.376147  HALLMARK_UNFOLDED_PROTEIN_RESPONSE
            1  Hallmark Mtorc1 Signaling                   0.466252  0.00289971  1.45713  0.120276   0.144986          3495   0.492063  HALLMARK_MTORC1_SIGNALING
           10  Hallmark Mtorc1 Signaling                   0.436109  0.0128987   1.36429  0.221197   0.644936          2563   0.380952  HALLMARK_MTORC1_SIGNALING
            1  Hallmark E2F Targets                        0.679257  0.00169983  2.08753  0.0150636  0.0849915         1767   0.616667  HALLMARK_E2F_TARGETS
            5  Hallmark E2F Targets                        0.67055   0.0017      2.0769   0.0153624  0.085             2384   0.683333  HALLMARK_E2F_TARGETS
           10  Hallmark E2F Targets                        0.665255  0.0021      2.05047  0.0205183  0.105             2002   0.638889  HALLMARK_E2F_TARGETS
            1  Hallmark Myc Targets V1                     0.546172  0.00269973  1.74489  0.0414776  0.134987          3199   0.611399  HALLMARK_MYC_TARGETS_V1
            5  Hallmark Myc Targets V1                     0.465722  0.0266973   1.49395  0.174041   1                 3162   0.466321  HALLMARK_MYC_TARGETS_V1
           10  Hallmark Myc Targets V1                     0.47658   0.019898    1.52473  0.146351   0.994901          3244   0.507772  HALLMARK_MYC_TARGETS_V1
            1  Hallmark Myc Targets V2                     0.611558  0.00460737  1.78412  0.0414776  0.230369          2746   0.672727  HALLMARK_MYC_TARGETS_V2
           10  Hallmark Myc Targets V2                     0.485087  0.0783253   1.41827  0.200093   1                 3809   0.636364  HALLMARK_MYC_TARGETS_V2
            1  Hallmark Epithelial Mesenchymal Transition  0.558487  0.0281      1.77131  0.0414776  1                 1912   0.505618  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
            5  Hallmark Epithelial Mesenchymal Transition  0.495334  0.0751925   1.56911  0.162131   1                 4256   0.657303  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
           10  Hallmark Epithelial Mesenchymal Transition  0.555664  0.0286      1.76413  0.061655   1                 3644   0.651685  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
           10  Hallmark Glycolysis                         0.43409   0.00319968  1.33593  0.221197   0.159984          3128   0.427711  HALLMARK_GLYCOLYSIS
            1  Hallmark P53 Pathway                        0.38869   0.0419958   1.20745  0.245631   1                 3050   0.363128  HALLMARK_P53_PATHWAY
            1  Hallmark Uv Response Up                     0.419914  0.0135986   1.28253  0.243126   0.679932          4130   0.5       HALLMARK_UV_RESPONSE_UP
            1  Hallmark Uv Response Dn                     0.390866  0.0776922   1.21106  0.245631   1                 5047   0.583333  HALLMARK_UV_RESPONSE_DN
            1  Hallmark Angiogenesis                       0.536293  0.045716    1.46426  0.120276   1                 4216   0.724138  HALLMARK_ANGIOGENESIS
           10  Hallmark Angiogenesis                       0.492858  0.0977943   1.34905  0.221197   1                 2905   0.517241  HALLMARK_ANGIOGENESIS
            1  Hallmark Coagulation                        0.472178  0.0243976   1.40242  0.153459   1                 3171   0.472527  HALLMARK_COAGULATION
            1  Hallmark Il2 Stat5 Signaling                0.39038   0.0665933   1.20759  0.245631   1                 3815   0.44      HALLMARK_IL2_STAT5_SIGNALING
            1  Hallmark Peroxisome                         0.413277  0.0607939   1.21535  0.245631   1                 2195   0.302326  HALLMARK_PEROXISOME
            1  Hallmark Spermatogenesis                    0.530093  0.0173      1.48874  0.120276   0.865             2051   0.403509  HALLMARK_SPERMATOGENESIS
            5  Hallmark Spermatogenesis                    0.546129  0.00989901  1.54434  0.162131   0.494951          1291   0.350877  HALLMARK_SPERMATOGENESIS
           10  Hallmark Spermatogenesis                    0.556826  0.0050005   1.56456  0.144704   0.250025          2920   0.526316  HALLMARK_SPERMATOGENESIS
            1  Hallmark Kras Signaling Up                  0.407231  0.0832917   1.25199  0.244149   1                 3927   0.461538  HALLMARK_KRAS_SIGNALING_UP
            1  Hallmark Pancreas Beta Cells                0.536574  0.0856039   1.3333   0.189703   1                 4661   0.705882  HALLMARK_PANCREAS_BETA_CELLS

</div>

Gene signatures (c2.cgp) from MSigDB.


```python
plot_ds(df_cgp, fdr=0.05, le_prop=0.0, abs=False)
```

![](figures/gsea-fa_cgp-es-plot-a_1){#cgp-es-plot-a }\



```python
plot_ds(df_cgp, fdr=0.05, le_prop=0.85, abs=False)
```

![](figures/gsea-fa_cgp-es-plot-b_1){#cgp-es-plot-b }\



```python
table_ds(df_cgp, fdr=0.25)
```


<div class="datatable">  mri_feature  gene_set                                                           es            p       nes        fdr      fwer    max_es_at    le_prop  gene_set_code
-------------  ----------------------------------------------------------  ---------  -----------  --------  ---------  --------  -----------  ---------  ----------------------------------------------------------
            1  Nakamura Cancer Microenvironment Up                         -0.750374  0.00646073   -1.81728  0.119811   1                9725   0.533333  NAKAMURA_CANCER_MICROENVIRONMENT_UP
           10  Nakamura Cancer Microenvironment Up                         -0.853511  0.000203625  -2.06878  0.203325   0.541845        10310   0.533333  NAKAMURA_CANCER_MICROENVIRONMENT_UP
            1  Nakamura Cancer Microenvironment Dn                          0.744813  0.00162173    1.84912  0.0308973  1                1455   0.682927  NAKAMURA_CANCER_MICROENVIRONMENT_DN
            4  Nakamura Cancer Microenvironment Dn                         -0.694099  0.0159574    -1.7234   0.192508   1                9501   0.585366  NAKAMURA_CANCER_MICROENVIRONMENT_DN
            5  Nakamura Cancer Microenvironment Dn                          0.641987  0.0611798     1.58466  0.225622   1                2379   0.560976  NAKAMURA_CANCER_MICROENVIRONMENT_DN
           10  Nakamura Cancer Microenvironment Dn                          0.663669  0.0361156     1.65992  0.164561   1                1971   0.658537  NAKAMURA_CANCER_MICROENVIRONMENT_DN
            1  West Adrenocortical Tumor Markers Up                         0.769349  0.0121581     1.74657  0.0588666  1                1165   0.647059  WEST_ADRENOCORTICAL_TUMOR_MARKERS_UP
            4  West Adrenocortical Tumor Markers Up                        -0.755466  0.0243227    -1.7098   0.192508   1                9270   0.588235  WEST_ADRENOCORTICAL_TUMOR_MARKERS_UP
            5  West Adrenocortical Tumor Markers Up                         0.754825  0.0244489     1.69785  0.14161    1                 810   0.529412  WEST_ADRENOCORTICAL_TUMOR_MARKERS_UP
           10  West Adrenocortical Tumor Markers Up                         0.695911  0.0711844     1.58599  0.224283   1                1423   0.647059  WEST_ADRENOCORTICAL_TUMOR_MARKERS_UP
            1  West Adrenocortical Tumor Markers Dn                        -0.788914  0.0136763    -1.69942  0.148746   1                9476   0.8       WEST_ADRENOCORTICAL_TUMOR_MARKERS_DN
            1  Winter Hypoxia Up                                            0.443227  0.109485      1.43726  0.245033   1                1621   0.341463  WINTER_HYPOXIA_UP
            4  Winter Hypoxia Up                                           -0.526616  0.0223745    -1.70209  0.192508   1                8611   0.54878   WINTER_HYPOXIA_UP
            1  Parent Mtor Signaling Dn                                    -0.44189   0.0466332    -1.49529  0.232098   1                9198   0.4       PARENT_MTOR_SIGNALING_DN
            1  Pyeon Hpv Positive Tumors Up                                 0.682153  0.00121139    1.97609  0.013915   1                1182   0.507246  PYEON_HPV_POSITIVE_TUMORS_UP
            5  Pyeon Hpv Positive Tumors Up                                 0.766034  0.000196889   2.18391  0.0172341  0.523922         1302   0.637681  PYEON_HPV_POSITIVE_TUMORS_UP
           10  Pyeon Hpv Positive Tumors Up                                 0.70127   0.00121852    2.01936  0.056431   1                1651   0.608696  PYEON_HPV_POSITIVE_TUMORS_UP
            4  Nakamura Tumor Zone Peripheral Vs Central Up                -0.501092  0.00199005   -1.84619  0.174932   1                8937   0.475983  NAKAMURA_TUMOR_ZONE_PERIPHERAL_VS_CENTRAL_UP
            1  Piccaluga Angioimmunoblastic Lymphoma Up                    -0.578937  0.0554435    -1.66807  0.162484   1                9349   0.507853  PICCALUGA_ANGIOIMMUNOBLASTIC_LYMPHOMA_UP
            1  Piccaluga Angioimmunoblastic Lymphoma Dn                    -0.43011   0.0674493    -1.51898  0.221755   1                7467   0.539683  PICCALUGA_ANGIOIMMUNOBLASTIC_LYMPHOMA_DN
            8  Korkola Embryonal Carcinoma Up                              -0.729166  0.000998203  -2.0165   0.21769    1                9737   0.636364  KORKOLA_EMBRYONAL_CARCINOMA_UP
            8  Korkola Yolk Sac Tumor Up                                   -0.795962  0.000597372  -1.99274  0.21769    1                8962   0.842105  KORKOLA_YOLK_SAC_TUMOR_UP
            1  Liu Sox4 Targets Dn                                          0.385056  0.00481928    1.66422  0.0969555  1                1854   0.326241  LIU_SOX4_TARGETS_DN
           10  Liu Sox4 Targets Dn                                          0.376752  0.00946055    1.6254   0.188431   1                3502   0.48227   LIU_SOX4_TARGETS_DN
            1  Bertucci Medullary Vs Ductal Breast Cancer Dn               -0.576533  0.054328     -1.66559  0.162484   1                8154   0.605634  BERTUCCI_MEDULLARY_VS_DUCTAL_BREAST_CANCER_DN
            1  Davicioni Pax Foxo1 Signature In Arms Dn                    -0.643505  0.0815617    -1.50068  0.230381   1                7873   0.785714  DAVICIONI_PAX_FOXO1_SIGNATURE_IN_ARMS_DN
            1  Fournier Acinar Development Late Dn                          0.802211  0.00262732    1.82132  0.0363179  1                 165   0.555556  FOURNIER_ACINAR_DEVELOPMENT_LATE_DN
            5  Fournier Acinar Development Late Dn                          0.762347  0.0148163     1.72484  0.126645   1                 942   0.611111  FOURNIER_ACINAR_DEVELOPMENT_LATE_DN
           10  Fournier Acinar Development Late Dn                          0.785896  0.00564061    1.78114  0.0904405  1                 714   0.611111  FOURNIER_ACINAR_DEVELOPMENT_LATE_DN
            1  Sengupta Nasopharyngeal Carcinoma Up                         0.48243   0.04885       1.61811  0.120857   1                 999   0.314516  SENGUPTA_NASOPHARYNGEAL_CARCINOMA_UP
            5  Sengupta Nasopharyngeal Carcinoma Up                         0.570672  0.00358066    1.90124  0.053926   1                1784   0.512097  SENGUPTA_NASOPHARYNGEAL_CARCINOMA_UP
            1  Kobayashi Egfr Signaling 6Hr Dn                             -0.66145   0.0485789    -1.57678  0.194122   1                8948   0.615385  KOBAYASHI_EGFR_SIGNALING_6HR_DN
            1  Sotiriou Breast Cancer Grade 1 Vs 3 Up                       0.879381  0.000201939   1.83869  0.0319879  0.537359          636   0.838983  SOTIRIOU_BREAST_CANCER_GRADE_1_VS_3_UP
            4  Sotiriou Breast Cancer Grade 1 Vs 3 Up                      -0.850474  0.00137147   -1.79385  0.192508   1                9865   0.872881  SOTIRIOU_BREAST_CANCER_GRADE_1_VS_3_UP
            5  Sotiriou Breast Cancer Grade 1 Vs 3 Up                       0.80978   0.0165405     1.69511  0.14296    1                 597   0.652542  SOTIRIOU_BREAST_CANCER_GRADE_1_VS_3_UP
           10  Sotiriou Breast Cancer Grade 1 Vs 3 Up                       0.832155  0.00439122    1.75002  0.108321   1                1169   0.830508  SOTIRIOU_BREAST_CANCER_GRADE_1_VS_3_UP
            1  Sotiriou Breast Cancer Grade 1 Vs 3 Dn                      -0.711536  0.00627023   -1.75074  0.139914   1                8708   0.818182  SOTIRIOU_BREAST_CANCER_GRADE_1_VS_3_DN
            1  Gazda Diamond Blackfan Anemia Myeloid Up                    -0.554389  0.0164758    -1.70447  0.147212   1                8713   0.5       GAZDA_DIAMOND_BLACKFAN_ANEMIA_MYELOID_UP
            1  Chemnitz Response To Prostaglandin E2 Up                     0.685735  0.00100685    1.98418  0.0134874  1                 812   0.475     CHEMNITZ_RESPONSE_TO_PROSTAGLANDIN_E2_UP
            4  Chemnitz Response To Prostaglandin E2 Up                    -0.610734  0.0238374    -1.77829  0.192508   1                9364   0.516667  CHEMNITZ_RESPONSE_TO_PROSTAGLANDIN_E2_UP
            5  Chemnitz Response To Prostaglandin E2 Up                     0.597055  0.035461      1.72372  0.126959   1                1940   0.508333  CHEMNITZ_RESPONSE_TO_PROSTAGLANDIN_E2_UP
           10  Chemnitz Response To Prostaglandin E2 Up                     0.638848  0.00940564    1.85515  0.0652612  1                1590   0.491667  CHEMNITZ_RESPONSE_TO_PROSTAGLANDIN_E2_UP
            1  Zhong Response To Azacitidine And Tsa Up                    -0.387384  0.0293238    -1.54017  0.209703   1                9311   0.299213  ZHONG_RESPONSE_TO_AZACITIDINE_AND_TSA_UP
            5  Zhong Response To Azacitidine And Tsa Dn                     0.451732  0.0253215     1.60251  0.209048   1                1553   0.306452  ZHONG_RESPONSE_TO_AZACITIDINE_AND_TSA_DN
            1  Davicioni Molecular Arms Vs Erms Dn                         -0.503683  0.0173966    -1.76357  0.139375   1                8920   0.439189  DAVICIONI_MOLECULAR_ARMS_VS_ERMS_DN
            1  Davicioni Targets Of Pax Foxo1 Fusions Up                   -0.463941  0.041432     -1.63327  0.170938   1                8458   0.490741  DAVICIONI_TARGETS_OF_PAX_FOXO1_FUSIONS_UP
            8  Turashvili Breast Ductal Carcinoma Vs Ductal Normal Up      -0.709013  0.00221953   -2.01877  0.21769    1               10563   0.439024  TURASHVILI_BREAST_DUCTAL_CARCINOMA_VS_DUCTAL_NORMAL_UP
            1  Turashvili Breast Ductal Carcinoma Vs Lobular Normal Dn     -0.534861  0.104968     -1.48485  0.241156   1                9066   0.428571  TURASHVILI_BREAST_DUCTAL_CARCINOMA_VS_LOBULAR_NORMAL_DN
            1  Turashvili Breast Lobular Carcinoma Vs Ductal Normal Up     -0.614809  0.125702     -1.53862  0.210383   1                9718   0.522388  TURASHVILI_BREAST_LOBULAR_CARCINOMA_VS_DUCTAL_NORMAL_UP
            1  Turashvili Breast Lobular Carcinoma Vs Lobular Normal Up    -0.451804  0.0818       -1.46961  0.248243   1                8836   0.395349  TURASHVILI_BREAST_LOBULAR_CARCINOMA_VS_LOBULAR_NORMAL_UP
            1  Turashvili Breast Lobular Carcinoma Vs Lobular Normal Dn    -0.567977  0.125177     -1.53275  0.214421   1                9752   0.450704  TURASHVILI_BREAST_LOBULAR_CARCINOMA_VS_LOBULAR_NORMAL_DN
            1  Chandran Metastasis Top50 Dn                                -0.462042  0.0251879    -1.60359  0.181764   1                7790   0.522727  CHANDRAN_METASTASIS_TOP50_DN
            1  Wilcox Response To Progesterone Up                           0.527384  0.00790274    1.83982  0.0319879  1                1207   0.381356  WILCOX_RESPONSE_TO_PROGESTERONE_UP
            5  Wilcox Response To Progesterone Up                           0.580054  0.000398248   2.01991  0.0336414  1                 945   0.364407  WILCOX_RESPONSE_TO_PROGESTERONE_UP
           10  Wilcox Response To Progesterone Up                           0.498565  0.0201532     1.73059  0.115709   1                1160   0.355932  WILCOX_RESPONSE_TO_PROGESTERONE_UP
            1  Zhou Inflammatory Response Live Up                          -0.369393  0.058006     -1.46816  0.249633   1                8711   0.33452   ZHOU_INFLAMMATORY_RESPONSE_LIVE_UP
            1  Samols Targets Of Khsv Mirnas Dn                            -0.458041  0.0194076    -1.63416  0.170938   1                8681   0.461538  SAMOLS_TARGETS_OF_KHSV_MIRNAS_DN
            1  Hooi St7 Targets Dn                                         -0.414933  0.0291223    -1.54921  0.202734   1                8467   0.413793  HOOI_ST7_TARGETS_DN
            1  Pramoonjago Sox4 Targets Dn                                  0.537904  0.0266296     1.65203  0.103487   1                2609   0.571429  PRAMOONJAGO_SOX4_TARGETS_DN
            4  Pramoonjago Sox4 Targets Dn                                 -0.547054  0.0168986    -1.67884  0.204423   1                9191   0.47619   PRAMOONJAGO_SOX4_TARGETS_DN
            1  Puiffe Invasion Inhibited By Ascites Up                      0.432473  0.0347041     1.5669   0.158229   1                1847   0.405797  PUIFFE_INVASION_INHIBITED_BY_ASCITES_UP
            1  Casorelli Acute Promyelocytic Leukemia Up                   -0.424676  0.0193939    -1.62951  0.170938   1                7372   0.555556  CASORELLI_ACUTE_PROMYELOCYTIC_LEUKEMIA_UP
            1  Huttmann B Cll Poor Survival Up                             -0.358406  0.0363563    -1.50248  0.23021    1                8399   0.367925  HUTTMANN_B_CLL_POOR_SURVIVAL_UP
            1  Casorelli Apl Secondary Vs De Novo Up                       -0.515482  0.0321537    -1.60558  0.180997   1                8512   0.5625    CASORELLI_APL_SECONDARY_VS_DE_NOVO_UP
            1  Liu Targets Of Vmyb Vs Cmyb Dn                              -0.541536  0.0177598    -1.67249  0.161802   1                7732   0.586207  LIU_TARGETS_OF_VMYB_VS_CMYB_DN
            1  Charafe Breast Cancer Basal Vs Mesenchymal Dn               -0.668777  0.0156031    -1.7524   0.139914   1                8790   0.702703  CHARAFE_BREAST_CANCER_BASAL_VS_MESENCHYMAL_DN
            1  Doane Breast Cancer Esr1 Up                                 -0.654283  0.0429298    -1.6836   0.156945   1                8190   0.733945  DOANE_BREAST_CANCER_ESR1_UP
            1  Borczuk Malignant Mesothelioma Up                            0.51997   0.0092873     1.81186  0.0385687  1                2567   0.515254  BORCZUK_MALIGNANT_MESOTHELIOMA_UP
            4  Borczuk Malignant Mesothelioma Up                           -0.490557  0.0235623    -1.70586  0.192508   1                8508   0.505085  BORCZUK_MALIGNANT_MESOTHELIOMA_UP
            5  Borczuk Malignant Mesothelioma Up                            0.466644  0.0469718     1.62701  0.19014    1                2416   0.416949  BORCZUK_MALIGNANT_MESOTHELIOMA_UP
           10  Borczuk Malignant Mesothelioma Up                            0.453275  0.0594258     1.57306  0.23642    1                2982   0.467797  BORCZUK_MALIGNANT_MESOTHELIOMA_UP
            1  Roy Wound Blood Vessel Up                                   -0.596268  0.0694804    -1.60019  0.181764   1                9718   0.511111  ROY_WOUND_BLOOD_VESSEL_UP
            1  Newman Ercc6 Targets Dn                                     -0.661347  0.0560691    -1.60619  0.180997   1                7911   0.76      NEWMAN_ERCC6_TARGETS_DN
            1  Horiuchi Wtap Targets Up                                    -0.455294  0.00581512   -1.75734  0.139914   1                8094   0.461847  HORIUCHI_WTAP_TARGETS_UP
            1  Horiuchi Wtap Targets Dn                                     0.607432  0.00161649    1.98566  0.0134874  1                1742   0.461832  HORIUCHI_WTAP_TARGETS_DN
            4  Horiuchi Wtap Targets Dn                                    -0.511193  0.0499603    -1.66334  0.214503   1                8895   0.473282  HORIUCHI_WTAP_TARGETS_DN
            5  Horiuchi Wtap Targets Dn                                     0.580652  0.00533808    1.89931  0.0539993  1                2350   0.526718  HORIUCHI_WTAP_TARGETS_DN
           10  Horiuchi Wtap Targets Dn                                     0.578088  0.00619628    1.88758  0.0627622  1                1555   0.416031  HORIUCHI_WTAP_TARGETS_DN
            1  Gal Leukemic Stem Cell Dn                                    0.45822   0.0428803     1.61504  0.123113   1                 547   0.233766  GAL_LEUKEMIC_STEM_CELL_DN
            1  Basaki Ybx1 Targets Up                                       0.637148  0.000805964   2.00779  0.0127991  1                1181   0.515284  BASAKI_YBX1_TARGETS_UP
            4  Basaki Ybx1 Targets Up                                      -0.54061   0.0428232    -1.69687  0.19736    1                9684   0.436681  BASAKI_YBX1_TARGETS_UP
            5  Basaki Ybx1 Targets Up                                       0.567451  0.0202824     1.7809   0.0972596  1                1603   0.471616  BASAKI_YBX1_TARGETS_UP
           10  Basaki Ybx1 Targets Up                                       0.597617  0.00606551    1.8807   0.0627622  1                1223   0.480349  BASAKI_YBX1_TARGETS_UP
            1  Nojima Sfrp2 Targets Up                                     -0.520429  0.0925775    -1.47048  0.24821    1                8283   0.565217  NOJIMA_SFRP2_TARGETS_UP
            1  Rodrigues Dcc Targets Dn                                    -0.519687  0.000200321  -1.92544  0.116729   0.533053         8372   0.5       RODRIGUES_DCC_TARGETS_DN
            1  Vecchi Gastric Cancer Advanced Vs Early Up                  -0.600111  0.093234     -1.58929  0.187154   1                9077   0.530612  VECCHI_GASTRIC_CANCER_ADVANCED_VS_EARLY_UP
            1  Vecchi Gastric Cancer Early Dn                              -0.578137  0.00141044   -1.93729  0.116729   1                9012   0.528409  VECCHI_GASTRIC_CANCER_EARLY_DN
            1  Slebos Head And Neck Cancer With Hpv Up                      0.685272  0.00205592    1.90521  0.0204626  1                1579   0.563636  SLEBOS_HEAD_AND_NECK_CANCER_WITH_HPV_UP
            5  Slebos Head And Neck Cancer With Hpv Up                      0.674996  0.00437724    1.8715   0.0603934  1                2090   0.618182  SLEBOS_HEAD_AND_NECK_CANCER_WITH_HPV_UP
           10  Slebos Head And Neck Cancer With Hpv Up                      0.668403  0.00555441    1.8562   0.0652612  1                1781   0.581818  SLEBOS_HEAD_AND_NECK_CANCER_WITH_HPV_UP
            1  Jaeger Metastasis Up                                         0.625615  0.00752338    1.80306  0.0411858  1                 501   0.382353  JAEGER_METASTASIS_UP
            4  Jaeger Metastasis Up                                        -0.662344  0.00339592   -1.90849  0.155097   1               10308   0.411765  JAEGER_METASTASIS_UP
            1  Ginestier Breast Cancer 20Q13 Amplification Dn               0.450426  0.12819       1.44877  0.238115   1                2839   0.43871   GINESTIER_BREAST_CANCER_20Q13_AMPLIFICATION_DN
            5  Ginestier Breast Cancer 20Q13 Amplification Dn               0.497489  0.0650631     1.58345  0.225686   1                3050   0.541935  GINESTIER_BREAST_CANCER_20Q13_AMPLIFICATION_DN
           10  Ginestier Breast Cancer 20Q13 Amplification Dn               0.493897  0.0701611     1.57789  0.231959   1                3106   0.483871  GINESTIER_BREAST_CANCER_20Q13_AMPLIFICATION_DN
            1  Gargalovic Response To Oxidized Phospholipids Turquoise Up  -0.569023  0.00440882   -1.86903  0.119811   1                8735   0.513889  GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_TURQUOISE_UP
            1  Gargalovic Response To Oxidized Phospholipids Turquoise Dn   0.698401  0.00781563    1.81583  0.0377212  1                 847   0.461538  GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_TURQUOISE_DN
            5  Gargalovic Response To Oxidized Phospholipids Turquoise Dn   0.728084  0.00299341    1.87734  0.0582517  1                1149   0.564103  GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_TURQUOISE_DN
           10  Gargalovic Response To Oxidized Phospholipids Turquoise Dn   0.703826  0.00620124    1.82661  0.0688673  1                1423   0.512821  GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_TURQUOISE_DN
            1  Gargalovic Response To Oxidized Phospholipids Red Up         0.577819  0.0952572     1.45794  0.233099   1                1484   0.5       GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_RED_UP
            1  Gargalovic Response To Oxidized Phospholipids Red Dn        -0.638873  0.0300524    -1.61623  0.175815   1                8892   0.5625    GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_RED_DN
            1  Gargalovic Response To Oxidized Phospholipids Magenta Up    -0.568096  0.0530693    -1.57107  0.195726   1               10059   0.37037   GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_MAGENTA_UP
            1  Gargalovic Response To Oxidized Phospholipids Grey Up       -0.609811  0.0568204    -1.54623  0.203357   1                9558   0.533333  GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_GREY_UP
            1  Gargalovic Response To Oxidized Phospholipids Grey Dn       -0.445593  0.0367913    -1.54959  0.202734   1                9124   0.358491  GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_GREY_DN
            1  Takeda Targets Of Nup98 Hoxa9 Fusion 6Hr Up                 -0.472316  0.0368049    -1.55617  0.202734   1                7763   0.588235  TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_6HR_UP
            1  Takeda Targets Of Nup98 Hoxa9 Fusion 6Hr Dn                 -0.576435  0.010101     -1.77593  0.134462   1                9831   0.40625   TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_6HR_DN
            1  Lopez Mesotelioma Survival Time Up                           0.815591  0.0118545     1.68271  0.0870062  1                1062   0.75      LOPEZ_MESOTELIOMA_SURVIVAL_TIME_UP
            5  Lopez Mesotelioma Survival Time Up                           0.789094  0.0283773     1.61405  0.199637   1                 600   0.5       LOPEZ_MESOTELIOMA_SURVIVAL_TIME_UP
           10  Lopez Mesotelioma Survival Time Up                           0.828806  0.00635804    1.71063  0.127651   1                1469   0.833333  LOPEZ_MESOTELIOMA_SURVIVAL_TIME_UP
            1  Tonks Targets Of Runx1 Runx1T1 Fusion Hsc Up                -0.477485  0.00919118   -1.74962  0.139914   1                7592   0.6       TONKS_TARGETS_OF_RUNX1_RUNX1T1_FUSION_HSC_UP
            1  Tonks Targets Of Runx1 Runx1T1 Fusion Granulocyte Up        -0.463983  0.0471118    -1.55301  0.202734   1                9395   0.309524  TONKS_TARGETS_OF_RUNX1_RUNX1T1_FUSION_GRANULOCYTE_UP
            4  Papaspyridonos Unstable Aterosclerotic Plaque Up            -0.589416  0.0558648    -1.64566  0.224407   1                9352   0.456522  PAPASPYRIDONOS_UNSTABLE_ATEROSCLEROTIC_PLAQUE_UP
            1  Papaspyridonos Unstable Aterosclerotic Plaque Dn            -0.69794   0.00543588   -1.87286  0.119811   1                8138   0.868421  PAPASPYRIDONOS_UNSTABLE_ATEROSCLEROTIC_PLAQUE_DN
            1  Dittmer Pthlh Targets Up                                     0.380977  0.0192854     1.56238  0.160717   1                2168   0.37037   DITTMER_PTHLH_TARGETS_UP
           10  Dittmer Pthlh Targets Dn                                     0.406327  0.0204661     1.56512  0.244424   1                2109   0.3125    DITTMER_PTHLH_TARGETS_DN
            1  Udayakumar Med1 Targets Up                                   0.423956  0.00688538    1.65301  0.103136   1                1906   0.333333  UDAYAKUMAR_MED1_TARGETS_UP
           10  Udayakumar Med1 Targets Up                                   0.475731  0.000601926   1.85657  0.0652612  1                2665   0.455285  UDAYAKUMAR_MED1_TARGETS_UP
            1  Udayakumar Med1 Targets Dn                                  -0.367394  0.0125812    -1.5821   0.190354   1                8922   0.334951  UDAYAKUMAR_MED1_TARGETS_DN
            1  Odonnell Tfrc Targets Up                                    -0.429481  0.00301144   -1.78995  0.129551   1                8361   0.431907  ODONNELL_TFRC_TARGETS_UP
            1  Odonnell Tfrc Targets Dn                                     0.76481   0.00100543    2.08236  0.0127991  1                 299   0.522727  ODONNELL_TFRC_TARGETS_DN
            5  Odonnell Tfrc Targets Dn                                     0.748729  0.00237436    2.02307  0.0336414  1                 937   0.625     ODONNELL_TFRC_TARGETS_DN
           10  Odonnell Tfrc Targets Dn                                     0.722104  0.00604717    1.94911  0.056431   1                 698   0.511364  ODONNELL_TFRC_TARGETS_DN
            1  Odonnell Targets Of Myc And Tfrc Up                         -0.554628  0.0150512    -1.76802  0.138176   1                8638   0.5       ODONNELL_TARGETS_OF_MYC_AND_TFRC_UP
            1  Odonnell Targets Of Myc And Tfrc Dn                          0.890416  0.000200723   1.97764  0.013915   0.534123          299   0.705882  ODONNELL_TARGETS_OF_MYC_AND_TFRC_DN
            5  Odonnell Targets Of Myc And Tfrc Dn                          0.843415  0.00158165    1.85947  0.0620503  1                 823   0.764706  ODONNELL_TARGETS_OF_MYC_AND_TFRC_DN
           10  Odonnell Targets Of Myc And Tfrc Dn                          0.798216  0.0107935     1.76541  0.100873   1                1160   0.705882  ODONNELL_TARGETS_OF_MYC_AND_TFRC_DN
            1  Senese Hdac1 Targets Dn                                     -0.466555  0.0155515    -1.7371   0.141912   1                8846   0.453125  SENESE_HDAC1_TARGETS_DN
            1  Senese Hdac1 And Hdac2 Targets Dn                           -0.462909  0.0427732    -1.63086  0.170938   1                8587   0.48503   SENESE_HDAC1_AND_HDAC2_TARGETS_DN
            1  Senese Hdac2 Targets Dn                                     -0.500434  0.0616776    -1.59014  0.187154   1                8587   0.54878   SENESE_HDAC2_TARGETS_DN
            1  Lee Neural Crest Stem Cell Up                               -0.497804  0.0728892    -1.56348  0.199134   1                8660   0.5       LEE_NEURAL_CREST_STEM_CELL_UP
            1  Lee Neural Crest Stem Cell Dn                               -0.441673  0.025641     -1.59042  0.187154   1                8948   0.358974  LEE_NEURAL_CREST_STEM_CELL_DN
            2  Tien Intestine Probiotics 6Hr Up                            -0.705459  0.0042228    -1.98754  0.202024   1                9959   0.509434  TIEN_INTESTINE_PROBIOTICS_6HR_UP
            1  Zhou Inflammatory Response Lps Up                           -0.363863  0.0267786    -1.51542  0.222565   1                8716   0.334842  ZHOU_INFLAMMATORY_RESPONSE_LPS_UP
            1  Kinsey Targets Of Ewsr1 Flii Fusion Dn                      -0.541107  0.00222222   -1.92795  0.116729   1                8406   0.527273  KINSEY_TARGETS_OF_EWSR1_FLII_FUSION_DN
            1  Sabates Colorectal Adenoma Dn                               -0.583493  0.00785657   -1.8235   0.119811   1                9012   0.495798  SABATES_COLORECTAL_ADENOMA_DN
            1  Scibetta Kdm5B Targets Dn                                    0.491682  0.0165456     1.71937  0.0715452  1                 614   0.253968  SCIBETTA_KDM5B_TARGETS_DN
            5  Scibetta Kdm5B Targets Dn                                    0.503839  0.0102533     1.74864  0.115744   1                 430   0.269841  SCIBETTA_KDM5B_TARGETS_DN
           10  Scibetta Kdm5B Targets Dn                                    0.466042  0.0347513     1.62221  0.18979    1                 872   0.269841  SCIBETTA_KDM5B_TARGETS_DN
            1  Nagashima Nrg1 Signaling Up                                 -0.494969  0.0486974    -1.64011  0.170938   1                9477   0.38      NAGASHIMA_NRG1_SIGNALING_UP
            1  Nagashima Egf Signaling Up                                  -0.710589  0.0135903    -1.82424  0.119811   1                8652   0.745098  NAGASHIMA_EGF_SIGNALING_UP
            1  Kim Wt1 Targets Up                                          -0.415182  0.0672807    -1.5073   0.22774    1                9419   0.326316  KIM_WT1_TARGETS_UP
            1  Kim Wt1 Targets 12Hr Up                                     -0.403891  0.0482633    -1.51771  0.222374   1                7756   0.489209  KIM_WT1_TARGETS_12HR_UP
            1  Elvidge Hypoxia Dn                                           0.480354  0.0129766     1.70197  0.078785   1                2654   0.471014  ELVIDGE_HYPOXIA_DN
            4  Elvidge Hypoxia Dn                                          -0.488682  0.0129818    -1.72624  0.192508   1                8808   0.492754  ELVIDGE_HYPOXIA_DN
            1  Elvidge Hif1A And Hif2A Targets Up                           0.498231  0.0550959     1.52824  0.183841   1                2309   0.473684  ELVIDGE_HIF1A_AND_HIF2A_TARGETS_UP
            1  Morosetti Facioscapulohumeral Muscular Distrophy Up         -0.741352  0.0378236    -1.62858  0.170938   1                9351   0.733333  MOROSETTI_FACIOSCAPULOHUMERAL_MUSCULAR_DISTROPHY_UP
            1  Graham Cml Quiescent Vs Normal Quiescent Up                  0.572904  0.0219581     1.7611   0.0538537  1                1687   0.523077  GRAHAM_CML_QUIESCENT_VS_NORMAL_QUIESCENT_UP
            4  Graham Cml Quiescent Vs Normal Quiescent Up                 -0.638799  0.00256511   -1.95316  0.155097   1               10290   0.446154  GRAHAM_CML_QUIESCENT_VS_NORMAL_QUIESCENT_UP
            5  Graham Cml Quiescent Vs Normal Quiescent Up                  0.5215    0.0670635     1.5807   0.226952   1                 838   0.353846  GRAHAM_CML_QUIESCENT_VS_NORMAL_QUIESCENT_UP
           10  Graham Cml Quiescent Vs Normal Quiescent Up                  0.596991  0.0105431     1.82833  0.0688673  1                 700   0.369231  GRAHAM_CML_QUIESCENT_VS_NORMAL_QUIESCENT_UP
            1  Graham Cml Dividing Vs Normal Quiescent Up                   0.723085  0.00140817    1.994    0.0129733  1                 931   0.6       GRAHAM_CML_DIVIDING_VS_NORMAL_QUIESCENT_UP
            4  Graham Cml Dividing Vs Normal Quiescent Up                  -0.673579  0.0134202    -1.8494   0.174932   1                9637   0.6       GRAHAM_CML_DIVIDING_VS_NORMAL_QUIESCENT_UP
            5  Graham Cml Dividing Vs Normal Quiescent Up                   0.679232  0.0123506     1.85658  0.0622723  1                1107   0.525926  GRAHAM_CML_DIVIDING_VS_NORMAL_QUIESCENT_UP
           10  Graham Cml Dividing Vs Normal Quiescent Up                   0.718243  0.00240096    1.96996  0.056431   1                 766   0.555556  GRAHAM_CML_DIVIDING_VS_NORMAL_QUIESCENT_UP
            1  Graham Cml Dividing Vs Normal Quiescent Dn                  -0.539793  0.0224787    -1.71332  0.146665   1                8562   0.466667  GRAHAM_CML_DIVIDING_VS_NORMAL_QUIESCENT_DN
            1  Graham Normal Quiescent Vs Normal Dividing Up               -0.545979  0.0100361    -1.74165  0.140775   1                9255   0.4       GRAHAM_NORMAL_QUIESCENT_VS_NORMAL_DIVIDING_UP
            1  Graham Normal Quiescent Vs Normal Dividing Dn                0.861045  0.000404449   1.91285  0.0193532  1                 855   0.830769  GRAHAM_NORMAL_QUIESCENT_VS_NORMAL_DIVIDING_DN
            4  Graham Normal Quiescent Vs Normal Dividing Dn               -0.768587  0.0232101    -1.70911  0.192508   1               10178   0.692308  GRAHAM_NORMAL_QUIESCENT_VS_NORMAL_DIVIDING_DN
            5  Graham Normal Quiescent Vs Normal Dividing Dn                0.802291  0.00993049    1.77891  0.0972596  1                1078   0.769231  GRAHAM_NORMAL_QUIESCENT_VS_NORMAL_DIVIDING_DN
           10  Graham Normal Quiescent Vs Normal Dividing Dn                0.839157  0.00140704    1.86278  0.0652612  1                1194   0.830769  GRAHAM_NORMAL_QUIESCENT_VS_NORMAL_DIVIDING_DN
            1  Bidus Metastasis Up                                          0.642535  0.00201694    1.93207  0.0172116  1                2183   0.585492  BIDUS_METASTASIS_UP
            5  Bidus Metastasis Up                                          0.628238  0.00641283    1.87849  0.0582517  1                2081   0.544041  BIDUS_METASTASIS_UP
           10  Bidus Metastasis Up                                          0.605674  0.0114688     1.81124  0.0756726  1                2602   0.595855  BIDUS_METASTASIS_UP
            1  Bidus Metastasis Dn                                         -0.453917  0.0175055    -1.70302  0.147387   1                8187   0.467153  BIDUS_METASTASIS_DN
            1  Wamunyokoli Ovarian Cancer Lmp Dn                           -0.51373   0.0215465    -1.72902  0.144333   1                8290   0.544944  WAMUNYOKOLI_OVARIAN_CANCER_LMP_DN
            1  Wamunyokoli Ovarian Cancer Grades 1 2 Dn                    -0.617426  0.00343712   -1.85157  0.119811   1                9472   0.534483  WAMUNYOKOLI_OVARIAN_CANCER_GRADES_1_2_DN
            1  Mayburd Response To L663536 Dn                               0.538478  0.00725514    1.77326  0.0503723  1                1896   0.469388  MAYBURD_RESPONSE_TO_L663536_DN
            5  Mayburd Response To L663536 Dn                               0.501712  0.0215272     1.65111  0.177566   1                2365   0.44898   MAYBURD_RESPONSE_TO_L663536_DN
            1  Watanabe Ulcerative Colitis With Cancer Up                  -0.620956  0.0497271    -1.56291  0.199134   1               10641   0.4       WATANABE_ULCERATIVE_COLITIS_WITH_CANCER_UP
            1  Rodrigues Thyroid Carcinoma Dn                              -0.453884  0.0397076    -1.54078  0.209394   1                7514   0.485714  RODRIGUES_THYROID_CARCINOMA_DN
            1  Hahtola Sezary Syndrom Up                                    0.464344  0.0195801     1.64116  0.107415   1                1000   0.283784  HAHTOLA_SEZARY_SYNDROM_UP
            1  Hahtola Mycosis Fungoides Cd4 Up                            -0.508113  0.0365854    -1.63455  0.170938   1               10059   0.333333  HAHTOLA_MYCOSIS_FUNGOIDES_CD4_UP
            5  Hahtola Mycosis Fungoides Cd4 Dn                             0.411375  0.009342      1.64652  0.178715   1                2571   0.372727  HAHTOLA_MYCOSIS_FUNGOIDES_CD4_DN
            1  Provenzani Metastasis Up                                     0.387803  0.0530318     1.51337  0.19807    1                1991   0.378049  PROVENZANI_METASTASIS_UP
           10  Provenzani Metastasis Up                                     0.426065  0.0232097     1.65248  0.169206   1                2605   0.439024  PROVENZANI_METASTASIS_UP
            1  Liu Cdx2 Targets Up                                         -0.544151  0.054319     -1.49317  0.233901   1                8713   0.444444  LIU_CDX2_TARGETS_UP
            1  Kokkinakis Methionine Deprivation 96Hr Up                   -0.459553  0.0122926    -1.70865  0.146665   1                8717   0.471154  KOKKINAKIS_METHIONINE_DEPRIVATION_96HR_UP
            1  Kokkinakis Methionine Deprivation 96Hr Dn                    0.439654  0.0296754     1.57347  0.155488   1                2864   0.47619   KOKKINAKIS_METHIONINE_DEPRIVATION_96HR_DN
           10  Kokkinakis Methionine Deprivation 96Hr Dn                    0.462416  0.00983541    1.65808  0.164959   1                1886   0.349206  KOKKINAKIS_METHIONINE_DEPRIVATION_96HR_DN
            1  Chow Rassf1 Targets Up                                       0.551289  0.0205745     1.64917  0.103998   1                 960   0.307692  CHOW_RASSF1_TARGETS_UP
            5  Chow Rassf1 Targets Up                                       0.527269  0.0357002     1.58038  0.226952   1                3518   0.730769  CHOW_RASSF1_TARGETS_UP
            1  Landis Breast Cancer Progression Dn                         -0.660298  0.000821187  -2.00033  0.116729   1               10179   0.459016  LANDIS_BREAST_CANCER_PROGRESSION_DN
            1  Geserick Tert Targets Dn                                    -0.779725  0.00646073   -1.83363  0.119811   1                9760   0.7       GESERICK_TERT_TARGETS_DN
            1  Delys Thyroid Cancer Dn                                     -0.56895   0.00503423   -1.88271  0.119811   1                9236   0.5       DELYS_THYROID_CANCER_DN
            1  Chiaradonna Neoplastic Transformation Kras Cdc25 Dn         -0.703369  0.00305748   -1.89016  0.119811   1                9239   0.634146  CHIARADONNA_NEOPLASTIC_TRANSFORMATION_KRAS_CDC25_DN
            1  Chiaradonna Neoplastic Transformation Kras Up                0.450146  0.0360762     1.56938  0.157043   1                1507   0.364486  CHIARADONNA_NEOPLASTIC_TRANSFORMATION_KRAS_UP
            4  Chiaradonna Neoplastic Transformation Kras Up               -0.50302   0.00733109   -1.74724  0.192508   1                8413   0.514019  CHIARADONNA_NEOPLASTIC_TRANSFORMATION_KRAS_UP
            1  Chiaradonna Neoplastic Transformation Kras Dn               -0.541951  0.00263051   -1.91922  0.116729   1                9476   0.414634  CHIARADONNA_NEOPLASTIC_TRANSFORMATION_KRAS_DN
            1  Chiaradonna Neoplastic Transformation Cdc25 Up              -0.41062   0.0724784    -1.47851  0.242484   1                9271   0.371429  CHIARADONNA_NEOPLASTIC_TRANSFORMATION_CDC25_UP
            1  Chiaradonna Neoplastic Transformation Cdc25 Dn              -0.39423   0.0404959    -1.51673  0.222374   1                9900   0.261905  CHIARADONNA_NEOPLASTIC_TRANSFORMATION_CDC25_DN
            1  Chebotaev Gr Targets Up                                     -0.487448  0.0166533    -1.66032  0.162756   1                8226   0.517857  CHEBOTAEV_GR_TARGETS_UP
            1  Chebotaev Gr Targets Dn                                     -0.566757  0.0210569    -1.74819  0.139914   1                8994   0.535211  CHEBOTAEV_GR_TARGETS_DN
            1  Berenjeno Transformed By Rhoa Reversibly Dn                 -0.688577  0.00362319   -1.85454  0.119811   1                9051   0.666667  BERENJENO_TRANSFORMED_BY_RHOA_REVERSIBLY_DN
            1  Lindgren Bladder Cancer Cluster 2A Dn                       -0.472217  0.0257536    -1.66433  0.162484   1                8538   0.436508  LINDGREN_BLADDER_CANCER_CLUSTER_2A_DN
            1  Lindgren Bladder Cancer Cluster 3 Up                         0.671432  0.000404122   2.04764  0.0127991  1                1467   0.517986  LINDGREN_BLADDER_CANCER_CLUSTER_3_UP
            4  Lindgren Bladder Cancer Cluster 3 Up                        -0.543586  0.0555446    -1.65523  0.217692   1                8840   0.510791  LINDGREN_BLADDER_CANCER_CLUSTER_3_UP
            5  Lindgren Bladder Cancer Cluster 3 Up                         0.612071  0.00916335    1.86174  0.0620503  1                1484   0.438849  LINDGREN_BLADDER_CANCER_CLUSTER_3_UP
           10  Lindgren Bladder Cancer Cluster 3 Up                         0.661924  0.000796972   2.03185  0.056431   1                1575   0.52518   LINDGREN_BLADDER_CANCER_CLUSTER_3_UP
            1  Lindgren Bladder Cancer Cluster 3 Dn                        -0.411676  0.0151945    -1.63379  0.170938   1                8622   0.4       LINDGREN_BLADDER_CANCER_CLUSTER_3_DN
            1  Creighton Akt1 Signaling Via Mtor Dn                         0.603609  0.0724782     1.54172  0.174598   1                2371   0.5       CREIGHTON_AKT1_SIGNALING_VIA_MTOR_DN
            4  Creighton Akt1 Signaling Via Mtor Dn                        -0.678178  0.0167271    -1.73054  0.192508   1                9064   0.636364  CREIGHTON_AKT1_SIGNALING_VIA_MTOR_DN
            1  Naderi Breast Cancer Prognosis Up                            0.829326  0.00101215    1.84518  0.0311256  1                 894   0.71875   NADERI_BREAST_CANCER_PROGNOSIS_UP
            5  Naderi Breast Cancer Prognosis Up                            0.772972  0.0112894     1.71932  0.12811    1                 633   0.59375   NADERI_BREAST_CANCER_PROGNOSIS_UP
           10  Naderi Breast Cancer Prognosis Up                            0.810478  0.00161878    1.7999   0.0809518  1                1706   0.84375   NADERI_BREAST_CANCER_PROGNOSIS_UP
            1  Naderi Breast Cancer Prognosis Dn                           -0.64659   0.0434171    -1.59825  0.181764   1                8576   0.6875    NADERI_BREAST_CANCER_PROGNOSIS_DN
            1  Markey Rb1 Chronic Lof Up                                    0.454296  0.0178715     1.6646   0.0969555  1                 983   0.369565  MARKEY_RB1_CHRONIC_LOF_UP
            5  Markey Rb1 Chronic Lof Up                                    0.445063  0.025874      1.62044  0.193542   1                1295   0.336957  MARKEY_RB1_CHRONIC_LOF_UP
           10  Markey Rb1 Chronic Lof Up                                    0.45485   0.0186037     1.65497  0.166946   1                1220   0.369565  MARKEY_RB1_CHRONIC_LOF_UP
            1  Markey Rb1 Acute Lof Up                                      0.652024  0.000805477   2.02963  0.0127991  1                1196   0.530612  MARKEY_RB1_ACUTE_LOF_UP
            5  Markey Rb1 Acute Lof Up                                      0.645322  0.00159776    1.99889  0.0354642  1                1737   0.571429  MARKEY_RB1_ACUTE_LOF_UP
           10  Markey Rb1 Acute Lof Up                                      0.622691  0.00435471    1.95205  0.056431   1                1588   0.535714  MARKEY_RB1_ACUTE_LOF_UP
            1  Chin Breast Cancer Copy Number Up                            0.510818  0.0928387     1.44541  0.238115   1                1804   0.363636  CHIN_BREAST_CANCER_COPY_NUMBER_UP
            1  Vanharanta Uterine Fibroid Dn                               -0.552693  0.0172691    -1.7466   0.139914   1                8665   0.576271  VANHARANTA_UTERINE_FIBROID_DN
            1  Vanharanta Uterine Fibroid With 7Q Deletion Up               0.524178  0.00243161    1.78314  0.0476102  1                2051   0.435484  VANHARANTA_UTERINE_FIBROID_WITH_7Q_DELETION_UP
            5  Vanharanta Uterine Fibroid With 7Q Deletion Up               0.476608  0.0195999     1.63146  0.186277   1                3721   0.596774  VANHARANTA_UTERINE_FIBROID_WITH_7Q_DELETION_UP
           10  Vanharanta Uterine Fibroid With 7Q Deletion Up               0.486502  0.0157151     1.66433  0.162408   1                2409   0.403226  VANHARANTA_UTERINE_FIBROID_WITH_7Q_DELETION_UP
            1  Landis Erbb2 Breast Tumors 65 Dn                            -0.547985  0.0256721    -1.69806  0.148746   1                9888   0.363636  LANDIS_ERBB2_BREAST_TUMORS_65_DN
            1  Landis Erbb2 Breast Preneoplastic Up                         0.508404  0.070186      1.45092  0.238115   1                 997   0.333333  LANDIS_ERBB2_BREAST_PRENEOPLASTIC_UP
            1  Concannon Apoptosis By Epoxomicin Up                        -0.35586   0.0364919    -1.50052  0.230381   1                9745   0.248804  CONCANNON_APOPTOSIS_BY_EPOXOMICIN_UP
            1  Concannon Apoptosis By Epoxomicin Dn                         0.411179  0.044         1.54733  0.171988   1                1321   0.286822  CONCANNON_APOPTOSIS_BY_EPOXOMICIN_DN
            1  Gaussmann Mll Af4 Fusion Targets A Dn                       -0.428588  0.0162009    -1.60128  0.181764   1                8344   0.423729  GAUSSMANN_MLL_AF4_FUSION_TARGETS_A_DN
            1  Gaussmann Mll Af4 Fusion Targets C Up                       -0.349406  0.0322646    -1.47202  0.246864   1                7668   0.492647  GAUSSMANN_MLL_AF4_FUSION_TARGETS_C_UP
            1  Gaussmann Mll Af4 Fusion Targets E Up                       -0.557414  0.00886203   -1.80033  0.124176   1                9012   0.479452  GAUSSMANN_MLL_AF4_FUSION_TARGETS_E_UP
            1  Gaussmann Mll Af4 Fusion Targets F Up                       -0.519915  0.0102471    -1.7997   0.124176   1                9002   0.437086  GAUSSMANN_MLL_AF4_FUSION_TARGETS_F_UP
            1  Berenjeno Transformed By Rhoa Forever Up                    -0.666392  0.0325779    -1.62967  0.170938   1                8601   0.666667  BERENJENO_TRANSFORMED_BY_RHOA_FOREVER_UP
            1  Berenjeno Transformed By Rhoa Forever Dn                    -0.51677   0.0561752    -1.52846  0.214812   1                8713   0.56      BERENJENO_TRANSFORMED_BY_RHOA_FOREVER_DN
            5  Cairo Pml Targets Bound By Myc Up                            0.595286  0.0245656     1.6382   0.184004   1                2330   0.7       CAIRO_PML_TARGETS_BOUND_BY_MYC_UP
            1  Cairo Pml Targets Bound By Myc Dn                           -0.598512  0.0689585    -1.5068   0.22774    1               10477   0.333333  CAIRO_PML_TARGETS_BOUND_BY_MYC_DN
            1  Ouellet Ovarian Cancer Invasive Vs Lmp Up                    0.671483  0.000402334   2.00475  0.0127991  1                2060   0.630631  OUELLET_OVARIAN_CANCER_INVASIVE_VS_LMP_UP
            4  Ouellet Ovarian Cancer Invasive Vs Lmp Up                   -0.632919  0.00178962   -1.88515  0.162281   1                9114   0.531532  OUELLET_OVARIAN_CANCER_INVASIVE_VS_LMP_UP
           10  Ouellet Ovarian Cancer Invasive Vs Lmp Up                    0.576855  0.021544      1.72452  0.117315   1                3004   0.648649  OUELLET_OVARIAN_CANCER_INVASIVE_VS_LMP_UP
            1  Landis Erbb2 Breast Tumors 324 Dn                           -0.462528  0.0183673    -1.71265  0.146665   1                9888   0.328244  LANDIS_ERBB2_BREAST_TUMORS_324_DN
            1  Baris Thyroid Cancer Up                                      0.529913  0.0755094     1.47007  0.22633    1                1929   0.473684  BARIS_THYROID_CANCER_UP
            1  Wang Barretts Esophagus And Esophagus Cancer Up             -0.553253  0.0615694    -1.51634  0.222374   1                8759   0.5       WANG_BARRETTS_ESOPHAGUS_AND_ESOPHAGUS_CANCER_UP
            1  Kerley Response To Cisplatin Up                             -0.670164  0.00672235   -1.85131  0.119811   1                8961   0.648649  KERLEY_RESPONSE_TO_CISPLATIN_UP
            1  Missiaglia Regulated By Methylation Dn                       0.759336  0.000404204   2.02836  0.0127991  1                 665   0.59434   MISSIAGLIA_REGULATED_BY_METHYLATION_DN
            5  Missiaglia Regulated By Methylation Dn                       0.720603  0.00218993    1.92022  0.0475578  1                1614   0.688679  MISSIAGLIA_REGULATED_BY_METHYLATION_DN
           10  Missiaglia Regulated By Methylation Dn                       0.718668  0.00159585    1.9256   0.056431   1                1189   0.603774  MISSIAGLIA_REGULATED_BY_METHYLATION_DN
            1  Tang Senescence Tp53 Targets Up                             -0.557106  0.037089     -1.58312  0.190354   1                9104   0.55      TANG_SENESCENCE_TP53_TARGETS_UP
            1  Tang Senescence Tp53 Targets Dn                              0.833028  0.00241838    1.90283  0.0208664  1                 272   0.666667  TANG_SENESCENCE_TP53_TARGETS_DN
            4  Tang Senescence Tp53 Targets Dn                             -0.751396  0.0377844    -1.70717  0.192508   1               10338   0.622222  TANG_SENESCENCE_TP53_TARGETS_DN
            5  Tang Senescence Tp53 Targets Dn                              0.857057  0.00059512    1.93724  0.0454213  1                 575   0.688889  TANG_SENESCENCE_TP53_TARGETS_DN
           10  Tang Senescence Tp53 Targets Dn                              0.796754  0.0138861     1.80937  0.0757489  1                 774   0.688889  TANG_SENESCENCE_TP53_TARGETS_DN
            1  Roylance Breast Cancer 16Q Copy Number Up                    0.62922   0.0882709     1.55718  0.164311   1                1684   0.527778  ROYLANCE_BREAST_CANCER_16Q_COPY_NUMBER_UP
           10  Roylance Breast Cancer 16Q Copy Number Up                    0.708424  0.0266213     1.74965  0.108321   1                 635   0.5       ROYLANCE_BREAST_CANCER_16Q_COPY_NUMBER_UP
            1  Johansson Gliomagenesis By Pdgfb Up                          0.443725  0.0974039     1.43805  0.244682   1                1355   0.327273  JOHANSSON_GLIOMAGENESIS_BY_PDGFB_UP
            4  Johansson Gliomagenesis By Pdgfb Up                         -0.540245  0.00856062   -1.74462  0.192508   1                9625   0.418182  JOHANSSON_GLIOMAGENESIS_BY_PDGFB_UP
            1  Dunne Targets Of Aml1 Mtg8 Fusion Dn                        -0.694276  0.0196635    -1.69462  0.150339   1                8638   0.6875    DUNNE_TARGETS_OF_AML1_MTG8_FUSION_DN
            1  Ouellet Cultured Ovarian Cancer Invasive Vs Lmp Dn          -0.575225  0.0247258    -1.66315  0.162484   1                9502   0.428571  OUELLET_CULTURED_OVARIAN_CANCER_INVASIVE_VS_LMP_DN
            1  Mcbryan Pubertal Tgfb1 Targets Up                           -0.430269  0.0642497    -1.55016  0.202734   1                8642   0.438272  MCBRYAN_PUBERTAL_TGFB1_TARGETS_UP
            1  Zirn Tretinoin Response Up                                  -0.563357  0.0699414    -1.51177  0.225344   1                8779   0.611111  ZIRN_TRETINOIN_RESPONSE_UP
            1  Darwiche Papilloma Risk Low Up                               0.382349  0.0528652     1.45285  0.237504   1                1935   0.364486  DARWICHE_PAPILLOMA_RISK_LOW_UP
            4  Darwiche Papilloma Risk High Up                             -0.438829  0.014875     -1.61617  0.238897   1                9474   0.352941  DARWICHE_PAPILLOMA_RISK_HIGH_UP
            4  Darwiche Squamous Cell Carcinoma Up                         -0.428444  0.00851317   -1.65069  0.220523   1                9984   0.32381   DARWICHE_SQUAMOUS_CELL_CARCINOMA_UP
            1  Tomida Metastasis Up                                         0.605712  0.053694      1.57153  0.156894   1                1606   0.45      TOMIDA_METASTASIS_UP
            4  Tomida Metastasis Up                                        -0.669883  0.010429     -1.73633  0.192508   1                8254   0.75      TOMIDA_METASTASIS_UP
            1  Hu Angiogenesis Up                                          -0.603789  0.00888171   -1.74185  0.140775   1               10308   0.333333  HU_ANGIOGENESIS_UP
            1  Hu Angiogenesis Dn                                           0.575777  0.0194371     1.69031  0.0843017  1                2546   0.628571  HU_ANGIOGENESIS_DN
            4  Hu Angiogenesis Dn                                          -0.554753  0.0307692    -1.63496  0.229428   1                9756   0.4       HU_ANGIOGENESIS_DN
            5  Hu Angiogenesis Dn                                           0.549798  0.0343088     1.63375  0.18509    1                1951   0.428571  HU_ANGIOGENESIS_DN
            1  Begum Targets Of Pax3 Foxo1 Fusion Up                       -0.566292  0.0275061    -1.67209  0.161802   1                7984   0.680851  BEGUM_TARGETS_OF_PAX3_FOXO1_FUSION_UP
            1  Begum Targets Of Pax3 Foxo1 Fusion Dn                       -0.573382  0.0436637    -1.62829  0.170938   1                8382   0.567568  BEGUM_TARGETS_OF_PAX3_FOXO1_FUSION_DN
            1  Olsson E2F3 Targets Dn                                       0.563951  0.0658238     1.59675  0.138478   1                 589   0.4       OLSSON_E2F3_TARGETS_DN
            4  Olsson E2F3 Targets Dn                                      -0.620563  0.0168117    -1.75178  0.192508   1               10551   0.342857  OLSSON_E2F3_TARGETS_DN
            5  Olsson E2F3 Targets Dn                                       0.557236  0.0645417     1.57746  0.229377   1                 910   0.428571  OLSSON_E2F3_TARGETS_DN
            5  Marks Hdac Targets Dn                                        0.609347  0.0326547     1.58144  0.226952   1                1925   0.454545  MARKS_HDAC_TARGETS_DN
            5  Choi Atl Stage Predictor                                     0.51741   0.0074178     1.73942  0.117295   1                2702   0.4       CHOI_ATL_STAGE_PREDICTOR
            1  Wong Endmetrium Cancer Dn                                   -0.767239  0.0124795    -1.78053  0.134462   1                9014   0.877551  WONG_ENDMETRIUM_CANCER_DN
            1  Kong E2F3 Targets                                            0.890171  0.000406091   1.97498  0.013915   1                 504   0.777778  KONG_E2F3_TARGETS
            4  Kong E2F3 Targets                                           -0.773936  0.0268786    -1.72112  0.192508   1                9689   0.722222  KONG_E2F3_TARGETS
            5  Kong E2F3 Targets                                            0.872585  0.000594295   1.9354   0.0454213  1                 625   0.819444  KONG_E2F3_TARGETS
           10  Kong E2F3 Targets                                            0.826391  0.00383528    1.83477  0.0688673  1                 921   0.763889  KONG_E2F3_TARGETS
            1  Eguchi Cell Cycle Rb1 Targets                                0.959741  0.000406669   1.75508  0.055805   1                 467   1         EGUCHI_CELL_CYCLE_RB1_TARGETS
            5  Eguchi Cell Cycle Rb1 Targets                                0.939836  0.00199402    1.7192   0.12811    1                 688   1         EGUCHI_CELL_CYCLE_RB1_TARGETS
           10  Eguchi Cell Cycle Rb1 Targets                                0.938755  0.00100422    1.72744  0.116755   1                 700   1         EGUCHI_CELL_CYCLE_RB1_TARGETS
            1  Ivanov Mutated In Colon Cancer                               0.582663  0.0539007     1.52842  0.183841   1                2459   0.461538  IVANOV_MUTATED_IN_COLON_CANCER
            1  Riz Erythroid Differentiation                                0.567144  0.00548669    1.86733  0.0267786  1                1887   0.444444  RIZ_ERYTHROID_DIFFERENTIATION
            5  Riz Erythroid Differentiation                                0.606241  0.00120072    1.98959  0.037673   1                2688   0.587302  RIZ_ERYTHROID_DIFFERENTIATION
           10  Riz Erythroid Differentiation                                0.545867  0.0099553     1.78316  0.0904153  1                2602   0.52381   RIZ_ERYTHROID_DIFFERENTIATION
            5  Riz Erythroid Differentiation Ccne1                          0.515231  0.00603743    1.735    0.121243   1                1818   0.433333  RIZ_ERYTHROID_DIFFERENTIATION_CCNE1
            1  Chassot Skin Wound                                          -0.888863  0.00623241   -1.7147   0.146665   1               10158   0.9       CHASSOT_SKIN_WOUND
            1  Riz Erythroid Differentiation Apobec2                       -0.568206  0.0701195    -1.48311  0.241988   1               10316   0.357143  RIZ_ERYTHROID_DIFFERENTIATION_APOBEC2
            1  Riz Erythroid Differentiation 12Hr                          -0.68199   0.000403633  -1.99767  0.116729   1                9502   0.580645  RIZ_ERYTHROID_DIFFERENTIATION_12HR
            1  Ebauer Targets Of Pax3 Foxo1 Fusion Up                      -0.420628  0.0202634    -1.6236   0.170938   1                8671   0.395349  EBAUER_TARGETS_OF_PAX3_FOXO1_FUSION_UP
            1  Perez Tp63 Targets                                          -0.410362  0.0131658    -1.64711  0.169344   1                7973   0.432314  PEREZ_TP63_TARGETS
            1  Perez Tp53 And Tp63 Targets                                 -0.487288  0.00359569   -1.83086  0.119811   1                7927   0.496296  PEREZ_TP53_AND_TP63_TARGETS
            1  Wakasugi Have Znf143 Binding Sites                           0.522243  0.027164      1.64027  0.107592   1                1546   0.375     WAKASUGI_HAVE_ZNF143_BINDING_SITES
            5  Wakasugi Have Znf143 Binding Sites                           0.606587  0.00098912    1.91058  0.0488819  1                2462   0.520833  WAKASUGI_HAVE_ZNF143_BINDING_SITES
           10  Wakasugi Have Znf143 Binding Sites                           0.519525  0.0289064     1.63145  0.181196   1                2710   0.479167  WAKASUGI_HAVE_ZNF143_BINDING_SITES
            1  Rodrigues Ntn1 And Dcc Targets                              -0.57208   0.0101358    -1.74839  0.139914   1                8540   0.586207  RODRIGUES_NTN1_AND_DCC_TARGETS
            1  Dacosta Uv Response Via Ercc3 Ttd Up                         0.462148  0.0680258     1.49182  0.213453   1                2030   0.436364  DACOSTA_UV_RESPONSE_VIA_ERCC3_TTD_UP
            4  Dacosta Uv Response Via Ercc3 Ttd Up                        -0.497871  0.0257183    -1.62902  0.230529   1                8334   0.563636  DACOSTA_UV_RESPONSE_VIA_ERCC3_TTD_UP
            1  Mahadevan Response To Mp470 Dn                              -0.772722  0.00634596   -1.86017  0.119811   1                9744   0.529412  MAHADEVAN_RESPONSE_TO_MP470_DN
            1  Appierto Response To Fenretinide Dn                          0.432053  0.0782713     1.4485   0.238115   1                1394   0.395349  APPIERTO_RESPONSE_TO_FENRETINIDE_DN
            1  Aung Gastric Cancer                                          0.442226  0.0227453     1.55805  0.16398    1                 634   0.25      AUNG_GASTRIC_CANCER
            5  Aung Gastric Cancer                                          0.478652  0.00854361    1.68099  0.151662   1                1844   0.361111  AUNG_GASTRIC_CANCER
            1  Breuhahn Growth Factor Signaling In Liver Cancer            -0.667314  0.00243605   -1.83829  0.119811   1                9983   0.47619   BREUHAHN_GROWTH_FACTOR_SIGNALING_IN_LIVER_CANCER
            1  Iwanaga E2F1 Targets Induced By Serum                        0.67607   0.00428397    1.82976  0.0345609  1                 625   0.375     IWANAGA_E2F1_TARGETS_INDUCED_BY_SERUM
            5  Iwanaga E2F1 Targets Induced By Serum                        0.728876  0.00139526    1.9667   0.0383868  1                 959   0.583333  IWANAGA_E2F1_TARGETS_INDUCED_BY_SERUM
           10  Iwanaga E2F1 Targets Induced By Serum                        0.671708  0.00650803    1.81079  0.0756726  1                1459   0.5       IWANAGA_E2F1_TARGETS_INDUCED_BY_SERUM
            1  Lastowska Coamplified With Mycn                              0.536636  0.147154      1.43189  0.249302   1                2034   0.448276  LASTOWSKA_COAMPLIFIED_WITH_MYCN
            5  Lastowska Coamplified With Mycn                              0.608352  0.0504669     1.6206   0.193542   1                2958   0.655172  LASTOWSKA_COAMPLIFIED_WITH_MYCN
            1  Li Lung Cancer                                               0.525867  0.0354496     1.60551  0.131322   1                1179   0.333333  LI_LUNG_CANCER
            1  Li Amplified In Lung Cancer                                  0.482065  0.0770903     1.54157  0.174598   1                2353   0.448052  LI_AMPLIFIED_IN_LUNG_CANCER
            4  Li Amplified In Lung Cancer                                 -0.581354  0.00680953   -1.86819  0.162281   1                8868   0.519481  LI_AMPLIFIED_IN_LUNG_CANCER
            1  Ebauer Myogenic Targets Of Pax3 Foxo1 Fusion                -0.545372  0.0566153    -1.55841  0.201403   1                9715   0.444444  EBAUER_MYOGENIC_TARGETS_OF_PAX3_FOXO1_FUSION
            1  Mattioli Mgus Vs Pcl                                         0.605691  0.00201572    1.92679  0.0180334  1                1674   0.526882  MATTIOLI_MGUS_VS_PCL
            4  Mattioli Mgus Vs Pcl                                        -0.527846  0.0225358    -1.6814   0.203728   1                9359   0.462366  MATTIOLI_MGUS_VS_PCL
            5  Mattioli Mgus Vs Pcl                                         0.510732  0.0352564     1.62543  0.19014    1                2470   0.494624  MATTIOLI_MGUS_VS_PCL
           10  Mattioli Mgus Vs Pcl                                         0.490114  0.0542996     1.5665   0.243334   1                2654   0.505376  MATTIOLI_MGUS_VS_PCL
            1  Lui Thyroid Cancer Cluster 1                                -0.419016  0.0469465    -1.49545  0.232098   1                8642   0.404255  LUI_THYROID_CANCER_CLUSTER_1
            1  Lui Thyroid Cancer Cluster 2                                -0.560327  0.0116851    -1.73209  0.142218   1                9780   0.371429  LUI_THYROID_CANCER_CLUSTER_2
            2  Lui Thyroid Cancer Cluster 3                                -0.754822  0.00240336   -2.00051  0.202024   1               10184   0.555556  LUI_THYROID_CANCER_CLUSTER_3
            1  Grasemann Retinoblastoma With 6P Amplification               0.722611  0.0272524     1.62905  0.112797   1                1283   0.666667  GRASEMANN_RETINOBLASTOMA_WITH_6P_AMPLIFICATION
            5  Grasemann Retinoblastoma With 6P Amplification               0.791499  0.00501605    1.77584  0.0983477  1                1095   0.583333  GRASEMANN_RETINOBLASTOMA_WITH_6P_AMPLIFICATION
            1  Schlosser Myc And Serum Response Synergy                     0.563286  0.0589528     1.56873  0.157043   1                3368   0.75      SCHLOSSER_MYC_AND_SERUM_RESPONSE_SYNERGY
            4  Schlosser Myc And Serum Response Synergy                    -0.587333  0.0366667    -1.63236  0.230376   1                9676   0.535714  SCHLOSSER_MYC_AND_SERUM_RESPONSE_SYNERGY
            1  Schlosser Myc Targets Repressed By Serum                     0.586818  0.00542605    1.87625  0.0247325  1                2413   0.557047  SCHLOSSER_MYC_TARGETS_REPRESSED_BY_SERUM
            5  Schlosser Myc Targets Repressed By Serum                     0.529806  0.0284398     1.70668  0.135933   1                3063   0.583893  SCHLOSSER_MYC_TARGETS_REPRESSED_BY_SERUM
           10  Schlosser Myc Targets Repressed By Serum                     0.52047   0.0327149     1.66856  0.159599   1                2622   0.52349   SCHLOSSER_MYC_TARGETS_REPRESSED_BY_SERUM
            1  Schlosser Myc Targets And Serum Response Dn                  0.693757  0.00485928    1.81991  0.0363179  1                1885   0.613636  SCHLOSSER_MYC_TARGETS_AND_SERUM_RESPONSE_DN
           10  Schlosser Myc Targets And Serum Response Dn                  0.628007  0.0393449     1.66395  0.162408   1                3060   0.681818  SCHLOSSER_MYC_TARGETS_AND_SERUM_RESPONSE_DN
            1  Schlosser Myc Targets And Serum Response Up                  0.647976  0.00279051    1.87461  0.024989   1                1661   0.545455  SCHLOSSER_MYC_TARGETS_AND_SERUM_RESPONSE_UP
            5  Schlosser Myc Targets And Serum Response Up                  0.556786  0.0525898     1.59956  0.209579   1                2761   0.522727  SCHLOSSER_MYC_TARGETS_AND_SERUM_RESPONSE_UP
            1  Mattioli Multiple Myeloma Subgroups                          0.54888   0.0507797     1.49542  0.212964   1                 894   0.454545  MATTIOLI_MULTIPLE_MYELOMA_SUBGROUPS
            1  Farmer Breast Cancer Cluster 2                               0.863154  0.0114137     1.63192  0.111343   1                1351   0.931034  FARMER_BREAST_CANCER_CLUSTER_2
           10  Farmer Breast Cancer Cluster 2                               0.847023  0.0224269     1.59721  0.211912   1                 327   0.724138  FARMER_BREAST_CANCER_CLUSTER_2
            1  Rosty Cervical Cancer Proliferation Cluster                  0.895266  0.000403877   1.92047  0.0188436  1                 682   0.843478  ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER
            4  Rosty Cervical Cancer Proliferation Cluster                 -0.832021  0.00488854   -1.79279  0.192508   1               10308   0.747826  ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER
            5  Rosty Cervical Cancer Proliferation Cluster                  0.837549  0.00656716    1.78993  0.0965026  1                 997   0.773913  ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER
           10  Rosty Cervical Cancer Proliferation Cluster                  0.850878  0.00300782    1.82872  0.0688673  1                 717   0.782609  ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER
            1  Kobayashi Response To Romidepsin                            -0.660305  0.0105725    -1.72492  0.146219   1               10096   0.357143  KOBAYASHI_RESPONSE_TO_ROMIDEPSIN
            1  La Men1 Targets                                             -0.549993  0.0687313    -1.49749  0.231875   1                8948   0.526316  LA_MEN1_TARGETS
            1  Honrado Breast Cancer Brca1 Vs Brca2                         0.626799  0.0208123     1.7032   0.0785025  1                 180   0.352941  HONRADO_BREAST_CANCER_BRCA1_VS_BRCA2
            5  Honrado Breast Cancer Brca1 Vs Brca2                         0.673661  0.00597253    1.82794  0.0744562  1                 547   0.411765  HONRADO_BREAST_CANCER_BRCA1_VS_BRCA2
           10  Honrado Breast Cancer Brca1 Vs Brca2                         0.606634  0.0284284     1.64479  0.17341    1                 622   0.411765  HONRADO_BREAST_CANCER_BRCA1_VS_BRCA2
            1  Hernandez Aberrant Mitosis By Docetacel 4Nm Up              -0.663499  0.0112744    -1.75553  0.139914   1               10198   0.444444  HERNANDEZ_ABERRANT_MITOSIS_BY_DOCETACEL_4NM_UP
            1  Hernandez Mitotic Arrest By Docetaxel 1 Up                  -0.509467  0.0744298    -1.47978  0.242366   1                9936   0.407407  HERNANDEZ_MITOTIC_ARREST_BY_DOCETAXEL_1_UP
            1  Hernandez Aberrant Mitosis By Docetacel 2Nm Up              -0.479564  0.0185334    -1.65122  0.168175   1                8564   0.5       HERNANDEZ_ABERRANT_MITOSIS_BY_DOCETACEL_2NM_UP
            1  Hernandez Mitotic Arrest By Docetaxel 2 Dn                  -0.581032  0.0577976    -1.51533  0.222565   1                9589   0.466667  HERNANDEZ_MITOTIC_ARREST_BY_DOCETAXEL_2_DN
            1  Xu Hgf Signaling Not Via Akt1 48Hr Dn                        0.672659  0.0125328     1.76186  0.0537552  1                 102   0.333333  XU_HGF_SIGNALING_NOT_VIA_AKT1_48HR_DN
            4  Xu Hgf Signaling Not Via Akt1 48Hr Dn                       -0.729088  0.00255956   -1.90809  0.155097   1               10457   0.388889  XU_HGF_SIGNALING_NOT_VIA_AKT1_48HR_DN
            5  Xu Hgf Signaling Not Via Akt1 48Hr Dn                        0.648691  0.0237288     1.69339  0.143083   1                 425   0.388889  XU_HGF_SIGNALING_NOT_VIA_AKT1_48HR_DN
           10  Xu Hgf Signaling Not Via Akt1 48Hr Dn                        0.630356  0.0358155     1.64062  0.174876   1                 212   0.333333  XU_HGF_SIGNALING_NOT_VIA_AKT1_48HR_DN
            1  Xu Hgf Targets Induced By Akt1 48Hr Dn                       0.614831  0.0315233     1.6244   0.115882   1                 536   0.333333  XU_HGF_TARGETS_INDUCED_BY_AKT1_48HR_DN
            5  Xu Hgf Targets Induced By Akt1 48Hr Dn                       0.676795  0.0061374     1.78246  0.0972596  1                1005   0.444444  XU_HGF_TARGETS_INDUCED_BY_AKT1_48HR_DN
            1  Shi Sparc Targets Up                                        -0.541839  0.0288132    -1.58939  0.187154   1                8232   0.571429  SHI_SPARC_TARGETS_UP
            1  Waesch Anaphase Promoting Complex                            0.626554  0.0406732     1.53023  0.183101   1                2052   0.5       WAESCH_ANAPHASE_PROMOTING_COMPLEX
            5  Waesch Anaphase Promoting Complex                            0.666388  0.0238381     1.61594  0.198959   1                 596   0.3       WAESCH_ANAPHASE_PROMOTING_COMPLEX
           10  Waesch Anaphase Promoting Complex                            0.694757  0.00881587    1.70513  0.13013    1                2382   0.6       WAESCH_ANAPHASE_PROMOTING_COMPLEX
            1  Neben Aml With Flt3 Or Nras Dn                               0.704086  0.0198596     1.64448  0.106001   1                  61   0.272727  NEBEN_AML_WITH_FLT3_OR_NRAS_DN
            1  Amundson Dna Damage Response Tp53                           -0.721171  0.00983936   -1.79067  0.129551   1               10694   0.5       AMUNDSON_DNA_DAMAGE_RESPONSE_TP53
            1  Dauer Stat3 Targets Up                                      -0.644463  0.011876     -1.8031   0.124176   1                9421   0.560976  DAUER_STAT3_TARGETS_UP
            1  Teramoto Opn Targets Cluster 6                              -0.611005  0.00222267   -1.85359  0.119811   1                9859   0.391304  TERAMOTO_OPN_TARGETS_CLUSTER_6
            1  Kang Gist With Pdgfra Up                                    -0.571749  0.0464695    -1.62813  0.170938   1                9010   0.487179  KANG_GIST_WITH_PDGFRA_UP
            1  Meinhold Ovarian Cancer Low Grade Dn                         0.66866   0.0192926     1.6837   0.0868165  1                2049   0.684211  MEINHOLD_OVARIAN_CANCER_LOW_GRADE_DN
            1  Dacosta Ercc3 Allele Xpcs Vs Ttd Dn                         -0.604902  0.0142199    -1.75548  0.139914   1                9825   0.464286  DACOSTA_ERCC3_ALLELE_XPCS_VS_TTD_DN
            1  Weinmann Adaptation To Hypoxia Up                           -0.540618  0.0760139    -1.50713  0.22774    1                8896   0.478261  WEINMANN_ADAPTATION_TO_HYPOXIA_UP
            1  Weinmann Adaptation To Hypoxia Dn                           -0.559896  0.0613646    -1.56907  0.197731   1                9061   0.451613  WEINMANN_ADAPTATION_TO_HYPOXIA_DN
            1  Tsunoda Cisplatin Resistance Up                             -0.721454  0.0150572    -1.73759  0.141912   1                9652   0.642857  TSUNODA_CISPLATIN_RESISTANCE_UP
            1  Tsunoda Cisplatin Resistance Dn                             -0.455676  0.0564533    -1.50354  0.229309   1                9611   0.348837  TSUNODA_CISPLATIN_RESISTANCE_DN
            1  Li Wilms Tumor Anaplastic Up                                 0.902637  0.000405268   1.89282  0.0220871  1                 154   0.666667  LI_WILMS_TUMOR_ANAPLASTIC_UP
            4  Li Wilms Tumor Anaplastic Up                                -0.780747  0.0333268    -1.62607  0.230529   1               10178   0.666667  LI_WILMS_TUMOR_ANAPLASTIC_UP
            5  Li Wilms Tumor Anaplastic Up                                 0.834676  0.00717846    1.74106  0.117295   1                1063   0.8       LI_WILMS_TUMOR_ANAPLASTIC_UP
           10  Li Wilms Tumor Anaplastic Up                                 0.805572  0.0186149     1.68765  0.141043   1                 766   0.666667  LI_WILMS_TUMOR_ANAPLASTIC_UP
            1  Mann Response To Amifostine Dn                               0.829019  0.00143325    1.84146  0.0319879  1                 251   0.4       MANN_RESPONSE_TO_AMIFOSTINE_DN
            5  Mann Response To Amifostine Dn                               0.750524  0.0146371     1.66728  0.158116   1                 810   0.4       MANN_RESPONSE_TO_AMIFOSTINE_DN
           10  Mann Response To Amifostine Dn                               0.758917  0.0133866     1.70269  0.130517   1                1189   0.5       MANN_RESPONSE_TO_AMIFOSTINE_DN
            1  Oxford Rala Or Ralb Targets Up                               0.840496  0.000398883   1.96048  0.0147521  1                1236   0.813953  OXFORD_RALA_OR_RALB_TARGETS_UP
            4  Oxford Rala Or Ralb Targets Up                              -0.753344  0.00858872   -1.76536  0.192508   1                9771   0.697674  OXFORD_RALA_OR_RALB_TARGETS_UP
            5  Oxford Rala Or Ralb Targets Up                               0.738079  0.0198098     1.71719  0.128531   1                 970   0.55814   OXFORD_RALA_OR_RALB_TARGETS_UP
           10  Oxford Rala Or Ralb Targets Up                               0.786529  0.00238711    1.84398  0.0674455  1                1404   0.744186  OXFORD_RALA_OR_RALB_TARGETS_UP
            1  Mohankumar Tlx1 Targets Dn                                  -0.454312  0.0437158    -1.56059  0.200146   1                8163   0.495798  MOHANKUMAR_TLX1_TARGETS_DN
            1  Johansson Brain Cancer Early Vs Late Dn                     -0.496131  0.0123332    -1.70925  0.146665   1               10372   0.263158  JOHANSSON_BRAIN_CANCER_EARLY_VS_LATE_DN
            1  Li Wilms Tumor Vs Fetal Kidney 2 Up                          0.534477  0.0483516     1.55062  0.169529   1                 829   0.428571  LI_WILMS_TUMOR_VS_FETAL_KIDNEY_2_UP
            5  Li Wilms Tumor Vs Fetal Kidney 2 Up                          0.594638  0.0180723     1.70797  0.135404   1                1013   0.52381   LI_WILMS_TUMOR_VS_FETAL_KIDNEY_2_UP
            1  Li Wilms Tumor Vs Fetal Kidney 2 Dn                         -0.620212  0.00949879   -1.84557  0.119811   1                8950   0.534884  LI_WILMS_TUMOR_VS_FETAL_KIDNEY_2_DN
            1  Inamura Lung Cancer Scc Up                                   0.631997  0.0785215     1.48682  0.215694   1                1346   0.461538  INAMURA_LUNG_CANCER_SCC_UP
            5  Inamura Lung Cancer Scc Up                                   0.716146  0.0163967     1.67173  0.15633    1                 422   0.384615  INAMURA_LUNG_CANCER_SCC_UP
           10  Inamura Lung Cancer Scc Up                                   0.669247  0.0433735     1.5599   0.249257   1                1514   0.538462  INAMURA_LUNG_CANCER_SCC_UP
            1  Inamura Lung Cancer Scc Subtypes Up                          0.688658  0.0251915     1.64935  0.103998   1                1813   0.538462  INAMURA_LUNG_CANCER_SCC_SUBTYPES_UP
            5  Inamura Lung Cancer Scc Subtypes Up                          0.762083  0.00339389    1.82299  0.0764752  1                1216   0.538462  INAMURA_LUNG_CANCER_SCC_SUBTYPES_UP
           10  Inamura Lung Cancer Scc Subtypes Up                          0.689609  0.0238806     1.64678  0.17341    1                2379   0.615385  INAMURA_LUNG_CANCER_SCC_SUBTYPES_UP
            1  Wattel Autonomous Thyroid Adenoma Dn                        -0.590229  0.0122564    -1.77807  0.134462   1                9376   0.522727  WATTEL_AUTONOMOUS_THYROID_ADENOMA_DN
            1  Sasai Resistance To Neoplastic Transfromation               -0.720255  0.00325203   -1.90595  0.116729   1                9281   0.680851  SASAI_RESISTANCE_TO_NEOPLASTIC_TRANSFROMATION
            1  Kim Myc Amplification Targets Up                             0.417771  0.020004      1.63511  0.110528   1                1578   0.343949  KIM_MYC_AMPLIFICATION_TARGETS_UP
            1  Yang Breast Cancer Esr1 Up                                  -0.747861  0.0280203    -1.66106  0.162756   1                9348   0.8       YANG_BREAST_CANCER_ESR1_UP
            1  Yang Breast Cancer Esr1 Laser Up                            -0.628462  0.0951413    -1.5128   0.224944   1                7628   0.766667  YANG_BREAST_CANCER_ESR1_LASER_UP
            1  Yang Breast Cancer Esr1 Bulk Up                             -0.659405  0.0891331    -1.51838  0.222033   1                7595   0.833333  YANG_BREAST_CANCER_ESR1_BULK_UP
            1  Yang Breast Cancer Esr1 Bulk Dn                              0.671958  0.092773      1.49965  0.209696   1                1433   0.526316  YANG_BREAST_CANCER_ESR1_BULK_DN
            5  Yang Breast Cancer Esr1 Bulk Dn                              0.717247  0.0316957     1.59243  0.215809   1                1466   0.578947  YANG_BREAST_CANCER_ESR1_BULK_DN
            1  Howlin Pubertal Mammary Gland                               -0.496972  0.0873806    -1.49987  0.230768   1                9730   0.38      HOWLIN_PUBERTAL_MAMMARY_GLAND
            1  Hatada Methylated In Lung Cancer Up                         -0.444852  0.00245349   -1.79776  0.124176   1                9005   0.382653  HATADA_METHYLATED_IN_LUNG_CANCER_UP
            1  Furukawa Dusp6 Targets Pci35 Dn                              0.761552  0.00121017    1.9956   0.0129733  1                 903   0.588235  FURUKAWA_DUSP6_TARGETS_PCI35_DN
            4  Furukawa Dusp6 Targets Pci35 Dn                             -0.669     0.0297713    -1.74785  0.192508   1               10030   0.54902   FURUKAWA_DUSP6_TARGETS_PCI35_DN
            5  Furukawa Dusp6 Targets Pci35 Dn                              0.699949  0.0118047     1.8134   0.0809656  1                 949   0.54902   FURUKAWA_DUSP6_TARGETS_PCI35_DN
           10  Furukawa Dusp6 Targets Pci35 Dn                              0.753336  0.00141357    1.9713   0.056431   1                 836   0.588235  FURUKAWA_DUSP6_TARGETS_PCI35_DN
            1  Sheth Liver Cancer Vs Txnip Loss Pam6                       -0.463506  0.0527043    -1.50831  0.22774    1                9682   0.366667  SHETH_LIVER_CANCER_VS_TXNIP_LOSS_PAM6
            1  Lien Breast Carcinoma Metaplastic                           -0.700764  0.0667614    -1.60157  0.181764   1                9642   0.7       LIEN_BREAST_CARCINOMA_METAPLASTIC
            1  Tomlins Prostate Cancer Dn                                  -0.639295  0.0371262    -1.65749  0.16369    1                8644   0.675676  TOMLINS_PROSTATE_CANCER_DN
            1  Tomlins Metastasis Dn                                       -0.580785  0.0269013    -1.63172  0.170938   1                8431   0.473684  TOMLINS_METASTASIS_DN
            1  Liang Hematopoiesis Stem Cell Number Large Vs Tiny Dn        0.463092  0.0369775     1.53739  0.178544   1                3140   0.472222  LIANG_HEMATOPOIESIS_STEM_CELL_NUMBER_LARGE_VS_TINY_DN
            1  Hui Mapk14 Targets Up                                       -0.829067  0.000603136  -1.98394  0.116729   1               10354   0.545455  HUI_MAPK14_TARGETS_UP
            1  Theodorou Mammary Tumorigenesis                             -0.606902  0.0567062    -1.55217  0.202734   1               10342   0.411765  THEODOROU_MAMMARY_TUMORIGENESIS
            1  Pujana Breast Cancer Lit Int Network                         0.529335  0.0169972     1.74801  0.0583926  1                1372   0.395604  PUJANA_BREAST_CANCER_LIT_INT_NETWORK
            5  Pujana Breast Cancer Lit Int Network                         0.551277  0.00860516    1.81641  0.0793115  1                1571   0.406593  PUJANA_BREAST_CANCER_LIT_INT_NETWORK
           10  Pujana Breast Cancer Lit Int Network                         0.496316  0.0390672     1.63732  0.177303   1                1425   0.373626  PUJANA_BREAST_CANCER_LIT_INT_NETWORK
            1  Pujana Xprss Int Network                                     0.729681  0.000409584   2.01727  0.0127991  1                1920   0.679739  PUJANA_XPRSS_INT_NETWORK
            5  Pujana Xprss Int Network                                     0.757037  0.000197824   2.09475  0.0336414  0.526409         1572   0.686275  PUJANA_XPRSS_INT_NETWORK
           10  Pujana Xprss Int Network                                     0.685018  0.00406256    1.89542  0.0620676  1                1419   0.568627  PUJANA_XPRSS_INT_NETWORK
            1  Pujana Brca Centered Network                                 0.73481   0.000611496   1.96465  0.0147521  1                1890   0.696429  PUJANA_BRCA_CENTERED_NETWORK
            5  Pujana Brca Centered Network                                 0.742924  0.000394711   1.99689  0.0354642  1                1572   0.669643  PUJANA_BRCA_CENTERED_NETWORK
           10  Pujana Brca Centered Network                                 0.694478  0.00466059    1.85677  0.0652612  1                1832   0.625     PUJANA_BRCA_CENTERED_NETWORK
            1  Schlesinger Methylated De Novo In Cancer                    -0.604359  0.00364225   -1.85733  0.119811   1                9355   0.534884  SCHLESINGER_METHYLATED_DE_NOVO_IN_CANCER
            1  Goering Blood Hdl Cholesterol Qtl Trans                     -0.62164   0.0759544    -1.47865  0.242484   1                8468   0.7       GOERING_BLOOD_HDL_CHOLESTEROL_QTL_TRANS
            1  Khetchoumian Trim24 Targets Up                              -0.656146  0.0132424    -1.80767  0.124176   1               10384   0.409091  KHETCHOUMIAN_TRIM24_TARGETS_UP
            1  Alcala Apoptosis                                             0.463621  0.0144289     1.6698   0.0948624  1                2245   0.417722  ALCALA_APOPTOSIS
            1  Cowling Mycn Targets                                        -0.615857  0.0664371    -1.56759  0.197985   1                7264   0.769231  COWLING_MYCN_TARGETS
            1  Kauffmann Melanoma Relapse Up                                0.793572  0.0010101     1.934    0.0171414  1                 660   0.630435  KAUFFMANN_MELANOMA_RELAPSE_UP
            5  Kauffmann Melanoma Relapse Up                                0.77036   0.00358066    1.87568  0.0582517  1                1271   0.695652  KAUFFMANN_MELANOMA_RELAPSE_UP
           10  Kauffmann Melanoma Relapse Up                                0.758354  0.00259585    1.85617  0.0652612  1                1160   0.673913  KAUFFMANN_MELANOMA_RELAPSE_UP
            1  Williams Esr1 Targets Up                                     0.567144  0.0584206     1.54454  0.173732   1                 481   0.380952  WILLIAMS_ESR1_TARGETS_UP
            1  Gross Elk3 Targets Up                                       -0.590946  0.00796894   -1.77493  0.134462   1                9248   0.5       GROSS_ELK3_TARGETS_UP
            1  Gross Hypoxia Via Hif1A Up                                   0.474172  0.0133573     1.69862  0.0802671  1                1951   0.405405  GROSS_HYPOXIA_VIA_HIF1A_UP
            1  Gross Hypoxia Via Hif1A Dn                                  -0.432367  0.0230507    -1.62373  0.170938   1                9307   0.364583  GROSS_HYPOXIA_VIA_HIF1A_DN
            1  Ingram Shh Targets Up                                       -0.565292  0.00101585   -1.96685  0.116729   1                9146   0.483516  INGRAM_SHH_TARGETS_UP
            1  Ingram Shh Targets Dn                                       -0.448908  0.0686747    -1.48275  0.241988   1                9348   0.395833  INGRAM_SHH_TARGETS_DN
            1  Mantovani Nfkb Targets Up                                   -0.511381  0.0578381    -1.5366   0.211252   1                9534   0.40625   MANTOVANI_NFKB_TARGETS_UP
            1  Ferreira Ewings Sarcoma Unstable Vs Stable Up                0.693581  0.000608149   2.06606  0.0127991  1                1078   0.511628  FERREIRA_EWINGS_SARCOMA_UNSTABLE_VS_STABLE_UP
            4  Ferreira Ewings Sarcoma Unstable Vs Stable Up               -0.547951  0.0661837    -1.63006  0.230529   1                9591   0.48062   FERREIRA_EWINGS_SARCOMA_UNSTABLE_VS_STABLE_UP
            5  Ferreira Ewings Sarcoma Unstable Vs Stable Up                0.674991  0.00117786    2.01466  0.0336414  1                1186   0.511628  FERREIRA_EWINGS_SARCOMA_UNSTABLE_VS_STABLE_UP
           10  Ferreira Ewings Sarcoma Unstable Vs Stable Up                0.610659  0.0148833     1.81968  0.0734069  1                1189   0.449612  FERREIRA_EWINGS_SARCOMA_UNSTABLE_VS_STABLE_UP
            1  Ferreira Ewings Sarcoma Unstable Vs Stable Dn               -0.445214  0.0251328    -1.59879  0.181764   1                7940   0.478873  FERREIRA_EWINGS_SARCOMA_UNSTABLE_VS_STABLE_DN
            1  Nakayama Fra2 Targets                                       -0.44935   0.0568272    -1.51175  0.225344   1                8404   0.461538  NAKAYAMA_FRA2_TARGETS
            1  Lopez Mbd Targets Imprinted And X Linked                     0.599731  0.0494949     1.57041  0.156915   1                1467   0.533333  LOPEZ_MBD_TARGETS_IMPRINTED_AND_X_LINKED
           10  Lopez Mbd Targets Imprinted And X Linked                     0.72755   0.000998602   1.89655  0.0620676  1                1606   0.666667  LOPEZ_MBD_TARGETS_IMPRINTED_AND_X_LINKED
            1  Lockwood Amplified In Lung Cancer                            0.377933  0.0927363     1.4317   0.249302   1                2765   0.403141  LOCKWOOD_AMPLIFIED_IN_LUNG_CANCER
            1  Rahman Tp53 Targets Phosphorylated                           0.61165   0.0577617     1.55417  0.166867   1                2146   0.52381   RAHMAN_TP53_TARGETS_PHOSPHORYLATED
            4  Rahman Tp53 Targets Phosphorylated                          -0.69929   0.00513732   -1.78083  0.192508   1                8218   0.809524  RAHMAN_TP53_TARGETS_PHOSPHORYLATED
            1  Scian Cell Cycle Targets Of Tp53 And Tp73 Dn                 0.868145  0.00161031    1.76955  0.0513559  1                 612   0.789474  SCIAN_CELL_CYCLE_TARGETS_OF_TP53_AND_TP73_DN
            4  Scian Cell Cycle Targets Of Tp53 And Tp73 Dn                -0.841167  0.00589971   -1.73283  0.192508   1               10154   0.736842  SCIAN_CELL_CYCLE_TARGETS_OF_TP53_AND_TP73_DN
            5  Scian Cell Cycle Targets Of Tp53 And Tp73 Dn                 0.803934  0.0211672     1.63787  0.184004   1                1031   0.789474  SCIAN_CELL_CYCLE_TARGETS_OF_TP53_AND_TP73_DN
           10  Scian Cell Cycle Targets Of Tp53 And Tp73 Dn                 0.845204  0.00380609    1.72443  0.117315   1                 921   0.789474  SCIAN_CELL_CYCLE_TARGETS_OF_TP53_AND_TP73_DN
            1  Guenther Growth Spherical Vs Adherent Dn                    -0.625666  0.018845     -1.71321  0.146665   1                8204   0.56      GUENTHER_GROWTH_SPHERICAL_VS_ADHERENT_DN
            1  Caffarel Response To Thc Dn                                  0.68506   0.0140449     1.72884  0.0663823  1                2550   0.833333  CAFFAREL_RESPONSE_TO_THC_DN
            5  Caffarel Response To Thc Dn                                  0.693114  0.00998004    1.74425  0.117295   1                1839   0.583333  CAFFAREL_RESPONSE_TO_THC_DN
            1  Caffarel Response To Thc 8Hr 5 Dn                            0.646302  0.0648991     1.50181  0.208475   1                2357   0.7       CAFFAREL_RESPONSE_TO_THC_8HR_5_DN
            1  Caffarel Response To Thc 24Hr 3 Dn                           0.648327  0.0606485     1.51925  0.192032   1                  96   0.2       CAFFAREL_RESPONSE_TO_THC_24HR_3_DN
            1  Caffarel Response To Thc 24Hr 5 Up                           0.545942  0.027739      1.62929  0.112797   1                2198   0.533333  CAFFAREL_RESPONSE_TO_THC_24HR_5_UP
            4  Caffarel Response To Thc 24Hr 5 Up                          -0.56781   0.0143644    -1.68597  0.201047   1                9422   0.466667  CAFFAREL_RESPONSE_TO_THC_24HR_5_UP
            1  Caffarel Response To Thc 24Hr 5 Dn                           0.584933  0.00659341    1.80071  0.0417664  1                2060   0.489796  CAFFAREL_RESPONSE_TO_THC_24HR_5_DN
           10  Caffarel Response To Thc 24Hr 5 Dn                           0.506233  0.0528818     1.5672   0.243334   1                1837   0.408163  CAFFAREL_RESPONSE_TO_THC_24HR_5_DN
            1  Bertucci Invasive Carcinoma Ductal Vs Lobular Up             0.5393    0.0424482     1.54549  0.173284   1                3611   0.684211  BERTUCCI_INVASIVE_CARCINOMA_DUCTAL_VS_LOBULAR_UP
            1  Bertucci Invasive Carcinoma Ductal Vs Lobular Dn            -0.756075  0.0106747    -1.78557  0.131106   1                9090   0.714286  BERTUCCI_INVASIVE_CARCINOMA_DUCTAL_VS_LOBULAR_DN
            1  Wotton Runx Targets Up                                      -0.531184  0.0756217    -1.48103  0.242366   1               10684   0.3125    WOTTON_RUNX_TARGETS_UP
            1  Wotton Runx Targets Dn                                      -0.572544  0.0203959    -1.66978  0.162484   1                8761   0.458333  WOTTON_RUNX_TARGETS_DN
            1  Garcia Targets Of Fli1 And Dax1 Up                          -0.521112  0.0148515    -1.68573  0.154958   1                7592   0.622222  GARCIA_TARGETS_OF_FLI1_AND_DAX1_UP
            1  Garcia Targets Of Fli1 And Dax1 Dn                           0.614429  0.000606551   1.97253  0.0140167  1                1854   0.504     GARCIA_TARGETS_OF_FLI1_AND_DAX1_DN
            4  Garcia Targets Of Fli1 And Dax1 Dn                          -0.509873  0.0347346    -1.63528  0.229428   1                9333   0.416     GARCIA_TARGETS_OF_FLI1_AND_DAX1_DN
            5  Garcia Targets Of Fli1 And Dax1 Dn                           0.542069  0.0149224     1.7569   0.109573   1                2091   0.472     GARCIA_TARGETS_OF_FLI1_AND_DAX1_DN
           10  Garcia Targets Of Fli1 And Dax1 Dn                           0.552897  0.00842021    1.78125  0.0904405  1                2547   0.496     GARCIA_TARGETS_OF_FLI1_AND_DAX1_DN
            1  Rickman Tumor Differentiated Well Vs Poorly Up               0.442576  0.0458308     1.58055  0.151135   1                2219   0.392157  RICKMAN_TUMOR_DIFFERENTIATED_WELL_VS_POORLY_UP
            5  Rickman Tumor Differentiated Well Vs Poorly Up               0.523735  0.00221551    1.85823  0.0620503  1                2039   0.431373  RICKMAN_TUMOR_DIFFERENTIATED_WELL_VS_POORLY_UP
            1  Rickman Tumor Differentiated Moderately Vs Poorly Dn         0.749255  0.00798722    1.71385  0.0743741  1                1919   0.692308  RICKMAN_TUMOR_DIFFERENTIATED_MODERATELY_VS_POORLY_DN
            5  Rickman Tumor Differentiated Moderately Vs Poorly Dn         0.756587  0.00895879    1.7146   0.128776   1                1383   0.538462  RICKMAN_TUMOR_DIFFERENTIATED_MODERATELY_VS_POORLY_DN
            1  Fridman Immortalization Dn                                  -0.582144  0.0447187    -1.62177  0.170938   1                9010   0.53125   FRIDMAN_IMMORTALIZATION_DN
            1  Schaeffer Prostate Development 12Hr Up                      -0.407874  0.0627403    -1.47311  0.246295   1                9526   0.333333  SCHAEFFER_PROSTATE_DEVELOPMENT_12HR_UP
            1  Schaeffer Prostate Development 48Hr Dn                      -0.451338  0.0409968    -1.62157  0.170938   1                8091   0.505455  SCHAEFFER_PROSTATE_DEVELOPMENT_48HR_DN
            1  Schaeffer Prostate Development And Cancer Box5 Up           -0.693274  0.0882646    -1.47333  0.246295   1               10686   0.4       SCHAEFFER_PROSTATE_DEVELOPMENT_AND_CANCER_BOX5_UP
            1  Koyama Sema3B Targets Up                                    -0.352407  0.0217173    -1.51284  0.224944   1                8457   0.381643  KOYAMA_SEMA3B_TARGETS_UP
            1  Rickman Head And Neck Cancer F                              -0.639446  0.0202648    -1.66546  0.162484   1               10236   0.48      RICKMAN_HEAD_AND_NECK_CANCER_F
            1  Wu Cell Migration                                           -0.436961  0.0882711    -1.4953   0.232098   1                7874   0.532051  WU_CELL_MIGRATION
            1  Colin Pilocytic Astrocytoma Vs Glioblastoma Up              -0.543165  0.0479115    -1.55063  0.202734   1                9959   0.391304  COLIN_PILOCYTIC_ASTROCYTOMA_VS_GLIOBLASTOMA_UP
            1  Amit Delayed Early Genes                                    -0.754088  0.00982949   -1.76708  0.138176   1                9677   0.625     AMIT_DELAYED_EARLY_GENES
            1  Den Interact With Lca5                                       0.614484  0.0346752     1.63248  0.111331   1                3249   0.8       DEN_INTERACT_WITH_LCA5
            4  Den Interact With Lca5                                      -0.716636  0.000593237  -1.90489  0.155097   1                8795   0.84      DEN_INTERACT_WITH_LCA5
            1  Benporath Es 1                                               0.581536  0.00121877    1.97094  0.0140888  1                1354   0.425606  BENPORATH_ES_1
            5  Benporath Es 1                                               0.5955    0.000396275   2.02673  0.0336414  1                1947   0.512111  BENPORATH_ES_1
           10  Benporath Es 1                                               0.566817  0.00140252    1.92964  0.056431   1                2170   0.515571  BENPORATH_ES_1
            1  Benporath Es 2                                               0.685916  0.035287      1.58077  0.151135   1                1822   0.6       BENPORATH_ES_2
            4  Benporath Es 2                                              -0.707951  0.0215983    -1.64467  0.224407   1                8903   0.666667  BENPORATH_ES_2
            5  Benporath Es 2                                               0.751878  0.00530452    1.7498   0.115744   1                1592   0.6       BENPORATH_ES_2
            1  Benporath Prc2 Targets                                      -0.516905  0.00320834   -1.83841  0.119811   1                8862   0.461957  BENPORATH_PRC2_TARGETS
            1  Benporath Myc Targets With Ebox                              0.395618  0.0226498     1.57736  0.153224   1                1928   0.368132  BENPORATH_MYC_TARGETS_WITH_EBOX
            1  Benporath Proliferation                                      0.802618  0.000202511   1.94012  0.0167144  0.538882         1550   0.803279  BENPORATH_PROLIFERATION
            4  Benporath Proliferation                                     -0.774421  0.00138176   -1.87385  0.162281   1                9659   0.737705  BENPORATH_PROLIFERATION
            5  Benporath Proliferation                                      0.720042  0.0165141     1.73935  0.117295   1                2536   0.795082  BENPORATH_PROLIFERATION
           10  Benporath Proliferation                                      0.764295  0.002414      1.85679  0.0652612  1                1282   0.696721  BENPORATH_PROLIFERATION
            1  Petretto Cardiac Hypertrophy                                -0.61881   0.0402534    -1.67429  0.161802   1                8204   0.787879  PETRETTO_CARDIAC_HYPERTROPHY
            1  Stark Hyppocampus 22Q11 Deletion Up                         -0.512729  0.0313874    -1.61148  0.176986   1                8702   0.5       STARK_HYPPOCAMPUS_22Q11_DELETION_UP
            1  Benporath Es Core Nine Correlated                            0.648961  0.0281951     1.68987  0.0843017  1                2728   0.658824  BENPORATH_ES_CORE_NINE_CORRELATED
            5  Benporath Es Core Nine Correlated                            0.677465  0.00851991    1.76527  0.102916   1                2203   0.6       BENPORATH_ES_CORE_NINE_CORRELATED
           10  Benporath Es Core Nine Correlated                            0.630978  0.0436394     1.63939  0.175537   1                3222   0.705882  BENPORATH_ES_CORE_NINE_CORRELATED
            1  Amit Egf Response 40 Hela                                   -0.665215  0.013089     -1.85199  0.119811   1                9061   0.583333  AMIT_EGF_RESPONSE_40_HELA
            1  Amit Egf Response 60 Hela                                   -0.50336   0.0478392    -1.59765  0.181764   1                9081   0.465116  AMIT_EGF_RESPONSE_60_HELA
            1  Amit Egf Response 120 Hela                                  -0.472236  0.0548878    -1.56305  0.199134   1                7497   0.534483  AMIT_EGF_RESPONSE_120_HELA
            1  Amit Egf Response 40 Mcf10A                                 -0.544273  0.0843687    -1.50827  0.22774    1                9642   0.526316  AMIT_EGF_RESPONSE_40_MCF10A
            1  Amit Egf Response 60 Mcf10A                                 -0.605341  0.022343     -1.71394  0.146665   1                9051   0.6       AMIT_EGF_RESPONSE_60_MCF10A
            1  Amit Serum Response 40 Mcf10A                               -0.624443  0.0539078    -1.62272  0.170938   1                9106   0.642857  AMIT_SERUM_RESPONSE_40_MCF10A
            1  Amit Serum Response 60 Mcf10A                               -0.602361  0.0137395    -1.81455  0.121024   1                7488   0.759259  AMIT_SERUM_RESPONSE_60_MCF10A
            1  Amit Serum Response 240 Mcf10A                              -0.490938  0.0383282    -1.60085  0.181764   1                8950   0.413043  AMIT_SERUM_RESPONSE_240_MCF10A
            1  Georges Cell Cycle Mir192 Targets                            0.47879   0.0967478     1.47858  0.219334   1                2036   0.462963  GEORGES_CELL_CYCLE_MIR192_TARGETS
            5  Georges Cell Cycle Mir192 Targets                            0.558757  0.0154062     1.71557  0.128776   1                1406   0.407407  GEORGES_CELL_CYCLE_MIR192_TARGETS
            1  Sakai Tumor Infiltrating Monocytes Dn                        0.527092  0.00756671    1.76287  0.0535254  1                2625   0.536232  SAKAI_TUMOR_INFILTRATING_MONOCYTES_DN
            1  Sakai Chronic Hepatitis Vs Liver Cancer Dn                  -0.588477  0.00859484   -1.77602  0.134462   1                8670   0.7       SAKAI_CHRONIC_HEPATITIS_VS_LIVER_CANCER_DN
            1  Jeon Smad6 Targets Up                                       -0.675789  0.033086     -1.64323  0.170938   1                8759   0.65      JEON_SMAD6_TARGETS_UP
            1  Jeon Smad6 Targets Dn                                        0.799058  0.00338578    1.82233  0.0363179  1                 352   0.583333  JEON_SMAD6_TARGETS_DN
           10  Jeon Smad6 Targets Dn                                        0.788261  0.00495835    1.80254  0.0795671  1                1747   0.833333  JEON_SMAD6_TARGETS_DN
            1  Sung Metastasis Stroma Up                                   -0.437796  0.0777912    -1.51624  0.222374   1                8744   0.43      SUNG_METASTASIS_STROMA_UP
            1  Sung Metastasis Stroma Dn                                    0.598808  0.0164234     1.75856  0.0545451  1                1224   0.465116  SUNG_METASTASIS_STROMA_DN
            5  Sung Metastasis Stroma Dn                                    0.591012  0.015121      1.73161  0.123266   1                1056   0.465116  SUNG_METASTASIS_STROMA_DN
           10  Sung Metastasis Stroma Dn                                    0.596261  0.0149671     1.75399  0.107638   1                1537   0.465116  SUNG_METASTASIS_STROMA_DN
            1  Shin B Cell Lymphoma Cluster 3                              -0.546083  0.0930462    -1.48373  0.241988   1                9275   0.375     SHIN_B_CELL_LYMPHOMA_CLUSTER_3
            5  Zhang Response To Ikk Inhibitor And Tnf Dn                   0.430269  0.0185449     1.60117  0.209048   1                1848   0.329114  ZHANG_RESPONSE_TO_IKK_INHIBITOR_AND_TNF_DN
            1  Mori Pre Bi Lymphocyte Up                                    0.736515  0.0014093     1.93607  0.0171414  1                1316   0.651515  MORI_PRE_BI_LYMPHOCYTE_UP
            4  Mori Pre Bi Lymphocyte Up                                   -0.693421  0.00903023   -1.82228  0.187832   1               10229   0.515152  MORI_PRE_BI_LYMPHOCYTE_UP
            5  Mori Pre Bi Lymphocyte Up                                    0.640083  0.0424765     1.66947  0.156755   1                1295   0.515152  MORI_PRE_BI_LYMPHOCYTE_UP
           10  Mori Pre Bi Lymphocyte Up                                    0.685208  0.0150973     1.8044   0.0789389  1                1255   0.575758  MORI_PRE_BI_LYMPHOCYTE_UP
            1  Mori Large Pre Bii Lymphocyte Up                             0.799753  0.000201816   1.9577   0.0149997  0.537033         1111   0.717949  MORI_LARGE_PRE_BII_LYMPHOCYTE_UP
            4  Mori Large Pre Bii Lymphocyte Up                            -0.724683  0.0167488    -1.77225  0.192508   1               10185   0.564103  MORI_LARGE_PRE_BII_LYMPHOCYTE_UP
            5  Mori Large Pre Bii Lymphocyte Up                             0.741553  0.00897666    1.81081  0.08228    1                1295   0.666667  MORI_LARGE_PRE_BII_LYMPHOCYTE_UP
           10  Mori Large Pre Bii Lymphocyte Up                             0.755948  0.00341091    1.8478   0.0662351  1                1417   0.692308  MORI_LARGE_PRE_BII_LYMPHOCYTE_UP
            1  Mori Immature B Lymphocyte Dn                                0.797096  0.00020141    1.96145  0.0147521  0.535952          894   0.701299  MORI_IMMATURE_B_LYMPHOCYTE_DN
            4  Mori Immature B Lymphocyte Dn                               -0.691896  0.0369352    -1.70234  0.192508   1                9392   0.649351  MORI_IMMATURE_B_LYMPHOCYTE_DN
            5  Mori Immature B Lymphocyte Dn                                0.754564  0.00539353    1.84905  0.0650475  1                 935   0.623377  MORI_IMMATURE_B_LYMPHOCYTE_DN
           10  Mori Immature B Lymphocyte Dn                                0.751815  0.00560224    1.85144  0.0652612  1                1189   0.649351  MORI_IMMATURE_B_LYMPHOCYTE_DN
            1  Mori Mature B Lymphocyte Dn                                  0.597691  0.00690496    1.80328  0.0411858  1                1205   0.42623   MORI_MATURE_B_LYMPHOCYTE_DN
            4  Mori Mature B Lymphocyte Dn                                 -0.632899  0.00256158   -1.91807  0.155097   1               10104   0.508197  MORI_MATURE_B_LYMPHOCYTE_DN
            1  Mori Emu Myc Lymphoma By Onset Time Up                       0.565347  0.00301205    1.89091  0.0220871  1                1216   0.395833  MORI_EMU_MYC_LYMPHOMA_BY_ONSET_TIME_UP
            5  Mori Emu Myc Lymphoma By Onset Time Up                       0.48639   0.0343245     1.62348  0.191336   1                1988   0.416667  MORI_EMU_MYC_LYMPHOMA_BY_ONSET_TIME_UP
           10  Mori Emu Myc Lymphoma By Onset Time Up                       0.545313  0.00522403    1.82647  0.0688673  1                1033   0.354167  MORI_EMU_MYC_LYMPHOMA_BY_ONSET_TIME_UP
            1  Collis Prkdc Regulators                                      0.542401  0.0708518     1.47859  0.219334   1                2469   0.5       COLLIS_PRKDC_REGULATORS
            1  Kauffmann Dna Repair Genes                                   0.496275  0.00909642    1.78282  0.0476102  1                2173   0.438503  KAUFFMANN_DNA_REPAIR_GENES
            5  Kauffmann Dna Repair Genes                                   0.526083  0.00200361    1.88472  0.0572692  1                2481   0.486631  KAUFFMANN_DNA_REPAIR_GENES
           10  Kauffmann Dna Repair Genes                                   0.523384  0.00302176    1.87627  0.0627622  1                2322   0.459893  KAUFFMANN_DNA_REPAIR_GENES
            1  Kauffmann Dna Replication Genes                              0.598834  0.00244848    1.9546   0.015348   1                1994   0.516393  KAUFFMANN_DNA_REPLICATION_GENES
            5  Kauffmann Dna Replication Genes                              0.589204  0.00178359    1.92652  0.0475578  1                1905   0.459016  KAUFFMANN_DNA_REPLICATION_GENES
           10  Kauffmann Dna Replication Genes                              0.600849  0.00201857    1.96832  0.056431   1                1934   0.516393  KAUFFMANN_DNA_REPLICATION_GENES
           10  Nikolsky Breast Cancer 6P24 P22 Amplicon                     0.67104   0.0567288     1.60269  0.207217   1                1530   0.444444  NIKOLSKY_BREAST_CANCER_6P24_P22_AMPLICON
            1  Nikolsky Breast Cancer 16Q24 Amplicon                        0.720098  0.0478981     1.63357  0.110867   1                2518   0.75      NIKOLSKY_BREAST_CANCER_16Q24_AMPLICON
           10  Nikolsky Breast Cancer 16Q24 Amplicon                        0.812968  0.0048048     1.85218  0.0652612  1                1142   0.78125   NIKOLSKY_BREAST_CANCER_16Q24_AMPLICON
            1  Stegmeier Pre-Mitotic Cell Cycle Regulators                  0.634473  0.0369622     1.56227  0.160717   1                 319   0.3       STEGMEIER_PRE-MITOTIC_CELL_CYCLE_REGULATORS
            1  Ji Metastasis Repressed By Stk11                            -0.576679  0.0467083    -1.60078  0.181764   1                7747   0.636364  JI_METASTASIS_REPRESSED_BY_STK11
            1  Onder Cdh1 Targets 2 Up                                     -0.635721  0.00565885   -1.92809  0.116729   1                8791   0.633333  ONDER_CDH1_TARGETS_2_UP
            1  Onder Cdh1 Targets 1 Up                                     -0.373342  0.0146074    -1.55354  0.202734   1                7960   0.403509  ONDER_CDH1_TARGETS_1_UP
            1  Onder Cdh1 Signaling Via Ctnnb1                             -0.601726  0.0129345    -1.81691  0.119811   1                8921   0.56      ONDER_CDH1_SIGNALING_VIA_CTNNB1
            1  Cervera Sdhb Targets 2                                      -0.463475  0.0317009    -1.61146  0.176986   1                8766   0.43038   CERVERA_SDHB_TARGETS_2
            1  Rozanov Mmp14 Targets Subset                                -0.669211  0.106082     -1.52981  0.214421   1                8759   0.642857  ROZANOV_MMP14_TARGETS_SUBSET
            1  Rozanov Mmp14 Targets Up                                    -0.424537  0.0319886    -1.61333  0.176986   1                8625   0.436275  ROZANOV_MMP14_TARGETS_UP
            5  Landemaine Lung Metastasis                                   0.690487  0.0218111     1.64964  0.177575   1                1546   0.571429  LANDEMAINE_LUNG_METASTASIS
            1  Cervera Sdhb Targets 1 Up                                   -0.520881  0.0186084    -1.71078  0.146665   1                8681   0.465753  CERVERA_SDHB_TARGETS_1_UP
            1  Ross Acute Myeloid Leukemia Cbf                             -0.465853  0.0382127    -1.58201  0.190354   1                8694   0.411765  ROSS_ACUTE_MYELOID_LEUKEMIA_CBF
            1  Kenny Ctnnb1 Targets Up                                      0.474824  0.00965018    1.68354  0.0868165  1                1264   0.297872  KENNY_CTNNB1_TARGETS_UP
            1  Ross Aml With Pml Rara Fusion                               -0.433865  0.029707     -1.57365  0.195122   1                7472   0.571429  ROSS_AML_WITH_PML_RARA_FUSION
            1  Park Hsc Vs Multipotent Progenitors Dn                       0.548159  0.0882591     1.45473  0.235824   1                1706   0.428571  PARK_HSC_VS_MULTIPOTENT_PROGENITORS_DN
            1  Shipp Dlbcl Cured Vs Fatal Dn                                0.444761  0.0501107     1.49357  0.213453   1                3111   0.555556  SHIPP_DLBCL_CURED_VS_FATAL_DN
           10  Shipp Dlbcl Cured Vs Fatal Dn                                0.470461  0.0265575     1.57858  0.231959   1                2839   0.527778  SHIPP_DLBCL_CURED_VS_FATAL_DN
            4  Houstis Ros                                                 -0.620256  0.00562023   -1.83252  0.187832   1                9940   0.533333  HOUSTIS_ROS
            1  Basso B Lymphocyte Network                                   0.569544  0.000604351   1.9498   0.015456   1                1517   0.476562  BASSO_B_LYMPHOCYTE_NETWORK
            4  Basso B Lymphocyte Network                                  -0.486546  0.0250939    -1.66427  0.214503   1                9431   0.492188  BASSO_B_LYMPHOCYTE_NETWORK
           10  Basso B Lymphocyte Network                                   0.469035  0.0398798     1.60929  0.203759   1                2193   0.398438  BASSO_B_LYMPHOCYTE_NETWORK
            1  Iizuka Liver Cancer Progression L1 G1 Dn                     0.596274  0.0610442     1.50845  0.201936   1                 967   0.363636  IIZUKA_LIVER_CANCER_PROGRESSION_L1_G1_DN
            1  Iizuka Liver Cancer Progression G1 G2 Dn                     0.5189    0.0840234     1.47626  0.22054    1                1935   0.44      IIZUKA_LIVER_CANCER_PROGRESSION_G1_G2_DN
            5  Iizuka Liver Cancer Progression G1 G2 Dn                     0.578115  0.0346575     1.62508  0.19014    1                2693   0.56      IIZUKA_LIVER_CANCER_PROGRESSION_G1_G2_DN
            1  Shepard Bmyb Morpholino Dn                                   0.579235  0.00101112    2.02893  0.0127991  1                 879   0.363636  SHEPARD_BMYB_MORPHOLINO_DN
            5  Shepard Bmyb Morpholino Dn                                   0.526956  0.0105096     1.82601  0.074948   1                2212   0.5       SHEPARD_BMYB_MORPHOLINO_DN
           10  Shepard Bmyb Morpholino Dn                                   0.553246  0.00242915    1.92562  0.056431   1                1201   0.371212  SHEPARD_BMYB_MORPHOLINO_DN
            1  Smith Tert Targets Up                                        0.372212  0.0518324     1.46601  0.227776   1                2126   0.385827  SMITH_TERT_TARGETS_UP
            1  Menssen Myc Targets                                          0.647656  0.00965406    1.77193  0.0504442  1                2436   0.686275  MENSSEN_MYC_TARGETS
            4  Menssen Myc Targets                                         -0.62442   0.0179593    -1.72685  0.192508   1                8218   0.784314  MENSSEN_MYC_TARGETS
           10  Menssen Myc Targets                                          0.5766    0.0585452     1.59239  0.215799   1                3089   0.647059  MENSSEN_MYC_TARGETS
            1  Kannan Tp53 Targets Up                                      -0.523611  0.00656545   -1.79347  0.127603   1                9754   0.42      KANNAN_TP53_TARGETS_UP
            1  Shepard Bmyb Targets                                         0.751241  0.00162173    2.00551  0.0127991  1                 529   0.543478  SHEPARD_BMYB_TARGETS
            4  Shepard Bmyb Targets                                        -0.632835  0.0500591    -1.67621  0.207177   1               10271   0.521739  SHEPARD_BMYB_TARGETS
            5  Shepard Bmyb Targets                                         0.681278  0.0191284     1.8004   0.089859   1                1010   0.543478  SHEPARD_BMYB_TARGETS
           10  Shepard Bmyb Targets                                         0.684005  0.0166497     1.81485  0.0752823  1                 751   0.5       SHEPARD_BMYB_TARGETS
            1  Sana Tnf Signaling Dn                                       -0.516161  0.0248574    -1.70683  0.146665   1                8619   0.544304  SANA_TNF_SIGNALING_DN
            1  Tarte Plasma Cell Vs Plasmablast Up                         -0.411314  0.0206771    -1.62141  0.170938   1                9196   0.296     TARTE_PLASMA_CELL_VS_PLASMABLAST_UP
            1  Peng Leucine Deprivation Dn                                  0.599936  0.00461199    1.86001  0.0279962  1                2151   0.573864  PENG_LEUCINE_DEPRIVATION_DN
            4  Peng Leucine Deprivation Dn                                 -0.630884  0.000596184  -1.94941  0.155097   1                9133   0.573864  PENG_LEUCINE_DEPRIVATION_DN
           10  Peng Leucine Deprivation Dn                                  0.524823  0.0467271     1.63504  0.178312   1                3099   0.590909  PENG_LEUCINE_DEPRIVATION_DN
            1  Shepard Crush And Burn Mutant Dn                             0.57622   0.00687424    1.90197  0.0208664  1                 660   0.356522  SHEPARD_CRUSH_AND_BURN_MUTANT_DN
            5  Shepard Crush And Burn Mutant Dn                             0.589156  0.00353843    1.93712  0.0454213  1                1010   0.408696  SHEPARD_CRUSH_AND_BURN_MUTANT_DN
           10  Shepard Crush And Burn Mutant Dn                             0.571781  0.00766747    1.87527  0.0627622  1                1655   0.434783  SHEPARD_CRUSH_AND_BURN_MUTANT_DN
            1  Kim Germinal Center T Helper Dn                             -0.596731  0.0506178    -1.55951  0.201037   1                7143   0.777778  KIM_GERMINAL_CENTER_T_HELPER_DN
            1  Yagi Aml Relapse Prognosis                                   0.486788  0.0636563     1.48878  0.215694   1                1944   0.392857  YAGI_AML_RELAPSE_PROGNOSIS
            4  Yagi Aml Relapse Prognosis                                  -0.567348  0.00811239   -1.74497  0.192508   1               10595   0.321429  YAGI_AML_RELAPSE_PROGNOSIS
            1  Shipp Dlbcl Cured Vs Fatal Up                               -0.565487  0.0718799    -1.47477  0.245513   1                8747   0.428571  SHIPP_DLBCL_CURED_VS_FATAL_UP
            1  Chen Lung Cancer Survival                                    0.587136  0.0784901     1.51285  0.19807    1                1134   0.35      CHEN_LUNG_CANCER_SURVIVAL
            1  Manalo Hypoxia Dn                                            0.727784  0.000403551   2.07722  0.0127991  1                1930   0.677686  MANALO_HYPOXIA_DN
            4  Manalo Hypoxia Dn                                           -0.582592  0.0431612    -1.66296  0.214503   1                8938   0.553719  MANALO_HYPOXIA_DN
            5  Manalo Hypoxia Dn                                            0.611744  0.0234499     1.74832  0.115744   1                2691   0.603306  MANALO_HYPOXIA_DN
           10  Manalo Hypoxia Dn                                            0.678106  0.0016        1.95289  0.056431   1                2078   0.615702  MANALO_HYPOXIA_DN
            1  Bystroem Correlated With Il5 Up                              0.423307  0.0542789     1.46898  0.22633    1                2234   0.410256  BYSTROEM_CORRELATED_WITH_IL5_UP
            4  Zucchi Metastasis Up                                        -0.497206  0.0123239    -1.68016  0.204016   1                9711   0.388889  ZUCCHI_METASTASIS_UP
            1  Vantveer Breast Cancer Metastasis Dn                         0.711961  0.00426309    1.83938  0.0319879  1                1732   0.663158  VANTVEER_BREAST_CANCER_METASTASIS_DN
            4  Vantveer Breast Cancer Metastasis Dn                        -0.6675    0.0265861    -1.70238  0.192508   1                8994   0.610526  VANTVEER_BREAST_CANCER_METASTASIS_DN
            5  Vantveer Breast Cancer Metastasis Dn                         0.653036  0.0333202     1.68086  0.151662   1                2088   0.673684  VANTVEER_BREAST_CANCER_METASTASIS_DN
           10  Vantveer Breast Cancer Metastasis Dn                         0.749883  0.00122624    1.90937  0.056431   1                1546   0.684211  VANTVEER_BREAST_CANCER_METASTASIS_DN
            1  Brocke Apoptosis Reversed By Il6                            -0.435079  0.0437941    -1.57875  0.193455   1                8638   0.377049  BROCKE_APOPTOSIS_REVERSED_BY_IL6
            1  Pomeroy Medulloblastoma Desmoplasic Vs Classic Dn           -0.491292  0.0855369    -1.50512  0.228916   1                9474   0.347826  POMEROY_MEDULLOBLASTOMA_DESMOPLASIC_VS_CLASSIC_DN
            1  Frasor Response To Estradiol Dn                             -0.460732  0.0154503    -1.66952  0.162484   1                8696   0.438356  FRASOR_RESPONSE_TO_ESTRADIOL_DN
            1  Tarte Plasma Cell Vs Plasmablast Dn                          0.641168  0.000403796   2.01072  0.0127991  1                1674   0.544484  TARTE_PLASMA_CELL_VS_PLASMABLAST_DN
            4  Tarte Plasma Cell Vs Plasmablast Dn                         -0.614694  0.00138149   -1.93551  0.155097   1                8933   0.587189  TARTE_PLASMA_CELL_VS_PLASMABLAST_DN
            5  Tarte Plasma Cell Vs Plasmablast Dn                          0.533477  0.0379393     1.67353  0.155982   1                2435   0.483986  TARTE_PLASMA_CELL_VS_PLASMABLAST_DN
           10  Tarte Plasma Cell Vs Plasmablast Dn                          0.533726  0.0359712     1.67995  0.149251   1                2538   0.523132  TARTE_PLASMA_CELL_VS_PLASMABLAST_DN
            1  Smith Liver Cancer                                           0.54732   0.0163444     1.71182  0.075365   1                2297   0.513514  SMITH_LIVER_CANCER
            1  Jechlinger Epithelial To Mesenchymal Transition Dn          -0.529871  0.00742823   -1.82037  0.119811   1                8948   0.525424  JECHLINGER_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_DN
            1  Le Egr2 Targets Up                                           0.741426  0.000406009   2.04057  0.0127991  1                 504   0.540816  LE_EGR2_TARGETS_UP
            4  Le Egr2 Targets Up                                          -0.647941  0.0232604    -1.7937   0.192508   1                9532   0.591837  LE_EGR2_TARGETS_UP
            5  Le Egr2 Targets Up                                           0.695197  0.00495442    1.9172   0.0475945  1                1031   0.530612  LE_EGR2_TARGETS_UP
           10  Le Egr2 Targets Up                                           0.679035  0.00724492    1.87723  0.0627622  1                 968   0.530612  LE_EGR2_TARGETS_UP
            1  Fernandez Bound By Myc                                       0.357474  0.0408937     1.45809  0.233099   1                2084   0.357143  FERNANDEZ_BOUND_BY_MYC
            1  Mootha Voxphos                                               0.580819  0.1332        1.48337  0.217553   1                3051   0.5875    MOOTHA_VOXPHOS
            4  Mootha Voxphos                                              -0.677522  0.0324195    -1.72078  0.192508   1                8664   0.75      MOOTHA_VOXPHOS
            1  Schuhmacher Myc Targets Up                                   0.666472  0.00319298    1.89135  0.0220871  1                2271   0.69863   SCHUHMACHER_MYC_TARGETS_UP
            4  Schuhmacher Myc Targets Up                                  -0.604246  0.0230342    -1.71217  0.192508   1                9169   0.561644  SCHUHMACHER_MYC_TARGETS_UP
           10  Schuhmacher Myc Targets Up                                   0.601185  0.0260956     1.71088  0.127651   1                2990   0.684932  SCHUHMACHER_MYC_TARGETS_UP
            1  Sweet Kras Targets Up                                       -0.624621  0.0235982    -1.76952  0.138176   1                8291   0.64      SWEET_KRAS_TARGETS_UP
            1  Shipp Dlbcl Vs Follicular Lymphoma Dn                       -0.572342  0.0197921    -1.70747  0.146665   1                8530   0.5       SHIPP_DLBCL_VS_FOLLICULAR_LYMPHOMA_DN
            1  Sana Response To Ifng Dn                                     0.469739  0.0443389     1.58391  0.150048   1                1593   0.421053  SANA_RESPONSE_TO_IFNG_DN
            4  Sana Response To Ifng Dn                                    -0.570434  0.00139804   -1.92891  0.155097   1                9501   0.447368  SANA_RESPONSE_TO_IFNG_DN
            1  Lee Liver Cancer Myc Dn                                     -0.45611   0.0304068    -1.55249  0.202734   1               10166   0.305556  LEE_LIVER_CANCER_MYC_DN
            1  Iizuka Liver Cancer Progression L1 G1 Up                    -0.654743  0.0155419    -1.69809  0.148746   1                8934   0.5625    IIZUKA_LIVER_CANCER_PROGRESSION_L1_G1_UP
            1  Pomeroy Medulloblastoma Prognosis Dn                         0.656622  0.00913423    1.79985  0.0417664  1                2343   0.589744  POMEROY_MEDULLOBLASTOMA_PROGNOSIS_DN
            4  Pomeroy Medulloblastoma Prognosis Dn                        -0.609861  0.0305958    -1.66659  0.214503   1               10113   0.358974  POMEROY_MEDULLOBLASTOMA_PROGNOSIS_DN
            1  Ferrando Lyl1 Neighbors                                     -0.71437   0.020028     -1.67208  0.161802   1                9275   0.666667  FERRANDO_LYL1_NEIGHBORS
            1  Peng Glucose Deprivation Dn                                  0.393394  0.0264228     1.56373  0.160717   1                2626   0.419118  PENG_GLUCOSE_DEPRIVATION_DN
            1  Vantveer Breast Cancer Metastasis Up                        -0.610715  0.0398122    -1.64225  0.170938   1                7298   0.877551  VANTVEER_BREAST_CANCER_METASTASIS_UP
            1  Coller Myc Targets Up                                        0.623549  0.0529894     1.56863  0.157043   1                2198   0.708333  COLLER_MYC_TARGETS_UP
            4  Coller Myc Targets Up                                       -0.67737   0.0164683    -1.71199  0.192508   1                9223   0.625     COLLER_MYC_TARGETS_UP
            1  Manalo Hypoxia Up                                           -0.493103  0.04375      -1.64543  0.170556   1                8790   0.451429  MANALO_HYPOXIA_UP
            1  Ren Bound By E2F                                             0.808134  0.000807265   1.94261  0.0164279  1                 745   0.673077  REN_BOUND_BY_E2F
            5  Ren Bound By E2F                                             0.765055  0.00441324    1.83275  0.0723024  1                1102   0.673077  REN_BOUND_BY_E2F
           10  Ren Bound By E2F                                             0.754431  0.00616425    1.83356  0.0688673  1                1423   0.730769  REN_BOUND_BY_E2F
            1  Shipp Dlbcl Vs Follicular Lymphoma Up                        0.73441   0.00301508    1.84624  0.0309982  1                2114   0.756098  SHIPP_DLBCL_VS_FOLLICULAR_LYMPHOMA_UP
            4  Shipp Dlbcl Vs Follicular Lymphoma Up                       -0.729967  0.00275645   -1.85341  0.173896   1                9266   0.707317  SHIPP_DLBCL_VS_FOLLICULAR_LYMPHOMA_UP
           10  Shipp Dlbcl Vs Follicular Lymphoma Up                        0.657769  0.0340366     1.66413  0.162408   1                3297   0.804878  SHIPP_DLBCL_VS_FOLLICULAR_LYMPHOMA_UP
            1  Peart Hdac Proliferation Cluster Dn                          0.552501  0.0125608     1.74834  0.0583926  1                1662   0.394366  PEART_HDAC_PROLIFERATION_CLUSTER_DN
            5  Peart Hdac Proliferation Cluster Dn                          0.581405  0.005002      1.8278   0.0744562  1                1772   0.43662   PEART_HDAC_PROLIFERATION_CLUSTER_DN
           10  Peart Hdac Proliferation Cluster Dn                          0.560607  0.00973433    1.77075  0.0985143  1                2042   0.464789  PEART_HDAC_PROLIFERATION_CLUSTER_DN
            1  Sasaki Adult T Cell Leukemia                                 0.476153  0.00843373    1.77204  0.0504442  1                1550   0.344371  SASAKI_ADULT_T_CELL_LEUKEMIA
            5  Sasaki Adult T Cell Leukemia                                 0.463213  0.0109627     1.71839  0.12811    1                1484   0.317881  SASAKI_ADULT_T_CELL_LEUKEMIA
           10  Sasaki Adult T Cell Leukemia                                 0.450672  0.0228778     1.66931  0.159594   1                2839   0.463576  SASAKI_ADULT_T_CELL_LEUKEMIA
            1  Peng Rapamycin Response Dn                                   0.577976  0.00360649    1.87954  0.0239806  1                2376   0.594714  PENG_RAPAMYCIN_RESPONSE_DN
            4  Peng Rapamycin Response Dn                                  -0.589732  0.00158228   -1.91898  0.155097   1                8933   0.572687  PENG_RAPAMYCIN_RESPONSE_DN
           10  Peng Rapamycin Response Dn                                   0.493744  0.0506784     1.60536  0.207217   1                3330   0.603524  PENG_RAPAMYCIN_RESPONSE_DN
            1  Bhattacharya Embryonic Stem Cell                             0.62372   0.00260887    1.89426  0.0220431  1                1507   0.485714  BHATTACHARYA_EMBRYONIC_STEM_CELL
            4  Bhattacharya Embryonic Stem Cell                            -0.571279  0.0211337    -1.72123  0.192508   1                9313   0.514286  BHATTACHARYA_EMBRYONIC_STEM_CELL
            5  Bhattacharya Embryonic Stem Cell                             0.569392  0.0223553     1.72249  0.127357   1                1504   0.442857  BHATTACHARYA_EMBRYONIC_STEM_CELL
           10  Bhattacharya Embryonic Stem Cell                             0.537053  0.0473067     1.62273  0.18979    1                2654   0.528571  BHATTACHARYA_EMBRYONIC_STEM_CELL
            1  Yao Hoxa10 Targets Via Progesterone Up                      -0.552568  0.00808898   -1.79881  0.124176   1                8352   0.589286  YAO_HOXA10_TARGETS_VIA_PROGESTERONE_UP
            1  Sansom Apc Targets Up                                        0.35424   0.0453286     1.45605  0.234831   1                 514   0.195652  SANSOM_APC_TARGETS_UP
            5  Sansom Apc Targets Up                                        0.38811   0.0143569     1.60597  0.20614    1                2415   0.402174  SANSOM_APC_TARGETS_UP
            1  Zhan Multiple Myeloma Subgroups                              0.795285  0.000405762   1.91787  0.0189688  1                1807   0.827586  ZHAN_MULTIPLE_MYELOMA_SUBGROUPS
            4  Zhan Multiple Myeloma Subgroups                             -0.720023  0.0102767    -1.7342   0.192508   1                8986   0.724138  ZHAN_MULTIPLE_MYELOMA_SUBGROUPS
            5  Zhan Multiple Myeloma Subgroups                              0.701775  0.016456      1.69595  0.142948   1                2796   0.862069  ZHAN_MULTIPLE_MYELOMA_SUBGROUPS
           10  Zhan Multiple Myeloma Subgroups                              0.727145  0.00521251    1.75964  0.104724   1                2062   0.724138  ZHAN_MULTIPLE_MYELOMA_SUBGROUPS
            1  Vernell Retinoblastoma Pathway Dn                           -0.612905  0.00869917   -1.75107  0.139914   1                9106   0.55      VERNELL_RETINOBLASTOMA_PATHWAY_DN
            1  Wang Immortalized By Hoxa9 And Meis1 Dn                     -0.657824  0.0538727    -1.59969  0.181764   1                9537   0.6875    WANG_IMMORTALIZED_BY_HOXA9_AND_MEIS1_DN
            1  Stossi Response To Estradiol                                -0.546475  0.0309922    -1.63254  0.170938   1                8829   0.405405  STOSSI_RESPONSE_TO_ESTRADIOL
            1  Yu Myc Targets Up                                            0.875337  0.000203874   1.9163   0.0189688  0.542508          656   0.783784  YU_MYC_TARGETS_UP
            4  Yu Myc Targets Up                                           -0.783478  0.0154508    -1.73381  0.192508   1               10138   0.702703  YU_MYC_TARGETS_UP
            5  Yu Myc Targets Up                                            0.837838  0.000994233   1.85097  0.0645773  1                 945   0.756757  YU_MYC_TARGETS_UP
           10  Yu Myc Targets Up                                            0.834911  0.00120144    1.8504   0.0652612  1                 986   0.756757  YU_MYC_TARGETS_UP
            1  Moreaux B Lymphocyte Maturation By Taci Dn                   0.55678   0.0204573     1.70323  0.0785025  1                1634   0.393939  MOREAUX_B_LYMPHOCYTE_MATURATION_BY_TACI_DN
            5  Moreaux B Lymphocyte Maturation By Taci Dn                   0.55572   0.0205466     1.70586  0.135933   1                1107   0.30303   MOREAUX_B_LYMPHOCYTE_MATURATION_BY_TACI_DN
            1  Haslinger B Cll With 13Q14 Deletion                         -0.517763  0.0421687    -1.52283  0.21968    1                8165   0.631579  HASLINGER_B_CLL_WITH_13Q14_DELETION
            5  Hofmann Cell Lymphoma Up                                     0.499785  0.0144956     1.66915  0.156755   1                1592   0.418605  HOFMANN_CELL_LYMPHOMA_UP
            1  Gildea Metastasis                                           -0.579175  0.0605998    -1.57191  0.195207   1                7889   0.689655  GILDEA_METASTASIS
            1  Zhan Multiple Myeloma Lb Up                                 -0.447054  0.0568767    -1.46693  0.249633   1                7339   0.606061  ZHAN_MULTIPLE_MYELOMA_LB_UP
            1  Brown Myeloid Cell Development Dn                            0.455597  0.0132805     1.69649  0.0811864  1                2291   0.494845  BROWN_MYELOID_CELL_DEVELOPMENT_DN
            5  Brown Myeloid Cell Development Dn                            0.440149  0.0211619     1.64298  0.181542   1                2330   0.443299  BROWN_MYELOID_CELL_DEVELOPMENT_DN
           10  Brown Myeloid Cell Development Dn                            0.429889  0.0259457     1.6018   0.207217   1                2434   0.412371  BROWN_MYELOID_CELL_DEVELOPMENT_DN
            1  Li Wilms Tumor Vs Fetal Kidney 1 Dn                          0.718335  0.000405351   2.10113  0.0127991  1                1088   0.555556  LI_WILMS_TUMOR_VS_FETAL_KIDNEY_1_DN
            5  Li Wilms Tumor Vs Fetal Kidney 1 Dn                          0.686603  0.000790202   2.01279  0.0336414  1                1102   0.534722  LI_WILMS_TUMOR_VS_FETAL_KIDNEY_1_DN
           10  Li Wilms Tumor Vs Fetal Kidney 1 Dn                          0.668863  0.00220132    1.97008  0.056431   1                1339   0.520833  LI_WILMS_TUMOR_VS_FETAL_KIDNEY_1_DN
            1  Halmos Cebpa Targets Dn                                     -0.543514  0.0216282    -1.69476  0.150339   1                7559   0.6       HALMOS_CEBPA_TARGETS_DN
            1  Ferrando T All With Mll Enl Fusion Dn                        0.643413  0.00141243    1.95323  0.015348   1                2160   0.586667  FERRANDO_T_ALL_WITH_MLL_ENL_FUSION_DN
            4  Ferrando T All With Mll Enl Fusion Dn                       -0.581924  0.0088845    -1.77391  0.192508   1                8971   0.533333  FERRANDO_T_ALL_WITH_MLL_ENL_FUSION_DN
            5  Ferrando T All With Mll Enl Fusion Dn                        0.524119  0.0502793     1.59269  0.215809   1                1690   0.4       FERRANDO_T_ALL_WITH_MLL_ENL_FUSION_DN
           10  Ferrando T All With Mll Enl Fusion Dn                        0.559635  0.0207202     1.69835  0.131604   1                1764   0.453333  FERRANDO_T_ALL_WITH_MLL_ENL_FUSION_DN
            1  Nemeth Inflammatory Response Lps Dn                          0.507168  0.0531332     1.53318  0.181681   1                1495   0.37037   NEMETH_INFLAMMATORY_RESPONSE_LPS_DN
            5  Nemeth Inflammatory Response Lps Dn                          0.530449  0.0304417     1.59183  0.215809   1                2325   0.555556  NEMETH_INFLAMMATORY_RESPONSE_LPS_DN
            1  Rorie Targets Of Ewsr1 Fli1 Fusion Up                       -0.604745  0.019774     -1.71295  0.146665   1                8660   0.541667  RORIE_TARGETS_OF_EWSR1_FLI1_FUSION_UP
            1  Zhan Multiple Myeloma Pr Up                                  0.971158  0.000201005   1.82853  0.0346151  0.534874          180   0.9375    ZHAN_MULTIPLE_MYELOMA_PR_UP
            4  Zhan Multiple Myeloma Pr Up                                 -0.894667  0.00908911   -1.6735   0.209628   1               10178   0.9375    ZHAN_MULTIPLE_MYELOMA_PR_UP
            5  Zhan Multiple Myeloma Pr Up                                  0.897177  0.0109562     1.67617  0.153585   1                 904   0.90625   ZHAN_MULTIPLE_MYELOMA_PR_UP
           10  Zhan Multiple Myeloma Pr Up                                  0.918636  0.00341297    1.73003  0.115709   1                 664   0.90625   ZHAN_MULTIPLE_MYELOMA_PR_UP
            1  Halmos Cebpa Targets Up                                     -0.505898  0.0663681    -1.52913  0.214421   1                7873   0.578947  HALMOS_CEBPA_TARGETS_UP
            1  Tavor Cebpa Targets Dn                                      -0.546324  0.0304487    -1.6156   0.175815   1                8921   0.454545  TAVOR_CEBPA_TARGETS_DN
            1  Wang Targets Of Mll Cbp Fusion Dn                            0.442312  0.0652303     1.46846  0.22633    1                2857   0.534884  WANG_TARGETS_OF_MLL_CBP_FUSION_DN
            4  Wang Targets Of Mll Cbp Fusion Dn                           -0.497707  0.0184554    -1.64901  0.221842   1                8936   0.44186   WANG_TARGETS_OF_MLL_CBP_FUSION_DN
            1  Nguyen Notch1 Targets Up                                    -0.471908  0.0292175    -1.54914  0.202734   1                8842   0.5       NGUYEN_NOTCH1_TARGETS_UP
            1  Abbud Lif Signaling 1 Up                                    -0.511454  0.0619164    -1.5498   0.202734   1                9526   0.368421  ABBUD_LIF_SIGNALING_1_UP
            1  Abraham Alpc Vs Multiple Myeloma Up                         -0.593211  0.0150997    -1.72457  0.146219   1                8011   0.8       ABRAHAM_ALPC_VS_MULTIPLE_MYELOMA_UP
            1  Gery Cebp Targets                                           -0.419682  0.0665449    -1.50449  0.229013   1                9898   0.287234  GERY_CEBP_TARGETS
            1  Kang Immortalized By Tert Up                                -0.437275  0.0490872    -1.52716  0.215992   1                9643   0.333333  KANG_IMMORTALIZED_BY_TERT_UP
            1  Petrova Endothelium Lymphatic Vs Blood Dn                   -0.484042  0.0500101    -1.61515  0.175815   1                8145   0.521739  PETROVA_ENDOTHELIUM_LYMPHATIC_VS_BLOOD_DN
            1  Jain Nfkb Signaling                                          0.512014  0.00101338    1.84445  0.0311256  1                2240   0.416667  JAIN_NFKB_SIGNALING
           10  Jain Nfkb Signaling                                          0.483699  0.00580581    1.74406  0.112179   1                2622   0.458333  JAIN_NFKB_SIGNALING
            1  Petrova Endothelium Lymphatic Vs Blood Up                    0.448901  0.0477713     1.56064  0.161827   1                 180   0.229358  PETROVA_ENDOTHELIUM_LYMPHATIC_VS_BLOOD_UP
            1  Li Wilms Tumor Vs Fetal Kidney 1 Up                         -0.508916  0.00305499   -1.86696  0.119811   1                8854   0.416667  LI_WILMS_TUMOR_VS_FETAL_KIDNEY_1_UP
            1  Staege Ewing Family Tumor                                   -0.608958  0.024027     -1.62506  0.170938   1                9014   0.692308  STAEGE_EWING_FAMILY_TUMOR
            5  Gale Apl With Flt3 Mutated Up                                0.464425  0.0278173     1.58323  0.225686   1                2870   0.462963  GALE_APL_WITH_FLT3_MUTATED_UP
            1  Affar Yy1 Targets Up                                        -0.462385  0.00122175   -1.86293  0.119811   1                8803   0.417323  AFFAR_YY1_TARGETS_UP
            1  Affar Yy1 Targets Dn                                         0.558149  0.00200401    1.95094  0.015456   1                1116   0.407643  AFFAR_YY1_TARGETS_DN
            4  Affar Yy1 Targets Dn                                        -0.478173  0.0428062    -1.66     0.21546    1                9381   0.401274  AFFAR_YY1_TARGETS_DN
            5  Affar Yy1 Targets Dn                                         0.514001  0.0137395     1.77924  0.0972596  1                1352   0.388535  AFFAR_YY1_TARGETS_DN
           10  Affar Yy1 Targets Dn                                         0.523837  0.0120895     1.81609  0.0752823  1                1236   0.382166  AFFAR_YY1_TARGETS_DN
            1  Vernell Retinoblastoma Pathway Up                            0.687446  0.00671687    1.86963  0.0262785  1                 879   0.578947  VERNELL_RETINOBLASTOMA_PATHWAY_UP
            5  Vernell Retinoblastoma Pathway Up                            0.734512  0.00158983    2.00108  0.0354642  1                1322   0.684211  VERNELL_RETINOBLASTOMA_PATHWAY_UP
           10  Vernell Retinoblastoma Pathway Up                            0.692288  0.00583971    1.88403  0.0627622  1                1082   0.578947  VERNELL_RETINOBLASTOMA_PATHWAY_UP
            1  Alcalay Aml By Npm1 Localization Dn                          0.460101  0.0226226     1.69228  0.0834695  1                1365   0.331081  ALCALAY_AML_BY_NPM1_LOCALIZATION_DN
            5  Alcalay Aml By Npm1 Localization Dn                          0.461716  0.0223807     1.68646  0.147682   1                1663   0.351351  ALCALAY_AML_BY_NPM1_LOCALIZATION_DN
           10  Alcalay Aml By Npm1 Localization Dn                          0.46536   0.0172241     1.70732  0.129418   1                1770   0.391892  ALCALAY_AML_BY_NPM1_LOCALIZATION_DN
            1  Kang Immortalized By Tert Dn                                -0.479336  0.030467     -1.62383  0.170938   1                8888   0.435484  KANG_IMMORTALIZED_BY_TERT_DN
            1  Hess Targets Of Hoxa9 And Meis1 Up                           0.570873  0.00360144    1.84042  0.0319879  1                1449   0.471698  HESS_TARGETS_OF_HOXA9_AND_MEIS1_UP
            5  Hess Targets Of Hoxa9 And Meis1 Up                           0.513519  0.025888      1.64067  0.182655   1                2036   0.45283   HESS_TARGETS_OF_HOXA9_AND_MEIS1_UP
           10  Hess Targets Of Hoxa9 And Meis1 Up                           0.549738  0.00992908    1.76877  0.0988582  1                2173   0.45283   HESS_TARGETS_OF_HOXA9_AND_MEIS1_UP
            1  Petrova Prox1 Targets Dn                                    -0.53983   0.0678478    -1.56641  0.198428   1                9556   0.407407  PETROVA_PROX1_TARGETS_DN
            1  Dorsey Gab2 Targets                                         -0.730996  0.0111925    -1.75126  0.139914   1                8888   0.722222  DORSEY_GAB2_TARGETS
            1  Zhan Multiple Myeloma Ms Up                                 -0.570219  0.0245966    -1.68903  0.153869   1                8603   0.542857  ZHAN_MULTIPLE_MYELOMA_MS_UP
            1  Radmacher Aml Prognosis                                     -0.46611   0.0166938    -1.64786  0.169344   1                8500   0.4375    RADMACHER_AML_PROGNOSIS
            1  Xu Crebbp Targets Up                                         0.548998  0.0186        1.65125  0.103686   1                1158   0.409091  XU_CREBBP_TARGETS_UP
           10  Xu Crebbp Targets Up                                         0.543777  0.0227181     1.64567  0.17341    1                1567   0.454545  XU_CREBBP_TARGETS_UP
            1  Verhaak Aml With Npm1 Mutated Dn                            -0.481051  0.0227685    -1.70379  0.147212   1                8291   0.518919  VERHAAK_AML_WITH_NPM1_MUTATED_DN
            1  Moreaux Multiple Myeloma By Taci Dn                          0.448512  0.0352042     1.60086  0.135148   1                2615   0.43871   MOREAUX_MULTIPLE_MYELOMA_BY_TACI_DN
            5  Moreaux Multiple Myeloma By Taci Dn                          0.499572  0.00504541    1.77939  0.0972596  1                2170   0.374194  MOREAUX_MULTIPLE_MYELOMA_BY_TACI_DN
            1  Schuringa Stat5A Targets Dn                                 -0.608833  0.0636492    -1.52166  0.220019   1                8626   0.416667  SCHURINGA_STAT5A_TARGETS_DN
            1  Kamminga Ezh2 Targets                                        0.86654   0.000201654   1.94709  0.0158608  0.5366            635   0.75      KAMMINGA_EZH2_TARGETS
            5  Kamminga Ezh2 Targets                                        0.793531  0.00797448    1.77356  0.0996009  1                1013   0.75      KAMMINGA_EZH2_TARGETS
           10  Kamminga Ezh2 Targets                                        0.800969  0.00420505    1.79806  0.0817333  1                 805   0.722222  KAMMINGA_EZH2_TARGETS
            1  Greenbaum E2A Targets Up                                     0.864045  0.000402982   1.96177  0.0147521  1                 535   0.71875   GREENBAUM_E2A_TARGETS_UP
            4  Greenbaum E2A Targets Up                                    -0.804301  0.00614105   -1.81877  0.188729   1               10351   0.6875    GREENBAUM_E2A_TARGETS_UP
            5  Greenbaum E2A Targets Up                                     0.825101  0.00295916    1.86197  0.0620503  1                 624   0.65625   GREENBAUM_E2A_TARGETS_UP
           10  Greenbaum E2A Targets Up                                     0.868721  0.000400721   1.9733   0.056431   1                 742   0.75      GREENBAUM_E2A_TARGETS_UP
            1  Xu Response To Tretinoin And Nsc682994 Dn                    0.751578  0.00815434    1.70354  0.0785025  1                2172   0.785714  XU_RESPONSE_TO_TRETINOIN_AND_NSC682994_DN
            1  Petrova Prox1 Targets Up                                     0.718227  0.00976205    1.76362  0.0535254  1                 180   0.458333  PETROVA_PROX1_TARGETS_UP
            4  Petrova Prox1 Targets Up                                    -0.72693   0.0084997    -1.80451  0.192508   1                9627   0.625     PETROVA_PROX1_TARGETS_UP
            5  Petrova Prox1 Targets Up                                     0.703079  0.0169158     1.73136  0.123266   1                1532   0.583333  PETROVA_PROX1_TARGETS_UP
           10  Petrova Prox1 Targets Up                                     0.673124  0.0310808     1.66071  0.164558   1                1382   0.625     PETROVA_PROX1_TARGETS_UP
            1  Faelt B Cll With Vh3 21 Up                                   0.558506  0.0312374     1.64919  0.103998   1                1377   0.395349  FAELT_B_CLL_WITH_VH3_21_UP
            1  Ly Aging Premature Dn                                        0.688436  0.0249695     1.67704  0.0902429  1                1111   0.62963   LY_AGING_PREMATURE_DN
           10  Ly Aging Premature Dn                                        0.708103  0.0142972     1.72933  0.115709   1                1423   0.666667  LY_AGING_PREMATURE_DN
            4  Takao Response To Uvb Radiation Up                          -0.564086  0.0547076    -1.61945  0.237236   1                8917   0.547945  TAKAO_RESPONSE_TO_UVB_RADIATION_UP
            1  Urs Adipocyte Differentiation Up                            -0.448162  0.0360794    -1.57633  0.194122   1                9042   0.470588  URS_ADIPOCYTE_DIFFERENTIATION_UP
            1  Traynor Rett Syndrom Up                                     -0.577644  0.0416751    -1.63119  0.170938   1                9459   0.484848  TRAYNOR_RETT_SYNDROM_UP
            5  Browne Hcmv Infection 16Hr Up                                0.394385  0.014998      1.60475  0.206924   1                2898   0.427083  BROWNE_HCMV_INFECTION_16HR_UP
            1  Jazaeri Breast Cancer Brca1 Vs Brca2 Dn                     -0.490296  0.0338346    -1.57719  0.194122   1                8510   0.485714  JAZAERI_BREAST_CANCER_BRCA1_VS_BRCA2_DN
            1  Zamora Nos2 Targets Up                                       0.638965  0.000995223   1.9336   0.0171414  1                2116   0.539683  ZAMORA_NOS2_TARGETS_UP
            4  Zamora Nos2 Targets Up                                      -0.623949  0.0017843    -1.88092  0.162281   1                9648   0.507937  ZAMORA_NOS2_TARGETS_UP
            4  Weigel Oxidative Stress Response                            -0.50678   0.0167541    -1.64254  0.225074   1                9512   0.451613  WEIGEL_OXIDATIVE_STRESS_RESPONSE
            1  Hedenfalk Breast Cancer Hereditary Vs Sporadic               0.560881  0.00622365    1.79989  0.0417664  1                2436   0.571429  HEDENFALK_BREAST_CANCER_HEREDITARY_VS_SPORADIC
            5  Hedenfalk Breast Cancer Hereditary Vs Sporadic               0.506377  0.0289281     1.63509  0.184221   1                2459   0.55102   HEDENFALK_BREAST_CANCER_HEREDITARY_VS_SPORADIC
           10  Hedenfalk Breast Cancer Hereditary Vs Sporadic               0.533756  0.0147147     1.72976  0.115709   1                2224   0.510204  HEDENFALK_BREAST_CANCER_HEREDITARY_VS_SPORADIC
            1  Brachat Response To Cisplatin                               -0.534753  0.0266585    -1.62335  0.170938   1                9090   0.4       BRACHAT_RESPONSE_TO_CISPLATIN
            1  Joseph Response To Sodium Butyrate Dn                       -0.443319  0.0358861    -1.54738  0.203357   1                8930   0.477273  JOSEPH_RESPONSE_TO_SODIUM_BUTYRATE_DN
            1  Browne Hcmv Infection 16Hr Dn                               -0.454717  0.0477837    -1.55684  0.202405   1                9591   0.412698  BROWNE_HCMV_INFECTION_16HR_DN
            1  Weston Vegfa Targets 6Hr                                    -0.67601   0.0104637    -1.8362   0.119811   1                9146   0.697674  WESTON_VEGFA_TARGETS_6HR
            1  Rhodes Cancer Meta Signature                                 0.661831  0.00782661    1.82065  0.0363179  1                2245   0.694915  RHODES_CANCER_META_SIGNATURE
            4  Rhodes Cancer Meta Signature                                -0.724536  0.000197355  -1.9927   0.155097   0.525163         9259   0.677966  RHODES_CANCER_META_SIGNATURE
            5  Rhodes Cancer Meta Signature                                 0.589742  0.0556779     1.62596  0.19014    1                3171   0.661017  RHODES_CANCER_META_SIGNATURE
            1  Nielsen Malignat Fibrous Histiocytoma Up                    -0.783718  0.0177563    -1.67905  0.157308   1                9172   0.8       NIELSEN_MALIGNAT_FIBROUS_HISTIOCYTOMA_UP
            1  Mcdowell Acute Lung Injury Up                               -0.591021  0.0433023    -1.64692  0.169344   1                9481   0.5       MCDOWELL_ACUTE_LUNG_INJURY_UP
            4  Sesto Response To Uv C4                                     -0.564617  0.022673     -1.63467  0.229428   1                8931   0.588235  SESTO_RESPONSE_TO_UV_C4
            1  Su Testis                                                    0.70717   0.00202265    1.90816  0.0202154  1                 635   0.409091  SU_TESTIS
            5  Su Testis                                                    0.703706  0.00454007    1.88816  0.0559509  1                2140   0.613636  SU_TESTIS
           10  Su Testis                                                    0.703521  0.0026252     1.89214  0.0625516  1                1724   0.568182  SU_TESTIS
            1  Mcclung Creb1 Targets Dn                                    -0.440725  0.0280821    -1.53576  0.211844   1                9353   0.425     MCCLUNG_CREB1_TARGETS_DN
            1  Kaab Heart Atrium Vs Ventricle Up                           -0.431286  0.0260799    -1.65002  0.168175   1                8711   0.435233  KAAB_HEART_ATRIUM_VS_VENTRICLE_UP
            1  Nielsen Gist Vs Synovial Sarcoma Up                         -0.676169  0.0682186    -1.57487  0.195122   1               10082   0.5       NIELSEN_GIST_VS_SYNOVIAL_SARCOMA_UP
            1  Vietor Ifrd1 Targets                                         0.623428  0.0130267     1.70306  0.0785025  1                 769   0.333333  VIETOR_IFRD1_TARGETS
            4  Vietor Ifrd1 Targets                                        -0.668576  0.004        -1.83486  0.187832   1               10457   0.388889  VIETOR_IFRD1_TARGETS
            1  Pal Prmt5 Targets Up                                         0.493648  0.00459908    1.83732  0.0321985  1                1460   0.377778  PAL_PRMT5_TARGETS_UP
            5  Pal Prmt5 Targets Up                                         0.443763  0.0324779     1.64444  0.180576   1                 942   0.294444  PAL_PRMT5_TARGETS_UP
           10  Pal Prmt5 Targets Up                                         0.438383  0.0354354     1.62331  0.18979    1                2137   0.405556  PAL_PRMT5_TARGETS_UP
            1  Zhang Response To Cantharidin Up                            -0.58639   0.052579     -1.54451  0.205051   1                7718   0.8125    ZHANG_RESPONSE_TO_CANTHARIDIN_UP
            1  Browne Hcmv Infection 24Hr Dn                               -0.509643  0.0195162    -1.74471  0.139914   1                8165   0.58871   BROWNE_HCMV_INFECTION_24HR_DN
            1  Chiba Response To Tsa Up                                    -0.585114  0.0118174    -1.80549  0.124176   1                8713   0.511111  CHIBA_RESPONSE_TO_TSA_UP
            1  Mody Hippocampus Neonatal                                    0.527331  0.0917961     1.46823  0.22633    1                1353   0.428571  MODY_HIPPOCAMPUS_NEONATAL
            4  Mody Hippocampus Neonatal                                   -0.625372  0.0157811    -1.7229   0.192508   1                9316   0.535714  MODY_HIPPOCAMPUS_NEONATAL
            1  Mcclung Cocaine Reward 5D                                   -0.403412  0.0204244    -1.53092  0.214421   1                9810   0.285714  MCCLUNG_COCAINE_REWARD_5D
            1  Chuang Oxidative Stress Response Dn                          0.794807  0.00583266    1.75095  0.0577278  1                 304   0.5       CHUANG_OXIDATIVE_STRESS_RESPONSE_DN
            5  Chuang Oxidative Stress Response Dn                          0.734716  0.0296754     1.60829  0.203818   1                1102   0.5       CHUANG_OXIDATIVE_STRESS_RESPONSE_DN
           10  Chuang Oxidative Stress Response Dn                          0.730779  0.0328751     1.61033  0.203347   1                 246   0.5       CHUANG_OXIDATIVE_STRESS_RESPONSE_DN
            1  Delaserna Myod Targets Dn                                   -0.47065   0.0359845    -1.60081  0.181764   1                8531   0.510204  DELASERNA_MYOD_TARGETS_DN
            1  Jiang Vhl Targets                                            0.352339  0.0455631     1.44499  0.238115   1                2654   0.415385  JIANG_VHL_TARGETS
            1  Ehrlich Icf Syndrom Up                                      -0.697827  0.0238386    -1.64016  0.170938   1                8681   0.6       EHRLICH_ICF_SYNDROM_UP
            1  Weston Vegfa Targets 12Hr                                   -0.682791  0.0175546    -1.76815  0.138176   1                9862   0.538462  WESTON_VEGFA_TARGETS_12HR
            1  Chuang Oxidative Stress Response Up                         -0.499438  0.0401867    -1.56837  0.197985   1                9617   0.416667  CHUANG_OXIDATIVE_STRESS_RESPONSE_UP
            4  Natsume Response To Interferon Beta Dn                      -0.541961  0.0288156    -1.65947  0.21546    1                8633   0.604651  NATSUME_RESPONSE_TO_INTERFERON_BETA_DN
            1  Baelde Diabetic Nephropathy Up                              -0.596951  0.0162009    -1.83206  0.119811   1                9012   0.6       BAELDE_DIABETIC_NEPHROPATHY_UP
            1  Lindvall Immortalized By Tert Up                            -0.463005  0.0216194    -1.63247  0.170938   1                9509   0.369231  LINDVALL_IMMORTALIZED_BY_TERT_UP
            1  Wang Cisplatin Response And Xpc Dn                          -0.357027  0.0120676    -1.57996  0.192439   1                9823   0.25      WANG_CISPLATIN_RESPONSE_AND_XPC_DN
            1  Browne Hcmv Infection 20Hr Dn                               -0.470052  0.0394658    -1.59687  0.181921   1                9606   0.352941  BROWNE_HCMV_INFECTION_20HR_DN
            1  Simbulan Parp1 Targets Up                                   -0.710019  0.00855049   -1.8082   0.124176   1                9171   0.666667  SIMBULAN_PARP1_TARGETS_UP
            1  Kuninger Igf1 Vs Pdgfb Targets Dn                           -0.513476  0.0332736    -1.61597  0.175815   1                9509   0.542857  KUNINGER_IGF1_VS_PDGFB_TARGETS_DN
            8  Ruan Response To Tnf Dn                                      0.585335  0.00079952    2.05543  0.104006   1                1555   0.439394  RUAN_RESPONSE_TO_TNF_DN
            1  Kayo Aging Muscle Dn                                         0.364361  0.0466051     1.46672  0.227528   1                1508   0.28      KAYO_AGING_MUSCLE_DN
            1  Nielsen Liposarcoma Up                                      -0.603645  0.0884885    -1.4662   0.249989   1                9580   0.454545  NIELSEN_LIPOSARCOMA_UP
            1  Inga Tp53 Targets                                           -0.690718  0.024351     -1.71006  0.146665   1                9061   0.8       INGA_TP53_TARGETS
            1  Mariadason Regulated By Histone Acetylation Up              -0.415473  0.0128308    -1.61222  0.176986   1                8584   0.411765  MARIADASON_REGULATED_BY_HISTONE_ACETYLATION_UP
            1  Lindvall Immortalized By Tert Dn                            -0.577741  0.037297     -1.69468  0.150339   1                9992   0.405797  LINDVALL_IMMORTALIZED_BY_TERT_DN
            1  Kang Doxorubicin Resistance Up                               0.903697  0.000805477   1.80947  0.0392887  1                 745   0.930233  KANG_DOXORUBICIN_RESISTANCE_UP
            4  Kang Doxorubicin Resistance Up                              -0.825198  0.0262898    -1.64382  0.224407   1               10178   0.767442  KANG_DOXORUBICIN_RESISTANCE_UP
            5  Kang Doxorubicin Resistance Up                               0.885605  0.00356577    1.76548  0.102916   1                 885   0.906977  KANG_DOXORUBICIN_RESISTANCE_UP
           10  Kang Doxorubicin Resistance Up                               0.874739  0.00481541    1.7459   0.111255   1                 678   0.860465  KANG_DOXORUBICIN_RESISTANCE_UP
            1  Browne Hcmv Infection 48Hr Up                                0.346555  0.0420858     1.44723  0.238115   1                1347   0.25      BROWNE_HCMV_INFECTION_48HR_UP
            5  Browne Hcmv Infection 48Hr Up                                0.412489  0.00282771    1.71448  0.128776   1                2183   0.378788  BROWNE_HCMV_INFECTION_48HR_UP
            4  Burton Adipogenesis 5                                       -0.504285  0.0148356    -1.73306  0.192508   1                9161   0.478261  BURTON_ADIPOGENESIS_5
            1  Weston Vegfa Targets                                        -0.551763  0.0191915    -1.76595  0.138176   1                9556   0.493827  WESTON_VEGFA_TARGETS
            1  Chen Lvad Support Of Failing Heart Up                       -0.586554  0.00765049   -1.83267  0.119811   1                9012   0.535714  CHEN_LVAD_SUPPORT_OF_FAILING_HEART_UP
            1  Verrecchia Response To Tgfb1 C2                             -0.631635  0.0742834    -1.56759  0.197985   1                9747   0.478261  VERRECCHIA_RESPONSE_TO_TGFB1_C2
            1  Brachat Response To Methotrexate Up                         -0.492375  0.0424354    -1.52142  0.220019   1                8587   0.47619   BRACHAT_RESPONSE_TO_METHOTREXATE_UP
            1  Rhodes Undifferentiated Cancer                               0.796703  0.00182815    1.77618  0.0492564  1                1484   0.836066  RHODES_UNDIFFERENTIATED_CANCER
            4  Rhodes Undifferentiated Cancer                              -0.806568  0.00117005   -1.82201  0.187832   1                9977   0.803279  RHODES_UNDIFFERENTIATED_CANCER
            5  Rhodes Undifferentiated Cancer                               0.715927  0.0545309     1.5987   0.20961    1                1096   0.606557  RHODES_UNDIFFERENTIATED_CANCER
           10  Rhodes Undifferentiated Cancer                               0.723972  0.0447252     1.61703  0.195005   1                1727   0.721311  RHODES_UNDIFFERENTIATED_CANCER
            1  Browne Hcmv Infection 14Hr Up                                0.364915  0.0552743     1.44091  0.242467   1                2009   0.3       BROWNE_HCMV_INFECTION_14HR_UP
            1  Yamazaki Tceb3 Targets Dn                                    0.435451  0.00967547    1.67942  0.0889374  1                2469   0.465116  YAMAZAKI_TCEB3_TARGETS_DN
            5  Yamazaki Tceb3 Targets Dn                                    0.442829  0.00876145    1.72527  0.126645   1                2297   0.412791  YAMAZAKI_TCEB3_TARGETS_DN
           10  Yamazaki Tceb3 Targets Dn                                    0.424346  0.0157853     1.64476  0.17341    1                2695   0.459302  YAMAZAKI_TCEB3_TARGETS_DN
            1  Lu Aging Brain Up                                           -0.414009  0.0274028    -1.60807  0.180545   1                8216   0.419913  LU_AGING_BRAIN_UP
            1  Burton Adipogenesis 3                                        0.785686  0.000604839   1.97665  0.013915   1                 635   0.677778  BURTON_ADIPOGENESIS_3
            4  Burton Adipogenesis 3                                       -0.673444  0.0493852    -1.68859  0.201047   1                9689   0.588889  BURTON_ADIPOGENESIS_3
            5  Burton Adipogenesis 3                                        0.761338  0.00378637    1.89733  0.0541218  1                1031   0.677778  BURTON_ADIPOGENESIS_3
           10  Burton Adipogenesis 3                                        0.744891  0.00859484    1.87825  0.0627622  1                1259   0.688889  BURTON_ADIPOGENESIS_3
            1  Burton Adipogenesis 8                                       -0.459543  0.0919726    -1.48845  0.237801   1                9491   0.342105  BURTON_ADIPOGENESIS_8
            1  Verrecchia Early Response To Tgfb1                          -0.540769  0.0771704    -1.58549  0.189217   1                9707   0.418182  VERRECCHIA_EARLY_RESPONSE_TO_TGFB1
            1  Zhu Cmv 8 Hr Dn                                             -0.513494  0.0106511    -1.73314  0.141971   1                9215   0.490196  ZHU_CMV_8_HR_DN
            1  Zhu Cmv 24 Hr Dn                                            -0.600796  0.0371876    -1.71676  0.146665   1                9596   0.55      ZHU_CMV_24_HR_DN
            1  Pal Prmt5 Targets Dn                                        -0.587574  0.00968473   -1.7443   0.139914   1               10192   0.333333  PAL_PRMT5_TARGETS_DN
            1  Browne Hcmv Infection 18Hr Dn                               -0.488547  0.0183561    -1.71505  0.146665   1                9586   0.393701  BROWNE_HCMV_INFECTION_18HR_DN
            4  Moreira Response To Tsa Up                                  -0.64021   0.00739408   -1.79406  0.192508   1                9226   0.555556  MOREIRA_RESPONSE_TO_TSA_UP
            1  Gentile Uv Response Cluster D6                              -0.473816  0.0355922    -1.57226  0.195207   1                9692   0.363636  GENTILE_UV_RESPONSE_CLUSTER_D6
            1  Zhang Response To Cantharidin Dn                             0.650813  0.00179784    1.91514  0.0190217  1                1661   0.546875  ZHANG_RESPONSE_TO_CANTHARIDIN_DN
           10  Zhang Response To Cantharidin Dn                             0.554415  0.0361711     1.63633  0.177567   1                2102   0.484375  ZHANG_RESPONSE_TO_CANTHARIDIN_DN
            1  Burton Adipogenesis Peak At 24Hr                             0.807616  0.00200884    1.84732  0.0308973  1                1180   0.789474  BURTON_ADIPOGENESIS_PEAK_AT_24HR
            4  Burton Adipogenesis Peak At 24Hr                            -0.784411  0.00706853   -1.78845  0.192508   1               10138   0.763158  BURTON_ADIPOGENESIS_PEAK_AT_24HR
            5  Burton Adipogenesis Peak At 24Hr                             0.710673  0.0562        1.61079  0.201237   1                1484   0.605263  BURTON_ADIPOGENESIS_PEAK_AT_24HR
           10  Burton Adipogenesis Peak At 24Hr                             0.782235  0.00859828    1.7944   0.0824479  1                 974   0.684211  BURTON_ADIPOGENESIS_PEAK_AT_24HR
            1  Brachat Response To Camptothecin Up                         -0.501366  0.0218724    -1.6238   0.170938   1               10594   0.25      BRACHAT_RESPONSE_TO_CAMPTOTHECIN_UP
            1  Georgantas Hsc Markers                                      -0.599095  0.00240722   -1.89437  0.119811   1                9081   0.54902   GEORGANTAS_HSC_MARKERS
            1  Han Jnk Singaling Up                                        -0.594218  0.0267337    -1.6795   0.157308   1                9799   0.37931   HAN_JNK_SINGALING_UP
            1  Mclachlan Dental Caries Dn                                  -0.546921  0.007242     -1.76544  0.138176   1                7927   0.603448  MCLACHLAN_DENTAL_CARIES_DN
            1  Burton Adipogenesis Peak At 2Hr                             -0.61807   0.0352486    -1.70699  0.146665   1                9602   0.521739  BURTON_ADIPOGENESIS_PEAK_AT_2HR
            1  Weigel Oxidative Stress By Hne And Tbh                      -0.409239  0.0187925    -1.5664   0.198428   1                9518   0.320755  WEIGEL_OXIDATIVE_STRESS_BY_HNE_AND_TBH
            1  Tseng Irs1 Targets Dn                                       -0.519286  0.00163432   -1.92921  0.116729   1                9747   0.39      TSENG_IRS1_TARGETS_DN
            1  Rodwell Aging Kidney No Blood Up                            -0.474652  0.0645557    -1.58559  0.189217   1                9198   0.395722  RODWELL_AGING_KIDNEY_NO_BLOOD_UP
            1  Burton Adipogenesis 4                                        0.457483  0.0709206     1.47782  0.219334   1                1951   0.454545  BURTON_ADIPOGENESIS_4
            4  Burton Adipogenesis 4                                       -0.544242  0.00674068   -1.76656  0.192508   1                8408   0.636364  BURTON_ADIPOGENESIS_4
            1  Zhu Cmv All Dn                                              -0.554097  0.0298507    -1.73554  0.141912   1                9083   0.551724  ZHU_CMV_ALL_DN
            1  Ivanova Hematopoiesis Intermediate Progenitor                0.369721  0.0478183     1.46543  0.227839   1                3014   0.463415  IVANOVA_HEMATOPOIESIS_INTERMEDIATE_PROGENITOR
            1  Blalock Alzheimers Disease Incipient Dn                      0.352422  0.047505      1.46089  0.23075    1                1848   0.302158  BLALOCK_ALZHEIMERS_DISEASE_INCIPIENT_DN
            1  Song Targets Of Ie86 Cmv Protein                             0.819504  0.000805153   1.90706  0.0202154  1                1314   0.807692  SONG_TARGETS_OF_IE86_CMV_PROTEIN
            5  Song Targets Of Ie86 Cmv Protein                             0.763671  0.0104146     1.76527  0.102916   1                1980   0.826923  SONG_TARGETS_OF_IE86_CMV_PROTEIN
           10  Song Targets Of Ie86 Cmv Protein                             0.794647  0.00281917    1.8505   0.0652612  1                1423   0.75      SONG_TARGETS_OF_IE86_CMV_PROTEIN
            1  Mms Mouse Lymph High 4Hrs Up                                 0.583022  0.0115584     1.72928  0.0663823  1                2356   0.558824  MMS_MOUSE_LYMPH_HIGH_4HRS_UP
           10  Mms Mouse Lymph High 4Hrs Up                                 0.58465   0.0105945     1.75157  0.108321   1                1562   0.441176  MMS_MOUSE_LYMPH_HIGH_4HRS_UP
            1  Burton Adipogenesis 1                                       -0.580719  0.0754679    -1.55801  0.201425   1                8283   0.586207  BURTON_ADIPOGENESIS_1
            1  Liu Smarca4 Targets                                         -0.572661  0.012477     -1.75498  0.139914   1                8122   0.565217  LIU_SMARCA4_TARGETS
            1  Hedenfalk Breast Cancer Brca1 Vs Brca2                       0.466464  0.00708933    1.7593   0.0544691  1                2198   0.425806  HEDENFALK_BREAST_CANCER_BRCA1_VS_BRCA2
            5  Hedenfalk Breast Cancer Brca1 Vs Brca2                       0.445223  0.0119355     1.6901   0.145171   1                2500   0.458065  HEDENFALK_BREAST_CANCER_BRCA1_VS_BRCA2
           10  Hedenfalk Breast Cancer Brca1 Vs Brca2                       0.468079  0.0070014     1.77432  0.0959699  1                2224   0.451613  HEDENFALK_BREAST_CANCER_BRCA1_VS_BRCA2
            1  Brachat Response To Camptothecin Dn                          0.470175  0.0687575     1.49523  0.212964   1                 943   0.333333  BRACHAT_RESPONSE_TO_CAMPTOTHECIN_DN
            1  Burton Adipogenesis Peak At 0Hr                             -0.517102  0.0572134    -1.59943  0.181764   1                8842   0.509091  BURTON_ADIPOGENESIS_PEAK_AT_0HR
            1  Jazaeri Breast Cancer Brca1 Vs Brca2 Up                      0.441017  0.0627874     1.47894  0.219334   1                2416   0.372093  JAZAERI_BREAST_CANCER_BRCA1_VS_BRCA2_UP
            4  Jazaeri Breast Cancer Brca1 Vs Brca2 Up                     -0.48948   0.0212724    -1.63361  0.22975    1                9170   0.418605  JAZAERI_BREAST_CANCER_BRCA1_VS_BRCA2_UP
            1  Kalma E2F1 Targets                                           0.90237   0.000811853   1.72888  0.0663823  1                 270   0.818182  KALMA_E2F1_TARGETS
            5  Kalma E2F1 Targets                                           0.879155  0.00501605    1.68876  0.145831   1                 910   0.818182  KALMA_E2F1_TARGETS
           10  Kalma E2F1 Targets                                           0.854752  0.0102513     1.65496  0.166946   1                1625   1         KALMA_E2F1_TARGETS
            1  Verrecchia Response To Tgfb1 C4                             -0.631429  0.0884654    -1.47616  0.244495   1                9642   0.416667  VERRECCHIA_RESPONSE_TO_TGFB1_C4
            1  Ramaswamy Metastasis Up                                      0.497788  0.025005      1.66695  0.0958915  1                1267   0.365079  RAMASWAMY_METASTASIS_UP
            1  Kayo Calorie Restriction Muscle Up                          -0.501062  0.040282     -1.62384  0.170938   1                9506   0.366197  KAYO_CALORIE_RESTRICTION_MUSCLE_UP
            1  Ly Aging Old Dn                                              0.704567  0.00978394    1.7812   0.0477546  1                1514   0.705882  LY_AGING_OLD_DN
            4  Ly Aging Old Dn                                             -0.694601  0.0163934    -1.75858  0.192508   1                9565   0.627451  LY_AGING_OLD_DN
            5  Ly Aging Old Dn                                              0.66408   0.0328882     1.67211  0.15633    1                2095   0.666667  LY_AGING_OLD_DN
           10  Ly Aging Old Dn                                              0.680343  0.0208584     1.72521  0.117315   1                1189   0.627451  LY_AGING_OLD_DN
            1  Wang Cisplatin Response And Xpc Up                           0.547543  0.000400641   2.00796  0.0127991  1                1431   0.385135  WANG_CISPLATIN_RESPONSE_AND_XPC_UP
            5  Wang Cisplatin Response And Xpc Up                           0.491053  0.00575854    1.78993  0.0965026  1                1031   0.324324  WANG_CISPLATIN_RESPONSE_AND_XPC_UP
           10  Wang Cisplatin Response And Xpc Up                           0.529663  0.000798563   1.93698  0.056431   1                1546   0.391892  WANG_CISPLATIN_RESPONSE_AND_XPC_UP
            1  Varela Zmpste24 Targets Up                                  -0.590182  0.00141901   -1.9195   0.116729   1                9090   0.5       VARELA_ZMPSTE24_TARGETS_UP
            1  Cheng Response To Nickel Acetate                             0.422847  0.0293651     1.51103  0.19957    1                1924   0.4       CHENG_RESPONSE_TO_NICKEL_ACETATE
            1  Tseng Adipogenic Potential Dn                               -0.702173  0.000203998  -2.0745   0.116729   0.54284          9888   0.575758  TSENG_ADIPOGENIC_POTENTIAL_DN
            5  Dazard Uv Response Cluster G6                                0.469459  0.0261491     1.63744  0.184004   1                2525   0.443662  DAZARD_UV_RESPONSE_CLUSTER_G6
           10  Sesto Response To Uv C7                                      0.470743  0.0297189     1.58264  0.228294   1                2108   0.431034  SESTO_RESPONSE_TO_UV_C7
            1  Simbulan Parp1 Targets Dn                                    0.825554  0.00483871    1.74943  0.058262   1                1162   0.785714  SIMBULAN_PARP1_TARGETS_DN
            4  Simbulan Parp1 Targets Dn                                   -0.826917  0.00453381   -1.74713  0.192508   1                9689   0.857143  SIMBULAN_PARP1_TARGETS_DN
            5  Simbulan Parp1 Targets Dn                                    0.791261  0.0157822     1.67917  0.152752   1                 910   0.642857  SIMBULAN_PARP1_TARGETS_DN
           10  Simbulan Parp1 Targets Dn                                    0.760081  0.0294476     1.6035   0.207217   1                1330   0.714286  SIMBULAN_PARP1_TARGETS_DN
            1  Safford T Lymphocyte Anergy                                 -0.454068  0.0109671    -1.6791   0.157308   1                8744   0.409836  SAFFORD_T_LYMPHOCYTE_ANERGY
            1  Ly Aging Middle Dn                                           0.937594  0.00120289    1.7454   0.0591965  1                 158   0.866667  LY_AGING_MIDDLE_DN
            4  Ly Aging Middle Dn                                          -0.924844  0.00294522   -1.71589  0.192508   1               10537   0.866667  LY_AGING_MIDDLE_DN
            5  Ly Aging Middle Dn                                           0.849444  0.0366047     1.56704  0.242478   1                 425   0.733333  LY_AGING_MIDDLE_DN
           10  Ly Aging Middle Dn                                           0.869974  0.0182667     1.60823  0.204197   1                1189   0.933333  LY_AGING_MIDDLE_DN
            1  Bacolod Resistance To Alkylating Agents Up                  -0.638826  0.0175227    -1.71466  0.146665   1               10461   0.388889  BACOLOD_RESISTANCE_TO_ALKYLATING_AGENTS_UP
            1  Hu Genotoxic Damage 4Hr                                      0.739342  0.00245048    1.88275  0.0235955  1                1495   0.666667  HU_GENOTOXIC_DAMAGE_4HR
            4  Hu Genotoxic Damage 4Hr                                     -0.676794  0.0201064    -1.72288  0.192508   1                8883   0.733333  HU_GENOTOXIC_DAMAGE_4HR
            5  Hu Genotoxic Damage 4Hr                                      0.687219  0.0156469     1.74648  0.116739   1                1096   0.566667  HU_GENOTOXIC_DAMAGE_4HR
           10  Hu Genotoxic Damage 4Hr                                      0.689719  0.0146381     1.75792  0.105571   1                1561   0.633333  HU_GENOTOXIC_DAMAGE_4HR
            1  Mcdowell Acute Lung Injury Dn                               -0.514864  0.0324083    -1.62344  0.170938   1                9217   0.424242  MCDOWELL_ACUTE_LUNG_INJURY_DN
            1  Wang Smarce1 Targets Up                                     -0.609061  0.0106597    -1.86788  0.119811   1                8192   0.651982  WANG_SMARCE1_TARGETS_UP
            1  Han Jnk Singaling Dn                                        -0.522851  0.0509061    -1.56398  0.199134   1                9606   0.387097  HAN_JNK_SINGALING_DN
            1  Gentile Response Cluster D3                                  0.5016    0.0221996     1.65922  0.0991974  1                 610   0.232143  GENTILE_RESPONSE_CLUSTER_D3
            1  Lee Aging Muscle Up                                         -0.483454  0.0180601    -1.65866  0.16369    1                9172   0.4       LEE_AGING_MUSCLE_UP
            1  Mulligan Ntf3 Signaling Via Insr And Igf1R Dn               -0.643472  0.0306143    -1.60773  0.180545   1                9318   0.5       MULLIGAN_NTF3_SIGNALING_VIA_INSR_AND_IGF1R_DN
            1  Chen Etv5 Targets Testis                                     0.798704  0.00738007    1.73179  0.0657044  1                1588   0.857143  CHEN_ETV5_TARGETS_TESTIS
            5  Chen Etv5 Targets Testis                                     0.882038  0.000197511   1.92359  0.0475578  0.525578          535   0.714286  CHEN_ETV5_TARGETS_TESTIS
           10  Chen Etv5 Targets Testis                                     0.780515  0.0122416     1.69313  0.136131   1                1487   0.714286  CHEN_ETV5_TARGETS_TESTIS
            1  Welcsh Brca1 Targets Dn                                      0.591993  0.0020016     1.93397  0.0171414  1                2046   0.507692  WELCSH_BRCA1_TARGETS_DN
           10  Welcsh Brca1 Targets Dn                                      0.516593  0.0285886     1.69826  0.131604   1                3103   0.553846  WELCSH_BRCA1_TARGETS_DN
            1  Dasu Il6 Signaling Scar Dn                                  -0.657541  0.081869     -1.52152  0.220019   1                7870   0.866667  DASU_IL6_SIGNALING_SCAR_DN
            1  Martinez Response To Trabectedin Dn                          0.384149  0.0558878     1.48815  0.215694   1                2291   0.352227  MARTINEZ_RESPONSE_TO_TRABECTEDIN_DN
            5  Martinez Response To Trabectedin Dn                          0.444083  0.00957815    1.73035  0.123397   1                2909   0.469636  MARTINEZ_RESPONSE_TO_TRABECTEDIN_DN
            1  Burton Adipogenesis 9                                       -0.536291  0.00564289   -1.82943  0.119811   1                8660   0.4875    BURTON_ADIPOGENESIS_9
            1  Urs Adipocyte Differentiation Dn                            -0.560925  0.101778     -1.49146  0.235557   1                9967   0.36      URS_ADIPOCYTE_DIFFERENTIATION_DN
            1  Visala Response To Heat Shock And Aging Dn                   0.661929  0.0267539     1.63452  0.110528   1                2084   0.615385  VISALA_RESPONSE_TO_HEAT_SHOCK_AND_AGING_DN
            4  Visala Response To Heat Shock And Aging Dn                  -0.700023  0.00962503   -1.74076  0.192508   1                9990   0.538462  VISALA_RESPONSE_TO_HEAT_SHOCK_AND_AGING_DN
            1  Browne Hcmv Infection 12Hr Dn                               -0.426209  0.0470821    -1.52279  0.21968    1                9448   0.364706  BROWNE_HCMV_INFECTION_12HR_DN
            1  Mahajan Response To Il1A Up                                 -0.484644  0.0241624    -1.6308   0.170938   1                8131   0.490566  MAHAJAN_RESPONSE_TO_IL1A_UP
            1  Jiang Aging Hypothalamus Up                                  0.497123  0.0713281     1.54723  0.171988   1                2719   0.595238  JIANG_AGING_HYPOTHALAMUS_UP
            4  Jiang Aging Hypothalamus Up                                 -0.549437  0.0263635    -1.69374  0.19736    1                9017   0.595238  JIANG_AGING_HYPOTHALAMUS_UP
            1  Weston Vegfa Targets 3Hr                                    -0.46095   0.0596527    -1.53013  0.214421   1                9279   0.428571  WESTON_VEGFA_TARGETS_3HR
            1  Gerhold Adipogenesis Up                                     -0.456562  0.0648692    -1.50115  0.230381   1                9003   0.431818  GERHOLD_ADIPOGENESIS_UP
            1  Burton Adipogenesis Peak At 16Hr                             0.715155  0.00463336    1.88657  0.0228257  1                 917   0.540541  BURTON_ADIPOGENESIS_PEAK_AT_16HR
            4  Burton Adipogenesis Peak At 16Hr                            -0.619409  0.0538993    -1.62525  0.230529   1                8834   0.594595  BURTON_ADIPOGENESIS_PEAK_AT_16HR
            5  Burton Adipogenesis Peak At 16Hr                             0.734985  0.000798403   1.93668  0.0454213  1                1095   0.567568  BURTON_ADIPOGENESIS_PEAK_AT_16HR
           10  Burton Adipogenesis Peak At 16Hr                             0.66353   0.0197526     1.75087  0.108321   1                1201   0.513514  BURTON_ADIPOGENESIS_PEAK_AT_16HR
            1  Bild E2F3 Oncogenic Signature                                0.353417  0.0193443     1.53289  0.181681   1                1508   0.248588  BILD_E2F3_ONCOGENIC_SIGNATURE
           10  Bild E2F3 Oncogenic Signature                                0.369708  0.0109003     1.59247  0.215799   1                2317   0.338983  BILD_E2F3_ONCOGENIC_SIGNATURE
            1  Durchdewald Skin Carcinogenesis Up                           0.374544  0.0532213     1.43363  0.248209   1                1053   0.236364  DURCHDEWALD_SKIN_CARCINOGENESIS_UP
            1  Li Cytidine Analogs Cyctotoxicity                            0.54514   0.0349398     1.57092  0.156915   1                2499   0.428571  LI_CYTIDINE_ANALOGS_CYCTOTOXICITY
            1  Li Cytidine Analog Pathway                                   0.574685  0.057577      1.49339  0.213453   1                2144   0.583333  LI_CYTIDINE_ANALOG_PATHWAY
            1  Zhong Secretome Of Lung Cancer And Macrophage                0.461034  0.0817183     1.47878  0.219334   1                2047   0.41791   ZHONG_SECRETOME_OF_LUNG_CANCER_AND_MACROPHAGE
            4  Zhong Secretome Of Lung Cancer And Macrophage               -0.548413  0.0111955    -1.76774  0.192508   1                8890   0.522388  ZHONG_SECRETOME_OF_LUNG_CANCER_AND_MACROPHAGE
            4  Zhong Secretome Of Lung Cancer And Endothelium              -0.502161  0.0177256    -1.69472  0.19736    1                9348   0.433333  ZHONG_SECRETOME_OF_LUNG_CANCER_AND_ENDOTHELIUM
            4  Zhong Secretome Of Lung Cancer And Fibroblast               -0.515887  0.00841346   -1.80572  0.192508   1                9064   0.440678  ZHONG_SECRETOME_OF_LUNG_CANCER_AND_FIBROBLAST
            1  Clasper Lymphatic Vessels During Metastasis Up              -0.670056  0.0362319    -1.63684  0.170938   1                9261   0.736842  CLASPER_LYMPHATIC_VESSELS_DURING_METASTASIS_UP
            1  Clasper Lymphatic Vessels During Metastasis Dn              -0.825529  0.00584089   -1.80389  0.124176   1               10204   0.7       CLASPER_LYMPHATIC_VESSELS_DURING_METASTASIS_DN
            4  Singh Nfe2L2 Targets                                        -0.700369  0.0161968    -1.70294  0.192508   1                9523   0.642857  SINGH_NFE2L2_TARGETS
            4  Pellicciotta Hdac In Antigen Presentation Up                -0.528004  0.0396904    -1.63877  0.229428   1                9605   0.409836  PELLICCIOTTA_HDAC_IN_ANTIGEN_PRESENTATION_UP
            1  Pellicciotta Hdac In Antigen Presentation Dn                 0.555846  0.106301      1.48431  0.217372   1                3487   0.73913   PELLICCIOTTA_HDAC_IN_ANTIGEN_PRESENTATION_DN
            1  Harris Brain Cancer Progenitors                             -0.758942  0.00121457   -1.91611  0.116729   1                9812   0.6       HARRIS_BRAIN_CANCER_PROGENITORS
            1  Howlin Cited1 Targets 1 Dn                                  -0.472572  0.0691447    -1.46681  0.249633   1                8277   0.451613  HOWLIN_CITED1_TARGETS_1_DN
            1  Lein Astrocyte Markers                                      -0.624322  0.00530179   -1.81671  0.119811   1               10117   0.52      LEIN_ASTROCYTE_MARKERS
            1  Lein Choroid Plexus Markers                                 -0.472267  0.0320097    -1.63061  0.170938   1                8331   0.507692  LEIN_CHOROID_PLEXUS_MARKERS
            1  Lein Midbrain Markers                                       -0.478903  0.0200884    -1.62136  0.170938   1                9248   0.386364  LEIN_MIDBRAIN_MARKERS
            2  Gavin Foxp3 Targets Cluster T7                              -0.561919  0.000800961  -2.00565  0.202024   1                8552   0.511111  GAVIN_FOXP3_TARGETS_CLUSTER_T7
            1  Gavin Foxp3 Targets Cluster P6                               0.772905  0.00120992    2.04759  0.0127991  1                 283   0.557377  GAVIN_FOXP3_TARGETS_CLUSTER_P6
            4  Gavin Foxp3 Targets Cluster P6                              -0.631561  0.0656287    -1.67273  0.209628   1               10261   0.52459   GAVIN_FOXP3_TARGETS_CLUSTER_P6
            5  Gavin Foxp3 Targets Cluster P6                               0.704113  0.0169826     1.85954  0.0620503  1                 942   0.557377  GAVIN_FOXP3_TARGETS_CLUSTER_P6
           10  Gavin Foxp3 Targets Cluster P6                               0.745484  0.00261569    1.97898  0.056431   1                 714   0.57377   GAVIN_FOXP3_TARGETS_CLUSTER_P6
            1  Gavin Foxp3 Targets Cluster P7                              -0.480397  0.0460553    -1.56527  0.199134   1                8042   0.476923  GAVIN_FOXP3_TARGETS_CLUSTER_P7
            5  Gavin Il2 Responsive Foxp3 Targets Up                        0.57187   0.0208459     1.63534  0.184221   1                1775   0.4       GAVIN_IL2_RESPONSIVE_FOXP3_TARGETS_UP
            1  Zheng Foxp3 Targets In T Lymphocyte Dn                      -0.599332  0.0348348    -1.61235  0.176986   1                8558   0.642857  ZHENG_FOXP3_TARGETS_IN_T_LYMPHOCYTE_DN
            1  Zheng Foxp3 Targets Up                                      -0.545523  0.083098     -1.47619  0.244495   1                7370   0.65      ZHENG_FOXP3_TARGETS_UP
            1  Sansom Apc Myc Targets                                       0.356351  0.0411148     1.47841  0.219334   1                2390   0.326203  SANSOM_APC_MYC_TARGETS
            1  Sansom Apc Targets Require Myc                               0.489572  0.00678778    1.79418  0.0440822  1                2075   0.404624  SANSOM_APC_TARGETS_REQUIRE_MYC
            1  Stein Esrra Targets Responsive To Estrogen Dn                0.731359  0.00159936    1.98834  0.0134874  1                 739   0.542857  STEIN_ESRRA_TARGETS_RESPONSIVE_TO_ESTROGEN_DN
            5  Stein Esrra Targets Responsive To Estrogen Dn                0.766117  0.000795545   2.06862  0.0336414  1                 910   0.542857  STEIN_ESRRA_TARGETS_RESPONSIVE_TO_ESTROGEN_DN
           10  Stein Esrra Targets Responsive To Estrogen Dn                0.733195  0.00261097    1.97665  0.056431   1                1160   0.571429  STEIN_ESRRA_TARGETS_RESPONSIVE_TO_ESTROGEN_DN
            1  Wang Lsd1 Targets Up                                        -0.60776   0.0554989    -1.55436  0.202734   1               10384   0.352941  WANG_LSD1_TARGETS_UP
            1  Wang Lsd1 Targets Dn                                        -0.561044  0.0466706    -1.47965  0.242366   1                9001   0.52381   WANG_LSD1_TARGETS_DN
            1  Ji Carcinogenesis By Kras And Stk11 Dn                      -0.746069  0.0674859    -1.55347  0.202734   1                9841   0.714286  JI_CARCINOGENESIS_BY_KRAS_AND_STK11_DN
            1  Kyng Dna Damage By 4Nqo Or Gamma Radiation                  -0.627095  0.00987505   -1.73305  0.141971   1                9864   0.384615  KYNG_DNA_DAMAGE_BY_4NQO_OR_GAMMA_RADIATION
            1  Kondo Prostate Cancer With H3K27Me3                         -0.527362  0.0507797    -1.54852  0.202971   1                7929   0.567568  KONDO_PROSTATE_CANCER_WITH_H3K27ME3
            1  Kondo Ezh2 Targets                                          -0.458023  0.0147352    -1.68011  0.157308   1                8562   0.429412  KONDO_EZH2_TARGETS
            1  Claus Pgr Positive Meningioma Dn                            -0.63501   0.0958849    -1.48012  0.242366   1               10552   0.416667  CLAUS_PGR_POSITIVE_MENINGIOMA_DN
            1  Wong Mitochondria Gene Module                                0.50624   0.0981522     1.52538  0.186438   1                3156   0.603015  WONG_MITOCHONDRIA_GENE_MODULE
            4  Wong Mitochondria Gene Module                               -0.56016   0.0413834    -1.68509  0.201047   1                8967   0.552764  WONG_MITOCHONDRIA_GENE_MODULE
            1  Wong Proteasome Gene Module                                  0.612814  0.00988701    1.77868  0.0483293  1                1094   0.413043  WONG_PROTEASOME_GENE_MODULE
            4  Wong Proteasome Gene Module                                 -0.588047  0.023344     -1.70689  0.192508   1                9767   0.456522  WONG_PROTEASOME_GENE_MODULE
            1  Amundson Gamma Radiation Response                            0.906139  0.000406421   1.82035  0.0363179  1                 206   0.757576  AMUNDSON_GAMMA_RADIATION_RESPONSE
            4  Amundson Gamma Radiation Response                           -0.869473  0.00296384   -1.74299  0.192508   1               10154   0.818182  AMUNDSON_GAMMA_RADIATION_RESPONSE
            5  Amundson Gamma Radiation Response                            0.846429  0.010917      1.69922  0.14091    1                 641   0.787879  AMUNDSON_GAMMA_RADIATION_RESPONSE
           10  Amundson Gamma Radiation Response                            0.85959   0.00482412    1.7394   0.113573   1                 921   0.818182  AMUNDSON_GAMMA_RADIATION_RESPONSE
            1  Amundson Poor Survival After Gamma Radiation 8G             -0.435571  0.0292234    -1.59307  0.18512    1                8576   0.384615  AMUNDSON_POOR_SURVIVAL_AFTER_GAMMA_RADIATION_8G
            1  Lee Metastasis And Rna Processing Up                         0.752898  0.00160901    1.86495  0.027325   1                1934   0.764706  LEE_METASTASIS_AND_RNA_PROCESSING_UP
            5  Lee Metastasis And Rna Processing Up                         0.651322  0.0326481     1.61196  0.200564   1                3015   0.705882  LEE_METASTASIS_AND_RNA_PROCESSING_UP
            1  Finetti Breast Cancer Kinome Red                             0.984247  0.000202799   1.70321  0.0785025  0.539647          189   1         FINETTI_BREAST_CANCER_KINOME_RED
            5  Finetti Breast Cancer Kinome Red                             0.971825  0.000394789   1.67723  0.153251   1                 327   1         FINETTI_BREAST_CANCER_KINOME_RED
           10  Finetti Breast Cancer Kinome Red                             0.957242  0.00140562    1.66286  0.162825   1                 489   1         FINETTI_BREAST_CANCER_KINOME_RED
            1  Finetti Breast Cancers Kinome Blue                          -0.670956  0.0334587    -1.61171  0.176986   1                8844   0.684211  FINETTI_BREAST_CANCERS_KINOME_BLUE
            1  Wallace Prostate Cancer Up                                   0.535606  0.0667056     1.49227  0.213453   1                2157   0.5       WALLACE_PROSTATE_CANCER_UP
            1  Sarrio Epithelial Mesenchymal Transition Up                  0.739178  0.000809061   1.9999   0.0129733  1                 940   0.58209   SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_UP
            4  Sarrio Epithelial Mesenchymal Transition Up                 -0.688986  0.00887224   -1.86663  0.162281   1                9701   0.574627  SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_UP
            5  Sarrio Epithelial Mesenchymal Transition Up                  0.711895  0.00378335    1.92079  0.0475578  1                1335   0.574627  SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_UP
           10  Sarrio Epithelial Mesenchymal Transition Up                  0.678489  0.0106276     1.84373  0.0674455  1                 805   0.492537  SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_UP
            5  Ouillette Cll 13Q14 Deletion Up                              0.464246  0.0182556     1.61387  0.199637   1                3281   0.490909  OUILLETTE_CLL_13Q14_DELETION_UP
            1  Yokoe Cancer Testis Antigens                                 0.51047   0.0363056     1.553    0.167527   1                1165   0.333333  YOKOE_CANCER_TESTIS_ANTIGENS
            1  Iwanaga Carcinogenesis By Kras Pten Dn                      -0.385134  0.0168051    -1.62341  0.170938   1                8702   0.368     IWANAGA_CARCINOGENESIS_BY_KRAS_PTEN_DN
            1  Iwanaga Carcinogenesis By Kras Dn                           -0.390474  0.0124899    -1.56258  0.199134   1                8922   0.356322  IWANAGA_CARCINOGENESIS_BY_KRAS_DN
            1  Fujii Ybx1 Targets Dn                                        0.671749  0.00101317    1.99355  0.0129733  1                1556   0.560976  FUJII_YBX1_TARGETS_DN
            4  Fujii Ybx1 Targets Dn                                       -0.563316  0.0542789    -1.66359  0.214503   1                9274   0.487805  FUJII_YBX1_TARGETS_DN
            5  Fujii Ybx1 Targets Dn                                        0.647075  0.00651273    1.91618  0.0475945  1                1592   0.530488  FUJII_YBX1_TARGETS_DN
           10  Fujii Ybx1 Targets Dn                                        0.64607   0.00504032    1.91292  0.056431   1                1404   0.506098  FUJII_YBX1_TARGETS_DN
            1  Osada Ascl1 Targets Dn                                      -0.622269  0.0444086    -1.58882  0.187174   1                8502   0.75      OSADA_ASCL1_TARGETS_DN
            1  Dairkee Cancer Prone Response Bpa                            0.555324  0.0048048     1.78962  0.0453485  1                1236   0.414634  DAIRKEE_CANCER_PRONE_RESPONSE_BPA
           10  Dairkee Cancer Prone Response Bpa                            0.527907  0.0143165     1.70904  0.12849    1                1988   0.463415  DAIRKEE_CANCER_PRONE_RESPONSE_BPA
            1  Dairkee Cancer Prone Response Bpa E2                         0.52919   0.00256006    1.86116  0.0278677  1                2055   0.483516  DAIRKEE_CANCER_PRONE_RESPONSE_BPA_E2
            5  Dairkee Cancer Prone Response Bpa E2                         0.445449  0.0410739     1.56344  0.246162   1                2448   0.362637  DAIRKEE_CANCER_PRONE_RESPONSE_BPA_E2
            1  Riggi Ewing Sarcoma Progenitor Dn                           -0.505112  0.0694473    -1.58998  0.187154   1                8881   0.632353  RIGGI_EWING_SARCOMA_PROGENITOR_DN
            1  Engelmann Cancer Progenitors Up                             -0.439053  0.0631127    -1.4821   0.241988   1                8142   0.525     ENGELMANN_CANCER_PROGENITORS_UP
            1  Engelmann Cancer Progenitors Dn                             -0.504895  0.0181443    -1.64972  0.168175   1                8504   0.511111  ENGELMANN_CANCER_PROGENITORS_DN
            1  Huang Foxa2 Targets Dn                                      -0.546304  0.10554      -1.47032  0.24821    1                8182   0.533333  HUANG_FOXA2_TARGETS_DN
            1  Mishra Carcinoma Associated Fibroblast Up                   -0.774087  0.00603986   -1.82005  0.119811   1                9342   0.75      MISHRA_CARCINOMA_ASSOCIATED_FIBROBLAST_UP
            1  Mccabe Hoxc6 Targets Dn                                     -0.66758   0.0261097    -1.63796  0.170938   1                9116   0.5       MCCABE_HOXC6_TARGETS_DN
            1  Ting Silenced By Dicer                                      -0.625762  0.0106898    -1.76263  0.139375   1                9580   0.44      TING_SILENCED_BY_DICER
            1  Molenaar Targets Of Ccnd1 And Cdk4 Up                       -0.518172  0.00624245   -1.73541  0.141912   1                7690   0.555556  MOLENAAR_TARGETS_OF_CCND1_AND_CDK4_UP
            1  Molenaar Targets Of Ccnd1 And Cdk4 Dn                        0.752872  0.00749291    1.82821  0.0346151  1                 711   0.575     MOLENAAR_TARGETS_OF_CCND1_AND_CDK4_DN
            4  Molenaar Targets Of Ccnd1 And Cdk4 Dn                       -0.676861  0.0576124    -1.63146  0.230376   1               10054   0.525     MOLENAAR_TARGETS_OF_CCND1_AND_CDK4_DN
            5  Molenaar Targets Of Ccnd1 And Cdk4 Dn                        0.79862   0.00118906    1.91412  0.0477898  1                 624   0.6       MOLENAAR_TARGETS_OF_CCND1_AND_CDK4_DN
           10  Molenaar Targets Of Ccnd1 And Cdk4 Dn                        0.738995  0.0156752     1.79415  0.0824479  1                 808   0.6       MOLENAAR_TARGETS_OF_CCND1_AND_CDK4_DN
            1  Acevedo Liver Tumor Vs Normal Adjacent Tissue Dn            -0.386648  0.00556439   -1.66303  0.162484   1                8966   0.333333  ACEVEDO_LIVER_TUMOR_VS_NORMAL_ADJACENT_TISSUE_DN
            1  Mitsiades Response To Aplidin Dn                             0.610813  0.00528993    1.88243  0.0235955  1                2003   0.530702  MITSIADES_RESPONSE_TO_APLIDIN_DN
            5  Mitsiades Response To Aplidin Dn                             0.639217  0.00219868    1.95899  0.0409356  1                1471   0.47807   MITSIADES_RESPONSE_TO_APLIDIN_DN
           10  Mitsiades Response To Aplidin Dn                             0.595015  0.0120192     1.83422  0.0688673  1                2173   0.513158  MITSIADES_RESPONSE_TO_APLIDIN_DN
            1  Hoque Methylated In Cancer                                  -0.521739  0.0141237    -1.70382  0.147212   1                9299   0.416667  HOQUE_METHYLATED_IN_CANCER
            1  Smid Breast Cancer Relapse In Brain Dn                      -0.61672   0.0837238    -1.57645  0.194122   1                8197   0.675676  SMID_BREAST_CANCER_RELAPSE_IN_BRAIN_DN
            1  Smid Breast Cancer Relapse In Bone Up                       -0.579838  0.0679631    -1.59581  0.18268    1                8065   0.5875    SMID_BREAST_CANCER_RELAPSE_IN_BONE_UP
            1  Smid Breast Cancer Relapse In Lung Dn                       -0.718395  0.0129319    -1.74093  0.140775   1                9205   0.7       SMID_BREAST_CANCER_RELAPSE_IN_LUNG_DN
            1  Smid Breast Cancer Relapse In Pleura Dn                      0.702564  0.0241086     1.62581  0.115119   1                1365   0.470588  SMID_BREAST_CANCER_RELAPSE_IN_PLEURA_DN
           10  Smid Breast Cancer Relapse In Pleura Dn                      0.721156  0.0137807     1.67199  0.157141   1                1331   0.470588  SMID_BREAST_CANCER_RELAPSE_IN_PLEURA_DN
            1  Smid Breast Cancer Luminal B Up                             -0.558385  0.079112     -1.56353  0.199134   1                7790   0.620438  SMID_BREAST_CANCER_LUMINAL_B_UP
            1  Smid Breast Cancer Luminal A Up                             -0.740111  0.00142363   -1.92296  0.116729   1                9205   0.780822  SMID_BREAST_CANCER_LUMINAL_A_UP
            1  Smid Breast Cancer Luminal A Dn                              0.931869  0.000201735   1.68929  0.0843288  0.536817          769   1         SMID_BREAST_CANCER_LUMINAL_A_DN
            4  Smid Breast Cancer Luminal A Dn                             -0.94402   0.000589391  -1.70451  0.192508   1               10641   0.75      SMID_BREAST_CANCER_LUMINAL_A_DN
            5  Smid Breast Cancer Luminal A Dn                              0.878247  0.010691      1.57929  0.227597   1                 418   0.833333  SMID_BREAST_CANCER_LUMINAL_A_DN
           10  Smid Breast Cancer Luminal A Dn                              0.936549  0.000605694   1.70235  0.130517   1                 717   1         SMID_BREAST_CANCER_LUMINAL_A_DN
            1  Smid Breast Cancer Normal Like Up                           -0.59231   0.0924152    -1.59432  0.184082   1                9010   0.461017  SMID_BREAST_CANCER_NORMAL_LIKE_UP
            1  Ito Pttg1 Targets Up                                        -0.623665  0.065587     -1.51922  0.221755   1                8638   0.636364  ITO_PTTG1_TARGETS_UP
            1  Le Ski Targets Up                                           -0.696037  0.0165984    -1.70217  0.147635   1                7898   0.857143  LE_SKI_TARGETS_UP
            1  Bonome Ovarian Cancer Poor Survival Up                      -0.570552  0.0472441    -1.58235  0.190354   1                7369   0.76      BONOME_OVARIAN_CANCER_POOR_SURVIVAL_UP
            1  Bonome Ovarian Cancer Poor Survival Dn                       0.713691  0.00264389    1.85414  0.029501   1                2487   0.722222  BONOME_OVARIAN_CANCER_POOR_SURVIVAL_DN
           10  Bonome Ovarian Cancer Poor Survival Dn                       0.598523  0.0462889     1.57667  0.231959   1                3830   0.777778  BONOME_OVARIAN_CANCER_POOR_SURVIVAL_DN
            1  Mcgarvey Silenced By Methylation In Colon Cancer            -0.544612  0.0772193    -1.48318  0.241988   1                9086   0.473684  MCGARVEY_SILENCED_BY_METHYLATION_IN_COLON_CANCER
            1  Basso Hairy Cell Leukemia Dn                                -0.595574  0.00542169   -1.86832  0.119811   1                8397   0.58209   BASSO_HAIRY_CELL_LEUKEMIA_DN
            1  Izadpanah Stem Cell Adipose Vs Bone Dn                      -0.593399  0.0104987    -1.84181  0.119811   1                8764   0.549451  IZADPANAH_STEM_CELL_ADIPOSE_VS_BONE_DN
            1  Bernard Ppapdc1B Targets Dn                                 -0.480854  0.0151822    -1.64118  0.170938   1                8493   0.413043  BERNARD_PPAPDC1B_TARGETS_DN
            1  Zheng Glioblastoma Plasticity Up                             0.56896   0.00121778    1.96018  0.0147521  1                 504   0.34104   ZHENG_GLIOBLASTOMA_PLASTICITY_UP
            5  Zheng Glioblastoma Plasticity Up                             0.5062    0.0193227     1.74118  0.117295   1                1096   0.352601  ZHENG_GLIOBLASTOMA_PLASTICITY_UP
           10  Zheng Glioblastoma Plasticity Up                             0.556446  0.00316832    1.92842  0.056431   1                1049   0.398844  ZHENG_GLIOBLASTOMA_PLASTICITY_UP
            1  Zheng Glioblastoma Plasticity Dn                            -0.609131  0.00944344   -1.7973   0.124176   1                9631   0.431818  ZHENG_GLIOBLASTOMA_PLASTICITY_DN
            8  Buckanovich T Lymphocyte Homing On Tumor Dn                  0.6456    0.000593942   2.01159  0.130059   1                2919   0.65      BUCKANOVICH_T_LYMPHOCYTE_HOMING_ON_TUMOR_DN
            1  Zheng Il22 Signaling Up                                     -0.515479  0.0576148    -1.52935  0.214421   1                9596   0.46875   ZHENG_IL22_SIGNALING_UP
            1  Honma Docetaxel Resistance                                   0.61082   0.0555556     1.59477  0.139853   1                2139   0.617647  HONMA_DOCETAXEL_RESISTANCE
            4  Honma Docetaxel Resistance                                  -0.698491  0.00639872   -1.8229   0.187832   1                8335   0.794118  HONMA_DOCETAXEL_RESISTANCE
            1  Matzuk Ovulation                                            -0.638923  0.0202582    -1.65664  0.16369    1                9299   0.5       MATZUK_OVULATION
            1  Matzuk Meiotic And Dna Repair                                0.601229  0.0156218     1.70296  0.0785025  1                1048   0.4       MATZUK_MEIOTIC_AND_DNA_REPAIR
           10  Matzuk Meiotic And Dna Repair                                0.715927  0.00040282    2.02782  0.056431   1                 959   0.48      MATZUK_MEIOTIC_AND_DNA_REPAIR
            1  Matzuk Central For Female Fertility                         -0.607689  0.0584455    -1.53794  0.210521   1               10174   0.583333  MATZUK_CENTRAL_FOR_FEMALE_FERTILITY
            1  Matzuk Embryonic Germ Cell                                  -0.71758   0.00402253   -1.74718  0.139914   1                8229   0.75      MATZUK_EMBRYONIC_GERM_CELL
            1  Matzuk Spermatocyte                                          0.448293  0.0336202     1.56305  0.160717   1                1198   0.295455  MATZUK_SPERMATOCYTE
            5  Matzuk Spermatocyte                                          0.470291  0.0156156     1.62764  0.19014    1                 535   0.227273  MATZUK_SPERMATOCYTE
           10  Matzuk Spermatocyte                                          0.46192   0.0219692     1.60189  0.207217   1                1107   0.295455  MATZUK_SPERMATOCYTE
            1  Chung Blister Cytotoxicity Up                                0.413601  0.0545052     1.50608  0.204051   1                1340   0.275229  CHUNG_BLISTER_CYTOTOXICITY_UP
            1  Whiteford Pediatric Cancer Markers                           0.889775  0.00020145    1.98682  0.0134874  0.53606           646   0.793103  WHITEFORD_PEDIATRIC_CANCER_MARKERS
            4  Whiteford Pediatric Cancer Markers                          -0.796356  0.0088426    -1.77576  0.192508   1               10104   0.735632  WHITEFORD_PEDIATRIC_CANCER_MARKERS
            5  Whiteford Pediatric Cancer Markers                           0.827429  0.00258501    1.83868  0.0694466  1                1006   0.747126  WHITEFORD_PEDIATRIC_CANCER_MARKERS
           10  Whiteford Pediatric Cancer Markers                           0.837092  0.00120555    1.86908  0.0627622  1                 805   0.758621  WHITEFORD_PEDIATRIC_CANCER_MARKERS
            1  Grade Metastasis Dn                                          0.644127  0.00899281    1.78799  0.0458337  1                3197   0.804878  GRADE_METASTASIS_DN
            4  Grade Metastasis Dn                                         -0.711317  0.00039604   -1.94776  0.155097   1                8531   0.731707  GRADE_METASTASIS_DN
            1  Grade Colon And Rectal Cancer Up                             0.611759  0.000201572   2.03732  0.0127991  0.536384         1671   0.503906  GRADE_COLON_AND_RECTAL_CANCER_UP
            4  Grade Colon And Rectal Cancer Up                            -0.581052  0.00138203   -1.93801  0.155097   1                8872   0.535156  GRADE_COLON_AND_RECTAL_CANCER_UP
            5  Grade Colon And Rectal Cancer Up                             0.502184  0.0354652     1.68182  0.151662   1                2105   0.441406  GRADE_COLON_AND_RECTAL_CANCER_UP
           10  Grade Colon And Rectal Cancer Up                             0.548695  0.00598325    1.82871  0.0688673  1                2163   0.492188  GRADE_COLON_AND_RECTAL_CANCER_UP
            1  Grade Colon And Rectal Cancer Dn                            -0.484634  0.00667611   -1.72312  0.146665   1                8897   0.453125  GRADE_COLON_AND_RECTAL_CANCER_DN
            1  Boquest Stem Cell Up                                        -0.694437  0.0105928    -1.8728   0.119811   1                9213   0.660194  BOQUEST_STEM_CELL_UP
            1  Boquest Stem Cell Dn                                        -0.509962  0.043879     -1.68617  0.154958   1                8713   0.450292  BOQUEST_STEM_CELL_DN
            1  Boquest Stem Cell Cultured Vs Fresh Dn                      -0.756819  0.0038736    -1.85046  0.119811   1                9051   0.785714  BOQUEST_STEM_CELL_CULTURED_VS_FRESH_DN
            1  Labbe Wnt3A Targets Up                                       0.410898  0.054883      1.48476  0.217372   1                1253   0.3       LABBE_WNT3A_TARGETS_UP
            1  Labbe Tgfb1 Targets Up                                      -0.486975  0.0538261    -1.57609  0.194122   1                8685   0.477612  LABBE_TGFB1_TARGETS_UP
            1  Nakamura Metastasis Model Dn                                -0.593752  0.00343782   -1.83403  0.119811   1                8888   0.535714  NAKAMURA_METASTASIS_MODEL_DN
            1  Gratias Retinoblastoma 16Q24                                 0.756223  0.0550532     1.57508  0.155023   1                1547   0.647059  GRATIAS_RETINOBLASTOMA_16Q24
           10  Gratias Retinoblastoma 16Q24                                 0.832985  0.00993246    1.73715  0.113573   1                1263   0.823529  GRATIAS_RETINOBLASTOMA_16Q24
            1  West Adrenocortical Tumor Up                                 0.579314  0.000402739   2.0178   0.0127991  1                1581   0.438462  WEST_ADRENOCORTICAL_TUMOR_UP
            4  West Adrenocortical Tumor Up                                -0.541456  0.00378486   -1.87027  0.162281   1                9203   0.473077  WEST_ADRENOCORTICAL_TUMOR_UP
            5  West Adrenocortical Tumor Up                                 0.499946  0.0212508     1.72889  0.124046   1                2240   0.419231  WEST_ADRENOCORTICAL_TUMOR_UP
           10  West Adrenocortical Tumor Up                                 0.48489   0.0277557     1.69273  0.136131   1                1423   0.330769  WEST_ADRENOCORTICAL_TUMOR_UP
            1  West Adrenocortical Carcinoma Vs Adenoma Up                  0.659243  0.0169287     1.66299  0.0974565  1                2512   0.666667  WEST_ADRENOCORTICAL_CARCINOMA_VS_ADENOMA_UP
            4  West Adrenocortical Carcinoma Vs Adenoma Up                 -0.674099  0.0110628    -1.70969  0.192508   1                9565   0.583333  WEST_ADRENOCORTICAL_CARCINOMA_VS_ADENOMA_UP
            1  Ray Tumorigenesis By Erbb2 Cdc25A Dn                        -0.391969  0.0282194    -1.54655  0.203357   1                9024   0.309524  RAY_TUMORIGENESIS_BY_ERBB2_CDC25A_DN
            1  Podar Response To Adaphostin Dn                              0.713514  0.00983936    1.76319  0.0535254  1                2416   0.8       PODAR_RESPONSE_TO_ADAPHOSTIN_DN
            4  Podar Response To Adaphostin Dn                             -0.709279  0.0113253    -1.74474  0.192508   1                9413   0.666667  PODAR_RESPONSE_TO_ADAPHOSTIN_DN
            5  Podar Response To Adaphostin Dn                              0.662151  0.0300318     1.6363   0.184221   1                1185   0.533333  PODAR_RESPONSE_TO_ADAPHOSTIN_DN
            1  Lin Tumor Escape From Immune Attack                         -0.658751  0.0876792    -1.4699   0.248243   1                9429   0.6       LIN_TUMOR_ESCAPE_FROM_IMMUNE_ATTACK
            4  Chen Hoxa5 Targets 9Hr Dn                                   -0.515726  0.00930141   -1.70365  0.192508   1                9231   0.5       CHEN_HOXA5_TARGETS_9HR_DN
            1  Cadwell Atg16L1 Targets Up                                  -0.493566  0.0242826    -1.63829  0.170938   1                8862   0.434783  CADWELL_ATG16L1_TARGETS_UP
            1  Bredemeyer Rag Signaling Via Atm Not Via Nfkb Up            -0.424165  0.0469812    -1.47593  0.244495   1                8312   0.432432  BREDEMEYER_RAG_SIGNALING_VIA_ATM_NOT_VIA_NFKB_UP
            1  Thum Mir21 Targets Heart Disease Up                         -0.914826  0.00121926   -1.74696  0.139914   1               10229   0.933333  THUM_MIR21_TARGETS_HEART_DISEASE_UP
           10  Taylor Methylated In Acute Lymphoblastic Leukemia            0.441607  0.0249551     1.59716  0.211912   1                1316   0.303571  TAYLOR_METHYLATED_IN_ACUTE_LYMPHOBLASTIC_LEUKEMIA
            1  Chng Multiple Myeloma Hyperploid Dn                          0.508213  0.0926966     1.46342  0.229649   1                3326   0.62963   CHNG_MULTIPLE_MYELOMA_HYPERPLOID_DN
            1  Huper Breast Basal Vs Luminal Dn                            -0.47907   0.0334967    -1.59743  0.181764   1                8324   0.52      HUPER_BREAST_BASAL_VS_LUMINAL_DN
            1  Maloney Response To 17Aag Dn                                 0.554825  0.036933      1.66673  0.0958915  1                1303   0.391892  MALONEY_RESPONSE_TO_17AAG_DN
            4  Maloney Response To 17Aag Dn                                -0.555169  0.0319842    -1.66931  0.212004   1                9498   0.445946  MALONEY_RESPONSE_TO_17AAG_DN
            1  Cui Glucose Deprivation                                      0.504469  0.0670479     1.5433   0.174459   1                2625   0.511111  CUI_GLUCOSE_DEPRIVATION
            1  Blum Response To Salirasib Dn                                0.622386  0.000808244   2.06383  0.0127991  1                1343   0.494983  BLUM_RESPONSE_TO_SALIRASIB_DN
            4  Blum Response To Salirasib Dn                               -0.499811  0.0583021    -1.65225  0.220221   1                9568   0.391304  BLUM_RESPONSE_TO_SALIRASIB_DN
            5  Blum Response To Salirasib Dn                                0.612875  0.00159363    2.02084  0.0336414  1                1644   0.491639  BLUM_RESPONSE_TO_SALIRASIB_DN
           10  Blum Response To Salirasib Dn                                0.583766  0.00619257    1.94031  0.056431   1                1423   0.458194  BLUM_RESPONSE_TO_SALIRASIB_DN
            1  Mueller Plurinet                                             0.636189  0.000598563   1.99749  0.0129733  1                1919   0.547718  MUELLER_PLURINET
            4  Mueller Plurinet                                            -0.520391  0.053357     -1.61753  0.237982   1                8919   0.506224  MUELLER_PLURINET
            5  Mueller Plurinet                                             0.572083  0.0136876     1.78697  0.0965026  1                1938   0.448133  MUELLER_PLURINET
           10  Mueller Plurinet                                             0.579078  0.00821643    1.8142   0.0752823  1                2995   0.609959  MUELLER_PLURINET
            5  Firestein Ctnnb1 Pathway                                     0.558671  0.00841009    1.70082  0.139959   1                2530   0.384615  FIRESTEIN_CTNNB1_PATHWAY
            1  Jiang Tip30 Targets Up                                      -0.476004  0.0466585    -1.55561  0.202734   1                9935   0.375     JIANG_TIP30_TARGETS_UP
            1  Wang Neoplastic Transformation By Ccnd1 Myc                 -0.591268  0.086515     -1.4955   0.232098   1                7831   0.785714  WANG_NEOPLASTIC_TRANSFORMATION_BY_CCND1_MYC
            1  Tsai Response To Radiation Therapy                          -0.641856  0.0247186    -1.7101   0.146665   1                9589   0.571429  TSAI_RESPONSE_TO_RADIATION_THERAPY
            1  Beier Glioma Stem Cell Dn                                   -0.461788  0.0287712    -1.60652  0.180997   1                8645   0.388889  BEIER_GLIOMA_STEM_CELL_DN
            1  Vart Kshv Infection Angiogenic Markers Up                   -0.586204  0.026252     -1.74873  0.139914   1                8861   0.566372  VART_KSHV_INFECTION_ANGIOGENIC_MARKERS_UP
            1  Vart Kshv Infection Angiogenic Markers Dn                   -0.519483  0.0344898    -1.66365  0.162484   1                8480   0.525     VART_KSHV_INFECTION_ANGIOGENIC_MARKERS_DN
            1  Mootha Tca                                                   0.664236  0.0256514     1.65856  0.0991974  1                2294   0.625     MOOTHA_TCA
            4  Mootha Tca                                                  -0.698817  0.0122122    -1.73814  0.192508   1                9514   0.625     MOOTHA_TCA
           10  Mootha Tca                                                   0.651271  0.0294        1.64195  0.174876   1                2243   0.6875    MOOTHA_TCA
            1  Lee Liver Cancer Survival Dn                                 0.652692  0.000604961   2.01531  0.0127991  1                1601   0.49359   LEE_LIVER_CANCER_SURVIVAL_DN
            4  Lee Liver Cancer Survival Dn                                -0.598346  0.00672335   -1.85607  0.173896   1                9223   0.519231  LEE_LIVER_CANCER_SURVIVAL_DN
            5  Lee Liver Cancer Survival Dn                                 0.570035  0.0142829     1.77777  0.0974471  1                2759   0.564103  LEE_LIVER_CANCER_SURVIVAL_DN
           10  Lee Liver Cancer Survival Dn                                 0.515733  0.0596865     1.59327  0.215799   1                3046   0.551282  LEE_LIVER_CANCER_SURVIVAL_DN
            8  Lee Liver Cancer Survival Up                                 0.54023   0.000198965   2.04195  0.104006   0.529447         2175   0.454545  LEE_LIVER_CANCER_SURVIVAL_UP
            1  Boylan Multiple Myeloma C Cluster Up                         0.529639  0.0237138     1.63489  0.110528   1                1718   0.518519  BOYLAN_MULTIPLE_MYELOMA_C_CLUSTER_UP
            5  Boylan Multiple Myeloma C Cluster Up                         0.55781   0.0142285     1.72006  0.12811    1                1053   0.37037   BOYLAN_MULTIPLE_MYELOMA_C_CLUSTER_UP
            1  Boylan Multiple Myeloma C Cluster Dn                        -0.506165  0.0404081    -1.57258  0.195207   1                9060   0.423077  BOYLAN_MULTIPLE_MYELOMA_C_CLUSTER_DN
            1  Boylan Multiple Myeloma C D Up                               0.368224  0.0539176     1.44667  0.238115   1                2045   0.37      BOYLAN_MULTIPLE_MYELOMA_C_D_UP
            5  Boylan Multiple Myeloma C D Up                               0.406198  0.0161226     1.60146  0.209048   1                2497   0.44      BOYLAN_MULTIPLE_MYELOMA_C_D_UP
            1  Wu Silenced By Methylation In Bladder Cancer                -0.678161  0.00263425   -1.91167  0.116729   1                8668   0.645833  WU_SILENCED_BY_METHYLATION_IN_BLADDER_CANCER
            1  Hellebrekers Silenced During Tumor Angiogenesis             -0.541003  0.0372319    -1.65684  0.16369    1                9012   0.444444  HELLEBREKERS_SILENCED_DURING_TUMOR_ANGIOGENESIS
            1  Gu Pdef Targets Up                                          -0.584112  0.0403877    -1.67259  0.161802   1                9070   0.546875  GU_PDEF_TARGETS_UP
            1  Bohn Primary Immunodeficiency Syndrom Up                     0.5458    0.0281068     1.64137  0.107415   1                1404   0.380952  BOHN_PRIMARY_IMMUNODEFICIENCY_SYNDROM_UP
            5  Bohn Primary Immunodeficiency Syndrom Up                     0.548528  0.028133      1.64737  0.178715   1                3311   0.595238  BOHN_PRIMARY_IMMUNODEFICIENCY_SYNDROM_UP
            1  Flotho Pediatric All Therapy Response Dn                     0.501702  0.0574249     1.50058  0.209258   1                 949   0.347826  FLOTHO_PEDIATRIC_ALL_THERAPY_RESPONSE_DN
           10  Flotho Pediatric All Therapy Response Dn                     0.575931  0.015437      1.7066   0.129418   1                1008   0.347826  FLOTHO_PEDIATRIC_ALL_THERAPY_RESPONSE_DN
            1  Zaidi Osteoblast Transcription Factors                      -0.738749  0.00403796   -1.79719  0.124176   1               10354   0.6       ZAIDI_OSTEOBLAST_TRANSCRIPTION_FACTORS
            1  Seki Inflammatory Response Lps Up                           -0.529866  0.102343     -1.50968  0.22708    1                8756   0.457627  SEKI_INFLAMMATORY_RESPONSE_LPS_UP
            1  Seki Inflammatory Response Lps Dn                           -0.641391  0.0242067    -1.63381  0.170938   1                8205   0.571429  SEKI_INFLAMMATORY_RESPONSE_LPS_DN
            1  Chauhan Response To Methoxyestradiol Up                      0.56573   0.0164493     1.71009  0.0761105  1                1480   0.534884  CHAUHAN_RESPONSE_TO_METHOXYESTRADIOL_UP
            1  Frasor Response To Serm Or Fulvestrant Dn                    0.855319  0.000401284   1.99536  0.0129733  1                 814   0.75      FRASOR_RESPONSE_TO_SERM_OR_FULVESTRANT_DN
            5  Frasor Response To Serm Or Fulvestrant Dn                    0.799002  0.00336367    1.85319  0.0639129  1                1006   0.704545  FRASOR_RESPONSE_TO_SERM_OR_FULVESTRANT_DN
           10  Frasor Response To Serm Or Fulvestrant Dn                    0.826999  0.000401768   1.92199  0.056431   1                1423   0.795455  FRASOR_RESPONSE_TO_SERM_OR_FULVESTRANT_DN
            1  Goldrath Antigen Response                                    0.548056  0.004828      1.88045  0.0239353  1                 924   0.330882  GOLDRATH_ANTIGEN_RESPONSE
            5  Goldrath Antigen Response                                    0.501115  0.023103      1.70463  0.136369   1                1484   0.360294  GOLDRATH_ANTIGEN_RESPONSE
           10  Goldrath Antigen Response                                    0.517177  0.0130785     1.76338  0.101982   1                1459   0.356618  GOLDRATH_ANTIGEN_RESPONSE
            1  Ichiba Graft Versus Host Disease 35D Dn                     -0.47002   0.048345     -1.49119  0.235557   1                9124   0.333333  ICHIBA_GRAFT_VERSUS_HOST_DISEASE_35D_DN
            1  Lindstedt Dendritic Cell Maturation C                       -0.491811  0.0255585    -1.66114  0.162756   1                8844   0.466667  LINDSTEDT_DENDRITIC_CELL_MATURATION_C
            1  Ruiz Tnc Targets Up                                         -0.542436  0.00141271   -1.90526  0.116729   1                8293   0.524823  RUIZ_TNC_TARGETS_UP
            1  Ruiz Tnc Targets Dn                                          0.6262    0.00402495    1.89794  0.0214742  1                1255   0.5       RUIZ_TNC_TARGETS_DN
            5  Ruiz Tnc Targets Dn                                          0.598827  0.0140873     1.80389  0.0875949  1                1551   0.490909  RUIZ_TNC_TARGETS_DN
           10  Ruiz Tnc Targets Dn                                          0.610648  0.0121756     1.84205  0.0677144  1                 714   0.418182  RUIZ_TNC_TARGETS_DN
            1  Schraets Mll Targets Dn                                     -0.664023  0.012087     -1.78474  0.131106   1                8245   0.666667  SCHRAETS_MLL_TARGETS_DN
            1  Chang Core Serum Response Up                                 0.61585   0.00100462    1.93944  0.0167144  1                1668   0.546448  CHANG_CORE_SERUM_RESPONSE_UP
            4  Chang Core Serum Response Up                                -0.605824  0.00118087   -1.90509  0.155097   1                8834   0.601093  CHANG_CORE_SERUM_RESPONSE_UP
           10  Chang Core Serum Response Up                                 0.550888  0.0151848     1.73687  0.113573   1                2488   0.535519  CHANG_CORE_SERUM_RESPONSE_UP
            1  Chang Core Serum Response Dn                                -0.506721  0.000606428  -1.92285  0.116729   1                8358   0.52551   CHANG_CORE_SERUM_RESPONSE_DN
            1  Chang Cycling Genes                                          0.812694  0.000405927   2.0362   0.0127991  1                 504   0.666667  CHANG_CYCLING_GENES
            4  Chang Cycling Genes                                         -0.68933   0.0342857    -1.72772  0.192508   1               10178   0.587719  CHANG_CYCLING_GENES
            5  Chang Cycling Genes                                          0.776515  0.00179569    1.93197  0.0454213  1                1031   0.675439  CHANG_CYCLING_GENES
           10  Chang Cycling Genes                                          0.777212  0.00161714    1.94761  0.056431   1                1223   0.701754  CHANG_CYCLING_GENES
            1  Vantveer Breast Cancer Esr1 Up                              -0.644271  0.0489623    -1.66864  0.162484   1                8213   0.740506  VANTVEER_BREAST_CANCER_ESR1_UP
            1  Vantveer Breast Cancer Esr1 Dn                               0.631772  0.083909      1.58289  0.150484   1                2602   0.575269  VANTVEER_BREAST_CANCER_ESR1_DN
           10  Vantveer Breast Cancer Esr1 Dn                               0.661934  0.0452854     1.6316   0.181196   1                2403   0.629032  VANTVEER_BREAST_CANCER_ESR1_DN
            1  Vantveer Breast Cancer Brca1 Up                              0.623451  0.0282787     1.64455  0.106001   1                1738   0.5       VANTVEER_BREAST_CANCER_BRCA1_UP
           10  Vantveer Breast Cancer Brca1 Up                              0.645569  0.0140091     1.71693  0.122525   1                2886   0.733333  VANTVEER_BREAST_CANCER_BRCA1_UP
            1  Vantveer Breast Cancer Brca1 Dn                             -0.582648  0.0203666    -1.6987   0.148746   1                9086   0.526316  VANTVEER_BREAST_CANCER_BRCA1_DN
            1  Wallace Jak2 Targets Up                                      0.660208  0.00833503    1.76558  0.0530966  1                2752   0.7       WALLACE_JAK2_TARGETS_UP
           10  Wallace Jak2 Targets Up                                      0.697142  0.00180542    1.87053  0.0627622  1                2052   0.6       WALLACE_JAK2_TARGETS_UP
            1  Hoffmann Pre Bi To Large Pre Bii Lymphocyte Up               0.454605  0.0258325     1.54239  0.174598   1                1356   0.333333  HOFFMANN_PRE_BI_TO_LARGE_PRE_BII_LYMPHOCYTE_UP
           10  Hoffmann Pre Bi To Large Pre Bii Lymphocyte Up               0.515135  0.00321802    1.75551  0.107034   1                2212   0.444444  HOFFMANN_PRE_BI_TO_LARGE_PRE_BII_LYMPHOCYTE_UP
            1  Hoffmann Large To Small Pre Bii Lymphocyte Up                0.765967  0.000405597   2.06899  0.0127991  1                1210   0.654412  HOFFMANN_LARGE_TO_SMALL_PRE_BII_LYMPHOCYTE_UP
            4  Hoffmann Large To Small Pre Bii Lymphocyte Up               -0.607274  0.0705836    -1.63755  0.229428   1                9472   0.529412  HOFFMANN_LARGE_TO_SMALL_PRE_BII_LYMPHOCYTE_UP
            5  Hoffmann Large To Small Pre Bii Lymphocyte Up                0.730845  0.000998403   1.96628  0.0383868  1                1096   0.580882  HOFFMANN_LARGE_TO_SMALL_PRE_BII_LYMPHOCYTE_UP
           10  Hoffmann Large To Small Pre Bii Lymphocyte Up                0.714051  0.00381986    1.93417  0.056431   1                 887   0.566176  HOFFMANN_LARGE_TO_SMALL_PRE_BII_LYMPHOCYTE_UP
            4  Hoffmann Immature To Mature B Lymphocyte Dn                 -0.494262  0.01618      -1.65661  0.216995   1                8958   0.5       HOFFMANN_IMMATURE_TO_MATURE_B_LYMPHOCYTE_DN
            1  Lee Early T Lymphocyte Up                                    0.835186  0.000400481   2.03648  0.0127991  1                 509   0.651515  LEE_EARLY_T_LYMPHOCYTE_UP
            4  Lee Early T Lymphocyte Up                                   -0.729008  0.0292548    -1.76092  0.192508   1               10270   0.621212  LEE_EARLY_T_LYMPHOCYTE_UP
            5  Lee Early T Lymphocyte Up                                    0.803044  0.0041866     1.93789  0.0454213  1                 931   0.666667  LEE_EARLY_T_LYMPHOCYTE_UP
           10  Lee Early T Lymphocyte Up                                    0.819189  0.000992654   1.99368  0.056431   1                 622   0.666667  LEE_EARLY_T_LYMPHOCYTE_UP
            1  Lee Double Polar Thymocyte                                  -0.620634  0.0362776    -1.58571  0.189217   1               10063   0.5       LEE_DOUBLE_POLAR_THYMOCYTE
            1  Finak Breast Cancer Sdpp Signature                          -0.631237  0.0631176    -1.53845  0.210383   1                9526   0.4375    FINAK_BREAST_CANCER_SDPP_SIGNATURE
            1  Shedden Lung Cancer Good Survival A4                        -0.428028  0.0152122    -1.64247  0.170938   1                7268   0.556338  SHEDDEN_LUNG_CANCER_GOOD_SURVIVAL_A4
            1  Zhang Tlx Targets 36Hr Up                                   -0.451356  0.00162437   -1.83736  0.119811   1                8110   0.475676  ZHANG_TLX_TARGETS_36HR_UP
            1  Zhang Tlx Targets 36Hr Dn                                    0.681096  0.00348289    1.93129  0.0172116  1                1213   0.464968  ZHANG_TLX_TARGETS_36HR_DN
            5  Zhang Tlx Targets 36Hr Dn                                    0.700799  0.000998403   1.97491  0.0383868  1                1882   0.630573  ZHANG_TLX_TARGETS_36HR_DN
           10  Zhang Tlx Targets 36Hr Dn                                    0.647751  0.00826446    1.8374   0.0688673  1                 915   0.414013  ZHANG_TLX_TARGETS_36HR_DN
            1  Zhang Tlx Targets 60Hr Up                                   -0.491274  0.00162371   -1.90687  0.116729   1                8930   0.433884  ZHANG_TLX_TARGETS_60HR_UP
            1  Zhang Tlx Targets 60Hr Dn                                    0.733543  0.000406669   2.07923  0.0127991  1                1013   0.568282  ZHANG_TLX_TARGETS_60HR_DN
            5  Zhang Tlx Targets 60Hr Dn                                    0.703006  0.00139054    1.98179  0.0380769  1                1571   0.594714  ZHANG_TLX_TARGETS_60HR_DN
           10  Zhang Tlx Targets 60Hr Dn                                    0.673568  0.00525996    1.91135  0.056431   1                1500   0.563877  ZHANG_TLX_TARGETS_60HR_DN
            1  Zhang Tlx Targets Dn                                         0.849986  0.000611995   2.01948  0.0127991  1                1000   0.757143  ZHANG_TLX_TARGETS_DN
            5  Zhang Tlx Targets Dn                                         0.833512  0.000991867   1.96712  0.0383868  1                 753   0.728571  ZHANG_TLX_TARGETS_DN
           10  Zhang Tlx Targets Dn                                         0.82618   0.000403633   1.96656  0.056431   1                 819   0.685714  ZHANG_TLX_TARGETS_DN
            1  Zhang Tlx Targets Up                                        -0.550445  0.00347719   -1.87882  0.119811   1                8573   0.541176  ZHANG_TLX_TARGETS_UP
            1  Mellman Tut1 Targets Up                                      0.51477   0.0761656     1.46127  0.23075    1                1795   0.444444  MELLMAN_TUT1_TARGETS_UP
            1  Shaffer Irf4 Targets In Activated B Lymphocyte               0.493281  0.0216786     1.68864  0.0843864  1                2421   0.540541  SHAFFER_IRF4_TARGETS_IN_ACTIVATED_B_LYMPHOCYTE
            4  Shaffer Irf4 Targets In Activated B Lymphocyte              -0.483237  0.025        -1.64784  0.222238   1                8947   0.445946  SHAFFER_IRF4_TARGETS_IN_ACTIVATED_B_LYMPHOCYTE
            1  Ferrando Hox11 Neighbors                                     0.783297  0.00661986    1.77994  0.0480368  1                 304   0.6       FERRANDO_HOX11_NEIGHBORS
            5  Ferrando Hox11 Neighbors                                     0.709513  0.0429363     1.61328  0.199637   1                1484   0.666667  FERRANDO_HOX11_NEIGHBORS
           10  Ferrando Hox11 Neighbors                                     0.805292  0.0023933     1.83563  0.0688673  1                 880   0.666667  FERRANDO_HOX11_NEIGHBORS
            1  Zhan Early Differentiation Genes Dn                          0.465175  0.101253      1.43631  0.245569   1                1305   0.277778  ZHAN_EARLY_DIFFERENTIATION_GENES_DN
            5  Zhan Early Differentiation Genes Dn                          0.555485  0.0171486     1.69109  0.144939   1                2091   0.444444  ZHAN_EARLY_DIFFERENTIATION_GENES_DN
            1  Liu Vav3 Prostate Carcinogenesis Up                         -0.516093  0.0994975    -1.49525  0.232098   1                9088   0.460317  LIU_VAV3_PROSTATE_CARCINOGENESIS_UP
            1  Liu Vav3 Prostate Carcinogenesis Dn                         -0.669796  0.024484     -1.66008  0.162756   1                8220   0.636364  LIU_VAV3_PROSTATE_CARCINOGENESIS_DN
            1  Croonquist Il6 Deprivation Dn                                0.900914  0.00020141    2.02599  0.0127991  0.535952          509   0.797619  CROONQUIST_IL6_DEPRIVATION_DN
            4  Croonquist Il6 Deprivation Dn                               -0.780033  0.0166797    -1.7518   0.192508   1               10178   0.654762  CROONQUIST_IL6_DEPRIVATION_DN
            5  Croonquist Il6 Deprivation Dn                                0.82382   0.00457711    1.83586  0.0707367  1                1083   0.785714  CROONQUIST_IL6_DEPRIVATION_DN
           10  Croonquist Il6 Deprivation Dn                                0.84615   0.00140704    1.9081   0.056431   1                 921   0.785714  CROONQUIST_IL6_DEPRIVATION_DN
            1  Croonquist Stromal Stimulation Up                           -0.672098  0.0332657    -1.71924  0.146665   1                9961   0.510638  CROONQUIST_STROMAL_STIMULATION_UP
            5  Croonquist Stromal Stimulation Dn                            0.683318  0.0091017     1.74192  0.117295   1                 241   0.3       CROONQUIST_STROMAL_STIMULATION_DN
            1  Croonquist Nras Signaling Dn                                 0.902911  0.000201126   2.01243  0.0127991  0.535197          270   0.770492  CROONQUIST_NRAS_SIGNALING_DN
            4  Croonquist Nras Signaling Dn                                -0.763868  0.0267317    -1.70703  0.192508   1               10178   0.655738  CROONQUIST_NRAS_SIGNALING_DN
            5  Croonquist Nras Signaling Dn                                 0.820344  0.00656325    1.81646  0.0793115  1                1031   0.770492  CROONQUIST_NRAS_SIGNALING_DN
           10  Croonquist Nras Signaling Dn                                 0.844733  0.00160224    1.89078  0.0625516  1                 854   0.770492  CROONQUIST_NRAS_SIGNALING_DN
            1  Croonquist Nras Vs Stromal Stimulation Up                   -0.481141  0.0648334    -1.49726  0.231875   1                7483   0.545455  CROONQUIST_NRAS_VS_STROMAL_STIMULATION_UP
            1  Croonquist Nras Vs Stromal Stimulation Dn                    0.519588  0.0378109     1.64103  0.107415   1                 262   0.353659  CROONQUIST_NRAS_VS_STROMAL_STIMULATION_DN
            1  Iritani Mad1 Targets Dn                                      0.600163  0.0236064     1.71602  0.0733361  1                2794   0.659091  IRITANI_MAD1_TARGETS_DN
            1  Zhan Late Differentiation Genes Up                          -0.496104  0.0572114    -1.53195  0.214421   1                8967   0.464286  ZHAN_LATE_DIFFERENTIATION_GENES_UP
            1  Zhan Variable Early Differentiation Genes Dn                 0.490017  0.0738        1.49688  0.212322   1                2072   0.517241  ZHAN_VARIABLE_EARLY_DIFFERENTIATION_GENES_DN
            1  Ishida E2F Targets                                           0.891125  0.000602531   1.85533  0.0293627  1                 612   0.857143  ISHIDA_E2F_TARGETS
            4  Ishida E2F Targets                                          -0.83143   0.0110149    -1.72606  0.192508   1               10138   0.738095  ISHIDA_E2F_TARGETS
            5  Ishida E2F Targets                                           0.868364  0.00278884    1.79409  0.0948409  1                 935   0.857143  ISHIDA_E2F_TARGETS
           10  Ishida E2F Targets                                           0.862252  0.00199641    1.79577  0.0824479  1                 805   0.857143  ISHIDA_E2F_TARGETS
            1  Valk Aml Cluster 1                                          -0.561522  0.0730058    -1.50706  0.22774    1                8638   0.65      VALK_AML_CLUSTER_1
            1  Valk Aml Cluster 4                                          -0.64966   0.00260312   -1.85068  0.119811   1                9477   0.5       VALK_AML_CLUSTER_4
            1  Valk Aml Cluster 10                                         -0.585577  0.0143637    -1.7045   0.147212   1                9596   0.464286  VALK_AML_CLUSTER_10
            1  Valk Aml Cluster 11                                         -0.507962  0.0862872    -1.4908   0.2356     1                8395   0.571429  VALK_AML_CLUSTER_11
            1  Valk Aml With Evi1                                          -0.536828  0.0666258    -1.50494  0.228916   1                8480   0.55      VALK_AML_WITH_EVI1
            1  Valk Aml With Cebpa                                         -0.543637  0.0183099    -1.64944  0.168175   1                9477   0.458333  VALK_AML_WITH_CEBPA
            4  Poola Invasive Breast Cancer Up                             -0.599657  0.0540487    -1.66054  0.21546    1                8414   0.565217  POOLA_INVASIVE_BREAST_CANCER_UP
            1  Poola Invasive Breast Cancer Dn                             -0.513476  0.0195328    -1.6879   0.154546   1                8562   0.509615  POOLA_INVASIVE_BREAST_CANCER_DN
            1  Ding Lung Cancer By Mutation Rate                           -0.545157  0.0673387    -1.49833  0.231855   1                8515   0.529412  DING_LUNG_CANCER_BY_MUTATION_RATE
            1  Boyault Liver Cancer Subclass G3 Up                          0.628306  0.00159394    1.9172   0.0189688  1                1868   0.519553  BOYAULT_LIVER_CANCER_SUBCLASS_G3_UP
            4  Boyault Liver Cancer Subclass G3 Up                         -0.56232   0.0284407    -1.70941  0.192508   1                8603   0.530726  BOYAULT_LIVER_CANCER_SUBCLASS_G3_UP
            5  Boyault Liver Cancer Subclass G3 Up                          0.545863  0.0444801     1.66521  0.159809   1                3194   0.575419  BOYAULT_LIVER_CANCER_SUBCLASS_G3_UP
           10  Boyault Liver Cancer Subclass G3 Up                          0.549791  0.0337797     1.67421  0.155378   1                3532   0.653631  BOYAULT_LIVER_CANCER_SUBCLASS_G3_UP
            1  Boyault Liver Cancer Subclass G3 Dn                         -0.540166  0.0266398    -1.63035  0.170938   1                9734   0.40625   BOYAULT_LIVER_CANCER_SUBCLASS_G3_DN
            1  Boyault Liver Cancer Subclass G23 Up                         0.724612  0.00162635    1.89617  0.021732   1                 584   0.5       BOYAULT_LIVER_CANCER_SUBCLASS_G23_UP
            5  Boyault Liver Cancer Subclass G23 Up                         0.705897  0.00576197    1.85879  0.0620503  1                1273   0.568182  BOYAULT_LIVER_CANCER_SUBCLASS_G23_UP
           10  Boyault Liver Cancer Subclass G23 Up                         0.741484  0.000808081   1.94601  0.056431   1                 738   0.545455  BOYAULT_LIVER_CANCER_SUBCLASS_G23_UP
            1  Boyault Liver Cancer Subclass G56 Dn                        -0.562006  0.0677052    -1.50066  0.230381   1                8480   0.5       BOYAULT_LIVER_CANCER_SUBCLASS_G56_DN
            1  Boyault Liver Cancer Subclass G123 Up                        0.702407  0.00302847    1.91892  0.0189688  1                1374   0.571429  BOYAULT_LIVER_CANCER_SUBCLASS_G123_UP
            5  Boyault Liver Cancer Subclass G123 Up                        0.597689  0.0440854     1.6326   0.185717   1                1752   0.5       BOYAULT_LIVER_CANCER_SUBCLASS_G123_UP
           10  Boyault Liver Cancer Subclass G123 Up                        0.637108  0.0167406     1.7365   0.113573   1                1783   0.547619  BOYAULT_LIVER_CANCER_SUBCLASS_G123_UP
            1  Chiang Liver Cancer Subclass Ctnnb1 Up                      -0.382927  0.038948     -1.4863   0.239688   1                8020   0.442308  CHIANG_LIVER_CANCER_SUBCLASS_CTNNB1_UP
            1  Chiang Liver Cancer Subclass Ctnnb1 Dn                      -0.530649  0.0324139    -1.67945  0.157308   1                9177   0.435115  CHIANG_LIVER_CANCER_SUBCLASS_CTNNB1_DN
            1  Chiang Liver Cancer Subclass Proliferation Up                0.679804  0.00221908    1.95268  0.015348   1                 745   0.4375    CHIANG_LIVER_CANCER_SUBCLASS_PROLIFERATION_UP
            4  Chiang Liver Cancer Subclass Proliferation Up               -0.658875  0.00892857   -1.89071  0.162281   1                9603   0.576389  CHIANG_LIVER_CANCER_SUBCLASS_PROLIFERATION_UP
            5  Chiang Liver Cancer Subclass Proliferation Up                0.654922  0.00752624    1.86593  0.0620503  1                 999   0.458333  CHIANG_LIVER_CANCER_SUBCLASS_PROLIFERATION_UP
           10  Chiang Liver Cancer Subclass Proliferation Up                0.658813  0.00665323    1.8778   0.0627622  1                1262   0.486111  CHIANG_LIVER_CANCER_SUBCLASS_PROLIFERATION_UP
            1  Chiang Liver Cancer Subclass Unannotated Up                 -0.432415  0.0210463    -1.57412  0.195122   1                7560   0.551724  CHIANG_LIVER_CANCER_SUBCLASS_UNANNOTATED_UP
            1  Chiang Liver Cancer Subclass Unannotated Dn                  0.563403  0.0113818     1.8035   0.0411858  1                2548   0.527778  CHIANG_LIVER_CANCER_SUBCLASS_UNANNOTATED_DN
            4  Chiang Liver Cancer Subclass Unannotated Dn                 -0.544974  0.0188867    -1.74333  0.192508   1                8860   0.533333  CHIANG_LIVER_CANCER_SUBCLASS_UNANNOTATED_DN
           10  Chiang Liver Cancer Subclass Unannotated Dn                  0.51322   0.0414178     1.64273  0.174876   1                2443   0.45      CHIANG_LIVER_CANCER_SUBCLASS_UNANNOTATED_DN
            1  Coulouarn Temporal Tgfb1 Signature Dn                       -0.380516  0.0104018    -1.58671  0.189217   1                9104   0.330189  COULOUARN_TEMPORAL_TGFB1_SIGNATURE_DN
            4  Kaposi Liver Cancer Met Up                                  -0.626192  0.00997208   -1.7357   0.192508   1                8103   0.722222  KAPOSI_LIVER_CANCER_MET_UP
            1  Lee Liver Cancer Hepatoblast                                -0.684377  0.0267122    -1.66326  0.162484   1                8888   0.6       LEE_LIVER_CANCER_HEPATOBLAST
            1  Wang Recurrent Liver Cancer Up                               0.505339  0.0824721     1.45112  0.238115   1                3161   0.631579  WANG_RECURRENT_LIVER_CANCER_UP
           10  Wang Recurrent Liver Cancer Up                               0.608244  0.0101379     1.73901  0.113573   1                2204   0.526316  WANG_RECURRENT_LIVER_CANCER_UP
            1  Woo Liver Cancer Recurrence Up                              -0.457969  0.0586207    -1.56153  0.199949   1                9832   0.34375   WOO_LIVER_CANCER_RECURRENCE_UP
            1  Montero Thyroid Cancer Poor Survival Up                      0.952205  0.00461662    1.65911  0.0991974  1                 165   0.9       MONTERO_THYROID_CANCER_POOR_SURVIVAL_UP
           10  Montero Thyroid Cancer Poor Survival Up                      0.91158   0.0232932     1.58141  0.229077   1                 327   0.9       MONTERO_THYROID_CANCER_POOR_SURVIVAL_UP
            1  Meissner Npc Hcp With H3K4Me2 And H3K27Me3                  -0.440376  0.024577     -1.59769  0.181764   1                8411   0.414894  MEISSNER_NPC_HCP_WITH_H3K4ME2_AND_H3K27ME3
            1  Yoshioka Liver Cancer Early Recurrence Up                    0.460418  0.0349848     1.52386  0.18751    1                1956   0.37931   YOSHIOKA_LIVER_CANCER_EARLY_RECURRENCE_UP
            1  Mikkelsen Mcv6 Hcp With H3K27Me3                            -0.486494  0.00881764   -1.73053  0.143251   1                8950   0.417391  MIKKELSEN_MCV6_HCP_WITH_H3K27ME3
            1  Mikkelsen Mcv6 Icp With H3K4Me3 And H3K27Me3                -0.679063  0.0272507    -1.63268  0.170938   1               10336   0.416667  MIKKELSEN_MCV6_ICP_WITH_H3K4ME3_AND_H3K27ME3
            1  Mikkelsen Mef Lcp With H3K4Me3                              -0.555824  0.0130078    -1.78449  0.131106   1                9865   0.382022  MIKKELSEN_MEF_LCP_WITH_H3K4ME3
            1  Mikkelsen Ips With Hcp H3K27Me3                             -0.669271  0.0245048    -1.64034  0.170938   1                8128   0.75      MIKKELSEN_IPS_WITH_HCP_H3K27ME3
            1  Mikkelsen Ips Icp With H3K4Me3 And H327Me3                  -0.51542   0.0392835    -1.59916  0.181764   1                8246   0.456522  MIKKELSEN_IPS_ICP_WITH_H3K4ME3_AND_H327ME3
            1  Mikkelsen Ips Lcp With H3K4Me3                              -0.453732  0.0511945    -1.5534   0.202734   1                8931   0.3625    MIKKELSEN_IPS_LCP_WITH_H3K4ME3
            1  Ono Foxp3 Targets Up                                        -0.706046  0.00987505   -1.75593  0.139914   1                9854   0.357143  ONO_FOXP3_TARGETS_UP
            1  Schoen Nfkb Signaling                                       -0.665059  0.017607     -1.7359   0.141912   1                9051   0.666667  SCHOEN_NFKB_SIGNALING
            1  Ngo Malignant Glioma 1P Loh                                  0.5861    0.106237      1.4462   0.238115   1                1242   0.4375    NGO_MALIGNANT_GLIOMA_1P_LOH
            4  Ngo Malignant Glioma 1P Loh                                 -0.712523  0.0108911    -1.7662   0.192508   1                9557   0.6875    NGO_MALIGNANT_GLIOMA_1P_LOH
            1  Winnepenninckx Melanoma Metastasis Dn                       -0.605268  0.0174916    -1.71866  0.146665   1                8645   0.576923  WINNEPENNINCKX_MELANOMA_METASTASIS_DN
            1  Pujana Breast Cancer With Brca1 Mutated Up                   0.801187  0.00183711    1.82187  0.0363179  1                1674   0.854167  PUJANA_BREAST_CANCER_WITH_BRCA1_MUTATED_UP
            5  Pujana Breast Cancer With Brca1 Mutated Up                   0.829147  0.000392696   1.88829  0.0559509  1                1572   0.916667  PUJANA_BREAST_CANCER_WITH_BRCA1_MUTATED_UP
           10  Pujana Breast Cancer With Brca1 Mutated Up                   0.805179  0.00163767    1.83082  0.0688673  1                1251   0.791667  PUJANA_BREAST_CANCER_WITH_BRCA1_MUTATED_UP
            1  Syed Estradiol Response                                      0.494662  0.0527054     1.49222  0.213453   1                1853   0.411765  SYED_ESTRADIOL_RESPONSE
            1  Stein Esr1 Targets                                           0.482871  0.0124821     1.73419  0.0649777  1                 888   0.309859  STEIN_ESR1_TARGETS
            5  Stein Esr1 Targets                                           0.485503  0.0121704     1.74183  0.117295   1                 910   0.295775  STEIN_ESR1_TARGETS
           10  Stein Esr1 Targets                                           0.450105  0.0274299     1.61985  0.192086   1                 941   0.295775  STEIN_ESR1_TARGETS
            1  Kobayashi Egfr Signaling 24Hr Up                            -0.46459   0.0107921    -1.69131  0.152216   1                7732   0.45977   KOBAYASHI_EGFR_SIGNALING_24HR_UP
            1  Kobayashi Egfr Signaling 24Hr Dn                             0.810706  0.000202716   2.09203  0.0127991  0.539428          924   0.699482  KOBAYASHI_EGFR_SIGNALING_24HR_DN
            4  Kobayashi Egfr Signaling 24Hr Dn                            -0.690025  0.0214019    -1.7841   0.192508   1                9498   0.61658   KOBAYASHI_EGFR_SIGNALING_24HR_DN
            5  Kobayashi Egfr Signaling 24Hr Dn                             0.772043  0.000798244   1.98167  0.0380769  1                1270   0.683938  KOBAYASHI_EGFR_SIGNALING_24HR_DN
           10  Kobayashi Egfr Signaling 24Hr Dn                             0.744107  0.00359425    1.93431  0.056431   1                1076   0.637306  KOBAYASHI_EGFR_SIGNALING_24HR_DN
            1  Fournier Acinar Development Late 2                           0.632754  0.000402334   2.02405  0.0127991  1                1663   0.522088  FOURNIER_ACINAR_DEVELOPMENT_LATE_2
            4  Fournier Acinar Development Late 2                          -0.530033  0.0379297    -1.6853   0.201047   1                8896   0.518072  FOURNIER_ACINAR_DEVELOPMENT_LATE_2
            5  Fournier Acinar Development Late 2                           0.579495  0.00795545    1.84316  0.0674856  1                1266   0.389558  FOURNIER_ACINAR_DEVELOPMENT_LATE_2
           10  Fournier Acinar Development Late 2                           0.589076  0.00462033    1.88212  0.0627622  1                2482   0.566265  FOURNIER_ACINAR_DEVELOPMENT_LATE_2
            1  Boylan Multiple Myeloma Pca3 Up                             -0.400654  0.0132877    -1.53876  0.210383   1                9060   0.37037   BOYLAN_MULTIPLE_MYELOMA_PCA3_UP
            5  Hoshida Liver Cancer Subclass S2                             0.420589  0.0336529     1.56477  0.24508    1                1988   0.373832  HOSHIDA_LIVER_CANCER_SUBCLASS_S2
            1  Pyeon Cancer Head And Neck Vs Cervical Up                    0.559192  0.0039151     1.89909  0.0213982  1                1358   0.379747  PYEON_CANCER_HEAD_AND_NECK_VS_CERVICAL_UP
            5  Pyeon Cancer Head And Neck Vs Cervical Up                    0.619115  0.000397298   2.09505  0.0336414  1                1565   0.474684  PYEON_CANCER_HEAD_AND_NECK_VS_CERVICAL_UP
           10  Pyeon Cancer Head And Neck Vs Cervical Up                    0.5689    0.00284669    1.93235  0.056431   1                2171   0.481013  PYEON_CANCER_HEAD_AND_NECK_VS_CERVICAL_UP
            1  Cairo Hepatoblastoma Classes Dn                             -0.585333  0.00119713   -1.99335  0.116729   1                8873   0.509677  CAIRO_HEPATOBLASTOMA_CLASSES_DN
            1  Gutierrez Chronic Lymphocytic Leukemia Dn                   -0.481621  0.0196928    -1.63501  0.170938   1                9104   0.413043  GUTIERREZ_CHRONIC_LYMPHOCYTIC_LEUKEMIA_DN
            1  Winnepenninckx Melanoma Metastasis Up                        0.764978  0.000201694   2.04922  0.0127991  0.536708          672   0.555556  WINNEPENNINCKX_MELANOMA_METASTASIS_UP
            4  Winnepenninckx Melanoma Metastasis Up                       -0.661414  0.021991     -1.77164  0.192508   1                9032   0.651852  WINNEPENNINCKX_MELANOMA_METASTASIS_UP
            5  Winnepenninckx Melanoma Metastasis Up                        0.689118  0.00835322    1.84231  0.0674856  1                1760   0.614815  WINNEPENNINCKX_MELANOMA_METASTASIS_UP
           10  Winnepenninckx Melanoma Metastasis Up                        0.711982  0.00279944    1.90928  0.056431   1                1558   0.622222  WINNEPENNINCKX_MELANOMA_METASTASIS_UP
            1  Uzonyi Response To Leukotriene And Thrombin                 -0.705281  0.0447195    -1.65191  0.168175   1                9961   0.59375   UZONYI_RESPONSE_TO_LEUKOTRIENE_AND_THROMBIN
            1  Vantveer Breast Cancer Poor Prognosis                        0.446087  0.0677932     1.48691  0.215694   1                1085   0.341463  VANTVEER_BREAST_CANCER_POOR_PROGNOSIS
           10  Vantveer Breast Cancer Poor Prognosis                        0.511316  0.0138973     1.69846  0.131604   1                 348   0.317073  VANTVEER_BREAST_CANCER_POOR_PROGNOSIS
            5  Kyng Werner Syndrom And Normal Aging Up                      0.39721   0.0120114     1.59941  0.209579   1                1835   0.30137   KYNG_WERNER_SYNDROM_AND_NORMAL_AGING_UP
            1  Chiaretti T All Relapse Prognosis                            0.750996  0.00222807    1.84725  0.0308973  1                 660   0.470588  CHIARETTI_T_ALL_RELAPSE_PROGNOSIS
            5  Chiaretti T All Relapse Prognosis                            0.73234   0.00414283    1.78893  0.0965026  1                 572   0.470588  CHIARETTI_T_ALL_RELAPSE_PROGNOSIS
           10  Chiaretti T All Relapse Prognosis                            0.793889  0.000201045   1.96076  0.056431   0.534982          475   0.470588  CHIARETTI_T_ALL_RELAPSE_PROGNOSIS
            1  Zhan Multiple Myeloma Cd1 Dn                                -0.431348  0.047619     -1.48025  0.242366   1                8562   0.416667  ZHAN_MULTIPLE_MYELOMA_CD1_DN
            1  Zhan Multiple Myeloma Pr Dn                                 -0.532111  0.0584634    -1.55505  0.202734   1                7146   0.72973   ZHAN_MULTIPLE_MYELOMA_PR_DN
            1  Dang Regulated By Myc Up                                     0.58917   0.0030036     1.83675  0.0321985  1                2271   0.6       DANG_REGULATED_BY_MYC_UP
           10  Dang Regulated By Myc Up                                     0.499316  0.0520564     1.5711   0.238345   1                2641   0.507692  DANG_REGULATED_BY_MYC_UP
            1  Dang Regulated By Myc Dn                                    -0.403122  0.0449619    -1.56415  0.199134   1                9060   0.327273  DANG_REGULATED_BY_MYC_DN
            1  Dang Myc Targets Up                                          0.605819  0.00778443    1.81357  0.0381776  1                2181   0.572581  DANG_MYC_TARGETS_UP
            4  Dang Myc Targets Up                                         -0.559075  0.0329519    -1.68647  0.201047   1                8552   0.564516  DANG_MYC_TARGETS_UP
            1  Tian Tnf Signaling Not Via Nfkb                             -0.688742  0.0362611    -1.69432  0.150339   1                9061   0.684211  TIAN_TNF_SIGNALING_NOT_VIA_NFKB
            1  Wong Embryonic Stem Cell Core                                0.702658  0.000602894   2.01837  0.0127991  1                1543   0.583893  WONG_EMBRYONIC_STEM_CELL_CORE
            4  Wong Embryonic Stem Cell Core                               -0.636175  0.00748621   -1.82796  0.187832   1                8903   0.59396   WONG_EMBRYONIC_STEM_CELL_CORE
            5  Wong Embryonic Stem Cell Core                                0.608089  0.0233109     1.76116  0.106138   1                1983   0.486577  WONG_EMBRYONIC_STEM_CELL_CORE
           10  Wong Embryonic Stem Cell Core                                0.619488  0.013581      1.78718  0.0877569  1                2271   0.580537  WONG_EMBRYONIC_STEM_CELL_CORE
            1  Nakayama Soft Tissue Tumors Pca2 Up                          0.752101  0.00406587    1.90731  0.0202154  1                 206   0.525424  NAKAYAMA_SOFT_TISSUE_TUMORS_PCA2_UP
            4  Nakayama Soft Tissue Tumors Pca2 Up                         -0.704811  0.0237525    -1.77864  0.192508   1               10272   0.576271  NAKAYAMA_SOFT_TISSUE_TUMORS_PCA2_UP
            5  Nakayama Soft Tissue Tumors Pca2 Up                          0.701492  0.0243471     1.78717  0.0965026  1                 425   0.491525  NAKAYAMA_SOFT_TISSUE_TUMORS_PCA2_UP
           10  Nakayama Soft Tissue Tumors Pca2 Up                          0.680157  0.0372362     1.72949  0.115709   1                 622   0.542373  NAKAYAMA_SOFT_TISSUE_TUMORS_PCA2_UP
            1  Nakayama Soft Tissue Tumors Pca2 Dn                         -0.690528  0.0145308    -1.82781  0.119811   1                9205   0.672727  NAKAYAMA_SOFT_TISSUE_TUMORS_PCA2_DN
            5  Chandran Metastasis Up                                       0.416414  0.0124298     1.67735  0.153251   1                2756   0.410112  CHANDRAN_METASTASIS_UP
            1  Chandran Metastasis Dn                                      -0.512602  0.000606428  -1.95048  0.116729   1                8758   0.475207  CHANDRAN_METASTASIS_DN
            1  Mikkelsen Es Icp With H3K27Me3                              -0.688137  0.0300421    -1.59858  0.181764   1                9572   0.5       MIKKELSEN_ES_ICP_WITH_H3K27ME3
            1  Mikkelsen Es Icp With H3K4Me3 And H3K27Me3                  -0.451894  0.0659475    -1.48932  0.237102   1                8722   0.381818  MIKKELSEN_ES_ICP_WITH_H3K4ME3_AND_H3K27ME3
            1  Mikkelsen Es Lcp With H3K4Me3                               -0.456968  0.0703518    -1.50093  0.230381   1                8369   0.423729  MIKKELSEN_ES_LCP_WITH_H3K4ME3
            1  Mikkelsen Npc Hcp With H3K4Me3 And H3K27Me3                 -0.478237  0.0077126    -1.72753  0.144779   1                8937   0.379747  MIKKELSEN_NPC_HCP_WITH_H3K4ME3_AND_H3K27ME3
            1  Yao Temporal Response To Progesterone Cluster 0             -0.608502  0.000806614  -2.00104  0.116729   1                8964   0.52381   YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_0
            1  Yao Temporal Response To Progesterone Cluster 1             -0.557391  0.000599161  -1.91201  0.116729   1                9090   0.491525  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_1
            1  Yao Temporal Response To Progesterone Cluster 3             -0.580605  0.0646204    -1.47961  0.242366   1                7779   0.727273  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_3
            1  Yao Temporal Response To Progesterone Cluster 10             0.462631  0.0975174     1.46871  0.22633    1                2266   0.416667  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_10
            4  Yao Temporal Response To Progesterone Cluster 10            -0.514789  0.0327056    -1.62582  0.230529   1                8918   0.533333  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_10
            1  Yao Temporal Response To Progesterone Cluster 11             0.523044  0.0286286     1.67065  0.0946504  1                1784   0.462366  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_11
           10  Yao Temporal Response To Progesterone Cluster 11             0.488028  0.0654857     1.56135  0.248199   1                2193   0.419355  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_11
            1  Yao Temporal Response To Progesterone Cluster 13             0.479252  0.139738      1.44697  0.238115   1                3565   0.607595  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_13
            4  Yao Temporal Response To Progesterone Cluster 13            -0.589702  0.0185148    -1.77298  0.192508   1                8719   0.594937  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_13
            1  Yao Temporal Response To Progesterone Cluster 14             0.546848  0.0116209     1.78125  0.0477546  1                2687   0.603053  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_14
            4  Yao Temporal Response To Progesterone Cluster 14            -0.641678  0.000196734  -2.09129  0.112621   0.52351          8593   0.648855  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_14
            1  Yao Temporal Response To Progesterone Cluster 17             0.439341  0.110175      1.44962  0.238115   1                3088   0.553571  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_17
            4  Yao Temporal Response To Progesterone Cluster 17            -0.563104  0.00436508   -1.8662   0.162281   1                8551   0.577381  YAO_TEMPORAL_RESPONSE_TO_PROGESTERONE_CLUSTER_17
            8  Mootha Ffa Oxydation                                         0.725195  0.000395023   2.07562  0.104006   1                1818   0.55      MOOTHA_FFA_OXYDATION
           10  Mootha Gluconeogenesis                                       0.61538   0.0274803     1.65093  0.170205   1                1411   0.545455  MOOTHA_GLUCONEOGENESIS
           10  Mootha Glycogen Metabolism                                   0.591869  0.0108805     1.73231  0.115709   1                3600   0.75      MOOTHA_GLYCOGEN_METABOLISM
           10  Mootha Glycolysis                                            0.653769  0.034421      1.64113  0.174876   1                1411   0.533333  MOOTHA_GLYCOLYSIS
            1  Nakamura Adipogenesis Early Up                              -0.565058  0.00730816   -1.80494  0.124176   1                9633   0.4375    NAKAMURA_ADIPOGENESIS_EARLY_UP
            1  Nakamura Adipogenesis Early Dn                              -0.695682  0.0119312    -1.77706  0.134462   1                9946   0.529412  NAKAMURA_ADIPOGENESIS_EARLY_DN
            1  Nakamura Adipogenesis Late Up                               -0.43965   0.00707929   -1.73762  0.141912   1                8573   0.467391  NAKAMURA_ADIPOGENESIS_LATE_UP
            1  Nakamura Adipogenesis Late Dn                               -0.62866   0.0413091    -1.66573  0.162484   1                7898   0.764706  NAKAMURA_ADIPOGENESIS_LATE_DN
            1  Gerhold Response To Tzd Dn                                  -0.680478  0.0208544    -1.69165  0.152216   1                9163   0.583333  GERHOLD_RESPONSE_TO_TZD_DN
            1  Zembutsu Sensitivity To Cisplatin                           -0.604962  0.012966     -1.67988  0.157308   1                9865   0.384615  ZEMBUTSU_SENSITIVITY_TO_CISPLATIN
            1  Zembutsu Sensitivity To Nimustine                            0.523658  0.0495768     1.47782  0.219334   1                 690   0.25      ZEMBUTSU_SENSITIVITY_TO_NIMUSTINE
            5  Zembutsu Sensitivity To Methotrexate                         0.68588   0.00280056    1.85829  0.0620503  1                 339   0.416667  ZEMBUTSU_SENSITIVITY_TO_METHOTREXATE
            1  Dazard Uv Response Cluster G5                                0.686274  0.0902736     1.47164  0.225424   1                1984   0.5       DAZARD_UV_RESPONSE_CLUSTER_G5
            4  Dazard Uv Response Cluster G5                               -0.793745  0.00841852   -1.70011  0.193938   1                9422   0.6       DAZARD_UV_RESPONSE_CLUSTER_G5
            1  Dazard Uv Response Cluster G24                              -0.539473  0.0750401    -1.48753  0.238525   1                9888   0.352941  DAZARD_UV_RESPONSE_CLUSTER_G24
            1  Dazard Uv Response Cluster G28                              -0.711657  0.0431887    -1.60596  0.180997   1                8501   0.75      DAZARD_UV_RESPONSE_CLUSTER_G28
            1  Nielsen Gist                                                -0.410269  0.0406835    -1.52529  0.217927   1                7851   0.457831  NIELSEN_GIST
            1  Nielsen Leiomyosarcoma Cnn1 Dn                              -0.665628  0.0697581    -1.532    0.214421   1                9862   0.470588  NIELSEN_LEIOMYOSARCOMA_CNN1_DN
            1  Nielsen Schwannoma Dn                                       -0.624368  0.0900739    -1.51052  0.226462   1                8298   0.705882  NIELSEN_SCHWANNOMA_DN
            1  Dasu Il6 Signaling Up                                       -0.537603  0.0318329    -1.68192  0.157308   1               10112   0.346939  DASU_IL6_SIGNALING_UP
            1  Whitfield Cell Cycle Literature                              0.895144  0.000201369   1.8916   0.0220871  0.535844          510   0.861111  WHITFIELD_CELL_CYCLE_LITERATURE
            4  Whitfield Cell Cycle Literature                             -0.768153  0.0448201    -1.62228  0.234192   1               10328   0.694444  WHITFIELD_CELL_CYCLE_LITERATURE
            5  Whitfield Cell Cycle Literature                              0.841483  0.00555666    1.76985  0.102316   1                 997   0.861111  WHITFIELD_CELL_CYCLE_LITERATURE
           10  Whitfield Cell Cycle Literature                              0.804091  0.0201516     1.70211  0.130517   1                 921   0.805556  WHITFIELD_CELL_CYCLE_LITERATURE
            1  Jiang Hypoxia Via Vhl                                        0.489086  0.0359712     1.58068  0.151135   1                2776   0.46875   JIANG_HYPOXIA_VIA_VHL
            1  Wu Apoptosis By Cdkn1A Via Tp53                              0.831075  0.000607287   1.92114  0.0188436  1                 499   0.729167  WU_APOPTOSIS_BY_CDKN1A_VIA_TP53
            5  Wu Apoptosis By Cdkn1A Via Tp53                              0.842083  0.000790826   1.94072  0.0454213  1                 882   0.833333  WU_APOPTOSIS_BY_CDKN1A_VIA_TP53
           10  Wu Apoptosis By Cdkn1A Via Tp53                              0.809611  0.00280955    1.8713   0.0627622  1                 910   0.770833  WU_APOPTOSIS_BY_CDKN1A_VIA_TP53
            1  Wu Apoptosis By Cdkn1A Not Via Tp53                          0.585016  0.0963139     1.44815  0.238115   1                 168   0.363636  WU_APOPTOSIS_BY_CDKN1A_NOT_VIA_TP53
            1  Browne Hcmv Infection 30Min Up                              -0.624858  0.00224581   -1.9486   0.116729   1               10291   0.40625   BROWNE_HCMV_INFECTION_30MIN_UP
            1  Browne Hcmv Infection 2Hr Up                                -0.564737  0.0423015    -1.62776  0.170938   1                9961   0.375     BROWNE_HCMV_INFECTION_2HR_UP
            1  Whitfield Cell Cycle G1 S                                    0.572708  0.0006035     2.02239  0.0127991  1                 975   0.366071  WHITFIELD_CELL_CYCLE_G1_S
            5  Whitfield Cell Cycle G1 S                                    0.582378  0.000394244   2.0385   0.0336414  1                1994   0.473214  WHITFIELD_CELL_CYCLE_G1_S
           10  Whitfield Cell Cycle G1 S                                    0.579435  0.000405844   2.02554  0.056431   1                1423   0.392857  WHITFIELD_CELL_CYCLE_G1_S
            1  Whitfield Cell Cycle S                                       0.560115  0.00101482    1.97435  0.013915   1                1361   0.361345  WHITFIELD_CELL_CYCLE_S
            5  Whitfield Cell Cycle S                                       0.566608  0.000398645   1.97339  0.0383868  1                1882   0.420168  WHITFIELD_CELL_CYCLE_S
           10  Whitfield Cell Cycle S                                       0.522829  0.00590752    1.82921  0.0688673  1                 622   0.285714  WHITFIELD_CELL_CYCLE_S
            1  Whitfield Cell Cycle G2                                      0.504707  0.00654798    1.84769  0.0308973  1                1351   0.342657  WHITFIELD_CELL_CYCLE_G2
            5  Whitfield Cell Cycle G2                                      0.559123  0.000197981   2.04506  0.0336414  0.526826         1041   0.356643  WHITFIELD_CELL_CYCLE_G2
           10  Whitfield Cell Cycle G2                                      0.519007  0.00380838    1.90678  0.056431   1                1513   0.363636  WHITFIELD_CELL_CYCLE_G2
            1  Whitfield Cell Cycle G2 M                                    0.542788  0.00389584    1.86108  0.0278677  1                1674   0.407609  WHITFIELD_CELL_CYCLE_G2_M
            5  Whitfield Cell Cycle G2 M                                    0.563266  0.00178112    1.93305  0.0454213  1                1753   0.407609  WHITFIELD_CELL_CYCLE_G2_M
           10  Whitfield Cell Cycle G2 M                                    0.526083  0.0104796     1.81046  0.0756726  1                1223   0.347826  WHITFIELD_CELL_CYCLE_G2_M
            1  Whitfield Cell Cycle M G1                                    0.374268  0.0481155     1.48728  0.215694   1                2183   0.356061  WHITFIELD_CELL_CYCLE_M_G1
            5  Whitfield Cell Cycle M G1                                    0.413158  0.013945      1.64126  0.182655   1                2783   0.409091  WHITFIELD_CELL_CYCLE_M_G1
            1  Karlsson Tgfb1 Targets Up                                    0.50494   0.0236995     1.69619  0.0811864  1                1784   0.464912  KARLSSON_TGFB1_TARGETS_UP
            4  Karlsson Tgfb1 Targets Up                                   -0.489393  0.0343738    -1.63571  0.229428   1                9267   0.438596  KARLSSON_TGFB1_TARGETS_UP
            1  Karlsson Tgfb1 Targets Dn                                   -0.394742  0.00419245   -1.72033  0.146665   1                8934   0.370968  KARLSSON_TGFB1_TARGETS_DN
            1  Wang Metastasis Of Breast Cancer Esr1 Up                     0.707027  0.0262683     1.66741  0.0958915  1                 164   0.444444  WANG_METASTASIS_OF_BREAST_CANCER_ESR1_UP
            4  Wang Metastasis Of Breast Cancer Esr1 Up                    -0.718517  0.0180484    -1.68161  0.203728   1                9364   0.611111  WANG_METASTASIS_OF_BREAST_CANCER_ESR1_UP
            5  Wang Metastasis Of Breast Cancer Esr1 Up                     0.711683  0.0257732     1.67045  0.156755   1                 904   0.555556  WANG_METASTASIS_OF_BREAST_CANCER_ESR1_UP
            1  Wang Metastasis Of Breast Cancer Esr1 Dn                    -0.503092  0.0678543    -1.49724  0.231875   1                8684   0.416667  WANG_METASTASIS_OF_BREAST_CANCER_ESR1_DN
            1  Martens Tretinoin Response Up                               -0.401991  0.0377126    -1.5477   0.203357   1                8701   0.358974  MARTENS_TRETINOIN_RESPONSE_UP
            1  Verhaak Glioblastoma Classical                              -0.418715  0.053146     -1.52415  0.218871   1                7960   0.45082   VERHAAK_GLIOBLASTOMA_CLASSICAL
            1  Chicas Rb1 Targets Low Serum                                 0.53562   0.0393495     1.64458  0.106001   1                2167   0.432432  CHICAS_RB1_TARGETS_LOW_SERUM
            4  Chicas Rb1 Targets Low Serum                                -0.544566  0.0311067    -1.67171  0.209664   1                9826   0.459459  CHICAS_RB1_TARGETS_LOW_SERUM
            1  Chicas Rb1 Targets Growing                                   0.469624  0.0286516     1.66114  0.0984775  1                1194   0.377551  CHICAS_RB1_TARGETS_GROWING
            5  Chicas Rb1 Targets Growing                                   0.534027  0.00434439    1.87576  0.0582517  1                 940   0.367347  CHICAS_RB1_TARGETS_GROWING
           10  Chicas Rb1 Targets Growing                                   0.492301  0.0155782     1.72124  0.118911   1                 886   0.352041  CHICAS_RB1_TARGETS_GROWING
            1  Wang Response To Gsk3 Inhibitor Sb216763 Dn                  0.600111  0.000806289   1.98337  0.0134874  1                1180   0.422383  WANG_RESPONSE_TO_GSK3_INHIBITOR_SB216763_DN
            4  Wang Response To Gsk3 Inhibitor Sb216763 Dn                 -0.502962  0.0412636    -1.65702  0.216995   1                8981   0.501805  WANG_RESPONSE_TO_GSK3_INHIBITOR_SB216763_DN
            5  Wang Response To Gsk3 Inhibitor Sb216763 Dn                  0.500852  0.0412022     1.64664  0.178715   1                2050   0.422383  WANG_RESPONSE_TO_GSK3_INHIBITOR_SB216763_DN
           10  Wang Response To Gsk3 Inhibitor Sb216763 Dn                  0.562997  0.00478278    1.86947  0.0627622  1                1286   0.393502  WANG_RESPONSE_TO_GSK3_INHIBITOR_SB216763_DN
            1  Hoelzel Nf1 Targets Up                                      -0.449139  0.0709976    -1.50657  0.22774    1                7714   0.551724  HOELZEL_NF1_TARGETS_UP
            1  Hoelzel Nf1 Targets Dn                                      -0.539489  0.0544327    -1.63949  0.170938   1                8204   0.626667  HOELZEL_NF1_TARGETS_DN
            1  Dutertre Estradiol Response 6Hr Up                           0.418884  0.0376226     1.57871  0.15239    1                1250   0.304569  DUTERTRE_ESTRADIOL_RESPONSE_6HR_UP
            5  Dutertre Estradiol Response 6Hr Up                           0.449757  0.0217479     1.69396  0.143083   1                1541   0.314721  DUTERTRE_ESTRADIOL_RESPONSE_6HR_UP
            1  Dutertre Estradiol Response 6Hr Dn                          -0.533332  0.00528133   -1.82763  0.119811   1                7817   0.576471  DUTERTRE_ESTRADIOL_RESPONSE_6HR_DN
            1  Dutertre Estradiol Response 24Hr Up                          0.752776  0.000400561   2.116    0.0127991  1                 745   0.6       DUTERTRE_ESTRADIOL_RESPONSE_24HR_UP
            4  Dutertre Estradiol Response 24Hr Up                         -0.581876  0.0815722    -1.6188   0.237236   1                9250   0.526923  DUTERTRE_ESTRADIOL_RESPONSE_24HR_UP
            5  Dutertre Estradiol Response 24Hr Up                          0.741784  0.00139637    2.05989  0.0336414  1                 935   0.6       DUTERTRE_ESTRADIOL_RESPONSE_24HR_UP
           10  Dutertre Estradiol Response 24Hr Up                          0.726462  0.00140028    2.03009  0.056431   1                1098   0.615385  DUTERTRE_ESTRADIOL_RESPONSE_24HR_UP
            1  Figueroa Aml Methylation Cluster 3 Dn                       -0.638106  0.00020016   -1.97641  0.116729   0.532626         9766   0.384615  FIGUEROA_AML_METHYLATION_CLUSTER_3_DN
            1  Wang Adipogenic Genes Repressed By Sirt1                     0.603796  0.134838      1.43853  0.244682   1                1577   0.578947  WANG_ADIPOGENIC_GENES_REPRESSED_BY_SIRT1
            1  Pangas Tumor Suppression By Smad1 And Smad5 Up              -0.530731  0.0237951    -1.71928  0.146665   1                8493   0.514286  PANGAS_TUMOR_SUPPRESSION_BY_SMAD1_AND_SMAD5_UP
            1  Li Dcp2 Bound Mrna                                           0.577779  0.0228311     1.70023  0.0795648  1                2985   0.653846  LI_DCP2_BOUND_MRNA
            4  Li Dcp2 Bound Mrna                                          -0.575916  0.025956     -1.68964  0.201047   1                8762   0.653846  LI_DCP2_BOUND_MRNA
            1  Ohguchi Liver Hnf4A Targets Up                              -0.593586  0.0075652    -1.77215  0.137045   1                8908   0.521739  OHGUCHI_LIVER_HNF4A_TARGETS_UP
            1  Thillainadesan Znf217 Targets Up                             0.501387  0.0201075     1.63838  0.108703   1                2334   0.5       THILLAINADESAN_ZNF217_TARGETS_UP
            5  Thillainadesan Znf217 Targets Up                             0.541613  0.00539245    1.76576  0.102916   1                1253   0.388889  THILLAINADESAN_ZNF217_TARGETS_UP
            1  Chyla Cbfa2T3 Targets Up                                    -0.451026  0.00645682   -1.75849  0.139914   1                8326   0.433333  CHYLA_CBFA2T3_TARGETS_UP
            1  Wierenga Stat5A Targets Group2                              -0.566422  0.0371039    -1.65773  0.16369    1                9135   0.428571  WIERENGA_STAT5A_TARGETS_GROUP2
            1  Wierenga Pml Interactome                                     0.612868  0.00140534    1.88796  0.0226571  1                1941   0.5       WIERENGA_PML_INTERACTOME
            1  Acosta Proliferation Independent Myc Targets Up              0.501137  0.0107463     1.73356  0.0650154  1                2546   0.513889  ACOSTA_PROLIFERATION_INDEPENDENT_MYC_TARGETS_UP
           10  Acosta Proliferation Independent Myc Targets Up              0.475031  0.0251447     1.64825  0.172777   1                2581   0.458333  ACOSTA_PROLIFERATION_INDEPENDENT_MYC_TARGETS_UP
            1  Kang Ar Targets Up                                          -0.820475  0.00180977   -1.8964   0.119811   1               10150   0.5       KANG_AR_TARGETS_UP
            1  Bhat Esr1 Targets Not Via Akt1 Dn                           -0.566652  0.00506997   -1.82505  0.119811   1                8404   0.589744  BHAT_ESR1_TARGETS_NOT_VIA_AKT1_DN
            1  Bhat Esr1 Targets Via Akt1 Dn                               -0.523426  0.00765203   -1.74394  0.139914   1                7534   0.662162  BHAT_ESR1_TARGETS_VIA_AKT1_DN
            1  Miyagawa Targets Of Ewsr1 Ets Fusions Up                    -0.387083  0.0503525    -1.50368  0.229309   1                8083   0.454545  MIYAGAWA_TARGETS_OF_EWSR1_ETS_FUSIONS_UP
            1  Miyagawa Targets Of Ewsr1 Ets Fusions Dn                    -0.533758  0.0302538    -1.69889  0.148746   1                8756   0.52071   MIYAGAWA_TARGETS_OF_EWSR1_ETS_FUSIONS_DN
            1  Kim Glis2 Targets Up                                        -0.650782  0.0430876    -1.66471  0.162484   1                9519   0.52      KIM_GLIS2_TARGETS_UP
            1  Steger Adipogenesis Dn                                      -0.770094  0.0320565    -1.6871   0.154762   1                8780   0.904762  STEGER_ADIPOGENESIS_DN
            1  Zhu Skil Targets Up                                         -0.687701  0.0277608    -1.66053  0.162756   1                8950   0.736842  ZHU_SKIL_TARGETS_UP
            1  Kim Pten Targets Up                                         -0.53362   0.0617932    -1.49376  0.233584   1                9946   0.333333  KIM_PTEN_TARGETS_UP
            1  Pasini Suz12 Targets Dn                                     -0.457175  0.0720048    -1.54999  0.202734   1                8108   0.510563  PASINI_SUZ12_TARGETS_DN
            1  Bakker Foxo3 Targets Up                                     -0.54274   0.00141729   -1.86698  0.119811   1                7871   0.588235  BAKKER_FOXO3_TARGETS_UP
            1  Vandesluis Normal Embryos Dn                                -0.566294  0.0334753    -1.57373  0.195122   1                9962   0.416667  VANDESLUIS_NORMAL_EMBRYOS_DN
            1  Azare Neoplastic Transformation By Stat3 Up                 -0.458957  0.0977323    -1.47282  0.246295   1                8573   0.46875   AZARE_NEOPLASTIC_TRANSFORMATION_BY_STAT3_UP
            1  Delpuech Foxo3 Targets Up                                   -0.426382  0.0294835    -1.55878  0.201403   1                8361   0.419355  DELPUECH_FOXO3_TARGETS_UP
            1  Delpuech Foxo3 Targets Dn                                    0.659947  0.0118856     1.79092  0.0450416  1                1852   0.605263  DELPUECH_FOXO3_TARGETS_DN
            5  Delpuech Foxo3 Targets Dn                                    0.609922  0.039992      1.64955  0.177575   1                2325   0.605263  DELPUECH_FOXO3_TARGETS_DN
           10  Delpuech Foxo3 Targets Dn                                    0.640701  0.0205138     1.73807  0.113573   1                2418   0.684211  DELPUECH_FOXO3_TARGETS_DN
            1  Wiederschain Targets Of Bmi1 And Pcgf2                      -0.565763  0.0403355    -1.62456  0.170938   1                8128   0.543478  WIEDERSCHAIN_TARGETS_OF_BMI1_AND_PCGF2
            1  Kasler Hdac7 Targets 2 Dn                                   -0.477951  0.0400814    -1.55288  0.202734   1                9900   0.37931   KASLER_HDAC7_TARGETS_2_DN
            1  Dormoy Elavl1 Targets                                        0.551238  0.0798787     1.44539  0.238115   1                2065   0.416667  DORMOY_ELAVL1_TARGETS
            1  Bilanges Rapamycin Sensitive Via Tsc1 And Tsc2               0.485559  0.123433      1.45019  0.238115   1                2309   0.484375  BILANGES_RAPAMYCIN_SENSITIVE_VIA_TSC1_AND_TSC2
            1  Bilanges Serum Sensitive Via Tsc1                           -0.560185  0.0533609    -1.53004  0.214421   1               10878   0.285714  BILANGES_SERUM_SENSITIVE_VIA_TSC1
            1  Kim Tial1 Targets                                            0.580634  0.0177814     1.68689  0.0852171  1                 873   0.428571  KIM_TIAL1_TARGETS
            4  Kim Tial1 Targets                                           -0.578032  0.0169053    -1.69392  0.19736    1                9725   0.428571  KIM_TIAL1_TARGETS
           10  Kim Tial1 Targets                                            0.567196  0.0240079     1.65875  0.164959   1                1131   0.428571  KIM_TIAL1_TARGETS
            1  Duan Prdm5 Targets                                          -0.439239  0.044991     -1.52071  0.220462   1                8515   0.380952  DUAN_PRDM5_TARGETS
            1  Dalessio Tsa Response                                       -0.668168  0.0844463    -1.51681  0.222374   1               10299   0.384615  DALESSIO_TSA_RESPONSE
            1  Winzen Degraded Via Khsrp                                   -0.608001  0.00423729   -1.89771  0.119811   1                8663   0.609375  WINZEN_DEGRADED_VIA_KHSRP
            1  Yuan Znf143 Partners                                         0.639028  0.0101833     1.7367   0.0638215  1                1354   0.333333  YUAN_ZNF143_PARTNERS
            5  Yuan Znf143 Partners                                         0.575909  0.049435      1.56938  0.239949   1                1960   0.52381   YUAN_ZNF143_PARTNERS
            1  Liu Topbp1 Targets                                          -0.527365  0.0520537    -1.47931  0.242366   1                9146   0.454545  LIU_TOPBP1_TARGETS
            1  Schmidt Por Targets In Limb Bud Up                           0.707484  0.0159713     1.75758  0.0547737  1                2416   0.761905  SCHMIDT_POR_TARGETS_IN_LIMB_BUD_UP
           10  Schmidt Por Targets In Limb Bud Up                           0.631974  0.0692292     1.56195  0.248199   1                1286   0.47619   SCHMIDT_POR_TARGETS_IN_LIMB_BUD_UP
            1  Katsanou Elavl1 Targets Up                                  -0.411478  0.0611043    -1.49917  0.231229   1                8317   0.416667  KATSANOU_ELAVL1_TARGETS_UP
            1  Juban Targets Of Spi1 And Fli1 Dn                            0.513234  0.00325071    1.81342  0.0381776  1                2084   0.475     JUBAN_TARGETS_OF_SPI1_AND_FLI1_DN
            1  Servitja Islet Hnf1A Targets Up                             -0.557471  0.043629     -1.6654   0.162484   1                9086   0.487179  SERVITJA_ISLET_HNF1A_TARGETS_UP
            1  Servitja Liver Hnf1A Targets Up                              0.416355  0.0486097     1.53102  0.183101   1                 326   0.212766  SERVITJA_LIVER_HNF1A_TARGETS_UP
            5  Servitja Liver Hnf1A Targets Up                              0.435947  0.031793      1.59765  0.21017    1                 564   0.234043  SERVITJA_LIVER_HNF1A_TARGETS_UP
           10  Servitja Liver Hnf1A Targets Up                              0.47102   0.0100301     1.72232  0.118647   1                 825   0.276596  SERVITJA_LIVER_HNF1A_TARGETS_UP
            1  Servitja Liver Hnf1A Targets Dn                             -0.440993  0.0232464    -1.58488  0.189446   1                8836   0.384615  SERVITJA_LIVER_HNF1A_TARGETS_DN
            1  Kohoutek Ccnt1 Targets                                      -0.473164  0.0563636    -1.53299  0.214421   1                7692   0.487805  KOHOUTEK_CCNT1_TARGETS
            1  Pedersen Metastasis By Erbb2 Isoform 1                      -0.553076  0.0881764    -1.52958  0.214421   1                8713   0.527778  PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_1
            1  Pedersen Metastasis By Erbb2 Isoform 6                      -0.653812  0.00184162   -1.83716  0.119811   1                9164   0.590909  PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_6
            1  Yang Bcl3 Targets Up                                         0.360512  0.0495069     1.48309  0.217553   1                 421   0.166008  YANG_BCL3_TARGETS_UP
            1  Le Neuronal Differentiation Dn                               0.784843  0.00742078    1.79364  0.0440822  1                 781   0.642857  LE_NEURONAL_DIFFERENTIATION_DN
            5  Le Neuronal Differentiation Dn                               0.826643  0.00258501    1.87778  0.0582517  1                 623   0.714286  LE_NEURONAL_DIFFERENTIATION_DN
           10  Le Neuronal Differentiation Dn                               0.691202  0.0564339     1.57737  0.231959   1                1584   0.642857  LE_NEURONAL_DIFFERENTIATION_DN
            1  Plasari Tgfb1 Targets 1Hr Up                                -0.713197  0.0151945    -1.78546  0.131106   1                8277   0.851852  PLASARI_TGFB1_TARGETS_1HR_UP
            1  Plasari Tgfb1 Targets 10Hr Up                               -0.519098  0.0246964    -1.71341  0.146665   1                9036   0.496454  PLASARI_TGFB1_TARGETS_10HR_UP
            1  Plasari Nfic Targets Basal Up                               -0.595249  0.0385086    -1.56059  0.200146   1                8520   0.642857  PLASARI_NFIC_TARGETS_BASAL_UP
            1  Plasari Tgfb1 Signaling Via Nfic 10Hr Up                    -0.475458  0.0834165    -1.46718  0.249633   1                8463   0.511628  PLASARI_TGFB1_SIGNALING_VIA_NFIC_10HR_UP
            1  Plasari Tgfb1 Signaling Via Nfic 10Hr Dn                    -0.572144  0.0366513    -1.62972  0.170938   1                9036   0.608696  PLASARI_TGFB1_SIGNALING_VIA_NFIC_10HR_DN
            1  Wang Mll Targets                                            -0.537007  0.00649219   -1.86252  0.119811   1                9016   0.5       WANG_MLL_TARGETS
            1  Kang Glis3 Targets                                          -0.610594  0.0453258    -1.5368   0.211252   1                9852   0.4       KANG_GLIS3_TARGETS
            1  Delacroix Rarg Bound Mef                                    -0.332244  0.0343769    -1.47447  0.245513   1                8738   0.339041  DELACROIX_RARG_BOUND_MEF
            1  Delacroix Rar Targets Up                                    -0.497754  0.0145313    -1.67459  0.161802   1                8801   0.461538  DELACROIX_RAR_TARGETS_UP
            1  Liu Il13 Priming Model                                      -0.792772  0.000408747  -2.00736  0.116729   1                9613   0.615385  LIU_IL13_PRIMING_MODEL
            1  Fu Interact With Alkbh8                                      0.730522  0.0230646     1.64507  0.106001   1                1834   0.769231  FU_INTERACT_WITH_ALKBH8
            4  Fu Interact With Alkbh8                                     -0.781067  0.00595356   -1.75518  0.192508   1                9809   0.769231  FU_INTERACT_WITH_ALKBH8
            1  Reichert Mitosis Lin9 Targets                                0.895955  0.000401284   1.92064  0.0188436  1                 218   0.730769  REICHERT_MITOSIS_LIN9_TARGETS
            4  Reichert Mitosis Lin9 Targets                               -0.780148  0.0407393    -1.65186  0.220221   1               10438   0.692308  REICHERT_MITOSIS_LIN9_TARGETS
            5  Reichert Mitosis Lin9 Targets                                0.890097  0.000789733   1.89202  0.0555334  1                 470   0.769231  REICHERT_MITOSIS_LIN9_TARGETS
           10  Reichert Mitosis Lin9 Targets                                0.814495  0.014979      1.74247  0.112869   1                1165   0.769231  REICHERT_MITOSIS_LIN9_TARGETS
            1  Pedrioli Mir31 Targets Dn                                   -0.411002  0.0103606    -1.65098  0.168175   1                7341   0.528517  PEDRIOLI_MIR31_TARGETS_DN
            1  Phong Tnf Targets Up                                        -0.633943  0.029203     -1.7271   0.144779   1                8888   0.509804  PHONG_TNF_TARGETS_UP
            1  Phong Tnf Response Via P38 Partial                          -0.488479  0.0333803    -1.66574  0.162484   1                8888   0.434109  PHONG_TNF_RESPONSE_VIA_P38_PARTIAL
            1  Yu Bap1 Targets                                              0.69255   0.00788835    1.86229  0.0278677  1                 510   0.48      YU_BAP1_TARGETS
            4  Yu Bap1 Targets                                             -0.611178  0.0432282    -1.62861  0.230529   1                9255   0.56      YU_BAP1_TARGETS
            5  Yu Bap1 Targets                                              0.716916  0.00280112    1.89351  0.0555334  1                 838   0.52      YU_BAP1_TARGETS
           10  Yu Bap1 Targets                                              0.585463  0.0719164     1.55927  0.249257   1                1620   0.52      YU_BAP1_TARGETS
            1  Vanoevelen Myogenesis Sin3A Targets                          0.347911  0.0129045     1.52998  0.183101   1                2295   0.331522  VANOEVELEN_MYOGENESIS_SIN3A_TARGETS
            1  Abramson Interact With Aire                                  0.752238  0.00101215    1.94436  0.0162816  1                1919   0.72093   ABRAMSON_INTERACT_WITH_AIRE
            5  Abramson Interact With Aire                                  0.69295   0.0076582     1.78237  0.0972596  1                2070   0.674419  ABRAMSON_INTERACT_WITH_AIRE
           10  Abramson Interact With Aire                                  0.683501  0.0110932     1.76842  0.0988582  1                2424   0.651163  ABRAMSON_INTERACT_WITH_AIRE
            1  Roessler Liver Cancer Metastasis Up                          0.372749  0.035461      1.46206  0.230639   1                1652   0.267442  ROESSLER_LIVER_CANCER_METASTASIS_UP
           10  Roessler Liver Cancer Metastasis Up                          0.435542  0.00427263    1.70156  0.130517   1                1998   0.313953  ROESSLER_LIVER_CANCER_METASTASIS_UP
            1  Holleman Prednisolone Resistance B All Up                    0.584444  0.0440073     1.57359  0.155488   1                3181   0.75      HOLLEMAN_PREDNISOLONE_RESISTANCE_B_ALL_UP
            1  Acevedo Fgfr1 Targets In Prostate Cancer Model Dn           -0.45362   0.0469502    -1.60166  0.181764   1                8138   0.532609  ACEVEDO_FGFR1_TARGETS_IN_PROSTATE_CANCER_MODEL_DN
            1  Anastassiou Cancer Mesenchymal Transition Signature         -0.765711  0.114794     -1.46724  0.249633   1                9840   0.770492  ANASTASSIOU_CANCER_MESENCHYMAL_TRANSITION_SIGNATURE
            1  Lim Mammary Luminal Progenitor Dn                           -0.640927  0.0414737    -1.58255  0.190354   1                8059   0.666667  LIM_MAMMARY_LUMINAL_PROGENITOR_DN
            1  Lim Mammary Luminal Mature Dn                               -0.604276  0.0847936    -1.57209  0.195207   1                8384   0.575     LIM_MAMMARY_LUMINAL_MATURE_DN
            1  Durand Stroma S Up                                          -0.506454  0.0101112    -1.81825  0.119811   1                8410   0.480952  DURAND_STROMA_S_UP
            1  Durand Stroma Ns Up                                         -0.390758  0.0270051    -1.54651  0.203357   1                7891   0.5       DURAND_STROMA_NS_UP
            1  Smirnov Response To Ir 2Hr Up                               -0.548116  0.0168836    -1.73347  0.141971   1                9051   0.555556  SMIRNOV_RESPONSE_TO_IR_2HR_UP
            1  Smirnov Response To Ir 6Hr Up                               -0.423105  0.00588593   -1.71862  0.146665   1                8569   0.404255  SMIRNOV_RESPONSE_TO_IR_6HR_UP
            1  Smirnov Response To Ir 6Hr Dn                                0.465019  0.0611055     1.53399  0.181681   1                 564   0.235955  SMIRNOV_RESPONSE_TO_IR_6HR_DN
            5  Smirnov Response To Ir 6Hr Dn                                0.478995  0.0519557     1.56933  0.239949   1                1830   0.359551  SMIRNOV_RESPONSE_TO_IR_6HR_DN
            1  Ghandhi Direct Irradiation Dn                               -0.549892  0.0294835    -1.61615  0.175815   1                9479   0.434783  GHANDHI_DIRECT_IRRADIATION_DN
            1  Warters Response To Ir Skin                                 -0.511391  0.0259109    -1.69095  0.152216   1                9071   0.462963  WARTERS_RESPONSE_TO_IR_SKIN
            1  Warters Ir Response 5Gy                                     -0.564787  0.019548     -1.71055  0.146665   1                9171   0.516129  WARTERS_IR_RESPONSE_5GY
            1  Zhou Cell Cycle Genes In Ir Response 6Hr                     0.849615  0.000201613   1.9672   0.0145844  0.536492          660   0.787879  ZHOU_CELL_CYCLE_GENES_IN_IR_RESPONSE_6HR
            4  Zhou Cell Cycle Genes In Ir Response 6Hr                    -0.703082  0.0585347    -1.62635  0.230529   1                9689   0.651515  ZHOU_CELL_CYCLE_GENES_IN_IR_RESPONSE_6HR
            5  Zhou Cell Cycle Genes In Ir Response 6Hr                     0.789637  0.00478469    1.81745  0.0793115  1                 776   0.712121  ZHOU_CELL_CYCLE_GENES_IN_IR_RESPONSE_6HR
           10  Zhou Cell Cycle Genes In Ir Response 6Hr                     0.826433  0.00120603    1.91201  0.056431   1                 705   0.757576  ZHOU_CELL_CYCLE_GENES_IN_IR_RESPONSE_6HR
            1  Zhou Cell Cycle Genes In Ir Response 24Hr                    0.802467  0.000402091   2.00656  0.0127991  1                 649   0.701923  ZHOU_CELL_CYCLE_GENES_IN_IR_RESPONSE_24HR
            4  Zhou Cell Cycle Genes In Ir Response 24Hr                   -0.709506  0.022763     -1.77102  0.192508   1               10027   0.576923  ZHOU_CELL_CYCLE_GENES_IN_IR_RESPONSE_24HR
            5  Zhou Cell Cycle Genes In Ir Response 24Hr                    0.744001  0.00676348    1.84728  0.0655189  1                 940   0.663462  ZHOU_CELL_CYCLE_GENES_IN_IR_RESPONSE_24HR
           10  Zhou Cell Cycle Genes In Ir Response 24Hr                    0.76165   0.00259119    1.91282  0.056431   1                 854   0.673077  ZHOU_CELL_CYCLE_GENES_IN_IR_RESPONSE_24HR
            1  Zwang Class 2 Transiently Induced By Egf                    -0.575127  0.0670855    -1.55084  0.202734   1                8948   0.53125   ZWANG_CLASS_2_TRANSIENTLY_INDUCED_BY_EGF
            1  Zwang Class 3 Transiently Induced By Egf                    -0.520241  0.0243608    -1.74449  0.139914   1                8391   0.502762  ZWANG_CLASS_3_TRANSIENTLY_INDUCED_BY_EGF
            1  Zwang Egf Interval Dn                                       -0.399561  0.00569453   -1.65563  0.164295   1                7582   0.505952  ZWANG_EGF_INTERVAL_DN
            1  Chemello Soleus Vs Edl Myofibers Up                         -0.615911  0.0766467    -1.48223  0.241988   1                8165   0.6       CHEMELLO_SOLEUS_VS_EDL_MYOFIBERS_UP
            1  Horton Srebf Targets                                         0.569     0.0904977     1.49142  0.213453   1                3570   0.727273  HORTON_SREBF_TARGETS
            1  Farmer Breast Cancer Cluster 4                              -0.897514  0.0310171    -1.53017  0.214421   1                9967   1         FARMER_BREAST_CANCER_CLUSTER_4
            1  Eppert Hsc R                                                -0.507481  0.0034219    -1.80132  0.124176   1                9044   0.421053  EPPERT_HSC_R
            1  Eppert Progenitor                                            0.524476  0.017551      1.72752  0.0667969  1                2005   0.451327  EPPERT_PROGENITOR
            5  Eppert Progenitor                                            0.542226  0.0122024     1.78417  0.0972596  1                1571   0.39823   EPPERT_PROGENITOR
           10  Eppert Progenitor                                            0.487038  0.0472582     1.60377  0.207217   1                2284   0.442478  EPPERT_PROGENITOR
            1  Eppert Ce Hsc Lsc                                           -0.623717  0.00724346   -1.77792  0.134462   1                9058   0.62069   EPPERT_CE_HSC_LSC

</div>

Pathways (c2.cp) from MSigDB.


```python
plot_ds(df_cp, fdr=0.05, le_prop=0.0)
```

![](figures/gsea-fa_cp-es-plot-a_1){#cp-es-plot-a }\



```python
plot_ds(df_cp, fdr=0.05, le_prop=0.8)
```

![](figures/gsea-fa_cp-es-plot-b_1){#cp-es-plot-b }\



```python
table_ds(df_cp, fdr=0.05)
```


<div class="datatable">  mri_feature  gene_set                                                              es            p      nes        fdr      fwer    max_es_at    le_prop  gene_set_code
-------------  --------------------------------------------------------------  --------  -----------  -------  ---------  --------  -----------  ---------  --------------------------------------------------------------
            1  Kegg Dna Replication                                            0.785845  0.000509892  2.18994  0.0396865  0.560881          442   0.625     KEGG_DNA_REPLICATION
            1  Kegg Mismatch Repair                                            0.781049  0.000102902  2.05903  0.0396865  0.113192          441   0.5       KEGG_MISMATCH_REPAIR
            1  Pid Plk1 Pathway                                                0.757247  0.000100513  2.07507  0.0396865  0.110564          199   0.394737  PID_PLK1_PATHWAY
            1  Pid Foxm1 Pathway                                               0.734502  0.00100796   1.97573  0.040293   1                1370   0.575758  PID_FOXM1_PATHWAY
            1  Reactome Activation Of The Pre Replicative Complex              0.755031  0.00324166   1.99102  0.040293   1                 497   0.590909  REACTOME_ACTIVATION_OF_THE_PRE_REPLICATIVE_COMPLEX
            1  Reactome Pol Switching                                          0.862758  0.000784138  2.01496  0.040293   0.862552          442   0.818182  REACTOME_POL_SWITCHING
            1  Reactome Repair Synthesis For Gap Filling By Dna Pol In Tc Ner  0.830559  0.000648368  1.98676  0.040293   0.713205          441   0.666667  REACTOME_REPAIR_SYNTHESIS_FOR_GAP_FILLING_BY_DNA_POL_IN_TC_NER
            1  Reactome Kinesins                                               0.812804  0.00074381   1.96325  0.0407759  0.818191          281   0.5       REACTOME_KINESINS
            1  Reactome Lagging Strand Synthesis                               0.836409  0.000318167  2.10387  0.0396865  0.349984          442   0.705882  REACTOME_LAGGING_STRAND_SYNTHESIS
            1  Reactome Dna Replication                                        0.637197  0.00059994   1.94291  0.0440208  0.659934         2021   0.558442  REACTOME_DNA_REPLICATION
            1  Reactome Activation Of Atr In Response To Replication Stress    0.751141  0.00132329   2.04387  0.0396865  1                 497   0.517241  REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS
            1  Reactome Unwinding Of Dna                                       0.956177  0.000504032  2.11311  0.0396865  0.554435          497   1         REACTOME_UNWINDING_OF_DNA
            1  Reactome Telomere Maintenance                                   0.710926  0.000705219  1.97457  0.040293   0.77574           665   0.515152  REACTOME_TELOMERE_MAINTENANCE
            1  Reactome Extension Of Telomeres                                 0.783725  0.00103896   2.06531  0.0396865  1                 665   0.681818  REACTOME_EXTENSION_OF_TELOMERES
            1  Reactome G2 M Checkpoints                                       0.771507  0.000607103  2.12886  0.0396865  0.667813          497   0.529412  REACTOME_G2_M_CHECKPOINTS
            1  Reactome Dna Strand Elongation                                  0.884704  0.000105485  2.38135  0.0396865  0.116034          497   0.814815  REACTOME_DNA_STRAND_ELONGATION

</div>
