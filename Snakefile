from functools import reduce
from operator import add
import os
from pathlib import Path

import click
import dotenv
import requests

# Configuration of executable names etc.
configfile: "config/snakemake.yaml"

# Credentials that need top be kept out of source control
dotenv.load_dotenv(str(Path("./.env").resolve()))
beehub_username = os.environ["BEEHUB_USERNAME"]
beehub_password = os.environ["BEEHUB_PASSWORD"]

# Allow importing of Python modules from src directory
os.environ["PYTHONPATH"] = str(Path("./src/").resolve())


all_targets = dict()

rule all:
    input: lambda _: reduce(add, all_targets.values())
rule all_from_make:
    input: lambda _: reduce(add, all_targets.values())

rule all_data:
    input: lambda _: all_targets['data']
rule all_reports:
    input: lambda _: all_targets['reports']
rule all_features:
    input: lambda _: all_targets['features']
rule all_models:
    input: lambda _: all_targets['models']
rule all_analyses:
    input: lambda _: all_targets['analyses']
rule all_notebooks:
    input: lambda _: all_targets['notebooks']


########################################################################
# DATA                                                                 #
########################################################################


all_targets['data'] = [
    "data/processed/gene-expression.nc",
    "data/processed/mri-features-all.nc",
    "data/processed/mri-features-er.nc",
    "data/processed/mri-eigenbreasts.nc",
    "data/processed/clinical.nc",
]


#------------------
# Download Raw Data

# Project data on remote share #

def download_scp(remotefile, dest):
    root = config['download_root'] 
    shell(f"""
        scp {root}{remotefile} {dest}
        touch {dest}
    """)

download_funcs = {'download_scp': download_scp}
download = download_funcs[config['download_func']]

rule download_gene_expression:
    output: "data/raw/gene-expression.nc"
    params:
        file="gene_expression/2017-02-01-gene-expression-imagene.nc",
    run:
        download(params.file, output[0])

rule download_sample_annotation:
    output: "data/raw/sample-tracking.tsv"
    params:
        file="gene_expression/2016-06-14-sample-tracking.tsv"
    run:
        download(params.file, output[0])

rule download_mri_features:
    output: "data/raw/mri-features.xlsx"
    params:
        file="mri/2016-03-31-Tumor_Parenchym_Features_"
             "variablenamesupdated.xlsx",
    run:
        download(params.file, output[0])

rule download_eigenbreasts:
    output: "data/raw/eigenbreasts_{subset}.xlsx"
    params:
        file="mri/2017-08-09-eigenbreasts/eigenbreasts_Elastix_{subset}.xlsx"
    run:
        download(params.file, output[0])

rule download_clinical_data:
    output: "data/raw/imagene_clinical.tsv"
    params:
        file="clinical/2016-01-19-imagene_clinical.tsv"
    run:
        download(params.file, output[0])

rule download_sfa_tcga_breast:
    output: "data/external/tcga-breast-gexp+rppa+cn-sfa-solution.h5"
    params:
        file="external/tcga-breast-gexp+rppa+cn-sfa-solution.h5"
    run:
        download(params.file, output[0])

rule download_sfa_tcga_breast_conf:
    output: "data/external/tcga-breast-gexp+rppa+cn-sfa-solution.yaml"
    params:
        file="external/tcga-breast-gexp+rppa+cn-sfa-solution.yaml"
    run:
        download(params.file, output[0])


# External resources #

rule download_ensembl_annotation:
    input:
        script="src/data/query_ensembl_reference.py"
    output:
        "data/external/ensembl_annotation.tsv"
    shell:
        "{config[python]} {input.script} {output}"

rule download_msigdb:
    output:
        "data/external/msigdb/msigdb_v5.2_files_to_download_locally.zip"
    message:
        "Please download msigdb_v5.2_files_to_download_locally.zip from "
        "http://software.broadinstitute.org/gsea/downloads.jsp#msigdb "
        "under Archived Releases and place it into data/external/msigdb."

rule unzip_msigdb:
    input:
        "data/external/msigdb/msigdb_v5.2_files_to_download_locally.zip"
    output:
        "data/external/msigdb/{gene_set}.v5.2.{gene_ids}.gmt"
    shell:
        "unzip -p {input} "
        "msigdb_v5.2_files_to_download_locally/msigdb_v5.2_GMTs/"
        "{wildcards.gene_set}.v5.2.{wildcards.gene_ids}.gmt > {output}\n"
        "touch {output}"


#-------------
# Process Data


rule tcga_factors_to_tsv:
    input:
        "data/external/tcga-breast-gexp+rppa+cn-sfa-solution.h5",
    output:
        "data/to_share/tcga-breast-gexp+rppa+cn-sfa-factors.tsv",
    run:
        import h5py
        import xarray as xr
        import pandas as pd
        import numpy as np

        with h5py.File(input[0]) as f:
            samples = [s.decode() for s in f['sample names']]
            factor_names = ['Factor {}'.format(i+1)
                            for i in range(f['factors'].shape[1])]
            tcga_factors = xr.DataArray(
                data=np.array(f['factors']),
                dims=['Patient', 'factor'],
                coords={
                    'Patient': np.array(samples, 'object'),
                    'factor': np.array(factor_names, 'object'),
                },
            )
            tcga_factors.to_pandas().to_csv(output[0], sep='\t')

rule process_mri_features:
    input:
        script="src/data/process_mri.py",
        xlsx="data/raw/mri-features.xlsx",
    output:
        "data/processed/mri-features-all.nc"
    shell:
        "{config[python]} {input.script} {input.xlsx} {output} "
        "--study-nr-col=MARGINSstudyNr"

rule process_eigenbreasts:
    input:
        "src/data/process_mri_eigenbreasts.py",
        "data/raw/eigenbreasts_ipsi_ds2.xlsx",
        "data/raw/eigenbreasts_ipsi_ds4.xlsx",
        "data/raw/eigenbreasts_ipsi_ds8.xlsx",
        "data/raw/eigenbreasts_ipsi_ds16.xlsx",
        "data/raw/eigenbreasts_contra_ds2.xlsx",
        "data/raw/eigenbreasts_contra_ds4.xlsx",
        "data/raw/eigenbreasts_contra_ds8.xlsx",
        "data/raw/eigenbreasts_contra_ds16.xlsx",
        "data/raw/eigenbreasts_both_ds2.xlsx",
        "data/raw/eigenbreasts_both_ds4.xlsx",
        "data/raw/eigenbreasts_both_ds8.xlsx",
        "data/raw/eigenbreasts_both_ds16.xlsx",
    output:
        "data/processed/mri-eigenbreasts.nc"
    shell:
        "{config[python]} {input} {output} "
        "--study-nr-col=StudyID"

rule process_gene_expression:
    input:
        script="src/data/process_gene_expression.py",
        gexp="data/raw/gene-expression.nc",
        sample_tracking="data/raw/sample-tracking.tsv",
        gene_annot="data/ensembl_annotation.tsv",
    output:
        "data/processed/gene-expression.nc"
    shell:
        "{config[python]} {input.script} {input.gexp} {input.sample_tracking} "
        "{input.gene_annot} {output}"

rule process_gene_expression_voom:
    input:
        script="src/data/process_gene_expression_voom.py",
        gexp="data/processed/gene-expression.nc"
    output:
        "data/processed/gene-expression-voom.nc"
    shell:
        "{config[python]} {input.script} {input.gexp} {output}"

rule process_clincal:
    input:
        script="src/data/process_clinical.py",
        tsv="data/raw/imagene_clinical.tsv",
    output:
        "data/processed/clinical.nc"
    shell:
        "{config[python]} {input.script} {input.tsv} {output} "


rule select_er:
    input:
        script="src/data/select_samples.py",
        mri="data/processed/mri-features-all.nc",
        clinical="data/processed/clinical.nc",
    output:
        "data/processed/mri-features-er.nc",
    shell:
        "{config[python]} {input.script} {input.mri} --er-positive "
        "{input.clinical} {output}"

rule concat_data_for_sfa:
    input:
        script="src/data/concat_data_sfa.py",
        gexp="data/processed/gene-expression-voom.nc",
        mri="data/processed/mri-features-all.nc",
    output:
        "data/processed/concat-data.nc",
    shell:
        "{config[python]} {input.script} {input.gexp} {input.mri} {output}"


########################################################################
# FEATURES                                                             #
########################################################################

all_targets['features'] = [
    "data/processed/mri-features-all.nc",
    "data/processed/mri-features-all-reg-volume.nc",
    "data/processed/mri-features-all-fa.nc",
    "data/processed/mri-features-er.nc",
    "data/processed/mri-features-er-reg-volume.nc",
    "data/processed/mri-features-er-fa.nc",
]

rule regress_out_volume:
    input:
        script="src/features/regress_out_mri_var.py",
        mri="data/processed/mri-features-{subset}.nc",
    output:
        "data/processed/mri-features-{subset}-reg-volume.nc"
    shell:
        "{config[python]} {input.script} {input.mri} {output}"

rule factor_analysis_mri_features:
    input:
        script="src/features/fa_mri_features.py",
        mri="data/processed/mri-features-{subset}.nc",
    output:
        "data/processed/mri-features-{subset}-fa.nc"
    shell:
        "{config[python]} {input.script} 10 {input.mri} {output}"


########################################################################
# MODELS                                                               #
########################################################################

all_targets['models'] = [
    "models/sfa_tcga/sfa.nc",
]


rule apply_tcga_sfa:
    input:
        script="src/models/apply_sfa.py",
        gexp="data/processed/gene-expression.nc",
        conf="data/external/tcga-breast-gexp+rppa+cn-sfa-solution.yaml",
        sfa_tcga="data/external/tcga-breast-gexp+rppa+cn-sfa-solution.h5",
    output:
        "models/sfa_tcga/sfa.nc"
    shell:
        "{config[python]} {input.script} {input.gexp} {input.sfa_tcga} "
        "{input.conf} {output}"

rule cross_validate_mri_from_factors:
    input:
        script="src/models/cv_mri_from_factors.py",
        sfa="models/sfa_tcga/sfa.nc",
        mri="data/processed/mri-features.nc",
    output:
        "models/mri_from_factors/performance.nc",
    shell:
        "{config[python]} {input.script} {input.sfa} {input.mri} "
        "{output}"

rule cross_validate_factors_from_mri:
    input:
        script="src/models/cv_factors_from_mri.py",
        mri="data/processed/mri-features.nc",
        sfa="models/sfa_tcga/sfa.nc",
    output:
        "models/factors_from_mri/performance.nc",
    shell:
        "{config[python]} {input.script} {input.mri} {input.sfa} "
        "{output}"


rule run_sfa_gexp_sweep:
    """Sweep penalty on gene expression data to determine range"""
    input:
        script="src/models/run_sfa.py",
        data="data/processed/concat-data.nc",
    threads:
        64
    output:
        "models/sfa_mri_cad/gexp_sweep.nc"
    shell:
        "{config[python]} {input.script} {input.data} {output} "
        "--k 5 --l-gexp=-15:5:.1 --l-mri=-15 --alpha=0.7 "
        "--eps 1e-3 --max-iter 10000 "
        "--threads {threads}"

rule run_sfa_mri_sweep:
    """Sweep penalty on MRI CAD data to determine range"""
    input:
        script="src/models/run_sfa.py",
        data="data/processed/concat-data.nc",
    threads:
        64
    output:
        "models/sfa_mri_cad/mri_sweep.nc"
    shell:
        "{config[python]} {input.script} {input.data} {output} "
        "--k 5 --l-gexp=-20 --l-mri=-25:10 --alpha=0.7 "
        "--eps 1e-3 --max-iter 10000 "
        "--threads {threads}"

rule run_sfa_grid_search_coarse:
    """Sweep penalty on MRI CAD and RNA-seq"""
    input:
        script="src/models/run_sfa.py",
        data="data/processed/concat-data.nc",
    threads:
        128
    output:
        "models/sfa_mri_cad/parameter_sweep_coarse.nc"
    shell:
        "{config[python]} {input.script} {input.data} {output} "
        "--k 2:10 --l-gexp=-13:1:1 --l-mri=-18:3:1 --alpha=0.1:.9:.2 "
        "--eps 1e-3 --max-iter 10000 "
        "--threads {threads}"

rule run_sfa_grid_search_fine:
    """Sweep penalty on MRI CAD and RNA-seq"""
    input:
        script="src/models/run_sfa.py",
        data="data/processed/concat-data.nc",
    threads:
        128
    output:
        "models/sfa_mri_cad/parameter_sweep_fine.nc"
    shell:
        "{config[python]} {input.script} {input.data} {output} "
        "--k 3 --l-gexp=-8:3:.5 --l-mri=-18:3:.5 --alpha=0.1:.9:.2 "
        "--eps 1e-3 --max-iter 10000 "
        "--threads {threads}"

rule eval_sfa:
    input:
        script="src/models/eval_sfa_bic.py",
        models="models/sfa_mri_cad/{sweep}.nc",
        data="data/processed/concat-data.nc",
    output:
        "models/sfa_mri_cad/{sweep}-bics.nc",
    shell:
        "{config[python]} {input.script} {input.models} {input.data} {output}"

rule select_best_sfa:
    input:
        script="src/models/select_best_sfa.py",
        models="models/sfa_mri_cad/{sweep}.nc",
        bics="models/sfa_mri_cad/{sweep}-bics.nc",
    output:
        "models/sfa_mri_cad/{sweep}-best.nc",
    shell:
        "{config[python]} {input.script} {input.models} {input.bics} {output}"


########################################################################
# ANALYSIS                                                             #
########################################################################

mri_features = [
    "mri-features-all", "mri-features-all-fa", "mri-features-all-reg-volume",
    "mri-features-er", "mri-features-er-fa", "mri-features-er-reg-volume",
]

all_targets['analyses'] = expand(
    "analyses/gsea/{features}_{gene_set_abs}.nc",
    features=mri_features,
    gene_set_abs=["c2.cgp_F", "c2.cp_T", "h.all_T"],
) + expand(
    "analyses/de/{features}.nc",
    features=mri_features,
)

rule differential_expression_analysis:
    input:
        script="src/analysis/differential-expression.R",
        gexp="data/processed/gene-expression.nc",
        mri="data/processed/{mri}.nc",
    output:
        "analyses/de/{mri}.nc"
    shell:
        "mkdir -p analyses/de; "
        "{config[r]} {input.script} {input.gexp} {input.mri} {output}"

rule analyse_gene_sets:
    input:
        script="src/analysis/analyse-gene-set-enrichment.R",
        gexp="data/processed/gene-expression.nc",
        mri="data/processed/{mri}.nc",
        gene_sets="data/external/msigdb/{gene_set}.v5.2.entrez.gmt",
    output:
        protected("analyses/gsea/{mri}_{gene_set}_{abs}.Rds"),
    threads:
        4 # Takes a lot of memory
    shell:
        "mkdir -p analyses/gsea; "
        "{config[r]} {input.script} {input.gexp} {input.mri} "
        "{input.gene_sets} {output} --abs {wildcards.abs} --threads {threads} "
        "--perms 10000"

rule gene_set_analysis_to_netcdf:
    input:
        script="src/analysis/gsea-rds-to-nc.R",
        rds="analyses/gsea{a}/{name}.Rds",
    output:
        "analyses/gsea{a,.*}/{name}.nc",
    shell:
        "{config[r]} {input.script} {input.rds} {output}"

rule gene_set_analysis_to_xlsx:
    input:
        script="src/analysis/gsea-rds-to-xlsx.R",
        rds="analyses/gsea{a}/{name}.Rds",
    output:
        "analyses/gsea{a,.*}/{name}.xlsx",
    shell:
        "{config[r]} {input.script} {input.rds} {output} .25 0.0"

rule analyse_gene_sets_sfa:
    input:
        script="src/analysis/analyse-gene-set-enrichment-sfa.R",
        gexp="data/processed/gene-expression.nc",
        model="models/sfa_mri_cad/{sweep}-best.nc",
        gene_sets="data/external/msigdb/{gene_set}.v5.2.entrez.gmt",
    output:
        protected("analyses/gsea-sfa/{sweep}_{gene_set}_{abs}.Rds"),
    threads:
        16
    shell:
        "mkdir -p analyses/gsea; "
        "{config[r]} {input.script} {input.gexp} {input.model} "
        "{input.gene_sets} {output} --abs {wildcards.abs} --threads {threads} "
        "--perms 10000"

n_pc_eigenbreasts = {
    'contra_ds8': 12,
}

rule analyse_gene_sets_eigenbreasts:
    input:
        script="src/analysis/analyse-gene-set-enrichment-eigenbreasts.R",
        gexp="data/processed/gene-expression.nc",
        mri="data/processed/mri-eigenbreasts.nc",
        gene_sets="data/external/msigdb/{gene_set}.v5.2.entrez.gmt",
    output:
        protected("analyses/gsea/eigenbreasts_{mri_var}_{gene_set}_{abs}.Rds"),
    params:
        n_pc=lambda w: n_pc_eigenbreasts[w.mri_var]
    threads:
        4 # Takes a lot of memory
    shell:
        "mkdir -p analyses/gsea; "
        "OPENBLAS_NUM_THREADS=1 {config[r]} {input.script} {input.gexp} "
        "{input.mri} {wildcards.mri_var} {params.n_pc} "
        "{input.gene_sets} {output} --abs {wildcards.abs} --threads {threads} "
        "--perms 10000"


########################################################################
# REPORTS                                                              #
########################################################################


all_targets['reports'] = [
    str(p.with_suffix(".html"))
    for p in Path('reports').glob("*.pmd")
]

report_deps = {
    "mri-remove-size": [
        "src/plot.py",
        "src/reports/setup-matplotlib.py",
        "data/processed/mri-features.nc",
    ],
    "cv-mri-from-factors": [
        "src/plot.py",
        "src/reports/setup-matplotlib.py",
        "models/mri_from_factors/performance.nc",
    ],
    "cv-factors-from-mri": [
        "src/plot.py",
        "src/reports/setup-matplotlib.py",
        "models/factors_from_mri/performance.nc",
    ],
    "er-factor-correlations": [
        "src/util.py",
        "src/plot.py",
        "src/reports/setup-matplotlib.py",
        "data/external/set-index.tsv",
    ],
    "sfa-mg-parameter-sweeps": [
        "src/plot.py",
        "src/reports/setup-matplotlib.py",
        "models/sfa_mri_cad/gexp_sweep-bics.nc",
        "models/sfa_mri_cad/mri_sweep-bics.nc",
    ]
}

features_to_report_name = {
    "mri-features-all": "gsea",
    "mri-features-all-fa": "gsea-fa",
    "mri-features-all-reg-volume": "gsea-reg",
    "mri-features-er": "gsea-er",
    "mri-features-er-fa": "gsea-er-fa",
    "mri-features-er-reg-volume": "gsea-er-reg",
}

for mri_f in mri_features:
    report_deps[features_to_report_name[mri_f]] = [
        f"analyses/gsea/{mri_f}_c2.cgp_F.nc",
        f"analyses/gsea/{mri_f}_h.all_T.nc",
        f"analyses/gsea/{mri_f}_c2.cp_T.nc",
        f"analyses/de/{mri_f}.nc",
        "src/plot.py",
        "src/reports/es-heatmap-fun.py",
        "src/reports/load-gsea-fun.py",
        "src/reports/setup-matplotlib.py",
    ]


rule weave_report:
    input:
        lambda w: report_deps.get(w['report'], []),
        pmd="reports/{report}.pmd",
    output:
        "reports/{report}.md",
    shell:
        "{config[pweave]} "
        "--kernel=python3 "
        "-f pandoc "
        "--input-format=markdown "
        "{input.pmd} -o {output}"

rule markdown_to_html:
    input:
        "reports/pandoc-template.html",
        md="reports/{report}.md",
    output:
        "reports/{report}.html"
    shell:
        "{config[pandoc]} "
        "-t html5 "
        "--smart "
        "--standalone "
        "--mathjax "
        "--template=reports/pandoc-template.html "
        "--dpi=300 "
        "--default-image-extension=png "
        "--toc "
        "--highlight-style pygments "
        "--section-divs "
        "--filter pandoc-sidenote "
        "{input.md} -o {output}"


all_targets['notebooks'] = [
    str(p.with_suffix(".html"))
    for p in Path('notebooks').glob("*.ipynb")
]


rule convert_notebook_to_html:
    input:
        "notebooks/{notebook}.ipynb",
    output:
        "notebooks/{notebook}.html",
    priority: -10
    shell:
        "{config[nbconvert]} "
        "--ExecutePreprocessor.timeout=1800 "
        "--to html "
        "--execute {input}"


########################################################################
# FIGURES                                                              #
########################################################################

rule figure_mri_cad_correlation:
    input:
        script="src/visualization/figure-mri-cad-correlation.py",
        cad_features="data/processed/mri-features.nc",
    output: "reports/figures/mri-cad-correlation.svg"
    shell:
        "{config[python]} {input.script} {input.cad_features} {output}"


########################################################################
# INTERACTIVE                                                          #
########################################################################


rule run_notebook:
    shell:
        "jupyter notebook"

rule run_lab:
    shell:
        "jupyter lab"

rule run_ipython:
    shell:
        "ipython3"
