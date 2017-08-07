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
    "data/processed/mri-features.nc",
    "data/processed/mri-eigenbreasts.nc",
    "data/processed/clinical.nc",
]


#------------------
# Download Raw Data

def download_beehub(file_location, out):
    file_url = 'https://beehub.nl/' + file_location
    out = Path(out).resolve()

    r = requests.get(file_url, stream=True,
                            auth=(beehub_username, beehub_password))
    r.raise_for_status()

    r_length = r.headers.get('content-length')
    if r_length is not None:
        r_length = int(r_length)
    chunk_size = 1024
    with out.open('wb') as file:
        with click.progressbar(r.iter_content(1024),
                               length=int(r_length/chunk_size)) as bar:
            for chunk in bar:
                file.write(chunk)


rule download_gene_expression:
    output: "data/raw/gene-expression.nc"
    run:
        download_beehub(
            "home/tychobismeijer/Imagene/gene_expression/"
            "2017-02-01-gene-expression-imagene.nc",
            output[0]
        )

rule download_sample_annotation:
    output: "data/raw/sample-tracking.tsv"
    run:
        download_beehub(
            "home/tychobismeijer/Imagene/"
            "gene_expression/2016-06-14-sample-tracking.tsv",
            output[0]
        )

rule download_mri_features:
    output: "data/raw/mri-features.xlsx"
    run:
        download_beehub(
            "home/tychobismeijer/Imagene/mri/"
            "2016-03-31-Tumor_Parenchym_Features_variablenamesupdated.xlsx",
            output[0]
        )

rule download_eigenbreasts:
    output: "data/raw/eigenbreasts_{subset}.xlsx"
    run:
        download_beehub(
            "home/tychobismeijer/Imagene/mri/2017-07-18-eigenbreasts/"
            "{}_ero4.xlsx".format(str(Path(output[0]).stem)),
            output[0],
        )

rule download_clinical_data:
    output: "data/raw/imagene_clinical.tsv"
    run:
        download_beehub(
            "home/tychobismeijer/Imagene/clinical/"
            "2016-01-19-imagene_clinical.tsv",
            output[0]
        )

rule download_sfa_tcga_breast:
    output: "data/external/tcga-breast-gexp+rppa+cn-sfa-solution.h5"
    run:
        download_beehub(
            "home/tychobismeijer/Imagene/external/"
            "tcga-breast-gexp+rppa+cn-sfa-solution.h5",
            output[0])


rule download_sfa_tcga_breast_conf:
    output: "data/external/tcga-breast-gexp+rppa+cn-sfa-solution.yaml"
    run:
        download_beehub(
            "home/tychobismeijer/Imagene/external/"
            "tcga-breast-gexp+rppa+cn-sfa-solution.yaml",
            output[0])


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
        "data/processed/mri-features.nc"
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
        gene_annot="data/external/ensembl_annotation.tsv",
    output:
        "data/processed/gene-expression.nc"
    shell:
        "{config[python]} {input.script} {input.gexp} {input.sample_tracking} "
        "{input.gene_annot} {output}"

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


########################################################################
# FEATURES                                                             #
########################################################################

all_targets['features'] = [
    "data/processed/mri-features.nc",
    "data/processed/mri-features-reg-volume.nc",
    "data/processed/mri-features-fa.nc",
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
    "models/sfa/sfa.nc",
]


rule apply_tcga_sfa:
    input:
        script="src/models/apply_sfa.py",
        gexp="data/processed/gene-expression.nc",
        conf="data/external/tcga-breast-gexp+rppa+cn-sfa-solution.yaml",
        sfa_tcga="data/external/tcga-breast-gexp+rppa+cn-sfa-solution.h5",
    output:
        "models/sfa/sfa.nc"
    shell:
        "{config[python]} {input.script} {input.gexp} {input.sfa_tcga} "
        "{input.conf} {output}"

rule cross_validate_mri_from_factors:
    input:
        script="src/models/cv_mri_from_factors.py",
        sfa="models/sfa/sfa.nc",
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
        sfa="models/sfa/sfa.nc",
    output:
        "models/factors_from_mri/performance.nc",
    shell:
        "{config[python]} {input.script} {input.mri} {input.sfa} "
        "{output}"



########################################################################
# ANALYSIS                                                             #
########################################################################

all_targets['analyses'] = expand(
    "analyses/gsea/{features}_{gene_set_abs}.nc",
    features=[
        "mri-features", "mri-features-fa", "mri-features-reg-volume",
        "mri-features-er", "mri-features-er-fa", "mri-features-er-reg-volume",
    ],
    gene_set_abs=["c2.cgp_F", "c2.cp_T", "h.all_T"],
)


rule analyse_gene_sets:
    input:
        script="src/analysis/analyse-gene-set-enrichment.R",
        gexp="data/processed/gene-expression.nc",
        mri="data/processed/{mri}.nc",
        gene_sets="data/external/msigdb/{gene_set}.v5.2.entrez.gmt",
    output:
        protected("analyses/gsea/{mri}_{gene_set}_{abs}.Rds"),
    threads:
        100
    shell:
        "mkdir -p analyses/gsea; "
        "{config[r]} {input.script} {input.gexp} {input.mri} "
        "{input.gene_sets} {output} --abs {wildcards.abs} --threads {threads} "
        "--perms 10000"

rule gene_set_analysis_to_netcdf:
    input:
        script="src/analysis/gsea-rds-to-nc.R",
        rds="analyses/gsea/{name}.Rds",
    output:
        "analyses/gsea/{name}.nc",
    shell:
        "{config[r]} {input.script} {input.rds} {output}"


########################################################################
# REPORTS                                                              #
########################################################################


all_targets['reports'] = ["reports/" + f for f in [
    str(p.with_suffix(".html"))
    for p in Path('notebooks').glob("*.pmd")
]]


gsea_deps = [
    "src/plot.py",
    "src/reports/es-heatmap-fun.py",
    "src/reports/load-gsea-fun.py",
    "src/reports/setup-matplotlib.py",
]

report_deps = {
    "gsea": [
        "analyses/gsea/mri-features_h.all_T.nc",
    ] + gsea_deps,
    "gsea-reg": [
        "analyses/gsea/mri-features-reg-volume_c2.cgp_F.nc",
        "analyses/gsea/mri-features-reg-volume_h.all_T.nc",
        "analyses/gsea/mri-features-reg-volume_c2.cp_T.nc",
    ] + gsea_deps,
    "gsea-fa": [
        "analyses/gsea/mri-features-fa_c2.cgp_F.nc",
        "analyses/gsea/mri-features-fa_h.all_T.nc",
        "analyses/gsea/mri-features-fa_c2.cp_T.nc",
    ] + gsea_deps,
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
}

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
    shell:
        "{config[nbconvert]} --to html --execute {input}"


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
