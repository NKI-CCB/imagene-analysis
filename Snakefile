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

# Other environment
os.environ["PYTHONPATH"] = str(Path("./src/").resolve())


########################################################################
# DATA                                                                 #
########################################################################

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


rule all_data:
    input:
        "data/processed/gene-expression.nc",
        "data/processed/mri-features.nc",
        "data/raw/imagene_clinical.tsv"

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
            "home/tychobismeijer/Imagene/mri/2017-07-12-eigenbreasts/" +
            str(Path(output[0]).name),
            output[0]
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
        "{wildcards.gene_set}.v5.2.{wildcards.gene_ids}.gmt > {output}"


#-------------
# Process Data


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

########################################################################
# FEATURES                                                             #
########################################################################

rule all_features:
    input:
        "data/processed/mri-features.nc",
        "data/processed/mri-features-reg-volume.nc",

rule regress_out_volume:
    input:
        script="src/features/regress_out_mri_var.py",
        mri="data/processed/mri-features.nc",
    output:
        "data/processed/mri-features-reg-volume.nc"
    shell:
        "{config[python]} {input.script} {input.mri} {output}"

########################################################################
# MODELS                                                               #
########################################################################

rule all_models:
    input:
        "models/sfa.nc"


rule apply_tcga_sfa:
    input:
        script="src/models/apply_sfa.py",
        gexp="data/processed/gene-expression.nc",
        sfa_tcga="data/external/tcga-breast-gexp+rppa+cn-sfa-solution.h5",
    output:
        "models/sfa.nc"
    shell:
        "{config[python]} {input.script} {input.gexp} {input.sfa_tcga} "
        "{output}"

########################################################################
# ANALYSIS                                                             #
########################################################################

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

report_deps = {
    "gsea": [
        "analyses/gsea/mri-features_h.all_T.nc",
        "analyses/gsea/mri-features-reg-volume_h.all_T.nc",
        "src/plot.py",
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
        "--toc "
        "--highlight-style pygments "
        "--section-divs "
        "{input.md} -o {output}"

rule run_notebook:
    shell:
        "jupyter notebook"

rule run_lab:
    shell:
        "jupyter lab"

rule run_ipython:
    shell:
        "ipython3"
