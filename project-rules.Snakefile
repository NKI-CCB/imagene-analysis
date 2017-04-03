import os
from pathlib import Path

import click
import dotenv
import requests

# Configuration of executable names etc.
configfile: "config/snakemake.yaml"

# Credentials that need top be kept out of source control
dotenv.load_dotenv(str(Path("./.env").resolve()))
beehub_username = os.environ.get("BEEHUB_USERNAME")
beehub_password = os.environ.get("BEEHUB_PASSWORD")

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
