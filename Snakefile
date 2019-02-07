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
rule all_figures:
    input: lambda _: all_targets['figures']
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

rule download_clinical_data:
    output: "data/raw/imagene_clinical.tsv"
    params:
        file="clinical/2016-01-19-imagene_clinical.tsv"
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

rule download_zwart2011_er_signature:
    output: "data/external/zwart2011/table_S2.xls"
    params:
        url="http://emboj.embopress.org/content/embojnl/30/23/4764/DC3/embed/"
            "inline-supplementary-material-3.xls?download=true"

    shell:
        "mkdir -p data/external/zwart2011/\n"
        "curl {params.url} -o {output}"

rule download_refseq_annotation:
    output: "data/external/refseq/rna.gbk.gz"
    params:
        url="ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/"
            "GRCh38.p10_interim_annotation/interim_GRCh38.p10_rna.gbk.gz",
    shell:
        "mkdir -p data/external/refseq/\n"
        "curl {params.url} -o {output}"


#-------------
# Process Data

rule process_genbank:
    input:
        genbank_flatfile="data/external/refseq/rna.gbk.gz",
        script="src/data/parse_genbank_flatfile.py",
    output: "data/external/refseq/rna_annotation.tsv"
    shell:
        "{config[python]} {input.script} {input.genbank_flatfile} {output}"

rule process_zwart2011_signature:
    input:
        xls="data/external/zwart2011/table_S2.xls",
        refseq="data/external/refseq/rna_annotation.tsv",
        ensembl="data/external/ensembl_annotation.tsv",
        script="src/data/map_genes_zwart2011.py",
    output: "data/external/zwart2011/er_responsive_genes.tsv",
    shell:
        "{config[python]} {input.script} {input.xls} {input.refseq} "
        "{input.ensembl} {output}"


rule process_mri_features:
    input:
        script="src/data/process_mri.py",
        xlsx="data/raw/mri-features.xlsx",
    output:
        "data/processed/mri-features-all.nc"
    shell:
        "{config[python]} {input.script} {input.xlsx} {output} "
        "--study-nr-col=MARGINSstudyNr"

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

rule process_clincal_all_patients:
    input:
        script="src/data/process_clinical_all-patients.py",
        tsv="data/raw/imagene_clinical_all-patients.tsv",
    output:
        "data/processed/clinical_all-patients.nc"
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
        "{config[python]} {input.script} 8 {input.mri} {output}"


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
        protected("analyses/gsea/{mri}_{gene_set}_{abs,T|F}.Rds"),
    threads:
        4 # Takes a lot of memory
    shell:
        "mkdir -p analyses/gsea; "
        "{config[r]} {input.script} {input.gexp} {input.mri} "
        "{input.gene_sets} {output} --abs {wildcards.abs} --threads {threads} "
        "--perms 10000"

rule analyse_gene_sets_keep_es:
    input:
        script="src/analysis/analyse-gene-set-enrichment.R",
        gexp="data/processed/gene-expression.nc",
        mri="data/processed/{mri}.nc",
        gene_sets="data/external/msigdb/{gene_set}.v5.2.entrez.gmt",
    output:
        protected("analyses/gsea/{mri}_{gene_set}_{abs}+es.Rds"),
    threads:
        4 # Takes a lot of memory
    shell:
        "mkdir -p analyses/gsea; "
        "{config[r]} {input.script} {input.gexp} {input.mri} "
        "{input.gene_sets} {output} --abs {wildcards.abs} --threads {threads} "
        "--perms 10000 --return es_null"

rule gene_set_analysis_to_netcdf:
    input:
        script="src/analysis/gsea-rds-to-nc.R",
        rds="analyses/gsea{a}/{variables}_{gene_set_collection}_{abs}.Rds",
    output:
        "analyses/gsea{a,.*}/"
        "{variables}_{gene_set_collection,[^_/]+}_{abs,[TF]}.nc",
    shell:
        "{config[r]} {input.script}  {input.rds} {output} "
        "--gene-set-collection={wildcards.gene_set_collection} "
        "--abs={wildcards.abs}"

rule gene_set_analysis_to_xlsx:
    input:
        script="src/analysis/gsea-rds-to-xlsx.R",
        rds="analyses/gsea{a}/{name}.Rds",
    output:
        "analyses/gsea{a,.*}/{name}.xlsx",
    shell:
        "{config[r]} {input.script} {input.rds} {output} .25 0.0"


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
        "data/processed/mri-features-all.nc",
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

all_targets['figures'] = expand(
    "figures/{fig}.{ext}",
    fig=[
        "mri-cad-correlation",
        "fa-variance-explained",
        "cad-factors-heatmap",
        "gsea-heatmap_all-fa_c2.cgp_F_1",
        "gsea-heatmap_all-fa_c2.cp_T_3",
        "gsea-heatmap_all-fa_c2.cp_T_7",
        "clin-boxplot-ihc_subtype-volume",
    ],
    ext=['svg', 'pdf', 'png'],
) + ['figures/figure1.pdf', 'figures/figure1.png']

rule svg_to_pdf:
    input: "figures/{fn}.svg"
    output: "figures/{fn}.pdf"
    shell: "inkscape --export-pdf {output} -D {input}"

rule svg_to_png:
    input: "figures/{fn}.svg"
    output: "figures/{fn}.png"
    shell: "inkscape --export-png {output} -D -d 300 {input}"

rule figure_mri_cad_correlation:
    input:
        script="src/visualization/figure-mri-cad-correlation.py",
        cad_features="data/processed/mri-features-all.nc",
    output: "figures/mri-cad-correlation.svg"
    shell:
        "{config[python]} {input.script} {input.cad_features} {output}"

rule figure_fa_variance_explained:
    input:
        script="src/visualization/figure-fa-variance-explained.py",
        cad_features="data/processed/mri-features-all.nc",
    output: "figures/fa-variance-explained.svg"
    shell:
        "{config[python]} {input.script} {input.cad_features} {output}"

rule figure_cad_factors_heatmap:
    input:
        script="src/visualization/figure-cad-factors-heatmap.py",
        cad_factors="data/processed/mri-features-{subset}-fa.nc",
    output: "figures/cad-factors-{subset}-heatmap.svg"
    shell:
        "{config[python]} {input.script} {input.cad_factors} {output}"

rule figure_clin_boxplot_factor:
    input:
        script="src/visualization/figure-mri-factor-clin-boxplot.py",
        cad_factors="data/processed/mri-features-all-fa.nc",
        factor_annotation="config/factor_annot_all.yaml",
        clinical_annotation="data/processed/clinical.nc",
    output:
        "figures/clin-boxplotf-{clin}-{factor}.svg",
        "figures/clin-boxplotf-{clin}-{factor}_stats.txt",
    shell:
        "{config[python]} {input.script} "
        "{input.cad_factors} {wildcards.factor} {input.factor_annotation} "
        "{input.clinical_annotation} {wildcards.clin} "
        "{output}"

rule figure_clin_boxplot_feature:
    input:
        script="src/visualization/figure-mri-feature-clin-boxplot.py",
        mri_features="data/processed/mri-features-all.nc",
        clinical_annotation="data/processed/clinical.nc",
    output:
        "figures/clin-boxplot-{clin}-{feature}.svg",
        "figures/clin-boxplot-{clin}-{feature}_stats.txt",
    shell:
        "{config[python]} {input.script} "
        "{input.mri_features} {wildcards.feature} "
        "{input.clinical_annotation} {wildcards.clin} "
        "{output}"

rule figure_gsea_heatmap_fa:
    input:
        script="src/visualization/figure-gsea-heatmap.py",
        gsea="analyses/gsea/mri-features-{subset}-fa_{gene_set}_{abs}.nc",
        sel_genesets="src/visualization/"
            "sel-gs_{subset}_{gene_set}_{abs}_{factor}.tsv"
    output: "figures/gsea-heatmap_{subset}-fa_{gene_set}_{abs}_{factor}.svg"
    shell:
        "{config[python]} {input.script} {input.gsea} {input.sel_genesets} "
        "{wildcards.factor} {output}"



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
