Integration of MRI and RNAseq data in the Imagene Project
=========================================================

Integration of MRI and RNAseq data of the Imagene project. Paper published in Radiology: "Radiogenomic analysis of breast cancer by linking MRI phenotypes with tumor gene expression." https://doi.org/10.1148/radiol.2020191453.

To make a virtual environment for R and Python usage:
```sh
make requirements
```

Then you can process the data, train models and run analyses:
```sh
make all
```

Remote files
------------

Data is not included in this repository. The download method and location of
the data can be specified in the `config/snakemake.yaml` file.

### Pathway Analysis Requirements ###

The pathway analyses requires the msigdb files release 5.2. These cannot be
downloaded automatically, as registration is required. Before running the
pathway analysis, download `msigdb_v5.2_files_to_download_locally.zip`
from the [msigdb site](http://software.broadinstitute.org/gsea/downloads.jsp#msigdb)
under Archived Releases and place it into `data/external/msigdb`. The scripts
will extract the required files from the archive.



Project Organization
------------

    ├── LICENSE
    │
    ├── Makefile           <- Makefile with commands like `make all` and `make venv`.
    │
    ├── Snakefile          <- Snakemake file with rules to run this project. These rules
    │                         should be run inside a virtual environment containing the packages
    │                         in requirements.txt (Python) and requirement.R (R). The Makefile can
    │                         set this environment up on some systems.
    │
    ├── README.md          <- The top-level README for developers using this project.
    │
    ├── docs               <- Project Documentation using Sphinx.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── config/snakemake.yaml <- Configuration of the workflow. Set version of Python and R used here.
    │
    ├── data
    │   ├── raw            <- The original, immutable data dump.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   ├── external       <- Data from third party sources, such as reference data.
    │   └── to_share       <- Data for third parties.
    │
    ├── notebooks          <- Jupyter notebooks for exploration. Naming convention is date (for
    │                         ordering), the creator's initials, and a short `-` delimited
    │                         description, e.g. `2017-01-01-tb-initial-data-exploration`.
    │
    ├── models             <- Trained models and model predictions.
    │
    ├── analyses           <- Results and summaries of statistical analyses.
    │
    ├── reports            <- Analysis reports, such as pweave or rmarkdown reports.
    │   └── figures        <- Generated graphics and figures to be used in reporting.
    ├── figures            <- Figures for in the manuscript.
    │
    ├── venv/              <- Suggested directory for virtual environment
    │
    ├── requirements.txt   <- The Python requirements file for reproducing the analysis
    │                         environment e.g. generated with `pip freeze > requirements.txt`.
    ├── requirements.R     <- Requirements and installation script for R.
    │
    └── src                <- Source code for use in this project.
        │
        ├── plot.py        <- Plotting library.
        ├── util.py        <- Utility function library.
        │
        ├── analysis       <- Scripts for statistical analysis of data or models.
        │
        ├── data           <- Scripts to download or generate data.
        │
        ├── features       <- Scripts to construct features from data.
        │
        ├── models         <- Scripts to train and apply models.
        │
        ├── reports        <- Source of reports, such as snippets used in multiple reports.
        │
        └── visualization  <- Scripts to visualize results.



--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
