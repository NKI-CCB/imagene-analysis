mri-rnaseq-integration
==============================

Integration of MRI and RNAseq data of the Imagene project.

To make a virtual environment for R and Python usage:
```sh
make requirements
```

Then you can process the data, train models and run analyses:
```sh
make all
```

Project Organization
------------

    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make all` and `make venv`.
    ├── project-rules.Makefile <- Makefile with rules to run this project. 
    |                                These rules are run inside a virtual environment.
    ├── README.md          <- The top-level README for developers using this project.
    │
    ├── docs               <- Project Documentation using Sphinx.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── data
    │   ├── raw            <- The original, immutable data dump.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── external       <- Data from third party sources, such as reference data.
    │
    ├── notebooks          <- Jupyter notebooks for exploration. Naming convention is date (for
    │                         ordering), the creator's initials, and a short `-` delimited
    │                         description, e.g. `2017-01-01-jan-initial-data-exploration`.
    │
    ├── models             <- Trained models and model predictions.
    │
    ├── analyses           <- Results and summaries of statistical analyses.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting.
    │
    ├── venv/              <- Suggested directory for virtual environment
    │
    ├── requirements.txt   <- The Python requirements file for reproducing the analysis 
    │                         environment e.g. generated with `pip freeze > requirements.txt`.
    ├── requirements.R     <- Requirements and installation script for R.
    │
    └── src                <- Source code for use in this project.
        │
        ├── lib            <- Modules with function useful for multiple scripts in this project.
        │
        ├── data           <- Scripts to download or generate data.
        │
        ├── models         <- Scripts to train and apply models.
        │
        ├── analysis       <- Scripts for statistical analysis of data or models.
        │
        ├── reports        <- Source of reports, such as pweave or rmarkdown reports.
        │
        └── visualization  <- Scripts to visualize results.



--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
