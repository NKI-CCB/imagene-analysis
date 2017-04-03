# Require R 3.3, don't forget to update Bioconductor as well when updating R.
stopifnot(version['major'] == '3')
stopifnot(substr(version['minor'], 1, 1) == '3')

cran_requirements = c(
    'tidyverse',
    'devtools',
    'ggplot2')

bioc_requirements = c(
    )

github_requirements = c(
    ggsea='NKI-CCB/ggsea'
)


# Installation script #

to_install_cran = setdiff(cran_requirements, installed.packages()[,1])
if (length(to_install_cran) > 0) {
    r <- getOption("repos")
    r["CRAN"] <- "https://cloud.r-project.org/"
    options(repos=r)
    install.packages(to_install_cran)
}

to_install_bioc = setdiff(bioc_requirements, installed.packages()[,1])
if (length(to_install_bioc) > 0) {
    if (!(requireNamespace('BiocInstaller'))) {
        install.packages('BiocInstaller',
                         repo='https://bioconductor.org/packages/3.4/bioc/')
        loadNamespace('BiocInstaller')
    }
    BiocInstaller::biocLite(to_install_bioc)
}

to_install_github = setdiff(names(github_requirements),
                            installed.packages()[,1])
if (length(to_install_github) > 0) {
    loadNamespace('devtools')
    devtools::install_github(github_requirements[to_install_github])
}
