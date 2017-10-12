# CRAN
library(ncdf4)
library(optparse)
library(tibble)

# Bioconductor
loadNamespace('limma')

# Other
library(ggsea)

source('src/analysis/gsea-common.R')

read_factors <- function(fn) {
    ds <- nc_open(fn)
    nc_read_matrix(ds, 'factor_value')
}

parse_args <- function() {
    args <- c('gene_expression', 'sfa_solution', 'gene_sets', 'out')
    option_list <- list(
        make_option("--threads", default=1L, help="Number of threads"),
        make_option("--perms", default=1000L, help="Number of permutations"),
        make_option("--abs", default=FALSE, help="Use absolute scores"))
    usage <- paste("%prog [options] ",  paste(args, collapse=" "), collapse="")
    parser <- OptionParser(usage=usage, option_list=option_list)
    arguments <- optparse::parse_args(parser,
                                      positional_arguments=length(args))
    names(arguments$args) <- args
    c(as.list(arguments$args), arguments$options)
}

main <- function(args) {

    if (args$threads > 1) {
        library(doMC)
        registerDoMC(args$threads)
    }

    gexp <- read_gexp(args$gene_expression)
    factors <- read_factors(args$sfa_solution)

    # Select data with no missing values
    sel_samples <- colnames(factors)
    sel_samples <- intersect(sel_samples, colnames(gexp$read_count))
    factors <- factors[, sel_samples]
    gexp$read_count <- gexp$read_count[, sel_samples]

    res <- run_gsea(gexp$read_count, factors, gexp$entrez_gene, args$gene_sets,
                    nperm=args$perms, abs=args$abs, n_threads=args$threads)
    saveRDS(res, args$out)

}

args <- parse_args()
main(args)
