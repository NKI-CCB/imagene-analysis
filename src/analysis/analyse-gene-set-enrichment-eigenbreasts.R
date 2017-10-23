# CRAN
library(ncdf4)
library(optparse)
library(tibble)

# Bioconductor
loadNamespace('limma')

# Other
library(ggsea)

source('src/analysis/gsea-common.R')

read_mri <- function(fn, variable, k) {
    ds <- nc_open(fn)
    on.exit({nc_close(ds)})

    ds_mat <- nc_read_matrix(ds, variable)
    t(ds_mat[, 1:k])
}

parse_args <- function() {
    args <- c('gene_expression', 'mri', 'mri_variable', 'pc_number',
              'gene_sets', 'out')
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
    mri <- read_mri(args$mri, args$mri_variable, args$pc_number)
    str(mri)

    # Select data with no missing values
    sel_samples <- names(which(colSums(is.na(mri)) == 0)) 
    sel_samples <- intersect(sel_samples, colnames(gexp$read_count))
    mri <- mri[, sel_samples]
    gexp$read_count <- gexp$read_count[, sel_samples]

    res <- run_gsea(gexp$read_count, mri, gexp$entrez_gene, args$gene_sets,
                    nperm=args$perms, abs=args$abs, n_threads=args$threads,
                    gene_score_fn=score_genes_limma)
    saveRDS(res, args$out)

}

args <- parse_args()
main(args)
