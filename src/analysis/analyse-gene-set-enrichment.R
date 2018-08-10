# CRAN
library(ncdf4)
library(optparse)
library(tibble)

# Bioconductor
loadNamespace('limma')

# Other
library(ggsea)

source('src/analysis/gsea-common.R')

read_mri <- function(fn) {
    ds <- nc_open(fn)
    on.exit({nc_close(ds)})

    if ('factors' %in% names(ds$var)) {
        mri_mat = nc_read_matrix(ds, 'factors')
    } else {
        mri <- nc_read_data_frame(ds, 'case')
        # Integer variables as index are not funny in R, so we make it character.
        mri$case <- as.character(mri$case)
        mri$MultiFocal <- NULL
        mri_mat <- t(as.matrix(mri[, sapply(mri, is.numeric)]))
        dimnames(mri_mat) <- list(
            mri_feature=rownames(mri_mat),
            case=mri$case)
    }
    mri_mat
}

parse_args <- function() {
    args <- c('gene_expression', 'mri', 'gene_sets', 'out')
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
    mri <- read_mri(args$mri)

    # Select data with no missing values
    sel_samples <- names(which(colSums(is.na(mri)) == 0)) 
    sel_samples <- intersect(sel_samples, colnames(gexp$read_count))
    mri <- mri[, sel_samples]
    gexp$read_count <- gexp$read_count[, sel_samples]

    res <- run_gsea(gexp$read_count, mri, gexp$entrez_gene, args$gene_sets,
                    nperm=args$perms, abs=args$abs, n_threads=args$threads)
    saveRDS(res, args$out)

}

args <- parse_args()
main(args)
