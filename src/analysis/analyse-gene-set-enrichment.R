# CRAN
library(ncdf4)
library(optparse)
library(tibble)

# Bioconductor
loadNamespace('limma')

# Other
library(ggsea)

# Read a matrix from a NetCDF4 file
nc_read_matrix <- function(ds, var_name) {
    mat = ncvar_get(ds, var_name)

    row_dim <- ds[['var']][[var_name]][['dim']][[1]]
    col_dim <- ds[['var']][[var_name]][['dim']][[2]]
    mat_dim <- list(c(row_dim[['vals']]), c(col_dim[['vals']]))
    names(mat_dim) <- c(row_dim[['name']], col_dim[['name']])
    dimnames(mat) <- mat_dim
    mat
}

read_gexp <- function(fn) {
    ds <- nc_open(fn)
    on.exit({nc_close(ds)})

    read_count <- nc_read_matrix(ds, 'read_count')
    entrez_gene <- ncvar_get(ds, 'entrez_gene_id')
    entrez_gene <- as.character(entrez_gene)

    list(
         read_count=read_count,
         entrez_gene=entrez_gene)
}

# Read a data frame from a NetCDF4 file
nc_read_data_frame <- function(ds, dim) {
    v <- ds$dim[[dim]]$vals
    df <- tibble(v)
    names(df) <- c(dim)

    for (var_name in names(ds[['var']])) {
        var <- ds[['var']][[var_name]]
        if ((length(var[['dim']]) == 1) &&
                (var[['dim']][[1]][['name']] == dim)) {
            df[[var_name]] <- ncvar_get(ds, var_name)
        }
    }
    df
}

read_mri <- function(fn) {
    ds <- nc_open(fn)
    on.exit({nc_close(ds)})

    mri <- nc_read_data_frame(ds, 'case')
    # Integer variables as index are not funny in R, so we make it character.
    mri$case <- as.character(mri$case)
    mri$MultiFocal <- NULL
    mri_mat <- t(as.matrix(mri[, sapply(mri, is.numeric)]))
    dimnames(mri_mat) <- list(
        mri_feature=rownames(mri_mat),
        case=mri$case)
    mri_mat
}

# Gene scoring with limma for use in the gene set enrichment analysis
score_genes_limma <- function (x, y, abs=F) {
    design <- stats::model.matrix(~ ., as.data.frame(y))
    fit <- limma::lmFit(x, design)
    fit <- limma::eBayes(fit)
    res <- fit$t[,2:ncol(fit$t), drop=F]
    if (abs) {
        abs(res)
    } else {
        res
    }
}

run_gsea <- function(gexp_counts, mri, gene_ids, gs_fn, nperm, abs,
                     n_threads) {

    stopifnot(colnames(gexp_counts) == colnames(mri))

    gene_sel <- rowSums(gexp_counts) > nrow(gexp_counts)
    gexp_counts <- gexp_counts[gene_sel, ]
    gene_ids <- gene_ids[gene_sel]
    gexp_dge <- edgeR::DGEList(gexp_counts)
    gexp_dge <- edgeR::calcNormFactors(gexp_dge, method='TMM')
    design <- stats::model.matrix(~ ., as.data.frame(t(mri)))
    gexp <- limma::voom(gexp_counts, design, plot=F)

    ggsea(gexp, t(mri), gs_fn,
        gene.score.fn=score_genes_limma,
        gene.names=gene_ids,
        es.fn=ggsea_weighted_ks, sig.fun=ggsea_calc_sig,
        verbose=T, nperm=nperm, block.size=64, abs=abs,
        parallel=(n_threads > 1))
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

main <- function() {
    args <- parse_args()

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

main()
