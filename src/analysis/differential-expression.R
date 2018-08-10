library(ncdf4)
library(optparse)
library(tibble)

source('src/analysis/gsea-common.R')

ncstr_def <- function(name, dim, vals) {
    str_dim <- max(stringr::str_length(vals))
    str_dim <- ncdf4::ncdim_def(paste0('string', str_dim), '',
                                  seq_len(str_dim), F, F)
    ncdf4::ncvar_def(name, '', list(str_dim, dim), prec='char')
}

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

write_fit_to_nc <- function(fit, out, mri_features, genes, gene_symbols) {
    nc_vars <- list()
    gene_dim <-  ncdim_def('gene', '', seq_along(genes), F, F)
    nc_vars$gene <- ncstr_def('gene', gene_dim, genes)
    mri_dim <-  ncdim_def('mri_feature', '', seq_along(mri_features), F, F)
    nc_vars$mri <- ncstr_def('mri_feature', mri_dim, mri_features)
    nc_vars$coefficient=ncvar_def(
        'coefficient', '', list(gene_dim, mri_dim), prec='double')
    nc_vars$t=ncvar_def(
        't', '', list(gene_dim, mri_dim), prec='double')
    nc_vars$hgnc_symbol <- ncstr_def('hgnc_symbol', gene_dim, gene_symbols)

    out_f <- nc_create(out, nc_vars, force_v4=T)
    ncvar_put(out_f, 'gene', genes)
    ncvar_put(out_f, 'hgnc_symbol', gene_symbols)
    ncvar_put(out_f, 'mri_feature', mri_features)
    ncvar_put(out_f, 'coefficient',
              fit$coefficients[, 2:ncol(fit$coefficients)])
    ncvar_put(out_f, 't',
              fit$t[, 2:ncol(fit$t)])
}


parse_args <- function() {
    args <- c('gene_expression', 'mri', 'out')
    option_list <- list()
    usage <- paste("%prog [options] ",  paste(args, collapse=" "), collapse="")
    parser <- OptionParser(usage=usage, option_list=option_list)
    arguments <- optparse::parse_args(parser,
                                      positional_arguments=length(args))
    names(arguments$args) <- args
    c(as.list(arguments$args), arguments$options)
}

main <- function(args) {
    gexp <- read_gexp(args$gene_expression)
    mri <- read_mri(args$mri)

    # Select data with no missing values
    sel_samples <- names(which(colSums(is.na(mri)) == 0))
    sel_samples <- intersect(sel_samples, colnames(gexp$read_count))
    mri <- mri[, sel_samples]
    gexp_counts <- gexp$read_count[, sel_samples]
    gene_ids <- gexp$entrez_gene
    gene_symbols <- gexp$hgnc_symbol
    y = mri


    gene_sel <- rowSums(gexp_counts) > nrow(gexp_counts)
    gexp_counts <- gexp_counts[gene_sel, ]
    gene_ids <- gene_ids[gene_sel]
    gene_symbols <- gene_symbols[gene_sel]
    gexp_dge <- edgeR::DGEList(gexp_counts)
    gexp_dge <- edgeR::calcNormFactors(gexp_dge, method='TMM')
    design <- stats::model.matrix(~ ., as.data.frame(t(y)))
    gexp <- limma::voom(gexp_counts, design, plot=F)

    fit <- limma::lmFit(gexp, design)
    fit <- limma::eBayes(fit)

    write_fit_to_nc(fit, args$out, rownames(mri), rownames(gexp_counts),
                    gene_symbols=gene_symbols)
}

args <- parse_args()
main(args)
