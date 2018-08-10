library('ncdf4')
library('optparse')

ncstr_def <- function(name, dim, vals) {
    str_dim <- max(stringr::str_length(vals))
    str_dim <- ncdf4::ncdim_def(paste0('string', str_dim), '',
                                  seq_len(str_dim), F, F)
    ncdf4::ncvar_def(name, '', list(str_dim, dim), prec='char')
}

long_names <- c(
    es='enrichment statistic',
    p='nominal p-value of enrichment statistic',
    nes='normalized enrichment statistic',
    fdr='false discovery rate of normalized enrichment statistic',
    fwer='Bonferroni adjusted p-value of enrichment statistic',
    max_es_at='gene rank where enrichment statistic is reached',
    le_prop='proportion of genes in gene set in leading edge',
    GeneSet='gene set identifier')

parse_args <- function() {
    args <- c('gsea_rds', 'out')
    option_list <- list(
        make_option('--gene-set-collection', default=""),
        make_option('--abs', default=""))
    usage <- paste("%prog [options] ",  paste(args, collapse=" "), collapse="")
    parser <- OptionParser(usage=usage, option_list=option_list)
    arguments <- optparse::parse_args(parser,
                                      positional_arguments=length(args))
    names(arguments$args) <- args
    c(as.list(arguments$args), arguments$options)
}

main <- function() {
    args <- parse_args()

    res <- readRDS(args$gsea_rds)
    vars <- setdiff(colnames(res$table[[1]]), 'GeneSet')
    gene_sets <- res$table[[1]]$GeneSet
    mri_features <- names(res$table)

    gene_set_dim <-  ncdim_def('gene_set', '', seq_along(gene_sets), F, F)
    gene_set_var <- ncstr_def('gene_set', gene_set_dim, gene_sets)
    mri_dim <-  ncdim_def('mri_feature', '', seq_along(res$table), F, F)
    mri_var <- ncstr_def('mri_feature', mri_dim, mri_features)
    nc_vars <- lapply(vars, function (var_name) {
        v <- res$table[[1]][[var_name]]
        if (is.integer(v)) {
            prec <- 'integer'
        } else if (is.double(v)) {
            prec <- 'double'
        } else {
            stop("Could not determine data type")
        }
        ncvar_def(
            var_name, '', list(gene_set_dim, mri_dim),
            prec=prec, longname=long_names[var_name])
    })
    nc_vars <- append(nc_vars, list(mri_var, gene_set_var))

    out_f <- nc_create(args$out, nc_vars, force_v4=T)
    ncvar_put(out_f, 'gene_set', gene_sets)
    ncvar_put(out_f, 'mri_feature', mri_features)
    for (feature_i in seq_along(mri_features)) {
        feature_name <- mri_features[[feature_i]]
        stopifnot(res$table[[feature_name]]$GeneSet == gene_sets)
        for (var_name in vars) {
            ncvar_put(out_f, var_name, res$table[[feature_name]][[var_name]],
                      start=c(1, feature_i), count=c(length(gene_sets), 1))
        }
    }
    ncatt_put(out_f, 0, 'gene_set_collection', args$`gene-set-collection`,
              'text')
    ncatt_put(out_f, 0, 'absolute', as.logical(args$abs), 'short')
    nc_close(out_f)
}

main()
