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

# Gene scoring with limma for use in the gene set enrichment analysis
score_genes_limma_independend <- function (x, y, abs=F) {
    res = matrix(NA, nrow(x), ncol(y))
    for (i in seq_len(ncol(y))) {
        design <- stats::model.matrix(~ y[, i])
        fit <- limma::lmFit(x, design)
        fit <- limma::eBayes(fit)
        res[, i] <- fit$t[, 2]
    }
    if (abs) {
        abs(res)
    } else {
        res
    }
}

run_gsea <- function(gexp_counts, y, gene_ids, gs_fn, nperm, abs,
                     n_threads, gene_score_fn=score_genes_limma) {

    stopifnot(colnames(gexp_counts) == colnames(y))

    gene_sel <- rowSums(gexp_counts) > nrow(gexp_counts)
    gexp_counts <- gexp_counts[gene_sel, ]
    gene_ids <- gene_ids[gene_sel]
    gexp_dge <- edgeR::DGEList(gexp_counts)
    gexp_dge <- edgeR::calcNormFactors(gexp_dge, method='TMM')
    design <- stats::model.matrix(~ ., as.data.frame(t(y)))
    gexp <- limma::voom(gexp_counts, design, plot=F)

    ggsea(gexp, t(y), gs_fn,
        gene.score.fn=gene_score_fn,
        gene.names=gene_ids,
        es.fn=ggsea_weighted_ks, sig.fun=ggsea_calc_sig,
        verbose=T, nperm=nperm, block.size=64, abs=abs,
        parallel=(n_threads > 1))
}
