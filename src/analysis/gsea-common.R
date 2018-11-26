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
    hgnc_symbol <- ncvar_get(ds, 'hgnc_symbol')

    list(
         read_count=read_count,
         hgnc_symbol=hgnc_symbol,
         entrez_gene=entrez_gene)
}

run_gsea <- function(gexp_counts, y, gene_ids, gs_fn, nperm, abs,
                     n_threads, gene_score_fn=score_genes_limma,
                     return_values) {

    stopifnot(colnames(gexp_counts) == colnames(y))

    gene_sel <- rowSums(gexp_counts) > nrow(gexp_counts)
    gexp_counts <- gexp_counts[gene_sel, ]
    gene_ids <- gene_ids[gene_sel]
    gexp_dge <- edgeR::DGEList(gexp_counts)
    gexp_dge <- edgeR::calcNormFactors(gexp_dge, method='TMM')
    design <- stats::model.matrix(~ ., as.data.frame(t(y)))
    gexp <- limma::voom(gexp_dge, design, plot=F)

    ggsea(gexp, design, gs_fn,
        gene.score.fn=ggsea_limma,
        gene.names=gene_ids,
        es.fn=ggsea_weighted_ks, sig.fun=ggsea_calc_sig,
        verbose=T, nperm=nperm, block.size=64, abs=abs,
        parallel=(n_threads > 1),
        return_values=return_values)
}
