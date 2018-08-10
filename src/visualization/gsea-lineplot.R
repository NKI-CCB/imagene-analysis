suppressPackageStartupMessages(library(tidyverse))
library(forcats)
library(ggplot2)

arg_spec = list(
    gsea_results=list(parse_fun=as.character),
    response_name=list(parse_fun=as.character),
    pdf_out=list(parse_fun=as.character))

parse_arguments <- function(arguments) {
    stopifnot(length(arguments) == length(arg_spec))
    arguments <- map(
        transpose(list(arg=arguments, spec=arg_spec)),
        function (a) {
            a$spec$parse_fun(a$arg)
        })
    names(arguments) <- names(arg_spec)
    arguments
}

gsea_lineplot <- function (gsea_results, response_name) {
    max_genesets <- 20

    gsea_res <- readRDS(gsea_results)
    gsea_df <- gsea_res$table[[response_name]]
    plot_df <- gsea_df %>%
        filter(fdr < .25) %>%
        arrange(-le_prop) %>%
        mutate(GeneSet=fct_rev(as_factor(GeneSet))) %>%
        slice(1:max_genesets)
    ggplot(plot_df, aes(x=max_es_at, y=GeneSet)) +
        geom_segment(aes(xend=0, yend=GeneSet, color=nes, size=le_prop)) +
        scale_color_gradient(guide=guide_colorbar(raster=T))

}

if (!exists('arguments')) {
    arguments <- parse_arguments(commandArgs(T))
    pdf(arguments$pdf_out, 14, 4)
    p <- gsea_lineplot(arguments$gsea_results, arguments$response_name)
    print(p)
    dev.off()
    invisible()
}
