library('optparse')
library('dplyr')
library('openxlsx')

msigdb_prefix = paste0('http://software.broadinstitute.org/gsea/msigdb/',
                       'geneset_page.jsp?geneSetName=')
num_style = openxlsx::createStyle(numFmt='NUMBER')

parse_args <- function() {
    args <- c('gsea_rds', 'out', 'fdr', 'le_prop')
    option_list <- list()
    usage <- paste("%prog [options] ",  paste(args, collapse=" "), collapse="")
    parser <- OptionParser(usage=usage, option_list=option_list)
    arguments <- optparse::parse_args(parser,
                                      positional_arguments=length(args))
    names(arguments$args) <- args
    res <- c(as.list(arguments$args), arguments$options)
    res$fdr <- as.numeric(res$fdr)
    res$le_prop <- as.numeric(res$le_prop)
    res
}

gsea_rds_to_xlsx <- function(gsea_rds, out, max_fdr, min_le_prop) {
    gsea_table <- readr::read_rds(gsea_rds)$table
    workbook <- openxlsx::createWorkbook('GGSEA')
    for (factor_name in names(gsea_table)) {
        worksheet_index <- openxlsx::addWorksheet(workbook, factor_name)
        filtered_table <- gsea_table[[factor_name]] %>%
            filter(fdr <= max_fdr, le_prop >= min_le_prop) %>%
            arrange(-le_prop) %>%
            select(-fwer, -es)
        if (nrow(filtered_table) > 0) {
            filtered_table$link = paste0(
                '=HYPERLINK("', msigdb_prefix, filtered_table$GeneSet,
                '", "link")')
            class(filtered_table$link) <- 'formula'
        } else {
            filtered_table$link = character()
        }
        openxlsx::writeDataTable(workbook, worksheet_index, filtered_table)
        openxlsx::addStyle(workbook, worksheet_index, num_style,
                           cols=which(colnames(filtered_table) %in%
                                      c('es', 'p', 'nes', 'fdr', 'le_prop',
                                        'fwer')),
                           rows=2:(nrow(filtered_table)+1), gridExpand=T)
        openxlsx::setColWidths(workbook, worksheet_index,
                               1:ncol(filtered_table), 'auto')
    }
    openxlsx::saveWorkbook(workbook, args$out, T)
}

args <- parse_args()
gsea_rds_to_xlsx(args$gsea_rds, args$out, args$fdr, args$le_prop)
