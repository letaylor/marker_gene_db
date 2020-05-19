# disable strings as factors, but re-enable upon exit
old <- options(stringsAsFactors = FALSE)
on.exit(options(old), add = TRUE)

require(dplyr)

#' mgi ids 2 human ensembl ids
#'
#' \code{hgnc2ensembl} converts HGNC symbols to ensembl gene ids. If no
#' ensembl then empty vector.
#'
#' @param hgnc_symbol Character vector.
#'     List of hgnc symbols.
#' @param host Character.
#'     Ensembl biomaRt host. "ensembl.org" points to the latest version.
#' @param ensembl_version
#'     Ensembl version.
#' @param drop_dot_ensembl_id Logical.
#'     Drop "." from ENSG00000072310.12
#' @param filter_chromosome_patches Logical.
#'     Drop genes on chromosome patches.
#'
#' @return Data frame
#'     With two columns: hgnc_symbol, ensembl_gene_id
#'
#' @examples
#' mgi2hgnc_ensembl(c("Tlr1", "Tmsb4x", "Ythdc2", "Zfp280d"))
#'
#' @importFrom biomaRt useMart
#' @importFrom biomaRt getLDS
#' @importFrom plyr rename
#' @importFrom dplyr left_join
#' @export
mgi2hgnc_ensembl <- function(
    mgi_symbols,
    host = "jan2020.archive.ensembl.org",
    ensembl_version = 99,
    drop_dot_ensembl_id = TRUE,
    filter_chromosome_patches = TRUE
) {
    # If this function fails, it is likely because version does not match up
    # with available versions. To check, run biomaRt::listMarts(host = host)
    human <- biomaRt::useMart(
        biomart = "ENSEMBL_MART_ENSEMBL",
        host = host,
        path = "/biomart/martservice",
        dataset = "hsapiens_gene_ensembl",
        version = paste("Ensembl Genes", ensembl_version)
    )
    mouse <- biomaRt::useMart(
        biomart = "ENSEMBL_MART_ENSEMBL",
        host = host,
        path = "/biomart/martservice",
        dataset = "mmusculus_gene_ensembl",
        version = paste("Ensembl Genes", ensembl_version)
    )
    
    query_results <- biomaRt::getLDS(
        attributes = c("mgi_symbol", "ensembl_gene_id"), 
        filters = "mgi_symbol", 
        values = mgi_symbols, 
        mart = mouse, 
        attributesL = c("hgnc_symbol", "chromosome_name", "ensembl_gene_id"), 
        martL = human, 
        uniqueRows = TRUE
    )
    query_results <- plyr::rename(query_results, warn_missing = FALSE, c(
        "MGI.symbol" = "mgi_symbol",
        "Gene.stable.ID" = "ensembl_mouse_gene_id",
        "HGNC.symbol" = "hgnc_symbol",
        "Gene.name" = "external_gene_name",
        "Chromosome.scaffold.name" = "chromosome_name",
        "Gene.stable.ID.1" = "ensembl_gene_id"
    ))
    if (filter_chromosome_patches) {
        filt <- grep("CHR", query_results$chromosome_name, value = TRUE)
        message("Dropping chromosome patches:\n", paste(filt, collapse = ","))
        query_results <- subset(
            query_results,
            !(query_results$chromosome_name %in% filt)
        )
    }
    query_results$chromosome_name <- NULL
    results <- dplyr::left_join(
        data.frame("mgi_symbol" = mgi_symbols), 
        query_results,
        by = c("mgi_symbol")
    )
    
    if (drop_dot_ensembl_id) {
        results$ensembl_mouse_gene_id <- unlist(lapply(
            results$ensembl_mouse_gene_id,
            FUN = function(x) { return(strsplit(x, '\\.')[[1]][1]) }
        ))
        results$ensembl_gene_id <- unlist(lapply(
            results$ensembl_gene_id,
            FUN = function(x) { return(strsplit(x, '\\.')[[1]][1]) }
        ))
    }
    
    # add the ensembl version
    results$ensembl_version <- ensembl_version

    return(results)
}


main <- function() {
    df <- read.csv("data_unprocessed-mouse.tsv", sep = "\t")
    if ("ensembl_gene_id" %in% colnames(df)) {
        df$ensembl_gene_id <- NULL
    }
    df <- unique(df)
    res <- mgi2hgnc_ensembl(
        df$mgi_symbol
    )
    results <- unique(merge(
        df,
        res,
        by = c("mgi_symbol"),
        all.x = TRUE,
        all.y = FALSE
    ))
    results <- results[order(results$grouping_id, results$study_id),]
    tmp <- table(unique(results[c("ensembl_gene_id", "mgi_symbol")]))
    if (any(colSums(tmp) > 1)) {
        message("WARNING: more than one ensembl ids for some genes:")
        tmp2 <- colSums(tmp)
        dup_genes <- names(tmp2[tmp2 > 1])
        print(subset(results, mgi_symbol %in% dup_genes))
    }
    results <- results %>%
        dplyr::select(
            hgnc_symbol, 
            grouping_id,
            study_id,
            description,
            ensembl_gene_id,
            ensembl_version,
            dplyr::everything()
        )
    write.table(
        results,
        "data_processed-mouse.tsv",
        row.names = FALSE,
        col.names = TRUE,
        quote = TRUE,
        sep = "\t",
        na = ""
    )
}


# like python if __name__ == '__main__'
if (sys.nframe() == 0) {
    main()
}
