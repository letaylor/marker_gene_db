# disable strings as factors, but re-enable upon exit
old <- options(stringsAsFactors = FALSE)
on.exit(options(old), add = TRUE)

suppressPackageStartupMessages(require(tidyverse))

#' Gene synonyms 2 HGNC symbols
#'
#' \code{gene_synonym2hgnc} converts gene synonyms to HGNC symbols HGNC symbols.
#'
#' @param gene_synonyms Character vector.
#'     List of gene names sum of which are hgnc symbols and some of which are
#'     gene synonyms.
#' @param select_first_synonyms_that_map_to_multiple_hgnc Logical.
#'     If TRUE, then sets the gene synonym to the first hgnc_symbol considered.
#'     If FALSE, sets gene synonyms that map to multiple hgnc_symbols to the
#'     original gene synonym.
#'
#' @return List
#'     List of hgnc symbols.
#'
#' @examples
#' gene_synonym2hgnc(c("CASP16P", "COL1A2", "KRT18", "COL6A2", "VWF", "HROB"))
#'
#' @importFrom org.Hs.eg.db org.Hs.eg_dbconn
#' @importFrom DBI dbGetQuery
#' @importFrom dplyr group_by
#' @importFrom dplyr filter
#' @importFrom dplyr n
#' @importFrom dplyr left_join
#' @importFrom magrittr %>%
#' @export
gene_synonym2hgnc <- function(
        gene_synonyms,
        select_first_synonyms_that_map_to_multiple_hgnc = FALSE
    ) {
    # use sql to get alias table and gene_info table (contains the symbols)
    # first open the database connection
    dbCon <- org.Hs.eg.db::org.Hs.eg_dbconn()
    # write your SQL query
    sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
    # execute the query on the database
    aliasSymbol <- data.frame(DBI::dbGetQuery(dbCon, sqlQuery))
    # delete extra columns
    for (col in c("gene_name", "X_id", "X_id.1")) {
        aliasSymbol[[col]] <- NULL
    }
    # get duplicate keys
    # duplicate_keys <- names(table(aliasSymbol$alias_symbol)[
    #   table(aliasSymbol$alias_symbol)>1
    # ])
    # print(subset(aliasSymbol, alias_symbol %in% duplicate_keys[[1]]))
    if (select_first_synonyms_that_map_to_multiple_hgnc) {
        # drop the >1st occurance  of duplicate keys
        aliasSymbol <- aliasSymbol[!duplicated(aliasSymbol$alias_symbol),]
    } else {
        # drop the all instances of duplicate keys
        aliasSymbol <- aliasSymbol %>%
            dplyr::group_by(alias_symbol) %>%
            dplyr::filter(dplyr::n() == 1)
    }
    # subset to get your results
    # result <- unique(merge(
    #     data.frame("alias_symbol" = gene_synonyms),
    #     aliasSymbol,
    #     by = c("alias_symbol"),
    #     all.x = TRUE,
    #     all.y = FALSE
    # ))
    # rownames(result) <- result$alias_symbol
    # result <- result[gene_synonyms,]$symbol
    # join does not change row order like marge
    result <- dplyr::left_join(
        data.frame("alias_symbol" = gene_synonyms),
        aliasSymbol,
        by = c("alias_symbol")
    )
    filt <- is.na(result$symbol)
    result$symbol[filt] <- result$alias_symbol[filt]
    result <- result$symbol
    if (length(result) != length(gene_synonyms)) {
        stop(paste0(
          "ERROR in gene_synonym2hgnc:\t",
          length(result),
          "\t",
          length(hgnc_symbols)
        ))
    }
    return(result)
}


#' Ensembl ids 2 hgnc
#'
#' \code{ensembl2hgnc} converts Ensembl gene ids to hgnc symbols. If no
#' hgnc symbol then uses external_gene_name. If no external_gene_name then
#' uses Ensembl gene id.
#'
#' @param ensembl_gene_ids Character vector.
#'     List of Ensembl gene ids to get hgnc symbols for.
#' @param host Character.
#'     Ensembl biomaRt host.
#' @param drop_dot_ensembl_id Logical.
#'     Drop "." from ENSG00000072310.12
#'
#' @return Character.
#'     The corresponding hgnc symbols to the Ensembl ids.
#'
#' @examples
#' ensembl2hgnc(c("ENSG00000109501", "ENSG00000058453", "ENSG00000030066"))
#'
#' @importFrom biomaRt useMart
#' @importFrom biomaRt getBM
#' @export
ensembl2hgnc <- function(
        ensembl_gene_ids,
        host = "grch37.ensembl.org",
        drop_dot_ensembl_id = TRUE
    ) {

    ensembl_gene_ids <- as.character(ensembl_gene_ids)
    ensembl_gene_ids_original <- ensembl_gene_ids
    if (drop_dot_ensembl_id) {
        ensembl_gene_ids <- unlist(lapply(ensembl_gene_ids,
            FUN = function(x) { return(strsplit(x, '\\.')[[1]][1]) }
        ))
    }
    # put this in a tryCatch block in case biomaRt is down!
    gene_info_tmp <- tryCatch({
        ensembl <- biomaRt::useMart(
            biomart = "ENSEMBL_MART_ENSEMBL",
            host = host,
            path = "/biomart/martservice",
            dataset = "hsapiens_gene_ensembl"
        )
        #biomaRt::listAttributes(ensembl)
        #cols <- c("ensembl_gene_id", "hgnc_symbol", "gene_biotype",
        #   "description")
        cols <- c("ensembl_gene_id", "hgnc_symbol", "external_gene_name")
        gene_info_tmp <- biomaRt::getBM(
            attributes = cols,
            mart = ensembl,
            filters = c("ensembl_gene_id"),
            values = unique(ensembl_gene_ids)
        )
    }, error = function(cond) {
        warning(paste0("Problems with biomaRt, probably due to a down",
                       " website. Returning ensembl_ids rather than hgnc_id."))
        return(data.frame())
    })

    # case where no results
    if (nrow(gene_info_tmp) == 0) {
        return_vector <- ensembl_gene_ids_original
        names(return_vector) <- ensembl_gene_ids_original
    } else {
        # make the gene_info_tmp the length of ensembl_gene_ids... so that
        # the return vector is the same length
        # this command also expands cases when there are duplicates.
        gene_info_tmp <- merge(
            data.frame(
                "ensembl_gene_id" = ensembl_gene_ids,
                "ensembl_gene_ids_original" = ensembl_gene_ids_original
            ),
            gene_info_tmp,
            by = c("ensembl_gene_id"),
            all.x = T,
            all.y = T
        )
        # if hgnc_symbol is empty, replace it with external_gene_name
        filt <- is.na(gene_info_tmp$hgnc_symbol) |
            gene_info_tmp$hgnc_symbol == ""
        if (any(filt)) {
            gene_info_tmp$hgnc_symbol[filt] <- as.character(
                gene_info_tmp$external_gene_name[filt]
            )
        }
        gene_info_tmp$external_gene_name <- NULL  # del external_gene_name
        # if hgnc_symbol is *still* empty, replace it with ensembl_gene_id
        filt <- is.na(gene_info_tmp$hgnc_symbol) |
            gene_info_tmp$hgnc_symbol == ""
        if (any(filt)) {
            gene_info_tmp$hgnc_symbol[filt] <- as.character(
                gene_info_tmp$ensembl_gene_id[filt]
            )
        }
        # get the return vector
        return_vector <- gene_info_tmp$hgnc_symbol
        names(return_vector) <- gene_info_tmp$ensembl_gene_ids_original
    }

    return(return_vector[ensembl_gene_ids_original])
}


#' Hgnc ids 2 ensembl ids
#'
#' \code{hgnc2ensembl} converts HGNC symbols to ensembl gene ids. If no
#' ensembl then empty vector.
#'
#' @param hgnc_symbol Character vector.
#'     List of hgnc symbols.
#' @param hgnc_symbol_mappings List.
#'     List ensembl mapping database to look up hgnc_symbol.
#' @param host Character.
#'     Ensembl biomaRt host. "ensembl.org" points to the latest version.
#' @param ensembl_version
#'     Ensembl version.
#' @param drop_dot_ensembl_id Logical.
#'     Drop "." from ENSG00000072310.12
#' @param filter_chromosome_patches Logical.
#'     Drop genes on chromosome patches.
#' @param verbose Logical.
#'     Print extra things.
#'
#' @return Data frame
#'     With two columns: hgnc_symbol, ensembl_gene_id
#'
#' @examples
#' hgnc2ensembl(c("INS", "PCSK1","PCSK2", "GLP1R"))
#'
#' @importFrom biomaRt useMart
#' @importFrom biomaRt getBM
#' @importFrom plyr rename
#' @export
hgnc2ensembl <- function(
    hgnc_symbols,
    hgnc_symbol_mappings = list(
        "hgnc_symbol", "external_gene_name", "external_synonym",
        "wikigene_name"
    ),
    host = "jan2020.archive.ensembl.org",
    ensembl_version = 99,
    drop_dot_ensembl_id = TRUE,
    filter_chromosome_patches = TRUE,
    verbose = TRUE
) {
    hgnc_symbols <- as.character(hgnc_symbols)
    if (verbose) {
        cat("Using hgnc_symbol_mappings of:\t", hgnc_symbol_mappings[[1]], "\n")
    }

    # If this function fails, it is likely because version does not match up
    # with available versions. To check, run biomaRt::listMarts(host = host)
    ensembl <- biomaRt::useMart(
        biomart = "ENSEMBL_MART_ENSEMBL",
        host = host,
        path = "/biomart/martservice",
        dataset = "hsapiens_gene_ensembl",
        version = paste("Ensembl Genes", ensembl_version)
    )
    #biomaRt::listAttributes(ensembl)
    query_results <- biomaRt::getLDS(
        attributes = c(hgnc_symbol_mappings[[1]]),
        filters = hgnc_symbol_mappings[[1]],
        values = hgnc_symbols,
        mart = ensembl,
        attributesL = c(
            #"external_gene_name",
            "chromosome_name",
            "ensembl_gene_id"
        ),
        martL = ensembl,
        uniqueRows = TRUE
    )
    colnames(query_results)[[1]] <- "HGNC.symbol"
    query_results <- plyr::rename(query_results, warn_missing = FALSE, c(
        "HGNC.symbol" = "hgnc_symbol",
        "Gene.name" = "external_gene_name",
        "Chromosome.scaffold.name" = "chromosome_name",
        "Gene.stable.ID" = "ensembl_gene_id"
    ))
    if (verbose) {
      message(
        "biomaRt query column example:\t",
        paste0(head(query_results[["hgnc_symbol"]]), collapse = ",")
      )
    }
    if (filter_chromosome_patches) {
        filt_chr_name <- grep(
            "CHR",
            query_results$chromosome_name,
            value = TRUE
        )
        if (length(filt_chr_name) > 0) {
            message(
              "Dropping chromosome patches:\n",
              paste(filt_chr_name, collapse = ",")
            )
            query_results <- subset(
                query_results,
                !(query_results$chromosome_name %in% filt_chr_name)
            )
        }
    }
    query_results$chromosome_name <- NULL
    results <- data.frame(
        "hgnc_symbol" = hgnc_symbols
    )
    results <- merge(
        results,
        query_results,
        by = c("hgnc_symbol"),
        all.x = TRUE,
        all.y = FALSE
    )

    if (drop_dot_ensembl_id) {
        results$ensembl_gene_id <- as.character(results$ensembl_gene_id)
        results$ensembl_gene_id <- unlist(lapply(results$ensembl_gene_id,
            FUN = function(x) { return(strsplit(x, '\\.')[[1]][1]) }
        ))
    }

    # add the ensembl version
    results$ensembl_version <- ensembl_version

    # apply recursive function to fill in NA
    na_filter <- is.na(results$ensembl_gene_id)
    if (any(na_filter) & (length(hgnc_symbol_mappings) > 1)) {
        if (verbose) {
            cat("recursive call (",
                paste0(hgnc_symbol_mappings, collapse = ","),
                "). number of ensembl_id NAs = ",
                sum(na_filter),
                "\n"
            )
        }
        hgnc_symbols_missing_ensembl <- results[na_filter,]$hgnc_symbol
        results_recursive_attempt_fill <- hgnc2ensembl(
            hgnc_symbols_missing_ensembl,
            hgnc_symbol_mappings = unlist(hgnc_symbol_mappings[-1]),
            host = host,
            ensembl_version = ensembl_version,
            drop_dot_ensembl_id = drop_dot_ensembl_id,
            filter_chromosome_patches = filter_chromosome_patches
        )
        results_no_na <- results[!na_filter,]
        results <- rbind(results_no_na, results_recursive_attempt_fill)
    }

    return(results)
}


main <- function() {
    df <- read.csv("data_unprocessed.tsv", sep = "\t")
    for (col in c("ensembl_gene_id", "ensembl_version")) {
        if (col %in% colnames(df)) {
            df[[col]] <- NULL
        }
    }
    df <- unique(df)
    res <- hgnc2ensembl(df$hgnc_symbol)
    results_temp <- unique(merge(
        df,
        res,
        by = c("hgnc_symbol"),
        all.x = TRUE,
        all.y = FALSE
    ))

    # try to fill in missing columns again by using another database to
    # identify hgnc
    use_gene_synonym2hgnc <- FALSE
    filt <- is.na(results_temp$ensembl_gene_id)
    if (use_gene_synonym2hgnc & (sum(filt) > 0)) {
        # try to map gene symbols and recall
        results_with_na <- results_temp[filt,]
        results_with_na$hgnc_symbol <- gene_synonym2hgnc(
            gene_synonyms = results_with_na$hgnc_symbol,
            select_first_synonyms_that_map_to_multiple_hgnc = FALSE
        )
        res2 <- hgnc2ensembl(unique(results_with_na$hgnc_symbol))
        results_temp2 <- unique(merge(
            results_with_na[colnames(df)],
            res2,
            by = c("hgnc_symbol"),
            all.x = TRUE,
            all.y = FALSE
        ))
        if (nrow(results_temp2) > nrow(results_with_na)) {
            warning(
                "Using gene_synonym2hgnc has resulted in genes with more",
                " than one ensembl id"
            )
        }
        results <- rbind(results_temp[!filt,], results_temp2)
    } else {
        results <- results_temp
    }

    results <- results[order(results$grouping_id, results$study_id),]
    tmp <- table(unique(results[c("ensembl_gene_id", "hgnc_symbol")]))
    if (any(colSums(tmp) > 1)) {
        message("WARNING: more than one ensembl ids for some genes:")
        tmp2 <- colSums(tmp)
        dup_genes <- names(tmp2[tmp2 > 1])
        print(subset(results, hgnc_symbol %in% dup_genes))
    }
    write.table(
        results,
        "data_processed.tsv",
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
