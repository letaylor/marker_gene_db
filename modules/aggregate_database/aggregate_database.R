# disable strings as factors, but re-enable upon exit
old <- options(stringsAsFactors = FALSE)
on.exit(options(old), add = TRUE)

suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(stringr))

main <- function() {
    # Creates an output file with the following columns:
    final_columns <- c(
        "score_id",
        "hgnc_symbol",
        "grouping_id",
        "study_id",
        "description",
        "ensembl_gene_id",
        "ensembl_version"
    )

    # List all of the directories within the database.
    database_files <- list.files(
        path = "../..", 
        pattern = "*database.tsv", 
        all.files = FALSE,
        full.names = TRUE, 
        recursive = TRUE
    )
    
    df_list <- list()

    # Process data in the celltypes dir. To make score_id, we add the dir
    # within celltypes to grouping_id.
    for(f in grep("celltypes", database_files, value = TRUE)) {
        df_list[[f]] <- read.csv(f, sep = "\t")
        subdir_celltypes <- strsplit(
            strsplit(f, "celltypes/")[[1]][2], 
            "/"
        )[[1]][1]
        df_list[[f]]$score_id <- paste0(
            subdir_celltypes, 
            "__", 
            unlist(lapply(
                df_list[[f]]$grouping_id, 
                FUN = function(x) {
                    #str_replace(tolower(x), subdir_celltypes, "")
                    tolower(x)
                }
            ))
        )
        # Clean up cases like b_cell__b_cell_activation
        replace_string <- paste0(
            c(subdir_celltypes, subdir_celltypes), 
            collapse = "__"
        )
        df_list[[f]]$score_id <- unlist(lapply(
            df_list[[f]]$score_id, 
            FUN = function(x) {
                stringr::str_replace(
                    x,
                    replace_string,
                    paste0(subdir_celltypes, "_")
                )
            }
        ))
    }
    
    # Process data in signatures dir. To make score_id, we use grouping_id.
    for(f in grep("signatures", database_files, value = TRUE)) {
        df_list[[f]] <- read.csv(f, sep = "\t")
        subdir_celltypes <- strsplit(
            strsplit(f, "signatures/")[[1]][2], 
            "/"
        )[[1]][1]
        if (subdir_celltypes == "cell_cycle") {
            df_list[[f]]$score_id <- "cell_cycle"
        } else {
            df_list[[f]]$score_id <- unlist(lapply(
                df_list[[f]]$grouping_id, 
                FUN = function(x) {
                  tolower(x)
                }
            ))
        }
    }
    
    # Merge the data.
    df <- data.frame(data.table::rbindlist(df_list, fill = TRUE))
    print(table(df$score_id))
    
    # Write the final data, subsetting to the columns.
    write.table(
        df[final_columns],
        "data_aggregated.tsv",
        row.names = FALSE,
        col.names = TRUE,
        quote = TRUE,
        sep = "\t",
        na = ""
    )
}
main()

# like python if __name__ == '__main__'
if (sys.nframe() == 0) {
    main()
}
