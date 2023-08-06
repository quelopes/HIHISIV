geo_query <- function(gse) {
  df_final <- NULL

  for (i in 1:length(gse)) {
    print(paste0("===== ", i, " ========"))
    query_result <- getGEO(gse[i], GSEMatrix = FALSE)

    df_final$geo_accession[i] <- Meta(query_result)$geo_accession
    df_final$title[i] <- ifelse(is.na(Meta(query_result)$title), NA, Meta(query_result)$title)
    df_final$pubmed_id[i] <- ifelse(is.null(Meta(query_result)$pubmed_id), NA, Meta(query_result)$pubmed_id)
    df_final$summary[i] <- ifelse(is.na(Meta(query_result)$summary), NA, Meta(query_result)$summary)
    df_final$overall_design[i] <- ifelse(is.na(Meta(query_result)$overall_design), "missing", Meta(query_result)$overall_design)

    df_final$keyword[i] <- ifelse(is.na(Meta(query_result)$summary), NA, Meta(query_result)$summary[2])
    df_final$file[i] <- ifelse(is.na(Meta(query_result)$supplementary_file), NA, Meta(query_result)$supplementary_file)
    df_final$type[i] <- ifelse(is.na(Meta(query_result)$type), NA, paste0(Meta(query_result)$type))
    df_final$taxid[i] <- ifelse(is.na(Meta(query_result)$taxid), NA, Meta(query_result)$taxid)
    df_final$sample_id[i] <- ifelse(is.na(Meta(query_result)$sample_id), NA, Meta(query_result)$sample_id)
    df_final$platform_id[i] <- ifelse(is.na(Meta(query_result)$platform_id), NA, Meta(query_result)$platform_id)

    rm(query_result)

    query_result_true <- getGEO(gse[i], GSEMatrix = TRUE)
    phenoData <- pData(phenoData(query_result_true[[1]]))


    pheno <- paste0(unique(phenoData$characteristics_ch1), unique(phenoData$source_name_ch1), collapse = " - ") # unique(phenoData$title)
    df_final$phenodata[i] <- pheno

    write.csv(phenoData, file = paste0("/results/phenoData_", df_final$geo_accession[i], ".csv"))


    rm(query_result_true)
    gc()
  }

  df_final
}