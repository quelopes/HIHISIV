extract_gse <- function(way_text) {
  
  texto <- readLines(way_text)


  # Create an empty dataframe with the desired column names
  df <- data.frame(
    title = character(),
    supplied = character(),
    organism = character(),
    type = character(),
    platform = character(),
    n_samples = character(),
    ftp = character(),
    series_accession = character(),
    stringsAsFactors = FALSE
  )

  # Iterate over the lines and extract the relevant information
  for (i in seq_along(texto)) {
    line <- texto[i]

    if (grepl("^\\d+\\. ", line)) {
      title <- gsub("^\\d+\\. ", "", line)

      j <- i + 1
      supplied <- ""
      while (j <= length(texto) && texto[j] != "") {
        supplied <- paste(supplied, texto[j])
        j <- j + 1
      }
      organism <- NA
      type <- NA
      platform <- NA
      n_samples <- NA
      ftp <- NA
      series_accession <- NA

      for (k in i:j) {
        if (grepl("^Organism:", texto[k])) {
          organism <- gsub("^Organism:\\s+", "", texto[k])
        } else if (grepl("^Type:", texto[k])) {
          type <- gsub("^Type:\\s+", "", texto[k])
        } else if (grepl("^Platform:", texto[k])) {
          platform <- gsub("^Platform:\\s+GPL(\\d+).*", "GPL\\1", texto[k])
          n_samples <- gsub("^Platform:\\s+.* (\\d+) Samples.*", "\\1", texto[k])
        } else if (grepl("^Series\\s+Accession:", texto[k])) {
          series_accession <- gsub("^Series\\s+Accession:\\s+(GSE\\d+).*", "\\1", texto[k])
        } else if (grepl("^FTP download:", texto[k])) {
          ftp <- gsub("^FTP download:\\s+GEO \\(TXT\\) ftp://(.*?)/.*", "ftp://\\1", texto[k])
          ftp <- gsub("ftp://.*?(geo/series/GSE\\d+[^/]*)/.*", "ftp://\\1", ftp)
        }
      }

      # Temporary dataframe with the extracted information
      temp_df <- data.frame(
        title = title,
        supplied = ifelse(nchar(trimws(supplied)) > 0, supplied, NA),
        organism = ifelse(nchar(trimws(organism)) > 0, organism, NA),
        type = ifelse(nchar(trimws(type)) > 0, type, NA),
        platform = ifelse(nchar(trimws(platform)) > 0, platform, NA),
        n_samples = ifelse(nchar(trimws(n_samples)) > 0, n_samples, NA),
        ftp = ifelse(nchar(trimws(ftp)) > 0, ftp, NA),
        series_accession = ifelse(nchar(trimws(series_accession)) > 0, series_accession, NA),
        stringsAsFactors = FALSE
      )

      df <- rbind(df, temp_df)
    }
  }

  df

  
}