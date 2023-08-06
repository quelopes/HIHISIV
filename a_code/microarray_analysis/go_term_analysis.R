go_terms_analysis <- function(df_annot, org_annot, name_exp, way_out) {
  print("*** ACTIVITY GO ANALYSIS ***")

  # --- parameters ---
  hgCutoff <- 0.001

  # remove probes without entrez id
  df_annot <- df_annot[-which(is.na(df_annot$entrezgene_id == TRUE)), ]

  # Create a vector of all genes in the dataset
  universe_genes <- unique(df_annot$entrezgene_id)

  sel_genes <- df_annot[df_annot$adj.P.Val < 0.05 & abs(df_annot$logFC) > 1, ]
  # check value of sel_genes

  if (nrow(sel_genes) < 3) {
    print("no probes selected, no go analysis!")
  } else {
    unique_sel_genes <- unique(sel_genes$entrezgene_id)
    params <- new("GOHyperGParams",
      geneIds = unique_sel_genes,
      universeGeneIds = universe_genes,
      annotation = org_annot,
      ontology = "BP",
      pvalueCutoff = hgCutoff,
      conditional = FALSE,
      testDirection = "over"
    )
    paramsCond <- params
    conditional(paramsCond) <- TRUE

    hgOver <- hyperGTest(params)
    hgCondOver <- hyperGTest(paramsCond)

    summary(hgOver)

    mat_GO <- NULL
    l_value <- NULL

    if (nrow(summary(hgOver)) >= 1) {
      # Get the entrez id for each GO term:
      hgOver_probes <- geneIds(hgOver)
      hgOver_probes <- geneIdsByCategory(hgOver, catids = summary(hgOver)$GOBPID)
      GOtable <- data.frame(
        GOterm = summary(hgOver)$GOBPID, term = summary(hgOver)$Term, Pvalue = format(summary(hgOver)$Pvalue, scientific = T),
        NoGOsize = paste(summary(hgOver)$Count, summary(hgOver)$Size, sep = ";")
      )

      # GOtable = cbind(GOtable,cluster = paste("C",k,nameExp,sep="_"))

      mat_agreg <- NULL

      for (j in 1:length(hgOver_probes)) {
        tmp <- hgOver_probes[[j]]
        tmp
        e_go <- unique_sel_genes[which(is.na(match(unique_sel_genes, tmp)) == FALSE)]
        l_symbol <- paste(e_go, ";", sep = "", collapse = "")
        l_symbol <- substr(l_symbol, 1, nchar(l_symbol) - 1)

        mat_agreg <- rbind(l_symbol, mat_agreg)
      }
      head(mat_agreg)
      row.names(mat_agreg) <- NULL

      colnames(mat_agreg) <- "entrez"
      GOtable <- cbind(GOtable, mat_agreg)
      mat_GO <- rbind(mat_GO, GOtable)
    }

    mat_GO <- cbind(experiment_id = name_exp, mat_GO)
    write.table(mat_GO,
      file = paste0(way_out, "/", name_exp, "/", name_exp, "_go_terms_analysis.csv"),
      sep = ",", col.names = T, row.names = F
    )
  }
}