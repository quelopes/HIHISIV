my_data <- read.csv("/data.csv", header = F)
head(my_data)
# setwd("../")
getwd()

my_data
# --- ways ---
raw_data <- "/data_microarray/"
way_result <- "/results"
way_code <- "/hihisiv_gitlab/a_code/microarray_analysis/"
source(paste0(way_code, "dependencies.R"))
dependencies(way_code)

i <- 1
for (i in 1:nrow(my_data)) {
  project <- my_data[i, 4]
  experiment_name <- my_data[i, 1]
  org_annot <- my_data[i, 2]
  
  # --- create a dir with the results ---
  if (file.exists(experiment_name)) {
    print("the directory already exists")
  } else {
    dir.create(file.path(way_result, experiment_name))
  }
  
  norm_method <- "mas5"

  # --- e-set ---
  pheno_data <- read.table(paste0(raw_data, project, "/", experiment_name, ".csv"), head = T, sep = ",")
  
  # ======================
  # === GERAL SETTING ===
  # ======================
  # --- limma ---
  fdr <- 0.05 # threshold
  log_fc <- 1 # threshold
  
  org_annot <- my_data$V2[i]
  
  print("==========================")
  print(experiment_name)
  print("==========================")
  
  # === NORMALIZED ===
  mat_normalized <- raw_activity(raw_data, way_result, pheno_data, norm_method, experiment_name)

  # === E-SET ===
  e_set <- e_set_activity(mat_normalized, pheno_data, org_annot)
  
  rm(mat_normalized)
  
  # === LIMMA ===
  result_limma <- limma_activity(e_set, fdr, log_fc, way_result, experiment_name)
  result_limma <- result_limma[, -2]
  
  # === ANNOTATION ===
  affyids <- result_limma$ID
  
  if (org_annot == "org.Mmu.eg.db") {
    ensembl <- useEnsembl(biomart = "genes", dataset = "mmulatta_gene_ensembl")
    result_annot <- getBM(
      attributes = c("affy_rhesus", "entrezgene_id", "hgnc_symbol"),
      filters = "affy_rhesus",
      values = affyids,
      mart = ensembl
    )
  } else {
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl") # , host="https://uswest.ensembl.org")
    result_annot <- getBM(
      attributes = c(my_data$V5[i], "entrezgene_id", "hgnc_symbol"),
      filters = my_data$V5[i],
      values = affyids,
      mart = ensembl
    )
  }
  
  colnames(result_annot)[1] <- "probe_id"
  
  with_annot <- merge(result_annot, result_limma, by.x = "probe_id", by.y = "ID")
  with_annot <- cbind(experiment_id = experiment_name, with_annot)
  write.table(with_annot,
              file = paste0(way_result, "/", experiment_name, "/", experiment_name, "_diffExpression.csv"),
              sep = ",", col.names = T, row.names = F
  )
  
  rm(result_annot)
  rm(e_set)
  rm(result_annot)
  
  # === GOstats ===
  go_terms_analysis(with_annot, org_annot, experiment_name, way_result)
  rm(go_terms)
  rm(with_annot)
}

# ===========================
# --- PLATFORM TRANSCRIPT ---
# ===========================
my_data <- read.csv("/data.csv", header = F)
dir_platforms <- "/hihisiv_gitlab/platforms/"
list.files(dir_platforms)
library(biomaRt)

my_data_plat <- my_data[, -c(1, 4)]
my_data_plat <- my_data_plat[!duplicated(my_data_plat), ]
platforms <- unique(my_data_plat$V3)

plat_transcript <- NULL
transcripts = NULL

j=5
for (j in 1:length(platforms)) {
  if (platforms[j] == "GPL21697") {
    df_plat <- read.csv(paste0(dir_platforms, platforms[j], ".csv"), header = T)
    head(df_plat)
    ucsc_ids <- df_plat$ucsc
    transcripts = c(transcripts, df_plat$ucsc)
    
    
    plat_transcript_temp <- getBM(
      attributes = c("ucsc", "entrezgene_id", "hgnc_symbol"),
      filters = "ucsc",
      values = ucsc_ids,
      mart = ensembl
    )
    head(plat_transcript_temp)
    # ucsc entrezgene_id hgnc_symbol
    # 1 uc001abu.1            NA
    # 2 uc001aef.2         54973      INTS11
    # 3 uc001aei.2         54973      INTS11
    # 4 uc001afp.3        441869     ANKRD65
    # 5 uc001aic.3         65220        NADK
    # 6 uc001aie.3         65220        NADK
    plat_transcript_temp$species <- "Homo sapiens"
  } else {
    df_plat <- read.table(paste0(dir_platforms, platforms[j], ".csv"), header = T, sep = "\t", quote = "", fill = T)
    head(df_plat)
    transcripts = c(transcripts, df_plat$ID)
    # df_plat[grep("213274_s_at",df_plat$ID),]
    
    probes_plat <- df_plat$ID
    
    if (platforms[j] == "GPL3535") {
      ensembl <- useEnsembl(biomart = "genes", dataset = "mmulatta_gene_ensembl")
      plat_transcript_temp <- getBM(
        attributes = c("affy_rhesus", "entrezgene_id", "hgnc_symbol"),
        filters = "affy_rhesus",
        values = probes_plat,
        mart = ensembl
      )
      plat_transcript_temp$species <- "Macaca mulatta"
      # --- map entrez macaca mulatta ---
      entrez_mulatta <- plat_transcript_temp$entrezgene_id
      entrez_mulatta <- entrez_mulatta[which(is.na(entrez_mulatta) == FALSE)]
    } else {
      ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl") # , host="https://uswest.ensembl.org")
      plat_transcript_temp <- getBM(
        attributes = c(my_data_plat[j, 3], "entrezgene_id", "hgnc_symbol"),
        filters = my_data_plat[j, 3],
        values = probes_plat,
        mart = ensembl
      )
      head(plat_transcript_temp)
      plat_transcript_temp$species <- "Homo sapiens"
    }
  }
  # grep("213274_s_at",plat_transcript_temp$affy_hg_u133a_2)
  
  plat_transcript_temp <- cbind(platform = platforms[j], plat_transcript_temp)
  colnames(plat_transcript_temp)[2] <- "probe"
  plat_transcript <- rbind(plat_transcript, plat_transcript_temp)
}

link_entrez_id <- "https://www.ncbi.nlm.nih.gov/gene/"
temp_link_entrez <- paste0(link_entrez_id, plat_transcript$entrezgene_id)

link_genecards <- "https://www.genecards.org/cgi-bin/carddisp.pl?gene="
temp_link_symbol <- paste0(link_genecards, plat_transcript$hgnc_symbol)

plat_transcript$entrez_id_link <- temp_link_entrez
plat_transcript$hgnc_id_link <- temp_link_symbol

# grep("213274_s_at",plat_transcript$probe)

head(plat_transcript)
dim(plat_transcript) # [1] 136189      7

# --- write results ---
head(plat_transcript)

plat_trans_final <- plat_transcript[, c(1, 2)]
# plat_trans_final[grep("201104_x_at", plat_trans_final[,2]),]
plat_trans_final <- plat_trans_final[!duplicated(plat_trans_final), ]
head(plat_trans_final)

head(plat_transcript)
transcript_table <- unique(plat_transcript$probe)
transcript_table <- unique(transcripts)
length(transcript_table) # [1] 77861 # [1] 185652



entrez_id_table <- plat_transcript[, c("entrezgene_id", "species", "entrez_id_link")]
entrez_id_table <- entrez_id_table[!duplicated(entrez_id_table), ]
entrez_id_table <- entrez_id_table[which(is.na(entrez_id_table$entrezgene_id) == FALSE), ]
head(entrez_id_table)


entrez_id_table[which(entrez_id_table$entrezgene_id==692213),]

symbol <- plat_transcript[, c("hgnc_symbol", "hgnc_id_link")]
symbol$hgnc_symbol[symbol$hgnc_symbol == ""] <- NA
symbol <- symbol[which(is.na(symbol$hgnc_symbol) == FALSE), ]
symbol <- symbol[!duplicated(symbol$hgnc_symbol), ]
# symbol[which(symbol$hgnc_symbol=='CCL3L1'),]


entrez_symbol <- plat_transcript[, c("entrezgene_id", "hgnc_symbol")]
entrez_symbol$hgnc_symbol[entrez_symbol$hgnc_symbol == ""] <- NA
entrez_symbol$entrezgene_id[entrez_symbol$entrezgene_id == ""] <- NA
entrez_symbol <- entrez_symbol[!duplicated(entrez_symbol), ]
entrez_symbol <- entrez_symbol[which(is.na(entrez_symbol$hgnc_symbol) == FALSE), ]
entrez_symbol <- entrez_symbol[which(is.na(entrez_symbol$entrezgene_id) == FALSE), ]

transc_gene_entrez <- plat_transcript[, c("probe", "entrezgene_id"), ]
transc_gene_entrez <- transc_gene_entrez[!duplicated(transc_gene_entrez), ]
transc_gene_entrez <- transc_gene_entrez[which(is.na(transc_gene_entrez$entrezgene_id) == FALSE), ]
head(transc_gene_entrez)

write.table(plat_trans_final, paste0("/b_database/tables/", "platform_transcript.csv"), row.names = F, col.names = T, sep = ",")
write.table(transcript_table, paste0("/b_database/tables/", "transcript.csv"), row.names = F, col.names = T, sep = ",")
write.table(entrez_id_table, paste0("/b_database/tables/", "gene_entrez.csv"), row.names = F, col.names = T, sep = ",")
write.table(symbol, paste0("/b_database/tables/", "gene_symbol.csv"), row.names = F, col.names = T, sep = ",")
write.table(entrez_symbol, paste0("/b_database/tables/", "entrez_symbol.csv"), row.names = F, col.names = T, sep = ",")
write.table(transc_gene_entrez, paste0("/b_database/tables/", "transcript_gene_entrez.csv"), row.names = F, col.names = T, sep = ",")



# grep("213274_s_at" ,transcript_table)

# --------------------------------------------------------
# --- ORTHOLOGS ---
human.mart <- biomaRt::useMart(host = "https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
macaca.mart <- biomaRt::useMart(host = "https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset = "mmulatta_gene_ensembl")

orthologs <- getLDS(
  attributes = "entrezgene_id",
  filters = "entrezgene_id",
  values = entrez_mulatta,
  mart = macaca.mart,
  # macaque_ensembl,
  attributesL = c("entrezgene_id"),
  martL = human.mart
  # ensembl_human
)

colnames(orthologs) <- c("entrez_id_mulatta", "entrez_id_human")
head(orthologs)

head(entrez_id_table)

non_matching_rows <- orthologs[!(orthologs$entrez_id_human %in% entrez_id_table$entrezgene_id), ]
non_matching_rows <- unique(non_matching_rows$entrez_id_human)

temp_link <- paste0(link_entrez_id, non_matching_rows)
non_matching_entrez <- cbind(entrez_gene_id = non_matching_rows, species = "Homo sapiens", entrez_id_link = temp_link)
non_matching_entrez <- as.data.frame(non_matching_entrez)
head(non_matching_entrez)

write.table(non_matching_entrez, paste0("/b_database/tables/", "additional_gene_entrez.csv"), row.names = F, col.names = T, sep = ",")
write.table(orthologs, paste0("/b_database/tables/", "orthologs.csv"), row.names = F, col.names = T, sep = ",")

# --------------------------------------------------------
# ver[grep("abi", ver$description),]

ver <- listAttributes(human.mart)
ver[grep("ulatta", ver$name), ]
# 1901                 mmulatta_homolog_ensembl_gene                           Macaque gene stable ID homologs
# 1902         mmulatta_homolog_associated_gene_name                                Macaque gene name homologs
# 1903              mmulatta_homolog_ensembl_peptide          Macaque protein or transcript stable ID homologs
# 1904                   mmulatta_homolog_chromosome                 Macaque chromosome/scaffold name homologs
# 1905                  mmulatta_homolog_chrom_start           Macaque chromosome/scaffold start (bp) homologs
# 1906                    mmulatta_homolog_chrom_end             Macaque chromosome/scaffold end (bp) homologs
# 1907 mmulatta_homolog_canonical_transcript_protein                   Query protein or transcript ID homologs
# 1908                      mmulatta_homolog_subtype                Last common ancestor with Macaque homologs
# 1909               mmulatta_homolog_orthology_type                            Macaque homology type homologs
# 1910                      mmulatta_homolog_perc_id %id. target Macaque gene identical to query gene homologs
# 1911                   mmulatta_homolog_perc_id_r1 %id. query gene identical to target Macaque gene homologs
# 1912                    mmulatta_homolog_goc_score            Macaque Gene-order conservation score homologs
# 1913                 mmulatta_homolog_wga_coverage          Macaque Whole-genome alignment coverage homologs
# 1914         mmulatta_homolog_orthology_confidence     Macaque orthology confidence [0 low, 1 high] homologs
