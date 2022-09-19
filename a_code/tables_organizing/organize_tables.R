# --- ways ---
way_table <- "../tables_hihisiv/"
way_results_microarray <- "../microarray_results/"
way_results_rna_seq <- "../rna-seq_results/"
way_ucsc_id <- "..s/ucsc_id/"


# --- libs ---
library("biomaRt")
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

# ==============
# === TISSUE ===
# ==============
tissue <- read.table(paste0(way_table, "saved_google_sheet/tissue.csv"), sep = ",", header = F)
tissue <- unique(tissue$V1)

write.table(tissue, file = paste0(way_table, "tissue.csv"), row.names = F, col.names = F)
rm(tissue)

# ==============
# === DESIGN ===
# ==============
design <- read.table(paste0(way_table, "saved_google_sheet/design.csv"), sep = ",", header = F)
design <- unique(design$V1)

write.table(design, file = paste0(way_table, "design.csv"), row.names = F, col.names = F)
rm(design)

# ========================
# === RNA-Seq GPL21697 ===
# ========================
ucsc_id <- read.table(file = paste0(way_ucsc_id, "ucsc_id_Hs.csv"), sep = ",", header = F)

# filters <- listFilters(ensembl)
# filters[grep("ucsc", filters$name), ]
# attributes <- listAttributes(ensembl)
# attributes[grep("ucsc", attributes$name), ]

entrez_genes <- getBM(
  filters = "ucsc",
  attributes = c("ucsc", "entrezgene_id"),
  values = ucsc_id,
  mart = mart
)

# --- map the ucsc found and not ---
unmapped <- ucsc_id[which(is.na(match(ucsc_id$V1, entrez_genes$ucsc)) == TRUE), ]
unmapped <- cbind(unmapped, entrezgene_id = NA)
colnames(unmapped)[1] <- "ucsc"

# --- unifying ---
ucsc_id <- rbind(entrez_genes, unmapped)

write.table(ucsc_id, file = paste0(way_table, "platform_transcript_id/GPL21697.csv"), sep = ",", col.names = T, row.names = F)
rm(list = c("ucsc_id", "unmapped", "entrez_genes"))

# ============
# === GENE ===
# ============
ncbi_base <- "https://www.ncbi.nlm.nih.gov/gene/"

# from ncbi and inserted in google sheet
# --- human 9606 ---
human <- read.table(paste0(way_table, "gene_species/gene_human_9606.csv"), sep = ",", head = T)
link <- rep(ncbi_base, nrow(human))
human$link <- paste(link, human$GeneID, sep = "")
gs_9606 <- cbind(human$GeneID, rep(9606, nrow(human)))

# --- macaca fascicularis 9541 ---
fascicularis <- read.table(paste0(way_table, "gene_species/gene_macaca_fascicularis_9541.csv"), sep = ",", head = T)
link <- rep(ncbi_base, nrow(fascicularis))
fascicularis$link <- paste(link, fascicularis$GeneID, sep = "")
gs_9541 <- cbind(fascicularis$GeneID, rep(9541, nrow(fascicularis)))

# --- macaca mulatta 9544 ---
mulatta <- read.table(paste0(way_table, "gene_species/gene_macaca_mulatta_9544.csv"), sep = ",", head = T)
link <- rep(ncbi_base, nrow(mulatta))
mulatta$link <- paste(link, mulatta$GeneID, sep = "")
gs_9544 <- cbind(mulatta$GeneID, rep(9544, nrow(mulatta)))

# --- cercocebus atys 9531 ---
atys <- read.table(paste0(way_table, "gene_species/gene_cercocebus_atys_9531.csv"), sep = ",", head = T)
link <- rep(ncbi_base, nrow(atys))
atys$link <- paste(link, atys$GeneID, sep = "")
gs_9531 <- cbind(atys$GeneID, rep(9531, nrow(atys)))

# --- chlorocebus sabaeus 60711 ---
sabaeus <- read.table(paste0(way_table, "gene_species/gene_chlorocebus_sabaeus_60711.csv"), sep = ",", head = T)
link <- rep(ncbi_base, nrow(sabaeus))
sabaeus$link <- paste(link, sabaeus$GeneID, sep = "")
gs_60711 <- cbind(sabaeus$GeneID, rep(60711, nrow(sabaeus)))

# --- macaca nemestrina 9545---
nemestrina <- read.table(paste0(way_table, "gene_species/gene_macaca_nemestrina_9545.csv"), sep = ",", head = T)
link <- rep(ncbi_base, nrow(nemestrina))
nemestrina$link <- paste(link, nemestrina$GeneID, sep = "")
gs_9545 <- cbind(nemestrina$GeneID, rep(9545, nrow(nemestrina)))

genes <- rbind(human, fascicularis, mulatta, atys, sabaeus, nemestrina)
genes <- genes[!duplicated(genes$GeneID), ]

# ====================
# === GENE_SPECIES ===
# ====================
gene_species <- rbind(gs_9606, gs_9531, gs_9541, gs_9544, gs_9545, gs_60711)
gene_species <- gene_species[!duplicated(gene_species), ]

rm(list = c(
  "gs_60711", "gs_9531", "gs_9541", "gs_9544", "gs_9545", "gs_9606", "human",
  "fascicularis", "mulatta", "nemestrina", "atys", "sabaeus"
))

# ===========================
# === PLATFORM_TRANSCRIPT ===
# ===========================
# --- GPL96 ---
gpl96 <- read.table(paste0(way_table, "platform_transcript_id/GPL96-57554.csv"), sep = "\t", header = T, quote = "")
plat_gpl96 <- rep("GPL96", nrow(gpl96))
gpl96 <- cbind(plat = plat_gpl96, gpl96)
gpl96 <- gpl96[, c(1, 2, 5)]

# --- GPL570 ---
gpl570 <- read.table(paste0(way_table, "platform_transcript_id/GPL570-55999.csv"), sep = "\t", header = T, quote = "")
plat_gpl570 <- rep("GPL570", nrow(gpl570))
gpl570 <- cbind(plat = plat_gpl570, gpl570)
gpl570 <- gpl570[, c(1, 2, 5)]

# --- GPL3535 ---
gpl3535 <- read.table(paste0(way_table, "platform_transcript_id/GPL3535-10024.csv"), sep = "\t", header = T, quote = "", fill = T)
plat_gpl3535 <- rep("GPL3535", nrow(gpl3535))
gpl3535 <- cbind(plat = plat_gpl3535, gpl3535)
gpl3535 <- gpl3535[, c(1, 2, 5)]

# --- GPL3921 ---
gpl3921 <- read.table(paste0(way_table, "platform_transcript_id/GPL3921-25447.csv"), sep = "\t", header = T, quote = "")
plat_gpl3921 <- rep("GPL3921", nrow(gpl3921))
gpl3921 <- cbind(plat = plat_gpl3921, gpl3921)
gpl3921 <- gpl3921[, c(1, 2, 5)]

# --- GPL2986 ---
gpl2986 <- read.table(paste0(way_table, "platform_transcript_id/GPL2986.csv"), sep = "\t", header = T, quote = "")
plat_gpl2986 <- rep("GPL2986", nrow(gpl2986))
gpl2986 <- cbind(plat = plat_gpl2986, gpl2986)
gpl2986 <- gpl2986[, c(1, 2, 5)]

# --- CREATE TABLES PLATFORM AND ID RNA-SEQ ---
gpl21697 <- read.table(paste0(way_table, "platform_transcript_id/GPL21697.csv"), sep = ",", header = T)
colnames(gpl21697) <- c("ID", "ENTREZ_GENE_ID")
plat_gpl21697 <- rep("GPL21697", nrow(gpl21697))
gpl21697 <- cbind(plat = plat_gpl21697, gpl21697)

# --- PLATFORM PROBE ---
plat_transc_id_full <- data.frame(rbind(gpl96, gpl570, gpl3535, gpl3921, gpl2986, gpl21697))
plat_transc_id <- plat_transc_id_full[, c(1, 2)]
plat_transc_id <- plat_transc_id[!duplicated(plat_transc_id), ]
write.table(plat_transc_id, file = paste0(way_table, "platform_transcript.csv"), sep = ",", col.names = F, row.names = F)

# ==================
# === TRANSCRIPT ===
# ==================
# --- UNIQUES TRANSCRIPTS IDS (PROBES AND RNA-SEQ ID) ---
transcript_ids <- unique(plat_transc_id$ID)
write.table(transcript_ids, file = paste0(way_table, "transcript.csv"), sep = ",", col.names = F, row.names = F)
rm(list = c(
  "gpl570", "gpl96", "gpl21697", "gpl2986", "gpl3535", "gpl3921", "plat_gpl21697",
  "plat_gpl2986", "plat_gpl3921", "plat_gpl3535", "plat_gpl570", "plat_gpl96"
))

# =======================
# === TRANSCRIPT_GENE ===
# =======================
# --- clean the data: removing probe without entrez_id ---
plat_transc_id_full <- plat_transc_id_full[!(is.na(plat_transc_id_full$ENTREZ_GENE_ID) | plat_transc_id_full$ENTREZ_GENE_ID == ""), ]

# --- split two df with unique entrez and not ---
plat_transc_id_multi_entrez <- plat_transc_id_full[grep("///", plat_transc_id_full$ENTREZ_GENE_ID), ]
plat_transc_id_unique_entrez <- plat_transc_id_full[-grep("///", plat_transc_id_full$ENTREZ_GENE_ID), ]

# --- split entrez, gene symbol and description gene ---
entrez_split <- strsplit(plat_transc_id_multi_entrez[, "ENTREZ_GENE_ID"], split = "///")
names(entrez_split) <- plat_transc_id_multi_entrez$ID

transcript_gene <- NULL
for (i in 1:length(entrez_split)) {
  ID <- rep(names(entrez_split[i]), length(entrez_split[[i]]))
  transcript_gene <- rbind(transcript_gene, cbind(plat = plat_transc_id_full$plat[i], ID, ENTREZ_GENE_ID = entrez_split[[i]]))
}
transcript_gene <- rbind(plat_transc_id_unique_entrez, transcript_gene)

# --- removing duplicates ---
transcript_gene$ID <- as.character(transcript_gene$ID)
transcript_gene$ENTREZ_GENE_ID <- as.numeric(transcript_gene$ENTREZ_GENE_ID)
# considero so o par de chaves probe-entrez duplicados
transcript_gene <- transcript_gene[!duplicated(transcript_gene[, c(2, 3)]), ]
# transcript_gene2 <- transcript_gene[-which(is.na(transcript_gene$ENTREZ_GENE_ID) == TRUE), ]

write.table(transcript_gene[, c(2, 3)], file = paste0(way_table, "transcript_gene.csv"), sep = ",", col.names = F, row.names = F)

rm(list = c(
  "plat_transc_id", "entrez_split",
  "plat_transc_id_multi_entrez", "plat_transc_id_unique_entrez", "plat_transc_id_full"
))




















# ================
# === ANALYSIS ===
# ================
# *** MICROARRAY DATA ***
list_experiments <- list.files(way_results_microarray)

all_data_micro <- NULL
i <- 1
for (i in 1:length(list_experiments)) {
  exp <- read.table(paste0(
    way_results_microarray, list_experiments[i], "/", list_experiments[i],
    "_diffExpression.csv"
  ), sep = ",", header = T)
  if (grepl("2932572286", list_experiments[i]) == TRUE || grepl("3070984318", list_experiments[i]) == TRUE){
    exp <- exp[, -2]
    temp <- cbind(exp, experiment_id = list_experiments[i], norm.method = "Log2", analysis_script = NA)
  }
  else {
    temp <- cbind(exp, experiment_id = list_experiments[i], norm.method = "MAS5", analysis_script = NA)
  }
  # head(temp)
  all_data_micro <- rbind(all_data_micro, temp)
  # tail(all_data_micro)
}

# --- reorganizing columns ---
all_data_micro <- all_data_micro[, c("experiment_id", "ID", "norm.method", "logFC", "adj.P.Val", "analysis_script")]

# --- check if the probe_id is ok ---
# probes <- read.table(paste0(way_table,"transcript.csv"), sep = ",", header = F)
all_data_micro <- all_data_micro[-which(is.na(match(all_data_micro$ID, transcript_ids)) == TRUE), ]

# *** RNA-SEQ DATA ***
list_experiments <- list.files(way_results_rna_seq)
list_experiments <- list_experiments[grep(".csv", list_experiments)]

all_data_rna_seq <- NULL
for (i in 1:length(list_experiments)) {
  exp <- read.table(paste0(way_results_rna_seq, list_experiments[i]), sep = ",", header = T)
  all_data_rna_seq <- rbind(all_data_rna_seq, exp)
}

# --- SAVING ---
all_data_analysis <- rbind(all_data_micro, all_data_rna_seq)
write.table(all_data_analysis, file = paste0(way_table, "analysis.csv"), sep = ",", col.names = F, row.names = F)

rm(list = c("all_data_analysis", "all_data_micro", "all_data_rna_seq", "exp", "temp"))

# ==================
# === GENE_TRAIT ===
# ==================
# --- CELL SURFACE  ---
cell_surface_ncbi <- read.table(paste0(way_table, "traits/cell_surface_NCBI.csv"), sep = ",", header = T)
cell_surface_ncbi <- cell_surface_ncbi[, c(2, 1)]
cell_surface_uniprot <- read.table(paste0(way_table, "traits/cell_surface_UniProtKB.csv"), sep = ",", header = T)
values <- cell_surface_uniprot$Symbol

filter_trait <- "hgnc_symbol"
attributes_traits <- c("entrezgene_id", "hgnc_symbol", "description")

entrez_cell_surface <- getBM(
  filters = filter_trait,
  attributes = attributes_traits,
  values = values,
  mart = mart
)

entrez_cell_surface <- entrez_cell_surface[-which(is.na(entrez_cell_surface$entrezgene_id) == TRUE), ]
cell_surface_uniprot <- cbind(
  GeneID = entrez_cell_surface$entrezgene_id,
  trait_id = rep("cell_surface_UniProtKB", nrow(entrez_cell_surface))
)

# --- CYTOKINE  ---
cytokine_ncbi <- read.table(paste0(way_table, "traits/cytokine_NCBI.csv"), sep = ",", header = T)
cytokine_ncbi <- cytokine_ncbi[, c(2, 1)]
cytokine_uniprot <- read.table(paste0(way_table, "traits/cytokine_UniProtKB.csv"), sep = ",", header = T)
values <- cytokine_uniprot$Symbol

entrez_cytokine <- getBM(
  filters = filter_trait,
  attributes = attributes_traits,
  values = values,
  mart = mart
)

entrez_cytokine <- entrez_cytokine[-which(is.na(entrez_cytokine$entrezgene_id) == TRUE), ]
cytokine_uniprot <- cbind(
  GeneID = entrez_cytokine$entrezgene_id,
  trait_id = rep("cytokine_UniProtKB", nrow(entrez_cytokine))
)

# --- RESTRICTION FACTORS  ---
restriction_factors_ncbi <- read.table(paste0(way_table, "traits/restriction_factor_NCBI.csv"), sep = ",", header = T)
restriction_factors_ncbi <- restriction_factors_ncbi[, c(2, 1)]

# --- TRANSCRIPTION FACTORS  ---
transcription_factors_ncbi <- read.table(paste0(way_table, "traits/transcription_factors_NCBI.csv"), sep = ",", header = T)
transcription_factors_ncbi <- transcription_factors_ncbi[, c(1, 2)]

# --- UNIFYING DATA ---
gene_trait <- rbind(
  cell_surface_ncbi,
  cell_surface_uniprot,
  cytokine_ncbi,
  cytokine_uniprot,
  restriction_factors_ncbi,
  transcription_factors_ncbi
)

# --- rsaving ---
gene_trait <- gene_trait[!duplicated(gene_trait), ]
write.table(gene_trait, file = paste0(way_table, "gene_trait.csv"), sep = ",", col.names = F, row.names = F)

rm(list = c(
  "cell_surface_ncbi", "cell_surface_uniprot", "cytokine_ncbi", "cytokine_uniprot",
  "transcription_factors_ncbi", "restriction_factors_ncbi", "values", "entrez_cytokine", "entrez_cell_surface"
))

# ==================
# === GENE_GENE  ===
# ==================
# *** PPI ***
# information from String database
ppi <- read.table(paste0(way_table, "gene_gene/9606.protein.links.v11.0.txt"), sep = " ", header = T)

# --- REMOVING COMBINATION WITH LOW SCORE COMBINED ---
ppi <- ppi[-which(ppi$combined_score < 300), ]

# --- get ENSP_id ---
only_id <- unique(ppi$protein1)
only_id_2 <- unique(ppi$protein2)
only_id <- unique(c(only_id, only_id_2))
only_id <- gsub("9606.", "", only_id)
rm(only_id_2)

# === BioMart ===
# ensembl <- useMart("ensembl")
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
#
# attributes <- listAttributes(ensembl)
# attributes[grep("gene", attributes$name), ]

entrez_genes <- getBM(
  filters = "ensembl_peptide_id",
  attributes = c("ensembl_peptide_id", "entrezgene_id"),
  values = only_id,
  mart = mart
)

ppi$protein1 <- gsub("9606.", "", ppi$protein1)
ppi$protein2 <- gsub("9606.", "", ppi$protein2)

ppi$protein1 <- entrez_genes[match(ppi$protein1, entrez_genes$ensembl_peptide_id), 2]
ppi$protein2 <- entrez_genes[match(ppi$protein2, entrez_genes$ensembl_peptide_id), 2]

# --- removing NA ---
ppi <- ppi[-which(is.na(ppi$protein1)), ]
ppi <- ppi[-which(is.na(ppi$protein2)), ]

# --- removing duplicated ---
dim(ppi)
# [1] 3057360       3
# [1] 3048161       5

ppi$relationship <- rep("PPI", nrow(ppi))
ppi$source <- rep("STRING", nrow(ppi))
ppi <- ppi[, c(1, 2, 4, 5, 3)]

colnames(ppi) <- c("V1", "V2", "V3", "V4", "V5")

# --- remover links duplicated ---
want <- which(ppi$V1 > ppi$V2)
ppi_change <- ppi[want, c(2, 1, 3, 4, 5)]
ppi <- ppi[-want, ]
ppi <- rbind(ppi, ppi_change)
ppi <- ppi[!duplicated(ppi[, c(1, 2)]), ]







# *** ORTHOLOGOUS AND HUMANS AND MACACA MULATTA ***

# new_config <- httr::config(ssl_verifypeer = FALSE)
# httr::set_config(new_config, override = FALSE)
# devtools::check(vignettes = FALSE, args = '--timings')

ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

genes_id_mm <- gene_species[which(gene_species[, 2] == 9544), 1]

ensembl_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
macaque_ensembl <- useDataset("mmulatta_gene_ensembl", mart = ensembl)


orthologos <- getLDS(
  attributes = "entrezgene_id",
  filters = "entrezgene_id",
  values = genes_id_mm,
  mart = macaque_ensembl,
  attributesL = c("hgnc_symbol", "entrezgene_id"),
  martL = ensembl_human
)


packageVersion("biomaRt") # [1] ‘2.52.0’

# BiocManager::install('grimbough/biomaRt', ref = 'RELEASE_3_14')
# library(biomaRt)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

listEnsembl()


ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
searchDatasets(mart = ensembl, pattern = "hsapiens")
# 81 hsapiens_gene_ensembl Human genes (GRCh38.p13) GRCh38.p13
searchDatasets(mart = ensembl, pattern = "mulatta")
# 105 mmulatta_gene_ensembl Macaque genes (Mmul_10) Mmul_10



orthologos <- orthologos[, c(3, 1)]
colnames(orthologos) <- c("entrez_humano", "entrez_mulatta")
orthologos <- orthologos[-which(is.na(orthologos$entrez_mulatta) == TRUE | is.na(orthologos$entrez_humano) == TRUE), ]

orthologos <- orthologos[!duplicated(orthologos), ]
orthologos <- cbind(orthologos, type = "orthologous", source = "biomaRt", obs = "homo sapiens to macaca mulatta")
colnames(orthologos) <- c("V1", "V2", "V3", "V4", "V5")

gene_gene <- rbind(ppi, orthologos)

write.table(gene_gene, file = paste0(way_table, "gene_gene.csv"), sep = ",", col.names = F, row.names = F)
rm(list = c("gene_gene", "genes_id_mm", "entrez_genes", "ppi_change", "only_id", "want"))

# --- other orthologou source ---<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<,
orthologos_2_source <- read.table(paste0(way_table, "gene_gene/gene_gene_ncbi.csv"), sep = ",", header = T)
# head(orthologos_2_source)
# head(orthologos)

# 8888888888888888888
# === REVIEW GENE ===
# 8888888888888888888
# --- transcript_id_genes ---
miss_genes_transcript_gene <- transcript_gene[which(is.na(match(transcript_gene$ENTREZ_GENE_ID, genes$GeneID)) == TRUE), ]

# --- trait ---
miss_genes_traits <- unique(gene_trait$GeneID)
miss_genes_traits <- miss_genes_traits[which(is.na(match(miss_genes_traits, genes$GeneID)) == TRUE)]

# --- ppi ---
miss_genes_ppi <- unique(c(unique(ppi$V1), unique(ppi$V2)))
miss_genes_ppi <- miss_genes_ppi[which(is.na(match(miss_genes_ppi, genes$GeneID)) == TRUE)]

# --- orthologous ---
miss_genes_orth <- unique(c(unique(orthologos$V1), unique(orthologos$V2)))
miss_genes_orth <- miss_genes_orth[which(is.na(match(miss_genes_orth, genes$GeneID)) == TRUE)]

miss_genes <- c(miss_genes_transcript_gene$ENTREZ_GENE_ID, miss_genes_traits, miss_genes_ppi, miss_genes_orth)
miss_genes <- unique(miss_genes)

miss_genes <- as.numeric(miss_genes)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes_info <- getBM(
  filters = "entrezgene_id",
  attributes = c("entrezgene_id", "description", "hgnc_symbol"),
  values = miss_genes,
  mart = mart
)

# which(genes_info$entrezgene_id==100129575)
genes_info <- genes_info[!duplicated(genes_info$entrezgene_id), ]

# miss_yet = miss_genes[which(is.na(match(genes_info$entrezgene_id,miss_genes))==FALSE)]
miss_yet <- miss_genes[which(is.na(match(miss_genes, genes_info$entrezgene_id)) == TRUE)]

ncbi_base <- "https://www.ncbi.nlm.nih.gov/gene/"
link <- paste(rep(ncbi_base, nrow(genes_info)), genes_info$entrezgene_id, sep = "")
genes_info_2 <- as.data.frame(cbind(
  entrezgene_id = genes_info$entrezgene_id,
  gene_symbol = genes_info$hgnc_symbol,
  link = link,
  description = genes_info$description,
  chr = NA
))

miss_yet_2 <- cbind(
  entrezgene_id = miss_yet,
  gene_symbol = NA,
  link = NA,
  description = NA,
  chr = NA
)

new_genes <- rbind(genes_info_2, miss_yet_2)
colnames(new_genes) <- c("GeneID", "Symbol", "link", "description", "chromosome")
new_genes <- as.data.frame(new_genes)

genes <- rbind(genes, new_genes)
genes <- genes[!duplicated(genes$GeneID), ]

write.table(genes, file = paste0(way_table, "gene.csv"), sep = ",", col.names = F, row.names = F)

# 888888888888888888888
# --- GENE_SPECIES ---
# 888888888888888888888
# --- probe_gene ---
miss_genes_transcript_gene <- miss_genes_transcript_gene[, c(3, 1)]

miss_genes_transcript_gene$plat <- gsub("GPL96", "9606", miss_genes_transcript_gene$plat)
miss_genes_transcript_gene$plat <- gsub("GPL570", "9606", miss_genes_transcript_gene$plat)
miss_genes_transcript_gene$plat <- gsub("GPL3535", "9544", miss_genes_transcript_gene$plat)
miss_genes_transcript_gene$plat <- gsub("GPL3921", "9606", miss_genes_transcript_gene$plat)
miss_genes_transcript_gene$plat <- gsub("GPL2986", "9606", miss_genes_transcript_gene$plat)
miss_genes_transcript_gene$plat <- gsub("GPL21697", "9606", miss_genes_transcript_gene$plat)
colnames(miss_genes_transcript_gene) <- c("V1", "V2")
gene_species <- rbind(gene_species, miss_genes_transcript_gene)
gene_species <- gene_species[!duplicated(gene_species), ]

# --- gene_sp_trait ---
new_genes_species_trait <- cbind(V1 = miss_genes_traits, V2 = "9606")
gene_species <- rbind(gene_species, new_genes_species_trait)
gene_species <- gene_species[!duplicated(gene_species), ]

# --- ortologos ---
new_genes_species_orth <- cbind(V1 = miss_genes_orth, V2 = "9606")
gene_species <- rbind(gene_species, new_genes_species_orth)
gene_species <- gene_species[!duplicated(gene_species), ]

write.table(gene_species, file = paste0(way_table, "gene_species.csv"), sep = ",", col.names = F, row.names = F)

# ==========
# === GO ===
# ==========
genes_id <- as.numeric(genes$GeneID)

ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

# --- GO ---
go_genes <- getBM(
  filters = "entrezgene_id",
  attributes = c("entrezgene_id", "go_id", "name_1006", "namespace_1003"),
  values = genes_id,
  mart = mart
)

go_genes <- go_genes[-which(go_genes$go_id == ""), ]
go <- go_genes[, 2:4]
go <- go[!duplicated(go), ]

# --- write tables ---
write.table(go, file = paste0(way_table, "go.csv"), sep = ",", col.names = F, row.names = F)

# ================
# === GENES_GO ===
# ================
write.table(go_genes[, 1:2], file = paste0(way_table, "gene_go.csv"), sep = ",", col.names = F, row.names = F)
