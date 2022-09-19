# =============================================
# ============ MEUS DADOS =====================
# =============================================
library(tximport)
library(tximportData)
library(readr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(limma)
library(edgeR)

name = "GSE119234_uninfected_infected"
# name = "GSE119234_B-cell-germinal_uninfected_infected"
# name = "GSE119234_B-cell-memory_uninfected_infected"
# name = "GSE119234_B-cell-naive_uninfected_infected"
# name = "GSE119234_B-cell-unswitched-memory_uninfected_infected"

# --- Transcript id to gene id ---
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# --- RNA-Seq results from RSEM ---
dir <- "../processed/rna-seq_results/SRR_results/"
list.files(dir)

# --- Metadata ---
metaData <- read.table(paste("../rna-seq_analysis/get_SRA/", name, ".csv", sep = ""), sep = ",", header = T)
head(metaData)
metaData$Cell_type
rownames(metaData) <- metaData$Run
metaData$infection <- gsub("HIV-infected", "Infected", metaData$infection)


files <- file.path(dir, paste0(metaData$Run, ".genes.results"))
names(files) <- metaData$Run

txi.rsem <- tximport(files, type = "rsem", tx2gene = tx2gene, txIn = FALSE, txOut = FALSE)
head(txi.rsem)
head(txi.rsem$abundance)
y <- DGEList(txi.rsem$counts)
y <- calcNormFactors(y)

f <- factor(metaData$infection)
# metaData$infection = as.factor(metaData$infection)

# Pq da diferença?
design <- model.matrix(~infection, data = metaData) # do manual do limma
# (Intercept) infectionUninfected
# SRR7768809           1                   1
# SRR7768810           1                   1
# SRR7768811           1                   1
# SRR7768812           1                   1
# SRR7768813           1                   1
# SRR7768814           1                   1
# SRR7768815           1                   1
# SRR7768816           1                   1
# SRR7768817           1                   1
# SRR7768818           1                   1
# SRR7768819           1                   1
# SRR7768820           1                   1
# SRR7768821           1                   1
# SRR7768822           1                   1
# SRR7768823           1                   1
# SRR7768824           1                   1
# SRR7768825           1                   1
# SRR7768826           1                   1
# SRR7768827           1                   1
# SRR7768828           1                   1
# SRR7768829           1                   0
# SRR7768830           1                   0
# SRR7768831           1                   0
# SRR7768832           1                   0
# SRR7768833           1                   0
# SRR7768834           1                   0
# SRR7768835           1                   0
# SRR7768836           1                   0
# SRR7768837           1                   0
# SRR7768838           1                   0
# SRR7768839           1                   0
# SRR7768840           1                   0
# SRR7768841           1                   0
# SRR7768842           1                   0
# SRR7768843           1                   0
# SRR7768844           1                   0
# SRR7768845           1                   0
# SRR7768846           1                   0
# SRR7768847           1                   0
# SRR7768848           1                   0
# SRR7768849           1                   0
# SRR7768850           1                   0
# SRR7768851           1                   0
# SRR7768852           1                   0
# SRR7768853           1                   0
# SRR7768854           1                   0
# SRR7768855           1                   0
# SRR7768856           1                   0
# SRR7768857           1                   0
# SRR7768858           1                   0
# SRR7768859           1                   0
# attr(,"assign")
# [1] 0 1
# attr(,"contrasts")
# attr(,"contrasts")$infection
# [1] "contr.treatment"

# design = model.matrix(~-1+f) # script antigo do Marcelo
# fInfected fUninfected
# 1          0           1
# 2          0           1
# 3          0           1
# 4          0           1
# 5          0           1
# 6          0           1
# 7          0           1
# 8          0           1
# 9          0           1
# 10         0           1
# 11         0           1
# 12         0           1
# 13         0           1
# 14         0           1
# 15         0           1
# 16         0           1
# 17         0           1
# 18         0           1
# 19         0           1
# 20         0           1
# 21         1           0
# 22         1           0
# 23         1           0
# 24         1           0
# 25         1           0
# 26         1           0
# 27         1           0
# 28         1           0
# 29         1           0
# 30         1           0
# 31         1           0
# 32         1           0
# 33         1           0
# 34         1           0
# 35         1           0
# 36         1           0
# 37         1           0
# 38         1           0
# 39         1           0
# 40         1           0
# 41         1           0
# 42         1           0
# 43         1           0
# 44         1           0
# 45         1           0
# 46         1           0
# 47         1           0
# 48         1           0
# 49         1           0
# 50         1           0
# 51         1           0
# attr(,"assign")
# [1] 1 1
# attr(,"contrasts")
# attr(,"contrasts")$f
# [1] "contr.treatment"

v <- voom(y, design, plot = FALSE)
# v is now ready for lmFit() see limma User's Guide

# save.image(file = "limma-voom.RData")


# --- Using limma ---
fit <- lmFit(v, design)
fit <- eBayes(fit, trend = TRUE)
mat <- topTable(fit, coef = ncol(design), number = nrow(y$counts))
dim(mat[which(mat$adj.P.Val < 0.05 & mat$logFC > 1.5), ])
# [1] 21  9 manual
# [1] 14552     6 Marcelo script

mat$ID <- rownames(mat)
mat$experiment_id <- name
mat$norm.method <- NA
mat$analysis_script <- NA
mat <- mat[, c(8, 7, 9, 1, 5, 10)]
head(mat)

write.table(mat, file = paste("../RNA-Seq_results/", name, ".csv", sep = ""), row.names = F, col.names = T, sep = ",")