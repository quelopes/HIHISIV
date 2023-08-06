# ====================
# --- Dependencies ---
# ====================
library(stringr)
library(readr)
library(GEOquery)


way <- "/extract/"
source(paste0(way, "extraction_gse.R"))
source(paste0(way, "filter_gse.R"))
source(paste0(way, "geo_query.R"))


my_data <- read.csv("/data.csv", header = F)

gse_hihisiv <- unique(my_data$V4)
gse_hihisiv <- gse_hihisiv[-c(14, 15)]

gse_hihisiv <- c(gse_hihisiv, "GSE119234")

df_gse_complete <- geo_query(gse_hihisiv)

df_gse_complete <- data.frame(df_gse_complete, stringsAsFactors = FALSE)
head(df_gse_complete)
write.csv(df_gse_complete, file = "/results/before_nlp.csv", row.names = F)