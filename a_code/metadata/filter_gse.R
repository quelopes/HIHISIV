filter_gse <- function(gse_data) {
  filter <- gse_data[which(is.na(gse_data$n_samples) == FALSE), ] # 16325
  filter <- filter[which(is.na(filter$type) == FALSE), ] # 764
  filter <- filter[which(is.na(filter$series_accession) == FALSE), ] # 753
  filter <- filter[which(filter$n_samples > 3), ] # 400

  filter


  # unique(filter$type)

  # [1] "Expression profiling by high throughput sequencing"
  # [2] "Methylation profiling by array"
  # [3] "Expression profiling by array"
  # [4] "Non-coding RNA profiling by array"
  # [5] "Other"
  # [6] "Genome binding/occupancy profiling by high throughput sequencing; Expression profiling by high throughput sequencing"
  # [7] "Non-coding RNA profiling by high throughput sequencing"
  # [8] "Genome binding/occupancy profiling by high throughput sequencing"
  # [9] "Non-coding RNA profiling by high throughput sequencing; Expression profiling by high throughput sequencing"
  # [10] "Expression profiling by high throughput sequencing; Genome binding/occupancy profiling by high throughput sequencing"
  # [11] "Other; Expression profiling by high throughput sequencing"
  # [12] "Expression profiling by high throughput sequencing; Non-coding RNA profiling by high throughput sequencing"
  # [13] "Methylation profiling by genome tiling array"
  # [14] "Expression profiling by high throughput sequencing; Other"
  # [15] "Methylation profiling by high throughput sequencing"
  # [16] "Expression profiling by array; Non-coding RNA profiling by array"
  # [17] "Expression profiling by high throughput sequencing; Methylation profiling by high throughput sequencing; Genome binding/occupancy profiling by high throughput sequencing"
  # [18] "Expression profiling by RT-PCR"
  # [19] "Genome variation profiling by SNP array"
  # [20] "Genome variation profiling by genome tiling array"
  # [21] "Methylation profiling by high throughput sequencing; Expression profiling by high throughput sequencing"
  # [22] "SNP genotyping by SNP array"
  # [23] "Expression profiling by SNP array"
  # [24] "Genome variation profiling by array"
  # [25] "Expression profiling by high throughput sequencing; Expression profiling by array"
  # [26] "Genome variation profiling by SNP array; SNP genotyping by SNP array"
  # [27] "Expression profiling by array; Genome variation profiling by array"
}
