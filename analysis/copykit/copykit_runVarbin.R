# Script to run CopyKit runVarbin() function on BAM 
# files and save the output as a .rds file

dir <- snakemake@params[['dir']]
resolution <- snakemake@params[['resolution']]
threads <- snakemake@threads
output <- snakemake@output[['rds']]

# Log the output
# https://github.com/kelly-sovacool/snakemake-Rscript-log-mwe/blob/master/script.R
log_smk <- function() {
  if (exists("snakemake") & length(snakemake@log) != 0) {
    log <- file(snakemake@log[1][[1]], open = "wt")
    sink(log, append = TRUE)
    sink(log, append = TRUE, type = "message")
  }
}

log_smk()

# Load required libraries
library(copykit)

# Parallel processing 
library(BiocParallel)
#register(MulticoreParam(progressbar = T, workers = threads), default = TRUE)
register(SerialParam(), default = TRUE)
BiocParallel::bpparam()
# https://navinlabcode.github.io/CopyKit-UserGuide/pre-processing.html#runvarbin-modules
scells <- runVarbin(
    dir,
    genome = "hg38",
    resolution = resolution, # Default is 220kb
    remove_Y = TRUE,
    is_paired_end = TRUE,
)

# QC metrics
scells <- runMetrics(scells)
scells <- findAneuploidCells(scells)
scells <- findOutliers(scells)

saveRDS(scells, file = output)

sessionInfo()