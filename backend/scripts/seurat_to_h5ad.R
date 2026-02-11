args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript seurat_to_h5ad.R <input.(rds|h5seurat)> <output.h5ad>")
}

input_path <- args[1]
output_path <- args[2]

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
})

if (grepl("\\.rds$", input_path, ignore.case = TRUE)) {
  obj <- readRDS(input_path)
  tmp_h5seurat <- tempfile(fileext = ".h5seurat")
  SaveH5Seurat(obj, filename = tmp_h5seurat, overwrite = TRUE)
  Convert(tmp_h5seurat, dest = "h5ad", overwrite = TRUE)
  # SeuratDisk writes base.h5ad (replaces .h5seurat), not base.h5seurat.h5ad
  converted <- sub("\\.h5seurat$", ".h5ad", tmp_h5seurat)
  if (!file.copy(converted, output_path, overwrite = TRUE)) {
    stop("file.copy failed: ", converted, " -> ", output_path)
  }
} else if (grepl("\\.h5seurat$", input_path, ignore.case = TRUE)) {
  Convert(input_path, dest = "h5ad", overwrite = TRUE)
  converted <- sub("\\.h5seurat$", ".h5ad", input_path)
  if (!file.copy(converted, output_path, overwrite = TRUE)) {
    stop("file.copy failed: ", converted, " -> ", output_path)
  }
} else {
  stop("Unsupported file extension. Use .rds or .h5seurat")
}

cat("Converted to:", output_path, "\n")

