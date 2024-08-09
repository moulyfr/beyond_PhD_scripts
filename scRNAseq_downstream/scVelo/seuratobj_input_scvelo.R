# load necessary libraries
library(Seurat)
library(Matrix)

# load your Seurat object
seurat_object <- readRDS("ILD_Lung_GSE135893.rds")

# extract raw count data
raw_counts <- GetAssayData(seurat_object, assay = "RNA", slot = "counts")
# extract  normalized data
normalized_counts <- GetAssayData(seurat_object, assay = "RNA", layer = "data")
# Extract metadata
metadata <- seurat_object@meta.data

# Save them to csv files
write.csv(as.matrix(raw_counts), "raw_counts.csv")
write.csv(as.matrix(normalized_counts), "normalized_counts.csv")
write.csv(metadata, "metadata.csv")
