# install necessary packages
# remotes::install_github("stuart-lab/signac", ref="develop")
# install.packages("Matrix", type = "source")
# install.packages("irlba", type = "source")
# BiocManager::install("EnsDb.Hsapiens.v75")

library(Signac)                # for single-cell ATAC-seq analysis
library(Seurat)                # for single-cell RNA-seq analysis and integration
library(EnsDb.Hsapiens.v75)    # human genome annotations
library(tidyverse)             # for data manipulation and visualization

# -----------------------------------------------------------------------------
# 1) read in data
# -----------------------------------------------------------------------------

# read fragment file 
frag_file_path <- 'data/atac_v1_pbmc_10k_fragments.tsv.gz'
frag_file <- read.delim(frag_file_path, header = FALSE, nrows = 10)
head(frag_file)

# read counts matrix from 10X Genomics HDF5 file
counts_file_path <- 'data/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5'
counts_matrix <- Read10X_h5(counts_file_path)
counts_matrix[1:10, 1:10]

# read metadata associated with cells
metadata_file_path <- 'data/atac_v1_pbmc_10k_singlecell.csv'
metadata <- read.csv(file = metadata_file_path, header = TRUE, row.names = 1)
View(metadata)

# -----------------------------------------------------------------------------
# 2) create objects and annotate
# -----------------------------------------------------------------------------

# create Chromatin Assay object using the counts matrix and fragment file
chrom_assay <- CreateChromatinAssay(
  counts = counts_matrix,
  sep = c(":", "-"),                     # separator for chromosome and position in the fragment file
  fragments = frag_file_path,
  min.cells = 10,                        # minimum number of cells a feature must be present in
  min.features = 200                     # minimum number of features a cell must have
)
str(chrom_assay)

# create Seurat object using Chromatin Assay and metadata
pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  meta.data = metadata,
  assay = 'ATAC'
)
str(pbmc)

# extract and add gene annotations to Seurat object
pbmc@assays$ATAC@annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# convert to UCSC-style chromosome names (data mapped to hg19)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))

# add gene annotations to Seurat object
Annotation(pbmc) <- annotations
pbmc@assays$ATAC@annotation

# -----------------------------------------------------------------------------
# 3) QC - check
# -----------------------------------------------------------------------------

# nucleosome signal score for each cell
pbmc <- NucleosomeSignal(pbmc)

# Transcription Start Site (TSS) enrichment score for each cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# calculate blacklist ratio and percentage of reads in peaks
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

# view updated metadata with QC metrics
View(pbmc@meta.data)

# visualize QC metrics with density scatter plots
qc_plot1 <- DensityScatter(pbmc, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
qc_plot2 <- DensityScatter(pbmc, x = 'nucleosome_signal', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
qc_plot1 | qc_plot2

# create violin plots for various QC metrics
VlnPlot(object = pbmc, 
        features = c('nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'nucleosome_signal', 'blacklist_ratio', 'pct_reads_in_peaks'),
        pt.size = 0.1,
        ncol = 6)

# -----------------------------------------------------------------------------
# 4) QC - filter
# -----------------------------------------------------------------------------

# filter out low-quality cells based on QC metrics
pbmc <- subset(x = pbmc,
               subset = nCount_ATAC > 3000 &
                 nCount_ATAC < 30000 &
                 pct_reads_in_peaks > 15 & 
                 blacklist_ratio < 0.05 &
                 nucleosome_signal < 4 &
                 TSS.enrichment > 3)

# -----------------------------------------------------------------------------
# 5) normalization and analysis
# -----------------------------------------------------------------------------

# normalize data using TF-IDF (Term Frequency-Inverse Document Frequency)
pbmc <- RunTFIDF(pbmc) 

# ID top features (peaks) based on variability
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0') 

# perform linear dimensionality reduction using Singular Value Decomposition (SVD)
pbmc <- RunSVD(pbmc) 

# check depth correlation (optional)
DepthCor(pbmc)

# perform non-linear dimensionality reduction (UMAP) and clustering
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, algorithm = 3)

# plot the UMAP visualization of the clusters
DimPlot(object = pbmc, label = TRUE) + NoLegend()
