# R 4.3.0
#BiocManager::install("AnnotationHub")
#library(AnnotationHub, quietly=T)
#BiocManager::install("ensembldb")
#library(ensembldb, quietly=T)
#BiocManager::install("multtest")
#library(multtest)
#BiocManager::install("glmGamPoi")
#library(glmGamPoi, quietly = T)
#library(RCurl)
#library(cowplot)
library(scales)
library(BiocManager, quietly=T)
library(Seurat, quietly=T)
library(dplyr, quietly=T)
library(tidyverse, quietly=T)
library(Matrix, quietly=T)
library(harmony, quietly=T)
library(DoubletFinder)
library(factoextra, quietly=T)
library(ggplot2)
sessionInfo()

# -----------------------------------------------------------------------------
################################################# read in counts --> seurat obj
# make seurat obj for each sample using unfiltered counts from cellranger
# list all folders in the current wd, which contains each samples' 3 musketeers
folders <- basename(list.dirs(path = ".", full.names = TRUE, recursive = FALSE))
# initialize a list to store Seurat objects
seurat_objects <- list()
# go thru each sample folder to make seurat obj and store in the list
for (file in folders){
  seurat_data <- Read10X(file)
  seurat_obj <- CreateSeuratObject(counts = seurat_data, project = file)
  #seurat_obj$sample <- file
  # add to list
  seurat_objects[[file]] <- seurat_obj
}
# merge seurat objects
first_object <- seurat_objects[[1]]
remaining_objects <- seurat_objects[-1]
merged_seurat <- merge(x = first_object, y = remaining_objects, add.cell.id = names(seurat_objects))
# concatenate count matrices of both samples
seurat_merged <- JoinLayers(merged_seurat)
head(seurat_merged@meta.data)
tail(seurat_merged@meta.data)

# -----------------------------------------------------------------------------
############################################################ add needed columns
# make genes/UMI for each cell column 
seurat_merged$log10GenesPerUMI <- log10(seurat_merged$nFeature_RNA) / 
  log10(seurat_merged$nCount_RNA)
# make mito ratio column
seurat_merged$mitoRatio <- PercentageFeatureSet(object = seurat_merged, pattern = "^MT-")
seurat_merged$mitoRatio <- seurat_merged@meta.data$mitoRatio / 100
# add some more columns to metadata separately 
metadata <-seurat_merged@meta.data
# rename some columns
metadata <- metadata %>%
  dplyr::rename(sample = orig.ident, umiPerCell = nCount_RNA, 
                genePerCell = nFeature_RNA)
# add metadata back to seurat obj
seurat_merged@meta.data <- metadata
head(seurat_merged@meta.data)
tail(seurat_merged@meta.data)

# -----------------------------------------------------------------------------
######################################################################## QC 
# visualize UMI/cell w histogram
metadata %>%
  ggplot(aes(color=sample, x=umiPerCell, fill=sample))+
  geom_density(alpha=0.2)+
  scale_x_log10()+
  theme_classic() + ylab('log10 cell density')+
  ylab('log10 cell density')+
  geom_vline(xintercept=750, linetype='dashed')
# visualize gene/cell w histogram
metadata %>% 
  ggplot(aes(color=sample, x=genePerCell, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() + ylab('log10 cell density')+
  scale_x_log10() + 
  geom_vline(xintercept = 250, linetype='dashed')
# visualize mito ratio w histogram
metadata %>%
  ggplot(aes(color=sample, x=mitoRatio, fill=sample))+
  geom_density(alpha=0.2)+
  scale_x_log10()+
  theme_classic() + ylab('log10 cell density')+
  geom_vline(xintercept = 0.1, linetype='dashed')
# visualize complexity w histogram
#metadata %>% 
#ggplot(aes(color=sample, x=log10GenesPerUMI, fill=sample))+
#geom_density(alpha=0.2)+
#theme_classic()+
#geom_vline(xintercept=0.8)

# filter with your desired params!
seurat_filtered <-subset(x=seurat_merged, subset=(umiPerCell>=500) &
                           (genePerCell>=250) & (mitoRatio<0.2))
# extract counts
counts <- GetAssayData(object = seurat_filtered, layer='counts')
# output logical matrix indicating per gene if it more than zero per cell
nonzero <- counts >0
# only keep genes expressed in 3 or more cells (rm cells with 0 counts)
# sums all true values and returns true if more than 3 true values per gene
keep_genes <-Matrix::rowSums(nonzero) >= 3
filtered_counts <-counts[keep_genes,]
seurat_filtered <-CreateSeuratObject(filtered_counts, 
                                     meta.data=seurat_filtered@meta.data)

# -----------------------------------------------------------------------------
##################### OPTIONAL #################### doublet removal w Scrubblet
# to be conducted with raw (unormalized UMI counts matrix) per sample w Python

# -----------------------------------------------------------------------------
############ normalize, FindVariableFeatures, ScaleData, run PCA, Harmony, UMAP
seurat_final <-NormalizeData(seurat_filtered)
# find variable features w defaults
seurat_final <-FindVariableFeatures(seurat_final, select.method='vst', 
                                    nfeatures =2000, verbose=F)
# scale the counts
seurat_final <-ScaleData(seurat_final)
# perform PCA
seurat_final <-RunPCA(seurat_final)
DimPlot(seurat_final, reduction='pca')
# SCTranform (if you wanted to regress something out)
#seurat_final <- SCTransform(seurat_final, vars.to.regress = c("mitoRatio"))
# do harmony on sample
seurat_final <- RunHarmony(seurat_final, group.by.vars = 'sample')

# -----------------------------------------------------------------------------
################# OPTIONAL #################### doublet removal w DoubletFinder 
##################### pK Identification (no ground-truth) 
# first, ID which PCs explain most variance with scree plot and elbow plot (altho usually PC1-10 used)
pca_variance <- seurat_final[["pca"]]@stdev^2
pca_variance_ratio <- pca_variance / sum(pca_variance)
ggplot(data.frame(PC = seq_along(pca_variance_ratio), Variance = pca_variance_ratio), 
       aes(x = PC, y = Variance)) +
  geom_point() +
  geom_line() +
  xlab("Principal Component") +
  ylab("Proportion of Variance Explained") +
  theme_minimal() 

sweep.res.list <- paramSweep(seurat_final, PCs = 1:10, sct = FALSE)
# can use ground truth data instead of mathematical if such wet lab info avail (check their github for instruc)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
# to deterimine optiimal pK value, which determines proportion of doublets
# bcmvn = mean-variance normalized bimodality coefficient
bcmvn <- find.pK(sweep.stats)
# extract optimal pK value that maximizes performance metric (BCmetric)
# meanBC = how well doublet scores are sep'd from singlet scores (lower the better)
# varBC = var in BC across cells (lower the more consistent)
# bcmatric = higher the better, more sep between singlets and doublets
print(optimal_pK <- bcmvn[which.max(bcmvn$BCmetric),])
optimal_pK_value <- as.numeric(optimal_pK$pK)
##################### homotypic Doublet Proportion Estimate 
# model homotypic doublets (doublets from same cluster) to adj for overest of doublets
annotations <- seurat_final@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
# Assuming 7.5% doublet formation rate - tailor for your dataset, can be 0.05-0.1
nExp_poi <- round(0.075*nrow(seurat_final@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

##################### run DoubletFinder with varying classification stringencies 
# pN (prior probability of doublet) typically ranges from 0.25-0.5, smaller the more stringent
seurat_final <- doubletFinder(seurat_final, PCs = 1:10, pN = 0.25, pK = optimal_pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# refine doublet detection from first run
pANN_colname <- grep("^pANN_", colnames(seurat_final@meta.data), value = TRUE)
seurat_final <- doubletFinder(seurat_final, PCs = 1:10, pN = 0.25, pK = optimal_pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_colname, sct = FALSE)
# filter out doublets
DF_colname <-grep('^DF.classifications_', colnames(seurat_final@meta.data), value=T)
is_singlet <- seurat_final@meta.data[[DF_colname]] == "Singlet"
seurat_final_sans_DFdoublets <- seurat_final[, is_singlet]
# check number of cells before and after filtering
n_before <- ncol(seurat_final)
n_after <- ncol(seurat_final_sans_DFdoublets)
print(paste("Number of cells before DF filter:", n_before))
print(paste("Number of cells after DF filter:", n_after))
seurat_final <- seurat_final_sans_DFdoublets

# -----------------------------------------------------------------------------
############################################ clustering, cluster ID, annotation
# make knn graph #diff
seurat_final <-FindNeighbors(object=seurat_final, 
                             #reduction = 'harmony', 
                             dims=1:30)
# find clusters
seurat_final <- FindClusters(object = seurat_final, resolution = 0.7)
#seurat_final@meta.data %>% 
# View()
# run UMAP
set.seed(123456)
seurat_final <-RunUMAP(seurat_final, 
                       #reduction='harmony',
                       dims=1:30)
# assign identity of clusters
Idents(object = seurat_final) <- "RNA_snn_res.0.7"
# plot umap
DimPlot(seurat_final, reduction='umap', label=T, label.size=6)

# -----------------------------------------------------------------------------
############################################################ SignacX annotation
library(SignacX)
# classify cells with SignacX prebuilt, already trained model
labels <-Signac(E= seurat_final)
celltypes = GenerateLabels(labels, E=seurat_final)
seurat_final$signacx_cell_states <- celltypes

################################################################################
##############################     la fin    ###################################
################################################################################