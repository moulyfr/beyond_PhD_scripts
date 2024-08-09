# R4.3.0
#devtools::install_github("jinworks/CellChat")
#devtools::install_github('immunogenomics/presto')
library(CellChat, quietly = TRUE)
library(patchwork, quietly = TRUE)
library(presto, quietly = TRUE)
options(stringsAsFactors = FALSE)

library(Seurat, quietly=TRUE)
packageVersion('Seurat') # should be 5.0.1
packageVersion('SeuratObject') # should be 5.0.0

library(purrr, quietly=TRUE)
library(Matrix, quietly=TRUE)
packageVersion('Matrix') # should be >= 1.5-0

rds_directory = '/cloud-data//mfr/scrnaseq/processed_seurat_rds/'
rds_file = 'ILD_Lung_GSE135893_cellbridge.rds'

seurat_obj <- readRDS(file.path(rds_directory, rds_file))

# rename the sample column to 'samples' so cellchat func works
seurat_obj@meta.data$samples <- seurat_obj@meta.data$sample

# subset the seurat_obj for the disease group to output 
seurat_obj <- seurat_obj[, seurat_obj@meta.data$disease == 'Control']

################################################################## run cellchat
cellChat <- createCellChat(object = seurat_obj, group.by = 'signacx_cellstates', assay = "RNA")
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB # use all CellChatDB for cell-cell communication analysis
cellChat@DB <- CellChatDB.use
future::plan("multisession", workers = 4) # do parallel
cellchat <- subsetData(cellChat) # the subset func subsets genes for just those known to be involved in signaling 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# *** parameter 1 that is subjective
#cellchat <- computeCommunProb(cellchat, type = "triMean") # type can also be truncatedMean, but then a value would need to be assigned to trim
# triMean produces fewer interactions that are stronger 
# population.size = TRUE can also be set, if you want to consider the effect of cell proportion in each cell group (let's not do that, which would make things more dissimilar to cpdb)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.05)
# *** parameter 2 that is subjective; min.cells
cellchat <- filterCommunication(cellchat, min.cells = 2)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

######################################################################## graphs 
ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

# Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level

netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair

# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# (1) show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)

# (2) show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c('TNF'), remove.isolate = FALSE)

#################################  la fin   ####################################
