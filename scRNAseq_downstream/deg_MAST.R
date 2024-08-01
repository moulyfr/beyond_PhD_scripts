# R v 4.3.0
# MAST: Model-based Analysis of Single-cell Transcriptomics to get DEG

#library(BiocManager)
#BiocManager::install("GGally", force=TRUE)

packages <- c('ggplot2','reshape2','data.table','knitr','stringr','rsvd','Biobase','MAST',"GGally","GSEABase",'limma',"SingleCellExperiment", "Seurat")
lapply(packages, library, character.only=TRUE, quietly = TRUE)  

################################################################################
###########################     THE FUNCTION  ##################################
################################################################################
DEG_func <- function(seurat_obj) {
  
  ##################    params for all datasets    ###############################
  disease <- 'disease'
  basegrp <- 'diseaseDiseased'
  freq_expressed <- 0.25 # how many cells should be expressing a particular gene, in the entire dataset per cell type (not by sample) - aka min.pct
  FCTHRESHOLD <- 0.25  # fold change of control (other norm is 0.5)
  min_cells <- 10 # for each sample, the min number of cells of a particular cell type
  min_samples <- 2 # min number of samples that passes the min_cells cutoff in each group
  
  # remove unnecessary cell states to be examined
  to_test <- sort(unique(seurat_obj$signacx_cellstates[!(seurat_obj$signacx_cellstates %in% c("Unclassified", "NonImmune"))]))
  #to_test <- c('Mon.Classical')
  
  #loop through each cell state
  get_degs <- lapply(to_test, function(cellst){
    cat('Now processing', cellst)
    
    # subset by cell state
    seurat_obj_subset <- subset(seurat_obj, signacx_cellstates == cellst)
    
    # subset by cell state for each sample (per disease group)
    samples_by_grp <- lapply(sort(unique(seurat_obj_subset@meta.data[[disease]])),
                             function(cellst) {
                               table(seurat_obj_subset@meta.data[[sample]][seurat_obj_subset@meta.data[[disease]] == cellst])
                             })
    # filter the samples by group based on the minimum number of cells needed per cell type
    samples_by_grp_filtered <- lapply(samples_by_grp,
                                      function(x){
                                        x[x >= min_cells ]
                                      })
    # determine min number of samples among filtered groups
    min_b_samples <- min(unlist(lapply(samples_by_grp_filtered, function(x) { length(x) })))
    
    # if min number of samples is less than the required minimum, then cancel
    if (min_b_samples < min_samples){
      return(NULL)
    }
    
    # count number of cells for each sample
    samples <- table(seurat_obj_subset[[sample]])
    # remove samples where they have less than the min number of cells
    keep_samples <- names(samples[samples > min_cells])
    seurat_obj_subset <- seurat_obj_subset[, seurat_obj_subset[[sample]][[sample]] %in% keep_samples]
    
    # Convert to SingleCellExperiment bc that's what MAST needs
    sce_obj <- as.SingleCellExperiment(seurat_obj_subset)
    sca <- SceToSingleCellAssay(sce_obj)
    # keep genes that are atleast expressed in whatever freq_expressed amount per cell state
    expressed_genes <- freq(sca) > freq_expressed
    sca <- sca[expressed_genes,]
    
    # Fit a hurdle (zero-inflated negative binomial) model 
    zlm_result <- zlm(formula, sca)
    # assess significance of the diagnosis factor
    summary_result <- summary(zlm_result, doLRT=basegrp)
    summaryDt <- summary_result$datatable
    
    # merge primerid, unadjusted pval(pr>chisq), coef, cis 
    fcHurdle <- merge(summaryDt[contrast==basegrp & component=='H',.(primerid, `Pr(>Chisq)`)], #p values
                      summaryDt[contrast==basegrp & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
    
    # compute FDR using Benjamini-Hochberg method (fdr) for the pvals  (Pr(>Chisq)) 
    fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    # filter results where p val <0.05 and abs(coef) is greater than preset for FCTHRESHOLD
    fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
    setorder(fcHurdleSig, fdr)
    return(fcHurdleSig)
    cat('Done processing', cellst)
  })
  
  # add name to each dataframe in the list output get_degs
  names(get_degs) <- to_test
  # make table
  table <- rbindlist(get_degs, idcol = "Case.CellStates")
  # rename some columns for compatibility with later meta analyses
  colnames(table)[colnames(table) == 'primerid'] <- 'gene'
  colnames(table)[colnames(table) == 'coef'] <- 'logFC'
  colnames(table)[colnames(table) == 'fdr'] <- 'FDR'
  table$Direction <- ifelse(table$logFC > 0, "UP", ifelse(table$logFC < 0, "DOWN", "UNCHANGED"))
  return(table)
}

################################################################################
###########################       INPUTS      ##################################
################################################################################
partial_path <- '.../CellphoneDB_analysis/cellbridge_rds/'
partial_export_path <-"/cloud-data/its-cmo-darwin-magellan-workspaces-folders/WS_PMCB/RAHMAN.Mouly/Analyses/Analysis_DEG_MAST/"

################################################################################
# Skin_Cell_Atlas - Eczema_non_lesion
# 1) rds
Skin_Cell_Atlas <- readRDS(paste0(partial_path, 'Skin_Cell_Atlas_cellbridge.rds'))
seurat_obj <- Skin_Cell_Atlas
seurat_obj@meta.data$disease <- paste(seurat_obj@meta.data$disease, seurat_obj@meta.data$site, sep = "_")
seurat_obj <- seurat_obj[, seurat_obj@meta.data$disease %in% c('Healthy_non_lesion', 'Eczema_non_lesion')]

# 2) renaming
#colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == 'diagnosis'] <- 'disease' 
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'Healthy_non_lesion'] <- 'Control'
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'Eczema_non_lesion'] <- 'Diseased'

# 3) individual params
sample <- 'sample'
formula <- ~ disease + nFeature_RNA + sample

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
write.csv(deg, file = paste0(partial_export_path, "Skin_Cell_Atlas_Eczema_nonlesion_deg.csv"), row.names = FALSE)

################################################################################
# Skin_Cell_Atlas - Eczema_lesion
# 1) rds
#seurat_obj <- readRDS(paste0(partial_path, 'Skin_Cell_Atlas_cellbridge.rds'))
seurat_obj <- Skin_Cell_Atlas
seurat_obj@meta.data$disease <- paste(seurat_obj@meta.data$disease, seurat_obj@meta.data$site, sep = "_")
seurat_obj <- seurat_obj[, seurat_obj@meta.data$disease %in% c('Healthy_non_lesion', 'Eczema_lesion')]

# 2) renaming
#colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == 'diagnosis'] <- 'disease' 
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'Healthy_non_lesion'] <- 'Control'
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'Eczema_lesion'] <- 'Diseased'

# 3) individual params
sample <- 'sample'
formula <- ~ disease + nFeature_RNA + sample

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
deg$dataset <- 'Skin_Cell_Atlas_Eczema_lesion'
write.csv(deg, file = paste0(partial_export_path, "Skin_Cell_Atlas_Eczema_lesion_deg.csv"), row.names = FALSE)

################################################################################
# Skin_Cell_Atlas - PSO_lesion
# 1) rds
#seurat_obj <- readRDS(paste0(partial_path, 'Skin_Cell_Atlas_cellbridge.rds'))
seurat_obj <- Skin_Cell_Atlas
seurat_obj@meta.data$disease <- paste(seurat_obj@meta.data$disease, seurat_obj@meta.data$site, sep = "_")
seurat_obj <- seurat_obj[, seurat_obj@meta.data$disease %in% c('Healthy_non_lesion', 'Psoriasis_lesion')]

# 2) renaming
#colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == 'diagnosis'] <- 'disease' 
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'Healthy_non_lesion'] <- 'Control'
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'Psoriasis_lesion'] <- 'Diseased'

# 3) individual params
sample <- 'sample'
formula <- ~ disease + nFeature_RNA + sample

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
write.csv(deg, file = paste0(partial_export_path, "Skin_Cell_Atlas_PSO_lesion_deg.csv"), row.names = FALSE)

################################################################################
# Skin_Cell_Atlas - PSO_nonlesion
# 1) rds
#seurat_obj <- readRDS(paste0(partial_path, 'Skin_Cell_Atlas_cellbridge.rds'))
seurat_obj <- Skin_Cell_Atlas
seurat_obj@meta.data$disease <- paste(seurat_obj@meta.data$disease, seurat_obj@meta.data$site, sep = "_")
seurat_obj <- seurat_obj[, seurat_obj@meta.data$disease %in% c('Healthy_non_lesion', 'Psoriasis_non_lesion')]

# 2) renaming
#colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == 'diagnosis'] <- 'disease' 
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'Healthy_non_lesion'] <- 'Control'
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'Psoriasis_non_lesion'] <- 'Diseased'

# 3) individual params
sample <- 'sample'
formula <- ~ disease + nFeature_RNA + sample

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
deg$dataset <- 'Skin_Cell_Atlas_PSO_nonlesion'
write.csv(deg, file = paste0(partial_export_path, "Skin_Cell_Atlas_PSO_nonlesion_deg.csv"), row.names = FALSE)

################################################################################
# AD_Skin_GSE147424 - nonlesional
# 1) rds
AD_Skin_GSE147424 <- readRDS(paste0(partial_path, 'AD_Skin_GSE147424_cellbridge.rds'))
seurat_obj <- AD_Skin_GSE147424
seurat_obj <- seurat_obj[, seurat_obj@meta.data$disease %in% c('Healthy', 'Non Lesional')]

# 2) renaming
#colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == 'diagnosis'] <- 'disease' 
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'Healthy'] <- 'Control'
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'Non Lesional'] <- 'Diseased'

# 3) individual params
sample <- 'sample'
formula <- ~ disease + nFeature_RNA + sample

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
write.csv(deg, file = paste0(partial_export_path, "AD_Skin_GSE147424_nonlesion_deg.csv"), row.names = FALSE)

################################################################################
# AD_Skin_GSE147424 - lesional
# 1) rds
#seurat_obj <- readRDS(paste0(partial_path, 'AD_Skin_GSE147424_cellbridge.rds'))
seurat_obj <- AD_Skin_GSE147424
seurat_obj <- seurat_obj[, seurat_obj@meta.data$disease %in% c('Healthy', 'Lesional')]

# 2) renaming
#colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == 'diagnosis'] <- 'disease' 
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'Healthy'] <- 'Control'
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'Lesional'] <- 'Diseased'

# 3) individual params
sample <- 'sample'
formula <- ~ disease + nFeature_RNA + sample

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
write.csv(deg, file = paste0(partial_export_path, "AD_Skin_GSE147424_lesion_deg.csv"), row.names = FALSE)

################################################################################
# AD_Skin_GSE153760 
# 1) rds
seurat_obj <- readRDS(paste0(partial_path, 'AD_Skin_GSE153760_cellbridge.rds'))

# 2) renaming
colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == 'diagnosis'] <- 'disease' 
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'HC'] <- 'Control'
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'AD'] <- 'Diseased'

# 3) individual params
sample <- 'sample'
formula <- ~ disease + nFeature_RNA + sample

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
write.csv(deg, file = paste0(partial_export_path, "AD_Skin_GSE153760_deg.csv"), row.names = FALSE)

################################################################################
# PSO_Skin_GSE173706 - nonlesional
# 1) rds
PSO_Skin_GSE173706 <- readRDS(paste0(partial_path, 'PSO_Skin_GSE173706_cellbridge.rds'))
seurat_obj <- PSO_Skin_GSE173706
seurat_obj <- seurat_obj[, seurat_obj@meta.data$disease %in% c('Control', 'PSO_non_lesion')]

# 2) renaming
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'PSO_non_lesion'] <- 'Diseased'

# 3) individual params
sample <- 'sample'
formula <- ~ disease + nFeature_RNA + sample

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
write.csv(deg, file = paste0(partial_export_path, "PSO_Skin_GSE173706_nonlesion_deg.csv"), row.names = FALSE)

################################################################################
# PSO_Skin_GSE173706 - lesional
# 1) rds
#seurat_obj <- readRDS(paste0(partial_path, 'PSO_Skin_GSE173706_cellbridge.rds'))
seurat_obj <- PSO_Skin_GSE173706
seurat_obj <- seurat_obj[, seurat_obj@meta.data$disease %in% c('Control', 'PSO_lesion')]

# 2) renaming
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'PSO_lesion'] <- 'Diseased'

# 3) individual params
sample <- 'sample'
formula <- ~ disease + nFeature_RNA + sample

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
write.csv(deg, file = paste0(partial_export_path, "PSO_Skin_GSE173706_lesion_deg.csv"), row.names = FALSE)

################################################################################
# PSO_Skin_GSE220116 - pso_pre_tx
# 1) rds
seurat_obj <- readRDS(paste0(partial_path, 'PSO_Skin_GSE220116_cellbridge.rds'))
seurat_obj@meta.data$disease <- paste(seurat_obj@meta.data$disease, seurat_obj@meta.data$treatment, sep = "_")
seurat_obj <- seurat_obj[, seurat_obj@meta.data$disease %in% c('Control_pre_tx', 'PSO_pre_tx')]

# 2) renaming
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'Control_pre_tx'] <- 'Control'
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'PSO_pre_tx'] <- 'Diseased'

# 3) individual params
sample <- 'sample'
formula <- ~ disease + nFeature_RNA + sample

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
write.csv(deg, file = paste0(partial_export_path, "PSO_Skin_GSE220116_deg.csv"), row.names = FALSE)

################################################################################
# COPD_Lung_GSE136831 - copd
# 1) rds
COPD_Lung_GSE136831 <- readRDS(paste0(partial_path, 'COPD_Lung_GSE136831_cellbridge.rds'))
seurat_obj <- COPD_Lung_GSE136831 
seurat_obj <- seurat_obj[, seurat_obj@meta.data$disease %in% c('Control', 'COPD')]

# 2) renaming
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'COPD'] <- 'Diseased'

# 3) individual params
sample <- 'sample'
formula <- ~ disease + nFeature_RNA + sample

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
deg$dataset <- 'COPD_Lung_GSE136831_copd' 
write.csv(deg, file = paste0(partial_export_path, "COPD_Lung_GSE136831_copd_deg.csv"), row.names = FALSE)

################################################################################
# COPD_Lung_GSE136831 - ipf
# 1) rds
#seurat_obj <- readRDS(paste0(partial_path, 'COPD_Lung_GSE136831_cellbridge.rds'))
seurat_obj <- COPD_Lung_GSE136831 
seurat_obj <- seurat_obj[, seurat_obj@meta.data$disease %in% c('Control', 'IPF')]

# 2) renaming
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'IPF'] <- 'Diseased'

# 3) individual params
sample <- 'sample'
formula <- ~ disease + nFeature_RNA + sample

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
deg$dataset <- 'COPD_Lung_GSE136831_ipf' 
write.csv(deg, file = paste0(partial_export_path, "COPD_Lung_GSE136831_IPF_deg.csv"), row.names = FALSE)

################################################################################
# COPD_Lung_GSE171541 
# 1) rds
seurat_obj <- readRDS(paste0(partial_path, 'COPD_Lung_GSE171541_cellbridge.rds'))

# 2) renaming
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'control'] <- 'Control'
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'copd'] <- 'Diseased'

# 3) individual params
sample <- 'sample'
formula <- ~ disease + nFeature_RNA + sample

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
deg$dataset <- 'COPD_Lung_GSE171541' 
write.csv(deg, file = paste0(partial_export_path, "COPD_Lung_GSE171541_deg.csv"), row.names = FALSE)

################################################################################
# # ILD_Lung_GSE1122960 
# 1) rds
seurat_obj <- readRDS(paste0(partial_path, 'ILD_Lung_GSE1122960_cellbridge.rds'))
seurat_obj <- seurat_obj[, seurat_obj@meta.data$disease %in% c('Control', 'IPF')]

# 2) renaming
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'IPF'] <- 'Diseased'

# 3) individual params
sample <- 'sample'
formula <- ~ disease + nFeature_RNA + sample

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
deg$dataset <- 'ILD_Lung_GSE1122960' 
write.csv(deg, file = paste0(partial_export_path, "ILD_Lung_GSE1122960_deg.csv"), row.names = FALSE)

################################################################################
# ILD_Lung_GSE135893 
# 1) rds
seurat_obj <- readRDS(paste0(partial_path, 'ILD_Lung_GSE135893_cellbridge.rds'))
seurat_obj <- seurat_obj[, seurat_obj@meta.data$disease %in% c('Control', 'IPF')]

# 2) renaming
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'IPF'] <- 'Diseased'

# 3) individual params
sample <- 'sample'
formula <- ~ disease + nFeature_RNA + sample

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
deg$dataset <- 'ILD_Lung_GSE135893' 
write.csv(deg, file = paste0(partial_export_path, "ILD_Lung_GSE135893_deg.csv"), row.names = FALSE)

################################################################################
# UC_Colon_SCP259 - noninflamed
# 1) rds
UC_Colon_SCP259 <- readRDS(paste0(partial_path, 'UC_Colon_SCP259_cellbridge.rds'))
seurat_obj <- UC_Colon_SCP259
seurat_obj <- seurat_obj[, seurat_obj@meta.data$disease %in% c('Healthy', 'UC_Non_Inflamed')]

# 2) renaming
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'Healthy'] <- 'Control'
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'UC_Non_Inflamed'] <- 'Diseased'

# 3) individual params
sample <- 'subject'
formula <- ~ disease + nFeature_RNA + sample

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
deg$dataset <- 'UC_Colon_SCP259_noninflamed' 
write.csv(deg, file = paste0(partial_export_path, "UC_Colon_SCP259_noninflamed_deg.csv"), row.names = FALSE)

################################################################################
# UC_Colon_SCP259 - inflamed
# 1) rds
#seurat_obj <- readRDS(paste0(partial_path, 'UC_Colon_SCP259_cellbridge.rds'))
seurat_obj <- UC_Colon_SCP259
seurat_obj <- seurat_obj[, seurat_obj@meta.data$disease %in% c('Healthy', 'UC_Inflamed')]

# 2) renaming
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'Healthy'] <- 'Control'
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'UC_Inflamed'] <- 'Diseased'

# 3) individual params
sample <- 'subject'
formula <- ~ disease + nFeature_RNA + sample

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
deg$dataset <- 'UC_Colon_SCP259_inflamed' 
write.csv(deg, file = paste0(partial_export_path, "UC_Colon_SCP259_inflamed_deg.csv"), row.names = FALSE)

################################################################################
# UC_Colon_GSE116222 -noninflamed
# 1) rds
UC_Colon_GSE116222 <- readRDS(paste0(partial_path, 'UC_Colon_GSE116222_cellbridge.rds'))
seurat_obj <- UC_Colon_GSE116222
seurat_obj <- seurat_obj[, seurat_obj@meta.data$disease %in% c('Healthy', 'UC_Non_Inflamed')]

# 2) renaming
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'Healthy'] <- 'Control'
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'UC_Non_Inflamed'] <- 'Diseased'

# 3) individual params
#sample <- 'sample'
formula <- ~ disease + nFeature_RNA 

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
deg$dataset <- 'UC_Colon_GSE116222_noninflamed' 
write.csv(deg, file = paste0(partial_export_path, "UC_Colon_GSE116222_noninflamed_deg.csv"), row.names = FALSE)

################################################################################
# UC_Colon_GSE116222 -inflamed
# 1) rds
#seurat_obj <- readRDS(paste0(partial_path, 'UC_Colon_GSE116222_cellbridge.rds'))
seurat_obj <- UC_Colon_GSE116222
seurat_obj <- seurat_obj[, seurat_obj@meta.data$disease %in% c('Healthy', 'UC_Inflamed')]

# 2) renaming
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'Healthy'] <- 'Control'
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'UC_Inflamed'] <- 'Diseased'

# 3) individual params
sample <- 'sample'
formula <- ~ disease + nFeature_RNA 

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
deg$dataset <- 'UC_Colon_GSE116222_inflamed' 
write.csv(deg, file = paste0(partial_export_path, "UC_Colon_GSE116222_inflamed_deg.csv"), row.names = FALSE)

################################################################################
# # UC_Colon_GSE231993 -noninflamed
# 1) rds
UC_Colon_GSE231993 <- readRDS(paste0(partial_path, 'UC_Colon_GSE231993_cellbridge.rds'))
seurat_obj <- UC_Colon_GSE231993
seurat_obj <- seurat_obj[, seurat_obj@meta.data$disease %in% c('Control', 'UC_uninflamed')]

# 2) renaming
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'UC_uninflamed'] <- 'Diseased'

# 3) individual params
sample <- 'sample'
formula <- ~ disease + nFeature_RNA + sample

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
deg$dataset <- 'UC_Colon_GSE231993_noninflamed' 
write.csv(deg, file = paste0(partial_export_path, "UC_Colon_GSE231993_noninflamed_deg.csv"), row.names = FALSE)

################################################################################
# # UC_Colon_GSE231993 -inflamed
# 1) rds
#seurat_obj <- readRDS(paste0(partial_path, 'UC_Colon_GSE231993_cellbridge.rds'))
seurat_obj <- UC_Colon_GSE231993
seurat_obj <- seurat_obj[, seurat_obj@meta.data$disease %in% c('Control', 'UC_inflamed')]

# 2) renaming
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'UC_inflamed'] <- 'Diseased'

# 3) individual params
sample <- 'sample'
formula <- ~ disease + nFeature_RNA + sample

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
deg$dataset <- 'UC_Colon_GSE231993_inflamed' 
write.csv(deg, file = paste0(partial_export_path, "UC_Colon_GSE231993_inflamed_deg.csv"), row.names = FALSE)

################################################################################
# IgAN_Cell_Reports_Zheng_Kidney
# 1) rds
seurat_obj <- readRDS(paste0(partial_path, 'IgAN_Cell_Reports_Zheng_Kidney_cellbridge.rds'))

# 2) renaming
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'normal control'] <- 'Control'
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'IgAN'] <- 'Diseased'

# 3) individual params
sample <- 'sample'
formula <- ~ disease + nFeature_RNA + sample

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
deg$dataset <- 'IgAN_Cell_Reports_Zheng_Kidney' 
write.csv(deg, file = paste0(partial_export_path, "IgAN_Cell_Reports_Zheng_Kidney_deg.csv"), row.names = FALSE)

################################################################################
# SLE_Phase2_scRNA 
# 1) rds
seurat_obj <- readRDS(paste0(partial_path, 'SLE_Phase2_scRNA_cellbridge.rds'))

# 2) renaming
seurat_obj@meta.data$disease[seurat_obj@meta.data$disease == 'SLE'] <-'Diseased'

# 3) individual params
sample <- 'sample'
formula <- ~ disease + nFeature_RNA + sample

# 4) run
deg <- DEG_func(seurat_obj)

# 5) export
deg$dataset <- 'SLE_Phase2_scRNA' 
write.csv(deg, file = paste0(partial_export_path, "SLE_Phase2_scRNA_deg.csv"), row.names = FALSE)

############################ la fin ################################
