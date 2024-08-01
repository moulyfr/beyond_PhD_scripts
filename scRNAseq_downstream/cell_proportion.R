# R v 4.3.0
# brunner munzel for cell proportion analysis

packages <- c("Seurat","tidyverse", "brunnermunzel", 'writexl')
lapply(packages, library, character.only=TRUE, quietly = TRUE)  
`%nin%` = Negate(`%in%`) # negate in operator

################################################################################
###########################     THE FUNCTION  ##################################
################################################################################
cell_prop_func <- function(seurat_obj) {
  
  ##################    params for all datasets    #############################
  exclude <- c('Unclassified', "NonImmune")
  anno <- 'signacx_cellstates'
  min_subject_nonzero_total <- 4 # both groups together must have >= 4 subjects with non-zero abundance measurements in order for analysis to be conducted (per cell state)
  min_subject_total <- 2 # min number of subjects in EACH group to hold a comparison (per cell state) - 2
  min_cells_per_subject <- 200 # TOTAL cells minus excluded cell types, original was 500
  
  ##############################  get pairwise groups    #######################
  # remove the cell types we don't need (e.g. unclassified)
  seurat_obj_test <- seurat_obj[ , seurat_obj[[anno]][[anno]] %nin% exclude]
  
  # get the unique diseases to pair up
  pairwiise <- combn(unique(seurat_obj_test[[disease]][[disease]]), 2)
  list_of_pairwiise <- lapply(seq(1,ncol(pairwiise)),
                              function(colname) { as.vector(pairwiise[,colname])})
  names(list_of_pairwiise) <- lapply(list_of_pairwiise, 
                                     function(x) { paste(x, collapse='_')})
  
  ############################   get abundance info   ###########################
  # df of all cells in the dataset, and which sample and disease they're from
  cellstates <- data.frame(cell=factor(seurat_obj_test[[anno]][[anno]]),
                           donor=factor(seurat_obj_test[[sample]][[sample]]),
                           disease=factor(seurat_obj_test[[disease]][[disease]]))
  # for each sample, get counts of each cell type
  cellstates_grp_1 <- cellstates %>% group_by(cell, donor) %>% summarize(n1=n()) %>% complete(donor, fill = list(n1 = 0))
  # for each sample, get total cell counts
  cellstates_grp_2 <- cellstates %>% group_by(donor) %>% summarize(total=n())
  # for each sample, get cell type percentage 
  cellstates_grp_fin <- left_join(cellstates_grp_1, cellstates_grp_2, by = "donor") %>% mutate(perc=n1/total*100)
  unique(cellstates_grp_fin$total)
  # filter out subjects that have lower than the threshold of tot cells
  cellstates_grp_fin <- cellstates_grp_fin[cellstates_grp_fin$total >= min_cells_per_subject, ]
  
  ############################   get abundance stats    ############################
  # iterate over each element of list_of_pairwiise, and for each pair, assigns to grp1_grp2
  abundance_stats <- lapply(list_of_pairwiise,
                            function(grp1_grp2) {
                              
                              # extract unique donors from group 1
                              grp1_subj <- unique(cellstates[cellstates$disease==grp1_grp2[1], ]$donor)
                              # extract unique donors from group 2
                              grp2_subj <- unique(cellstates[cellstates$disease==grp1_grp2[2], ]$donor)
                              
                              # for each cell state...
                              alls <- lapply(sort(unique(cellstates_grp_fin$cell)), function(state) {
                                
                                # per group, extract perc of cells corresponding to that particular cell state
                                grp1 = cellstates_grp_fin[cellstates_grp_fin$cell==state & cellstates_grp_fin$donor %in% grp1_subj, ]$perc
                                grp2 = cellstates_grp_fin[cellstates_grp_fin$cell==state & cellstates_grp_fin$donor %in% grp2_subj, ]$perc
                                
                                # do the stat test if the total number of subjects (per group and together meet thresholds)
                                combined <- c(grp1, grp2)
                                if (length(combined[combined != 0]) >= min_subject_nonzero_total
                                    & length(grp1) >= min_subject_total & length(grp2) >= min_subject_total) {
                                  brunnermunzel.test(grp1, grp2)
                                }
                              })
                              # add names of cell states to the alls result
                              names(alls) <- sort(unique(cellstates_grp_fin$cell))
                              
                              ######
                              print('all results no fdr yet~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ')
                              print(alls <- alls[lapply(alls, is.null)==FALSE])
                              combined <- as.data.frame(do.call(rbind, alls))
                              
                              # do BH fdr correction
                              combined$p.adjust <- p.adjust(combined$p.value)
                              
                              # this clunky code is all to indicate which group is higher
                              higher_freq_group <- as.vector(unlist(combined$statistic))
                              higher_freq_group[higher_freq_group > 0] <- grp1_grp2[2]
                              higher_freq_group[higher_freq_group != grp1_grp2[2]] <- grp1_grp2[1]
                              combined$higher_freq_group <- higher_freq_group
                              # columns to keep
                              combined <- combined[,c("higher_freq_group", "statistic", "p.value", "p.adjust", "conf.int")]
                              combined
                              
                            })
  # make into table
  table <- bind_rows(abundance_stats, .id = "Comparison")
  table$cell <- rownames(table)
  rownames(table) <- NULL
  table$cell <- gsub("\\.\\.\\..*\\d$", "", table$cell)
  print(table)
  print(' ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~ DONE  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~')
  return(table)
}

################################################################################
#################################     inputs   #################################
################################################################################
partial_path <- '.../CellphoneDB_analysis/cellbridge_rds/'

# AD_Skin_GSE147424 
seurat_obj <- readRDS(paste0(partial_path, 'AD_Skin_GSE147424_cellbridge.rds'))
disease <- 'disease'
sample <- 'sample'
AD_Skin_GSE147424_result <- cell_prop_func(seurat_obj)
AD_Skin_GSE147424_result$dataset <- 'AD_Skin_GSE147424'

# AD_Skin_GSE153760 
seurat_obj <- readRDS(paste0(partial_path, 'AD_Skin_GSE153760_cellbridge.rds'))
disease <- 'diagnosis'
sample <- 'sample'
AD_Skin_GSE153760_result <- cell_prop_func(seurat_obj)
AD_Skin_GSE153760_result$dataset <- 'AD_Skin_GSE153760'

# PSO_Skin_GSE173706 
seurat_obj <- readRDS(paste0(partial_path, 'PSO_Skin_GSE173706_cellbridge.rds'))
disease <- 'disease'
sample <- 'sample'
PSO_Skin_GSE173706_result <- cell_prop_func(seurat_obj)
PSO_Skin_GSE173706_result$dataset <- 'PSO_Skin_GSE173706'

# PSO_Skin_GSE220116 
seurat_obj <- readRDS(paste0(partial_path, 'PSO_Skin_GSE220116_cellbridge.rds'))
seurat_obj@meta.data$disease2 <- paste(seurat_obj@meta.data$disease, seurat_obj@meta.data$treatment, sep = "_")
disease <- 'disease2'
sample <- 'sample'
PSO_Skin_GSE220116_result <- cell_prop_func(seurat_obj)
PSO_Skin_GSE220116_result$dataset <- 'PSO_Skin_GSE220116'

# Skin_Cell_Atlas 
seurat_obj <- readRDS(paste0(partial_path, 'Skin_Cell_Atlas_cellbridge.rds'))
seurat_obj@meta.data$disease2 <- paste(seurat_obj@meta.data$disease, seurat_obj@meta.data$site, sep = "_")
disease <- 'disease2'
sample <- 'sample'
Skin_Cell_Atlas_result <- cell_prop_func(seurat_obj)
Skin_Cell_Atlas_result$dataset <- 'Skin_Cell_Atlas'

# COPD_Lung_GSE136831 
seurat_obj <- readRDS(paste0(partial_path, 'COPD_Lung_GSE136831_cellbridge.rds'))
disease <- 'disease'
sample <- 'sample'
COPD_Lung_GSE136831_result <- cell_prop_func(seurat_obj)
COPD_Lung_GSE136831_result$dataset <- 'COPD_Lung_GSE136831'

# COPD_Lung_GSE171541 
seurat_obj <- readRDS(paste0(partial_path, 'COPD_Lung_GSE171541_cellbridge.rds'))
disease <- 'disease'
sample <- 'sample'
COPD_Lung_GSE171541_result <- cell_prop_func(seurat_obj)
COPD_Lung_GSE171541_result$dataset <- 'COPD_Lung_GSE171541'

# ILD_Lung_GSE1122960 
seurat_obj <- readRDS(paste0(partial_path, 'ILD_Lung_GSE1122960_cellbridge.rds'))
disease <- 'disease'
sample <- 'sample'

# ILD_Lung_GSE135893 
seurat_obj <- readRDS(paste0(partial_path, 'ILD_Lung_GSE135893_cellbridge.rds'))
disease <- 'disease'
sample <- 'sample'

# UC_Colon_SCP259 
seurat_obj <- readRDS(paste0(partial_path, 'UC_Colon_SCP259_cellbridge.rds'))
disease <- 'disease'
sample <- 'subject'
UC_Colon_SCP259_result <- cell_prop_func(seurat_obj)
UC_Colon_SCP259_result$dataset <- 'UC_Colon_SCP259'

# UC_Colon_GSE231993 
seurat_obj <- readRDS(paste0(partial_path, 'UC_Colon_GSE231993_cellbridge.rds'))
disease <- 'disease'
sample <- 'sample'
UC_Colon_GSE231993_result <- cell_prop_func(seurat_obj)
UC_Colon_GSE231993_result$dataset <- 'UC_Colon_GSE231993'

# IgAN_Cell_Reports_Zheng_Kidney
seurat_obj <- readRDS(paste0(partial_path, 'IgAN_Cell_Reports_Zheng_Kidney_cellbridge.rds'))
disease <- 'disease'
sample <- 'patient' 
IgAN_Cell_Reports_Zheng_Kidney_result <- cell_prop_func(seurat_obj)

# SLE_Phase2_scRNA 
seurat_obj <- readRDS(paste0(partial_path, 'SLE_Phase2_scRNA_cellbridge.rds'))
disease <- 'disease'
sample <- 'sample'
SLE_Phase2_scRNA_result <- cell_prop_func(seurat_obj)
SLE_Phase2_scRNA_result$dataset <- 'SLE_Phase2_scRNA'

################################################################################
#################################     export   #################################
################################################################################
export <- bind_rows(AD_Skin_GSE147424_result, AD_Skin_GSE153760_result, PSO_Skin_GSE173706_result, PSO_Skin_GSE220116_result, Skin_Cell_Atlas_result, COPD_Lung_GSE136831_result, UC_Colon_SCP259_result, UC_Colon_GSE231993_result, SLE_Phase2_scRNA_result)

# take significance from the combined frame
export <- export[export$p.adjust < 0.05, ]

# fix some of the columns so they show up properly
export$statistic <- gsub(".*\\n", "", export$statistic)
export$p.value <- sapply(export$p.value, function(x) as.character(x))
CI_lower <- numeric(length(export$conf.int))
CI_upper <- numeric(length(export$conf.int))
for (i in seq_along(export$conf.int)) {
  CI_lower[i] <- export$conf.int[[i]]["lower"]
  CI_upper[i] <- export$conf.int[[i]]["upper"]
}
export$CI_lower <- CI_lower
export$CI_upper <- CI_upper
conf.level <- sapply(export$conf.int, attr, "conf.level")
export$conf.int <- conf.level

export_path <- ".../Analyses/Analysis_cell_props/analysis_cell_proportion.xlsx"
write_xlsx(export, export_path)

##### # # # #  #  #  #  #   #   #   #   # done 