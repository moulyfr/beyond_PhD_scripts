library(Seurat, quietly=TRUE)
packageVersion('Seurat') # should be 5.0.1
packageVersion('SeuratObject') # should be 5.0.0

library(purrr, quietly=TRUE)
library(Matrix, quietly=TRUE)
packageVersion('Matrix') # should be >= 1.5-0

############ DEFINE INPUTS #############
# where all the rds are from 7b
rds_directory <- '/cloud-home/mfr/scrnaseq/CellphoneDB_analysis/processed_seuratobj_rds/'
# where all the matrix.mtx, meta.txt, genes.tsv, barcodes.tsv will go, within each dataset folder
output_directory <-  '/cloud-home/mfr/scrnaseq/CellphoneDB_analysis/cpdb_inputs/'
select_annotation <- 'signacx_cellstates'

############ CHECK VERSION OF ALL RDS, should be 5.0.1 #############
# List all RDS files in the directory
rds_files <- list.files(rds_directory, pattern = "\\.rds$", full.names = TRUE)
# Loop through each file and print its version
for (file_path in rds_files) {
  # Read the RDS file
  rds_object <- readRDS(file_path)
  # Print file path to know which RDS being inspected
  cat("File:", file_path, "\n")
  cat("Version:\n")
  # Print version
  print(rds_object@version)
  cat("\n")
}

############ DEFINE INPUTS - B #############
# indicate rds file splittings
select_patterns <- list('Skin_Cell_Atlas_cellbridg.rds' = c('disease','site'),
                        'AD_Skin_GSE147424_cellbridge.rds' = 'disease',
                        'AD_Skin_GSE153760_cellbridge.rds' = 'diagnosis',
                        'PSO_Skin_GSE173706_cellbridge.rds' = 'disease',
                        'PSO_Skin_GSE220116_cellbridge.rds' = c('disease', 'treatment'), 
                        'Human_Lung_Cell_Atlas_cellbridge.rds' = 'disease', 
                        'COPD_Lung_GSE136831_cellbridge.rds' = 'disease',
                        'COPD_Lung_GSE171541_cellbridge.rds' = 'disease', 
                        'ILD_Lung_GSE1122960_cellbridge.rds' = 'disease',
                        'ILD_Lung_GSE135893_cellbridge.rds' = 'disease', 
                        'UC_Colon_SCP259_cellbridge.rds' = 'disease', 
                        'UC_Colon_GSE116222_cellbridge.rds' = 'disease', 
                        'UC_Colon_GSE231993_cellbridge.rds' = 'disease',
                        'IgAN_Cell_Reports_Zheng_Kidney_cellbridge.rds' = 'disease',
                        'SLE_Phase2_scRNA_cellbridge.rds' = 'disease' 
                        )

########## DEFINE FUNCTION #############
create_cpdb_input <- function(infile, patterns=select_patterns,
                              annotation=select_annotation,
                              outdir=output_directory) {
  
  if (dir.exists(file.path(outdir, infile))==FALSE){
    dir.create(file.path(outdir, infile))
  }
  # what are we doing
  cat('processing', ' ', infile, '\n')
  
  # get rds
  in_rds <- readRDS(infile)
  
  # get metadata
  df_meta <- data.frame(Cell=colnames(in_rds[['RNA']]),
                        cell_type=in_rds@meta.data[,annotation],
                        row.names=colnames(in_rds[['RNA']]))
  # extract subset of metadata from seurat obj based on select_patterns
  def_select <- in_rds@meta.data[,patterns[[infile]]]
  
  # combine multiple metadata columns into 1 string for each cell to 
  # create unique identifier based on combined metadata attributes
  if (!is.null(ncol(def_select))) {
    def_select <- as.vector(apply(def_select, 1, 
                                  function(row) paste(row, collapse = "__")))
  }
  
  # go through unique metadata categories, format category name
  # create a directory for each category if it does not already exist
  for (category in unique(def_select)){
    out_category <- gsub(" ", "_", category)
    if (dir.exists(file.path(outdir, infile, out_category))==FALSE){
      dir.create(file.path(outdir, infile, out_category))
    }
    
    # get names of the cells beloinging to that category 
    selected <- def_select == category
    get_names <- colnames(in_rds[['RNA']])[selected]
    
    # write the matrix file; counts is raw, data is normalized
    writeMM(GetAssayData(in_rds[,selected], assay = "RNA", layer = "data"),
            file.path(outdir, infile, out_category, "matrix.mtx"))
    
    # write the features file; add genes column and move it to front, then save file
    write.table(rownames(in_rds[['RNA']]), 
                file.path(outdir, infile, out_category, "features.tsv"), 
                col.names=F, row.names=F, quote=F)
    
    # write the barcodes file
    write.table(get_names, 
                file.path(outdir, infile, out_category, "barcodes.tsv"), 
                col.names=F, row.names=F, quote=F)
    
    # check if it is getting the correct disease
    print(unique(in_rds@meta.data[get_names,]$disease))
    
    # now metadata
    write.table(df_meta[get_names,], file.path(outdir, infile, out_category, "meta.txt"),
                row.names = FALSE, sep='\t', quote=F)
    
    # what are we done
    cat('Done processing', ' ', infile, '\n')
  }
}

########## RUN FUNCTION #############
setwd(rds_directory)
files <- list.files(pattern = "\\.rds$")
# simpler way than to do for loop for file in files running the function on file
purrr::walk(files, create_cpdb_input)

