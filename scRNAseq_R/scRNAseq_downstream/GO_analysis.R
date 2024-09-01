library(clusterProfiler, quietly = TRUE)
library(org.Hs.eg.db, quietly = TRUE)
library(BiocManager, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(reshape2, quietly= TRUE)
library(AnnotationDbi, quietly = TRUE)

# get list of RDS files in current directory
file_names <- list.files(pattern = "\\.rds$")

# create empty list to store the loaded data
data_list <- list()

# iterate over the file names and load each file
for (file_name in file_names) {
  # Load the RDS file
  data <- readRDS(file_name)
  
  # add the loaded data to the list
  data_list[[file_name]] <- data
}
##################################################################
####### Go through each disease dataset for GO analysis ##########
##################################################################

for (i in seq_along(data_list)) {
  # access the second element [[2]] within each RDS file
  data <- data_list[[i]][[2]]
  
  # format the data into proper data frames
  data <- lapply(names(data), function(name) { # get names of each list
    data_frame <- data[[name]]
    data_frame$cell <- name  # add "cell" column and fill it with the list name
    data_frame$gene <- rownames(data_frame)  # create a new "gene" column from row names
    rownames(data_frame) <- NULL  # remove row names
    data_frame
  })
  
  # bind each cell type's dataframe into one dataframe
  dataset <- bind_rows(data)
  
  # subset data for rows where AdjustedPValue is <= 0.05 & log2foldchange is greater than 0
  upregulated_data <-subset(dataset, FDR <=0.05 & logFC >0)
  
  # subset data for significant genes that are downregulated genes
  downregulated_data <-subset(dataset, FDR <=0.05 & logFC <0)
  
  # put the up and downregulated dataframes into a list
  up_and_down_list <- list(upregulated_data, downregulated_data)
  
  # set this up ahead of time
  combined_results <- data.frame()
  
  for (input_data in up_and_down_list) {
    
    logFC_direction <- ifelse(identical(input_data, upregulated_data), "upregulated", "downregulated")
    
    # create a new column called "Case.Cellstates" in upregulated_data from the "cell" column
    input_data$Case.CellStates <- substr(input_data$cell, 1, regexpr("_", input_data$cell) - 1)
    
    # use dcast to melt the data
    transformed_matrix <- dcast(input_data, gene ~ Case.CellStates, value.var = "logFC")
    
    # make each cell type into a separate table
    # create a list of data frames, each containing a single cell type column
    cell_type_columns <- colnames(transformed_matrix)[-1]  # get the column names except the first one (gene column)
    
    cell_type_list <- list()
    
    # loop through each cell type column and create a new data frame for each cell type
    for (cell_type in cell_type_columns) {
      # create a new data frame with gene and the current cell type column
      new_data_frame <- data.frame(Gene = transformed_matrix$gene, Expression = transformed_matrix[[cell_type]])
      
      # set a new name for the data frame based on the cell type
      new_data_frame_name <- cell_type
      
      # add the new data frame to the list with the appropriate name
      cell_type_list[[new_data_frame_name]] <- new_data_frame
    }
    
    # iterate over each data frame in the cell_type_list to remove nas
    for (CellState in names(cell_type_list)) {
      # get the current data frame
      current_df <- cell_type_list[[CellState]]
      
      # remove rows with missing values (NAs)
      complete_rows <- complete.cases(current_df)
      current_df <- current_df[complete_rows, ]
      
      # do GO analysis with gene symbols for each cell type
      gene_symbols <- current_df$Gene
      go_enrichment <- enrichGO(gene = gene_symbols, OrgDb = org.Hs.eg.db, keyType = "SYMBOL")
      
      if (nrow(go_enrichment) > 0) {
      
        # convert the go_enrichment object to a data frame
        go_enrichment_df <- as.data.frame(go_enrichment)
        
        # add new columns for CellState, logFC_direction, and dataset_name
        go_enrichment_df$CellState <- CellState
        go_enrichment_df$logFC_direction <- logFC_direction
        go_enrichment_df$dataset_name <- gsub(".rds", "", file_names[i])
        
        # append the go_enrichment_df to the combined_results data frame
        combined_results <- rbind(combined_results, go_enrichment_df)
  
      }
    }
  }
  
  # construct the output CSV file name for this dataset
  dataset_name <- gsub(".rds", "", file_names[i])
  output_file <- paste0("GO_analysis_", dataset_name, ".csv")
  
  # write the combined_results to the output CSV file
  write.csv(combined_results, file = output_file, row.names = FALSE)
  
}
