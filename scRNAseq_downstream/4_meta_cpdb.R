# SETWD OF WHERE YOU WANT OUTPUT    &    SELECT WHICH 'THE_FILE'   &   HIGHLIGHT ALL & RUN #
# run on R VERSON 4.2.0

packages <- c('readxl', 'ggplot2', 'dplyr', 'tidyverse', 'writexl', 'openxlsx')
lapply(packages, library, character.only=TRUE, quietly = TRUE) 
partial_path <- '.../CellphoneDB_analysis/cpdb_outputs/'
the_file <- '/statistical_analysis_significant_means_.txt_table_chemokine.csv'
#the_file <- '/statistical_analysis_significant_means_.txt_table_ice.csv'

####################################################################################################################
####################################################################################################################
###################################                   FUNCTIONS                   ##################################
####################################################################################################################
####################################################################################################################

# for processing multiple cpdb csvs for comparison
process_dfs <-function(data_frames) {
  # combine the dataframes
  combined_df <- bind_rows(data_frames, .id = "dataset")
  # filter for columns that will be assessed for matches and uniqueness (get rows of columns like cell1_cell2)
  df <- combined_df %>% select(dataset, Value, ligand, receptor, cell1, cell2)
  # combine dataset and its corresponding value for later needs
  df$dataset <- paste(df$dataset, df$Value, sep = "---")
  # indicate which columns will be checked for uniqueness
  columns_to_check <- c("ligand", "receptor", "cell1", "cell2")
  # make new column "match" with dataset names having matches
  df_matched <- df %>%
    group_by(across(all_of(columns_to_check))) %>%
    summarize(match = paste(dataset, collapse = ", ")) %>%
    ungroup()
  return(df_matched)
}

###########################################################################
# for comparing 3 disease datasets to ID common interactions
compare3Disease <- function(disease1, disease2, disease3) {
  # make list of dataframes
  data_frames <- list(disease1 = disease1, disease2 = disease2, disease3 = disease3)
  # combine the dataframes
  combined_df <- bind_rows(data_frames, .id = "dataset")
  # filter for columns that will be assessed for matches and uniqueness (get rows of columns like cell1_cell2)
  df <- combined_df %>% select(dataset, Value, ligand, receptor, cell1, cell2)
  # combine dataset and its corresponding value for later needs
  df$dataset <- paste(df$dataset, df$Value, sep = "---")
  # indicate which columns will be checked for uniqueness
  columns_to_check <- c("ligand", "receptor", "cell1", "cell2")
  # make new column "match" with dataset names having matches
  df_matched <- df %>%
    group_by(across(all_of(columns_to_check))) %>%
    summarize(match = paste(dataset, collapse = ", ")) %>%
    ungroup()
  # split the matches based on dataset name
  split_match <- strsplit(df_matched$match, ",")
  # make new columns called dataset1 or 2 or 3 and fill them with what's in split_match
  df_matched$dataset1 <- sapply(split_match, function(x) ifelse(length(x) >= 1, x[1], NA))
  df_matched$dataset2 <- sapply(split_match, function(x) ifelse(length(x) >= 2, x[2], NA))
  df_matched$dataset3 <- sapply(split_match, function(x) ifelse(length(x) >= 3, x[3], NA))
  # make new columns called value1, 2 or 3 and fill them with the cpdb values for dataset1, 2 or 3
  df_matched$value1 <- sapply(strsplit(df_matched$dataset1, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
  df_matched$value2 <- sapply(strsplit(df_matched$dataset2, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
  df_matched$value3 <- sapply(strsplit(df_matched$dataset3, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
  # if there is an entity in dataset3 column, then make a new column called average_value and take the average of the values in columns value 1, 2, and 3
  # (first values need to be numeric)
  df_matched <- df_matched %>% mutate(across(c(value1, value2, value3), as.numeric))
  df_matched$average_value <- ifelse(!is.na(df_matched$value3), rowMeans(df_matched[, c("value1", "value2", "value3")], na.rm = TRUE), NA)
  # make df with relevant columns
  df_matched<-df_matched[,c('ligand', 'receptor', 'cell1', 'cell2', 'average_value')]
  df_matched <- df_matched[!is.na(df_matched$average_value), ]
  disease_commons <- df_matched
  print('disease_commons is ready')
  return(disease_commons)
}

###########################################################################
# for comparing 2 disease datasets to ID common interactions
compare2Disease <- function(disease1, disease2) {
  # make list of dataframes
  data_frames <- list(disease1 = disease1, disease2 = disease2)
  # combine the dataframes
  combined_df <- bind_rows(data_frames, .id = "dataset")
  # filter for columns that will be assessed for matches and uniqueness (get rows of columns like cell1_cell2)
  df <- combined_df %>% select(dataset, Value, ligand, receptor, cell1, cell2)
  # combine dataset and its corresponding value for later needs
  df$dataset <- paste(df$dataset, df$Value, sep = "---")
  # indicate which columns will be checked for uniqueness
  columns_to_check <- c("ligand", "receptor", "cell1", "cell2")
  # make new column "match" with dataset names having matches
  df_matched <- df %>%
    group_by(across(all_of(columns_to_check))) %>%
    summarize(match = paste(dataset, collapse = ", ")) %>%
    ungroup()
  # split the matches based on dataset name
  split_match <- strsplit(df_matched$match, ",")
  # make new columns called dataset1 or 2 or 3 and fill them with what's in split_match
  df_matched$dataset1 <- sapply(split_match, function(x) ifelse(length(x) >= 1, x[1], NA))
  df_matched$dataset2 <- sapply(split_match, function(x) ifelse(length(x) >= 2, x[2], NA))
  # make new columns called value1, 2 or 3 and fill them with the cpdb values for dataset1, 2 or 3
  df_matched$value1 <- sapply(strsplit(df_matched$dataset1, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
  df_matched$value2 <- sapply(strsplit(df_matched$dataset2, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
  # if there is an entity in dataset3 column, then make a new column called average_value and take the average of the values in columns value 1, 2, and 3
  # (first values need to be numeric)
  df_matched <- df_matched %>% mutate(across(c(value1, value2), as.numeric))
  df_matched$average_value <- ifelse(!is.na(df_matched$value2), rowMeans(df_matched[, c("value1", "value2")], na.rm = TRUE), NA)
  # make df with relevant columns
  df_matched<-df_matched[,c('ligand', 'receptor', 'cell1', 'cell2', 'average_value')]
  df_matched <- df_matched[!is.na(df_matched$average_value), ]
  disease_commons <- df_matched
  print('disease_commons is ready')
  return(disease_commons)
}

###########################################################################
# for comparing 3 disease datasets to ID common interactions across 2 datasets 
compare3Disease_UC <- function(disease1, disease2, disease3) {
  # make list of dataframes
  data_frames <- list(disease1 = disease1, disease2 = disease2, disease3 = disease3)
  # combine the dataframes
  combined_df <- bind_rows(data_frames, .id = "dataset")
  # filter for columns that will be assessed for matches and uniqueness (get rows of columns like cell1_cell2)
  df <- combined_df %>% select(dataset, Value, ligand, receptor, cell1, cell2)
  # combine dataset and its corresponding value for later needs
  df$dataset <- paste(df$dataset, df$Value, sep = "---")
  # indicate which columns will be checked for uniqueness
  columns_to_check <- c("ligand", "receptor", "cell1", "cell2")
  # make new column "match" with dataset names having matches
  df_matched <- df %>%
    group_by(across(all_of(columns_to_check))) %>%
    summarize(match = paste(dataset, collapse = ", ")) %>%
    ungroup()
  # split the matches based on dataset name
  split_match <- strsplit(df_matched$match, ",")
  # make new columns called dataset1 or 2 or 3 and fill them with what's in split_match
  df_matched$dataset1 <- sapply(split_match, function(x) ifelse(length(x) >= 1, x[1], NA))
  df_matched$dataset2 <- sapply(split_match, function(x) ifelse(length(x) >= 2, x[2], NA))
  # make new columns called value1, 2 or 3 and fill them with the cpdb values for dataset1, 2 or 3
  df_matched$value1 <- sapply(strsplit(df_matched$dataset1, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
  df_matched$value2 <- sapply(strsplit(df_matched$dataset2, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
  # if there is an entity in dataset3 column, then make a new column called average_value and take the average of the values in columns value 1, 2, and 3
  # (first values need to be numeric)
  df_matched <- df_matched %>% mutate(across(c(value1, value2), as.numeric))
  df_matched$average_value <- ifelse(!is.na(df_matched$value2), rowMeans(df_matched[, c("value1", "value2")], na.rm = TRUE), NA)
  # make df with relevant columns
  df_matched<-df_matched[,c('ligand', 'receptor', 'cell1', 'cell2', 'average_value')]
  df_matched <- df_matched[!is.na(df_matched$average_value), ]
  disease_commons <- df_matched
  print('disease_commons is ready')
  return(disease_commons)
}

###########################################################################
# for comparing the common disease interactions to the common healthy interactions
DvH_compare <- function(healthy, disease) {
  healthy$tissue <- 'healthy'
  # make list of dataframes
  data_frames <- list(healthy = healthy, disease = disease)
  # combine the dataframes
  combined_df <- bind_rows(data_frames, .id = "tissue")
  # rename df as combined_df
  df <- combined_df
  # combine dataset and its corresponding value for later needs
  df$dataset <- paste(df$tissue, df$average_value, sep = "---")
  # indicate which columns will be checked for uniqueness
  columns_to_check <- c("ligand", "receptor", "cell1", "cell2")
  # make new column "match" with dataset names having matches
  df_matched <- df %>%
    group_by(across(all_of(columns_to_check))) %>%
    summarize(match = paste(dataset, collapse = ", ")) %>%
    ungroup()
  # split the matches based on dataset name
  split_match <- strsplit(df_matched$match, ",")
  # make new columns called dataset1 or 2 or 3 and fill them with what's in split_match
  df_matched$dataset1 <- sapply(split_match, function(x) ifelse(length(x) >= 1, x[1], NA))
  df_matched$dataset2 <- sapply(split_match, function(x) ifelse(length(x) >= 2, x[2], NA))
  # make new columns called value1, 2 or 3 and fill them with the cpdb values for dataset1, 2 or 3
  df_matched$value1 <- sapply(strsplit(df_matched$dataset1, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
  df_matched$value2 <- sapply(strsplit(df_matched$dataset2, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
  # in columns dataset1, dataset2, if not filled with na, just keep whatever is written before the string '---' 
  df_matched$dataset1 <- ifelse(!is.na(df_matched$dataset1), gsub("---.*", "", df_matched$dataset1), NA)
  df_matched$dataset2 <- ifelse(!is.na(df_matched$dataset2), gsub("---.*", "", df_matched$dataset2), NA)
  #if dataset1 is filled with the word disease and dataset2 is filled with NA, make a column called disease only and fill it with whatever is in value1
  df_matched$disease_only <- ifelse(df_matched$dataset1 == "disease" & is.na(df_matched$dataset2), df_matched$value1, NA)
  # make disease direction column if interaction is taking place both in healthy and disease
  df_matched <- df_matched %>% mutate(across(c(value1, value2), as.numeric))
  df_matched <- df_matched %>% mutate(disease_direction = ifelse(!is.na(dataset1) & !is.na(dataset2), value2 - value1, NA))
  # make df with relevant columns
  df_matched <-df_matched[,c('ligand', 'receptor', 'cell1', 'cell2', 'disease_only', 'disease_direction')]
  DvH_df <- df_matched
  print('DvH_df is ready')
  return(DvH_df)
}

###########################################################################
# for comparing across the common datasets for either healthy or disease
metaTissue <- function(data_frames) {
  # combine the dataframes
  combined_df <- bind_rows(data_frames, .id = "tissue")
  df <- combined_df 
  # combine dataset and its corresponding value for later needs
  df$tissue <- paste(df$tissue, df$average_value, sep = "---")
  # indicate which columns will be checked for uniqueness
  columns_to_check <- c("ligand", "receptor", "cell1", "cell2")
  # make new column "match" with dataset names having matches
  df_matched <- df %>%
    group_by(across(all_of(columns_to_check))) %>%
    summarize(match = paste(tissue, collapse = ", ")) %>%
    ungroup()
  # split the matches based on dataset name
  split_match <- strsplit(df_matched$match, ",")
  # make new columns called dataset1 or 2 or 3 and fill them with what's in split_match
  df_matched$tissue1 <- sapply(split_match, function(x) ifelse(length(x) >= 1, x[1], NA))
  df_matched$tissue2 <- sapply(split_match, function(x) ifelse(length(x) >= 2, x[2], NA))
  df_matched$tissue3 <- sapply(split_match, function(x) ifelse(length(x) >= 3, x[3], NA))
  df_matched$tissue4 <- sapply(split_match, function(x) ifelse(length(x) >= 4, x[4], NA))
  df_matched$tissue5 <- sapply(split_match, function(x) ifelse(length(x) >= 5, x[5], NA))
  df_matched$tissue6 <- sapply(split_match, function(x) ifelse(length(x) >= 6, x[6], NA))
  df_matched$tissue7 <- sapply(split_match, function(x) ifelse(length(x) >= 7, x[7], NA))
  df_matched$tissue8 <- sapply(split_match, function(x) ifelse(length(x) >= 8, x[8], NA))
  df_matched$tissue9 <- sapply(split_match, function(x) ifelse(length(x) >= 9, x[9], NA))
  df_matched$tissue10 <- sapply(split_match, function(x) ifelse(length(x) >= 10, x[10], NA))
  # make new columns called value1, 2 or 3 and fill them with the cpdb values for dataset1, 2 or 3
  # (if there are no matches along 4 datasets, there will be an error - that's fine)
  df_matched$value1 <- sapply(strsplit(df_matched$tissue1, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
  df_matched$value2 <- sapply(strsplit(df_matched$tissue2, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
  df_matched$value3 <- sapply(df_matched$tissue3, function(x) {
    if (is.na(x) || x == "") {
      NA
    } else {
      parts <- strsplit(x, "---")
      if (length(parts[[1]]) >= 2) {
        parts[[1]][2]
      } else {
        NA
      } } })
  df_matched$value4 <- sapply(df_matched$tissue4, function(x) {
    if (is.na(x) || x == "") {
      NA
    } else {
      parts <- strsplit(x, "---")
      if (length(parts[[1]]) >= 2) {
        parts[[1]][2]
      } else {
        NA
      } } })
  df_matched$value5 <- sapply(df_matched$tissue5, function(x) {
    if (is.na(x) || x == "") {
      NA
    } else {
      parts <- strsplit(x, "---")
      if (length(parts[[1]]) >= 2) {
        parts[[1]][2]
      } else {
        NA
      } } })
  df_matched$value6 <- sapply(df_matched$tissue6, function(x) {
    if (is.na(x) || x == "") {
      NA
    } else {
      parts <- strsplit(x, "---")
      if (length(parts[[1]]) >= 2) {
        parts[[1]][2]
      } else {
        NA
      } } })
  df_matched$value7 <- sapply(df_matched$tissue7, function(x) {
    if (is.na(x) || x == "") {
      NA
    } else {
      parts <- strsplit(x, "---")
      if (length(parts[[1]]) >= 2) {
        parts[[1]][2]
      } else {
        NA
      } } })
  df_matched$value8 <- sapply(df_matched$tissue8, function(x) {
    if (is.na(x) || x == "") {
      NA
    } else {
      parts <- strsplit(x, "---")
      if (length(parts[[1]]) >= 2) {
        parts[[1]][2]
      } else {
        NA
      } } })
  df_matched$value9 <- sapply(df_matched$tissue9, function(x) {
    if (is.na(x) || x == "") {
      NA
    } else {
      parts <- strsplit(x, "---")
      if (length(parts[[1]]) >= 2) {
        parts[[1]][2]
      } else {
        NA
      } } })
  df_matched$value10 <- sapply(df_matched$tissue10, function(x) {
    if (is.na(x) || x == "") {
      NA
    } else {
      parts <- strsplit(x, "---")
      if (length(parts[[1]]) >= 2) {
        parts[[1]][2]
      } else {
        NA
      } } })
  # get average_values
  # (first values need to be numeric)
  df_matched <- df_matched %>% mutate(across(c(value1, value2, value3, value4, value5, value6, value7, value8, value9, value10), as.numeric))
  df_matched <- df_matched %>% rowwise() %>%
    mutate( average_value = if (!is.na(value10)) {
      mean(c(value1, value2, value3, value4, value5, value6, value7, value8, value9, value10), na.rm = TRUE)
    } else if (!is.na(value9)) {
      mean(c(value1, value2, value3, value4, value5, value6, value7, value8, value9), na.rm = TRUE)
    } else if (!is.na(value8)) {
      mean(c(value1, value2, value3, value4, value5, value6, value7, value8), na.rm = TRUE)
    } else if (!is.na(value7)) {
      mean(c(value1, value2, value3, value4, value5, value6, value7), na.rm = TRUE)
    } else if (!is.na(value6)) {
      mean(c(value1, value2, value3, value4, value5, value6), na.rm = TRUE)
    } else if (!is.na(value5)) {
      mean(c(value1, value2, value3, value4, value5), na.rm = TRUE)
    } else if (!is.na(value4)) {
      mean(c(value1, value2, value3, value4), na.rm = TRUE)
    } else if (!is.na(value3)) {
      mean(c(value1, value2, value3), na.rm = TRUE)
    } else if (!is.na(value2)) {
      mean(c(value1, value2), na.rm = TRUE)
    } else {
      value1 })
  
  # prep for final df....
  # in columns dataset1, dataset2, dataset3 and dataset4, if not filled with na, just keep whatever is written before the string '---' 
  df_matched$tissue1 <- ifelse(!is.na(df_matched$tissue1), gsub("---.*", "", df_matched$tissue1), NA)
  df_matched$tissue2 <- ifelse(!is.na(df_matched$tissue2), gsub("---.*", "", df_matched$tissue2), NA)
  df_matched$tissue3 <- ifelse(!is.na(df_matched$tissue3), gsub("---.*", "", df_matched$tissue3), NA)
  df_matched$tissue4 <- ifelse(!is.na(df_matched$tissue4), gsub("---.*", "", df_matched$tissue4), NA)
  df_matched$tissue5 <- ifelse(!is.na(df_matched$tissue5), gsub("---.*", "", df_matched$tissue5), NA)
  df_matched$tissue6 <- ifelse(!is.na(df_matched$tissue6), gsub("---.*", "", df_matched$tissue6), NA)
  df_matched$tissue7 <- ifelse(!is.na(df_matched$tissue7), gsub("---.*", "", df_matched$tissue7), NA)
  df_matched$tissue8 <- ifelse(!is.na(df_matched$tissue8), gsub("---.*", "", df_matched$tissue8), NA)
  df_matched$tissue9 <- ifelse(!is.na(df_matched$tissue9), gsub("---.*", "", df_matched$tissue9), NA)
  df_matched$tissue10 <- ifelse(!is.na(df_matched$tissue10), gsub("---.*", "", df_matched$tissue10), NA)
  
  # grab needed columns
  columns_to_exclude <- c('match')
  df_matched <- df_matched[, !(names(df_matched) %in% columns_to_exclude)]
  return(df_matched)
}
####################################################################################################################
######     1      ############## CHEMOKINE L-R INTX IN HEALTHY SKIN ###############################################
####################################################################################################################

AD_Skin_GSE147424 <- read.csv(paste0(partial_path, 'AD_Skin_GSE147424/Healthy', the_file))
AD_Skin_GSE153760 <- read.csv(paste0(partial_path, 'AD_Skin_GSE153760/HC', the_file))
PSO_Skin_GSE173706 <- read.csv(paste0(partial_path, 'PSO_Skin_GSE173706/Control', the_file))
PSO_Skin_GSE220116 <- read.csv(paste0(partial_path, 'PSO_Skin_GSE220116/Control__pre_tx', the_file))
Skin_Cell_Atlas <- read.csv(paste0(partial_path, 'Skin_Cell_Atlas/Healthy__non_lesion', the_file))

# make list of dataframes
data_frames <- list(AD_Skin_GSE147424 = AD_Skin_GSE147424, AD_Skin_GSE153760 = AD_Skin_GSE153760,
                    PSO_Skin_GSE173706 = PSO_Skin_GSE173706, PSO_Skin_GSE220116 = PSO_Skin_GSE220116,
                    Skin_Cell_Atlas = Skin_Cell_Atlas)

# combine the dataframes
combined_df <- bind_rows(data_frames, .id = "dataset")
# perform processing function
df_matched <- process_dfs(data_frames)
# split the matches based on dataset name
split_match <- strsplit(df_matched$match, ",")
# make new columns called dataset1 or 2 or 3 and fill them with what's in split_match
df_matched$dataset1 <- sapply(split_match, function(x) ifelse(length(x) >= 1, x[1], NA))
df_matched$dataset2 <- sapply(split_match, function(x) ifelse(length(x) >= 2, x[2], NA))
df_matched$dataset3 <- sapply(split_match, function(x) ifelse(length(x) >= 3, x[3], NA))
df_matched$dataset4 <- sapply(split_match, function(x) ifelse(length(x) >= 4, x[4], NA))
df_matched$dataset5 <- sapply(split_match, function(x) ifelse(length(x) >= 5, x[5], NA))
# make new columns called value1, 2 or 3 and fill them with the cpdb values for dataset1, 2 or 3
df_matched$value1 <- sapply(strsplit(df_matched$dataset1, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
df_matched$value2 <- sapply(strsplit(df_matched$dataset2, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
df_matched$value3 <- sapply(strsplit(df_matched$dataset3, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
df_matched$value4 <- sapply(strsplit(df_matched$dataset4, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
df_matched$value5 <- sapply(strsplit(df_matched$dataset5, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
# if there is an entity in dataset3 column, then make a new column called average_value and take the average of the values in columns value 1, 2, and 3
# (first values need to be numeric)
df_matched <- df_matched %>% mutate(across(c(value1, value2, value3, value4, value5), as.numeric))
df_matched$average_value <- ifelse(!is.na(df_matched$value5), rowMeans(df_matched[, c("value1", "value2", "value3", 'value4', 'value5')], na.rm = TRUE), NA)
# make df with relevant columns
df_matched<-df_matched[,c('ligand', 'receptor', 'cell1', 'cell2', 'average_value')]
df_matched <- df_matched[!is.na(df_matched$average_value), ]

###### specify name of working file based on tissue and make a column named tissue
skin <-df_matched
skin$tissue <- 'skin'

####################################################################################################################
#########     2     ############# CHEMOKINE L-R INTX IN HEALTHY LUNG ###############################################
####################################################################################################################
COPD_Lung_GSE136831 <- read.csv(paste0(partial_path, 'COPD_Lung_GSE136831/Control', the_file))
COPD_Lung_GSE171541 <- read.csv(paste0(partial_path, 'COPD_Lung_GSE171541/control', the_file))
Human_Lung_Cell_Atlas <- read.csv(paste0(partial_path, 'Human_Lung_Cell_Atlas/normal', the_file))
ILD_Lung_GSE1122960 <- read.csv(paste0(partial_path, 'ILD_Lung_GSE1122960/Control', the_file))
ILD_Lung_GSE135893 <- read.csv(paste0(partial_path, 'ILD_Lung_GSE135893/Control', the_file))

# make list of dataframes
data_frames <- list(COPD_Lung_GSE136831 = COPD_Lung_GSE136831, COPD_Lung_GSE171541 = COPD_Lung_GSE171541,
                    ILD_Lung_GSE1122960 = ILD_Lung_GSE1122960, ILD_Lung_GSE135893 = ILD_Lung_GSE135893, 
                    Human_Lung_Cell_Atlas = Human_Lung_Cell_Atlas)
# combine the dataframes
combined_df <- bind_rows(data_frames, .id = "dataset")
# perform processing function
df_matched <- process_dfs(data_frames)
# split the matches based on dataset name
split_match <- strsplit(df_matched$match, ",")
# make new columns called dataset1 or 2 or 3 and fill them with what's in split_match
df_matched$dataset1 <- sapply(split_match, function(x) ifelse(length(x) >= 1, x[1], NA))
df_matched$dataset2 <- sapply(split_match, function(x) ifelse(length(x) >= 2, x[2], NA))
df_matched$dataset3 <- sapply(split_match, function(x) ifelse(length(x) >= 3, x[3], NA))
df_matched$dataset4 <- sapply(split_match, function(x) ifelse(length(x) >= 4, x[4], NA))
df_matched$dataset5 <- sapply(split_match, function(x) ifelse(length(x) >= 5, x[5], NA))
# make new columns called value1, 2 or 3 and fill them with the cpdb values for dataset1, 2 or 3
df_matched$value1 <- sapply(strsplit(df_matched$dataset1, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
df_matched$value2 <- sapply(strsplit(df_matched$dataset2, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
df_matched$value3 <- sapply(strsplit(df_matched$dataset3, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
df_matched$value4 <- sapply(strsplit(df_matched$dataset4, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
df_matched$value5 <- sapply(strsplit(df_matched$dataset5, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
# if there is an entity in dataset3 column, then make a new column called average_value and take the average of the values in columns value 1, 2, and 3
# (first values need to be numeric)
# (for lung dataset, I averaged value3 so that if 3 datasets had a commonality [rather than 4], it was averaged)
df_matched <- df_matched %>% mutate(across(c(value1, value2, value3, value4, value5), as.numeric))
df_matched$average_value <- ifelse(!is.na(df_matched$value5), rowMeans(df_matched[, c("value1", "value2", "value3", 'value4','value5')], na.rm = TRUE), NA)
# make df with relevant columns
df_matched <-df_matched[,c('ligand', 'receptor', 'cell1', 'cell2', 'average_value')]
df_matched <- df_matched[!is.na(df_matched$average_value), ]

###### specify name of working file based on tissue and make a column named tissue
lung <-df_matched
lung$tissue <- 'lung'

####################################################################################################################
##########   3  ############## CHEMOKINE L-R INTX IN HEALTHY COLON  ###############################################
####################################################################################################################
UC_Colon_GSE116222 <- read.csv(paste0(partial_path, 'UC_Colon_GSE116222/Healthy', the_file))
UC_Colon_GSE231993 <- read.csv(paste0(partial_path, 'UC_Colon_GSE231993/Control', the_file))
UC_Colon_SCP259 <- read.csv(paste0(partial_path, 'UC_Colon_SCP259/Healthy', the_file))

# make list of dataframes
data_frames <- list(UC_Colon_SCP259 = UC_Colon_SCP259, UC_Colon_GSE116222 = UC_Colon_GSE116222, UC_Colon_GSE231993 = UC_Colon_GSE231993)
# combine the dataframes
combined_df <- bind_rows(data_frames, .id = "dataset")
# perform processing function
df_matched <- process_dfs(data_frames)
# split the matches based on dataset name
split_match <- strsplit(df_matched$match, ",")
# make new columns called dataset1 or 2 or 3 and fill them with what's in split_match
df_matched$dataset1 <- sapply(split_match, function(x) ifelse(length(x) >= 1, x[1], NA))
df_matched$dataset2 <- sapply(split_match, function(x) ifelse(length(x) >= 2, x[2], NA))
df_matched$dataset3 <- sapply(split_match, function(x) ifelse(length(x) >= 3, x[2], NA))
# make new columns called value1, 2 or 3 and fill them with the cpdb values for dataset1, 2 or 3
df_matched$value1 <- sapply(strsplit(df_matched$dataset1, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
df_matched$value2 <- sapply(strsplit(df_matched$dataset2, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
df_matched$value3 <- sapply(strsplit(df_matched$dataset3, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
# if there is an entity in dataset3 column, then make a new column called average_value and take the average of the values in columns value 1, 2, and 3
# (first values need to be numeric)
df_matched <- df_matched %>% mutate(across(c(value1, value2, value3), as.numeric))
df_matched$average_value <- ifelse(!is.na(df_matched$value3), rowMeans(df_matched[, c("value1", "value2", 'value3')], na.rm = TRUE), NA)
# make df with relevant columns
df_matched <-df_matched[,c('ligand', 'receptor', 'cell1', 'cell2', 'average_value')]
df_matched <- df_matched[!is.na(df_matched$average_value), ]

###### specify name of working file based on tissue and make a column named tissue
colon <-df_matched
colon$tissue <- 'colon'

####################################################################################################################
########  4    ################ CHEMOKINE L-R INTX IN HEALTHY KIDNEY ###############################################
####################################################################################################################
IgAN_Cell_Reports_Zheng_Kidney <- read.csv(paste0(partial_path, 'IgAN_Cell_Reports_Zheng_Kidney/normal_control', the_file))
SLE_Phase2_scRNA <- read.csv(paste0(partial_path, 'SLE_Phase2_scRNA/Control', the_file))

# make list of dataframes
data_frames <- list(IgAN_Cell_Reports_Zheng_Kidney = IgAN_Cell_Reports_Zheng_Kidney, SLE_Phase2_scRNA = SLE_Phase2_scRNA)
# combine the dataframes
combined_df <- bind_rows(data_frames, .id = "dataset")
# perform processing function
df_matched <- process_dfs(data_frames)
# split the matches based on dataset name
split_match <- strsplit(df_matched$match, ",")
# make new columns called dataset1 or 2 or 3 and fill them with what's in split_match
df_matched$dataset1 <- sapply(split_match, function(x) ifelse(length(x) >= 1, x[1], NA))
df_matched$dataset2 <- sapply(split_match, function(x) ifelse(length(x) >= 2, x[2], NA))
# make new columns called value1, 2 or 3 and fill them with the cpdb values for dataset1, 2 or 3
df_matched$value1 <- sapply(strsplit(df_matched$dataset1, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
df_matched$value2 <- sapply(strsplit(df_matched$dataset2, "---"), function(x) ifelse(length(x) >= 2, x[2], NA))
# if there is an entity in dataset3 column, then make a new column called average_value and take the average of the values in columns value 1, 2, and 3
# (first values need to be numeric)
df_matched <- df_matched %>% mutate(across(c(value1, value2), as.numeric))
df_matched$average_value <- ifelse(!is.na(df_matched$value2), rowMeans(df_matched[, c("value1", "value2")], na.rm = TRUE), NA)
# make df with relevant columns
df_matched <-df_matched[,c('ligand', 'receptor', 'cell1', 'cell2', 'average_value')]
df_matched <- df_matched[!is.na(df_matched$average_value), ]

###### specify name of working file based on tissue and make a column named tissue
kidney <-df_matched
kidney$tissue <- 'kidney'

####################################################################################################################
####################################################################################################################
####################################################################################################################
############################## CHEMOKINE L-R INTX IN HEALTHY TISSUES ###############################################
####################################################################################################################
####################################################################################################################
####################################################################################################################

# make list of dataframes
data_frames <- list(skin = skin, lung = lung, colon = colon, kidney = kidney)
# run the function
meta_healthy <- metaTissue(data_frames)

# save as csv 
filename <- paste(getwd(), "cellphonedb_chemokine_healthytissue.csv", sep = "/")
write.csv(meta_healthy, file = filename, row.names = FALSE)

##### csv processing 
csv <- read.csv("cellphonedb_chemokine_healthytissue.csv")

# remove any cell1/2 columns that have NonImmune or Unclassified
csv <- csv[csv$cell1 != "NonImmune",]
csv <- csv[csv$cell2 != "NonImmune",]
csv <- csv[csv$cell1 != "Unclassified",]
csv <- csv[csv$cell2 != "Unclassified",]

# write the csv to excel
write_xlsx(csv, "cellphonedb_chemokine_healthytissue.xlsx")

# read in the excel file
excel <- read_excel('cellphonedb_chemokine_healthytissue.xlsx')

# indicate excel file for workbook function
wb <- loadWorkbook("cellphonedb_chemokine_healthytissue.xlsx")

# make needed columns per tissue and add its corresponding value
tissues <- c("kidney", "colon", "skin", "lung")
excel[, tissues] <- NA

# put tissue1 value in the appropriate tissue column
for (tissue in tissues) {
  value <- excel$tissue1 == tissue
  excel[value, tissue] <- excel$value1[value]}

# put tissue2 value in the appropriate tissue column
for (tissue in tissues) {
  value <- !is.na(excel$tissue2) & excel$tissue2 == tissue
  excel[value, tissue] <- excel$value2[value]}

# put tissue3 value in the appropriate tissue column
for (tissue in tissues) {
  value <- !is.na(excel$tissue3) & excel$tissue3 == tissue
  excel[value, tissue] <- excel$value3[value]}

# put tissue4 value in the appropriate tissue column
for (tissue in tissues) {
  value <- !is.na(excel$tissue4) & excel$tissue4 == tissue
  excel[value, tissue] <- excel$value4[value]}

# filter rows where 'value2' is not NA (= commons)
commons <- excel[!is.na(excel$value2), ]

# add filtered data  to new sheet called "commons"
addWorksheet(wb, "commons")
writeData(wb, "commons", commons)

# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_healthytissue.xlsx", overwrite = T)

# grab all the rows where there is no value2, meaning it is tissue specific
specifics <- excel[is.na(excel$value2),]

# put tissue1 value in the appropriate tissue column
for (tissue in tissues) {
  value <- specifics$tissue1 == tissue
  specifics[value, tissue] <- specifics$value1[value]}

# filter rows where column 'kidney' is not NA, which = kidney only
kidney <- specifics[!is.na(specifics$kidney), ]

# filter rows where column 'colon' is not NA, which = colon only
colon <- specifics[!is.na(specifics$colon), ]

# filter rows where column 'skin' is not NA, which = skin only
skin <- specifics[!is.na(specifics$skin), ]

# filter rows where column 'lung' is not NA, which = lung only
lung <- specifics[!is.na(specifics$lung), ]

# add to new sheet and save updated excel file
addWorksheet(wb, "skin")
writeData(wb, "skin", skin)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_healthytissue.xlsx", overwrite = T)

# add to new sheet and save updated excel file
addWorksheet(wb, "colon")
writeData(wb, "colon", colon)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_healthytissue.xlsx", overwrite = T)

# add to new sheet and save updated excel file
addWorksheet(wb, "kidney")
writeData(wb, "kidney", kidney)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_healthytissue.xlsx", overwrite = T)

# add to new sheet and save updated excel file
addWorksheet(wb, "lung")
writeData(wb, "lung", lung)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_healthytissue.xlsx", overwrite = T)


####################################################################################################################
######     5      ############## CHEMOKINE L-R INTX IN AD NONLESION SKIN ###########################################
####################################################################################################################
AD_Skin_GSE147424 <- read.csv(paste0(partial_path, 'AD_Skin_GSE147424/Non_Lesional', the_file))
Skin_Cell_Atlas <- read.csv(paste0(partial_path, 'Skin_Cell_Atlas/Eczema__non_lesion', the_file))

disease1 <- AD_Skin_GSE147424
disease2 <- Skin_Cell_Atlas

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### ##### COMPARE HEALTHY TO DISEASE TISSUE #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

####  get what's commmon accross the datasets
disease <- compare2Disease(disease1, disease2)
disease$tissue <- 'disease'

#### get intx present in this disease that is absent in all healthy tissues' intx (regardless of tissue type)
DvH_data <- DvH_compare(meta_healthy, disease)
# remove rows where either disease_only is empty
AD_nonlesion_only <- DvH_data [!is.na(DvH_data$disease_only), ]

#### get directional data comparing to the healthy tissue of this disease
# indicate healthy data
healthy <- skin
DvH_data <- DvH_compare(healthy, disease)
# remove rows where disease_direction empty
AD_nonlesion_direction <- DvH_data [!is.na(DvH_data$disease_direction), ]

####################################################################################################################
######     6      ############## CHEMOKINE L-R INTX IN AD LESION SKIN ##############################################
####################################################################################################################
AD_Skin_GSE147424 <- read.csv(paste0(partial_path, 'AD_Skin_GSE147424/Lesional', the_file))
AD_Skin_GSE153760 <- read.csv(paste0(partial_path, 'AD_Skin_GSE153760/AD', the_file))
Skin_Cell_Atlas <- read.csv(paste0(partial_path, 'Skin_Cell_Atlas/Eczema__lesion', the_file))


disease1 <- AD_Skin_GSE147424
disease2 <- AD_Skin_GSE153760
disease3 <- Skin_Cell_Atlas

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### ##### COMPARE HEALTHY TO DISEASE TISSUE #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

####  get what's commmon accross the datasets
disease <- compare3Disease(disease1, disease2, disease3)
disease$tissue <- 'disease'

#### get intx present in this disease that is absent in all healthy tissues' intx (regardless of tissue type)
DvH_data <- DvH_compare(meta_healthy, disease)
# remove rows where disease_only empty
AD_lesion_only <- DvH_data[!is.na(DvH_data$disease_only), ]

#### get directional data comparing to the healthy tissue of this disease
# indicate healthy data
healthy <- skin
# run the function
DvH_data <- DvH_compare(healthy, disease)
# remove rows where disease_direction empty
AD_lesion_direction <- DvH_data[!is.na(DvH_data$disease_direction), ]

####################################################################################################################
######     7      ##############  L-R INTX IN PSO NONLESION SKIN          ##########################################
####################################################################################################################
PSO_Skin_GSE173706 <- read.csv(paste0(partial_path, 'PSO_Skin_GSE173706/PSO_non_lesion', the_file))
Skin_Cell_Atlas <- read.csv(paste0(partial_path, 'Skin_Cell_Atlas/Psoriasis__non_lesion', the_file))

disease1 <- PSO_Skin_GSE173706
disease2 <- Skin_Cell_Atlas

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### ##### COMPARE HEALTHY TO DISEASE TISSUE #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

####  get what's commmon accross the datasets
disease <- compare2Disease(disease1, disease2)
disease$tissue <- 'disease'

#### get intx present in this disease that is absent in all healthy tissues' intx (regardless of tissue type)
DvH_data <- DvH_compare(meta_healthy, disease)
# remove rows where either disease_only is empty
PSO_nonlesion_only <- DvH_data [!is.na(DvH_data$disease_only), ]

#### get directional data comparing to the healthy tissue of this disease
# indicate healthy data
healthy <- skin
# run the function
DvH_data <- DvH_compare(healthy, disease)
# remove rows where disease_direction empty
PSO_nonlesion_direction <- DvH_data [!is.na(DvH_data$disease_direction), ]

####################################################################################################################
######     8      ############## CHEMOKINE L-R INTX IN PSO LESION SKIN #############################################
####################################################################################################################
PSO_Skin_GSE173706 <- read.csv(paste0(partial_path, 'PSO_Skin_GSE173706/PSO_lesion', the_file))
PSO_Skin_GSE220116 <- read.csv(paste0(partial_path, 'PSO_Skin_GSE220116/PSO__pre_tx', the_file))
Skin_Cell_Atlas <- read.csv(paste0(partial_path, 'Skin_Cell_Atlas/Psoriasis__lesion', the_file))

disease1 <- PSO_Skin_GSE173706
disease2 <- PSO_Skin_GSE220116
disease3 <- Skin_Cell_Atlas

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### ##### COMPARE HEALTHY TO DISEASE TISSUE #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

####  get what's commmon accross the datasets
disease <- compare3Disease(disease1, disease2, disease3)
disease$tissue <- 'disease'

#### get intx present in this disease that is absent in all healthy tissues' intx (regardless of tissue type)
DvH_data <- DvH_compare(meta_healthy, disease)
# remove rows where disease_direction empty
PSO_lesion_only <- DvH_data [!is.na(DvH_data$disease_only), ]

#### get directional data comparing to the healthy tissue of this disease
# indicate healthy data
healthy <- skin
DvH_data <- DvH_compare(healthy, disease)
# remove rows where disease_direction empty
PSO_lesion_direction <- DvH_data [!is.na(DvH_data$disease_direction), ]

####################################################################################################################
######     9      ############## CHEMOKINE L-R INTX IN COPD LUNG  ##################################################
####################################################################################################################
COPD_Lung_GSE136831 <- read.csv(paste0(partial_path, 'COPD_Lung_GSE136831/COPD', the_file))
COPD_Lung_GSE171541 <- read.csv(paste0(partial_path, 'COPD_Lung_GSE171541/copd', the_file))

disease1 <- COPD_Lung_GSE136831
disease2 <- COPD_Lung_GSE171541

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### ##### COMPARE HEALTHY TO DISEASE TISSUE #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

####  get what's commmon accross the datasets
disease <- compare2Disease(disease1, disease2)
disease$tissue <- 'disease'

#### get intx present in this disease that is absent in all healthy tissues' intx (regardless of tissue type)
DvH_data <- DvH_compare(meta_healthy, disease)
# remove rows where either disease_only is empty
COPD_only <- DvH_data [!is.na(DvH_data$disease_only), ]

#### get directional data comparing to the healthy tissue of this disease
# indicate healthy data
healthy <- lung
DvH_data <- DvH_compare(healthy, disease)
# remove rows where disease_direction empty
COPD_direction <- DvH_data [!is.na(DvH_data$disease_direction), ]

####################################################################################################################
######     10     ############## CHEMOKINE L-R INTX IN IPF LUNG  ##################################################
####################################################################################################################
COPD_Lung_GSE136831 <- read.csv(paste0(partial_path, 'COPD_Lung_GSE136831/IPF', the_file))
ILD_Lung_GSE1122960 <- read.csv(paste0(partial_path, 'ILD_Lung_GSE1122960/IPF', the_file))
ILD_Lung_GSE135893 <- read.csv(paste0(partial_path, 'ILD_Lung_GSE135893/IPF', the_file))

disease1 <- COPD_Lung_GSE136831
disease2 <- ILD_Lung_GSE1122960
disease3 <- ILD_Lung_GSE135893

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### ##### COMPARE HEALTHY TO DISEASE TISSUE #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

####  get what's commmon accross the datasets
disease <- compare3Disease(disease1, disease2, disease3)
disease$tissue <- 'disease'

#### get intx present in this disease that is absent in all healthy tissues' intx (regardless of tissue type)
DvH_data <- DvH_compare(meta_healthy, disease)
# remove rows where either disease_only is empty
IPF_only <- DvH_data [!is.na(DvH_data$disease_only), ]

#### get directional data comparing to the healthy tissue of this disease
# indicate healthy data
healthy <- lung
DvH_data <- DvH_compare(healthy, disease)
# remove rows where disease_direction empty
IPF_direction <- DvH_data [!is.na(DvH_data$disease_direction), ]

####################################################################################################################
######     11     ############## CHEMOKINE L-R INTX IN DISEASE COLON UNINFLAMED  ###################################
####################################################################################################################
UC_Colon_GSE116222 <- read.csv(paste0(partial_path, 'UC_Colon_GSE116222/UC_Non_Inflamed', the_file))
UC_Colon_GSE231993 <- read.csv(paste0(partial_path, 'UC_Colon_GSE231993/UC_uninflamed', the_file))
UC_Colon_SCP259 <- read.csv(paste0(partial_path, 'UC_Colon_SCP259/UC_Non_Inflamed', the_file))

disease1 <- UC_Colon_GSE116222
disease2 <- UC_Colon_GSE231993
disease3 <- UC_Colon_SCP259 

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### ##### COMPARE HEALTHY TO DISEASE TISSUE #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

####  get what's commmon accross the datasets
disease <- compare3Disease_UC(disease1, disease2, disease3)
disease$tissue <- 'disease'

#### get intx present in this disease that is absent in all healthy tissues' intx (regardless of tissue type)
DvH_data <- DvH_compare(meta_healthy, disease)
# remove rows where either disease_only is empty
UC_uninflamed_only <- DvH_data [!is.na(DvH_data$disease_only), ]

#### get directional data comparing to the healthy tissue of this disease
# indicate healthy data
healthy <- colon
DvH_data <- DvH_compare(healthy, disease)
# remove rows where disease_direction empty
UC_uninflamed_direction <- DvH_data [!is.na(DvH_data$disease_direction), ]

####################################################################################################################
######     12     ############## CHEMOKINE L-R INTX IN DISEASE COLON INFLAMED  #####################################
####################################################################################################################

UC_Colon_GSE116222 <- read.csv(paste0(partial_path, 'UC_Colon_GSE116222/UC_Inflamed', the_file))
UC_Colon_GSE231993 <- read.csv(paste0(partial_path, 'UC_Colon_GSE231993/UC_inflamed', the_file))
UC_Colon_SCP259 <- read.csv(paste0(partial_path, 'UC_Colon_SCP259/UC_Inflamed', the_file))

disease1 <- UC_Colon_GSE116222
disease2 <- UC_Colon_GSE231993
disease3 <- UC_Colon_SCP259 

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### ##### COMPARE HEALTHY TO DISEASE TISSUE #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

####  get what's commmon accross the datasets
disease <- compare3Disease_UC(disease1, disease2, disease3)
disease$tissue <- 'disease'

#### get intx present in this disease that is absent in all healthy tissues' intx (regardless of tissue type)
DvH_data <- DvH_compare(meta_healthy, disease)
# remove rows where either disease_only is empty
UC_inflamed_only <- DvH_data [!is.na(DvH_data$disease_only), ]

#### get directional data comparing to the healthy tissue of this disease
# indicate healthy data
healthy <- colon
DvH_data <- DvH_compare(healthy, disease)
# remove rows where disease_direction empty
UC_inflamed_direction <- DvH_data [!is.na(DvH_data$disease_direction), ]

####################################################################################################################
######     13     ############## CHEMOKINE L-R INTX IN IGAN  #######################################################
####################################################################################################################

IgAN_Cell_Reports_Zheng_Kidney <- read.csv(paste0(partial_path, 'IgAN_Cell_Reports_Zheng_Kidney/IgAN', the_file))

kidney_disease <- IgAN_Cell_Reports_Zheng_Kidney[,c('ligand', 'receptor', 'cell1', 'cell2', 'Value')]
# rename Value to average_value
kidney_disease <- kidney_disease %>% rename(average_value = Value)
# indicate disease file
disease <-kidney_disease
disease$tissue <- 'disease'

#### get intx present in this disease that is absent in all healthy tissues' intx (regardless of tissue type)
DvH_data <- DvH_compare(meta_healthy, disease)
# remove rows where either disease_only is empty
IgAN_only <- DvH_data [!is.na(DvH_data$disease_only), ]

#### get directional data comparing to the healthy tissue of this disease
# indicate healthy data
healthy <- kidney
DvH_data <- DvH_compare(healthy, disease)
# remove rows where disease_direction empty
IgAN_direction <- DvH_data [!is.na(DvH_data$disease_direction), ]

####################################################################################################################
######     14     ############## CHEMOKINE L-R INTX IN SLE   #######################################################
####################################################################################################################
SLE_Phase2_scRNA <- read.csv(paste0(partial_path, 'SLE_Phase2_scRNA/SLE', the_file))

kidney_disease <- SLE_Phase2_scRNA[,c('ligand', 'receptor', 'cell1', 'cell2', 'Value')]
# rename Value to average_value
kidney_disease <- kidney_disease %>% rename(average_value = Value)
# indicate disease file
disease <-kidney_disease
disease$tissue <- 'disease'

#### get intx present in this disease that is absent in all healthy tissues' intx (regardless of tissue type)
DvH_data <- DvH_compare(meta_healthy, disease)
# remove rows where either disease_only is empty
SLE_only <- DvH_data [!is.na(DvH_data$disease_only), ]

#### get directional data comparing to the healthy tissue of this disease
# indicate healthy data
healthy <- kidney
DvH_data <- DvH_compare(healthy, disease)
# remove rows where disease_direction empty
SLE_direction <- DvH_data [!is.na(DvH_data$disease_direction), ]

####################################################################################################################
####################################################################################################################
####################################################################################################################
############################## CHEMOKINE L-R INTX ONLY IN DISEASE TISSUES ########################################## 
####################################################################################################################
####################################################################################################################
####################################################################################################################

# in order to run metaTissue function;
# remove the disease_direction column
AD_nonlesion_only  <- AD_nonlesion_only[,c('ligand', 'receptor', 'cell1', 'cell2', 'disease_only')]
AD_lesion_only  <- AD_lesion_only[,c('ligand', 'receptor', 'cell1', 'cell2', 'disease_only')]
PSO_nonlesion_only  <- PSO_nonlesion_only[,c('ligand', 'receptor', 'cell1', 'cell2', 'disease_only')]
PSO_lesion_only  <- PSO_lesion_only[,c('ligand', 'receptor', 'cell1', 'cell2', 'disease_only')]
COPD_only  <- COPD_only[,c('ligand', 'receptor', 'cell1', 'cell2', 'disease_only')]
IPF_only  <- IPF_only[,c('ligand', 'receptor', 'cell1', 'cell2', 'disease_only')]
UC_uninflamed_only  <- UC_uninflamed_only[,c('ligand', 'receptor', 'cell1', 'cell2', 'disease_only')]
UC_inflamed_only  <- UC_inflamed_only[,c('ligand', 'receptor', 'cell1', 'cell2', 'disease_only')]
IgAN_only  <- IgAN_only[,c('ligand', 'receptor', 'cell1', 'cell2', 'disease_only')]
SLE_only  <- SLE_only[,c('ligand', 'receptor', 'cell1', 'cell2', 'disease_only')]

# rename disease_only column to average_value 
AD_nonlesion_only <- AD_nonlesion_only %>% rename(average_value = disease_only)
AD_lesion_only <- AD_lesion_only %>% rename(average_value = disease_only)
PSO_nonlesion_only <- PSO_nonlesion_only %>% rename(average_value = disease_only)
PSO_lesion_only <- PSO_lesion_only %>% rename(average_value = disease_only)
COPD_only <- COPD_only %>% rename(average_value = disease_only)
IPF_only <- IPF_only %>% rename(average_value = disease_only)
UC_uninflamed_only <- UC_uninflamed_only %>% rename(average_value = disease_only)
UC_inflamed_only <- UC_inflamed_only %>% rename(average_value = disease_only)
IgAN_only <- IgAN_only %>% rename(average_value = disease_only)
SLE_only <-  SLE_only %>% rename(average_value = disease_only)

# make a column called tissue and fill with AD_only, IPF_only, etc
AD_nonlesion_only$tissue <- 'AD_nonlesion_only'
AD_lesion_only$tissue <- 'AD_lesion_only'
PSO_nonlesion_only$tissue <- 'PSO_nonlesion_only'
PSO_lesion_only$tissue <- 'PSO_lesion_only'
COPD_only$tissue <- 'COPD_only'
IPF_only$tissue <- 'IPF_only'
UC_uninflamed_only$tissue <- 'UC_uninflamed_only'
UC_inflamed_only$tissue <- 'UC_inflamed_only'
IgAN_only$tissue <- 'IgAN_only'
SLE_only$tissue <- 'SLE_only'

# make list of dataframes
data_frames <- list(AD_nonlesion_only = AD_nonlesion_only, AD_lesion_only = AD_lesion_only, 
                    PSO_nonlesion_only = PSO_nonlesion_only, PSO_lesion_only = PSO_lesion_only, 
                    COPD_only = COPD_only, IPF_only = IPF_only, 
                    UC_uninflamed_only = UC_uninflamed_only, UC_inflamed_only = UC_inflamed_only,
                    IgAN_only = IgAN_only, SLE_only = SLE_only)
# run the function
meta_disease <- metaTissue(data_frames)

# save csv
filename <- paste(getwd(), "cellphonedb_chemokine_diseaseonly.csv", sep = "/")
write.csv(meta_disease, file = filename, row.names = FALSE)

##### csv processing 
csv <- read.csv("cellphonedb_chemokine_diseaseonly.csv")

# remove any cell1/2 columns that have NonImmune
csv <- csv[csv$cell1 != "NonImmune",]
csv <- csv[csv$cell2 != "NonImmune",]
csv <- csv[csv$cell1 != "Unclassified",]
csv <- csv[csv$cell2 != "Unclassified",]

# write the csv to excel
write_xlsx(csv, "cellphonedb_chemokine_diseaseonly.xlsx")

# read in the excel file
excel <- read_excel('cellphonedb_chemokine_diseaseonly.xlsx')

# indicate excel file for workbook function
wb <- loadWorkbook("cellphonedb_chemokine_diseaseonly.xlsx")

# make needed columns per tissue and add its corresponding value
tissues <- c('AD_nonlesion_only', 'AD_lesion_only', 'PSO_nonlesion_only', 'PSO_lesion_only', 'COPD_only', 'IPF_only',
             'UC_uninflamed_only', 'UC_inflamed_only', 'IgAN_only', 'SLE_only')
excel[, tissues] <- NA

# put tissue1 value in the appropriate tissue column
for (tissue in tissues) {
  value <- excel$tissue1 == tissue
  excel[value, tissue] <- excel$value1[value]}

# put tissue2 value in the appropriate tissue column
for (tissue in tissues) {
  value <- !is.na(excel$tissue2) & excel$tissue2 == tissue
  excel[value, tissue] <- excel$value2[value]}

# put tissue3 value in the appropriate tissue column
for (tissue in tissues) {
  value <- !is.na(excel$tissue3) & excel$tissue3 == tissue
  excel[value, tissue] <- excel$value3[value]}

# put tissue4 value in the appropriate tissue column
for (tissue in tissues) {
  value <- !is.na(excel$tissue4) & excel$tissue4 == tissue
  excel[value, tissue] <- excel$value4[value]}

# put tissue5 value in the appropriate tissue column
for (tissue in tissues) {
  value <- !is.na(excel$tissue5) & excel$tissue5 == tissue
  excel[value, tissue] <- excel$value5[value]}

# put tissue6 value in the appropriate tissue column
for (tissue in tissues) {
  value <- !is.na(excel$tissue6) & excel$tissue6 == tissue
  excel[value, tissue] <- excel$value6[value]}

# put tissue7 value in the appropriate tissue column
for (tissue in tissues) {
  value <- !is.na(excel$tissue7) & excel$tissue7 == tissue
  excel[value, tissue] <- excel$value7[value]}

# put tissue8 value in the appropriate tissue column
for (tissue in tissues) {
  value <- !is.na(excel$tissue8) & excel$tissue8 == tissue
  excel[value, tissue] <- excel$value8[value]}

# put tissue9 value in the appropriate tissue column
for (tissue in tissues) {
  value <- !is.na(excel$tissue9) & excel$tissue9 == tissue
  excel[value, tissue] <- excel$value9[value]}

# put tissue10 value in the appropriate tissue column
for (tissue in tissues) {
  value <- !is.na(excel$tissue10) & excel$tissue10 == tissue
  excel[value, tissue] <- excel$value10[value]}

# filter rows where 'value2' is not NA (= commons)
commons <- excel[!is.na(excel$value2), ]

# make ADPSO_commons tab if tissue5 is NA and tissue1+2+3+4 say either AD/PSO lesional/nonlesional (or NA for tissue3/4)
ADPSO_commons <- commons[(is.na(commons$tissue5)) &
                           (commons$tissue1 %in% c("AD_nonlesion_only", "AD_lesion_only", "PSO_nonlesion_only", "PSO_lesion_only"))
                         & (commons$tissue2 %in% c("AD_nonlesion_only", "AD_lesion_only", "PSO_nonlesion_only", "PSO_lesion_only"))
                         & (commons$tissue3 %in% c("AD_nonlesion_only", "AD_lesion_only", "PSO_nonlesion_only", "PSO_lesion_only", NA)) 
                         & (commons$tissue4 %in% c("AD_nonlesion_only", "AD_lesion_only", "PSO_nonlesion_only", "PSO_lesion_only", NA))
                         , ]
addWorksheet(wb, "ADPSO_commons")
writeData(wb, "ADPSO_commons", ADPSO_commons)
saveWorkbook(wb, "cellphonedb_chemokine_diseaseonly.xlsx", overwrite = T)

# make COPDIPF_commons tab if tissue3 is NA and tissue1+2 say COPD/IPF
COPDIPF_commons <- commons[(is.na(commons$tissue3)) & (commons$tissue1 == "COPD_only") & (commons$tissue2 == "IPF_only"), ]
addWorksheet(wb, "COPDIPF_commons")
writeData(wb, "COPDIPF_commons", COPDIPF_commons)
saveWorkbook(wb, "cellphonedb_chemokine_diseaseonly.xlsx", overwrite = T)

# make SLEIGAN_commons tab if tissue3 is NA and tissue1+2 say IGAN/SLE
SLEIGAN_commons <- commons[(is.na(commons$tissue3)) & (commons$tissue1 == "IgAN_only") & (commons$tissue2 == "SLE_only"), ]
addWorksheet(wb, "SLEIGAN_commons")
writeData(wb, "SLEIGAN_commons", SLEIGAN_commons)
saveWorkbook(wb, "cellphonedb_chemokine_diseaseonly.xlsx", overwrite = T)

# make UC_commons tab if tissue3 is NA and tissue1+2 say UC_inflamed/noninflamed
UC_commons <- commons[(is.na(commons$tissue3)) & (commons$tissue1 == "UC_uninflamed_only") & (commons$tissue2 == "UC_inflamed_only"), ]
addWorksheet(wb, "UC_commons")
writeData(wb, "UC_commons", UC_commons)
saveWorkbook(wb, "cellphonedb_chemokine_diseaseonly.xlsx", overwrite = T)

# update commons
commons <- anti_join(commons, ADPSO_commons)
commons <- anti_join(commons, COPDIPF_commons)
commons <- anti_join(commons, SLEIGAN_commons)
commons <- anti_join(commons, UC_commons)

# add filtered data  to new sheet called "commons"
addWorksheet(wb, "commons")
writeData(wb, "commons", commons)
saveWorkbook(wb, "cellphonedb_chemokine_diseaseonly.xlsx", overwrite = T)

# grab all the rows where there is no value2, meaning it is tissue specific
specifics <- excel[is.na(excel$value2),]

# put tissue1 value in the appropriate tissue column
for (tissue in tissues) {
  value <- specifics$tissue1 == tissue
  specifics[value, tissue] <- specifics$value1[value]}

# filter rows where column [each disease] is not NA, which = disease_only
AD_nonlesion <- specifics[!is.na(specifics$AD_nonlesion_only), ]
AD_lesion <- specifics[!is.na(specifics$AD_lesion_only), ]
PSO_nonlesion <- specifics[!is.na(specifics$PSO_nonlesion_only), ]
PSO_lesion <- specifics[!is.na(specifics$PSO_lesion_only), ]
COPD <- specifics[!is.na(specifics$COPD_only), ]
IPF <- specifics[!is.na(specifics$IPF_only), ]
UC_uninflamed <- specifics[!is.na(specifics$UC_uninflamed_only), ]
UC_inflamed <- specifics[!is.na(specifics$UC_inflamed_only), ]
IGAN <- specifics[!is.na(specifics$IgAN_only), ]
SLE <- specifics[!is.na(specifics$SLE_only), ]

# add to new sheet and save updated excel file
addWorksheet(wb, "AD_nonlesion")
writeData(wb, "AD_nonlesion", AD_nonlesion)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_diseaseonly.xlsx", overwrite = T)

# add to new sheet and save updated excel file
addWorksheet(wb, "AD_lesion")
writeData(wb, "AD_lesion", AD_lesion)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_diseaseonly.xlsx", overwrite = T)

# add to new sheet and save updated excel file
addWorksheet(wb, "PSO_nonlesion")
writeData(wb, "PSO_nonlesion", PSO_nonlesion)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_diseaseonly.xlsx", overwrite = T)

# add to new sheet and save updated excel file
addWorksheet(wb, "PSO_lesion")
writeData(wb, "PSO_lesion", PSO_lesion)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_diseaseonly.xlsx", overwrite = T)

# add to new sheet and save updated excel file
addWorksheet(wb, "COPD")
writeData(wb, "COPD", COPD)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_diseaseonly.xlsx", overwrite = T)

# add to new sheet and save updated excel file
addWorksheet(wb, "IPF")
writeData(wb, "IPF", IPF)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_diseaseonly.xlsx", overwrite = T)

# add to new sheet and save updated excel file
addWorksheet(wb, "UC_uninflamed")
writeData(wb, "UC_uninflamed", UC_uninflamed)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_diseaseonly.xlsx", overwrite = T)

# add to new sheet and save updated excel file
addWorksheet(wb, "UC_inflamed")
writeData(wb, "UC_inflamed", UC_inflamed)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_diseaseonly.xlsx", overwrite = T)

# add to new sheet and save updated excel file
addWorksheet(wb, "IGAN")
writeData(wb, "IGAN", IGAN)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_diseaseonly.xlsx", overwrite = T)

# make sheet of top 50 intx
IGAN <- IGAN %>% arrange(desc(value1))
IGAN_top50 <- head(IGAN, 50)
addWorksheet(wb, "IGAN_top50")
writeData(wb, "IGAN_top50", IGAN_top50)
saveWorkbook(wb, "cellphonedb_chemokine_diseaseonly.xlsx", overwrite = T)

# add to new sheet and save updated excel file
addWorksheet(wb, "SLE")
writeData(wb, "SLE", SLE)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_diseaseonly.xlsx", overwrite = T)

# make sheet of top 50 intx
SLE <- SLE %>% arrange(desc(value1))
SLE_top50 <- head(SLE, 50)
addWorksheet(wb, "SLE_top50")
writeData(wb, "SLE_top50", SLE_top50)
saveWorkbook(wb, "cellphonedb_chemokine_diseaseonly.xlsx", overwrite = T)


####################################################################################################################
####################################################################################################################
####################################################################################################################
############################## CHEMOKINE L-R INTX UP/DOWN IN DISEASE TISSUES #######################################
####################################################################################################################
####################################################################################################################
####################################################################################################################

# in order to run metaTissue function;
# remove the disease_only column
AD_nonlesion_direction  <- AD_nonlesion_direction[,c('ligand', 'receptor', 'cell1', 'cell2', 'disease_direction')]
AD_lesion_direction  <- AD_lesion_direction[,c('ligand', 'receptor', 'cell1', 'cell2', 'disease_direction')]
PSO_nonlesion_direction  <- PSO_nonlesion_direction[,c('ligand', 'receptor', 'cell1', 'cell2', 'disease_direction')]
PSO_lesion_direction  <- PSO_lesion_direction[,c('ligand', 'receptor', 'cell1', 'cell2', 'disease_direction')]
COPD_direction  <- COPD_direction[,c('ligand', 'receptor', 'cell1', 'cell2', 'disease_direction')]
IPF_direction  <- IPF_direction[,c('ligand', 'receptor', 'cell1', 'cell2', 'disease_direction')]
UC_uninflamed_direction  <- UC_uninflamed_direction[,c('ligand', 'receptor', 'cell1', 'cell2', 'disease_direction')]
UC_inflamed_direction  <- UC_inflamed_direction[,c('ligand', 'receptor', 'cell1', 'cell2', 'disease_direction')]
IgAN_direction  <- IgAN_direction [,c('ligand', 'receptor', 'cell1', 'cell2', 'disease_direction')]
SLE_direction  <- SLE_direction[,c('ligand', 'receptor', 'cell1', 'cell2', 'disease_direction')]

# rename disease_direction column to average_value 
AD_nonlesion_direction <- AD_nonlesion_direction %>% rename(average_value = disease_direction)
AD_lesion_direction <- AD_lesion_direction %>% rename(average_value = disease_direction)
PSO_nonlesion_direction <- PSO_nonlesion_direction %>% rename(average_value = disease_direction)
PSO_lesion_direction <- PSO_lesion_direction %>% rename(average_value = disease_direction)
COPD_direction <- COPD_direction %>% rename(average_value = disease_direction)
IPF_direction <- IPF_direction %>% rename(average_value = disease_direction)
UC_uninflamed_direction <- UC_uninflamed_direction %>% rename(average_value = disease_direction)
UC_inflamed_direction <- UC_inflamed_direction %>% rename(average_value = disease_direction)
IgAN_direction <- IgAN_direction %>% rename(average_value = disease_direction)
SLE_direction <- SLE_direction %>% rename(average_value = disease_direction)

# make a column called tissue and fill with AD_direction, IPF_direction, etc
AD_nonlesion_direction$tissue <- 'AD_nonlesion_direction'
AD_lesion_direction$tissue <- 'AD_lesion_direction'
PSO_nonlesion_direction$tissue <- 'PSO_nonlesion_direction'
PSO_lesion_direction$tissue <- 'PSO_lesion_direction'
COPD_direction$tissue <- 'COPD_direction'
IPF_direction$tissue <- 'IPF_direction'
UC_uninflamed_direction$tissue <- 'UC_uninflamed_direction'
UC_inflamed_direction$tissue <- 'UC_inflamed_direction'
IgAN_direction$tissue <- 'IgAN_direction'
SLE_direction$tissue <- 'SLE_direction'

# make list of dataframes
data_frames <- list(AD_nonlesion_direction = AD_nonlesion_direction, AD_lesion_direction = AD_lesion_direction,
                    PSO_nonlesion_direction = PSO_nonlesion_direction, PSO_lesion_direction = PSO_lesion_direction,
                    COPD_direction = COPD_direction, IPF_direction = IPF_direction, 
                    UC_uninflamed_direction = UC_uninflamed_direction, UC_inflamed_direction = UC_inflamed_direction, 
                    IgAN_direction = IgAN_direction, SLE_direction = SLE_direction)
# run the function
meta_disease_direction <- metaTissue(data_frames)
# save as csv to check out the entire file if desired 
filename <- paste(getwd(), "cellphonedb_chemokine_diseasedirection.csv", sep = "/")
write.csv(meta_disease_direction, file = filename, row.names = FALSE)

##### csv processing 
csv <- read.csv("cellphonedb_chemokine_diseasedirection.csv")

# remove any cell1/2 columns that have NonImmune
csv <- csv[csv$cell1 != "NonImmune",]
csv <- csv[csv$cell2 != "NonImmune",]
csv <- csv[csv$cell1 != "Unclassified",]
csv <- csv[csv$cell2 != "Unclassified",]

# write the csv to excel
write_xlsx(csv, "cellphonedb_chemokine_diseasedirection.xlsx")

# read in the excel file
excel <- read_excel('cellphonedb_chemokine_diseasedirection.xlsx')

# indicate excel file for workbook function
wb <- loadWorkbook("cellphonedb_chemokine_diseasedirection.xlsx")

# make needed columns per tissue and add its corresponding value
tissues <- c("AD_nonlesion_direction", "AD_lesion_direction", "PSO_nonlesion_direction", "PSO_lesion_direction",
             'COPD_direction', 'IPF_direction', 'UC_uninflamed_direction', 'UC_inflamed_direction',
             'IgAN_direction', 'SLE_direction')
excel[, tissues] <- NA

# put tissue1 value in the appropriate tissue column
for (tissue in tissues) {
  value <- excel$tissue1 == tissue
  excel[value, tissue] <- excel$value1[value]}

# put tissue2 value in the appropriate tissue column
for (tissue in tissues) {
  value <- !is.na(excel$tissue2) & excel$tissue2 == tissue
  excel[value, tissue] <- excel$value2[value]}

# put tissue3 value in the appropriate tissue column
for (tissue in tissues) {
  value <- !is.na(excel$tissue3) & excel$tissue3 == tissue
  excel[value, tissue] <- excel$value3[value]}

# put tissue4 value in the appropriate tissue column
for (tissue in tissues) {
  value <- !is.na(excel$tissue4) & excel$tissue4 == tissue
  excel[value, tissue] <- excel$value4[value]}

# put tissue5 value in the appropriate tissue column
for (tissue in tissues) {
  value <- !is.na(excel$tissue5) & excel$tissue5 == tissue
  excel[value, tissue] <- excel$value5[value]}

# put tissue6 value in the appropriate tissue column
for (tissue in tissues) {
  value <- !is.na(excel$tissue6) & excel$tissue6 == tissue
  excel[value, tissue] <- excel$value6[value]}

# put tissue7 value in the appropriate tissue column
for (tissue in tissues) {
  value <- !is.na(excel$tissue7) & excel$tissue7 == tissue
  excel[value, tissue] <- excel$value7[value]}

# put tissue8 value in the appropriate tissue column
for (tissue in tissues) {
  value <- !is.na(excel$tissue8) & excel$tissue8 == tissue
  excel[value, tissue] <- excel$value8[value]}

# put tissue9 value in the appropriate tissue column
for (tissue in tissues) {
  value <- !is.na(excel$tissue9) & excel$tissue9 == tissue
  excel[value, tissue] <- excel$value9[value]}

# put tissue10 value in the appropriate tissue column
for (tissue in tissues) {
  value <- !is.na(excel$tissue10) & excel$tissue10 == tissue
  excel[value, tissue] <- excel$value10[value]}

# filter rows where 'value2' is not NA (= commons)
commons <- excel[!is.na(excel$value2), ]

# make ADPSO_commons tab if tissue5 is NA and tissue1+2+3+4 say either AD/PSO lesional/nonlesional (or NA for tissue3/4)
ADPSO_commons <- commons[(is.na(commons$tissue5)) &
                           (commons$tissue1 %in% c("AD_nonlesion_direction", "AD_lesion_direction", "PSO_nonlesion_direction", "PSO_lesion_direction"))
                         & (commons$tissue2 %in% c("AD_nonlesion_direction", "AD_lesion_direction", "PSO_nonlesion_direction", "PSO_lesion_direction"))
                         & (commons$tissue3 %in% c("AD_nonlesion_direction", "AD_lesion_direction", "PSO_nonlesion_direction", "PSO_lesion_direction", NA)) 
                         & (commons$tissue4 %in% c("AD_nonlesion_direction", "AD_lesion_direction", "PSO_nonlesion_direction", "PSO_lesion_direction", NA))
                         , ]
addWorksheet(wb, "ADPSO_commons")
writeData(wb, "ADPSO_commons", ADPSO_commons)
saveWorkbook(wb, "cellphonedb_chemokine_diseasedirection.xlsx", overwrite = T)

# make COPDIPF_commons tab if tissue3 is NA and tissue1+2 say COPD/IPF
COPDIPF_commons <- commons[(is.na(commons$tissue3)) & (commons$tissue1 == "COPD_direction") & (commons$tissue2 == "IPF_direction"), ]
addWorksheet(wb, "COPDIPF_commons")
writeData(wb, "COPDIPF_commons", COPDIPF_commons)
saveWorkbook(wb, "cellphonedb_chemokine_diseasedirection.xlsx", overwrite = T)

# make SLEIGAN_commons tab if tissue3 is NA and tissue1+2 say IGAN/SLE
SLEIGAN_commons <- commons[(is.na(commons$tissue3)) & (commons$tissue1 == "IgAN_direction") & (commons$tissue2 == "SLE_direction"), ]
addWorksheet(wb, "SLEIGAN_commons")
writeData(wb, "SLEIGAN_commons", SLEIGAN_commons)
saveWorkbook(wb, "cellphonedb_chemokine_diseasedirection.xlsx", overwrite = T)

# make UC_commons tab if tissue3 is NA and tissue1+2 say UC_inflamed/noninflamed
UC_commons <- commons[(is.na(commons$tissue3)) & (commons$tissue1 == "UC_uninflamed_direction") & (commons$tissue2 == "UC_inflamed_direction"), ]
addWorksheet(wb, "UC_commons")
writeData(wb, "UC_commons", UC_commons)
saveWorkbook(wb, "cellphonedb_chemokine_diseasedirection.xlsx", overwrite = T)

# update commons
commons <- anti_join(commons, ADPSO_commons)
commons <- anti_join(commons, COPDIPF_commons)
commons <- anti_join(commons, SLEIGAN_commons)
commons <- anti_join(commons, UC_commons)

# add filtered data  to new sheet called "commons"
addWorksheet(wb, "commons")
writeData(wb, "commons", commons)
saveWorkbook(wb, "cellphonedb_chemokine_diseasedirection.xlsx", overwrite = T)

# grab all the rows where there is no value2, meaning it is tissue specific
specifics <- excel[is.na(excel$value2),]

# put tissue1 value in the appropriate tissue column
for (tissue in tissues) {
  value <- specifics$tissue1 == tissue
  specifics[value, tissue] <- specifics$value1[value]}

# filter rows where column [each disease] is not NA, which = disease_only
AD_nonlesion <- specifics[!is.na(specifics$AD_nonlesion_direction), ]
AD_lesion <- specifics[!is.na(specifics$AD_lesion_direction), ]
PSO_nonlesion <- specifics[!is.na(specifics$PSO_nonlesion_direction), ]
PSO_lesion <- specifics[!is.na(specifics$PSO_lesion_direction), ]
COPD <- specifics[!is.na(specifics$COPD_direction), ]
IPF <- specifics[!is.na(specifics$IPF_direction), ]
UC_uninflamed <- specifics[!is.na(specifics$UC_uninflamed_direction), ]
UC_inflamed <- specifics[!is.na(specifics$UC_inflamed_direction), ]
IGAN <- specifics[!is.na(specifics$IgAN_direction), ]
SLE <- specifics[!is.na(specifics$SLE_direction), ]

# add to new sheet and save updated excel file
addWorksheet(wb, "AD_nonlesion")
writeData(wb, "AD_nonlesion", AD_nonlesion)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_diseasedirection.xlsx", overwrite = T)

# add to new sheet and save updated excel file
addWorksheet(wb, "AD_lesion")
writeData(wb, "AD_lesion", AD_lesion)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_diseasedirection.xlsx", overwrite = T)

# add to new sheet and save updated excel file
addWorksheet(wb, "PSO_nonlesion")
writeData(wb, "PSO_nonlesion", PSO_nonlesion)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_diseasedirection.xlsx", overwrite = T)

# add to new sheet and save updated excel file
addWorksheet(wb, "PSO_lesion")
writeData(wb, "PSO_lesion", PSO_lesion)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_diseasedirection.xlsx", overwrite = T)

# add to new sheet and save updated excel file
addWorksheet(wb, "COPD")
writeData(wb, "COPD", COPD)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_diseasedirection.xlsx", overwrite = T)

# add to new sheet and save updated excel file
addWorksheet(wb, "IPF")
writeData(wb, "IPF", IPF)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_diseasedirection.xlsx", overwrite = T)

# add to new sheet and save updated excel file
addWorksheet(wb, "UC_uninflamed")
writeData(wb, "UC_uninflamed", UC_uninflamed)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_diseasedirection.xlsx", overwrite = T)

# add to new sheet and save updated excel file
addWorksheet(wb, "UC_inflamed")
writeData(wb, "UC_inflamed", UC_inflamed)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_diseasedirection.xlsx", overwrite = T)

# add to new sheet and save updated excel file
addWorksheet(wb, "IGAN")
writeData(wb, "IGAN", IGAN)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_diseasedirection.xlsx", overwrite = T)

# add to new sheet and save updated excel file
addWorksheet(wb, "SLE")
writeData(wb, "SLE", SLE)
# save updated excel file
saveWorkbook(wb, "cellphonedb_chemokine_diseasedirection.xlsx", overwrite = T)

######################## la fin ############################