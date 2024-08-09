library(tidyr)

# set your folder path
folder_path <- getwd()
# list all files in the folder
files <- list.files(folder_path, full.names = TRUE)

# extract folder names (everything before first dot)
folder_names <- gsub("^([^.]*)\\..*$", "\\1", basename(files))

# create unique folders for each name
unique_folder_names <- unique(folder_names)
for (name in unique_folder_names) {
  folder_name <- file.path(folder_path, name)
  # check if the folder already exists, if not, create it
  if (!file.exists(folder_name)) {
    dir.create(folder_name)
  }
}

# move files to their respective folders
for (i in seq_along(files)) {
  file_path <- files[i]
  file_name <- basename(file_path)
  folder_name <- folder_names[i]
  target_folder <- file.path(folder_path, folder_name)
  # move the file to the corresponding folder
  file.rename(file_path, file.path(target_folder, file_name))
}

# list folders in the directory
folders <- list.dirs(path = getwd(), full.names = TRUE, recursive = FALSE)

# rename each file in each folder
for (folder in folders) {
  files <- list.files(path = folder, full.names = TRUE)
  # loop through each file in the folder
  for (file in files) {
    # extract the filename and the extension
    file_parts <- strsplit(basename(file), "_")[[1]]
    # keep only what is after the last underscore
    new_name <- paste(tail(file_parts, 1), collapse = "_")
    # construct the new file path
    new_path <- file.path(dirname(file), new_name)
    # rename the file
    file.rename(file, new_path)
  }
}

########################### remove the GSM stuff from the folder names #########
# get a list of directories in the current directory
folders <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)

# rename each folder
for (folder in folders) {
  # extract the current folder name
  folder_name <- basename(folder)
  # remove everything before the first underscore, including the underscore
  new_folder_name <- gsub("^[^_]*_", "", folder_name)
  # construct the new folder path
  new_folder_path <- file.path(dirname(folder), new_folder_name)
  # rename the folder
  file.rename(folder, new_folder_path)
}

###################################################### make the csv file #######
# get the list of folders in the current working directory
folders <- list.dirs(path = getwd(), full.names = FALSE, recursive = FALSE)

# create a data frame with a column named "sample"
data <- data.frame(sample = folders)

# write the data frame to a CSV file
write.csv(data, file = "metadeta.csv", row.names = FALSE)

#################################### convert all csv.gz to txt.gz ##############
# list all folders in the current directory
folders <- list.dirs()
# iterate through each folder
for (folder in folders) {
  # construct the full path for each folder 
  full_folder_path <- file.path(getwd(), folder)
  # check if the folder contains a csv.gz file
  csv_files <- list.files(path = folder, pattern = "\\.csv\\.gz$", full.names = TRUE)

  if (length(csv_files) == 1) {
    # read compressed CSV file dynamically
    input_file <- csv_files[1]
    data <- read.csv(gzfile(input_file))
    # use system command to read compressed CSV file and save it to a data frame
    system(paste("zcat", input_file), intern = TRUE) -> csv_content
    data <- read.table(text = csv_content, header = TRUE, sep = ",")
    # remove the column name X
    col_index <- which(names(data) == 'X')
    colnames(data)[col_index] <- ""
    # flip the dataframe
    reshaped_data <- gather(data, key = "Sample", value = "Value", -1)
    # extract the "ENSG" part from the first column and create a new column
    reshaped_data$ENSG <- gsub("\\..*", "", as.character(data[[1]]))
    # select relevant columns
    reshaped_data <- reshaped_data[c("ENSG", "Sample", "Value")]
    reshaped_data <- as.data.frame(reshaped_data)
    reshaped_data <- spread(reshaped_data, key = Sample, value = Value)
    # remove the word Sample
    names(reshaped_data) <- gsub("ENSG", "", names(reshaped_data))
    # dynamically set the output file name based on the original csv.gz file
    output_file <- file.path(full_folder_path, paste0(folder, ".txt.gz"))
    # write the data to a tab-separated text file
    write.table(reshaped_data, file = gzfile(output_file), sep = "\t", row.names = FALSE, quote = FALSE)
    # print a message indicating the processing is done for this folder
    cat("Processed folder:", folder, "\n")
  }
}

############# remove the csv.gz files from all the folders ####################
# list all folders in the current directory
folders <- list.dirs()

# iterate through each folder
for (folder in folders) {
  # list all files with a .csv extension in the current folder
  csv_files <- list.files(path = folder, pattern = "\\.csv", full.names = TRUE)
  
  # remove the files
  file.remove(csv_files)
  
  # print a message indicating the removal
  cat("Removed .csv files from folder:", folder, "\n")
}


################################################################################
########### extra code: to convert csv.gz to txt.gz one sample at a time
################################################################################
output_file <-'output-4.txt'
# read compressed CSV file
data <- read.csv(gzfile("NS-AR001.csv.gz"))
# use system command to read compressed CSV file and save it to a data frame
system(paste("zcat", input_file), intern = TRUE) -> csv_content
data <- read.table(text = csv_content, header = TRUE, sep = ",")

# checkout the contents
print(data)
colnames(data)
str(data)

# remove the column name X
col_index <- which(names(data) == 'X')
colnames(data)[col_index] <- ""

# flip the dataframe
library(tidyr)
reshaped_data <- gather(data, key = "Sample", value = "Value", -1)

# extract the "ENSG" part from the first column and create a new column
reshaped_data$ENSG <- gsub("\\..*", "", as.character(data[[1]]))

# select relevant columns
reshaped_data <- reshaped_data[c("ENSG", "Sample", "Value")]
reshaped_data <- as.data.frame(reshaped_data)
reshaped_data <- spread(reshaped_data, key = Sample, value = Value)

str(reshaped_data)
# remove the word Sample
names(reshaped_data) <- gsub("ENSG", "", names(reshaped_data))

# write the data to a tab-separated text file
write.table(reshaped_data, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

