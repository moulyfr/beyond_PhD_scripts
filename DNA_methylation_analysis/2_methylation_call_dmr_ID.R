# this script is for methylation calling and DMR ID (with methylkit) 
library(methylKit)

# use processBismarkAln func to process bismark alignment files to filter based on
  # min cov and base qual score, and categorize according to the group
    # assemnly = rn6 (rat genome)
    # save.context = CpG (methylation context)
    # mincov = min coverage (how many times it was read) req.d to include a cytosine site 
    # minqual = min base quality score required (based on signal intensity in the machine)
# Read the file and sample ID information from a CSV file (sample id, treatment
info <- read.csv("sample_info.csv")

# extract the relevant columns
sample_id <- info$sample_id
treatment <- info$treatment

# call the function
data_bismark <- processBismarkAln(
  sample.id = sample_ids,
  assembly = "rn6",
  save.context = "CpG",
  read.context = "CpG",
  mincov = 10,
  minqual = 20,
  treatment = treatment
)               
# filter again 
  # lo.count: only keep CpG sites with at least 10 reads
  # hi.perc = 99.9: filter out highest coverage CpG sites, as these may be outliers
data_filtered=filterByCoverage(data_bismark, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)

# combine meth info from all samples into single object
  # destrand = F --> whether to combine methylation data from both DNA strands
data_methyl=unite(data_filtered, destrand=FALSE)
nrow(data_methyl)

# combine strands to get strand-agnostic view
data_methyl_agnos=unite(data_filtered, destrand=TRUE)
row(data_methyl_agnos)

# dmr analysis              
data_dmr =calculateDiffMeth(data_methyl,num.cores=2)
nrow(data_dmr)

# extract differentially methylated CpG sites from the DM analysis based on:
  # DMCs with 5% meth diff with p value 0.05
data_dmr_5=getMethylDiff(data_dmr,difference=5,qvalue=0.05)
nrow(data_dmr_5)

# save
write.table(data_dmr_5, file = "data_methyl_5.txt", sep = "\t", row.names = FALSE)

## QC
# correlation plot to see how similar/diff samples are to each other (batch effects?)
pdf("correlation.pdf")
getCorrelation(data_methyl,plot=T)
dev.off()

# sample clustering to see how samples group together based on meth profile
pdf("cluster.pdf")  
clusterSamples(data_methyl,dist="correlation", method="ward",plot=TRUE)
dev.off()

# PCA to visualize overall structure, variance in data, clustering, outliers
pdf("PCA.pdf")
PCASamples(data_methyl)
dev.off()

################################# Annotation
# rat gene annotation
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("compEpiTools")
# annotation package for rat genes
BiocManager::install("org.Rn.eg.db")
# annotation package containing rat transcript annotations 
BiocManager::install("TxDb.Rnorvegicus.UCSC.rn6.refGene")
# rat genome sequence
BiocManager::install("BSgenome.Rnorvegicus.UCSC.rn6")

library(compEpiTools)
library(org.Rn.eg.db)
library(TxDb.Rnorvegicus.UCSC.rn6.refGene)
library(BSgenome.Rnorvegicus.UCSC.rn6)
rattxdb <- TxDb.Rnorvegicus.UCSC.rn6.refGene

# read in the differential methylation cites 
data_dmr_5_input <- read.delim(file = "data_methyl_5.txt", stringsAsFactors = FALSE)

# convert to GRanges object, which represents genomic ranges useful for genomic interval anal
dmr_gr <- makeGRangesFromDataFrame(data_dmr_5_input, keep.extra.columns=TRUE)

# annotate the genomic ranges
  # object: use midpoint of GRanges for annotation
  # txdb specifies gene annotations
  # EG2GS provides mapping from gene IDs to gene symbols for annotation
  # upstream/ds extends annotation to include 2000 bps up/down stream of regions of interest
data_annotat <- GRannotate(Object = GRmidpoint(dmr_gr), txdb = rattxdb, EG2GS = org.Rn.eg.db,
                                  upstream = 2000, downstream = 2000)

# convert the GRanges object to data frame          
data_annotat_df <- as.data.frame(data_annotat)
colnames(data_annotat_df)

# take out "location_tx_id" and "location_gene_id" columns
data_annotat_df <- data_annotat_df[,-c(14,15)]

# function to get gene hits
func_gene_uniq <- function(x){ # let x be an entry of a character vector
  # use strsplit to separate the entries of the character vector
  a <- unlist(strsplit(x, split = ";"))
  # get unique entry of a
  b <- unique(a)
  # collapse results of b into a character that separates any unique entries by a semicolon
  c <- paste(b, collapse = ";")
  return(c)
}  

# apply func_gene_uniq func to location_gene_symbol column 
annot_loc_gene_uniq <- as.character(unlist(sapply(data_annotat_df[,"location_gene_symbol"], func_gene_uniq)))                                  

# apply func again to remove duplication location entries
annot_loc_unique <- as.character(unlist(sapply(data_annotat_df[,"location"], func_gene_uniq)))                                  

# update the df
data_annotat_df[,"location_gene_symbol"] <- annot_loc_gene_uniq
data_annotat_df[,"location"] <- annot_loc_unique    

# replace rows with NA values in the column "location_gene_symbol" with blank
dmc_annot_df_2 <- replace(data_annotat_df[,"location_gene_symbol"],grep("NA",
	data_annotat_df[,"location_gene_symbol"],perl=TRUE),"")
data_annotat_df [,"location_gene_symbol"] <- dmc_annot_df_2           

# count number of features per chromosome 
table(data_annotat_df$seqnames)

# count gene symbols 
gene_name <- table(data_annotat_df$location_gene_symbol)
write.table(data_annotat_df, file = "data_methyl_5_annot.txt", row.names = FALSE, sep = "\t")

################################## la fin #####################################
