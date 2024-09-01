**1) Preprocessing with Cellranger mkfastq and count commands.**
- For Pertub-seq
  - Use feature barcodes csv file, columns need to be called: id, name, read, pattern, sequence, feature_type, target_gene_id, target_gene_name
  - This will allow for extracting cells of particular feature (e.g. gRNA_A for gene_A) during down-stream analysis
- For reporter assays (e.g. GFP)
  - Cellranger mkref needs to be ran before running Cellranger mkfastq and count commands
  - Before running mkref, the reference genome fasta and GTF files needs to be updated with the reporter gene annotation
  - This will allow for extracting cells with the "GFP" feature during down-stream analysis

**2) Processing with Seurat pipeline (for the most part).**
