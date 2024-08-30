1) Preprocessing with Cellranger mkfastq and count commands.
- If using feature barcodes (e.g. for Perturb-seq), column with barcodes needs to be called: id, name, read, pattern, sequence, feature_type, target_gene_id, target_gene_name
- This will allow for extracting cells of particular feature (e.g. gRNA_A for gene_A) during down-stream analysis

2) Processing with Seurat pipeline (for the most part).
