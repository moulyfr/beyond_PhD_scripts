# scATAC-Seq *(single cell assay for transposase-accessible chromatin with seq)*

- Isolated nuclei treated with a transposes enzyme, usually Tn5, which cuts DNA and inserts seq adapters into regions of open/accessible chromatin (not bound by histones)
- DNA fragments with adapters are PCR-amp’d for sequencing, which will then inform genome regions that are accessible (=peaks)
- Open regions correlate with active regulatory elements like promoters, enhancer and TF binding sites
- If scATAC-seq is conducted simultaneous to scRNA-seq, then cell barcodes can be used to unify the cell annotations 

## Preprocessing Shell script:
- *With cellranger-atac mkfastq & cellranger-atac count*
- demultiplex
- adaptor trimming
- alignment

## Processing R script:
- read in raw data (h5), metadata, fragments file, fragments file index
- peak annotation
- QC (nucleosome signal, transcriptional start site enrichment score, pct of fragments in peaks, ratio reads in genomic blacklist regions)
  - tot no fragments in peaks = measure of cellular seq depth/complexity
  - fraction of fragments in peaks: cells with <15-20% fragments in peaks are low quality/technical artifacts
  - ratio reads in genomic blacklist regions: cells with many reads mapping to regions associated with artificial signals
- matrix formation
- batch correction
- dimension reduction / visualization / clustering
- cell ID annotation

## Down-stream script:
- differentially accessible regions


