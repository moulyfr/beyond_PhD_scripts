#!/bin/bash

# dirs………………………………………………………………………………………………………………………………………………………
# dir with fastq files      
FASTQ_DIR="/cloud-home/mfr/bulk_rnaseq/exp_mixed_brain/fastqs"  

# dir to store trimmed fastq files      
TRIMMED_DIR="/cloud-home/mfr/bulk_rnaseq/exp_mixed_brain/trimmed_fastqs" 

# dir with STAR index for human genome
STAR_INDEX="/cloud-home/mfr/pkg/trimmedstar_index/human" 

# dir for STAR-aligned BAMs
ALIGNMENT_DIR="/cloud-home/mfr/bulk_rnaseq/exp_mixed_brain/aligned_bams" 

 # dir for featureCounts output
FEATURECOUNTS_DIR="/cloud-home/mfr/bulk_rnaseq/exp_mixed_brain/featurecount_txts"

# dir for fastqc reports before trimming
FASTQC_PRE_TRIM_DIR=“/cloud-home/mfr/bulk_rnaseq/exp_mixed_brain/fastqc_before_trim" 

# dir for fastqc reports post trimming 
FASTQC_POST_TRIM_DIR="/cloud-home/mfr/bulk_rnaseq/exp_mixed_brain/fastqc_post_trim"  

# ………………………………………………………………………………………………………………………………………………………

# mk dirs if they do not exist
mkdir -p "$TRIMMED_DIR" "$ALIGNMENT_DIR" "$FEATURECOUNTS_DIR" "$FASTQC_PRE_TRIM_DIR" "$FASTQC_POST_TRIM_DIR"

# run fastQC before trimming
fastqc -o "$FASTQC_PRE_TRIM_DIR" "$FASTQ_DIR"/*.fastq.gz

# run trimmomatic 
for sample in $(ls "$FASTQ_DIR"/*_R1.fastq.gz | sed 's/_R1.fastq.gz//'); do
    base=$(basename "$sample")
    R1="${base}_R1.fastq.gz"
    R2="${base}_R2.fastq.gz"
    
    trimmomatic PE -phred33 \
        "$FASTQ_DIR/$R1" "$FASTQ_DIR/$R2" \
        "$TRIMMED_DIR/${base}_R1_paired.fastq.gz" "$TRIMMED_DIR/${base}_R1_unpaired.fastq.gz" \
        "$TRIMMED_DIR/${base}_R2_paired.fastq.gz" "$TRIMMED_DIR/${base}_R2_unpaired.fastq.gz" \
        ILLUMINACLIP:/path/to/adapter/seq/adapter.fa:2:30:10 \
        SLIDINGWINDOW:4:20 \
        MINLEN:36
done

# fastQC after trimming
fastqc -o "$FASTQC_POST_TRIM_DIR" "$TRIMMED_DIR"/*.fastq.gz

# STAR alignment 
for sample in $(ls "$TRIMMED_DIR"/*_R1_paired.fastq.gz | sed 's/_R1_paired.fastq.gz//'); do
    base=$(basename "$sample")
    R1="${base}_R1_paired.fastq.gz"
    R2="${base}_R2_paired.fastq.gz"
    
    STAR --runThreadN 8 \
         --genomeDir "$STAR_INDEX" \
         --readFilesIn "$TRIMMED_DIR/$R1" "$TRIMMED_DIR/$R2" \
         --readFilesCommand zcat \
         --outFileNamePrefix "$ALIGNMENT_DIR/$base" \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMstrandField intronMotif
done

# featureCounts for quant
featureCounts -a /path/to/annotation.gtf -o "$FEATURECOUNTS_DIR/counts.txt" "$ALIGNMENT_DIR"/*.bam
