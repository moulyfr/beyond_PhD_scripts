#!/bin/bash

#conda install -c bioconda velocyto

# dirs
# bams per sample
bam_dir="/cloud-home/mfr/scrnaseq/ingested/ILD_Lung_GSE135893/bam"  
# output dir per loom
loom_dir="/cloud-home/mfr/scrnaseq/ingested/ILD_Lung_GSE135893/looms"  

# mk loom_dir if it does not exist
mkdir -p "$loom_dir"

# loop over each sample directory
for sample_dir in "$bam_dir"/*; do
    if [ -d "$sample_dir" ]; then
        sample_name=$(basename "$sample_dir")
        bam_file="$sample_dir/bam/possorted_genome_bam.bam"
        
        if [ -f "$bam_file" ]; then
            loom_file="$loom_dir/${sample_name}.loom"
            
            # run Velocyto to convert BAM to Loom
            velocyto run -b "$bam_file" -o "$loom_file" --sample "$sample_name"
        else
            echo "BAM file not found in $sample_dir"
        fi
    fi
done
