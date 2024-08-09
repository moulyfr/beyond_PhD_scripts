#!/bin/bash

#conda install -c bioconda velocyto

# Define directories
bam_dir="path/to/bam_dir"  # Directory containing subdirectories for each sample
loom_dir="path/to/loom_dir"  # Directory to store output Loom files

# Create loom_dir if it does not exist
mkdir -p "$loom_dir"

# Loop over each sample directory
for sample_dir in "$bam_dir"/*; do
    if [ -d "$sample_dir" ]; then
        sample_name=$(basename "$sample_dir")
        
        # Assuming BAM files are in the "bam" subdirectory
        bam_file="$sample_dir/bam/possorted_genome_bam.bam"
        
        if [ -f "$bam_file" ]; then
            loom_file="$loom_dir/${sample_name}.loom"
            
            # Run Velocyto to convert BAM to Loom
            velocyto run -b "$bam_file" -o "$loom_file" --sample "$sample_name"
        else
            echo "BAM file not found in $sample_dir"
        fi
    fi
done
