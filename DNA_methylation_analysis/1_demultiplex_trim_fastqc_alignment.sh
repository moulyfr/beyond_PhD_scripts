#! /bin/bash
# this script performs demultiplexing (with Fumi tool), trimming and fasted (with trim galore), alignment (with Bismarck), per dedup (with fume tool)
# input = seq reads, output = bam files

# base paths
BASE_PATH="/cloud-home/mfr/DNAmeth/"

# paths to tools
TRIM_GALORE_PATH="${BASE_PATH}trim_galore_zip/trim_galore"
FUMITOOLS_PATH="${BASE_PATH}fumitools_zip/fumitools"
BISMARK_PATH="${BASE_PATH}bismark_zip/bismark"
REFERENCE_GENOME="${BASE_PATH}bismark_zip/bismark/reference/genome"

# input combined FastQ files for demultiplexing
COMBINED_R1="combined_R1.fastq.gz"
COMBINED_R2="combined_R2.fastq.gz"

# barcode file for demultiplexing
BARCODE_FILE="${BASE_PATH}barcode_file.txt"

# output directory
OUTPUT_DIR="${BASE_PATH}bam_outputs"

# create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# perform demultiplexing
echo "Starting demultiplexing..."
$FUMITOOLS_PATH/fumi demultiplex --barcode $BARCODE_FILE --input_r1 $COMBINED_R1 --input_r2 $COMBINED_R2 --output $OUTPUT_DIR

# perform trimming and FastQC
echo "Starting trimming and FastQC..."
for read1 in $OUTPUT_DIR/*_R1_*.fastq.gz; do
    # Extract the corresponding read2 file name
    read2=${read1/_R1_/_R2_}
    
    if [ -f "$read2" ]; then
        echo "Trimming files: $(basename $read1) and $(basename $read2)"
        $TRIM_GALORE_PATH --rrbs --non_directional --paired --length 15 --fastqc "$read1" "$read2" -o $OUTPUT_DIR
    else
        echo "Warning: No corresponding read2 file for $(basename $read1). Skipping trimming for this pair."
    fi
done

# perform alignment with Bismark
echo "Starting alignment with Bismark..."
for trimmed_read1 in $OUTPUT_DIR/*_R1_val_1.fq.gz; do
    trimmed_read2=${trimmed_read1/_R1_val_1.fq.gz/_R2_val_2.fq.gz}
    
    if [ -f "$trimmed_read2" ]; then
        echo "Aligning files: $(basename $trimmed_read1) and $(basename $trimmed_read2)"
        $BISMARK_PATH/bismark --genome $REFERENCE_GENOME --bowtie2 --paired "$trimmed_read1" "$trimmed_read2" -o $OUTPUT_DIR
    else
        echo "Warning: No corresponding read2 file for $(basename $trimmed_read1). Skipping alignment for this pair."
    fi
done

# perform PCR deduplication with FumiTools
echo "Starting PCR deduplication..."
for bam_file in $OUTPUT_DIR/*.bam; do
    echo "Deduplicating file: $(basename $bam_file)"
    $FUMITOOLS_PATH/fumi dedup --input $bam_file --output ${bam_file%.bam}_dedup.bam
done

echo "Processing complete!"
