#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@brief: script to convert BCL files to FASTQ using CellRanger mkfastq and run CellRanger count with  human reference genome

"""

import os
import subprocess

# Define paths and directories
CELLRANGER_PATH = '/cloud-home/mfr/pkg/cellranger-7.2.0/bin'      
REFERENCE_PATH = '/cloud-home/mfr/scrnaseq/admin/human_reference'           
FASTQ_OUTPUT_DIR = '/cloud-home/mfr/scrnaseq/raw_data/experiment_030508/output_fastq'            
COUNT_OUTPUT_DIR = '/cloud-home/mfr/scrnaseq/raw_data/experiment_030508/output_count'           
RUN_PATH = '/cloud-home/mfr/scrnaseq/raw_data/experiment_030508/folder_bcl’                   
SAMPLE_SHEET = '/cloud-home/mfr/scrnaseq/raw_data/experiment_030508/sample_sheet.csv'
FEATURE_BARCODES_FILE = '/cloud-home/mfr/scrnaseq/raw_data/experiment_030508/feature_barcodes.csv' 

def run_cellranger_mkfastq(cellranger_path, run_dir, sample_sheet, output_dir):
    “””convert BCL files to FASTQ using Cell Ranger mkfastq."""
    command = [
        os.path.join(cellranger_path, 'cellranger'),
        'mkfastq',
        '--run', run_dir,
        '--output-dir', output_dir,
        '--sample-sheet', sample_sheet
    ]
    print(“running Cell Ranger mkfastq with command:", ' '.join(command))
    subprocess.run(command, check=True)

def run_cellranger_count(cellranger_path, fastq_dir, reference_path, output_dir, feature_barcodes_file=None):
    """run Cell Ranger count."""
    command = [
        os.path.join(cellranger_path, 'cellranger'),
        'count',
        '--id=sample_count',
        '--transcriptome', reference_path,
        '--fastqs', fastq_dir,
        '--sample', 'sample_name'
        ‘—expected-cells’, ‘10000’
    ]
    if feature_barcodes_file and os.path.exists(feature_barcodes_file): 
        command += ['--feature-barcode', feature_barcodes_file] 

    print(“running Cell Ranger count with command:", ' '.join(command))
    subprocess.run(command, check=True)

def main():
    # run Cell Ranger mkfastq
    run_cellranger_mkfastq(CELLRANGER_PATH, RUN_PATH, SAMPLE_SHEET, FASTQ_OUTPUT_DIR)

    # run Cell Ranger count
    run_cellranger_count(CELLRANGER_PATH, FASTQ_OUTPUT_DIR, REFERENCE_PATH, COUNT_OUTPUT_DIR, FEATURE_BARCODES_FILE)

if __name__ == "__main__":
    main()
