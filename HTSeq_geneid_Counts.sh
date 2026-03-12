#!/bin/bash

# Define paths
GTF="/home/s2451842/GTF/Homo_sapiens.GRCh38.112.gtf"
OUTPUT_DIR="HTSeq_geneid_Results"
THREADS=8  # Adjust based on available resources

# Step 1: Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Step 2: Create a list of sorted BAM files in the current directory
BAM_FILES=$(ls ./*.sorted.bam)

# Step 3: Run HTSeq for each BAM file
for BAM_FILE in $BAM_FILES
do
    SAMPLE_NAME=$(basename "$BAM_FILE" .sorted.bam)  # Get sample name from BAM file name
    echo "Processing sample: $SAMPLE_NAME"
    
    # Run HTSeq-count for gene quantification with reverse strandedness
    echo "Running HTSeq-count for sample: $SAMPLE_NAME"
    htseq-count -f bam -r pos -s reverse -t exon -i gene_id "$BAM_FILE" "$GTF" > "$OUTPUT_DIR/${SAMPLE_NAME}_counts.txt"

    # Check if HTSeq output file was created successfully
    if [ ! -s "$OUTPUT_DIR/${SAMPLE_NAME}_counts.txt" ]; then
        echo "Error: HTSeq count file for $SAMPLE_NAME is empty. Check logs." | tee -a "$OUTPUT_DIR/htseq_errors.log"
    else
        echo "HTSeq-count completed for sample: $SAMPLE_NAME"
    fi
done

echo "HTSeq quantification completed for all samples."
