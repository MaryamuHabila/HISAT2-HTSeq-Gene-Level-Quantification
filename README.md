# HISAT2-HTSeq-Gene-Level-Quantification
This is a full working code to run HISAT2 for genome alignment followed by HTSeq for gene level quantification
Overall Summary of Workflow
Read sample names from Samplesheet_HAoSMCs.csv
Run FastQC on raw reads
Run Trimmomatic
Build HISAT2 index once
Align each sample with HISAT2
Sort and index BAMs
Run HTSeq-count on each sorted BAM
Run MultiQC at the end
