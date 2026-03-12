
# HISAT2–HTSeq RNA-seq Alignment and Quantification Pipeline

## Overview

This workflow performs RNA-seq read alignment and gene-level quantification using a reproducible **Snakemake pipeline**. It integrates several widely used tools to process paired-end RNA-seq data from raw FASTQ files through to gene count matrices.

The pipeline performs the following steps:

1. Quality control of raw reads using **FastQC**
2. Adapter trimming and quality filtering using **Trimmomatic**
3. Genome index generation using **HISAT2**
4. Alignment of trimmed reads to the reference genome using **HISAT2**
5. BAM sorting and indexing using **SAMtools**
6. Gene-level quantification using **HTSeq-count**
7. Aggregated quality reports using **MultiQC**

The workflow is designed for reproducibility, scalability, and automatic handling of intermediate files.

---

# Pipeline Structure

The workflow is implemented using **Snakemake**, which defines each processing step as a rule and automatically resolves dependencies between steps.

Processing steps:

```
FASTQ files
   ↓
FastQC
   ↓
Trimmomatic
   ↓
HISAT2 index generation
   ↓
HISAT2 alignment
   ↓
SAMtools sort & index
   ↓
HTSeq-count
   ↓
MultiQC report
```

---

# Directory Structure

After running the pipeline, the project directory will contain the following structure:

```
hisat2_htseq_snakemake/
│
├── Snakefile
├── config.yaml
├── Samplesheet_HAoSMCs.csv
├── README.md
│
├── adapters/
├── qc_reports/
├── trimmed_data/
├── hisat2_index/
├── logs/
├── HAoSMC_HISAT2_Trimmomatic_Output_Ensembl112/
├── HTSeq_geneid_Results/
└── multiqc_reports/
```

---

# Input Files

## 1. Samplesheet

The pipeline uses a CSV samplesheet describing each RNA-seq sample.

Example format:

```
sample,fastq_1,fastq_2,strandedness
Sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,reverse
Sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz,reverse
```

Columns:

| Column       | Description                     |
| ------------ | ------------------------------- |
| sample       | Unique sample identifier        |
| fastq_1      | Path to forward read FASTQ file |
| fastq_2      | Path to reverse read FASTQ file |
| strandedness | Library strandedness            |

---

## 2. Reference Genome

The pipeline requires:

```
FASTA reference genome
GTF annotation file
```

Example:

```
/home/s2451842/GTF/Homo_sapiens.GRCh38.dna.primary_assembly.fa
/home/s2451842/GTF/Homo_sapiens.GRCh38.112.gtf
```

---

# Configuration File

Pipeline parameters are defined in `config.yaml`.

Example configuration:

```yaml
fasta: "/home/s2451842/GTF/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
gtf: "/home/s2451842/GTF/Homo_sapiens.GRCh38.112.gtf"

index_prefix: "hisat2_index/Homo_sapiens.GRCh38_direct_index"

samplesheet: "Samplesheet_HAoSMCs.csv"

threads: 8

adapters: "adapters/TruSeq3-PE.fa"

trimdir: "trimmed_data"
qc_dir: "qc_reports"
align_dir: "HAoSMC_HISAT2_Trimmomatic_Output_Ensembl112"
log_dir: "logs"
multiqc_dir: "multiqc_reports"
htseq_dir: "HTSeq_geneid_Results"

hisat2_rna_strandness: "RF"
htseq_stranded: "reverse"
```

---

# Software Requirements

The following software must be installed and accessible in the environment:

| Software    | Purpose                  |
| ----------- | ------------------------ |
| Snakemake   | Workflow management      |
| FastQC      | Raw read quality control |
| Trimmomatic | Adapter trimming         |
| HISAT2      | RNA-seq alignment        |
| SAMtools    | BAM processing           |
| HTSeq       | Gene quantification      |
| MultiQC     | Combined QC reporting    |

Example environment check:

```
snakemake --version
fastqc --version
hisat2 --version
samtools --version
htseq-count --version
multiqc --version
```

---

# Running the Pipeline

Navigate to the workflow directory:

```
cd /home/s2451842/Scripts/hisat2_htseq_snakemake
```

## Dry Run

Before executing the workflow, perform a dry run to verify dependencies:

```
snakemake -n -p
```

This command prints the steps without running them.

---

## Execute the Pipeline

Run the pipeline using multiple cores:

```
snakemake --cores 8 -p
```

Snakemake will automatically execute steps in the correct order and skip completed files.

---

# Output Files

## Trimmed Reads

```
trimmed_data/
   sample_1P.fastq.gz
   sample_2P.fastq.gz
```

---

## Alignment Files

```
HAoSMC_HISAT2_Trimmomatic_Output_Ensembl112/
   sample.sorted.bam
   sample.sorted.bam.bai
```

---

## Gene Counts

```
HTSeq_geneid_Results/
   sample_counts.txt
```

Each file contains gene-level read counts derived from the GTF annotation.

---

## Quality Control Reports

FastQC reports:

```
qc_reports/
```

Combined QC summary:

```
multiqc_reports/multiqc_report.html
```

---

# Logs

Pipeline logs are stored in:

```
logs/
```

Including:

```
sample.trimmomatic.log
sample.htseq.log
hisat2_build.log
```

---

# Notes

• The workflow automatically skips completed steps.
• Intermediate files are reused if the pipeline is restarted.
• HISAT2 index generation is performed only once.

---

# Author

Maryam Uhabila Usman
RNA-seq Differential Expression Analysis Workflow
University of Edinburgh
