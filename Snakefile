import os
import csv

configfile: "config.yaml"

# ----------------------------
# Read samples from CSV
# ----------------------------
SAMPLES = {}
with open(config["samplesheet"]) as f:
    reader = csv.DictReader(f)
    for row in reader:
        SAMPLES[row["sample"]] = {
            "fastq_1": row["fastq_1"],
            "fastq_2": row["fastq_2"],
            "strandedness": row.get("strandedness", "reverse")
        }

SAMPLE_NAMES = list(SAMPLES.keys())

# ----------------------------
# Convenience variables
# ----------------------------
FASTA = config["fasta"]
GTF = config["gtf"]
INDEX_PREFIX = config["index_prefix"]
THREADS = config["threads"]

TRIMDIR = config["trimdir"]
QC_DIR = config["qc_dir"]
ALIGN_DIR = config["align_dir"]
LOG_DIR = config["log_dir"]
MULTIQC_DIR = config["multiqc_dir"]
HTSEQ_DIR = config["htseq_dir"]
ADAPTERS = config["adapters"]

INDEX_DIR = os.path.dirname(INDEX_PREFIX)

# ----------------------------
# Final rule
# ----------------------------
rule all:
    input:
        expand(f"{QC_DIR}" + "/{sample}_1_fastqc.html", sample=SAMPLE_NAMES),
        expand(f"{QC_DIR}" + "/{sample}_2_fastqc.html", sample=SAMPLE_NAMES),
        expand(f"{TRIMDIR}" + "/{sample}_1P.fastq.gz", sample=SAMPLE_NAMES),
        expand(f"{TRIMDIR}" + "/{sample}_2P.fastq.gz", sample=SAMPLE_NAMES),
        expand(f"{ALIGN_DIR}" + "/{sample}.sorted.bam", sample=SAMPLE_NAMES),
        expand(f"{ALIGN_DIR}" + "/{sample}.sorted.bam.bai", sample=SAMPLE_NAMES),
        expand(f"{HTSEQ_DIR}" + "/{sample}_counts.txt", sample=SAMPLE_NAMES),
        f"{MULTIQC_DIR}/multiqc_report.html"

# ----------------------------
# Download adapters if missing
# ----------------------------
rule download_adapters:
    output:
        ADAPTERS
    shell:
        """
        mkdir -p adapters
        wget -O {output} https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE.fa
        """

# ----------------------------
# FastQC on raw FASTQs
# ----------------------------
rule fastqc:
    input:
        r1=lambda wc: SAMPLES[wc.sample]["fastq_1"],
        r2=lambda wc: SAMPLES[wc.sample]["fastq_2"]
    output:
        html1=f"{QC_DIR}" + "/{sample}_1_fastqc.html",
        zip1=f"{QC_DIR}" + "/{sample}_1_fastqc.zip",
        html2=f"{QC_DIR}" + "/{sample}_2_fastqc.html",
        zip2=f"{QC_DIR}" + "/{sample}_2_fastqc.zip"
    threads: THREADS
    shell:
        """
        mkdir -p {QC_DIR}
        fastqc -t {threads} -o {QC_DIR} {input.r1} {input.r2}
        """

# ----------------------------
# Trimmomatic
# ----------------------------
rule trimmomatic:
    input:
        r1=lambda wc: SAMPLES[wc.sample]["fastq_1"],
        r2=lambda wc: SAMPLES[wc.sample]["fastq_2"],
        adapters=ADAPTERS
    output:
        p1=f"{TRIMDIR}" + "/{sample}_1P.fastq.gz",
        u1=f"{TRIMDIR}" + "/{sample}_1U.fastq.gz",
        p2=f"{TRIMDIR}" + "/{sample}_2P.fastq.gz",
        u2=f"{TRIMDIR}" + "/{sample}_2U.fastq.gz"
    log:
        f"{LOG_DIR}" + "/{sample}.trimmomatic.log"
    threads: THREADS
    shell:
        """
        mkdir -p {TRIMDIR} {LOG_DIR}
        trimmomatic PE -threads {threads} \
            {input.r1} {input.r2} \
            {output.p1} {output.u1} \
            {output.p2} {output.u2} \
            ILLUMINACLIP:{input.adapters}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50 \
            > {log} 2>&1
        """

# ----------------------------
# Extract splice sites
# ----------------------------
rule extract_splice_sites:
    input:
        GTF
    output:
        "hisat2_index/splicesites.txt"
    shell:
        """
        mkdir -p hisat2_index
        hisat2_extract_splice_sites.py {input} > {output}
        """

# ----------------------------
# Extract exons
# ----------------------------
rule extract_exons:
    input:
        GTF
    output:
        "hisat2_index/exons.txt"
    shell:
        """
        mkdir -p hisat2_index
        hisat2_extract_exons.py {input} > {output}
        """

# ----------------------------
# Build HISAT2 index
# ----------------------------
rule build_hisat2_index:
    input:
        fasta=FASTA,
        ss="hisat2_index/splicesites.txt",
        exon="hisat2_index/exons.txt"
    output:
        expand("hisat2_index/Homo_sapiens.GRCh38_direct_index.{i}.ht2", i=range(1, 9))
    log:
        f"{LOG_DIR}/hisat2_build.log"
    threads: THREADS
    shell:
        """
        mkdir -p {INDEX_DIR} {LOG_DIR}
        hisat2-build --ss {input.ss} --exon {input.exon} \
            {input.fasta} {INDEX_PREFIX} > {log} 2>&1
        """

# ----------------------------
# HISAT2 alignment -> unsorted BAM
# ----------------------------
rule hisat2_align:
    input:
        idx=expand("hisat2_index/Homo_sapiens.GRCh38_direct_index.{i}.ht2", i=range(1, 9)),
        r1=f"{TRIMDIR}" + "/{sample}_1P.fastq.gz",
        r2=f"{TRIMDIR}" + "/{sample}_2P.fastq.gz"
    output:
        bam=temp(f"{ALIGN_DIR}" + "/{sample}.bam"),
        summary=f"{ALIGN_DIR}" + "/{sample}.hisat2.summary.log"
    threads: THREADS
    params:
        rgid=lambda wc: wc.sample,
        rgsm=lambda wc: wc.sample,
        strandness=config["hisat2_rna_strandness"]
    shell:
        """
        mkdir -p {ALIGN_DIR}
        hisat2 -x {INDEX_PREFIX} \
            -1 {input.r1} \
            -2 {input.r2} \
            --rna-strandness {params.strandness} \
            --summary-file {output.summary} \
            --threads {threads} \
            --rg-id {params.rgid} --rg SM:{params.rgsm} \
            --no-mixed --no-discordant --sensitive -I 1 -X 1000 \
            | samtools view -@ {threads} -bS - > {output.bam}
        """

# ----------------------------
# Sort BAM
# ----------------------------
rule sort_bam:
    input:
        f"{ALIGN_DIR}" + "/{sample}.bam"
    output:
        f"{ALIGN_DIR}" + "/{sample}.sorted.bam"
    threads: THREADS
    shell:
        """
        samtools sort -@ {threads} -o {output} {input}
        """

# ----------------------------
# Index BAM
# ----------------------------
rule index_bam:
    input:
        f"{ALIGN_DIR}" + "/{sample}.sorted.bam"
    output:
        f"{ALIGN_DIR}" + "/{sample}.sorted.bam.bai"
    threads: THREADS
    shell:
        """
        samtools index {input}
        """

# ----------------------------
# HTSeq-count
# ----------------------------
rule htseq_count:
    input:
        bam=f"{ALIGN_DIR}" + "/{sample}.sorted.bam",
        bai=f"{ALIGN_DIR}" + "/{sample}.sorted.bam.bai",
        gtf=GTF
    output:
        f"{HTSEQ_DIR}" + "/{sample}_counts.txt"
    log:
        f"{LOG_DIR}" + "/{sample}.htseq.log"
    params:
        stranded=config["htseq_stranded"]
    shell:
        """
        mkdir -p {HTSEQ_DIR} {LOG_DIR}
        htseq-count -f bam -r pos -s {params.stranded} -t exon -i gene_id \
            {input.bam} {input.gtf} > {output} 2> {log}
        """

# ----------------------------
# MultiQC
# ----------------------------
rule multiqc:
    input:
        expand(f"{QC_DIR}" + "/{sample}_1_fastqc.zip", sample=SAMPLE_NAMES),
        expand(f"{QC_DIR}" + "/{sample}_2_fastqc.zip", sample=SAMPLE_NAMES),
        expand(f"{HTSEQ_DIR}" + "/{sample}_counts.txt", sample=SAMPLE_NAMES),
        expand(f"{ALIGN_DIR}" + "/{sample}.hisat2.summary.log", sample=SAMPLE_NAMES),
        expand(f"{LOG_DIR}" + "/{sample}.trimmomatic.log", sample=SAMPLE_NAMES)
    output:
        f"{MULTIQC_DIR}/multiqc_report.html"
    shell:
        """
        mkdir -p {MULTIQC_DIR}
        multiqc -o {MULTIQC_DIR} {QC_DIR} {ALIGN_DIR} {LOG_DIR} {HTSEQ_DIR}
        """
