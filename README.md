---
editor_options: 
  markdown: 
    wrap: sentence
---

# Molecular Modelling of Trachoma Transmission Dynamics

## Project Description

This project focuses on modeling the transmission dynamics of trachoma at the molecular level using genomic data.
The objectives are:

1.  **Identify genetic markers** and variants associated with transmission.

2.  **Analyze genomic data** to infer evolutionary and transmission patterns.

3.  **Simulate and model** molecular dynamics to understand trachoma transmission mechanisms.

### Significance

Trachoma remains a major public health concern, and understanding its transmission dynamics can inform better control strategies.
This project integrates genomic and computational approaches to provide new insights into the disease’s spread and evolution.

------------------------------------------------------------------------

## Repository Structure

```         
|-- data/                  # Input datasets (FASTQ files, reference genomes, etc.)
|-- scripts/               # Analysis and processing scripts
|-- results/               # Output files (e.g., BAM, VCF, reports, figures)
|-- docs/                  # Documentation and additional resources
|-- notebooks/             # Jupyter Notebooks for exploratory analysis
|-- requirements.txt       # Python/R dependencies
|-- README.md              # Project documentation
```

------------------------------------------------------------------------

## Requirements and Dependencies

The project requires the following tools, libraries, and databases:

### Tools and Software

-   **FastQC**: Quality control for raw reads.
-   **BWA**: Alignment to reference genomes.
-   **GATK**: Variant calling.
-   **Kraken2**: Taxonomic classification.
-   **Mash**: Identify the closest reference genome.

### Python Packages

-   numpy (v1.23+)
-   pandas (v1.5+)
-   biopython (v1.79+)
-   matplotlib (v3.6+)

### Databases

-   NCBI RefSeq or Ensembl for reference genomes.

### Installation

Install dependencies using Conda:

``` bash
conda env create -f environment.yml
conda activate trachoma_env
```

------------------------------------------------------------------------

## Data Sources

-   **Genomic Data**: Paired-end FASTQ files from trachoma samples. Samples were obtained from different countries including Gambia, Guinea and Tanzania and were adequate for transmission analysis as the they were designed for such studies.
-   **Reference Genomes**: Downloaded from NCBI RefSeq or Ensembl.
-   **Metadata**: Sample metadata including geographic and temporal information.

------------------------------------------------------------------------

## Workflow Overview

Whole Genome Sequencing (WGS) variant calling for *Chlamydia trachomatis* (the causative agent of trachoma) typically involves several key steps.
Below is a detailed guide tailored for *C. trachomatis*, which is a bacterium with a relatively small genome (\~1 Mb):

### 1. Data Preparation/ Preprocessing

-   Perform quality control on raw reads using FastQC.
-   Trim adapters and low-quality bases with Trimmomatic.

1.  **Obtain Raw Reads:**

    -   Get the raw FASTQ files (paired-end or single-end reads) from sequencing.

    -   Check the sequencing quality (Illumina, PacBio, Nanopore, etc.).

2.  **Quality Control (QC):**

    -   Use tools like **FastQC** or **MultiQC** to assess read quality (e.g., base quality scores, adapter contamination).

    -   Trim low-quality bases and adapter sequences using **Trimmomatic**, **Cutadapt, Fastp, Sickle, bbduk or Trim Galore**.

### 2. Closest Reference Identification

-   Use Kraken2 or Mash to identify the most appropriate reference genome for alignment.

<!-- -->

-   **Download the Reference Genome:**

    -   Obtain the reference genome of *C. trachomatis* from databases like [NCBI](https://www.ncbi.nlm.nih.gov/) or Ensembl Bacteria.

-   **Index the Genome:**

    -   Use tools like **bwa index**, **bowtie2-build**, or **samtools faidx** to prepare the genome for alignment.

### 3. Read Alignment

-   Align paired-end reads to the identified reference genome using BWA.

**Map Reads to the Reference Genome:**

-   Align the trimmed reads to the *C. trachomatis* reference genome using a high-performance aligner like **BWA-MEM** or **Bowtie2**.

``` bash
bwa mem ref_genome.fasta reads_R1.fastq reads_R2.fastq > aligned_reads.sam
```

**Convert and Sort BAM File:**

-   Convert SAM to BAM, sort, and index using **samtools**:

``` bash
samtools view -S -b aligned_reads.sam | samtools sort -o sorted_reads.bam
samtools index sorted_reads.bam
```

**Assess Alignment Quality:**

-   Use **Qualimap** or **samtools stats** to check alignment metrics (e.g., mapping rate, coverage depth).

------------------------------------------------------------------------

### 4. Variant Calling

-   Use GATK to call SNPs and indels.

1.  **Mark Duplicates:**

    -   Remove or mark PCR duplicates with **Picard Tools**:

    ``` bash
    picard MarkDuplicates I=sorted_reads.bam O=dedup_reads.bam M=metrics.txt
    ```

2.  **Call Variants:**

    -   Use a variant caller such as **GATK HaplotypeCaller**, **FreeBayes**, or **Bcftools**:

    ``` bash
    gatk HaplotypeCaller \
        -R ref_genome.fasta \
        -I dedup_reads.bam \
        -O raw_variants.vcf
    ```

3.  **Filter Variants:**

    -   Filter low-quality variants to retain high-confidence SNPs and indels:

    ``` bash
    gatk VariantFiltration \
        -R ref_genome.fasta \
        -V raw_variants.vcf \
        -O filtered_variants.vcf \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
        --filter-name "LowQuality"
    ```

------------------------------------------------------------------------

### **5. Downstream Analysis**

1.  **Annotate Variants:**

    -   Annotate the variants using tools like **SnpEff** or **ANNOVAR** to identify the impact on genes.

    ``` bash
    snpEff ann ref_genome_db filtered_variants.vcf > annotated_variants.vcf
    ```

2.  **Analyze Genomic Regions:**

    -   Use software like **BEDTools** to assess coverage and identify variants in specific genes.

3.  **Phylogenetic or Population Analysis (Optional):**

    -   Use variant data to study genetic diversity, construct phylogenetic trees (e.g., with **IQ-TREE** or **RAxML**), or analyze population structure.

------------------------------------------------------------------------

### **6. Quality Assessment and Troubleshooting**

-   Verify the following metrics:
    -   Coverage depth (e.g., \>30x for reliable variant calling).
    -   SNP/indel quality scores.
    -   Mapping rates (e.g., \>95% alignment to the reference genome).

------------------------------------------------------------------------

### **Tools Required**

| Step              | Recommended Tool                       |
|-------------------|----------------------------------------|
| Quality Control   | FastQC, MultiQC, Trimmomatic, Cutadapt |
| Alignment         | BWA, Bowtie2, SAMtools                 |
| Duplicate Marking | Picard                                 |
| Variant Calling   | GATK, FreeBayes, Bcftools              |
| Annotation        | SnpEff, ANNOVAR                        |
| Phylogenetics     | IQ-TREE, RAxML, MEGA                   |

------------------------------------------------------------------------

### **Considerations for Trachoma Data**

-   **Recombination Hotspots:** *C. trachomatis* exhibits recombination in certain genomic regions. Take these into account when interpreting results.
-   **Plasmid Sequences:** Be aware of plasmid-associated sequences, as these may require separate handling or analysis.
-   **Strain Variation:** Choose the most appropriate reference genome for the strain you're analyzing to minimize alignment biases.

### 5. Genetic Diversity Analysis

-   Compute nucleotide diversity, Tajima’s D, and other measures.

### 6. Transmission Modelling

-   Build phylogenetic trees and statistical models to infer transmission dynamics.

### 7. Simulation Studies

-   Use SLiM to simulate molecular evolution and transmission scenarios.

### 8. Visualization

-   Generate heatmaps, phylogenetic trees, and other plots.

------------------------------------------------------------------------

## Step-by-Step Instructions

### Format metadata

``` r
library(tidyverse)

df1 <- readxl::read_xlsx("data/metadata/SeneGamMetadata.xlsx", sheet = 1)
names(df1)[c(1,6)] <-  c("samples", "year")
df2 <- readxl::read_xlsx("data/metadata/SeneGamMetadata.xlsx", sheet = 2)
df2 <- df2 %>% distinct(samples, .keep_all = TRUE)

df <- df1 %>% full_join(., df2, by = "samples")

writexl::write_xlsx(df, "data/metadata/metadata.xlsx")

df <- read_tsv("data/metadata/Gambia_metadata.tsv")

df1 <- readxl::read_xlsx("data/metadata/metadata.xlsx")

df1 %>% 
   inner_join(., df, by = c("samples" = "SampleID")) %>% 
   write_tsv("data/metadata/metadata_to_use.tsv", col_names = TRUE)
```

### 
