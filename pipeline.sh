#!/bin/bash
# This script is used to run the pipeline for the project.
set -e  # Exit immediately if a command exits with a non-zero status.

project_root="your_project_root_folder"
cd $project_root

# project_root (the folder of the project should look like this)
#  |
#  +-- data
#  |
#  +-- raw
#  |
#  +-- report
#  |
#  +-- vcf
#  |
#  +-- ref
#  |  +-- Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
#  |  +-- r64
#  |     +-- genome.1.ht2
#  |     +-- ......



# Install softwares using conda
conda create -n project_env
conda activate project_env
conda install -c bioconda fastp
conda install -c bioconda hisat2
conda install -c bioconda samtools
conda install -c bioconda gatk4
conda install -c bioconda sra-tools

# Use prefetch to download the original sra files
cd ./raw
prefetch SRR5251643
prefetch SRR5251644

# Use fastq-dump to convert sra files to fastq files
fasterq-dump SRR5251644.sra
fasterq-dump SRR5251643.sra
cd ..

# Use fastp to trim the fastq files
fastp -i ./raw/SRR5251644.fastq \
      -o ./data/SRR5251644.clean.fastq.gz \
      -h ./report/SRR5251644_fastp_report.html \
      -j ./report/SRR5251644_fastp_report.json \
      -q 20 \
      -u 30

fastp -i ./raw/SRR5251643.fastq \
      -o ./data/SRR5251643.clean.fastq.gz \
      -h ./report/SRR5251643_fastp_report.html \
      -j ./report/SRR5251643_fastp_report.json \
      -q 20 \
      -u 30

# Use hisat2 to align the clean fastq files to the reference genome
hisat2 -p 8 \ 
-x ./ref/r64/genome \ 
-U ./data/SRR5251643.clean.fastq.gz \
-S ./data/SRR5251643_aligned.sam

hisat2 -p 8 \
-x ./ref/r64/genome \
-U ./data/SRR5251644.clean.fastq.gz \
-S ./data/SRR5251644_aligned.sam

# Use samtools to convert the sam files to bam files, creating index files
samtools view -b ./data/SRR5251644_aligned.sam | samtools sort -o ./data/SRR5251644.sorted.bam
samtools index ./data/SRR5251644.sorted.bam

samtools view -b ./data/SRR5251643_aligned.sam | samtools sort -o ./data/SRR5251643.sorted.bam
samtools index ./data/SRR5251643.sorted.bam

# Use GATK to mark duplicates and index the bam files
gatk MarkDuplicates \
     -I ./data/SRR5251644.sorted.bam \
     -O ./data/SRR5251644.dedup.bam \
     -M ./report/SRR5251644.dedup.metrics.txt
samtools index ./data/SRR5251644.dedup.bam

gatk MarkDuplicates \
     -I ./data/SRR5251643.sorted.bam \
     -O ./data/SRR5251643.dedup.bam \
     -M ./report/SRR5251643.dedup.metrics.txt
samtools index ./data/SRR5251643.dedup.bam

# Build the reference genome index using GATK and samtools
gatk CreateSequenceDictionary \
    -R ./ref/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
    -O ./ref/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.dict

samtools faidx ./ref/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa

# Add Read Group information
gatk AddOrReplaceReadGroups \
     -I ./data/SRR5251644.dedup.bam \
     -O ./data/SRR5251644.final.bam \
     --RGID SRR5251644 \
     --RGSM SRR5251644 \
     --RGPL ILLUMINA \
     --RGLB lib1 \
     --RGPU unit1

gatk AddOrReplaceReadGroups \
     -I ./data/SRR5251643.dedup.bam \
     -O ./data/SRR5251643.final.bam \
     --RGID SRR5251643 \
     --RGSM SRR5251643 \
     --RGPL ILLUMINA \
     --RGLB lib1 \
     --RGPU unit1

# Index the final bam files
samtools index ./data/SRR5251644.final.bam
samtools index ./data/SRR5251643.final.bam

# Use GATK HaplotypeCaller to call variants
gatk HaplotypeCaller \
     -R ./ref/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
     -I ./data/SRR5251644.final.bam \
     -O ./vcf/SRR5251644.vcf

gatk HaplotypeCaller \
     -R ./ref/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
     -I ./data/SRR5251643.final.bam \
     -O ./vcf/SRR5251643.vcf

# Filter the VCF files using GATK VariantFiltration
gatk VariantFiltration \
     -V ./vcf/SRR5251644.vcf \
     -O ./vcf/SRR5251644.filtered.vcf \
     --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
     --filter-name "filtration"

gatk VariantFiltration \
     -V ./vcf/SRR5251643.vcf \
     -O ./vcf/SRR5251643.filtered.vcf \
     --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
     --filter-name "filtration"

# Install SnpEff for annotation
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip

# Find and download the SnpEff database for Saccharomyces cerevisiae
java -jar snpEff/snpEff.jar databases | grep Saccharomyces_cerevisiae
java -jar snpEff/snpEff.jar download Saccharomyces_cerevisiae

# Run the annotation using SnpEff and generate statistical reports
java -Xmx8g -jar snpEff/snpEff.jar Saccharomyces_cerevisiae \
     ./vcf/SRR5251644.filtered.vcf > ./vcf/SRR5251644.annotated.vcf -stats ./vcf/SRR5251644_snpEff_summary.html

java -Xmx8g -jar snpEff/snpEff.jar Saccharomyces_cerevisiae \
     ./vcf/SRR5251643.filtered.vcf > ./vcf/SRR5251643.annotated.vcf -stats ./vcf/SRR5251643_snpEff_summary.html 