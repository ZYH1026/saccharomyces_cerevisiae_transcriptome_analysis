# saccharomyces_cerevisiae_transcriptome_analysis

> 基于RNA-Seq数据检测重组酵母菌株的遗传变异（SNP/Indel）

## 项目背景
利用李云城等(2019)的公开RNA-Seq数据（[SRR5251643](https://www.ncbi.nlm.nih.gov/sra/SRR5251643), [SRR5251644](https://www.ncbi.nlm.nih.gov/sra/SRR5251644)），分析携带不同木糖利用途径（XR-XDH vs XI）的工业酵母菌株的遗传变异。

## 快速开始

### 1. 环境配置
```bash
conda env create -f config/environment.yml
conda activate rna_seq_env

### 2.下载酵母基因组
wget -P ref/ ftp://ftp.ensembl.org/pub/release-98/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
gunzip ref/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz

### 3.建立索引
hisat2-build ref/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa ref/r64/genome
gatk CreateSequenceDictionary -R ref/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
samtools faidx ref/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
