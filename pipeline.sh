#!/bin/bash
# RNA-seq variant calling pipeline v1.1
set -e

# ====== 目录配置 ======
project_root=$(pwd)
data_dir="$project_root/data"
raw_dir="$project_root/raw"
report_dir="$project_root/report"
vcf_dir="$project_root/vcf"
ref_dir="$project_root/ref"

mkdir -p {$data_dir,$raw_dir,$report_dir,$vcf_dir}

# ====== 环境检查 ======
if ! command -v hisat2 &> /dev/null; then
    echo "错误：hisat2未安装！请先配置conda环境"
    exit 1
fi

# ====== 数据下载 ======
cd $raw_dir
prefetch SRR5251643 SRR5251644
fasterq-dump SRR5251643.sra
fasterq-dump SRR5251644.sra
cd $project_root

# ====== 质控 ======
fastp -i $raw_dir/SRR5251643.fastq -o $data_dir/SRR5251643.clean.fastq.gz \
      -h $report_dir/SRR5251643_fastp.html -j $report_dir/SRR5251643_fastp.json \
      -q 20 -u 30

fastp -i $raw_dir/SRR5251644.fastq -o $data_dir/SRR5251644.clean.fastq.gz \
      -h $report_dir/SRR5251644_fastp.html -j $report_dir/SRR5251644_fastp.json \
      -q 20 -u 30

# ====== 比对 ======
hisat2 -p 8 \
       -x $ref_dir/r64/genome \
       -U $data_dir/SRR5251643.clean.fastq.gz \
       -S $data_dir/SRR5251643_aligned.sam

hisat2 -p 8 \
       -x $ref_dir/r64/genome \
       -U $data_dir/SRR5251644.clean.fastq.gz \
       -S $data_dir/SRR5251644_aligned.sam

# ====== SAM转BAM ======
samtools view -b $data_dir/SRR5251643_aligned.sam | \
samtools sort -o $data_dir/SRR5251643.sorted.bam
samtools index $data_dir/SRR5251643.sorted.bam

samtools view -b $data_dir/SRR5251644_aligned.sam | \
samtools sort -o $data_dir/SRR5251644.sorted.bam
samtools index $data_dir/SRR5251644.sorted.bam

# ====== 标记重复 ======
gatk MarkDuplicates \
     -I $data_dir/SRR5251643.sorted.bam \
     -O $data_dir/SRR5251643.dedup.bam \
     -M $report_dir/SRR5251643.dedup.metrics.txt
samtools index $data_dir/SRR5251643.dedup.bam

gatk MarkDuplicates \
     -I $data_dir/SRR5251644.sorted.bam \
     -O $data_dir/SRR5251644.dedup.bam \
     -M $report_dir/SRR5251644.dedup.metrics.txt
samtools index $data_dir/SRR5251644.dedup.bam

# ====== 添加Read Group ======
gatk AddOrReplaceReadGroups \
     -I $data_dir/SRR5251643.dedup.bam \
     -O $data_dir/SRR5251643.final.bam \
     --RGID SRR5251643 \
     --RGSM SRR5251643 \
     --RGPL ILLUMINA \
     --RGLB lib1 \
     --RGPU unit1

gatk AddOrReplaceReadGroups \
     -I $data_dir/SRR5251644.dedup.bam \
     -O $data_dir/SRR5251644.final.bam \
     --RGID SRR5251644 \
     --RGSM SRR5251644 \
     --RGPL ILLUMINA \
     --RGLB lib1 \
     --RGPU unit1

# ====== 变异检测 ======
gatk HaplotypeCaller \
     -R $ref_dir/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
     -I $data_dir/SRR5251643.final.bam \
     -O $vcf_dir/SRR5251643.g.vcf.gz \
     -ERC GVCF \
     --dont-use-soft-clipped-bases \
     --native-pair-hmm-threads 4

gatk HaplotypeCaller \
     -R $ref_dir/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
     -I $data_dir/SRR5251644.final.bam \
     -O $vcf_dir/SRR5251644.g.vcf.gz \
     -ERC GVCF \
     --dont-use-soft-clipped-bases \
     --native-pair-hmm-threads 4

# ====== 变异过滤 ======
gatk VariantFiltration \
     -V $vcf_dir/SRR5251643.g.vcf.gz \
     -O $vcf_dir/SRR5251643.filtered.vcf \
     --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
     --filter-name "filtration"

gatk VariantFiltration \
     -V $vcf_dir/SRR5251644.g.vcf.gz \
     -O $vcf_dir/SRR5251644.filtered.vcf \
     --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
     --filter-name "filtration"

# ====== 变异注释 ======
java -Xmx8g -jar snpEff/snpEff.jar Saccharomyces_cerevisiae \
     $vcf_dir/SRR5251643.filtered.vcf > $vcf_dir/SRR5251643.annotated.vcf \
     -stats $vcf_dir/SRR5251643_snpEff_summary.html

java -Xmx8g -jar snpEff/snpEff.jar Saccharomyces_cerevisiae \
     $vcf_dir/SRR5251644.filtered.vcf > $vcf_dir/SRR5251644.annotated.vcf \
     -stats $vcf_dir/SRR5251644_snpEff_summary.html
