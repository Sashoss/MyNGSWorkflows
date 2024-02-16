#!/bin/sh
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name rnaseq_analysis
#SBATCH --time=2-00:00:00

module load fastqc

fastqc -o fastqc_out R1.fastq.gz
fastqc -o fastqc_out R2.fastq.gz

module purge
module load Bowtie2
module load SAMtools

mkdir aligned_bams

hisat2 -x Reference/GRCH38 -p 16 --no-discordant -1 R1.fastq.gz -2 R2.fastq.gz  \
    | samtools view -@ 16 -b -F 4 - \
    | samtools fixmate -@ 16 -O SAM -m - - \
    | samtools sort -@ 16 -O BAM \
    | samtools markdup -@ 16 -r - - -O SAM \
    | samtools view -b - \
    > aligned_bams/aligned.bam

samtools index -@ 16 aligned_bams/aligned.bam

java -jar picard.jar MarkDuplicates --INPUT aligned_bams/aligned.bam --OUTPUT aligned_bams/duplicates_removed.bam --METRICS_FILE metrics.out --REMOVE_DUPLICATES true
samtools index aligned_bams/duplicates_removed.bam

featureCounts -T 12 -O -t gene_id -g ID -a data/gencode.v38.annotation.gtf/gencode.v38.annotation.gtf -o counts.txt aligned_bams/duplicates_removed.bam
