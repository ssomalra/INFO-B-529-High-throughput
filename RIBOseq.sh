#!/bin/bash
#SBATCH --mail-user=ssomalra@iu.edu
#SBATCH --nodes=2
#SBATCH -p gpu
#SBATCH --ntasks-per-node=5
#SBATCH --gpus-per-node=1
#SBATCH --time=1-23:59:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --job-name=riboseq
#SBATCH -o riboseq.out
#SBATCH -A r00270

module load python
module load sra-toolkit
module load fastqc
#pip install umi_tools
module load trimgalore
module load star
module load subread

fasterq-dump SRR26942054
fasterq-dump SRR26942053
fasterq-dump SRR26942052
fasterq-dump SRR26942051	

mv SRR26942054.fastq riboseq_WT_1.fastq
mv SRR26942053.fastq riboseq_WT_2.fastq
mv SRR26942052.fastq riboseq_KO_1.fastq
mv SRR26942051.fastq riboseq_KO_2.fastq

# fastqc
fastqc -o fastqc/ *.fastq

# umi extraction
mkdir preprocessing
mkdir preprocessing/umi_extract

samples=("KO_1" "KO_2" "WT_1" "WT_2")

for sample in "${samples[@]}"
do
    umi_tools extract \
        --bc-pattern=XXXXXXNNNNNNNNN \
        -I /N/scratch/ssomalra/HighThroughput_Project/ribo-seq/riboseq_${sample}.fastq \
        -S /N/scratch/ssomalra/HighThroughput_Project/ribo-seq/preprocessing/umi_extract/riboseq_${sample}_umi_extract.fastq \
        --log=${sample}_umi_extract.log
done

# trimming
mkdir preprocessing/trimming

for sample in "${samples[@]}"
do
    trim_galore --phred33 --length 15 --fastqc \
        -o /N/scratch/ssomalra/HighThroughput_Project/ribo-seq/preprocessing/trimming \
        /N/scratch/ssomalra/HighThroughput_Project/ribo-seq/preprocessing/umi_extract/riboseq_${sample}_umi_extract.fastq
done

# alignment
mkdir alignment
for sample in "${samples[@]}"
do
    STAR --genomeDir /N/scratch/ssomalra/HighThroughput_Project/rna-seq/STAR_index \
         --readFilesIn /N/scratch/ssomalra/HighThroughput_Project/ribo-seq/preprocessing/trimming/riboseq_${sample}_umi_extract_trimmed.fq \
         --outFileNamePrefix /N/scratch/ssomalra/HighThroughput_Project/ribo-seq/alignment/riboseq_${sample} \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN 8 \
         --seedSearchStartLmax 12 \
         --outFilterMismatchNoverLmax 0.1 \
         --outFilterMultimapNmax 5
done

# deduplicate UMI
for sample in "${samples[@]}"
do
    umi_tools dedup --umi-separator="_" \
        --stdin=/N/scratch/ssomalra/HighThroughput_Project/ribo-seq/alignment/riboseq_${sample}Aligned.sortedByCoord.out.bam \
        --stdout=/N/scratch/ssomalra/HighThroughput_Project/ribo-seq/riboseq_${sample}Aligned_dedup.sortedByCoord.out.bam \
        --log=dedup_${sample}.log
done

# transcript quantification
featureCounts -a /N/scratch/ssomalra/HighThroughput_Project/references/GCF_000001405.40_GRCh38.p14_genomic.gtf -T 12 -o featureCounts_output.txt alignment/*.bam
