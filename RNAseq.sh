#!/bin/bash
#SBATCH --mail-user=ssomalra@iu.edu
#SBATCH --nodes=2
#SBATCH -p gpu
#SBATCH --ntasks-per-node=5
#SBATCH --gpus-per-node=1
#SBATCH --time=1-23:59:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --job-name=rnaseq
#SBATCH -o rnaseq.out
#SBATCH -A r00270

module load python
module load sra-toolkit
module load fastqc
module load trimgalore
module load star
module load subread

# load data
fasterq-dump SRR26942041
fasterq-dump SRR26942042
fasterq-dump SRR26942038
fasterq-dump SRR26942039

mv SRR26942041.fastq rnaseq_WT_1.fastq
mv SRR26942042.fastq rnaseq_WT_2.fastq
mv SRR26942038.fastq rnaseq_KO_1.fastq
mv SRR26942039.fastq rnaseq_KO_2.fastq

# fastqc
fastqc -o fastqc/ *.fastq

# trimming
for sample in rnaseq_WT_1 rnaseq_WT_2 rnaseq_KO_1 rnaseq_KO_2
do
    trim_galore --phred33 --fastqc -o trimming ${sample}.fastq
done

# download reference genome
wget -P /N/scratch/ssomalra/HighThroughput_Project/references/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz

# assembly
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /N/scratch/ssomalra/HighThroughput_Project/ribo-seq/STAR_index --genomeFastaFiles /N/scratch/ssomalra/HighThroughput_Project/references/GCF_000001405.40_GRCh38.p14_genomic.fna --sjdbGTFfile /N/scratch/ssomalra/HighThroughput_Project/references/GCF_000001405.40_GRCh38.p14_genomic.gtf --sjdbOverhang 28

# alignment
mkdir alignment

for sample in rnaseq_KO_1 rnaseq_KO_2 rnaseq_WT_1 rnaseq_WT_2
do
    STAR --genomeDir /N/scratch/ssomalra/HighThroughput_Project/ribo-seq/STAR_index \
         --readFilesIn /N/scratch/ssomalra/HighThroughput_Project/rna-seq/trimming/${sample}_trimmed.fq \
         --outFileNamePrefix /N/scratch/ssomalra/HighThroughput_Project/rna-seq/alignment/${sample} \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN 8 \
         --outFilterMultimapNmax 1
done

# download gene annotation files
wget -P /N/scratch/ssomalra/HighThroughput_Project/references/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.gtf.gz

# transcript quantification
featureCounts -a /N/scratch/ssomalra/HighThroughput_Project/references/GCF_000001405.40_GRCh38.p14_genomic.gtf -T 12 -o featureCounts_output.txt alignment/*.bam


