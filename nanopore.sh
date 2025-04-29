#!/bin/bash
#SBATCH --mail-user=ssomalra@iu.edu
#SBATCH --nodes=2
#SBATCH -p gpu
#SBATCH --ntasks-per-node=5
#SBATCH --gpus-per-node=1
#SBATCH --time=1-23:59:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --job-name=nanopore
#SBATCH -o nanopore.out
#SBATCH -A r00270

module load python
module load nanopolish
module load minimap
conda activate m6anet

# download data
wget https://sra-pub-src-2.s3.amazonaws.com/SRR9646141/HEK_WT_Fast5.tar.gz
tar -xzf HEK_WT_Fast5.tar.gz

# convert multi_fast5 to pod5
pod5 convert fast5 /N/scratch/ssomalra/HighThroughput_Project/nanopore/multi_read_fast5/*.fast5 --output /N/scratch/ssomalra/HighThroughput_Project/nanopore/pod5_files/ --one-to-one /N/scratch/ssomalra/HighThroughput_Project/nanopore/multi_read_fast5/

# dorado basecalling
mkdir input_files
/N/slate/ssomalra/Packages/dorado-0.7.0-linux-x64/bin/dorado basecaller --device cpu /N/slate/ssomalra/Packages/dorado_models/rna002_70bps_fast\@v3/ /N/scratch/ssomalra/HighThroughput_Project/nanopore/pod5_files/ --emit-fastq > /N/scratch/ssomalra/HighThroughput_Project/nanopore/input_files/dorado.fastq

# convert fastq > fasta
sed -n '1~4s/^@/>/p;2~4p' input_files/dorado.fastq > input_files/merged_reads.fasta

# index fasta file
nanopolish index -d /N/scratch/ssomalra/HighThroughput_Project/nanopore/multi_read_fast5 /N/scratch/ssomalra/HighThroughput_Project/nanopore/input_files/merged_reads.fasta

# alignment
minimap2 --secondary=no -a -x map-ont /N/scratch/ssomalra/HighThroughput_Project/references/GCF_000001405.40_GRCh38.p14_genomic.fna /N/scratch/ssomalra/HighThroughput_Project/nanopore/input_files/merged_reads.fasta > /N/scratch/ssomalra/HighThroughput_Project/nanopore/input_files/merged_reads.sam

# sam > bam
samtools view -S -b /N/scratch/ssomalra/HighThroughput_Project/nanopore/input_files/merged_reads.sam > /N/scratch/ssomalra/HighThroughput_Project/nanopore/input_files/merged_reads.bam

samtools sort /N/scratch/ssomalra/HighThroughput_Project/nanopore/input_files/merged_reads.bam -o /N/scratch/ssomalra/HighThroughput_Project/nanopore/input_files/merged_reads.sorted.bam

samtools index /N/scratch/ssomalra/HighThroughput_Project/nanopore/input_files/merged_reads.sorted.bam

# nanopolish eventalign
mkdir m6A
nanopolish eventalign --reads /N/scratch/ssomalra/HighThroughput_Project/nanopore/input_files/merged_reads.fasta --bam /N/scratch/ssomalra/HighThroughput_Project/nanopore/input_files/merged_reads.sorted.bam --genome /N/scratch/ssomalra/HighThroughput_Project/references/GCF_000001405.40_GRCh38.p14_genomic.fna --scale-events > /N/scratch/ssomalra/HighThroughput_Project/nanopore/m6A/eventalign_output.txt

# m6anet
m6anet dataprep --eventalign /N/scratch/ssomalra/HighThroughput_Project/nanopore/m6A/eventalign_output.txt --out_dir /N/scratch/ssomalra/HighThroughput_Project/nanopore/m6A/dataprep --n_processes 4

m6anet inference --input_dir /N/scratch/ssomalra/HighThroughput_Project/nanopore/m6A/dataprep --out_dir /N/scratch/ssomalra/HighThroughput_Project/nanopore/m6A/inference/ --n_processes 4 --num_iterations 1000








