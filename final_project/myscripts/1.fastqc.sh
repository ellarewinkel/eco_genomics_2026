#!/bin/bash
#SBATCH --job-name=QiimeQC
#SBATCH --account=project_2003522
#SBATCH --time=48:00:00
#SBATCH --mem=300G
#SBATCH --partition=small
#SBATCH --cpus-per-task=16
#SBATCH -e prepro_err_%j.txt 
#SBATCH -o prepro_out_%j.txt 
#SBATCH --mail-type=END
#SBATCH --mail-user daniel.sadler@uvm.edu

module load bio-apps
module load fastqc
module load multiqc

cd /users/desadlet/project_2003522/Danny/HypoxiaMicrobiome/all16S
OUT=/users/desadlet/project_2003522/Danny/HypoxiaMicrobiome/alldataHyp/Qiime/fastqc


####fastqc to check parameters for fastp

fastqc *_1.fastq.gz *_2.fastq.gz \
    --noextract \
    --outdir $OUT \
    -threads $SLURM_CPUS_PER_TASK

cd $OUT
multiqc fastqc -o fastqc