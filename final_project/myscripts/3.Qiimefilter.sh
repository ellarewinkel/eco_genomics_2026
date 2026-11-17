#!/bin/bash
#SBATCH --job-name=Qiimefilter
#SBATCH --account=project_2003522
#SBATCH --time=48:00:00
#SBATCH --mem=300G
#SBATCH --partition=small
#SBATCH --cpus-per-task=8
#SBATCH -e logs/prepro_err_%j.txt
#SBATCH -o logs/prepro_out_%j.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.sadler@uvm.edu

module load bio-apps
module load qiime2

cd /users/desadlet/project_2003522/Danny/HypoxiaMicrobiome/Qiime

MANIFEST=/users/desadlet/project_2003522/Danny/HypoxiaMicrobiome/Qiime/manifest.tsv

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $MANIFEST \
  --output-path demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux-summary.qzv

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trunc-len-f 220 \
  --p-trunc-len-r 200 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-n-threads 8 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza \
  --o-base-transition-stats base-transition-stats.qza

qiime feature-table summarize \
  --i-table table.qza \
  --o-feature-frequencies table-feature-frequencies.qza \
  --o-sample-frequencies table-sample-frequencies.qza \
  --o-summary table.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv