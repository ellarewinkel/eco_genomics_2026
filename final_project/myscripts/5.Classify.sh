#!/bin/bash
#SBATCH --job-name=QiimeTaxFilter
#SBATCH --account=project_2003522
#SBATCH --time=24:00:00
#SBATCH --mem=200G
#SBATCH --partition=small
#SBATCH --cpus-per-task=8
#SBATCH -e logs/taxfilter_err_%j.txt
#SBATCH -o logs/taxfilter_out_%j.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.sadler@uvm.edu

module load bio-apps
module load qiime2

cd /users/desadlet/project_2003522/Danny/HypoxiaMicrobiome/Qiime

# Taxonomy classification

qiime feature-classifier classify-sklearn \
  --i-classifier SILVA/silva-138-99-515-806-nb-classifier-local.qza \
  --i-reads rep-seqs.qza \
  --p-n-jobs 8 \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

#  Frequency filtering — remove rare ASVs

qiime feature-table filter-features \
  --i-table table.qza \
  --p-min-frequency 10 \
  --o-filtered-table table-f10.qza

# Taxonomy-based filtering

qiime taxa filter-table \
  --i-table table-f10.qza \
  --i-taxonomy taxonomy.qza \
  --p-include p__ \
  --p-exclude Mitochondria,Chloroplast,Eukaryota,Unassigned \
  --o-filtered-table table-f10-filtered.qza

qiime feature-table summarize \
   --i-table table-f10-filtered.qza \
   --o-feature-frequencies table-f10-filtered-frequencies.qza \
  --o-sample-frequencies table-f10-filtered-frequencies.qza \
  --o-summary table-f10-filtered.qzv


#Sample-level filtering — drop samples with too few reads to, 5000 too much?

# qiime feature-table filter-samples \
#   --i-table table-f10-filtered.qza \
#   --p-min-frequency 5000 \
#   --o-filtered-table table-final.qza

# # Re-summarize the final table 
# qiime feature-table summarize \
#   --i-table table-final.qza \
#   --o-feature-frequencies table-final-feature-frequencies.qza \
#   --o-sample-frequencies table-final-sample-frequencies.qza \
#   --o-summary table-final.qzv