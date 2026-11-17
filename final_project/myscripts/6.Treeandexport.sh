#!/bin/bash
#SBATCH --job-name=QiimeTreeExport
#SBATCH --account=project_2003522
#SBATCH --time=12:00:00
#SBATCH --mem=150G
#SBATCH --partition=small
#SBATCH --cpus-per-task=8
#SBATCH -e logs/treeexport_err_%j.txt
#SBATCH -o logs/treeexport_out_%j.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.sadler@uvm.edu

module load bio-apps
module load qiime2

cd /users/desadlet/project_2003522/Danny/HypoxiaMicrobiome/Qiime

# 1. Filter rep-seqs down to only the ASVs still present in table-final.qza
#    (mitochondria/chloroplast/Eukaryota/low-frequency ASVs were dropped
#    from the table earlier but rep-seqs.qza still has all of them)
qiime feature-table filter-seqs \
  --i-data rep-seqs.qza \
  --i-table table-f10-filtered.qza \
  --o-filtered-data rep-seqs-f10-filtered.qza

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-f10-filtered.qza \
  --o-visualization rep-seqs-f10-filtered.qzv

# Build the phylogenetic tree from the filtered sequences
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-f10-filtered.qza \
  --p-n-threads 8 \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

# Export everything for phyloseq 
mkdir -p exported
qiime tools export --input-path table-f10-filtered.qza --output-path exported/
qiime tools export --input-path taxonomy.qza --output-path exported/
qiime tools export --input-path rooted-tree.qza --output-path exported/
qiime tools export --input-path rep-seqs-f10-filtered.qza --output-path exported/
biom convert -i exported/feature-table.biom -o exported/feature-table.tsv --to-tsv

biom summarize-table \
  -i exported/feature-table.biom \
  -o exported/table-summary.txt

