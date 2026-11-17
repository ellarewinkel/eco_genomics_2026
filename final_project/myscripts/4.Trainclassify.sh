#!/bin/bash
#SBATCH --job-name=TrainClassifier
#SBATCH --account=project_2003522
#SBATCH --time=24:00:00
#SBATCH --mem=200G
#SBATCH --partition=small
#SBATCH --cpus-per-task=8
#SBATCH -e logs/trainclf_err_%j.txt
#SBATCH -o logs/trainclf_out_%j.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.sadler@uvm.edu
 
module load bio-apps
module load qiime2
 
cd /users/desadlet/project_2003522/Danny/HypoxiaMicrobiome/Qiime/SILVA
 
##needed to train own sequences in the latest version of Qiime on Roihu

# Import reference sequences + taxonomy 
qiime feature-classifier extract-reads \
  --i-sequences silva-138-2-nr99-seqs.qza \
  --p-f-primer GTGYCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACNVGGGTWTCTAAT \
  --p-n-jobs 8 \
  --o-reads silva-138-99-515-806-extracted.qza
 
# train the naive Bayes classifier on the extracted region
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-99-515-806-extracted.qza \
  --i-reference-taxonomy silva-138-2-nr99-tax.qza \
  --o-classifier silva-138-99-515-806-nb-classifier-local.qza
 



