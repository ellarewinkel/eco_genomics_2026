#!/bin/bash
#SBATCH --job-name=admixture
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=./logs/%x.%j.out
#SBATCH --error=./logs/%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=erewinke@uvm.edu

cd /gpfs1/home/e/r/erewinke/projects/eco_genomics_2026/population_genomics/myscripts

OUTDIR=/gpfs1/home/e/r/erewinke/projects/eco_genomics_2026/population_genomics/myresults

echo $OUTDIR

VCF=Lottia_gigantea_pl_flt3.renamed.vcf.gz

echo $VCF

module load plink

DATA_DIR=/gpfs1/home/e/r/erewinke/projects/eco_genomics_2026/population_genomics/myresults

VCF=$DATA_DIR/Lottia_gigantea_pl_flt3.renamed.vcf.gz

cd $DATA_DIR

plink --vcf $VCF --double-id --allow-extra-chr \
    --set-missing-var-ids @:# \
    --indep-pairwise 50 10 0.2 \
    --out $OUTDIR/Lg_LD_0.2_alldata

plink --vcf $VCF --double-id --allow-extra-chr \
    --set-missing-var-ids @:# \
    --extract $OUTDIR/Lg_LD_0.2_alldata.prune.in \
    --make-bed --pca --out $OUTDIR/Lg_LD_0.2_alldata

echo "my script is done!"
