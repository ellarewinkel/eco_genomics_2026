#!/bin/bash
#SBATCH --job-name=Qiimefastp
#SBATCH --account=project_2003522
#SBATCH --time=12:00:00
#SBATCH --mem=200G
#SBATCH --partition=small
#SBATCH --cpus-per-task=8
#SBATCH -e prepro_err_%j.txt
#SBATCH -o prepro_out_%j.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.sadler@uvm.edu
 
module load bio-apps
module load fastp
 
IN=/users/desadlet/project_2003522/Danny/HypoxiaMicrobiome/all16S
OUT=${IN}/fastp_trimmed
LOGS=${OUT}/logs
 
mkdir -p "$OUT" "$LOGS"


##loop fastp through the sequences

for r1 in "$IN"/*_1.fastq.gz; do
    sample=$(basename "$r1" _1.fastq.gz)
    r2="$IN/${sample}_2.fastq.gz"
 
    if [[ ! -f "$r2" ]]; then
        echo "WARNING: missing R2 for sample $sample (expected $r2) — skipping" >&2
        continue
    fi
 
    echo "Running fastp on sample: $sample"
 
    fastp \
        -i "$r1" -I "$r2" \
        -o "$OUT/${sample}_1.trim.fastq.gz" -O "$OUT/${sample}_2.trim.fastq.gz" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 15 \
        --unqualified_percent_limit 40 \
        --length_required 50 \
        --n_base_limit 5 \
        --thread 8 \
        --json "$LOGS/${sample}.fastp.json" \
        --html "$LOGS/${sample}.fastp.html" \
        2> "$LOGS/${sample}.fastp.log"
done
 
echo "Done. Trimmed reads in $OUT, per-sample reports in $LOGS"
 