#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 01:00:00
#SBATCH -J bam_sort_index
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ibra818@gmail.com
#SBATCH --output=bam_post_%j.out
#SBATCH --error=bam_post_%j.err

module load bioinfo-tools bwa/0.7.18 samtools/1.19

IN_DIR=/proj/uppmax2025-3-3/private/Efaecium_Project/alignment/RNA_alignment
OUT_DIR=${IN_DIR}/sorted

mkdir -p "$OUT_DIR"

for bam in "${IN_DIR}"/*.bam; do
  sample=$(basename "$bam" .bam)
  samtools sort -@ 4 -o "${OUT_DIR}/${sample}.sorted.bam" "$bam"
  samtools index    "${OUT_DIR}/${sample}.sorted.bam"
done

