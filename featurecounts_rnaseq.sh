#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 01:00:00
#SBATCH -J featureCounts
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ibra818@gmail.com
#SBATCH --output=featurecounts_%j.out
#SBATCH --error=featurecounts_%j.err

module load bioinfo-tools subread/2.0.3

BAM_DIR=/proj/uppmax2025-3-3/private/Efaecium_Project/alignment/RNA_alignment/sorted
GFF=/proj/uppmax2025-3-3/private/Efaecium_Project/annotation/efaecium.gff
OUT_DIR=/proj/uppmax2025-3-3/private/Efaecium_Project/alignment/RNA_alignment/counts

mkdir -p "$OUT_DIR"

featureCounts -T 4 -p -B \
  -F GFF \
  -t CDS \
  -g locus_tag \
  -a "$GFF" \
  -o "$OUT_DIR/efaecium_featureCounts.txt" \
  "$BAM_DIR"/*.sorted.bam

