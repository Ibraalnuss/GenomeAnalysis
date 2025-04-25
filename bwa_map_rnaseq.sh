#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4               # 4 threads total
#SBATCH -t 02:00:00
#SBATCH -J bwa_map_RNA
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ibra818@gmail.com
#SBATCH --output=logs/bwa_map_RNA_%j.out
#SBATCH --error=logs/bwa_map_RNA_%j.err

module load bioinfo-tools
module load bwa/0.7.18       # ← specify version
module load samtools/1.19    # ← specify version

####### VARIABLES (✏️ edit these!) #######
REF=/proj/uppmax2025-3-3/private/Efaecium_Project/assembly/canu/efaecium_assembly.contigs.fasta
RAW_BASE=/proj/uppmax2025-3-3/Genome_Analysis/1_Zhang_2017/transcriptomics_data
OUT_DIR=/proj/uppmax2025-3-3/private/Efaecium_Project/alignment/RNA_alignment
##########################################

set -euo pipefail

mkdir -p "${OUT_DIR}" logs

if [[ ! -f "${REF}.bwt" ]]; then
  bwa index "${REF}"
fi

for DIR in "${RAW_BASE}"/RNA-Seq_*; do
  for R1 in "${DIR}"/trim_paired_*_pass_1.fastq.gz; do
    SAMPLE=$(basename "${R1}" _pass_1.fastq.gz | sed 's/trim_paired_//')
    R2=${R1%_pass_1.fastq.gz}_pass_2.fastq.gz

    bwa mem -t 4 "${REF}" "${R1}" "${R2}" \
      | samtools view -b - \
      | samtools sort -@2 -o "${OUT_DIR}/${SAMPLE}.bam"

    samtools index "${OUT_DIR}/${SAMPLE}.bam"
    echo ">>> ${SAMPLE}.bam ready"
  done
done

