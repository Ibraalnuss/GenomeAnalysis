#!/bin/bash -l
#SBATCH -A uppmax2025-3-3            # Project/account
#SBATCH -M snowy                     # Cluster
#SBATCH -p core                      # Queue
#SBATCH -n 4                         # Cores
#SBATCH -t 00:30:00                  # Walltime
#SBATCH -J fastqc_rna                # Job name
#SBATCH --mail-type=END,FAIL         # Only mail on completion/failure
#SBATCH --mail-user=you@your.email   
#SBATCH --output=fastqc_rna.%j.out   # STDOUT (%j = jobID)

module purge
module load bioinfo-tools
module load FastQC

OUTPUT_DIR="/proj/uppmax2025-3-3/private/Efaecium_Project/analysis/qc/fastqc_rn$
mkdir -p "${OUTPUT_DIR}"

DIR1="/proj/uppmax2025-3-3/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RN$
DIR2="/proj/uppmax2025-3-3/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RN$

echo "[$(date)] FastQC on RNA‑Seq (Serum) …"
fastqc -t 4 -o "${OUTPUT_DIR}" \
  "${DIR1}"/trim_paired_*_pass_1.fastq.gz \
  "${DIR1}"/trim_paired_*_pass_2.fastq.gz \
  "${DIR1}"/trim_single_*_pass_1.fastq.gz \
  "${DIR1}"/trim_single_*_pass_2.fastq.gz

echo "[$(date)] FastQC on RNA‑Seq (BH) …"
fastqc -t 4 -o "${OUTPUT_DIR}" \
  "${DIR2}"/trim_paired_*_pass_1.fastq.gz \
  "${DIR2}"/trim_paired_*_pass_2.fastq.gz \
  "${DIR2}"/trim_single_*_pass_1.fastq.gz \
  "${DIR2}"/trim_single_*_pass_2.fastq.gz

echo "[$(date)] FastQC complete; reports in ${OUTPUT_DIR}"
