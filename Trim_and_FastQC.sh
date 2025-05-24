#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 02:00:00
#SBATCH -J trim_RNAseq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ibra818@gmail.com
#SBATCH --output=trim_and_qc.%j.out

set -euo pipefail

module load bioinfo-tools
module load FastQC/0.11.9
module load trimmomatic/0.39
module load MultiQC/1.22.2

RAW_DIRS=(
  "/proj/uppmax2025-3-3/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RNA-Seq_Serum"
  "/proj/uppmax2025-3-3/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RNA-Seq_BH"
)

OUT_ROOT="/proj/uppmax2025-3-3/private/Efaecium_Project/analysis/qc/trimmed"
mkdir -p "$OUT_ROOT"/paired "$OUT_ROOT"/unpaired "$OUT_ROOT"/trimmed_qc "$OUT_ROOT"/multiqc

ADAPTERS="$TRIMMOMATIC_ROOT/adapters/TruSeq3-PE.fa"

SAMPLES=(1797969 1797970 1797971 1797972 1797973 1797974)

echo ">>> Trimmomatic PE trimming"
for DIR in "${RAW_DIRS[@]}"; do
  for SAMPLE in "${SAMPLES[@]}"; do
    IN1="$DIR/trim_paired_ERR${SAMPLE}_pass_1.fastq.gz"
    IN2="$DIR/trim_paired_ERR${SAMPLE}_pass_2.fastq.gz"
    OUT1_P="$OUT_ROOT/paired/ERR${SAMPLE}_1.paired.fastq.gz"
    OUT1_U="$OUT_ROOT/unpaired/ERR${SAMPLE}_1.unpaired.fastq.gz"
    OUT2_P="$OUT_ROOT/paired/ERR${SAMPLE}_2.paired.fastq.gz"
    OUT2_U="$OUT_ROOT/unpaired/ERR${SAMPLE}_2.unpaired.fastq.gz"

    trimmomatic PE -threads 2 -phred33 \
      "$IN1" "$IN2" \
      "$OUT1_P" "$OUT1_U" \
      "$OUT2_P" "$OUT2_U" \
      ILLUMINACLIP:"$ADAPTERS":2:30:10 \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  done
done

echo ">>> FastQC on trimmed reads"
for SAMPLE in "${SAMPLES[@]}"; do
  fastqc -t 1 -o "$OUT_ROOT/trimmed_qc" \
    "$OUT_ROOT/paired/ERR${SAMPLE}_1.paired.fastq.gz" \
    "$OUT_ROOT/paired/ERR${SAMPLE}_2.paired.fastq.gz"
done

echo ">>> Generating MultiQC report"
cd "$OUT_ROOT"
multiqc trimmed_qc -o multiqc

echo "All done! See: $OUT_ROOT/multiqc/multiqc_report.html"
