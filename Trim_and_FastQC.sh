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

### Modules ###
module load bioinfo-tools
module load FastQC
module load trimmomatic/0.39
module load MultiQC

### Configuration ###
# Paired‐end raw data directories:
RAW_DIRS=( \
  "/proj/uppmax2025-3-3/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RNA-Seq_Serum" \
  "/proj/uppmax2025-3-3/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RNA-Seq_BH" \
)

# Output root directory:
OUT_ROOT="/proj/uppmax2025-3-3/private/Efaecium_Project/analysis/qc/trimmed"
mkdir -p "$OUT_ROOT"/{raw_qc,trimmed,trimmed_qc,multiqc}

# Adapter file:
ADAPTERS="/proj/uppmax2025-3-3/private/FastQC/adapters/TruSeq3-PE.fa"

# Sample IDs:
SAMPLES=(ERR1797969 ERR1797970 ERR1797971 ERR1797972 ERR1797973 ERR1797974)

### 1) FastQC on raw reads ###
echo ">>> FastQC on RAW reads"
for dir in "${RAW_DIRS[@]}"; do
  for sample in "${SAMPLES[@]}"; do
    fastqc -t 1 -o "$OUT_ROOT/raw_qc" \
      "$dir/${sample}_1.fastq.gz" \
      "$dir/${sample}_2.fastq.gz"
  done
done

### 2) Trimmomatic ###
echo ">>> Trimmomatic trimming"
for dir in "${RAW_DIRS[@]}"; do
  for sample in "${SAMPLES[@]}"; do
    trimmomatic PE -threads 2 -phred33 \
      "$dir/${sample}_1.fastq.gz" "$dir/${sample}_2.fastq.gz" \
      "$OUT_ROOT/trimmed/${sample}_1.paired.fastq.gz"  "$OUT_ROOT/trimmed/${sample}_1.unpaired.fastq.gz" \
      "$OUT_ROOT/trimmed/${sample}_2.paired.fastq.gz"  "$OUT_ROOT/trimmed/${sample}_2.unpaired.fastq.gz" \
      ILLUMINACLIP:"$ADAPTERS":2:30:10 \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  done
done

### 3) FastQC on trimmed pairs ###
echo ">>> FastQC on TRIMMED reads"
for sample in "${SAMPLES[@]}"; do
  fastqc -t 1 -o "$OUT_ROOT/trimmed_qc" \
    "$OUT_ROOT/trimmed/${sample}_1.paired.fastq.gz" \
    "$OUT_ROOT/trimmed/${sample}_2.paired.fastq.gz"
done

### 4) MultiQC ###
echo ">>> Generating MultiQC report"
cd "$OUT_ROOT"
multiqc raw_qc trimmed_qc -o multiqc

echo "All done! See: $OUT_ROOT/multiqc/multiqc_report.html"
