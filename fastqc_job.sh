#!/bin/bash -l
#SBATCH -A uppmax2025-3-3            # Project/Account name
#SBATCH -M snowy                     # Cluster name (if required)
#SBATCH -p core                      # Partition/queue
#SBATCH -n 2                         # Number of cores
#SBATCH -t 00:30:00                  # Time limit (HH:MM:SS)
#SBATCH -J fastqc_job                # Job name
#SBATCH --mail-type=ALL              # Mail notifications for job events
#SBATCH --mail-user=your_email@example.com  # Replace with your actual email
#SBATCH --output=fastqc_job.%j.out   # Output file name (%j will be replaced by the job id)

module load bioinfo-tools  
module load FastQC

FASTQC_EXEC="/proj/uppmax2025-3-3/private/FastQC/fastqc"
OUTPUT_DIR="/proj/uppmax2025-3-3/private/Efaecium_Project/analysis/qc/fastqc_raw"
ILLUMINA_DIR="/proj/uppmax2025-3-3/Genome_Analysis/1_Zhang_2017/genomics_data/Illumina"
PACBIO_DIR="/proj/uppmax2025-3-3/Genome_Analysis/1_Zhang_2017/genomics_data/PacBio"
NANOPORE_DIR="/proj/uppmax2025-3-3/Genome_Analysis/1_Zhang_2017/genomics_data/Nanopore"


echo "Starting FastQC on Illumina data..."
${FASTQC_EXEC} -o ${OUTPUT_DIR} ${ILLUMINA_DIR}/*.fq.gz

echo "Starting FastQC on PacBio data..."
${FASTQC_EXEC} -o ${OUTPUT_DIR} ${PACBIO_DIR}/*.fastq.gz

echo "Starting FastQC on Nanopore data (FASTA format; quality metrics will be limited)..."
${FASTQC_EXEC} -o ${OUTPUT_DIR} ${NANOPORE_DIR}/*.fasta.gz

echo "FastQC analysis complete."

