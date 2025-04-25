#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 06:00:00
#SBATCH -J canu_efaecium
#SBATCH --output=canu_%j.out
#SBATCH --error=canu_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ibra818@gmail.com

module load bioinfo-tools
module load canu/2.2

canu \
  -p efaecium_assembly \
  -d /proj/uppmax2025-3-3/private/Efaecium_Project/assembly/canu \
  genomeSize=3m \
  -pacbio /proj/uppmax2025-3-3/Genome_Analysis/1_Zhang_2017/genomics_data/PacBio/*.fastq.gz \
  useGrid=false \
  maxThreads=4

