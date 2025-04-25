#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 00:30:00
#SBATCH -J blastn_synteny
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ibra818@gmail.com
#SBATCH --output=synteny_%j.out
#SBATCH --error=synteny_%j.err

module load bioinfo-tools
module load blast/2.15.0+

set -euo pipefail

ASSEMBLY=/proj/uppmax2025-3-3/private/Efaecium_Project/assembly/canu/efaecium_assembly.contigs.fasta
REF_FASTA=/proj/uppmax2025-3-3/private/Efaecium_Project/assembly/data/reference/Efaecium_ref.fasta
REFDB=/proj/uppmax2025-3-3/private/Efaecium_Project/assembly/data/reference/Efaecium_ref_db
WORKDIR=/proj/uppmax2025-3-3/private/Efaecium_Project/assembly/synteny

mkdir -p "${WORKDIR}"
cd       "${WORKDIR}"

if [ ! -f "${REFDB}.nin" ]; then
  makeblastdb \
    -in "${REF_FASTA}" \
    -dbtype nucl \
    -parse_seqids \
    -out "${REFDB}"
fi

blastn \
  -task blastn \
  -query "${ASSEMBLY}" \
  -db    "${REFDB}" \
  -evalue 1e-5 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -num_threads ${SLURM_CPUS_ON_NODE} \
  > efaecium_vs_ref_synteny.tsv

