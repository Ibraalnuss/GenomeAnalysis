#!/bin/bash
#SBATCH -A uppmax2025-3-3   # your project account
#SBATCH -M snowy            # cluster
#SBATCH -p core             # single node
#SBATCH -n 1                # 1 core is enough for download
#SBATCH -t 00:10:00
#SBATCH -J get_ref
#SBATCH --output=get_ref_%j.out
#SBATCH --error =get_ref_%j.err

set -euo pipefail

# 1) Variables
REF_ACC=GCF_000007565.1_ASM756v1
FTP_BASE=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/565/${REF_ACC}
OUTDIR=/proj/uppmax2025-3-3/private/Efaecium_Project/assembly/data/reference

# 2) Make directory & cd
mkdir -p "${OUTDIR}"
cd       "${OUTDIR}"

# 3) Download the genomic FASTA
wget "${FTP_BASE}/${REF_ACC}_genomic.fna.gz"

# 4) Uncompress & rename
gunzip "${REF_ACC}_genomic.fna.gz"
mv   "${REF_ACC}_genomic.fna" Efaecium_ref.fasta

echo "Downloaded and prepared Efaecium_ref.fasta in ${OUTDIR}"

