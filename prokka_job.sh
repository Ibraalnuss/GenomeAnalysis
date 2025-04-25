#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4                  # Prokka will use 4 threads
#SBATCH -t 02:00:00
#SBATCH -J prokka_efaecium
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ibra818@gmail.com
#SBATCH --output=prokka_%j.out
#SBATCH --error=prokka_%j.err

module load bioinfo-tools
module load prokka/1.45-5b58020

# --------- paths ---------
ASSEMBLY=/proj/uppmax2025-3-3/private/Efaecium_Project/assembly/canu/efaecium_assembly.contigs.fasta
OUTDIR=/proj/uppmax2025-3-3/private/Efaecium_Project/annotation
# -------------------------

set -euo pipefail

prokka --outdir "${OUTDIR}" \
       --force \
       --prefix efaecium \
       --cpus 4 \
       "${ASSEMBLY}"

