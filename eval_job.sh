#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 02:00:00
#SBATCH -J assembly_eval_efa
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ibra818@gmail.com
#SBATCH --output=assembly_eval_%j.out
#SBATCH --error=assembly_eval_%j.err

module load bioinfo-tools
module load quast/5.0.2
module load BUSCO/5.7.1        # loads Augustus & sets env‑vars

ASSEMBLY=/proj/uppmax2025-3-3/private/Efaecium_Project/assembly/canu/efaecium_assembly.contigs.fasta
WORKDIR=/proj/uppmax2025-3-3/private/Efaecium_Project/assembly/eval
LINEAGE=$BUSCO_LINEAGE_SETS/bacteria_odb10   # <-- valid lineage folder

set -euo pipefail
mkdir -p "$WORKDIR"
cd       "$WORKDIR"

echo "=== $(date)   QUAST ==="
quast.py -t 4 -o quast_out "$ASSEMBLY"

echo "=== $(date)   BUSCO ==="
# create writable Augustus config
source "$AUGUSTUS_CONFIG_COPY"

busco -i "$ASSEMBLY" \
      -l "$LINEAGE" \
      -o efaecium_busco \
      -m genome \
      -c 4 \
      --offline

echo "=== $(date)   Done ==="

