#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 00:30:00
#SBATCH -J mummer_synteny
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ibra818@gmail.com
#SBATCH --output=synteny_%j.out
#SBATCH --error=synteny_%j.err

module load bioinfo-tools MUMmer/4.0.0rc1

REF_FASTA=/proj/uppmax2025-3-3/private/Efaecium_Project/assembly/data/sequence.fasta
ASM_FASTA=/proj/uppmax2025-3-3/private/Efaecium_Project/assembly/canu/efaecium_assembly.contigs.fasta
WORKDIR=/proj/uppmax2025-3-3/private/Efaecium_Project/analysis/synteny

mkdir -p "$WORKDIR"
cd "$WORKDIR"

# quick check for any alignments
nucmer --maxmatch -p quicktest "$REF_FASTA" "$ASM_FASTA"
show-coords -rcl quicktest.delta | head -n 10

# stringent one‑to‑one synteny
nucmer --maxmatch -l 100 -c 500 -p efaec_vs_ref "$REF_FASTA" "$ASM_FASTA"
delta-filter -1 efaec_vs_ref.delta > efaec_vs_ref.1delta
show-coords -rcl efaec_vs_ref.1delta > efaec_vs_ref.coords.tsv

# plot
mummerplot --png --layout -p synteny_plot efaec_vs_ref.1delta
