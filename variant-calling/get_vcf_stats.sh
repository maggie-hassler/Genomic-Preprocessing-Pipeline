#!/bin/sh
#SBATCH --job-name="chr1_stats"
#SBATCH --mail-user=mhassle1@asu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/slurm.%A_%a.out
#SBATCH -e logs/slurm.%A_%a.err
#SBATCH -n 1
#SBATCH --time=0-04:00:00
#SBATCH --partition=htc
#SBATCH --mem=16G

# modified on 5/28/25 for WGS chr1 

# exit on silent errors
set -euo pipefail

# load bcftools
module load bcftools-1.14-gcc-11.2.0

# define paths
INPUT_DIR="/path/to/vcfs"
OUTPUT_DIR="/path/to/output_directory"

# make sure output directory exists
mkdir -p "$OUTPUT_DIR"

# run bcftools stats 
for vcf in "$INPUT_DIR"/*.chr1.g.vcf.gz; do
	sample=$(basename "$vcf" .chr1.g.vcf.gz)
	echo "Processing $sample..."
	bcftools stats "$vcf" > "$OUTPUT_DIR/${sample}.vcf.stats"
done