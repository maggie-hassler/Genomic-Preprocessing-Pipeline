#!/bin/bash
#SBATCH --job-name="bam_stats"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mhassle1@asu.edu
#SBATCH -o logs/slurm.%A_%a.out
#SBATCH -e logs/slurm.%A_%a.err
#SBATCH --mem=16G
#SBATCH --time=0-01:00:00
#SBATCH --partition=general
#SBATCH --cpus-per-task=16
#SBATCH --array=0-84

# NOTE updated for WGS full directory on 5/15/25

# exit on silent errors
set -euo pipefail 

# load samtools
module load samtools-1.21-gcc-12.1.0

# define paths 
IN_DIR="/path/to/bams"
OUT_DIR="/path/to/output_directory"

# make sure output directory exists
mkdir -p "$OUT_DIR"

# define variables 
NAME=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))"p /path/to/manifest.txt)
TYPE="wgs"

samtools flagstat ${IN_DIR}/${NAME}.${TYPE}.*.dedup.bam > ${OUT_DIR}/${NAME}.samtools.txt
