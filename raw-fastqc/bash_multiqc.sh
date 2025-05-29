#!/bin/bash
#SBATCH --job-name="multiqc"
#SBATCH --mail-user=mhassle1@asu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/slurm.%A_%a.out
#SBATCH -e logs/slurm.%A_%a.err
#SBATCH -n 1
#SBATCH --time=0-01:00:00
#SBATCH --partition=htc
#SBATCH --mem=64G

# define variables 
JOB_NAME="my_job"
DATA_DIR="/path/to/data"
OUTPUT_DIR="/path/to/output"

# make sure output directory exists 
mkdir -p "$OUTPUT_DIR"

# load modules 
module load mamba/latest
source activate multiqc_env2 

# run multiqc
multiqc "$DATA_DIR" -o "$OUTPUT_DIR"

# rename report 
mv "$OUTPUT_DIR/multiqc_report.html" "$OUTPUT_DIR/${JOB_NAME}.multiqc.html"
