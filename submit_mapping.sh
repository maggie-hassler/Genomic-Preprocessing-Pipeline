#!/bin/bash
#SBATCH --job-name="mapping_dryrun"
#SBATCH -o logs/slurm.%j.out
#SBATCH -e logs/slurm.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mhassle1@asu.edu
#SBATCH --time=0-00:00:00
#SBATCH --partition=general
#SBATCH --export=NONE
#SBATCH --mem=80G
#SBATCH --cpus-per-task=16

# load mamba
module load mamba/latest

# activate snakemake environment
source activate snakemake_env

# run mapping Snakefile # NOTE DRY RUN FLAG PRESENT
snakemake \
	--snakefile mapping.smk \
	--configfile master_config_wgs_vcf.json \
	--dry-run \
	--cores ${SLURM_CPUS_PER_TASK} \
	--printshellcmds \
	--latency-wait 60 \
	--use-conda 
