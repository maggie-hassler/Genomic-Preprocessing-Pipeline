#!/bin/bash
#SBATCH --job-name="multiqc"
#SBATCH -o logs/slurm.%j.out
#SBATCH -e logs/slurm.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mhassle1@asu.edu
#SBATCH --time=0-03:00:00
#SBATCH --partition=htc
#SBATCH --export=NONE #isolate the job (start job from clean slate)
#SBATCH -c 8

# load mamba
module load mamba/latest

# load snakemake env
source activate snakemake_env

# run snakefile
snakemake \
	--snakefile multiqc.smk \
	--configfile master_config_wgs_vcf.json \
	--latency-wait 60 \
	--cores 8 \
	--use-conda \
	--conda-frontend mamba \
	--printshellcmds
