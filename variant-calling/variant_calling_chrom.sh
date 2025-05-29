#!/bin/bash
#SBATCH --job-name="variant_calling_chr1"
#SBATCH --mail-user=mhassle1@asu.edu
#SBATCH --mail-type=ALL
#SBATCH -o /data/CEM/wilsonlab/lab_generated/kenya/hassler/logs/%A_%a.out
#SBATCH -e /data/CEM/wilsonlab/lab_generated/kenya/hassler/logs/%A_%a.err
#SBATCH --mem=16G
#SBATCH --time=0-02:00:00
#SBATCH --partition=general
#SBATCH --cpus-per-task=16
#SBATCH --array=0-84

# TODO generate manifest textfile with sample names 
# TODO verify contig name is updated for reference genome (modified for T2T CHM13v2.0)

# exit on silent errors
set -euo pipefail

# load modules
module load mamba/latest
source activate gatk4_env

# define variables
SAMPLES="/path/to/manifest.txt"
REF="/path/to/ref.fna"
CONTIG="CP068270.2"
CHROM="chr8"

# store current sample name inside variable
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLES")

# define input and output
BAM=$(ls /path/to/bams/${SAMPLE}.wgs.*.dedup.bam) # NOTE you can remove dynamic syntax depending on naming scheme
OUTDIR="/path/to/outfile/${CHROM}"
mkdir -p "$OUTDIR"
OUTFILE="${OUTDIR}/${SAMPLE}.${CHROM}.g.vcf.gz"

# run gatk
gatk HaplotypeCaller \
	-R $REF \
	-I $BAM \
	-O $OUTFILE \
	-ERC GVCF \
	-L $CONTIG \
	--native-pair-hmm-threads $SLURM_CPUS_PER_TASK

# example output: ${OUTDIR}/chr8/A100.chr8.g.vcf.gz
