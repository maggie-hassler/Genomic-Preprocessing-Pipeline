#!/bin/bash
#SBATCH --job-name="trimming"
#SBATCH --mail-user=mhassle1@asu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/slurm.%A_%a.out
#SBATCH -e logs/slurm.%A_%a.err
#SBATCH --mem=64G
#SBATCH --time=0-02:00:00
#SBATCH --partition=general
#SBATCH --cpus-per-task=16
#SBATCH --array=0

# TODO generate a manfest textfile with sample names 
# TODO verify naming scheme in lines 35, 42, and 43 match input from $INPUT_DIR

# load mamba
module load mamba/latest

# activate the bbduk environment
source activate bbduk_env2

INPUT_DIR="/path/to/data"
OUTPUT_DIR="/path/to/output"
SAMPLES="/path/to/manifest.txt"

# make sure output directory exists
mkdir -p "$OUTPUT_DIR"

# store current file inside a variable
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLES")
echo "Running BBDuk on $SAMPLE..."

# sanity check - make sure the script found the file
# check that input files exist
if [[ ! -f "${INPUT_DIR}/${SAMPLE}_R1.fastq.gz" ]] || [[ ! -f "${INPUT_DIR}/${SAMPLE}_R2.fastq.gz" ]]; then
	echo "Missing input FASTQ files for $SAMPLE, skipping..."
	exit 1
fi

# run bbduk & adjust params as needed 
bbduk.sh \
		in1="${INPUT_DIR}/${SAMPLE}_R1.fastq.gz" \
		in2="${INPUT_DIR}/${SAMPLE}_R2.fastq.gz" \
		out1="${OUTPUT_DIR}/${SAMPLE}_R1_trimmed.fastq.gz" \
		out2="${OUTPUT_DIR}/${SAMPLE}_R2_trimmed.fastq.gz" \
		ref=/path/to/adapters.fa \
		ktrim=r \
		k=21 \
		mink=11 \
		hdist=2 \
		qtrim=rl \
		trimq=15 \
		minlen=75 \
		maq=20 \
		threads=16
