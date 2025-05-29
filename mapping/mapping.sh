#!/bin/bash
#SBATCH --job-name="mapping"
#SBATCH -o logs/slurm_%A_%a.out
#SBATCH -e logs/slurm_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mhassle1@asu.edu
#SBATCH --time=2-00:00:00
#SBATCH --partition=general
#SBATCH --mem=200G
#SBATCH --cpus-per-task=32
#SBATCH --array=0

# TODO generate manifest textfile with sample names before running 
# TODO make sure naming scheme in lines 42 and 43 match input from $DATA_DIR

# exit on silent errors
set -euo pipefail

#load mamba
module load mamba/latest

# load bwa/latest and samtools/latest
module load bwa-0.7.17-gcc-12.1.0
module load samtools-1.21-gcc-12.1.0
module load picard-2.26.2-gcc-12.1.0

# define data directory & reference sequence
DATA_DIR="/path/to/data"
OUTPUT_DIR="/path/to/output"
REFS="/path/to/ref.fa" 
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" path/to/manifest.txt)
TYPE="wgs"
REF_NAME="t2t_ypars"

# make sure output directory exists
mkdir -p "$OUTPUT_DIR"

# map to reference (16t), convert sam to bam (4t), then sort(12t) # TODO check and change input file names 
bwa mem -t 16 \
	-R "@RG\tID:1\tSM:${SAMPLE}\tLB:lib1\tPU:unit1\tPL:ILLUMINA" \
	$REFS \
	${DATA_DIR}/"${SAMPLE}_R1.fastq.gz" \
	${DATA_DIR}/"${SAMPLE}_R2.fastq.gz" \
	| samtools view -Sb -@ 8 - \
	| samtools sort -@ 8 -o ${OUTPUT_DIR}/${SAMPLE}.${TYPE}.${REF_NAME}.sorted.bam
	#2> /data/CEM/wilsonlab/lab_generated/kenya/imsad/logs/${NAME}_samtools.err	#custom log per sample (optional)


# index sorted bams
samtools index ${OUTPUT_DIR}/${SAMPLE}.${TYPE}.${REF_NAME}.sorted.bam
# example output: A11.sorted.bam

# mark duplicates - #NOTE forgot this the first time but in the pipeline this'll be included
picard MarkDuplicates \
	I=${OUTPUT_DIR}/${SAMPLE}.${TYPE}.${REF_NAME}.sorted.bam \
	O=${OUTPUT_DIR}/${SAMPLE}.${TYPE}.${REF_NAME}.dedup.bam \
	M=/data/CEM/wilsonlab/lab_generated/kenya/hassler/logs/${SAMPLE}_dedup_metrics.txt \
	REMOVE_DUPLICATES=true \
	VALIDATION_STRINGENCY=SILENT \
	CREATE_INDEX=true

# remove intermediate sorted BAM after MarkDuplicates completes (optional)
if [[ $? -eq 0 ]]; then
	rm ${OUTPUT_DIR}/${SAMPLE}.${TYPE}.${REF_NAME}.sorted.bam
	rm ${OUTPUT_DIR}/${SAMPLE}.${TYPE}.${REF_NAME}.sorted.bam.bai
else
	echo "Picard failed for ${SAMPLE}, not removing ${SAMPLE}.wes.sorted.bam" >&2
fi

# example output: A11.wgs.t2t_ypars.sorted.bam
