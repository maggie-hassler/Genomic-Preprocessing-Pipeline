#!/bin/bash
#SBATCH --job-name="concat_raw_fastq"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mhassle1@asu.edu
#SBATCH -o logs/slurm.%j.out
#SBATCH -e logs/slurm.%j.err
#SBATCH --time=5-00:00:00
#SBATCH --partition=general
#SBATCH --mem=64G

OUT_DIR="/path/to/ouput_directory"
mkdir -p "$OUT_DIR"

# get unique identifier
for sample in $(ls *_R1_001.fastq.gz | cut -d'_' -f1 | sort -u); do
	echo "Processing sample: $sample"

	# find and sort all R1 and R2 files for given unique identifier 
	R1_files=($(ls ${sample}_*_R1_001.fastq.gz | sort))
	R2_files=($(ls ${sample}_*_R2_001.fastq.gz | sort))

	# check the two variables
	# (1) that # of items in r1_files = r2_files
	if [ "${#R1_files[@]}" -ne "${#R2_files[@]}" ]; then
		echo "Mismatched number of R1 and R2 files for sample $sample! Skipping..."
		continue
	fi

	# (2) pairing check - make sure all the sample names in r1_files = r2_files 
	for i in "${!R1_files[@]}"; do
		r1="${R1_files[$i]}"
		r2="${R2_files[$i]}"
		# strip _R1_ vs _R2_ and compare rest of name
		if [[ "${r1/_R1_/_}" != "${r2/_R2_/_}" ]]; then
			echo "Unmatching names in pair: $r1 and $r2"
				continue 2
		fi
	done

	# sample matching = complete

	# merge R1s and R2s
	cat "${R1_files[@]}" > ${OUT_DIR}/${sample}_R1_merged.fastq.gz
	cat "${R2_files[@]}" > ${OUT_DIR}/${sample}_R2_merged.fastq.gz
done
