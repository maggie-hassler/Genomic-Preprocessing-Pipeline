#!/bin/bash
#SBATCH --job-name="generate_sex_lists"
#SBATCH --mail-user=mhassle1@asu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/slurm.%A_%a.out
#SBATCH -e logs/slurm.%A_%a.err
#SBATCH -n 1
#SBATCH --time=0-01:00:00
#SBATCH --partition=htc
#SBATCH --mem=64G

# NOTE set threshold based on R plot
THRESHOLD=5

# define variables 
INPUT="/path/to/concatenated.csv"
FEMALE_OUT="/path/to/females.txt"
MALE_OUT="/path/to/males.txt"

# make sure output exists
mkdir -p "FEMALE_OUT"
mkdir -p "MALE_OUT"

# Clear previous output files if they exist
> "$FEMALE_OUT"
> "$MALE_OUT"

# Read the input file line by line
while IFS=',' read -r SAMPLE YDEPTH; do
	# Skip header line (optional, if present)
	[[ "$SAMPLE" == "Sample" ]] && continue

	# Compare Y depth
	awk -v depth="$YDEPTH" -v name="$SAMPLE" -v t="$THRESHOLD" '
		BEGIN {
			if (depth < t) print name >> "wes_females.txt";
			 else print name >> "wes_males.txt";
		}'
done < "$INPUT"
