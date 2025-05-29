# import config
configfile: "config_trimming.json"

# store input data in variables 
SAMPLES = config["test_samples"]			# which samples to process
#SAMPLE_PATH = config["samples"] 			# where to find them 
OUTPUT_DIR = config["output_directory"] 	# path to /results 
ADAPTERS = config["adapters"] 

trimmed_r1 = expand("{output}/trimmed/{sample}_trimmed_R1.fastq.gz", output=OUTPUT_DIR, sample=SAMPLES)
trimmed_r2 = expand("{output}/trimmed/{sample}_trimmed_R2.fastq.gz", output=OUTPUT_DIR, sample=SAMPLES)

rule all:
	input:
		trimmed_r1,
		trimmed_r2

rule bbduk_trim:
	input:
		r1 = lambda wc: config["samples"][wc.sample]["fq_r1"],
		r2 = lambda wc: config["samples"][wc.sample]["fq_r2"],
		adapters = lambda wc: config["adapters"]
		#adapters = ADAPTERS
	output:
		# r1_trimmed = lambda wc: f"{OUTPUT_DIR}/trimmed/{wc.sample}_trimmed_R1.fastq.gz",
		# r2_trimmed = lambda wc: f"{OUTPUT_DIR}/trimmed/{wc.sample}_trimmed_R2.fastq.gz"
		r1_trimmed = "{OUTPUT_DIR}/trimmed/{sample}_trimmed_R1.fastq.gz",
		r2_trimmed = "{OUTPUT_DIR}/trimmed/{sample}_trimmed_R2.fastq.gz"
	threads: 2
	conda:
		"bbduk_env2"
	shell:
		"""
		bbduk.sh \
			in1={input.r1} \
			in2={input.r2} \
			out1={output.r1_trimmed} \
			out2={output.r2_trimmed} \
			ref={input.adapters} \
			ktrim=r k=23 mink=11 hdist=1 \
			qtrim=r trimq=20 minlen=30 minavgquality=20 maq=20 \
			threads={threads} \
						-eoom
		"""
