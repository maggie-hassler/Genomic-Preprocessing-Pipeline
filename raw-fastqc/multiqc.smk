import os

# import config #TODO change as needed 
configfile: "master_config_wgs_vcf.json"

# store input in variable, cleaner to read # TODO change as needed
SAMPLES = config["test_samples"]
READS = ["R1", "R2"]

rule all:
	input:
		expand("qc_reports/{sample}_{read}_001_fastqc.html", sample=SAMPLES, read=READS),
		"qc_reports/multiqc_report.html"

# generate fastqcs 
rule fastqc_analysis:
	input:
		lambda wc: config["samples"][wc.sample][f"fq_{wc.read.lower()}"]
	output:
		"qc_reports/{sample}_{read}_001_fastqc.html"
		#html = temp("qc_reports/{sample}_{read}_fastqc.html")
	threads: 4
	shell:
		"{config[fastqc_path]} -o qc_reports {input}"

# generate multiqc 
rule multiqc_analysis:
	input:
		expand("qc_reports/{sample}_{read}_001_fastqc.html", sample=SAMPLES, read=READS)
	output:
		"qc_reports/multiqc_report.html",
	conda:
		"multiqc_env2"
	shell: 
		"multiqc -o qc_reports qc_reports"
# done 
