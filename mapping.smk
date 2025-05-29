import os

# import config
configfile: "master_config_wgs_vcf.json"

SAMPLES = config["one_sample"]
REF = config["ypars_ref"]
OUTDIR = config["output_directory"]
DATADIR = config["trimmed_fastq_dir"]

# make sure logs and stats directories exist 
os.makedirs(f"{OUTDIR}/logs", exist_ok=True)
os.makedirs(f"{OUTDIR}/stats", exist_ok=True)
os.makedirs(f"{OUTDIR}/tmp", exist_ok=True)

# define output: bams + multiqc report 
rule all:
	input:
		expand("{out}/{sample}.wes.dedup.bam.bai", out=OUTDIR, sample=SAMPLES),
		f"{OUTDIR}/multiqc_report.html"

rule map_and_sort:
	input:
		r1 = lambda wc: f"{DATADIR}/{wc.sample}_R1_001.fastq.gz",
		r2 = lambda wc: f"{DATADIR}/{wc.sample}_R2_001.fastq.gz",
		ref = REF
	output:
		bam = temp(f"{OUTDIR}/{{sample}}.wes.sorted.bam")
	threads: 32
	shell:
		"""
		bwa mem -t 16 \
			-R "@RG\\tID:1\\tSM:{wildcards.sample}\\tLB:lib1\\tPU:unit1\\tPL:ILLUMINA" \
			{input.ref} \
			{input.r1} {input.r2} \
			| samtools view -Sb -@ 8 - \
			| samtools sort -@ 8 -o {output.bam}
		"""

rule mark_duplicates:
	input:
		bam = f"{OUTDIR}/{{sample}}.wes.sorted.bam"
	output:
		bam = f"{OUTDIR}/{{sample}}.wes.dedup.bam",
		bai = f"{OUTDIR}/{{sample}}.wes.dedup.bam.bai",
		metrics = f"{OUTDIR}/logs/{{sample}}_dedup_metrics.txt"
	params:
		tmpdir = f"{OUTDIR}/tmp/{{sample}}"
	shell:
		"picard MarkDuplicates I={input.bam} O={output.bam} M={output.metrics} REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true TMP_DIR={params.tmpdir}"

rule samtools_stats:
	input:
		bam = f"{OUTDIR}/{{sample}}.wes.dedup.bam"
	output:
		txt = temp(f"{OUTDIR}/stats/{{sample}}_stats.txt")
	shell:
		"samtools flagstat {input.bam} > {output.txt}"

rule multiqc:
	input:
		lambda wc: expand(rules.samtools_stats.output.txt, sample=SAMPLES)
	output:
		f"{OUTDIR}/multiqc_report.html"
	shell:
		f"multiqc {{input}} -o {OUTDIR}/stats"
