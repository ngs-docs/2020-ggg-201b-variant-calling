# list out samples
SAMPLES=['SRR2584403_1',
	 'SRR2584405_1',
         'SRR2584404_1',
         'SRR2584857_1']

rule all:
    input:
        # create a new filename for every entry in SAMPLES,
        # replacing {name} with each entry.
        expand("{name}-variants.vcf", name=SAMPLES)

# copy data from /home/ctbrown/data/ggg201b
rule copy_data:
    input:
        "/home/ctbrown/data/ggg201b/{sample}.fastq.gz"
    output: "{sample,[a-zA-Z0-9_]+}.fastq.gz"
    shell:
        "ln -fs {input} {output}"

rule download_genome:
    output:
        "ecoli-rel606.fa.gz"
    shell:
        "wget https://osf.io/8sm92/download -O ecoli-rel606.fa.gz"

rule uncompress_genome:
    input:
        "ecoli-rel606.fa.gz"
    output:
        "ecoli-rel606.fa"
    shell:
        "gunzip ecoli-rel606.fa.gz"

rule index_genome_bwa:
    input:
        "ecoli-rel606.fa"
    output:
        "ecoli-rel606.fa.amb",
        "ecoli-rel606.fa.ann",
        "ecoli-rel606.fa.sa"
    shell:
        "bwa index ecoli-rel606.fa"

rule map_reads:
    input:
        "ecoli-rel606.fa.amb",
        "{sample}.fastq.gz"
    output:
        "{sample}.sam"
    shell:
        "bwa mem -t 4 ecoli-rel606.fa {wildcards.sample}.fastq.gz > {wildcards.sample}.sam"

rule index_genome_samtools:
    input:
        "ecoli-rel606.fa"
    output:
        "ecoli-rel606.fa.fai"
    shell:
        "samtools faidx ecoli-rel606.fa"


rule samtools_import:
    input:
        "ecoli-rel606.fa.fai", "{sample}.sam"
    output:
        "{sample}.bam"
    shell:
        "samtools import ecoli-rel606.fa.fai {wildcards.sample}.sam {wildcards.sample}.bam"

rule samtools_sort:
    input:
        "{sample}.bam"
    output:
        "{sample}.sorted.bam"
    shell:
        "samtools sort {wildcards.sample}.bam -o {wildcards.sample}.sorted.bam"

rule samtools_index_sorted:
    input: "{sample}.sorted.bam"
    output: "{sample}.sorted.bam.bai"
    shell: "samtools index {wildcards.sample}.sorted.bam"


rule samtools_mpileup:
    input: "ecoli-rel606.fa", "{sample}.sorted.bam"
    output: "{sample}-variants.raw.bcf"
    shell:
        """samtools mpileup -u -t DP -f ecoli-rel606.fa {wildcards.sample}.sorted.bam | \
    bcftools call -mv -Ob -o - > {wildcards.sample}-variants.raw.bcf"""

rule make_vcf:
    input: "{sample}-variants.raw.bcf"
    output: "{sample}-variants.vcf"
    shell: "bcftools view {wildcards.sample}-variants.raw.bcf > {wildcards.sample}-variants.vcf"

## samtools tview -p ecoli:4202391 SRR2584857_1.sorted.bam ecoli-rel606.fa
