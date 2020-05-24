# download data from https://osf.io/vzfc6/

rule download_data:
    output: "SRR2584857_1.fastq.gz"
    shell:
        "wget https://osf.io/4rdza/download -O {output}"

rule download_genome:
    output:
        "ecoli-rel606.fa.gz"
    shell:
        "wget https://osf.io/8sm92/download -O {output}"

rule uncompress_genome:
    input: "ecoli-rel606.fa.gz"
    output: "ecoli-rel606.fa"
    shell:
        "gunzip -c {input} > {output}"

rule index_genome_bwa:
    input: "ecoli-rel606.fa"
    output:
        "ecoli-rel606.fa.sa",
        "ecoli-rel606.fa.amb",
        "ecoli-rel606.fa.ann",
        "ecoli-rel606.fa.pac",
        "ecoli-rel606.fa.bwt",
    shell:
        "bwa index {input}"

rule map_reads:
    input:
        genome="ecoli-rel606.fa",
        reads="SRR2584857_1.fastq.gz",
        index="ecoli-rel606.fa.sa",
    output: "SRR2584857.sam"
    shell:
        "bwa mem -t 4 {input.genome} {input.reads} > {output}"

rule index_genome_samtools:
    input: "ecoli-rel606.fa"
    output: "ecoli-rel606.fa.fai"
    shell:
        "samtools faidx {input}"

rule samtools_import:
    input:
        index="ecoli-rel606.fa.fai",
        samfile="SRR2584857.sam"
    output:
        "SRR2584857.bam"
    shell:
        "samtools import {input.index} {input.samfile} {output}"

rule samtools_sort:
    input: "SRR2584857.bam"
    output: "SRR2584857.sorted.bam"
    shell:
        "samtools sort {input} -o {output}"

rule samtools_index_sorted:
    input: "SRR2584857.sorted.bam"
    output: "SRR2584857.sorted.bam.bai"
    shell: "samtools index {input}"

rule samtools_mpileup:
    input:
        index="ecoli-rel606.fa",
        sorted="SRR2584857.sorted.bam",
        sorted_bai="SRR2584857.sorted.bam.bai"
    output:
        "variants.raw.bcf"
    shell: """
        samtools mpileup -u -t DP -f {input.index} {input.sorted} | \
             bcftools call -mv -Ob -o - > {output}
    """

rule make_vcf:
    input: "variants.raw.bcf",
    output: "variants.vcf"
    shell: "bcftools view {input} > {output}"

# at end, run:
## samtools tview -p ecoli:4202391 SRR2584857.sorted.bam ecoli-rel606.fa
