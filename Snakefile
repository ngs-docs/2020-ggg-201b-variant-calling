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
    shell:
        "samtools faidx ecoli-rel606.fa"

rule samtools_import:
    shell:
        "samtools import ecoli-rel606.fa.fai SRR2584857.sam SRR2584857.bam"

rule samtools_sort:
    shell:
        "samtools sort SRR2584857.bam -o SRR2584857.sorted.bam"

rule samtools_index_sorted:
    shell: "samtools index SRR2584857.sorted.bam"


rule samtools_mpileup:
    shell:
        """samtools mpileup -u -t DP -f ecoli-rel606.fa SRR2584857.sorted.bam | \
    bcftools call -mv -Ob -o - > variants.raw.bcf"""

rule make_vcf:
    shell: "bcftools view variants.raw.bcf > variants.vcf"

# at end, run:
## samtools tview -p ecoli:4202391 SRR2584857.sorted.bam ecoli-rel606.fa
