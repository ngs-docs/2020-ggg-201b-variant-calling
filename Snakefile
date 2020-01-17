rule download_data:
    shell:
        "wget https://osf.io/4rdza/download -O SRR2584857_1.fastq.gz"

rule download_genome:
    shell:
        "wget https://osf.io/8sm92/download -O ecoli-rel606.fa.gz"

rule uncompress_genome:
    shell:
        "gunzip ecoli-rel606.fa.gz"

rule index_genome_bwa:
    output:
        "ecoli-rel606.fa.amb",
        "ecoli-rel606.fa.ann",
        "ecoli-rel606.fa.bwt",
        "ecoli-rel606.fa.pac",
        "ecoli-rel606.fa.sa"
    shell:
        "bwa index ecoli-rel606.fa"

rule map_reads:
    input:
        "ecoli-rel606.fa.amb",
        "SRR2584857_1.fastq.gz"
    shell:
        "bwa mem -t 4 ecoli-rel606.fa SRR2584857_1.fastq.gz > SRR2584857.sam"

rule index_genome_samtools:
    input:
        "ecoli-rel606.fa"
    output:
        "ecoli-rel606.fa.fai"
    shell:
        "samtools faidx ecoli-rel606.fa"

rule samtools_import:
    input:
        "ecoli-rel606.fa.fai", "SRR2584857.sam"
    output:
        "SRR2584857.bam"
    shell:
        "samtools import ecoli-rel606.fa.fai SRR2584857.sam SRR2584857.bam"

rule samtools_sort:
    input:
        "SRR2584857.bam"
    output:
        "SRR2584857.sorted.bam"
    shell:
        "samtools sort SRR2584857.bam -o SRR2584857.sorted.bam"

rule samtools_index_sorted:
    shell: "samtools index SRR2584857.sorted.bam"

rule samtools_mpileup:
    input: "ecoli-rel606.fa", "SRR2584857.sorted.bam.bai"
    output: "variants.raw.bcf"
    shell:
        """samtools mpileup -u -t DP -f ecoli-rel606.fa SRR2584857.sorted.bam | \
    bcftools call -mv -Ob -o - > variants.raw.bcf"""

rule make_vcf:
    output: "variants.vcf"
    shell: "bcftools view variants.raw.bcf > variants.vcf"

## samtools tview -p ecoli:4202391 SRR2584857.sorted.bam ecoli-rel606.fa
