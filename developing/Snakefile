configfile: "config.yaml"

rule all:
    input:
        "config[outdir]/config[name].dmrpt"

rule fastp:
    input:
        i1="config[fastq1]",
        i2="config[fastq2]"
    output:
        o1="config[tempdir]/config[name]_1.fastp.gz",
        o2="config[tempdir]/config[name]_2.fastp.gz"
    run:
        shell("fastp -z 1 -i {input.i1} -I {input.i2} -o {output.o1} -O {output.o2}")

rule bowtie_align:
    input:
        o1="config[tempdir]/config[name]_1.fastp.gz",
        o2="config[tempdir]/config[name]_2.fastp.gz"
    output:
        "config[tempdir]/config[name].bam"
    shell:
        "bowtie2 -x hg19/hg19full -p {config[threads]} -1 {input.i1} -2 {input.i2} | samtools view -bS - > {output}"

rule fixmate:
    input:
        "config[tempdir]/config[name].bam"
    output:
        "config[tempdir]/config[name].fixmate"
    shell:
        "samtools fixmate -O bam {input} {output}"

rule sort_and_index:
    input:
        "config[tempdir]/config[name].fixmate"
    output:
        bam="config[outdir]/config[name].sorted",
        bai="config[outdir]/config[name].sorted.bai"
    run:
        shell("samtools sort {input} -m 20G -o {output.bam}"),
        shell("samtools index {output.bam}")

#rule chr12_split:
#    input:
#        "config[outdir]/config[name].sorted"
#    output:
#        "config[outdir]/config[name].chr12.sorted"
#    run:
#        shell("samtools view {input} chr12 -b > {output}"),
#        shell("samtools index {output}"),
#        shell("samtools flagstat {output}")

rule cnvkit_cns:
    input:
        "config[outdir]/config[name].sorted",
    output:
        config[tempdir]/cnvkit/config[name].cns"
    shell:
        "cnvkit.py batch {input} -r hg19/hg19_cnvkit_filtered_ref.cnn -p {config[threads]} -d {config[tempdir]}"

rule cns2bed:
    input:
        config[tempdir]/cnvkit/config[name].cns"
    output:
        "config[outdir]/config[name].bed"
    run:
        shell("convert_cns_to_bed.py --cns_file={input}"),
        shell("cp {config[tempdir]}/{name}_CNV_CALLS.bed {output}")

rule bam2cfg:
    input:
        "config[outdir]/config[name].sorted"
    output:
        "config[tempdir]/config[name].cfg"
    shell:
        "bam2cfg.pl {input} > {output}"

rule breakdancer:
    input:
        "config[tempdir]/config[name].cfg"
    output:
        "config[tempdir]/config[name].brk"
    shell:
        "breakdancer-max {input} | tail -n +5 > {output}"

rule breakdancer2vcf:
    input:
        "config[tempdir]/config[name].brk"
    output:
        "config[outdir]/config[name].vcf"
    shell:
        "breakdancer2vcf.py -i {input} -o {output}"

rule dmfinder:
    input:
        sv="config[outdir]/config[name].vcf",
        cn="config[outdir]/config[name].bed",
        bam="config[outdir]/config[name].sorted"
    output:
        "config[outdir]/config[name].dmrpt"
    shell:
        "dm_find.pl --input_bam {input.bam} --sv {input.sv}--cn {input.cn} --report_file {output}"

