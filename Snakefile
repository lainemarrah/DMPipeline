import argparse

parser = argparse.ArgumentParser(description='Input arguments necessary for running DMFinder.')
parser.add_argument('-f1', type=str, help='The full filepath to the first fastq file of paired-end reads.')
parser.add_argument('-f2', type=str, help='The full filepath to the second fastq file of paired-end reads.')
parser.add_argument('-n', type=str, help='The name of your sample. Output files will have this name.')
parser.add_argument('-o', type=str, help='Choose a directory for output files created by this pipeline. This directory should contain an indexed and sorted BAM file.')
parser.add_argument('-p', nargs='?', const=0, type=int, help='[Optional] If you want to multithread, choose number of threads here.')
args = parser.parse_args()

fastq1 = args.1
fastq2 = args.f2
name = args.n
outdir = args.o
threads = args.p
tempdir = outdir+"/temp"

#output so it knows how to backtrack
rule all:
    input:
        outdir+"/"+name+".dmrpt"

rule fastp:
    input:
        i1=fastq1,
        i2=fastq2
    output:
        o1=tempdir+"/"+name+"_1.fastp.gz",
        o2=tempdir+"/"+name+"_2.fastp.gz"
    run:
        shell("fastp -z 1 -i {input.i1} -I {input.i2} -o {output.o1} -O {output.o2}")

rule bowtie_align:
    input:
        i1=tempdir+"/"+name+"_1.fastp.gz",
        i2=tempdir+"/"+name+"_2.fastp.gz"
    output:
        tempdir+"/"+name+".bam"
    shell:
        "bowtie2 -x hg19/hg19full -p "+threads+" -1 {input.i1} -2 {input.i2} | samtools view -bS - > {output}"

rule fixmate:
    input:
        tempdir+"/"+name+".bam"
    output:
        tempdir+"/"+name+".fixmate"
    shell:
        "samtools fixmate -O bam {input} {output}"

rule sort_and_index:
    input:
        tempdir+"/"+name+".fixmate"
    output:
        bam=outdir+"/"+name+".sorted",
        bai=outdir+"/"+name+".sorted.bai"
    run:
        shell("samtools sort {input} -m 20G -o {output.bam}"),
        shell("samtools index {output.bam}")

#rule chr12_split:
#    input:
#        outdir+"/"+name+".sorted"
#    output:
#        outdir+"/"+name+"chr12.sorted"
#    run:
#        shell("samtools view {input} chr12 -b > {output}"),
#        shell("samtools index {output}"),
#        shell("samtools flagstat {output}"),

rule cnvkit_cns:
    input:
        outdir+"/"+name+".sorted"
    output:
        tempdir+"/cnvkit/"+name+".cns"
    shell:
        "cnvkit.py batch {input} -r {cnvkitref} -p {threads}-d {tempdir}/cnvkit/{name}"

rule cns2bed:
    input:
        tempdir+"/cnvkit/"+name+".cns"
    output:
        outdir+"/"+name+".bed"
    run:
        shell("convert_cns_to_bed.py --cns_file={input}"),
        shell("cp {tempdir}/cnvkit/{name}_CNV_CALLS.bed {output}")

rule bam2cfg:
    input:
        tempdir+"/"+name+".sorted"
    output:
        tempdir+"/"+name+".cfg"
    shell:
        "bam2cfg.pl {input} > {output}"

rule breakdancer:
    input:
        tempdir+"/"+name+".cfg"
    output:
        tempdir+"/"+name+".brk"
    shell:
        "breakdancer-max {input} | tail -n +5 > {output}"

rule breakdancer2vcf:
    input:
        tempdir+"/"+name+".brk"
    output:
        outdir+"/"+name+".vcf"
    shell:
        "breakdancer2vcf.py -i {input} -o {output}"

rule dmfinder:
    input:
        sv=outdir+"/"+name+".vcf",
        cn=outdir+"/"+name+".bed",
        bam=outdir+"/"+name+".sorted",
    output:
        outdir+"/"+name+".dmrpt"
    shell:
        "dm_find.pl --input_bam {input.bam} --sv {input.sv}--cn {input.cn} --report_file {output}"

