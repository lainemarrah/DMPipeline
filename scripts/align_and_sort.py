import argparse

#manipulate ref so we can take off .fa suffix to use for bowtie index
#if statement to see if bowtie indexing is necessary

parser = argparse.ArgumentParser(description='Input arguments necessary for running DMFinder.')
parser.add_argument('-i', type=str, help='The SRA ID of your sample. This should be something like SRRXXXXXXX.')
parser.add_argument('-n', type=str, help='The name of your sample. Output files will have this name.')
parser.add_argument('-o', type=str, help='Choose a directory for output files created by this pipeline. This directory should contain an indexed and sorted BAM file.')
parser.add_argument('-p', nargs='?', const=0, type=int, help='[Optional] If you want to multithread, choose number of threads here.')
parser.add_argument('-r', nargs='?', type=str, help='[Optional for BED] Reference genome .cnn file. This is for the CNVkit portion of the pipeline.')
parser.add_argument('-t', type=str, help='Choose a directory for temporary files to be downloaded to.')
parser.add_argument('-x', nargs='?', type=str, help='[Optional for aligning] The Bowtie2 index name. This should be everything but the suffix for your Bowtie2 reference index files.'
args = parser.parse_args()

id = args.i
sample = args.n
outdir = args.o
threads = args.p
cnvkitref = args.r
tempdir = args.t
bowtieidx = args.x

#output so it knows how to backtrack
rule all:
    input:
        "{outdir}/{name}.sorted"

rule fastp:
    input:
        i1="{tempdir}/{id}_1.fastq.gz"
        i2="{tempdir}/{id}_2.fastq.gz"
    output:
        o1="{tempdir}/{id}_1.fastp.gz"
        o2="{tempdir}/{id}_2.fastp.gz"
    params:
        tempdir=tempdir
        id=id
    run:
        shell("fastp -z 1 -i {input.i1} -I {input.i2} -o {output.o1} -O {output.o2}")

rule bowtie_index:
    input:
    output:
    params:
    shell:
      "bowtie2-build --large-index /home/lm2ku/AA_data_repo/hg19/hg19full.fa /home/lm2ku/AA_data_repo/hg19/hg19full"

rule bowtie_align:
    input:
        i1="{tempdir}/{id}_1.fastp.gz"
        i2="{tempdir}/{id}_2.fastp.gz"
    output:
        "{tempdir}/{name}.bam"
    params:
        tempdir=tempdir
        id=id
        name=name
        bowtieidx=bowtieidx
        threads=threads
    shell:
        "bowtie2 -x {bowtieidx} -p {threads} -1 {input.i1} -2 {input.i2} | samtools view -bS - > {output}"

rule fixmate:
    input:
        "{tempdir}/{name}.bam"
    output:
        "{tempdir}/{name}.fixmate"
    params:
        tempdir=tempdir
        name=name
    shell:
        "samtools fixmate -O bam {input} {output}"

rule sort_and_index:
    input:
        "{tempdir}.fixmate"
    output:
        bam="{outdir}/{name}.sorted"
        bai="{outdir}/{name}.sorted.bai"
    params:
        tempdir=tempdir
        outdir=outdir
        name=name
    run:
        shell("samtools sort {input} -m 20G -o {output.bam}")
        shell("samtools index {output.bam}")

