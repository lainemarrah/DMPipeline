import argparse

#todo: add info on where to get ref genome for cnvkit and how to make bowtie2 index to readme
#todo: make conda environment with everything necessary in it for download - this would make prereqs very easy
#todo: figure out how to skip aligning if there is already a sorted bam, could make separate py scripts but there is probably a more elegant way
#todo: add cleanup function to delete tempdir + files, find how to call it in workflow
#todo: test on my machine

parser = argparse.ArgumentParser(description='Input arguments necessary for running DMFinder.')
parser.add_argument('-i', type=str, help='The SRA ID of your sample. This should be something like SRRXXXXXXX.')
parser.add_argument('-n', type=str, help='The name of your sample. Output files will have this name.')
parser.add_argument('-o', type=str, help='Choose a directory for output files created by this pipeline. This directory should contain an indexed and sorted BAM file.')
parser.add_argument('-p', nargs='?', const=0, type=int, help='[Optional] If you want to multithread, choose number of threads here.')
parser.add_argument('-r', nargs='?', type=str, help='[Optional for BED] Reference genome .cnn file. This is for the CNVkit portion of the pipeline.')
parser.add_argument('-t', type=str, help='Choose a directory for temporary files to be downloaded to.')
args = parser.parse_args()

id = args.i
sample = args.n
outdir = args.o
threads = args.p
cnvkitref = args.r
tempdir = args.t

#output so it knows how to backtrack
rule all:
    input:
        "{sample}.dmrpt"

#rule chr12_split:
#    input:
#        "{outdir}/{name}.sorted"
#    output:
#        "{outdir}/{name}.chr12.sorted"
#    params:
#        outdir=outdir
#        name=name
#    run:
#        shell("samtools view {input} chr12 -b > {output}")
#        shell("samtools index {output}")
#        shell("samtools flagstat {output}")

rule cnvkit_cns:
    input:
        "{outdir}/{name}.sorted"
    output:
        "{tempdir}/cnvkit/{name}.cns"
    params:
        tempdir=tempdir
        outdir=outdir
        name=name
        cnvkitref=cnvkitref
        threads=threads
    shell:
        "cnvkit.py batch {input} -r {cnvkitref} -p {threads}-d {tempdir}/cnvkit/{name}"

rule cns2bed:
    input:
        "{tempdir}/cnvkit/{name}.cns"
    output:
        "{outdir}/{name}.bed"
    params:
        tempdir=tempdir
        outdir=outdir
        name=name
    run:
        shell("convert_cns_to_bed.py --cns_file={input}")
        shell("cp {tempdir}/cnvkit/{name}_CNV_CALLS.bed {output}")

rule bam2cfg:
    input:
        "{tempdir}/{name}.sorted"
    output:
        "{tempdir}/{name}.cfg"
    params:
        tempdir=tempdir
        name=name
    shell:
        "bam2cfg.pl {input} > {output}"

rule breakdancer:
    input:
        "{tempdir}/{name}.cfg"
    output:
        "{tempdir}/{name}.brk"
    params:
        tempdir=tempdir
        name=name
    shell:
        "breakdancer-max {input} | tail -n +5 > {output}"

rule breakdancer2vcf:
    input:
        "{tempdir}/{name}.brk"
    output:
        "{outdir}/{name}.vcf"
    params:
        tempdir=tempdir
        outdir=outdir
        name=name
    shell:
        "breakdancer2vcf.py -i {input} -o {output}"

rule dmfinder:
    input:
        sv="{outdir}/{name}.vcf"
        cn="{outdir}/{name}.bed"
        bam="{outdir}/{name}.sorted"
    output:
        "{outdir}/{name}.dmrpt"
    params:
        outdir=outdir
        name=name
    shell:
        "dm_find.pl --input_bam {input.bam} --sv {input.sv}--cn {input.cn} --report_file {output}"

