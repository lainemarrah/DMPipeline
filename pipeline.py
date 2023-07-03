sample = x
id = x
outdir = x
tempdir = x
bowtieidx = x
cnvkitref = x

#todo: assign variables by parsing from command line

#output so it knows how to backtrack
rule all:
    input:
        "{sample}.dmrpt"

#add a python if statement to see if aligning and sorting is needed
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
    shell:
        "bowtie2 -x {bowtieidx} -1 {input.i1} -2 {input.i2} | samtools view -bS - > {output}"

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

#make sure the indexing part works? not explicitly called by other rules
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

rule chr12_split:
    input:
        "{outdir}/{name}.sorted"
    output:
        "{outdir}/{name}.chr12.sorted"
    params:
        outdir=outdir
        name=name
    run:
        shell("samtools view {input} chr12 -b > {output}")
        shell("samtools index {output}")
        shell("samtools flagstat {output}")

rule cnvkit_cns:
    input:
        "{outdir}/{name}.chr12.sorted"
    output:
        "{tempdir}/cnvkit/{name}.chr12.cns"
    params:
        tempdir=tempdir
        outdir=outdir
        name=name
        cnvkitref=cnvkitref
    shell:
        "cnvkit.py batch {input} -r {cnvkitref} -d {tempdir}/cnvkit/{name}.chr12"

rule cns2bed:
    input:
        "{tempdir}/cnvkit/{name}.chr12.cns"
    output:
        "{outdir}/{name}.chr12.bed"
    params:
        tempdir=tempdir
        outdir=outdir
        name=name
    run:
        shell("convert_cns_to_bed.py --cns_file={input}")
        shell("cp {tempdir}/cnvkit/{name}_CNV_CALLS.bed {output}")

rule bam2cfg:
    input:
        "{tempdir}/{name}.chr12.sorted"
    output:
        "{tempdir}/{name}.chr12.cfg"
    params:
        tempdir=tempdir
        name=name
    shell:
        "bam2cfg.pl {input} > {output}"

rule breakdancer:
    input:
        "{tempdir}/{name}.chr12.cfg"
    output:
        "{tempdir}/{name}.chr12.brk"
    params:
        tempdir=tempdir
        name=name
    shell:
        "breakdancer-max {input} | tail -n +5 > {output}"

rule breakdancer2vcf:
    input:
        "{tempdir}/{name}.chr12.brk"
    output:
        "{outdir}/{name}.chr12.vcf"
    params:
        tempdir=tempdir
        outdir=outdir
        name=name
    shell:
        "breakdancer2vcf.py -i {input} -o {output}"

rule dmfinder:
    input:
        sv="{outdir}/{name}.chr12.vcf"
        cn="{outdir}/{name}.chr12.bed"
        bam="{outdir}/{name}.chr12.sorted"
    output:
        "{outdir}/{name}.chr12.dmrpt"
    params:
        outdir=outdir
        name=name
    shell:
        "dm_find.pl --input_bam {input.bam} --sv {input.sv}--cn {input.cn} --report_file {output}"

