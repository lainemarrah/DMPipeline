sample = x
id = x
outdir = x
tempdir = x
bowtieidx = x
cnvkitref = x

#double check assignments and variable use
# make sure directories are good throughout

#just output so it knows how to backtrack
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
    run:
        shell("prefetch --max-size 100G {id}")
        shell("fastq-dump --outdir {tempdir} --gzip --skip-technical  -F  --split-files --clip ./{id}")
        shell("rm -rf {id}")
        shell("fastp -z 1 -i {input.i1} -I {input.i2} -o {output.o1} -O {output.o2}")

# figure out how to delete fastq if successful

rule bowtie_align:
    input:
        i1="{tempdir}/{id}_1.fastp.gz"
        i2="{tempdir}/{id}_2.fastp.gz"
    output:
        "{tempdir}/{name}.bam"
    shell:
        "bowtie2 -x {bowtieidx} -1 {input.i1} -2 {input.i2} | samtools view -bS - > {output}"

rule fixmate:
    input:
        "{tempdir}/{name}.bam"
    output:
        "{tempdir}/{name}.fixmate"
    shell:
        "samtools fixmate -O bam {input} {output}"

#make sure the indexing part works? not explicitly called by other rules
rule sort_and_index:
    input:
        "{tempdir}.fixmate"
    output:
        bam="{outdir}/{name}.sorted"
        bai="{outdir}/{name}.sorted.bai"
    run:
        shell("samtools sort {input} -m 20G -o {output.bam}")
        shell("samtools index {output.bam}")

rule chr12_split:
    input:
        "{outdir}/{name}.sorted"
    output:
        "{outdir}/{name}.chr12.sorted"
    run:
        shell("samtools view {input} chr12 -b > {output}")
        shell("samtools index {output}")
        shell("samtools flagstat {output}")

rule cnvkit_cns:
    input:
        "{outdir}/{name}.chr12.sorted"
    output:
        "{tempdir}/cnvkit/{name}.chr12.cns"
    shell:
        "cnvkit.py batch {input} -r {cnvkitref} -d {tempdir}/cnvkit/{name}.chr12"

rule cns2bed:
    input:
        "{tempdir}/cnvkit/{name}.chr12.cns"
    output:
        "{outdir}/{name}.chr12.bed"
    run:
        shell("convert_cns_to_bed.py --cns_file={input}")
        shell("cp {tempdir}/cnvkit/{name}_CNV_CALLS.bed {output}")

rule bam2cfg:
    input:
        "{tempdir}/{name}.chr12.sorted"
    output:
        "{tempdir}/{name}.chr12.cfg"
    shell:
        "bam2cfg.pl {input} > {output}"

rule breakdancer:
    input:
        "{tempdir}/{name}.chr12.cfg"
    output:
        "{tempdir}/{name}.chr12.brk"
    shell:
        "breakdancer-max {input} | tail -n +5 > {output}
#make sure above line works with pipe

rule breakdancer2vcf:
    input:
        "{tempdir}/{name}.chr12.brk"
    output:
        "{outdir}/{name}.chr12.vcf"
    shell:
        "breakdancer2vcf.py -i {input} -o {output}"

rule dmfinder:
    input:
        sv="{outdir}/{name}.chr12.vcf"
        cn="{outdir}/{name}.chr12.bed"
        bam="{outdir}/{name}.chr12.sorted"
    output:
        "{outdir}/{name}.chr12.dmrpt"
    shell:
        "dm_find.pl --input_bam {input.bam} --sv {input.sv}--cn {input.cn} --report_file {output}

