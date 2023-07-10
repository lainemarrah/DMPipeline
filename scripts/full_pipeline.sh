#!/bin/bash

Help()
{  
   #show help
   echo
   echo "Syntax: full_pipeline.sh [-h] [-i SRA_ID] [-n SAMPLE_NAME] [-o OUTPUT_DIRECTORY] [-p THREADS]"
   echo "options:"  
   echo "f1	First paired-end fastq file filepath."
   echo "f2	Second paired-end fastq file filepath."
   echo "h     Print this Help."
   echo "i     The SRA ID of your sample. This should be something like SRRXXXXXXX."
   echo "n     The name of your sample. Output files will have this name."
   echo "o     Choose a directory for output files created by this pipeline. This directory should contain an indexed and sorted BAM file."
   echo "p     [Optional] If you want to multithread, choose the number of threads." 
}

fastq1=''
fastq2=''
name=''
outdir=''
threads=0

while getopts ":hi:n:o:t:p:r:t" option; do
   case "$option" in
      f1) #assign fastq1
	 fastq1="$OPTARG";;
      f2) #assign fastq2
      	 fastq2="$OPTARG";;
      h) #display Help
         Help
         exit;;
      n) #assign name
         name="$OPTARG";;
      o) #assign output dir
         outdir="$OPTARG";;
      p) #number of threads
         threads="$OPTARG";;
     \?) #invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done

id="$(echo ${fastq1} | grep -oP "SRR.......")"

#making a dir for new files
if [ ! -d ${outdir} ]; then
        mkdir ${outdir}
fi

echo "Making and sorting BAM file..."

#create sam -> bam -> sorted bam
if [ -f ${name}.sorted ]; then
        echo "Sorted BAM already exists"
else
        source /home/lm2ku/.bashrc
        conda activate bowtie2
        bowtie2  -p ${threads} -x hg19/hg19full -1 ${fastq1} -2 ${fastq2} | samtools view -bS > ${outdir}/${name}.bam
        conda deactivate

        samtools fixmate -O bam ${outdir}/${name}.bam ${outdir}/${name}.fixmate
        if [ -s ${outdir}/${name}.fixmate ]; then
                rm ${outdir}/${name}.bam
        fi

        samtools sort -m 20G -o ${outdir}/${name}.sorted ${outdir}/${name}.fixmate
        if [ -s ${outdir}/${name}.sorted ]; then
                rm ${outdir}/${name}.fixmate
        fi

        samtools index ${outdir}/${name}.sorted
fi


if [ -f ${outdir}/${name}.chr12.sorted ]; then
        echo "Sorted and split BAM already exists"
else
        samtools view ${outdir}/${name}.sorted chr12 -b > ${outdir}/${name}.chr12
        java -jar picard.jar AddOrReplaceReadGroups I=${outdir}/${name}.chr12 O=${outdir}/${name}.chr12.sorted SORT_ORDER=coordinate RGID=${id} RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=${name}
        rm ${outdir}/${name}.chr12
        samtools index ${outdir}/${name}.chr12.sorted
        samtools flagstat ${outdir}/${name}.chr12.sorted
fi

if [ ! -f ${outdir}/${name}.chr12.sorted.bai ]; then
        samtools index ${outdir}/${name}.chr12.sorted
        samtools flagstat ${outdir}/${name}.chr12.sorted
fi

echo "Making BED file..."
#create cn file
if [ -f ${outdir}/${name}.chr12.bed ]; then
	echo "BED file already exists"
else
	if [ ! -d ${outdir}/cnvkit ]; then
      		mkdir ${outdir}/cnvkit
	fi
	cnvkit.py batch ${outdir}/${name}.chr12.sorted -r hg19/hg19_cnvkit_filtered_ref.cnn -p ${threads} -d ${outdir}/cnvkit
	scripts/convert_cns_to_bed.py --cns_file=${outdir}/${name}.chr12.cns
	cp ${outdir}/cnvkit/${name}_CNV_CALLS.bed ${outdir}/${name}.chr12.bed
fi

if [ -s ${outdir}/${name}.chr12.bed ]; then
	rm -r ${outdir}/cnvkit
fi

echo "Making VCF file..."
#create sv file
if [ -f ${outdir}/${name}.chr12.vcf ]; then
  echo "VCF file already exists"
else
	scripts/bam2cfg.pl ${outdir}/${name}.chr12.sorted > ${outdir}/${name}.chr12.cfg
	breakdancer-max ${outdir}/${name}.chr12.cfg | tail -n +5 > ${outdir}/${name}.chr12.brk
	scripts/breakdancer2vcf.py -i ${outdir}/${name}.chr12.brk -o ${outdir}/${name}.chr12.vcf
fi

if [ -s ${outdir}/${name}.chr12.vcf ]; then
	rm ${outdir}/${name}.chr12.cfg
 	rm ${outdir}/${name}.chr12.brk
fi

echo "DMFinder prerequisite files complete, running DMFinder..."

#actually run dmfinder
if [ -s ${outdir}/${name}.chr12.vcf -a -s ${outdir}/${name}.chr12.bed ]; then
	dm_find.pl --input_bam ${outdir}/${name}.chr12.sorted --sv ${outdir}/${name}.chr12.vcf --cn ${outdir}/${name}.chr12.bed --report_file ${outdir}/${name}.chr12.dmrpt --graph_file ${outdir}/${name}.chr12.dmgraph --verbose
else
        echo "Your VCF file, BED file, or both are empty."
fi

echo "Pipeline complete"
