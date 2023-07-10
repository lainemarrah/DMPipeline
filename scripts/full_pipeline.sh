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

f1=''
f2=''
id=''
name=''
outdir=''
threads=0

while getopts ":hi:n:o:t:p:r:t" option; do
   case "$option" in
      h) #display Help
         Help
         exit;;
      i) #assign id
         id="$OPTARG";;
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

#making a dir for new files
if [ ! -d ${outdir} ]; then
        mkdir ${outdir}
fi
cd ${outdir}

echo "Making and sorting BAM file..."

#create sam -> bam -> sorted bam
if [ -f ${name}.sorted ]; then
        echo "Sorted BAM already exists"
else
        source /home/lm2ku/.bashrc
        conda activate bowtie2
        bowtie2  -p ${threads} -x hg19/hg19full -1 ${fastq1} -2 ${fastq2} | samtools view -bS > ${name}.bam
        conda deactivate

        samtools fixmate -O bam ${name}.bam ${name}.fixmate
        if [ -s ${name}.fixmate ]; then
                rm ${name}.bam
        fi

        samtools sort -m 20G -o ${name}.sorted ${name}.fixmate
        if [ -s ${name}.sorted ]; then
                rm ${name}.fixmate
        fi

        samtools index ${name}.sorted
fi


if [ -f ${name}.chr12.sorted ]; then
        echo "Sorted and split BAM already exists"
else
        samtools view ${name}.sorted chr12 -b > ${name}.chr12
        java -jar picard.jar AddOrReplaceReadGroups I=${name}.chr12 O=${name}.chr12.sorted SORT_ORDER=coordinate RGID=${id} RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=${name}
        rm ${name}.chr12
        samtools index ${name}.chr12.sorted
        samtools flagstat ${name}.chr12.sorted
fi

if [ ! -f ${name}.chr12.sorted.bai ]; then
        samtools index ${name}.chr12.sorted
        samtools flagstat ${name}.chr12.sorted
fi

echo "Making BED file..."
#create cn file
if [ -f ${name}.chr12.bed ]; then
	echo "BED file already exists"
else
	if [ ! -d ${outdir}/cnvkit ]; then
      		mkdir ${outdir}/cnvkit
	fi
	cd ${outdir}/cnvkit
	cnvkit.py batch ${outdir}/${name}.chr12.sorted -r hg19/hg19_cnvkit_filtered_ref.cnn -p ${threads} -d ${outdir}/cnvkit
	scripts/convert_cns_to_bed.py --cns_file=${name}.chr12.cns
	cp ${name}_CNV_CALLS.bed ${outdir}/${name}.chr12.bed
fi

if [ -s ${outdir}/${name}.chr12.bed ]; then
	rm -r ${outdir}/cnvkit
fi

echo "Making VCF file..."
#create sv file
if [ -f ${name}.chr12.vcf ]; then
  echo "VCF file already exists"
else
	scripts/bam2cfg.pl ${outdir}/${name}.chr12.sorted > ${name}.chr12.cfg
	breakdancer-max ${name}.chr12.cfg | tail -n +5 > ${name}.chr12.brk
	breakdancer2vcf.py -i ${name}.chr12.brk -o ${outdir}/${name}.chr12.vcf
fi

echo "DMFinder prerequisite files complete, running DMFinder..."

cd ${outdir}
#actually run dmfinder
if [ -s ${name}.chr12.vcf -a -s ${name}.chr12.bed ]; then
	/home/lm2ku/DMFinder/dm_find.pl --input_bam ${name}.chr12.sorted --sv ${name}.chr12.vcf --cn ${name}.chr12.bed --report_file ${name}.chr12.dmrpt --graph_file ${name}.chr12.dmgraph --verbose
else
        echo "Your VCF file, BED file, or both are empty."
fi

echo "Pipeline complete"
