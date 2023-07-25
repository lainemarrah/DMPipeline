#!/bin/bash

Help()
{  
   #show help
   echo
   echo "Syntax: pipeline.sh [-h] -1 [FASTQ_FILE1] -2 [FASTQ_FILE2] -n [SAMPLE_NAME] -o [OUTPUT_DIRECTORY] [-p THREADS] [-c CHROMOSOME]"
   echo "options:"  
   echo "1	First paired-end fastq file filepath."
   echo "2	Second paired-end fastq file filepath."
   echo "h     Print this Help."
   echo "n     The name of your sample. Output files will have this name."
   echo "o     Choose a directory for output files created by this pipeline. If you are providing your own sorted and indexed BAM file, it should be in this directory."
   echo "p     [Optional] If you want to multithread, choose the number of threads." 
   echo "c     [Optional] If you would like to limit your DM search to a single chromosome, specify here. The format should look like 'chr12'."
}

fastq1=''
fastq2=''
name=''
outdir=''
chr=""
threads=0
pldir=$PWD

while getopts ":h1:2:n:o:p:c:" option; do
   case "$option" in
      1) #assign fastq1
	 fastq1="$OPTARG";;
      2) #assign fastq2
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
      c) #chromosome
      	 chr="$OPTARG";;
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

source ~/.bashrc

echo "Making and sorting BAM file..."

#create sam -> bam -> sorted bam
if [ -f ${outdir}/${name}.sort ]; then
        echo "Sorted BAM already exists"
else
        conda activate bowtie2
        bowtie2  -p ${threads} -x hg19/hg19full -1 ${fastq1} -2 ${fastq2} | samtools view -bS > ${outdir}/${name}.bam
        conda deactivate

        samtools fixmate -O bam ${outdir}/${name}.bam ${outdir}/${name}.fixmate
        if [ -s ${outdir}/${name}.fixmate ]; then
                rm ${outdir}/${name}.bam
        fi

        samtools sort -m 20G -o ${outdir}/${name}.sort ${outdir}/${name}.fixmate
        if [ -s ${outdir}/${name}.sort ]; then
                rm ${outdir}/${name}.fixmate
        fi

        samtools index ${outdir}/${name}.sort
fi

if [ "${chr}" != "" ]; then
	if [ -f ${outdir}/${name}.${chr}.sorted ]; then
        	echo "Split BAM already exists"
	 	chr=".${chr}"
	else
        	samtools view ${outdir}/${name}.sort ${chr} -b > ${outdir}/${name}.${chr}
        	java -jar bin/picard.jar AddOrReplaceReadGroups I=${outdir}/${name}.${chr} O=${outdir}/${name}.${chr}.sorted SORT_ORDER=coordinate RGID=${id} RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=${name}
        	rm ${outdir}/${name}.${chr}
        	samtools index ${outdir}/${name}.${chr}.sorted
        	samtools flagstat ${outdir}/${name}.${chr}.sorted
	 	chr=".${chr}"
	fi
else
	java -jar bin/picard.jar AddOrReplaceReadGroups I=${outdir}/${name}.sort O=${outdir}/${name}.sorted SORT_ORDER=coordinate RGID=${id} RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=${name}
	samtools index ${outdir}/${name}.sorted
	samtools flagstat ${outdir}/${name}.sorted

fi

if [ -s ${outdir}/${name}.sorted ]; then
        rm ${outdir}/${name}.sort
fi

if [ ! -f ${outdir}/${name}${chr}.sorted.bai ]; then
        samtools index ${outdir}/${name}${chr}.sorted
        samtools flagstat ${outdir}/${name}${chr}.sorted
fi

echo "Making BED file..."
#create cn file
if [ -f ${outdir}/${name}${chr}.bed ]; then
	echo "BED file already exists"
else
	if [ ! -d ${outdir}/cnvkit ]; then
      		mkdir ${outdir}/cnvkit
	fi
	conda activate cnvkit
	cnvkit.py batch ${outdir}/${name}${chr}.sorted -r hg19/hg19_cnvkit_filtered_ref.cnn -p ${threads} -d ${outdir}/cnvkit
        cd ${outdir}/cnvkit
        ${pldir}/bin/convert_cns_to_bed.py --cns_file=${outdir}/cnvkit/${name}${chr}.cns
        cd ${pldir}
	cp ${outdir}/cnvkit/${name}${chr}_CNV_CALLS.bed ${outdir}/${name}${chr}.bed
	conda deactivate
fi

if [ -s ${outdir}/${name}${chr}.bed ]; then
	rm -r ${outdir}/cnvkit
fi

echo "Making VCF file..."
#create sv file
if [ -f ${outdir}/${name}${chr}.vcf ]; then
  echo "VCF file already exists"
else
	bin/bam2cfg.pl ${outdir}/${name}${chr}.sorted > ${outdir}/${name}${chr}.cfg
	breakdancer-max ${outdir}/${name}${chr}.cfg | tail -n +5 > ${outdir}/${name}${chr}.brk
	bin/breakdancer2vcf.py -i ${outdir}/${name}${chr}.brk -o ${outdir}/${name}${chr}.vcf
fi

if [ -s ${outdir}/${name}${chr}.vcf ]; then
	rm ${outdir}/${name}${chr}.cfg
 	rm ${outdir}/${name}${chr}.brk
fi

echo "DMFinder prerequisite files complete, running DMFinder..."

#actually run dmfinder
if [ -s ${outdir}/${name}${chr}.vcf -a -s ${outdir}/${name}${chr}.bed ]; then
	dm_find.pl --input_bam ${outdir}/${name}${chr}.sorted --sv ${outdir}/${name}${chr}.vcf --cn ${outdir}/${name}${chr}.bed --report_file ${outdir}/${name}${chr}.dmrpt --graph_file ${outdir}/${name}${chr}.dmgraph --verbose
else
        echo "Your VCF file, BED file, or both are empty."
fi

echo "Pipeline complete"
