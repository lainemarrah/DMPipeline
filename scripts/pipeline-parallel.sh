#!/bin/bash

Help()
{  
   #show help
   echo
   echo "Syntax: pipeline.sh [-h] -i [INPUT_FILE] -m [MEMORY] -o [OUTPUT_DIRECTORY] [-p THREADS] [-c CHROMOSOME]"
   echo "options:"  
   echo "h     Print this Help."
   echo "i     Input text file with desired samples formatted as (SRA accession number)/(desired sample name), e.g. SRRXXXX/SAMPLE."
   echo "m	   Desired total memory for this operation in MB. Default is 25000 (or 25GB)."
   echo "o     Choose a directory for output files created by this pipeline."
   echo "p     Number of threads to use. The memory will automatically divide to maximize efficiency." 
   echo "c     If you would like to limit your DM search to a single chromosome, specify here. The format should look like 'chr12' This is recommended if interested in a specific gene."
}

infile=''
mem='20000'
outdir=''
chr=""
threads=0
pldir=$PWD

while getopts ":hi:m:o:p:c:" option; do
   case "$option" in
      h) #display Help
         Help
         exit;;
      i) #assign input file
      	 infile="$OPTARG";;
      m) #assign total memory
         mem="$OPTARG";;
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


#add error catching for infile format

source ~/.bashrc

#making variables for later
rows=$(wc -l $infile | cut -d " " -f 1)
cpuperfile=$(($cpu / $rows))
memperfile=$(($mem / $rows - 2000))
sortmem=$(($memperfile / 1000))"G"
#check chr period consistency
samples=$(awk -v chr=${chr} -v dir=${outdir} -F "/" '{print dir $2 "/" $2 chr ".sorted"}' ${infile})

#making directories for each sample
while read -r accession name; do
       if [ ! -d ${outdir}/${name} ]; then
               mkdir ${outdir}/${name}
       fi
       if [ ! -d ${outdir}/${name}/fastq ]; then
               mkdir ${outdir}/${name}/fastq
        fi     
done < ${infile}

echo "Retrieving fastqs and aligning a BAM"

#add something for searching for fastqs so we can skip download step if there is already a directory/available fastq

##downloading fastq
#fastp step? adds like an hour so im not sure
while read -r accession name; do
       if [ ! -f ${outdir}/${name}/fastq/${id} ]; then
               module load sratoolkit/2.10.5
               prefetch --max-size 100G $id
               fastq-dump --outdir fastq --gzip --skip-technical  -F  --split-files --clip ./$id
               rm -rf $id       
               fastp -z 1 -i  /project/UVACC_public/CCLE/CCLE_WGS/${id}_1.fq.gz -I  /project/UVACC_public/CCLE/CCLE_WGS/${id}_2.fq.gz -o /scratch/lm2ku/fastp/${id}_1.fastq.gz -O /scratch/lm2ku/fastp/${id}_2.fastq.gz
               if [ -d fastq -a -f /scratch/lm2ku/fastp/${id}_1.fastq.gz -a -f /scratch/lm2ku/fastp/${id}_2.fastq.gz ]; then
                      rm -r fastq
               fi
       fi  
done < ${infile}


conda activate bowtie2
bowtie2  -p ${threads} -x hg19/hg19full -1 ${fastq1} -2 ${fastq2} | samtools view -bS > ${outdir}/${name}.bam
conda deactivate

parallel -a ${infile} 'name={/}; cd ${outdir}/${name}; samtools fixmate -O bam ${name}.bam ${name}.fixmate; samtools sort -m 20G -o ${name}.sort ${name}.fixmate; samtools index ${name}.sort'
parallel -a ${infile} 'name={/}; cd ${outdir}/${name}; samtools fixmate -O bam ${name}.bam ${name}.fixmate'
#parallel -a ${infile} 'name={/}; cd ${outdir}/${name}; samtools sort -m 10G -o ${name}.sort ${name}.fixmate'
#have not tested below line but it should work
parallel -a ${infile} 'name={/}; cd ${outdir}/${name}; samtools sort -m '${sortmem}' -@ '${cpuperfile}' -o ${name}.sort ${name}.fixmate'
parallel -a ${infile} 'name={/}; cd ${outdir}/${name}; samtools index ${name}.sort'

#adjust chr parts post-split
parallel -a ${infile} 'name={/}; cd ${outdir}/${name}; samtools view ${name}.sort '${chr}' -b > ${name}.'${chr}
parallel -a ${infile} 'name={/}; cd ${outdir}/${name}; java -jar ~/picard.jar AddOrReplaceReadGroups I=${name}.'${chr}' O=${name}.'${chr}'.sorted SORT_ORDER=coordinate RGID='${id}' RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=${name}; rm ${name}.'${chr}'; samtools index ${name}.'${chr}.'sorted'
parallel -a ${infile} 'name={/}; cd ${outdir}/${name}; samtools flagstat ${name}.${chr}.sorted'

echo "Making BED file"
mkdir ${outdir}/cnvkit
conda activate cnvkit
cnvkit.py batch ${samples} -r hg19/hg19_cnvkit_filtered_ref.cnn -p ${threads} -d ${outdir}/cnvkit
conda deactivate
while IFS="/" read accession name; do
  scripts/convert_cns_to_bed.py --cns_file=${outdir}/cnvkit/${name}.chr12.cns
  cp ${outdir}/cnvkit/${name}.${chr}_CNV_CALLS.bed ${outdir}/${name}/${name}.${chr}.bed 
done < ${infile}

echo "Making VCF file/s"
while IFS="/" read accession name; do
       ~/scripts/bam2cfg.pl ${outdir}/${name}/${name}.${chr}.sorted > ${outdir}/${name}/${name}.${chr}.cfg
done < ${infile}
parallel -a ${infile} 'name={/}; cd ${outdir}/${name}; breakdancer-max ${outdir}/${name}/${name}.'${chr}'.cfg | tail -n +5 > ${outdir}/${name}/${name}.'${chr}'.brk'
while IFS="/" read accession name; do
       bin/breakdancer2vcf.py -i ${outdir}/${name}/${name}.${chr}.brk -o ${outdir}/${name}/${name}.${chr}.vcf
done < ${infile}

if [ -s ${outdir}/${name}/${name}.${chr}.vcf -a -s ${outdir}/${name}/${name}.${chr}.bed ]; then
	dm_find.pl --input_bam ${outdir}/${name}/${name}.${chr}.sorted --sv ${outdir}/${name}/${name}.${chr}.vcf --cn ${outdir}/${name}/${name}.${chr}.bed --report_file ${outdir}/${name}/${name}.${chr}.dmrpt --graph_file ${outdir}/${name}/${name}.${chr}.dmgraph --verbose
else
        echo "Your VCF file, BED file, or both are empty."
fi

echo "Pipeline complete"
