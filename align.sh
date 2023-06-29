#!/bin/bash

#not sure if help/usage is necessary if i'm calling this from another script

Help()
{  
   #show help
   echo
   echo "Syntax: align.sh [-h] [-i SRA_ID] [-n SAMPLE_NAME] [-o OUTPUT_DIRECTORY] [-t TEMP_DIRECTORY]"
   echo "options:"   
   echo "h     Print this Help."
   echo "i     The SRA ID of your sample. This should be something like SRRXXXXXXX."
   echo "n     The name of your sample."
   echo "o     Choose a directory for output files created by this pipeline."
   echo "t     Choose a directory for temporary files to be downloaded to."
   echo
}

bt=''
id=''
name=''
outdir=''
tempdir=''
threads=1

while getopts "b:hi:n:o:t:p" option; do
   case "$option" in
      b) #bowtie index name
         bt="$OPTARG";;
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
      t) #assign temp dir
         tempdir="$OPTARG";;
     \?) #invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done

cd $tempdir

echo "Downloading and quality checking fastqs..."
if [ -f ${id}_1.fastq.gz -a -f ${id}_2.fastq.gz ]; then
	echo "Fastqs already exist"
else
  prefetch --max-size 100G $id
  fastq-dump --outdir $tempdir --gzip --skip-technical  -F  --split-files --clip ./$id
  rm -rf $id
  
	fastp -z 1 -i ${id}_1.fastq.gz -I ${id}_2.fastq.gz -o ${id}_1.fastp.gz -O ${id}_2.fastp.gz
fi
#delete unnecessary files
if [ -f $id_1.fastp.gz -a -f ${id}_2.fastp.gz ]; then
	rm /scratch/lm2ku/fastq/${id}_1.fastq.gz
  rm /scratch/lm2ku/fastq/${id}_2.fastq.gz
fi
echo "Fastq download and quality control complete"

echo "Making and sorting BAM file..."
#create sam -> bam -> sorted bam
#can make some changes to catch various file issues
if [ -f $name.sorted ]; then
	echo "Sorted BAM already exists"
else
  bowtie2  -p $threads -x /home/lm2ku/AA_data_repo/hg19/hg19full -1 /scratch/lm2ku/fastp/${id}_1.fastq.gz -2 /scratch/lm2ku/fastp/${id}_2.fastq.gz > ${name}.sam
	samtools fixmate -O bam $name.sam $name.fixmate	
	samtools sort $name.fixmate -m 20G -o $outdir/$name.sorted
fi

cd $outdir
if [ -f $name.sorted.bai ]; then
	echo "BAM index already exists"
else
  samtools index $name.sorted
fi

echo "Complete"
