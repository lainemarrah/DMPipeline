#!/bin/bash

Help()
{  
   #show help
   echo
   echo "Syntax: pipeline.sh [-h] [-i SRA_ID] [-n SAMPLE_NAME] [-o OUTPUT_DIRECTORY] [-p THREADS] [-p REFERENCE] [-t TEMP_DIRECTORY]"
   echo "options:"   
   echo "h     Print this Help."
   echo "i     The SRA ID of your sample. This should be something like SRRXXXXXXX."
   echo "n     The name of your sample. Output files will have this name."
   echo "o     Choose a directory for output files created by this pipeline. This directory should contain an indexed and sorted BAM file."
   echo "p     [Optional for BED] If you want to multithread, choose number of threads here. This is for the CNVkit portion of the pipeline." 
   echo "r     [Optional for BED] Reference genome .cnn file. This is for the CNVkit portion of the pipeline."
   #make sure i confirm where to get .cnn from cnvkit
   echo "t     Choose a directory for temporary files to be downloaded to."
   echo
}

bt=''
id=''
name=''
outdir=''
tempdir=''
threads=1

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
      r) #reference genome file
         ref="$OPTARG";;
      t) #assign temp dir
         tempdir="$OPTARG";;
     \?) #invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done

#error catch here: check for existence of sorted and indexed bam file

cd $outdir
samtools view $name.sorted chr12 -b > $name.chr12.sorted
samtools index $name.chr12.sorted
samtools flagstat $name.chr12.sorted

echo "Making BED file..."
#create cn file
if [ -f $name.chr12.bed ]; then
	echo "BED file already exists"
else
	if [ ! -d $tempdir/cnvkit/$name.chr12 ]; then
    mkdir $tempdir/cnvkit/$name.chr12
	fi
	cd $tempdir/cnvkit/$name.chr12
	cnvkit.py batch $outdir/$name.chr12.sorted -r $ref -p $threads -d $tempdir/cnvkit/$name.chr12
	convert_cns_to_bed.py --cns_file=$name.chr12.cns
	cp ${name}_CNV_CALLS.bed $outdir/$name.chr12.bed
fi

echo "Making VCF file..."
#create sv file
if [ -f ${name}.chr12.vcf ]; then
  echo "VCF file already exists"
else
  cd $tempdir 
	bam2cfg.pl $outdir/$name.chr12.sorted > $name.chr12.cfg
	breakdancer-max $name.chr12.cfg > $name.chr12.bd
	tail -n +5 $name.chr12.bd > $name.chr12.brk
	breakdancer2vcf.py -i $name.chr12.brk -o $outdir/$name.chr12.vcf
fi

echo "DMFinder prerequisite files complete, running DMFinder..."

cd $outdir
#actually run dmfinder
if [ -s $name.chr12.vcf -a -s $name.chr12.bed ]; then
	/home/lm2ku/DMFinder/dm_find.pl --input_bam $name.chr12.sorted --sv $name.chr12.vcf --cn $name.chr12.bed --report_file $name.chr12.dmrpt --graph_file $name.chr12.dmgraph --verbose
else
        echo "Your VCF file, BED file, or both are empty."
fi

echo "Pipeline complete"
