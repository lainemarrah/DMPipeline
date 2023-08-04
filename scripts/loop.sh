#!/bin/bash

infile=$1

outdir=''
threads=0
chr=""

while getopts ":o:p:" option; do
   case "$option" in
      o)
  	     outdir="$OPTARG";;
      p)
        threads="$OPTARG";;
      c) 
        chr="$OPTARG";;
   esac
done

while read -r field1 field2 field3; do mkdir ${outdir}/${field3}; sbatch --error ${field3}.err --output ${field3}.out --job-name=${field3} pipeline.sh -1 ${field1} -2 ${field2} -n ${field3} -o ${outdir}/${field3} -p ${threads} -c ${chr}; done < ${infile}
#echo "Process complete"
