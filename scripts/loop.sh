#!/bin/bash

infile=$1

outdir=''
threads=0

while getopts ":o:p:" option; do
   case "$option" in
      o)
  	     outdir="$OPTARG";;
      t)
      	 threads="$OPTARG";;
   esac
done

while read -r field1 field2 field3; do sbatch --error ${field3}.err --output ${field3}.out --job-name=${field3} pipeline.sh -f1 ${field1} -f2 ${field2} -n ${field3} -o ${outdir} -p ${threads}; done < ${infile}
#echo "Process complete"
