#!/bin/bash

infile=''
outdir=''
threads=0
chr=''

while getopts ":i:o:p:c:" option; do
   case "$option" in
      i)
        infile="$OPTARG";;
      o)
  	     outdir="$OPTARG";;
      p)
        threads="$OPTARG";;
      c) 
        chr="$OPTARG";;
     \?)
        echo "Invalid option"
        exit;;
   esac
done

while read -r field1 field2 field3; do sbatch --error ${outdir}/${field3}/${field3}.err --output ${outdir}/${field3}/${field3}.out --job-name=${field3} scripts/pipeline.sh -1 ${field1} -2 ${field2} -n ${field3} -o ${outdir}/${field3} -p ${threads} -c ${chr}; done < ${infile}
#echo "Process complete"
