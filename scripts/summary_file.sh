#!/bin/bash

infile=''
outfile=''
outdir=''
chr=''

while getopts ":i:f:o:c:" option; do
   case "$option" in
      i)
        infile="$OPTARG";;
      f)
        outfile="$OPTARG";;
      o)
  	    outdir="$OPTARG";;
      c) 
        chr="$OPTARG";;
   esac
done

if [ "${chr}" != '' ]; then
	 	chr=".${chr}"
fi

awk '{print $3}' ${infile} | while read p;  do echo $p >> ${outfile}; cat ${outdir}/${p}/${p}${chr}.dmrpt >> ${outfile}; done
