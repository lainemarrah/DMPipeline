#!/bin/bash

infile=$1
outfile=$2

awk '{print $2}' ${infile} |  while read p;  do echo $p >> ${outfile}; cat /scratch/lm2ku/dmfinder/$p/$p.chr12.dmrpt >> ${outfile}; done
