#!/bin/bash

infile=$1
outfile=$2

awk '{print $3}' ${infile} |  while read p;  do echo $p >> ${outfile}; cat $p.chr12.dmrpt >> ${outfile}; done
