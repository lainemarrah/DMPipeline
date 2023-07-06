#!/bin/bash

infile=$1
while read -r field1 field2; do sbatch --error ${field2}.err --output ${field2}.out --job-name=${field2} pipeline.sh  ${field2} ${field1}; done < ${infile}
#echo "Process complete"
