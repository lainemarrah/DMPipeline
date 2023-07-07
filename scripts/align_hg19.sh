#!/bin/bash

conda activate bowtie2

cd DMPipeline
bowtie2-build --large-index hg19/hg19full.fa hg19/hg19full

conda deactivate
