# DMPipeline
User-friendly pipeline to predict double minutes with DMFinder

Prerequisites:
* GCC (9.2.0+)
* SRAToolkit (2.10.5+)
* SAMTools (1.12+)
* Perl (5+)
* Python (3.7.16+)
* Picard: https://github.com/broadinstitute/picard/releases/tag/3.0.0
* Breakdancer: https://github.com/genome/breakdancer
* CNVkit: https://github.com/etal/cnvkit
* DMFinder: https://github.com/rmarduga/DMFinder
* Conda environment found

Scripts to download:
* convert_cns_to_bed.py: https://github.com/AmpliconSuite/AmpliconSuite-pipeline/blob/master/scripts/convert_cns_to_bed.py
* bam2cfg.pl: https://github.com/genome/breakdancer/blob/master/perl/bam2cfg.pl
* breakdancer2vcf.py: https://github.com/rmarduga/DMFinder/blob/master/tools/breakdancer2vcf.py

Quickstart:

Once all prerequisites are installed as well as DMPipeline, configure a Conda environment by running the code below to use the rest of the prerequisites. If you don't already have conda, see here: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
```
conda env create -n dmpipeline --file environment.yml
conda activate dmpipeline
```
Run the pipeline by running:
```
#example code
```
This will yield several files in your chosen output directory. These include: a sorted and indexed BAM file, a sorted and indexed BAM file split to only the 12th chromosome, a BED copy number file, a VCF structural variant file, and the DMFinder output files (which have suffixes .dmrpt and .dmgraph). 

You can then run the following to create a summary file of all desired samples, where the input file is the same as the sample input file used for the previous pipeline:
```
summary_file.sh [INPUT_FILE] [OUTPUT_FILE]
```
