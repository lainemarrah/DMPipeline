# DMPipeline
User-friendly pipeline to predict double minutes with DMFinder

**Prerequisites:**
* GCC (9.2.0+)
* SRAToolkit (2.10.5+)
* SAMTools (1.12+)
* Perl (5+)
* Python (3.7.16+)
* Conda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
* Picard: https://github.com/broadinstitute/picard/releases/tag/3.0.0
* Breakdancer: https://github.com/genome/breakdancer
* DMFinder: https://github.com/rmarduga/DMFinder

Scripts to download:
* convert_cns_to_bed.py: https://github.com/AmpliconSuite/AmpliconSuite-pipeline/blob/master/scripts/convert_cns_to_bed.py
* bam2cfg.pl: https://github.com/genome/breakdancer/blob/master/perl/bam2cfg.pl
* breakdancer2vcf.py: https://github.com/rmarduga/DMFinder/blob/master/tools/breakdancer2vcf.py

Make sure these scripts are located in the directory DMPipeline/scripts.

**Setup:**

Once all prerequisites listed above are downloaded, install DMPipeline by running:
```
git clone https://github.com/lainemarrah/DMPipeline
```

Next, configure the conda environments by running the code below. 
```
cd DMPipeline
conda env create -n bowtie2 --file envs/bowtie2-env.yml
conda env create -n cnvkit --file envs/cnvkit-env.yml
conda env create -n snakemake --file envs/snakemake-env.yml
```

**Basic Usage**

Run the general pipeline by running:
```
#example code
```
This will yield several files in your chosen output directory. These include: a sorted and indexed BAM file, a BED copy number file, a VCF structural variant file, and the DMFinder output files (which have suffixes .dmrpt and .dmgraph). 

You can then run the following to create a summary file of all desired samples, where the input file is the same as the sample input file used for the previous pipeline:
```
summary_file.sh [INPUT_FILE] [OUTPUT_FILE]
```
