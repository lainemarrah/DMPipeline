# DMPipeline
User-friendly pipeline to predict double minutes with DMFinder

**Prerequisites**

Your sequencing data should be paired-end, whole genome sequencing.

Tools to install:
* GCC (9.2.0+)
* SAMTools (1.12+)
* Perl (5+)
* Python (3.7.16+)
* Java (17+)
* Conda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
* Breakdancer: https://github.com/genome/breakdancer
* DMFinder: https://github.com/rmarduga/DMFinder

All of the above should be located in your $PATH.

Scripts to download:
* picard.jar: https://github.com/broadinstitute/picard/releases/tag/3.0.0
* convert_cns_to_bed.py: https://github.com/AmpliconSuite/AmpliconSuite-pipeline/blob/master/scripts/convert_cns_to_bed.py
* bam2cfg.pl: https://github.com/genome/breakdancer/blob/master/perl/bam2cfg.pl
* breakdancer2vcf.py: https://github.com/rmarduga/DMFinder/blob/master/tools/breakdancer2vcf.py

All of the above scripts should be located in the directory DMPipeline/bin and should be executable (may need to run chmod u+x [SCRIPT NAME]). An easy way to do this for the GitHub-based scripts is to copy and paste with vim, or simply download and move. 

**Setup**

Once all prerequisites listed above are downloaded, install DMPipeline by running:
```
git clone https://github.com/lainemarrah/DMPipeline
```

Next, you will need to find the file libstdc++.so.6 in your conda setup. It is usually in /home/user/conda/lib/, but change the filepath in the below code if it is elsewhere. Then, run the following code to configure the conda environments and set the LD_LIBRARY_PATH environment variable:
```
cd DMPipeline
conda env create -n bowtie2 --file envs/bowtie2-env.yml
conda env create -n cnvkit --file envs/cnvkit-env.yml
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/user/conda/lib/
```

Download the hg19 reference genome by using the below code. The bowtie indexing script (align_hg19.sh) will also create a bowtie index that will help you align files later, but note this takes a significant amount of time so you may prefer to run this on HPC if possible.
```
cd DMPipeline
wget https://datasets.genepattern.org/data/module_support_files/AmpliconArchitect/hg19.tar.gz
tar -zxf hg19.tar.gz
rm hg19.tar.gz
scripts/align_hg19.sh
```

**Usage**

To run an individual sample, run: 
```
scripts/pipeline.sh [-h] -1 [FASTQ_FILE1] -2 [FASTQ_FILE2] -n [SAMPLE_NAME] -o [OUTPUT_DIRECTORY] [-p THREADS] [-c CHROMOSOME]
```
Options:
* 1:	First paired-end fastq file filepath.
* 2:	Second paired-end fastq file filepath.
* h:     Print Help.
* n:     The name of your sample. Output files will have this name.
* o:     Choose a directory for output files created by this pipeline. **If you are providing your own sorted and indexed BAM file**, it (as well as its index) should be in this directory.
* p:     [Optional] If you want to multithread, choose the number of threads. This will default to 0, but multithreading is recommended, especially if you are not providing a sorted BAM file.
* c:     [Optional] If you would like to limit your DM search to a single chromosome, specify here. The format should look like 'chr12'. This is recommended if you are targeting a specific gene, as it makes the operation much faster.


To run several samples on UVA Rivanna, run:
```
scripts/loop.sh [INPUT_FILE] -o [OUTPUT_DIRECTORY] [-p THREADS] [-c CHROMOSOME]
```
This input file should be tab-separated, with the first column being the first fastq file, the second column being the second fastq file, and the third column being the sample name. Within the output directory, a separate directory will be created for each sample. This will submit a Rivanna job for each of the provided samples. Threads and Chromosome fields are the same as the above.

**Output**

In your chosen output directory, the pipeline will output eight files. These include: a sorted and indexed BAM file, a sorted and indexed BAM file split to only chromosome 12, a BED copy number file for chromosome 12, a VCF structural variant file for chromosome 12, and the DMFinder output files for chromosome 12 (which have suffixes .dmrpt and .dmgraph). 

Once all the samples have run successfully, you can then run the following to create a basic summary file, where the input file is the same as the sample input file used for the previous script loop.sh:
```
cd [OUTPUT_DIRECTORY]
scripts/summary_file.sh [INPUT_FILE] [OUTPUT_FILE]
```
