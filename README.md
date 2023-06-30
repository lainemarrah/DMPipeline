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
* Bowtie2: https://github.com/BenLangmead/bowtie2 (necessary if you're not supplying your own aligned BAM file)
* Fastp: https://github.com/OpenGene/fastp (optional but recommended)

Scripts to download:
* convert_cns_to_bed.py: https://github.com/AmpliconSuite/AmpliconSuite-pipeline/blob/master/scripts/convert_cns_to_bed.py
* bam2cfg.pl: https://github.com/genome/breakdancer/blob/master/perl/bam2cfg.pl
* breakdancer2vcf.py: https://github.com/rmarduga/DMFinder/blob/master/tools/breakdancer2vcf.py

todo: figure out how to configure scripts
