# agotron_detector

## Introduction

Under development....

## Prerequisites

Ago-CLIP dataset, e.g. [PRJNA312501](https://www.ebi.ac.uk/ena/data/view/PRJNA312501)
[trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) (tested with v0.4.1) to process adaptor sequences off the raw reads

reference genome, e.g. [hg19](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/)
[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to align reads with reference genome

python (tested with v2.7.11) with modules: [pysam](https://github.com/pysam-developers/pysam), [numpy](https://github.com/numpy/numpy), and [mySQL](https://pypi.python.org/pypi/MySQL-python/)
[R](https://www.r-project.org) (tested with v3.2.5) with libraries: ggplot2, ggregex, dplyr, and optparse


## Scripts / Usage

The agotron_detector repository basically contains three scripts:

 1. agotron_coordinates.py, a python script to extract short introns from the UCSC mySQL server
    TODO...
 
 2. agotron_detector.py, a python script to intersect mapped reads with coordinates of interest (COI) and output agotron-relevant features
    TODO...
 
 3. agotron_visualizer.R, an R script that takes the output from detector.py and makes scatter and density plot
    TODO...
 



## Example

1. Download dataset, e.g. by using the following bash-script:
```bash
#!/bin/bash

array=( \
ftp.sra.ebi.ac.uk/vol1/fastq/SRR317/008/SRR3177718/SRR3177718.fastq.gz \
ftp.sra.ebi.ac.uk/vol1/fastq/SRR317/009/SRR3177719/SRR3177719.fastq.gz \
ftp.sra.ebi.ac.uk/vol1/fastq/SRR317/000/SRR3177720/SRR3177720.fastq.gz \
ftp.sra.ebi.ac.uk/vol1/fastq/SRR317/001/SRR3177721/SRR3177721.fastq.gz \
ftp.sra.ebi.ac.uk/vol1/fastq/SRR317/002/SRR3177722/SRR3177722.fastq.gz \
ftp.sra.ebi.ac.uk/vol1/fastq/SRR317/003/SRR3177723/SRR3177723.fastq.gz)

for i in "${array[@]}"
do
   wget $i   
done
```

2. Trim adapters:
```bash
trim_galore -A TCAGTCACTTCCAGC –length 18 *.fastq.gz
```

3. Map trimmed reads to reference genome

```bash
TODO
```

4. Extract short introns from UCSC mySQL database, and intersect coordinates with mapped bam-files

```bash
python agotron_coordinates.py | python agotron_detector -c agotron_coverage.txt -g /path/to/genome.fa -f *.bam > agotrons.bed
```

5. Visualize Extract short introns from UCSC mySQL database, and intersect coordinates with mapped bam-files

```bash
Rscript agotron_visualize.R 
```


## Citation

**Hansen TB. Agotron_detector. MiMB, 2017, in press**


## License

Copyright (C) 2017 ncRNALab.  See the [LICENSE](LICENSE.txt)
file for license rights and limitations (MIT).

