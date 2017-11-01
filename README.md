# agotron_detector

## Introduction

The *agotron_detector* toolkit identifies and quantifies agotrons in Ago CLIPseq datasets.



## Prerequisites


[bowtie2](https://github.com/BenLangmead/bowtie2) (tested with v2.2.8) 

[samtools](https://github.com/samtools/) (tested with v1.3)

[python](https://www.python.org/) (tested with v2.7.11) including modules: [pysam](https://github.com/pysam-developers/pysam), [numpy](https://github.com/numpy/numpy), and [mySQL](https://github.com/arnaudsj/mysql-python)

[R](https://www.r-project.org) (tested with v3.2.5) including packages: ggplot2, ggregex, dplyr, tidyr, and optparse



## Scripts / Usage

The agotron_detector repository basically contains three scripts that should be run sequencially:



(1) **UCSC_intron_retriever.py**, a python script to extract short introns from the UCSC mySQL server    
 
``` python
Usage:
  UCSC_intron_retriever.py [ARGUMENTS] > [OUTPUT]
  
Options:
  -db <string>      MySQL database (default='hg19')
  -table <string>   MySQL table (default='refGene')
  -max <int>        Maximum intron length (default=150)
  -min <int>        Minimum intron length (default=50)
```

Alternatively, use **UCSC_mirna_retriever.py** or **tophat_intron_retriever.py** to retrieve RefSeq annotated miRNA coordinates or intron coordinates from tophat-produced *junctions.bed*.



 
(2) **analyzer.py**, a python script to intersect mapped reads with coordinates of interest and output agotron-relevant features    

``` python
Usage:
  analyzer.py [ARGUMENTS] < [INPUT] > [OUTPUT]
  
Required arguments:	
  -g <file>         Path to reference genome fastafile, must be indexed with samtools faidx
Optional arguments:	
  -f <files…>       Input bam-files (default=*.bam)
  -c <file>         Filename for coverage output (if empty, no coverage file is produced) 
                    (default='coverage.txt')
  -tr <int>         RPMM (reads per mapped million) expression threshold for output
                    (default = 5)
  -ts <int/'all'>   How many samples to meet RPMM expression threshold. For all samples, type ‘all’ 
                    (default = 2)
  -m <float>        Tolerance for reads mapping partly outside locus 
                    (default=0.1, e.g. 10% of the reads is allowed to map outside locus)
  -q <int>          Threshold for mapping quality 
                    (default=13)
  -a <int>          Add <int> flanking nucleotides (up and downstream) to the loci sequence output 
                    (default=10) 
``` 
 
 
 
(3) **annotater.R**, an R script that annotates agotron and outputs a few different plots
 
``` bash
Usage: 
  annotater.R [ARGUMENTS] < [INPUT]	

optional arguments:	
  -c <file>         Input coverage file (the -c output from analyzer.py) 
                    (default='coverage.txt')
  -p <string>       Prefix used in output files 
                    (default='agotron')
agotron definition:	
  -m <int>          Threshold for median read length 
                    (default=30)
  -h <float>        Minimum fraction of reads with distance (-d) between 5’end of read and 5’end of locus 
                    (default=0.7)
  -d <int>          Maximum distance allowed from predominant 5’end of reads to 5’end of locus 
                    (default=1)
```

## Example

Clone the repository and run the example script, *all_commands.sh*. 
This:
1. Downloads and prepares reference genome (hg19)
2. Downloads and trims Ago CLIPseq dataset (GSE78059)
3. Maps dataset to the reference genome using bowtie2
4. Intersects the mapped reads (bam-files) with annotated short introns to detect and annotate agotrons.

```bash
sh all_commands.sh
```

To annotate agotrons in mapped data (use -db and -g options to specify the reference genome used):

```bash
python UCSC_intron_retriever.py -db hg19 | python analyzer.py -g /path/to/hg19.fa -f /path/to/*.bam | Rscript annotater.R 
```




## Citation


**Hansen TB. Detecting agotrons in Ago CLIPseq data. MiMB, 2017, submitted**


## License

Copyright (C) 2017 ncRNALab.  See the [LICENSE](LICENSE.txt)
file for license rights and limitations (MIT).

