#!/bin/bash




date_start=$(date +"%s")
echo "========= Starting: `date` =========="

# -------------------------------------------
# Section 3.1: Download genome (hg19)
# -------------------------------------------

echo "Section 3.1..."

date_section=$(date +"%s")

#1
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz

#2
tar -zxvf chromFa.tar.gz &>/dev/null

#3
cat *.fa > hg19.fa 

#4
samtools faidx hg19.fa

#5
echo "Section 3.1... bowtie2-build (have patience)"
bowtie2-build hg19.fa hg19 &>/dev/null

duration=$(($(date +"%s")-$date_section))
echo "Section 3.1... $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."


# -------------------------------------------
# Section 3.2: 
# -------------------------------------------

echo "Section 3.2..."

date_section=$(date +"%s")

#6

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

duration=$(($(date +"%s")-$date_section))
echo "Section 3.2... $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."



# -------------------------------------------
# Section 3.3: 
# -------------------------------------------

echo "Section 3.3..."

date_section=$(date +"%s")

#8

trim_galore -A TCAGTCACTTCCAGC -length 18 *.fastq.gz

duration=$(($(date +"%s")-$date_section))
echo "Section 3.3... duration: $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."


# -------------------------------------------
# Section 3.4: 
# -------------------------------------------

echo "Section 3.4..."

#9-#12, mapping and indexing all trimmed fastq files

date_section=$(date +"%s")

for i in *_trimmed.fq.gz
do
    echo $i
    bowtie2 -q --local -x hg19 -U $i | samtools sort - > $i.sort.bam    
    samtools index $i.sort.bam
    
done

duration=$(($(date +"%s")-$date_section))
echo "Section 3.4... duration: $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."


# -------------------------------------------
# Section 3.5-3.7: 
# -------------------------------------------

echo "Section 3.5 - 3.7..."

date_section=$(date +"%s")

#13-15

python UCSC_intron_retriever.py | python analyzer.py -g hg19.fa | Rscript annotater.R 


duration=$(($(date +"%s")-$date_section))
echo "Section 3.5 - 3.7... duration: $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."


echo "========= Job finished at `date` =========="

duration=$(($(date +"%s")-$date_start))
echo "Overall duration: $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

