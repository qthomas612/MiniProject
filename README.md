# MiniProject
This project is for Comp483

# Description:

This program will retrieve the reads from the resequencing of the K-12 E. coli strain and assemble its genome using SPAdes. After some minor filtering, this program then uses GeneMarkS2 to predict coding regions from our assembly. These coding regions are blasted against an established E coli database to annotate the genome. Coding regions are compared to the original E. coli sequencing project.  




# Prerequisites: This python wrapper uses publicly available software tools; you must download and install these tools in your home directory and export them to your path. 

SRA-toolkit 
-https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump

SPADES 
-https://cab.spbu.ru/files/release3.15.4/manual.html#sec3.1

GeneMarkS-2 
-http://exon.gatech.edu/GeneMark/license_download.cgi

TopHat2
-https://ccb.jhu.edu/software/tophat/manual.shtml

Bowtie2
-http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

CuffLinks
-http://cole-trapnell-lab.github.io/cufflinks/install/

Export path to each binary in your bash_profile.

Ex. `export PATH="your-dir:$PATH"`

See https://gist.github.com/nex3/c395b2f8fd4b02068be37c961301caa7 for exporting path in different operating systems.


# Output:
