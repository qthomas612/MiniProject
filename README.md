# MiniProject
This project is for Comp483

# Description:

This program will retrieve the reads from the resequencing of the K-12 E. coli strain and assemble its genome using SPAdes. After some minor filtering, this program then uses GeneMarkS2 to predict coding regions from our assembly. These coding regions are blasted against an established E coli database to annotate the genome. Coding regions are compared to the original E. coli sequencing project.  



# Prerequisites: 

This python wrapper uses publicly available software tools; you must download and install these tools in your home directory and export them to your path. 

SRA-toolkit 
-https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump

SPADES 
-https://cab.spbu.ru/files/release3.15.4/manual.html#sec3.1

GeneMarkS-2 
-http://exon.gatech.edu/GeneMark/license_download.cgi

TopHat2
-https://ccb.jhu.edu/software/tophat/manual.shtml
* It is worth noting that when I ran this on my Mac: 
* I had to edit the first line in `tophat` from `#!/usr/bin/env python` to `#!/usr/bin/env python2`
* The most recent version of tophat2 did NOT work, I downloaded and installed version 2.0.11

Bowtie2
-http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

CuffLinks
-http://cole-trapnell-lab.github.io/cufflinks/install/

Export path to each binary in your bash_profile.

Ex. `export PATH="your-dir:$PATH"`

See https://gist.github.com/nex3/c395b2f8fd4b02068be37c961301caa7 for exporting path in different operating systems.


# Output:
A results folder will be created in whichever directory you are running the program from. Below is a list of files you should see in your output directory:
1.  annotated.fasta (original Ecoli genome)
2.  annotatedRef.*.bt2 (bowtie2 index)
3.  ecoliRef files (generated from making blast database)
4.  all SPAdes output, most notably scaffolds.fasta and contigs.fasta
5.  largeContigs.fasta (filtered contigs)
6.  proSeq (geneMarkS2 coding regions)
7.  results.log (all relevant and important information regarding the pipeline will be output here)
8.  tophat_out (directory for outputfiles of tophat2)

# Sample Data
Sample SRA files were too large to be upladed to github. Instead run the sample.py file to fetch SRA files and download them into a folder in your working directory called sampleData. Ecoli.fasta is a required dataset for the program put into the results folder. The program will detect the results folder being present and will not overwrite if you have it in your working directory.
