#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 20 13:10:44 2022

@author: QuinnThomas
"""
'''
1. Retrieve the Illumina reads for the resequencing of K-12 project:
https://www.ncbi.nlm.nih.gov/sra/SRX5005282. These are single-end Illumina reads. 
Your code will retrieve this file; i.e., you cannot just retrieve it.
'''

#This allows us to run the commandline from python
import os
import sys
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastpCommandline



home_dir = os.system("cd ~")
currDir = os.getcwd()
#print(currDir)
#workingDirectory  = os.path.realpath(sys.argv[0])

sraPath = currDir + '/sratoolkit.3.0.0-mac64/bin/'

#The first tool we are going to use is SRAtoolkit
#fetchCommand = 'prefetch -v SRR8185310'
import subprocess
#list_files = subprocess.run(["ls", "-l"])
#print("The exit code was: %d" % list_files.returncode)

#prefetch = subprocess.run(["prefetch", "-v", "SRR8185310"])
fqD = sraPath + 'fasterq-dump'
#fasterqDump = subprocess.run(fqD)
#fasterqDump2 = subprocess.run([fqD, '--split-files', 'SRR8185310', '-O', '/results'])

#subprocess.call('prefetch')
#convertSRA = subprocess.run(['fastq-dump', '--outdir', home_dir, '--split-files', '/home/QuinnThomas/ncbi/public/sra/SRR925811.sra'])


'''
2. Using SPAdes, assemble the genome. Write the SPAdes command to the log file.
This step takes a long time. My personal computer has 8 processors so when I thread these commands I will use 7 of them

'''
spadesPath = currDir + '/SPAdes-3.15.4-Darwin/bin'
spadesCom = spadesPath + '/spades.py'
#assembly = subprocess.run([spadesCom, '-o', './results', '-s', 'SRR8185310.fastq', '-t', ' 7'])


'''
3. Calculate the number of contigs with a length > 1000 and write the # out to the log file as follows:
There are # contigs > 1000 in the assembly.
From here on out, you will only consider those contigs > 1000 bp in length.
One of our outputs from spades is the contigs.fasta file. This is what we will read into python to learn more about the results.
'''

seqs=[]
records=[]

#save each fasta record into lists identified above
for seq_record in SeqIO.parse("results/contigs.fasta", "fasta"):
    #headers.append(seq_record.id)
    if len(seq_record.seq) > 1000:
        seqs.append(str(seq_record.seq))
        records.append(seq_record)

#we need the fasta version of this so that we can use it for geneMarkS-2
SeqIO.write(records, "results/largeContigs.fasta", "fasta")
    
        
#print(len(seqs))

largeContigs = len(seqs)
output3 = 'There are ' + str(largeContigs) + ' contigs > 1000 bp in the assembly.'
print(output3)

'''
4. Calculate the length of the assembly (the total number of bp in all of the contigs > 1000 bp in length) 
and write this # out to the log file as follows:
There are # bp in the assembly.
'''

totalLength = 0
for seq in seqs:
    totalLength+=len(seq)

output4 = 'There are ' +str(totalLength)+ ' bp in the assembly.'  
print(output4)  
    
'''
5. Use GeneMarkS-2 to predict coding regions.
Using GeneMarkS-2, output the predicted protein sequences for the identified genes.
I'm assuming output results in results folder?
Look at readME file for parameter options.
'''

gmsPath = currDir+'/gms2_macos/'
#Use largeContigs.fasta
predictionCommand = gmsPath + 'gms2.pl'
#prediction = subprocess.run([predictionCommand,'-seq', './results/largeContigs.fasta', 
     #                        '--genome-type', 'auto', '--gcode', 'auto', '--faa', 'results/proSeq'])

#we will use the proSeq file in the results folder for the next question

'''
6. Query your predicted amino acid sequences (from #5) against these sequences to predict their function.    
Query via BLAST your predicted amino acid sequences in order to identify the predicted function of each. 
BLAST+ (and all its functions) are accessible to you already.
Produce an output file called “predicted_functionality.csv”. It should be a csv file formatted as follows: 
    Query sequence ID, Subject sequence ID, % Identity, % Query Coverage
    
Note, only report the best hit for each of your predicted sequences.
'''
#help(NCBIWWW.qblast)

#first make a database with the ecoliFASTA file
blastDBPath = currDir+'/ncbi-blast-2.12.0+/bin/makeblastdb'
#makeblastdb -in Ecoli.fasta -out ecoliRef -title ecoliRef -dbtype prot

#makeDB = subprocess.run([blastDBPath, '-in', 'Ecoli.fasta', '-out', './results/ecoliRef', 
#                         '-title', 'ecoliRef', '-dbtype', 'prot'])


blastQpath = currDir+'/ncbi-blast-2.12.0+/bin/'
NcbiblastpCommandline()
'''
blastQuery = NcbiblastpCommandline(query="./results/proSeq", db="./results/ecoliRef", 
                                   out="./results/predicted_functionality", 
                                   outfmt='"10 qseqid sseqid pident qcovs"')
'''
query = 'blastp -query ./results/proSeq -db ./results/ecoliRef -out ./results/predicted_functionality.csv -outfmt "10 qseqid sseqid pident qcovs" -max_target_seqs 1 '
#blastQuery2 = os.system(blastQpath+query)
#print(blastQuery2)


'''
7. The assembled genome in RefSeq for E. coli K-12 (NC_000913) has 4140 CDS and 89 tRNAs annotated. 
Write to the log file the discrepancy (if any) found. 
For instance, if my GeneMark annotation predicted 4315 CDS, I would write,
GeneMarkS found 175 additional CDS than the RefSeq.
'''

count=0
for seq_record in SeqIO.parse("results/proSeq", "fasta"):
    count+=1
    
#print(count)

print("GeneMarkS found "+ str(count-4140)+ " additional CDS than the RefSeq.")



'''
8. Use TopHat & Cufflinks to map the reads of the E. coli transcriptome project of a 
K-12 derivative BW38028 (https://www.ncbi.nlm.nih.gov/sra/SRX604287) and quantify 
their expression, respectively. [You can use wget to retrieve this sequence; 
look at #1 to figure out how you construct the ftp address.] We’ll map these reads to 
the complete annotated genome NC_000913. [You do not need to do novel splice site discovery; 
include flag ‘--no-novel-juncs’.] Details re: the processes of TopHat and Cufflinks 
can be found in the class slides and Trapnell et al. 2013 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3334321/.
'''

#SRR1411276
#fetchSRA4 = os.system(fqD+ ' --split-files SRR1411276 -O ./results')
#fetchSRA2 = os.system('./fasterq-dump --split-files SRR1411276 -O ./results')
#fetchSRA3 = subprocess.run(['./fasterq-dump', '--split-files', 'SRR1411276', '-O', './results'])
#'SRR1411276.fastq'

#Here we fetch the older E. coli reference genome (NC_000913) and save it to a fasta file to be used later in tophat2
from Bio import Entrez
Entrez.email = 'qthomas@luc.edu'    #let ncbi know who is accessing
handle = Entrez.efetch(db = 'nucleotide', id = 'NC_000913', rettype="fasta")
#print(handle)
header = handle.readline()
#seqs = handle.read().strip()
#print(seqs)

'''
outfile2 = open(currDir + '/results/annotated.fasta', 'w')
outfile2.write(header)
for line in handle:
    outfile2.write(line)
outfile2.close()
'''
#build the index 
#os.system('bowtie2-build ./results/annotated.fasta ./results/annotatedRef')
#hardcode
#os.system('./bowtie2-2.4.5-macos-arm64/bowtie2-build ./results/annotated.fasta ./results/annotatedRef')


#run tophat2
#os.system('tophat2 ./results/annotatedRef ./results/SRR1411276.fq')
#hardcode
os.system('./tophat-2.1.1.OSX_x86_64/tophat2 --no-novel-juncs ./results/annotatedRef ./results/SRR1411276.fastq')
#realCode
#os.system('tophat2 --no-novel-juncs -p 8 ./results/annotatedRef ./results/SRR1411276.fastq')


