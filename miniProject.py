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

#PS. Sorry Dr. Putonti! Busy past few weeks and I just didn't have time to troubleshoot problems 8 and 9. 

#This allows us to run the commandline from python
import os
import sys
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastpCommandline
import csv


#identify the current working directory
currDir = os.getcwd()
#initialize the results folder and log file
makeResults = os.system('mkdir results')
outfile = open(currDir + '/results/miniproject.log', 'w')

#print(currDir)
#workingDirectory  = os.path.realpath(sys.argv[0])

sraPath = currDir + '/sratoolkit.3.0.0-mac64/bin/'

#The first tool we are going to use is SRAtoolkit
fasterQDump = os.system('fasterq-dump --split-files SRR8185310 -O ./results')
#hardcode
#fasterQDump = os.system(sraPath+'fasterq-dump --split-files SRR8185310 -O ./results')

'''
2. Using SPAdes, assemble the genome. Write the SPAdes command to the log file.
This step takes a long time. My personal computer has 8 processors so when I thread these commands I will use 7 of them

'''
spadesPath = currDir + '/SPAdes-3.15.4-Darwin/bin'
spadesCom = spadesPath + '/spades.py'
mySpadesCom = 'spades.py -o ./results -s SRR8185310.fastq -t 7'
#hardcode
#assembly = subprocess.run([spadesCom, '-o', './results', '-s', 'SRR8185310.fastq', '-t', ' 7'])
#realcode 
assembly = os.system('spades.py -o ./results -s SRR8185310.fastq -t 7')
outfile.write('Running SPAdes.... ')
outfile.write(mySpadesCom)
outfile.write('\n')


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
outfile.write(output3)
outfile.write('\n')


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
outfile.write(output4)
outfile.write('\n')  

'''
5. Use GeneMarkS-2 to predict coding regions.
Using GeneMarkS-2, output the predicted protein sequences for the identified genes.
I'm assuming output results in results folder?
Look at readME file for parameter options.
'''

gmsPath = currDir+'/gms2_macos/'
#Use largeContigs.fasta
predictionCommand = gmsPath + 'gms2.pl'

#hardcode
#prediction = subprocess.run([predictionCommand,'-seq', './results/largeContigs.fasta', 
     #                        '--genome-type', 'auto', '--gcode', 'auto', '--faa', 'results/proSeq'])

#realcode
prdCom = 'gms2.pl -seq ./results/largeContigs.fasta --genome-type auto --gcode auto --faa results/proSeq'
prediction = os.system(prdCom)
#we will use the proSeq file in the results folder for the next question
outfile.write('Running GeneMarkS2.... ')
outfile.write(prdCom)
outfile.write('\n') 

'''
6. Query your predicted amino acid sequences (from #5) against these sequences to predict their function.    
Query via BLAST your predicted amino acid sequences in order to identify the predicted function of each. 
BLAST+ (and all its functions) are accessible to you already.
Produce an output file called “predicted_functionality.csv”. It should be a csv file formatted as follows: 
    Query sequence ID, Subject sequence ID, % Identity, % Query Coverage
    
Note, only report the best hit for each of your predicted sequences.
'''
#help(NCBIWWW.qblast)


blastDBPath = currDir+'/ncbi-blast-2.12.0+/bin/makeblastdb'

#first make a database with the ecoliFASTA file
#hardcode
#makeDB = subprocess.run([blastDBPath, '-in', 'Ecoli.fasta', '-out', './results/ecoliRef', 
#                         '-title', 'ecoliRef', '-dbtype', 'prot'])
#realcode
makeDB = os.system('makeblastdb -in Ecoli.fasta -out ./results/ecoliRef -title ecoliRef -dbtype prot')


blastQpath = currDir+'/ncbi-blast-2.12.0+/bin/'



#Blast our proSeq results against the database that we just built
query = 'blastp -query ./results/proSeq -db ./results/ecoliRef -out ./results/predicted_functionality.csv -outfmt "10 qseqid sseqid pident qcovs" -max_target_seqs 1 '

outfile.write('Running BLAST.... ')
outfile.write(query)
outfile.write('\n') 
#hardcode
#blastQuery2 = os.system(blastQpath+query)
#realCode
blastQuery2 = os.system(query)
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

output7 = "GeneMarkS found "+ str(count-4140)+ " additional CDS than the RefSeq."
print(output7)
outfile.write(output7)
outfile.write('\n')



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
#fetch the K-12 derivative

#hardcode
#fetchSRA4 = os.system(fqD+ ' --split-files SRR1411276 -O ./results')
#realcode
fetchSRA2 = os.system('fasterq-dump --split-files SRR1411276 -O ./results')



#Here we fetch the older E. coli reference genome (NC_000913) and save it to a fasta file to be used later in tophat2
from Bio import Entrez

Entrez.email = 'qthomas@luc.edu'    #let ncbi know who is accessing
handle = Entrez.efetch(db = 'nucleotide', id = 'NC_000913', rettype="fasta")
#print(handle)
header = handle.readline()
#seqs = handle.read().strip()
#print(seqs)


outfile2 = open(currDir + '/results/annotated.fasta', 'w')
outfile2.write(header)
for line in handle:
    outfile2.write(line)
outfile2.close()

#build the index 
os.system('bowtie2-build ./results/annotated.fasta ./results/annotatedRef')
#hardcode
#os.system('./bowtie2-2.4.5-macos-arm64/bowtie2-build ./results/annotated.fasta ./results/annotatedRef')


#run tophat2
#hardcode
#os.system('./tophat-2.1.1.OSX_x86_64/tophat2 --no-novel-juncs -p 8 ./results/annotatedRef ./results/SRR1411276.fastq')
#realCode
os.system('tophat2 --no-novel-juncs -p 8 ./results/annotatedRef ./results/SRR1411276.fastq')

#The output for tophat2 is accepted_hits.bam (a sam file) which we will use as input for our cufflinks call
os.system('cp ./tophat_out/accepted_hits.bam ./results/')
#We will now run cufflinks to generate a transcriptome assembly for each condition
os.system('cufflinks ./results/accepted_hits.bam -p 8')

#Let's copy this output to our results folder
os.system('cp ./transcripts.gtf ./results/')

'''
9. Parse through the Cufflinks output to create a file called “transcriptome_data.fpkm”. 
In this CSV format file, write the seqname, start, end, strand, and FPKM for each record 
in the Cufflinks output file
'''

#FPKM: Fragments Per kb of transcript per Million mapped reads

#create output file
outfile2 = open(currDir +'/results/transcriptome_data.fpkm', 'w')
readGTF = open(currDir+ '/results/transcripts.gtf').read().split('\n')
'''
from BCBio import GFF

limit_info = dict(gff_id=["chr1"], gff_source=["Coding_transcript"])


for rec in GFF.parse(readGTF):
    print(rec.features[0])
readGTF.close()
'''


with open(currDir +'/results/transcriptome_data.fpkm', 'w') as f:
    # create the csv writer
    writer = csv.writer(f)
    count=2
    
    for read in readGTF:
        #print(read)
        line = read.split("\t")
        #print(len(line))
        if len(line) <9:
            break
        lastBit = line[8].split(' ')
        #print(lastBit)
        #newLine = [line[0]+line[3]+line[4]+line[6]+lastBit[7]+lastBit[8]]
        if count%2== 0:
            newLine = [line[0], line[3],line[4],line[6],lastBit[6],lastBit[7]]
        elif count%2==1:
            newLine = [line[0], line[3],line[4],line[6],lastBit[3],lastBit[4]]
            print('y')
        count+=1
        writer.writerow(newLine)


outfile.close()
