#Sample SRA files were too large to be upladed to github. Instead run the sample.py file to fetch SRA files 
#and download them into a folder in your working directory called sampleData.

import os
import sys

os.system('mkdir sampleData')
fasterqDump2 = os.system('fasterq-dump --split-files SRR8185310 -O ./sampleData'])
fasterqDump3 = os.system('fasterq-dump --split-files SRR1411276 -O ./sampleData'])
