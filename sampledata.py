import os
import sys

os.system('mkdir sampleData')
fasterqDump2 = os.system('fasterq-dump --split-files SRR8185310 -O ./sampleData'])
fasterqDump3 = os.system('fasterq-dump --split-files SRR1411276 -O ./sampleData'])
