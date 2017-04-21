import os
import sys
import csv
import configparser
from active.scripts import *

filebase=os.path.dirname(os.path.realpath(__file__))
consensusfilebase=filebase+'/consensus'
outfilebase=filebase+'/output'
if not os.path.exists(consensusfilebase):
    os.makedirs(consensusfilebase)
#Defining lists
info=[];pinfo=[]; csvlist=[]; hitnum_list=[]; hitting=[]
#Reading info
with open('Input.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        if len (row) != 0:
            csvlist = csvlist + [row]
csvfile.close()
select = csvlist[0][0]

config = ConfigParser.ConfigParser()
config.read("./uploads/CONFIGFILE")



#Opening files
with open(consensusfilebase+'//'+select+'.fasta.txt_mutations.txt', 'r') as consensus_file:
    for line in consensus_file:
        info.append(line.split())
with open(outfilebase+'//'+'simple_output.csv', 'r') as pdbresults:
    for line in pdbresults:
        pinfo.append(line.split(','))
#Generating search hit
for i in range(0, len(pinfo)):
    hitnum_list.append(pinfo[i][1])
#Writing output file
f = open(consensusfilebase+'//'+select+'.fasta.txt_mutations_trimmed.txt','w')
for r in range(0, len(info)):
    if info[r][2] not in hitnum_list:
        f.write(' '.join(info[r])+'\n')
    else:
        hitting.append(info[r][2])
f.write('Residue(s) removed: ')
for i in range(0, len(hitting)):
    f.write(hitting[i]+' ')
#Closing files
f.close()
consensus_file.close()
