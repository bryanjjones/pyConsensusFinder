#Program is written and developed by Enoch Kan
#Kazlauskas Lab, Department of Biochemistry, Molecular Biology and Biophysics
#University of Minnesota, Twin Cities; Fall 2016
#Copyright (c) Enoch Kan 2016
#File OS
import os
import csv
import numpy as np
import pdb_res as pr
from Bio.PDB import PDBList
from Bio.PDB import PDBParser
filebase=os.path.dirname(os.path.realpath(__file__))
outfilebase=filebase+'/output'
if not os.path.exists(outfilebase):
    os.makedirs(outfilebase)
#Arrays
atr=[];val=[];csvlist=[];ares=[];resname=[];rescode=[];heteroname=[];heterocode=[];atomlist=[];atomnamelist=[];hitnum_list=[];hitting=[]
sort1=[];sort2=[];sort3=[];sort4=[];sort5=[]
#Reading info
with open('Input.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        if len (row) != 0:
            csvlist = csvlist + [row]
csvfile.close()


select = csvlist[0][0];achn=csvlist[1][0];mk=csvlist[3][0]
for i in range(0, len(csvlist[2])):
    ares.append(csvlist[2][i])

#Opening the file
pdbl = PDBList()
pdbl.retrieve_pdb_file(select,pdir='pdb')
file_path=filebase+'/pdb/pdb'+select+'.ent'
#Read the file
parser=PDBParser(QUIET=1)
structure=parser.get_structure('test', file_path)

reslist()
resinfo()



