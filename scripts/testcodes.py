#!/usr/bin/python
import sys
import _mypath
import numpy as np # need numpy for arrays and the like
#import subprocess # need to run binaries (blast, CD-hit, clustal)
import ConfigParser# need to read Config file
import time
from runbin import Command
from trimmer import trimmer
import os
from Bio import SeqIO, Entrez

def quitsy(given):
    sys.exit(given)
quitsy(0)
'''
RUNCLUSTAL = 'sleep aeu'  
print "Aligning sequences using Clustal Omega. This can take a few minutes"
start = time.time()
command = Command(RUNCLUSTAL)
CLUSTALOUT=command.run(timeout=3)
#command.run(timeout=900)
print command.status
if command.status == -15:
    print'BLAST TOOK TOO LONG! Try with a smaller evalue, or fewer maximum sequences.'
elif command.status != 0:
	print 'BLAST FAILED. I do not know why, check your inputs, internet connection, etc.'
print 'still running'
'''
