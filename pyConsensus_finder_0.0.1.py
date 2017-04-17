#!/usr/bin/python

import sys #need sys to use system variables
import numpy as np # need numpy for arrays and the like
import os # need os to be able to pull in variables like "SOURCE" and "THRESHOLD"

import subprocess # need to run binaries (blast, CD-hit, clustal)
import ConfigParser# need to read config file
# args = ['./binaries/blastp', '-db', 'nr', '-query', './uploads/Gc_hydrolase.fst', '-out', './processing/1BLAST.out', '-evalue', '0.1', '-max_target_seqs', '20', '-outfmt', '6 sacc sseq pident', '-remote']

Config = ConfigParser.ConfigParser()
Config

Config.read("./uploads/CONFIGFILE")

Config.sections()
Config.options('BlastSettings')
Config.get('BlastSettings', 'maximumsequences')

SOURCE=Config.get('BasicSettings', 'FileName')
MAXSEQ=Config.get('BlastSettings', 'maximumsequences')
EVALUE=Config.get('BlastSettings', 'blastevalue')
THRESHOLD=Config.get('AlignmentSettings', 'ConsensusThreshold')
USE_COMPLETE_SEQS=Config.getboolean('AlignmentSettings', 'UseCompleteSequences')
ALIGNMENTITER=Config.get('AlignmentSettings', 'AlignmentIterations')
ID_REDUNDANCY=Config.get('AlignmentSettings', 'MaximumRedundancyThreshold')
LOGGING=Config.getboolean('TroubleShooting', 'Logging')
KEEPTEMP=Config.getboolean('TroubleShooting', 'KeepTempFiles')

print '$BLAST -db nr -query ./processing/' + SOURCE + ' -out ./processing/1BLAST.out -evalue ' + EVALUE +' -max_target_seqs ' + MAXSEQ + ' -outfmt "6 sacc sseq pident" -remote'

#Run Blast
#popen.wait()
#args = ['./binaries/blastp', '-db', 'nr', '-query', './uploads/Gc_hydrolase.fst', '-evalue', '0.1', '-max_target_seqs', '20', '-outfmt', '6 sacc sseq pident', '-remote']
#popen.wait()
#output = popen.stdout.read()
#print output

