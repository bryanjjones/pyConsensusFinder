#!/usr/bin/python
import sys
import _mypath
import numpy as np # need numpy for arrays and the like
#import subprocess # need to run binaries (blast, CD-hit, clustal)
import ConfigParser# need to read Config file
import time
from runbin import Command
from config import updatefromconfig, updatebooleanconfig

#defaults
FILENAME= None
MAXIMUMSEQUENCES="5"#2000
BLASTEVALUE="1e-3"	
CONSENSUSTHRESHOLD="0.6"
USECOMPLETESEQUENCES="1"
ALIGNMENTITERATIONS="1"
MAXIMUMREDUNDANCYTHRESHOLD="0.9"
LOGGING="0"
KEEPTEMPFILES="0"

#def printf(x):
#	print x
'''
def updatefromconfig(cat, opt):
	if Config.has_option(cat, opt):
		return Config.get(cat, opt)
	else:
		return globals()[opt]
def updatebooleanconfig(cat, opt):
	if Config.has_option(cat, opt):
		return Config.getboolean(cat, opt)
	else:
		return globals()[opt]
'''
#read from the config file
Config = ConfigParser.ConfigParser()
#Config.read("./config/config.cfg")
configfile="./config/config.cfg"
#need to add a check to verify that the config.cfg file is there, otherwise results in a KeyError

FILENAME=updatefromconfig('BasicSettings', 'filename', configfile)
if not FILENAME:
#	print FILENAME
	print "No query file specified in CONFIGFILE. EXITING"
	sys.exit(0)
print FILENAME
#overwrite defaults if settings specified in CONFIGFILE
test = updatefromconfig('BlastSettings', 'MAXIMUMSEQUENCES', configfile)
print test
print globals()[globals()['test']]
'''
MAXIMUMSEQUENCES=updatefromconfig('BlastSettings', 'MAXIMUMSEQUENCES', configfile)
BLASTEVALUE=updatefromconfig('BlastSettings', 'BLASTEVALUE', configfile)
CONSENSUSTHRESHOLD=updatefromconfig('AlignmentSettings', 'ConsensusThreshold', configfile)
USECOMPLETESEQUENCES=updatebooleanconfig('AlignmentSettings', 'UseCompleteSequences', configfile)
ALIGNMENTITERATIONS=updatefromconfig('AlignmentSettings', 'AlignmentIterations', configfile)
MAXIMUMREDUNDANCYTHRESHOLD=updatefromconfig('AlignmentSettings', 'MaximumRedundancyThreshold', configfile)
LOGGING=updatebooleanconfig('TroubleShooting', 'Logging', configfile)
KEEPTEMPFILES=updatebooleanconfig('TroubleShooting', 'KeepTempFiles', configfile)

#location of binaries to call
BLAST="./binaries/blastp" #change this if another version is installed locally
CDHIT="./binaries/cd-hit" #change this if another version is installed locally
CLUSTAL="./binaries/clustalo-1.2.0-Ubuntu-x86_64" #change this if another version is installed locally


#Run Blast
RUNBLAST = BLAST+' -db nr -query ./uploads/'+FILENAME+' -evalue '+BLASTEVALUE+' -max_target_seqs '+MAXIMUMSEQUENCES+' -outfmt "6 sacc sseq pident" -remote' 
print RUNBLAST
'''
'''
start = time.time()
command = Command(RUNBLAST)
BLASTOUT=command.run(timeout=900)
end = time.time()
print 'Time elapsed (seconds)'
print(end - start)
print 'BLAST Results'
print BLASTOUT[1]
'''
# BLASTOUT[1] needs to become an array
# when complete sequences =1, curl BLASTOUT[1] VERSIONS (first entry on line), ideally replace existing sequences with curled (complete sequences)



