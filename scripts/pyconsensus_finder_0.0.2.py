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
Entrez.email = "bryanjjones@gmail.com"

#make temporary processing directory
#os.mkdir('./processing')

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

def cleanexit():
	if KEEPTEMPFILES:
		print "Exiting due to failure"
		print "Leaving temporary files"
		sys.exit("Exiting due to failure")
	else:
		print 'Exiting due to failure'
		print 'Deleting temporary files'
		os.remove('./processing/'+FILENAME+'_all_sequences.tmp')
		os.remove('./processing/cdhit.tmp')
		os.remove('./processing/cdhit.tmp.clstr')
		os.remove('./processing/clustal.tmp')
		sys.exit("Exiting due to failure")
		
#read from the config file
Config = ConfigParser.ConfigParser()
Config.read("./config/default.cfg") #read defaults from here
Config.read("./config/config.cfg") #overwrite defaults from here

#need to add a check to verify that the config.cfg file is there, otherwise results in a KeyError

FILENAME=Config.get('BasicSettings', 'filename')
if not FILENAME:
#	print FILENAME
	print "No query file specified in CONFIGFILE. EXITING"
	sys.exit("No query file specified in CONFIGFILE. EXITING")

#need to add a check to veriyf that file is protein sequence

MAXIMUMSEQUENCES=Config.get('BlastSettings', 'MAXIMUMSEQUENCES')
BLASTEVALUE=Config.get('BlastSettings', 'BLASTEVALUE')
CONSENSUSTHRESHOLD=Config.get('AlignmentSettings', 'ConsensusThreshold')
USECOMPLETESEQUENCES=Config.getboolean('AlignmentSettings', 'UseCompleteSequences')
ALIGNMENTITERATIONS=Config.get('AlignmentSettings', 'AlignmentIterations')
MAXIMUMREDUNDANCYTHRESHOLD=Config.get('AlignmentSettings', 'MaximumRedundancyThreshold')
LOGGING=Config.getboolean('TroubleShooting', 'Logging')
KEEPTEMPFILES=Config.getboolean('TroubleShooting', 'KeepTempFiles')

#location of binaries to call
BLAST="./binaries/blastp" #change this if another version is installed locally
CDHIT="./binaries/cd-hit" #change this if another version is installed locally
CLUSTAL="./binaries/clustalo-1.2.0-Ubuntu-x86_64" #change this if another version is installed locally


#Run Blast
RUNBLAST = BLAST+' -db nr -query ./uploads/'+FILENAME+' -evalue '+BLASTEVALUE+' -max_target_seqs '+MAXIMUMSEQUENCES+' -outfmt "6 sacc sseq pident" -remote' 
print "Begining BLAST search of NCBI. This will take a few minutes."
start = time.time()
command = Command(RUNBLAST)
BLASTOUT=command.run(timeout=900)

#add instructions for timeout (i.e. try again w/ fewer sequences)

#temporary data set to skip actual blast search
'''
BLASTOUT=[0,"""WP_083418319    PCPSFLIEHDRGLVLFDAGFDPRGLDDMAAYYPEISKALPMAGNRDLGIDRQLDGLGYRPSQISYVIPSHLHFDHAGGLYLFPDSTFLMGSGEMAFALQAHDKPQ---AGFFRVEDLLPTRHFDW--IETAHDFDLFGDGSVVLLFSPGHTPGSLALFVRLPNQ-NIILSGDVCHFPLEVDMGIIATSSFSPSYATF-ALRRLRMISKAWDARIWIQHEEDHWNEWPHAPE 26.840
WP_073456977    PMPAYLIEHPKGLVLFDTMLVPDAADDPERVYGPLAEHLGLKYTREQRVDNQIKALGYRLEDVTHVIASHTHFDHSGGLYLFPHAKFYVGEGELRFAL--WPDPAGAGFFRQADIEA--TRSFNW--VQVGFDHDLFGDGSVVVLHTPGHTPGELSLLVRL-KSRNFILTGDTVHLRQALEDEIPMP---YDSNTELAIRSIRRLKLLRESADATVWITHDPEDWAEFQHAPYCY       29.362
WP_049268273    ESYEIPVPWFLLTHPDGFTLIDGGLAVEGLKDPFAYWGGAVEQFKPVMPEEQGCL-EQLKRIGVAPEDIRYVILSHLHSDHTGAIGRFPHATHVVQRREYEYAFA--PDWFTSGAYCRYDYD---HPELNWFFLNGLSEDNYDLYGDGTLQCIFTPGHSPGHQSFLIRLSSGTNFTLAIDAAYTLDHYHERAL-PGLMTSATDVAQSVQKLRQLTERYNAILIPGHDPEEWEKIRLAPAWY 29.876
WP_025778256    LAVPIPTFLIQHEGGLLVFDTGLATDAAGDPARAYGPLAEAFDMSFPPEARIDTQLESLGFSTSDVTDVVLSHMHFDHTGGLELFPTARGFIGEGEL-----AYSRSPRRLDAAMYREEDIAAAGQIDWLEIPQGVDHDIFGDGSVVVLSMPGHTHGTLSLKLSPPDHRTIILTSDAAHLQSNIDETTGMPLD-----VDTRNKERSLRRLRLLASQPNTTVWANHDPDHWKQFRR      27.966
OCL02351        KLWLLHVGNLECDEAWFKRGGGTSTLSNPHPNRERRKLIMVSVLIEHPVEGLILFETGSGKDY--PE-IWGAPINDIFARVDYTEEQELDVQIKKTGHDIKDVKMVVIGHLHLDHAGGLEYFRNTGIPIYVHEKELKHAFYSVATKSDLGVY----LPHYLTFDLNW--VPFSGAFLEIAQGLN-LHHAPGHTPGLCILQVNMPKSGAWIFTTDMYHVSENFEESVPQGWLAREHDDWVRSNHMIHMLQRRTGARMVFGHCTKALEGLNMAPHAY       26.545"""]
'''
end = time.time()
print 'BLAST search took '+str(int(end - start))+' seconds'

#print 'BLAST Results'
BLASTOUT=BLASTOUT[1].splitlines()
for index in range(len(BLASTOUT[:])):
	BLASTOUT[index]=BLASTOUT[index].split()
VERSIONS = [item[0] for item in BLASTOUT]
print 'BLAST returned'+str(len(VERSIONS))+' sequences.'
#if USECOMPLETESEQUENCES:

handle = Entrez.efetch(db="sequences", id=",".join(VERSIONS), retmax=MAXIMUMSEQUENCES, rettype="fasta", retmode="text")

record = (handle.read()).splitlines()
np.savetxt(('./processing/'+FILENAME+'_all_sequences.tmp'),record,delimiter="",fmt="%s") 

RUNCDHIT = CDHIT+' -i ./processing/'+FILENAME+'_all_sequences.tmp -o ./processing/cdhit.tmp -c '+MAXIMUMREDUNDANCYTHRESHOLD+' -M 0 -T 0'  
print "Removing redundant sequences using CD-HIT"
start = time.time()
command = Command(RUNCDHIT)
command.run(timeout=900)
end = time.time()
print 'Removing redundant sequences took '+str(int(end - start))+' seconds'

#add query sequence to list
querry_sequences = [] # Setup an empty list
querry_sequences.append(SeqIO.read('./uploads/'+FILENAME, "fasta"))
for record in SeqIO.parse('./processing/cdhit.tmp', "fasta"):
    querry_sequences.append(record)
SeqIO.write(querry_sequences, './processing/cdhit.tmp', "fasta")

RUNCLUSTAL = CLUSTAL+' --iter='+ALIGNMENTITERATIONS+' -i ./processing/cdhit.tmp  -o ./processing/clustal.tmp --outfmt=fa --force -v -v'  
print RUNCLUSTAL
print "Aligning sequences using Clustal Omega. This can take a few minutes"
start = time.time()
command = Command(RUNCLUSTAL)
#CLUSTALOUT=command.run(timeout=900)
command.run(timeout=900)
end = time.time()
print 'Aligning sequences took '+str(int(end - start))+' seconds'

trimmer('./processing/clustal.tmp', FILENAME, MAXIMUMREDUNDANCYTHRESHOLD)


#trimmer('test')
#cleanexit()
# when complete sequences =1, curl BLASTOUT[1] VERSIONS (first entry on line), ideally replace existing sequences with curled (complete sequences)
