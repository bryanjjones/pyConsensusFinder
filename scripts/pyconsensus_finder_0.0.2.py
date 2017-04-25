#!/usr/bin/python
import sys
import _mypath
import numpy as np # need numpy for arrays and the like
#import subprocess # need to run binaries (blast, CD-hit, clustal)
import ConfigParser# need to read Config file
import time
from runbin import Command
import os

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

#read from the config file
Config = ConfigParser.ConfigParser()
Config.read("./config/default.cfg") #read defaults from here
Config.read("./config/config.cfg") #overwrite defaults from here

#need to add a check to verify that the config.cfg file is there, otherwise results in a KeyError

FILENAME=Config.get('BasicSettings', 'filename')
if not FILENAME:
#	print FILENAME
	print "No query file specified in CONFIGFILE. EXITING"
	sys.exit(0)

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
print RUNBLAST

start = time.time()
#command = Command(RUNBLAST)
#BLASTOUT=command.run(timeout=900)

#add instructions for timeout (i.e. try again w/ fewer sequences)

BLASTOUT=[0,"""WP_083418319    PCPSFLIEHDRGLVLFDAGFDPRGLDDMAAYYPEISKALPMAGNRDLGIDRQLDGLGYRPSQISYVIPSHLHFDHAGGLYLFPDSTFLMGSGEMAFALQAHDKPQ---AGFFRVEDLLPTRHFDW--IETAHDFDLFGDGSVVLLFSPGHTPGSLALFVRLPNQ-NIILSGDVCHFPLEVDMGIIATSSFSPSYATF-ALRRLRMISKAWDARIWIQHEEDHWNEWPHAPE 26.840
WP_073456977    PMPAYLIEHPKGLVLFDTMLVPDAADDPERVYGPLAEHLGLKYTREQRVDNQIKALGYRLEDVTHVIASHTHFDHSGGLYLFPHAKFYVGEGELRFAL--WPDPAGAGFFRQADIEA--TRSFNW--VQVGFDHDLFGDGSVVVLHTPGHTPGELSLLVRL-KSRNFILTGDTVHLRQALEDEIPMP---YDSNTELAIRSIRRLKLLRESADATVWITHDPEDWAEFQHAPYCY       29.362
WP_049268273    ESYEIPVPWFLLTHPDGFTLIDGGLAVEGLKDPFAYWGGAVEQFKPVMPEEQGCL-EQLKRIGVAPEDIRYVILSHLHSDHTGAIGRFPHATHVVQRREYEYAFA--PDWFTSGAYCRYDYD---HPELNWFFLNGLSEDNYDLYGDGTLQCIFTPGHSPGHQSFLIRLSSGTNFTLAIDAAYTLDHYHERAL-PGLMTSATDVAQSVQKLRQLTERYNAILIPGHDPEEWEKIRLAPAWY 29.876
WP_025778256    LAVPIPTFLIQHEGGLLVFDTGLATDAAGDPARAYGPLAEAFDMSFPPEARIDTQLESLGFSTSDVTDVVLSHMHFDHTGGLELFPTARGFIGEGEL-----AYSRSPRRLDAAMYREEDIAAAGQIDWLEIPQGVDHDIFGDGSVVVLSMPGHTHGTLSLKLSPPDHRTIILTSDAAHLQSNIDETTGMPLD-----VDTRNKERSLRRLRLLASQPNTTVWANHDPDHWKQFRR      27.966
OCL02351        KLWLLHVGNLECDEAWFKRGGGTSTLSNPHPNRERRKLIMVSVLIEHPVEGLILFETGSGKDY--PE-IWGAPINDIFARVDYTEEQELDVQIKKTGHDIKDVKMVVIGHLHLDHAGGLEYFRNTGIPIYVHEKELKHAFYSVATKSDLGVY----LPHYLTFDLNW--VPFSGAFLEIAQGLN-LHHAPGHTPGLCILQVNMPKSGAWIFTTDMYHVSENFEESVPQGWLAREHDDWVRSNHMIHMLQRRTGARMVFGHCTKALEGLNMAPHAY       26.545"""]
end = time.time()
print 'Time elapsed (seconds)'
print(end - start)
#print 'BLAST Results'
BLASTOUT=BLASTOUT[1].splitlines()
for index in range(len(BLASTOUT[:])):
	BLASTOUT[index]=BLASTOUT[index].split()
VERSIONS = [item[0] for item in BLASTOUT]

#if USECOMPLETESEQUENCES:
	
#	curl -s  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&amp;id=${VERSIONS[i]}&amp;retmode=text&amp;rettype=fasta" | sed -r 's:'\>.+'::g' | sed ':a;N;$!ba;s/\n//g') 

#SEQS=SEQS.split()
#SEQS = []
#np.asarray(SEQS)
#print SEQS[0,1]
#print SEQS[1,1]
#	np.asarray(a)

# when complete sequences =1, curl BLASTOUT[1] VERSIONS (first entry on line), ideally replace existing sequences with curled (complete sequences)