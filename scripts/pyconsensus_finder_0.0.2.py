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
from Bio import SeqIO, Entrez, SeqRecord, Seq

def cleanexit(message): #function to exit cleanly by deleting any temporary files (unless indicated to save them)
    if KEEPTEMPFILES:
        print 'Leaving temporary files'
    else:
        print 'Deleting temporary files'
        for file in ['./processing/'+FILENAME+'_all_sequences.tmp','./processing/cdhit.tmp','./processing/cdhit.tmp.clstr','./processing/clustal.tmp']:
            if os.path.isfile(file):
                os.remove(file)
    print 'Exiting'
    sys.exit(message)

#location of binaries to call
BLAST="./binaries/blastp" #change this if another version is installed locally
CDHIT="./binaries/cd-hit" #change this if another version is installed locally
CLUSTAL="./binaries/clustalo-1.2.0-Ubuntu-x86_64" #change this if another version is installed locally

programstart = time.time()
# check directory structure
dirs = filter(os.path.isdir, ['./config', './uploads', './processing', './completed', './modules']) #check for presence of directory structure
if len(dirs) < 5:
     cleanexit('Missing directories. Expected to find /config, /uploads, /processing, /completed, and /modules in the working directory.')
bins = filter(os.path.isfile, [BLAST, CDHIT, CLUSTAL])
if len(bins) < 3:
     cleanexit('Missing binaries. Expected to find '+BLAST+', '+CDHIT+', and '+CLUSTAL+'.')
#check for presence of config files
if not os.path.isfile('./config/default.cfg'):
    cleanexit('Defaults config file missing. Expected to find it at /config/default.cfg')
if not os.path.isfile('./config/config.cfg'):
    cleanexit('Config file missing. Expected to find it at /config/config.cfg')
#read from the config file
Config = ConfigParser.ConfigParser()
Config.read("./config/default.cfg") #read defaults from here
Config.read("./config/config.cfg") #overwrite defaults from here

FILENAME=Config.get('BasicSettings', 'filename')
if not FILENAME:
    print "No query file specified in CONFIGFILE. EXITING"
    sys.exit("No query file specified in CONFIGFILE. EXITING")
if not 1 == (len(list(SeqIO.parse('./uploads/'+FILENAME, "fasta")))):
    print "WARNING, MORE THAN ONE SEQUENCE FOUND IN "+FILENAME+". Proceding using only the first sequence."
#check to verify the file is a protein sequence
protonly = ('F','L','I','M','V','S','P','Y','H','Q','K','D','E','W','R')
hasprotonly=0
for letter in protonly:
    if (next(SeqIO.parse('./uploads/'+FILENAME, "fasta"))).seq.count(letter):
        hasprotonly=1
if not hasprotonly:
    print "WARNING, IT LOOKS LIKE YOUR SEQUENCE IS NOT PROTEIN. Consensus Finder works on protein sequences. Proceding anyway."
#Get other variables from config.cfg file
Entrez.email = Config.get('BasicSettings', 'Email')
MAXIMUMSEQUENCES=Config.get('BlastSettings', 'MAXIMUMSEQUENCES')
BLASTEVALUE=Config.get('BlastSettings', 'BLASTEVALUE')
CONSENSUSTHRESHOLD=Config.get('AlignmentSettings', 'ConsensusThreshold')
USECOMPLETESEQUENCES=Config.getboolean('AlignmentSettings', 'UseCompleteSequences')
ALIGNMENTITERATIONS=Config.get('AlignmentSettings', 'AlignmentIterations')
MAXIMUMREDUNDANCYTHRESHOLD=Config.get('AlignmentSettings', 'MaximumRedundancyThreshold')
LOGGING=Config.getboolean('TroubleShooting', 'Logging')
KEEPTEMPFILES=Config.getboolean('TroubleShooting', 'KeepTempFiles')

#Run Blast
RUNBLAST = BLAST+' -db nr -query ./uploads/'+FILENAME+' -evalue '+BLASTEVALUE+' -max_target_seqs '+MAXIMUMSEQUENCES+' -outfmt "6 sacc sseq pident" -remote' 
print "Begining BLAST search of NCBI. This will take a few minutes."
start = time.time()
REDUCED_MAX_SEQS=0 #ticker to indicate if maximum sequences had to be reduced due to blast timeout
'''
command = Command(RUNBLAST)
BLASTOUT=command.run(timeout=2000)
if command.status == -15: #error code from Command indicating timeout
    print'BLAST TOOK TOO LONG! Tring with 200 maximum sequences.'
    RUNBLAST = BLAST+' -db nr -query ./uploads/'+FILENAME+' -evalue '+BLASTEVALUE+' -max_target_seqs 200 -outfmt "6 saccver sseq pident" -remote' #repeat blast search with only 200 max sequences
    print "Begining BLAST search of NCBI. This will take a few minutes."
    start = time.time()
    REDUCED_MAX_SEQS=1 #ticker to indicate that we had to reduce the maximum sequences to get results
    command = Command(RUNBLAST)
    BLASTOUT=command.run(timeout=2000)
    if command.status == -15: #error code from Command indicating timeout
        cleanexit('BLAST STILL TOOK TOO LONG! Giving up.')
if command.status != 0: #any other error code
    cleanexit('BLAST FAILED. I do not know why, check your inputs, internet connection, etc.')
'''
#temporary data set to skip actual blast search

BLASTOUT=[0,"""WP_083418319    PCPSFLIEHDRGLVLFDAGFDPRGLDDMAAYYPEISKALPMAGNRDLGIDRQLDGLGYRPSQISYVIPSHLHFDHAGGLYLFPDSTFLMGSGEMAFALQAHDKPQ---AGFFRVEDLLPTRHFDW--IETAHDFDLFGDGSVVLLFSPGHTPGSLALFVRLPNQ-NIILSGDVCHFPLEVDMGIIATSSFSPSYATF-ALRRLRMISKAWDARIWIQHEEDHWNEWPHAPE 26.840
WP_073456977    PMPAYLIEHPKGLVLFDTMLVPDAADDPERVYGPLAEHLGLKYTREQRVDNQIKALGYRLEDVTHVIASHTHFDHSGGLYLFPHAKFYVGEGELRFAL--WPDPAGAGFFRQADIEA--TRSFNW--VQVGFDHDLFGDGSVVVLHTPGHTPGELSLLVRL-KSRNFILTGDTVHLRQALEDEIPMP---YDSNTELAIRSIRRLKLLRESADATVWITHDPEDWAEFQHAPYCY       29.362
WP_049268273    ESYEIPVPWFLLTHPDGFTLIDGGLAVEGLKDPFAYWGGAVEQFKPVMPEEQGCL-EQLKRIGVAPEDIRYVILSHLHSDHTGAIGRFPHATHVVQRREYEYAFA--PDWFTSGAYCRYDYD---HPELNWFFLNGLSEDNYDLYGDGTLQCIFTPGHSPGHQSFLIRLSSGTNFTLAIDAAYTLDHYHERAL-PGLMTSATDVAQSVQKLRQLTERYNAILIPGHDPEEWEKIRLAPAWY 29.876
WP_025778256    LAVPIPTFLIQHEGGLLVFDTGLATDAAGDPARAYGPLAEAFDMSFPPEARIDTQLESLGFSTSDVTDVVLSHMHFDHTGGLELFPTARGFIGEGEL-----AYSRSPRRLDAAMYREEDIAAAGQIDWLEIPQGVDHDIFGDGSVVVLSMPGHTHGTLSLKLSPPDHRTIILTSDAAHLQSNIDETTGMPLD-----VDTRNKERSLRRLRLLASQPNTTVWANHDPDHWKQFRR      27.966
OCL02351        KLWLLHVGNLECDEAWFKRGGGTSTLSNPHPNRERRKLIMVSVLIEHPVEGLILFETGSGKDY--PE-IWGAPINDIFARVDYTEEQELDVQIKKTGHDIKDVKMVVIGHLHLDHAGGLEYFRNTGIPIYVHEKELKHAFYSVATKSDLGVY----LPHYLTFDLNW--VPFSGAFLEIAQGLN-LHHAPGHTPGLCILQVNMPKSGAWIFTTDMYHVSENFEESVPQGWLAREHDDWVRSNHMIHMLQRRTGARMVFGHCTKALEGLNMAPHAY       26.545"""]

end = time.time()
print 'BLAST search took '+str(int(end - start))+' seconds'

#print 'BLAST Results'
BLASTOUT=BLASTOUT[1].splitlines() #redefines BLASTOUT as the second item in the list, since BLASTOUT is the stdout from the blast search, all the data is in position 1
for index in range(len(BLASTOUT[:])): #for each sequence split at the white space (between VERSION and sequence)
    BLASTOUT[index]=BLASTOUT[index].split()
VERSIONS = [item[0] for item in BLASTOUT]
print 'BLAST returned '+str(len(VERSIONS))+' sequences.'

record=[] #initialize list for recording sequences
if USECOMPLETESEQUENCES: #if use complete sequences is true, dowload sequences from Entrez, otherwise just use returned BLAST sequences
    start = time.time() #start timer for downloads
    print 'Downloading complete sequences from NCBI'
    handle = Entrez.efetch(db="sequences", id=",".join(VERSIONS), retmax=MAXIMUMSEQUENCES, rettype="fasta", retmode="text") #retrieving all sequences from Entrez
    record = list(SeqIO.parse(handle, "fasta"))
    end = time.time()
    print 'Downloading sequences took '+str(int(end - start))+' seconds'
else:
    for item in BLASTOUT:
        item[1]=item[1].translate(None, '-') #remove gaps (-) since they would mess up CD-HIT later
        record.append(SeqRecord.SeqRecord(Seq.Seq(item[1]),id=item[0],description="")) # add each sequence to the record

SeqIO.write(record, './processing/'+FILENAME+'_all_sequences.tmp', "fasta") # Save all sequences identified by blast as fasta file for processing with CD-HIT

RUNCDHIT = CDHIT+' -i ./processing/'+FILENAME+'_all_sequences.tmp -o ./processing/cdhit.tmp -c '+MAXIMUMREDUNDANCYTHRESHOLD+' -M 0 -T 0'  
print "Removing redundant sequences using CD-HIT"
start = time.time()
command = Command(RUNCDHIT)
CDHITOUT=command.run(timeout=900)
if command.status == -15: #error code for Command timeout
    cleanexit('CDHIT TOOK TOO LONG! Try with fewer sequences.')
elif command.status != 0: #any other error code
    cleanexit('CDHIT FAILED. I do not know why, check your inputs.')
end = time.time()
print 'Removing redundant sequences took '+str(int(end - start))+' seconds'

#add query sequence to list
querry_sequences = [] # Setup an empty list
querry_sequences.append(SeqIO.read('./uploads/'+FILENAME, "fasta"))
for record in SeqIO.parse('./processing/cdhit.tmp', "fasta"):
    querry_sequences.append(record)
SeqIO.write(querry_sequences, './processing/cdhit.tmp', "fasta")

#Run clustal omega alingment
RUNCLUSTAL = CLUSTAL+' --iter='+ALIGNMENTITERATIONS+' -i ./processing/cdhit.tmp  -o ./processing/clustal.tmp --outfmt=fa --force -v -v'  
print "Aligning sequences using Clustal Omega. This can take a few minutes, espicially with many sequences"
start = time.time()
command = Command(RUNCLUSTAL)
CLUSTALOUT=command.run(timeout=7200)
#command.run(timeout=900)
print CLUSTALOUT
if command.status == -15:
    cleanexit('CLUSTAL TOOK TOO LONG! Try with fewer sequences.')
elif command.status != 0:
    cleanexit('CLUSTAL FAILED. I do not know why, check your inputs')
end = time.time()
print 'Aligning sequences took '+str(int(end - start))+' seconds'

trimmer('./processing/clustal.tmp', FILENAME, CONSENSUSTHRESHOLD) #Trim alingment to query length, and write output files.
programend = time.time()
print 'Consensus Finder Completed.'
print 'Your results are in the /completed/ directory.'
print 'Process took '+str(int(programend - programstart))+' seconds'
cleanexit(0)