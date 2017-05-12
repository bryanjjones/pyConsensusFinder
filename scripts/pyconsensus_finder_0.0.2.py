#!/usr/bin/python
import sys
import _mypath
import numpy as np # need numpy for arrays and the like
import ConfigParser# need to read Config file
import time # need to time various steps
from runbin import Command # need to run binaries (e.g. BLAST, CDHIT, CLUSTAL O) with timeout function
import analyze # for processing alingments and suggesting mutations
import os
from Bio import SeqIO, Entrez, SeqRecord, Seq

programstart = time.time()
warnings = []

#DEFAULT SETTINGS, only used if not specified in config.cfg, else they get replaced
FILENAME = None
Entrez.email = None
MAXIMUMSEQUENCES = 2000
BLASTEVALUE = 1e-3
CONSENSUSTHRESHOLD = 0
RATIO = 7
USECOMPLETESEQUENCES = 1
ALIGNMENTITERATIONS = 1
MAXIMUMREDUNDANCYTHRESHOLD = 0.9
LOGGING = 0
KEEPTEMPFILES = 0

#function to exit cleanly by deleting any temporary files (unless indicated to save them)
def cleanexit(message):
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

#read from the config file
Config = ConfigParser.SafeConfigParser()

#Config.read("./config/default.cfg") #read defaults from here
Config.read("./config/config.cfg") #overwrite defaults from here

def get_with_default(Config,section,name,default):
    if Config.has_option(section,name):
        return Config.get(section,name)
    else:
        return default

# check directory structure
dirs = filter(os.path.isdir, ['./config', './uploads', './processing', './completed', './modules']) #check for presence of directory structure
if len(dirs) < 5:
     cleanexit('Missing directories. Expected to find /config, /uploads, /processing, /completed, and /modules in the working directory.')
bins = filter(os.path.isfile, [BLAST, CDHIT, CLUSTAL])
if len(bins) < 3:
     cleanexit('Missing binaries. Expected to find '+BLAST+', '+CDHIT+', and '+CLUSTAL+'.')

#check for presence of config file
if not os.path.isfile('./config/config.cfg'):
    cleanexit('Config file missing. Expected to find it at ./config/config.cfg')

#get filename for query file and check if the file is valid.
FILENAME=get_with_default(Config, 'BasicSettings', 'filename', FILENAME)
#check if there is a specified file
if not FILENAME:
    print "No query file specified in CONFIGFILE. EXITING"
    sys.exit("No query file specified in CONFIGFILE. EXITING")
#check if that specified file actually exists
if not os.path.isfile('./uploads/'+FILENAME):
    if os.path.isfile('./'+FILENAME):
        os.rename('./'+FILENAME, './uploads/'+FILENAME)
    else:
        cleanexit('Query file specified in the config file (config.cfg), '+FILENAME+' not found. Expected to find it at ./uploads/'+FILENAME)
#check if that specified file has only one sequence
if not 1 == (len(list(SeqIO.parse('./uploads/'+FILENAME, "fasta")))):
    message = "WARNING, MORE THAN ONE SEQUENCE FOUND IN "+FILENAME+". Using only the first sequence. "
    print message
    warnings.append(message)
#check to verify the file is a protein sequence
protonly = ('F','L','I','M','V','S','P','Y','H','Q','K','D','E','W','R')
hasprotonly=0
for letter in protonly:
    if (next(SeqIO.parse('./uploads/'+FILENAME, "fasta"))).seq.count(letter):
        hasprotonly=1
if not hasprotonly:
    message = "WARNING, IT LOOKS LIKE YOUR SEQUENCE IS NOT PROTEIN. Consensus Finder works on protein sequences. "
    print message
    warnings.append(message)
    print "Attempting to procede anyway."

#Get other variables from config.cfg file
Entrez.email = get_with_default(Config, 'BasicSettings', 'Email', Entrez.email)
MAXIMUMSEQUENCES = int(get_with_default(Config, 'BlastSettings', 'MAXIMUMSEQUENCES', MAXIMUMSEQUENCES))
BLASTEVALUE = float(get_with_default(Config, 'BlastSettings', 'BLASTEVALUE', BLASTEVALUE))
CONSENSUSTHRESHOLD = float(get_with_default(Config, 'AlignmentSettings', 'ConsensusThreshold', CONSENSUSTHRESHOLD))
RATIO = float(get_with_default(Config, 'AlignmentSettings', 'ConsensusRatio', RATIO))
USECOMPLETESEQUENCES = bool('yes' == get_with_default(Config,'AlignmentSettings', 'UseCompleteSequences', USECOMPLETESEQUENCES))
ALIGNMENTITERATIONS = int(get_with_default(Config,'AlignmentSettings', 'AlignmentIterations', ALIGNMENTITERATIONS))
MAXIMUMREDUNDANCYTHRESHOLD = float(get_with_default(Config, 'AlignmentSettings', 'MaximumRedundancyThreshold', MAXIMUMREDUNDANCYTHRESHOLD))
LOGGING = bool('yes' == get_with_default(Config,'TroubleShooting', 'Logging', LOGGING))
KEEPTEMPFILES = bool('yes' == get_with_default(Config,'TroubleShooting', 'KeepTempFiles',KEEPTEMPFILES))

#Need to check if user specified ratio and threshold.
#If both are specified results from both will be suggested. Otherwise, only the one specified will be used.
#So if user specified threshold, but no ratio, don't use default ratio.
#Default ratio should only be used when user specified neither ratio nor thershold in config file.
hasratio = Config.has_option('AlignmentSettings', 'ConsensusRatio')
hasthreshold = Config.has_option('AlignmentSettings', 'ConsensusThreshold')
if hasthreshold and not hasratio:
    RATIO=None

#Run Blast
RUNBLAST = BLAST+' -db nr -query ./uploads/'+FILENAME+' -evalue '+str(BLASTEVALUE)+' -max_target_seqs '+str(MAXIMUMSEQUENCES)+' -outfmt "6 sacc sseq pident" -remote' 
print "Begining BLAST search of NCBI. This will take a few minutes."
start = time.time()

command = Command(RUNBLAST)
BLASTOUT=command.run(timeout=1800)
if command.status == -15: #error code from Command indicating timeout
    message = 'BLAST TOOK TOO LONG. Maximum sequences reduced to 200 BLAST results. '
    print message
    RUNBLAST = BLAST+' -db nr -query ./uploads/'+FILENAME+' -evalue '+str(BLASTEVALUE)+' -max_target_seqs 200 -outfmt "6 saccver sseq pident" -remote' #repeat blast search with only 200 max sequences. Thus even if user (or default) specifies too many sequences, they can still get some useful output.
    print "Begining BLAST search of NCBI. This will take a few minutes."
    start = time.time()
    warnings.append(message) #add warning that blast took too long to the warnings list so the user knows what happened.
    command = Command(RUNBLAST)
    BLASTOUT=command.run(timeout=1800)
    if command.status == -15: #error code from Command indicating timeout
        cleanexit('BLAST STILL TOOK TOO LONG! Giving up.')
if command.status != 0: #any other error code
    cleanexit('BLAST FAILED. I do not know why, check your inputs, internet connection, etc.')

#temporary data set to skip actual blast search just to speed troubleshooting
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
BLASTOUT=BLASTOUT[1].splitlines() #redefines BLASTOUT as the second item in the list, since BLASTOUT is the stdout from the blast search, all the data is in position 1
for index in range(len(BLASTOUT[:])): #for each sequence split at the white space (between VERSION and sequence)
    BLASTOUT[index]=BLASTOUT[index].split()
VERSIONS = [item[0] for item in BLASTOUT]
print 'BLAST returned '+str(len(VERSIONS))+' sequences.'
if not len(VERSIONS):
    cleanexit('BLAST did not return any hits. Check your input sequence and blast parameters. \nTry increasing BlastEValue. If it still does not work, you may have an orphan sequence with no known homologs and Consensus Finder is unable to help you.')

record=[] #initialize list for recording sequences
if USECOMPLETESEQUENCES: #if use complete sequences is true, dowload sequences from Entrez, otherwise just use returned BLAST sequences
    start = time.time() #start timer for downloads
    print 'Downloading complete sequences from NCBI'
    handle = Entrez.efetch(db="protein", id=",".join(VERSIONS), retmax=MAXIMUMSEQUENCES, rettype="fasta", retmode="text") #retrieving all sequences from Entrez, db="sequences"
    record = list(SeqIO.parse(handle, "fasta"))
    end = time.time()
    print 'Downloading sequences took '+str(int(end - start))+' seconds'
else:
    for item in BLASTOUT:
        item[1]=item[1].translate(None, '-') #remove gaps (-) since they would mess up CD-HIT later
        record.append(SeqRecord.SeqRecord(Seq.Seq(item[1]),id=item[0],description="")) # add each sequence to the record
SeqIO.write(record, './processing/'+FILENAME+'_all_sequences.tmp', "fasta") # Save all sequences identified by blast as fasta file for processing with CD-HIT

#process sequences with CDHIT to cluster similar sequences and take only one sequence from each cluster. This helps eliminated oversampling bias.
RUNCDHIT = CDHIT+' -i ./processing/'+FILENAME+'_all_sequences.tmp -o ./processing/cdhit.tmp -c '+str(MAXIMUMREDUNDANCYTHRESHOLD)+' -M 0 -T 0'  
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

#add query sequence to the top of the list of sequences so it can be aligned too and be used for trimming.
querry_sequences = [] # Setup an empty list
querry_sequences.append(SeqIO.read('./uploads/'+FILENAME, "fasta"))
for record in SeqIO.parse('./processing/cdhit.tmp', "fasta"):
    querry_sequences.append(record)
SeqIO.write(querry_sequences, './processing/cdhit.tmp', "fasta")

#Run clustal omega alingment
RUNCLUSTAL = CLUSTAL+' --iter='+str(ALIGNMENTITERATIONS)+' -i ./processing/cdhit.tmp  -o ./processing/clustal.tmp --outfmt=fa --force -v -v'  
print "Aligning sequences using Clustal Omega. This can take a few minutes, espicially with many sequences"
start = time.time()
command = Command(RUNCLUSTAL)
CLUSTALOUT=command.run(timeout=7200)
print CLUSTALOUT[1] #output from Clustal
if command.status == -15:
    cleanexit('CLUSTAL TOOK TOO LONG! Try with fewer sequences.')
elif command.status != 0:
    cleanexit('CLUSTAL FAILED. I do not know why, check your inputs')
end = time.time()
print 'Aligning sequences took '+str(int(end - start))+' seconds'

#processing steps from analyze module.
#Trimmer will trim the alingment down to the size of the query sequence.
#aacounts will count the ocourances of each amino acid at each position and save as csv
#aafrequencies will calculate amino acid frequencies from the counts at each position and save as csv.
#consensus will calculate the overall consensus sequence from the amino acid frequencies, and save it as fasta
trimmed = analyze.trimmer('./processing/clustal.tmp',filename='./completed/'+FILENAME+'_trimmed_alignment.fst')
counts = analyze.aacounts(trimmed, filename='./completed/'+FILENAME+'_counts.csv')
freqs = analyze.aafrequencies(counts, filename='./completed/'+FILENAME+'_frequencies.csv')
consensus = analyze.consensus(freqs, filename='./completed/'+FILENAME+'_consensus.fst')

#get suggested mutations based either upon specified ratio, cutoff, or both (depending on what was specified)
if RATIO:
    ratiomutations = analyze.ratioconsensus(trimmed[0], freqs, RATIO)
if CONSENSUSTHRESHOLD:
    thresholdmutations = analyze.cutoffconsensus(trimmed[0], freqs, CONSENSUSTHRESHOLD)
    mutations = thresholdmutations
    if RATIO:
        mutations = np.append(mutations, ratiomutations)
else:
	mutations = ratiomutations
#format mutations list by eliminating duplicates and adding descriptor text
mutations, output = analyze.formatmutations(mutations)

#Save output suggestions file with any warnings that have been added
file = open('./completed/'+FILENAME+'_mutations.txt','wb')
file.write(''.join(warnings) + '\n')# 's16\n')
np.savetxt(file, output, delimiter=",", fmt='%s') 

#move query file from uploads directory to completed to keep uploads clean
os.rename('./uploads/'+FILENAME, './completed/'+FILENAME)
programend = time.time()
print 'Consensus Finder Completed.'
print 'Your results are in the ./completed/ directory.'
print ''.join(warnings)
print 'Process took '+str(int(programend - programstart))+' seconds'
cleanexit(0)
