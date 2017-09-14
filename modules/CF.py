#!/usr/bin/python
import os
import Bio.SeqIO
import Bio.Entrez
import Bio.SeqRecord
import Bio.Seq
import sys
import runbin
import time
import StringIO
import _mypath
import analyze
import ConfigParser

HOME = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
#Main CF program, takes settings object and runs blast, cdhit, and clustalo. Then trims gaps, calculates counts, frequencies, and consensus at each position. Finally produces an output file with suggested mutations.
class CF(object):
    def __init__(self, defaults, configfile, write=True):
        settings=setsettings(defaults,configfile)
        #do checks of settings
        warnings=[]
        #Run CF checks to check all settings, return any warnings to warnings variable
        warnings=warnings+checks(settings).warnings

        if write:
            trimmed_output=HOME+'/completed/'+settings.FILENAME+'_trimmed_alignment.fst'
            counts_output=HOME+'/completed/'+settings.FILENAME+'_counts.csv'
            freqs_output=HOME+'/completed/'+settings.FILENAME+'_frequencies.csv'
            consensus_output=HOME+'/completed/'+settings.FILENAME+'_consensus.fst'
            summary_output=HOME+'/completed/'+settings.FILENAME+'_mutations.txt'
        else:
            trimmed_output=None
            counts_output=None
            freqs_output=None
            consensus_output=None
            summary_output=None
        #run binaries, pass results from blast to cdhit, and from cdhit to clustalo
        self.blast=runblast(settings)
        self.cdhit=runcdhit(settings,self.blast.out)
        self.clustalo=runclustalo(settings,self.cdhit.out)
        #processing steps from analyze module.
        self.aaaray = analyze.aaaray(self.clustalo.out)
        self.trimmed = analyze.trimmer(self.aaaray, sequenceids=self.clustalo.out, filename=trimmed_output)
        self.counts = analyze.aacounts(self.trimmed, filename=counts_output)
        self.freqs = analyze.aafrequencies(self.counts, filename=freqs_output)
        self.consensus = analyze.consensus(self.freqs, filename=consensus_output)
        self.mutations, self.output = analyze.mutations(settings, self.trimmed, self.freqs)
        self.warnings=warnings + self.blast.warnings
        self.settings = settings
        if write:
            analyze.saveoutput(settings, self.warnings, self.output, summary_output)

#clean up any files in /processing and exit the program displaying any error messages.
def cleanexit(message=None, keeptemp=False): #function to exit cleanly by deleting any temporary files (unless indicated to save them)
    if keeptemp:
        print('Leaving temporary files')
    else:
        print('Deleting temporary files...')
        for file in os.listdir(HOME+"/processing"):
            os.remove(HOME+"/processing/"+file)
            print(file+" deleted")
    print('Exiting')
    sys.exit(message)

#assign settings based upon specified defaults and values read from configfile
class setsettings(object):
    def __init__(self, defs, configfile=HOME+"/config/config.cfg"):
        #check for presence of config file
        if not os.path.isfile(configfile):
            cleanexit('Config file missing. Expected to find it at '+configfile)
        #check threshold only
        NoDefaultConfig=ConfigParser.SafeConfigParser()
        NoDefaultConfig.read(configfile)
        self.hasthreshold = NoDefaultConfig.has_option('Options', 'Consensus_Threshold')
        #default ratio should only be used when user specified neither ratio nor thershold in config file
        if self.hasthreshold:
            defs['Consensus_Ratio'] ='0'
        self.hasratio = NoDefaultConfig.has_option('Options', 'Consensus_Ratio')
        #set defaults and read the config file
        Config = ConfigParser.SafeConfigParser(defaults=defs)
        Config.read(configfile) 
        #get variables from config file
        self.FILENAME = Config.get('Options', 'File_Name')
        self.EMAIL = Config.get('Options', 'Email')
        self.MAXIMUMSEQUENCES = Config.getint('Options', 'Maximum_Sequences')
        self.BLASTEVALUE = Config.getfloat('Options', 'Blast_E_Value')
        self.CONSENSUSTHRESHOLD = Config.getfloat('Options', 'Consensus_Threshold')
        self.RATIO = Config.getfloat('Options', 'Consensus_Ratio')
        self.USECOMPLETESEQUENCES = Config.getboolean('Options', 'Use_Complete_Sequences')
        self.ALIGNMENTITERATIONS = Config.getint('Options', 'Alignment_Iterations')
        self.MAXIMUMREDUNDANCYTHRESHOLD = Config.getfloat('Options', 'Maximum_Redundancy_Threshold')
        #logging and keeptemfiles don't do anything as of now.
        self.LOGGING = Config.getboolean('Options', 'Logging')
        self.KEEPTEMPFILES = Config.getboolean('Options', 'Keep_Temp_Files')
        self.PDB = Config.get('Options', 'PDB_Name')
        self.CHAIN = Config.get('Options', 'Chain')
        self.RESIDUE = Config.getint('Options', 'Residue')
        self.ANG = Config.getfloat('Options', 'Angstrom')
        #location of binaries to call
        self.BLAST = Config.get('Options', 'Blast_binary')
        self.CDHIT = Config.get('Options', 'CDHIT_binary')
        self.CLUSTAL = Config.get('Options', 'ClustalO_binary')
        
# check file structure, query is given, sequence is protein, if binaries are present.
# exits if fatal problems are present, returns 'warnings' if non-fatal problems are present
class checks(object):
    def __init__(self, settings):
        self.warnings=[]
        #check directory structure
        dirs = filter(os.path.isdir, [HOME+'/config', HOME+'/uploads', HOME+'/processing', HOME+'/completed', HOME+'/modules']) #check for presence of directory structure
        if len(dirs) < 5:
             cleanexit('Missing directories. Expected to find /config, /uploads, /processing, /completed, and /modules in the working directory.')
        #check query file
        if not settings.FILENAME:
            cleanexit("No query file specified in CONFIGFILE. EXITING")
        if not os.path.isfile(HOME+'/uploads/'+settings.FILENAME):
            cleanexit('No query file. File '+settings.FILENAME+' does not exist in uploads directory.')
        if not 1 == (len(list(Bio.SeqIO.parse(HOME+'/uploads/'+settings.FILENAME, "fasta")))):
            message = "WARNING, MORE THAN ONE SEQUENCE FOUND IN "+settings.FILENAME+". Using only the first sequence. "
            print message
            self.warnings.append(message)
        #check to verify the file is a protein sequence
        protonly = ('F','L','I','M','V','S','P','Y','H','Q','K','D','E','W','R','f','l','i','m','v','s','p','y','h','q','k','d','e','w','r')
        hasprotonly=0
        for letter in protonly:
            if (next(Bio.SeqIO.parse(HOME+'/uploads/'+settings.FILENAME, "fasta"))).seq.count(letter):
                hasprotonly=1
        if not hasprotonly:
            message = "WARNING, IT LOOKS LIKE YOUR SEQUENCE IS NOT PROTEIN. Consensus Finder works on protein sequences. "
            print message
            self.warnings.append(message)
            print "Attempting to procede anyway."
        for binary in [settings.BLAST, settings.CLUSTAL, settings.CDHIT]:
            if not (os.path.isfile(binary)):
                cleanexit('Missing binaries. Expected to find '+binary+'.')

#Run Blast
#returns 'warnings'= warning if blast settings were changed, 'blastout'=blast results as Bio SeqRecord objects, 'versions'=list of version numbers, 'out'=list of complete sequences as Bio SeqRecord objects, and (if 'filename' is given) writes output to fasta formatted file
class runblast(object):
    def __init__(self, settings, filename=None):
        RUNBLAST = settings.BLAST+' -db nr -query '+HOME+'/uploads/'+settings.FILENAME+' -evalue '+str(settings.BLASTEVALUE)+' -max_target_seqs '+str(settings.MAXIMUMSEQUENCES)+' -outfmt "6 sacc sseq pident" -remote' 
        print "\nBegining BLAST search of NCBI. This will take a few minutes."
        start = time.time()
        self.warnings=[]
        command = runbin.Command(RUNBLAST)
        self.blastout=command.run(timeout=1800) #died by itself at 1884, succeeded at  1225 seconds
        if command.status == -15: #error code from Command indicating timeout
            message = 'BLAST TOOK TOO LONG. Maximum sequences reduced to 200 BLAST results. '
            print message
            RUNBLAST = settings.BLAST+' -db nr -query '+HOME+'/uploads/'+settings.FILENAME+' -evalue '+str(settings.BLASTEVALUE)+' -max_target_seqs 200 -outfmt "6 saccver sseq pident" -remote' #repeat blast search with only 200 max sequences
            print "Begining BLAST search of NCBI. This will take a few minutes."
            start = time.time()
            self.warnings=message
            command = runbin.Command(RUNBLAST)
            self.blastout=command.run(timeout=2000)
            if command.status == -15: #error code from Command indicating timeout
                cleanexit('BLAST STILL TOOK LONG! Even after reducing maximum sequences. Giving up.',keeptemp=settings.KEEPTEMPFILES)
        if command.status != 0: #any other error code
            cleanexit('BLAST FAILED. I do not know why, check your inputs, internet connection, etc.',keeptemp=settings.KEEPTEMPFILES)

        #temporary data set to skip actual blast search for testing
        '''
        self.blastout=[0,"""WP_083418319    PCPSFLIEHDRGLVLFDAGFDPRGLDDMAAYYPEISKALPMAGNRDLGIDRQLDGLGYRPSQISYVIPSHLHFDHAGGLYLFPDSTFLMGSGEMAFALQAHDKPQ---AGFFRVEDLLPTRHFDW--IETAHDFDLFGDGSVVLLFSPGHTPGSLALFVRLPNQ-NIILSGDVCHFPLEVDMGIIATSSFSPSYATF-ALRRLRMISKAWDARIWIQHEEDHWNEWPHAPE 26.840
        WP_073456977    PMPAYLIEHPKGLVLFDTMLVPDAADDPERVYGPLAEHLGLKYTREQRVDNQIKALGYRLEDVTHVIASHTHFDHSGGLYLFPHAKFYVGEGELRFAL--WPDPAGAGFFRQADIEA--TRSFNW--VQVGFDHDLFGDGSVVVLHTPGHTPGELSLLVRL-KSRNFILTGDTVHLRQALEDEIPMP---YDSNTELAIRSIRRLKLLRESADATVWITHDPEDWAEFQHAPYCY       29.362
        WP_049268273    ESYEIPVPWFLLTHPDGFTLIDGGLAVEGLKDPFAYWGGAVEQFKPVMPEEQGCL-EQLKRIGVAPEDIRYVILSHLHSDHTGAIGRFPHATHVVQRREYEYAFA--PDWFTSGAYCRYDYD---HPELNWFFLNGLSEDNYDLYGDGTLQCIFTPGHSPGHQSFLIRLSSGTNFTLAIDAAYTLDHYHERAL-PGLMTSATDVAQSVQKLRQLTERYNAILIPGHDPEEWEKIRLAPAWY 29.876
        WP_025778256    LAVPIPTFLIQHEGGLLVFDTGLATDAAGDPARAYGPLAEAFDMSFPPEARIDTQLESLGFSTSDVTDVVLSHMHFDHTGGLELFPTARGFIGEGEL-----AYSRSPRRLDAAMYREEDIAAAGQIDWLEIPQGVDHDIFGDGSVVVLSMPGHTHGTLSLKLSPPDHRTIILTSDAAHLQSNIDETTGMPLD-----VDTRNKERSLRRLRLLASQPNTTVWANHDPDHWKQFRR      27.966
        OCL02351        KLWLLHVGNLECDEAWFKRGGGTSTLSNPHPNRERRKLIMVSVLIEHPVEGLILFETGSGKDY--PE-IWGAPINDIFARVDYTEEQELDVQIKKTGHDIKDVKMVVIGHLHLDHAGGLEYFRNTGIPIYVHEKELKHAFYSVATKSDLGVY----LPHYLTFDLNW--VPFSGAFLEIAQGLN-LHHAPGHTPGLCILQVNMPKSGAWIFTTDMYHVSENFEESVPQGWLAREHDDWVRSNHMIHMLQRRTGARMVFGHCTKALEGLNMAPHAY       26.545"""]
        '''
        end = time.time()
        print('BLAST search took '+str(int(end - start))+' seconds')
        self.blastout=self.blastout[1].splitlines() #redefines BLASTOUT as the second item in the list, since BLASTOUT is the stdout from the blast search, all the data is in position 1
        for index in range(len(self.blastout[:])): #for each sequence split at the white space (between VERSION and sequence)
            self.blastout[index]=self.blastout[index].split()
        self.versions = [item[0] for item in self.blastout]
        print('BLAST returned '+str(len(self.versions))+' sequences.')
        
        self.out=[] #initialize list for recording sequences
        if settings.USECOMPLETESEQUENCES: #if use complete sequences is true, dowload sequences from Entrez, otherwise just use returned BLAST sequences
            start = time.time() #start timer for downloads
            print('\nDownloading complete sequences from NCBI')
            Bio.Entrez.email=settings.EMAIL
            handle = Bio.Entrez.efetch(db="protein", id=",".join(self.versions), retmax=settings.MAXIMUMSEQUENCES, rettype="fasta", retmode="text") #retrieving all sequences from Entrez, db="sequences"
            self.out = list(Bio.SeqIO.parse(handle, "fasta"))
            end = time.time()
            print('Downloading sequences took '+str(int(end - start))+' seconds')
        else:
            for item in self.blastout:
                item[1]=item[1].translate(None, '-') #remove gaps (-) since they would mess up CD-HIT later
                self.out.append(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(item[1]),id=item[0],description="")) # add each sequence to the record
        if filename is not None:
            Bio.SeqIO.write(self.out, filename, "fasta") # Save all sequences identified by blast as fasta file

#Run CD-HIT
#needs settings and sequences as a list of Bio SeqRecord objects
#returns 'out'= sequences as a list of Biopython SeqRecord objects, and (if write=True) also writes fasta formatted sequences to file.
class runcdhit(object):
    def __init__(self,settings,sequences,filename=None):
        cdhitinput=HOME+'/processing/'+settings.FILENAME+'_cdhit_input.tmp'
        cdhitoutput=HOME+'/processing/'+settings.FILENAME+'_cdhit_output.tmp'
        Bio.SeqIO.write(sequences, cdhitinput , "fasta") #converts Bio SeqRecord objects from self.Bioout into fasta strings
        RUNCDHIT = settings.CDHIT+' -i '+cdhitinput+' -o '+cdhitoutput+' -c '+str(settings.MAXIMUMREDUNDANCYTHRESHOLD)+' -M 0 -T 0'  
        print("\nRemoving redundant sequences using CD-HIT")
        start = time.time()
        command = runbin.Command(RUNCDHIT)
        #self.out=command.run(timeout=900) Do I need this?
        command.run(timeout=900)
        if command.status == -15: #error code for Command timeout
            cleanexit('CDHIT TOOK TOO LONG! Try with fewer sequences.',keeptemp=settings.KEEPTEMPFILES)
        elif command.status != 0: #any other error code
            cleanexit('CDHIT FAILED. I do not know why, check your inputs.',keeptemp=settings.KEEPTEMPFILES)
        end = time.time()
        print 'Removing redundant sequences took '+str(int(end - start))+' seconds'
        #add query sequence to list
        self.out = [] # Setup an empty list
        self.out.append(next(Bio.SeqIO.parse(HOME+'/uploads/'+settings.FILENAME, "fasta"))) #add query sequence to list as Bio SeqRecord object
        for record in Bio.SeqIO.parse(cdhitoutput, "fasta"): #add other sequences to list as Bio SeqRecord objects
            self.out.append(record)
        if filename is not None:
            Bio.SeqIO.write(self.out, filename, "fasta")
        os.remove(cdhitinput)
        os.remove(cdhitoutput)
        os.remove(cdhitoutput+'.clstr')

#Run clustal omega alingment
#needs settings and sequences (as Bio SeqRecord objects)
#returns 'out'= aligned sequences Bio SeqRecord object
class runclustalo(object):
    def __init__(self, settings, sequences, filename=None):
        out_handle = StringIO.StringIO()
        Bio.SeqIO.write(sequences, out_handle, "fasta") #converts Bio SeqRecord objects from self.Bioout into fasta strings
        fastaseqs = out_handle.getvalue() #using StringIO to replace writting to a file
        RUNCLUSTAL = settings.CLUSTAL+' --iter='+str(settings.ALIGNMENTITERATIONS)+' -i - --outfmt=fa --force'
        print("\nAligning sequences using Clustal Omega. This can take a few minutes, espicially with many sequences")
        start = time.time()
        command = runbin.Command(RUNCLUSTAL)
        self.out = command.run(timeout=7200, stdin=fastaseqs)
        handle=StringIO.StringIO(self.out[1]) #turns the fasta string into a 'file like' object
        self.out = list(Bio.SeqIO.parse(handle, "fasta")) #reads the file-like object intoa list of Bio SeqRecord objects
        if command.status == -15:
            cleanexit('CLUSTAL TOOK TOO LONG! Try with fewer sequences.',keeptemp=settings.KEEPTEMPFILES)
        elif command.status != 0:
            cleanexit('CLUSTAL FAILED. I do not know why, check your inputs',keeptemp=settings.KEEPTEMPFILES)
        if filename is not None:
            print(filename)
            Bio.SeqIO.write(self.out, filename, "fasta")
        end = time.time()
        print 'Aligning sequences took '+str(int(end - start))+' seconds'