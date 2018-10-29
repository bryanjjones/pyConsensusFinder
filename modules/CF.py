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
import httplib

HOME = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
#Main CF program, takes settings object and runs blast, cdhit, and clustalo. Then trims gaps, calculates counts, frequencies, and consensus at each position. Finally produces an output file with suggested mutations.
class CF(object):
    def __init__(self,defaults=None,configfile=None,settings=None,write=True):
        if settings is None:
            settings=setsettings(defaults,configfile)
        #do checks of settings
        print'query file: '+str(settings.FILENAME)
        print'email address: '+str(settings.EMAIL)
        print'BLAST maximum sequences: '+str(settings.MAXIMUMSEQUENCES)
        print'BLAST maximum e value: '+str(settings.BLASTEVALUE)
        print'threshold value for consensus: '+str(settings.CONSENSUSTHRESHOLD)
        print'ratio value for consensus: '+str(settings.RATIO)
        print'download complete sequences from Entrez: '+str(settings.USECOMPLETESEQUENCES)
        print'Clustal omega alignment itterations: '+str(settings.ALIGNMENTITERATIONS)
        print'threshold for removing redundant sequences: '+str(settings.MAXIMUMREDUNDANCYTHRESHOLD)
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
        if len(self.blast.versions)==0: #sanity check, make sure we have data
            cleanexit('Blast Returned zero hits', keeptemp=settings.KEEPTEMPFILES)
        self.cdhit=runcdhit(settings,self.blast.out)
        if len(self.cdhit.out)==0: #sanity check, make sure we have data
            cleanexit('CD-HIT returned zero sequences. Cannot continue without sequences.', keeptemp=settings.KEEPTEMPFILES)
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
        if configfile is not None:
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
            cleanexit("No query file specified. EXITING")
        if not os.path.isfile(HOME+'/uploads/'+settings.FILENAME):
            cleanexit('No query file. File '+settings.FILENAME+' does not exist in uploads directory.')
        if not 1 == (len(list(Bio.SeqIO.parse(HOME+'/uploads/'+settings.FILENAME, "fasta")))):
            #print(list(Bio.SeqIO.parse(HOME+'/uploads/'+settings.FILENAME, "fasta")))
            outputBase = '{}'.format(settings.PDB) # output.1.txt, output.2.txt, etc.
            input = open(HOME+'/uploads/'+settings.FILENAME, 'r').read().replace('>','#>').split('#')
            at = 1
            for lines in range(1, len(input), 1):
             # First, get the list slice
                outputData = input[lines:lines+1]
                # Now open the output file, join the new slice with newlines
                # and write it out. Then close the file.
                output = open(HOME+'/uploads/' + outputBase + '-' + str(at) + '.fasta', 'w+')
                output.write('\n'.join(outputData))
                output.seek(6)
                chainLetter = output.read(1)
                os.rename(HOME+'/uploads/'+ outputBase + '-' + str(at) + '.fasta',HOME+'/uploads/'+ outputBase + '-' + '{}'.format(chainLetter) + '.fasta')
                output.close()
                # Increment the counter
            if settings.CHAIN:
                settings.FILENAME=settings.PDB+'-{}'.format(settings.CHAIN)+'.fasta'
                message = "Using Chain {}".format(settings.CHAIN)
                print message
            else:
                message = "WARNING, MORE THAN ONE SEQUENCE FOUND IN "+settings.FILENAME+". Using only the first sequence."
                print message
                settings.FILENAME=settings.PDB+'-A'+'.fasta'
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
        numberofhits = len(self.versions)
        print('BLAST returned '+str(numberofhits)+' sequences.')        
        self.out=[] #initialize list for recording sequences
        if settings.USECOMPLETESEQUENCES is "True" and numberofhits != 0: #if use complete sequences is true, download sequences from Entrez, otherwise just use returned BLAST sequences
            start = time.time() #start timer for downloads
            print('\nDownloading '+str(numberofhits)+' complete sequences from NCBI.')
            Bio.Entrez.email=settings.EMAIL
            Bio.Entrez.tool = "Consensus Finder"
            def fetchseqs(ids, maxtries):
                retmax=len(ids)
                tries=1
                while tries <= maxtries:
                    try:
                        handle = Bio.Entrez.efetch(db="protein", id=ids, retmax=retmax, rettype="fasta", retmode="text") #retrieving all sequences from Entrez, db="sequences"
                        newseqs = list(Bio.SeqIO.parse(handle, "fasta"))
                        handle.close()
                        #assert len(newseqs)==retmax
                        time.sleep(1) #delay to avoid flooding entrez server
                        return newseqs
                    except httplib.HTTPException, e:
                        print('Network problem on try %s of %s: %s')  % (tries+1,maxtries,e)
                        tries += 1
                        if tries <= maxtries:
                            time.sleep(3**tries) #delay to avoid flooding entrez server
                            print('Trying again...')
                        else:
                            print('Giving up')
                            print('Try at a less busy time of day, or request fewer sequences')
                            cleanexit(message='Network problem trying to communicate with Entrez. Try again at a less busy time of day, or request fewer sequences')
                    except AssertionError:
                        print('Entrez failed to return all the sequences asked for on try %s of %s') % (tries+1,maxtries)
                        tries += 1
                        if tries <= maxtries:
                            time.sleep(3**tries) #delay to avoid flooding entrez server
                            print('Trying again...')
                        else:
                            print('Giving up')
                            print('Try at a less busy time of day, or request fewer sequences')
                            cleanexit(message='Network problem trying to communicate with Entrez. Try again at a less busy time of day, or request fewer sequences')
                    except:
                        print('Entrez faild to return the requsted sequences asked for on try %s of %s') % (tries+1,maxtries)
                        tries += 1
                        if tries <= maxtries:
                            time.sleep(3**tries) #delay to avoid flooding entrez server
                            print('Trying again...')
                        else:
                            print('Giving up')
                            print('Try at a less busy time of day, or request fewer sequences')
                            cleanexit(message='Network problem trying to communicate with Entrez. Try again at a less busy time of day, or request fewer sequences')
            #break the list of sequences into chunks
            chunksize=500
            chunks=[ self.versions[x:x+chunksize] for x in xrange(0, len(self.versions), chunksize)]
            for index, chunk in enumerate(chunks):
                versionlist=",".join(chunk)
                print('Downloading sequences %s to %s of %s total') % (index*chunksize+1, index*chunksize+len(chunk), len(self.versions))
                self.out=self.out+fetchseqs(versionlist, 3) #fetch a chunk of sequences and add it on to the list
            print('Downloaded '+str(len(self.out))+' sequences')
            #sanity check to make sure we got all the requested sequences.
            try:
                assert numberofhits==len(self.out)
            except:
                print('Did not get the expected '+str(numberofhits)+' sequences from Entrez.')
                print('Continuing with only '+str(len(self.out))+' sequences.')
                self.warnings.append('WARNING! Entrez only returned '+len(self.out)+' sequences of '+str(numberofhits)+' requested. Results based on '+len(self.out)+' sequences.')
            end = time.time()
            print('Downloading sequences took '+str(int(end - start))+' seconds.')
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
        print("\nRemoving redundant sequences using CD-HIT.")
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
        print(str(len(self.out))+' sequences after removind redundants.')
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
        print("\nAligning "+str(len(sequences))+" sequences using Clustal Omega. This can take a few minutes, espicially with many sequences")
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
