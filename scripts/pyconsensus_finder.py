#!/usr/bin/python
import _mypath
import os
import CF
import time

HOME = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))

configfile = HOME+"/config/config.cfg"
defaults = {
            'File_Name' : "None",
            'Email' : "None",
            'Maximum_Sequences' : "2000",
            'Blast_E_Value' : "1e-3",
            'Consensus_Threshold' : "0",#need to check if zero will give desired results. Might need to be 1
            'Consensus_Ratio' : "7",
            'Use_Complete_sequences' : "True",
            'Alignment_Iterations' : "1",
            'Maximum_Redundancy_Threshold' : "0.9",
            'Logging' : "False",
            'Keep_Temp_Files' : "False",
            'Chain' : 'A',
            'Angstrom' : "1",
            'Residue' : "0",
            'PDB_Name' : "None",
            }

settings=CF.setsettings(defaults,configfile)
#location of binaries to call
settings.BLAST=HOME+"/binaries/blastp"
settings.CDHIT=HOME+"/binaries/cd-hit"
settings.CLUSTAL=HOME+"/binaries/clustalo-1.2.4-Ubuntu-x86_64"

programstart = time.time()
#do checks of settings
warnings=[]
#Run CF checks to check all settings, return any warnings to warnings variable
warnings=warnings+CF.checks(settings).warnings
MainProgram=CF.CF(settings,warnings=warnings)
programend = time.time()
print '\nConsensus Finder Completed.'
os.rename(HOME+'/uploads/'+settings.FILENAME, HOME+'/completed/'+settings.FILENAME)
for i in MainProgram.output[:]:
    print(i)
print 'Your results are in the ./completed/ directory.'
print ''.join(MainProgram.warnings)
print 'Process took '+str(int(programend - programstart))+' seconds'
CF.cleanexit(0)