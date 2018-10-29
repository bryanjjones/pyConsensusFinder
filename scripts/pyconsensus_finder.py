#!/usr/bin/python
import _mypath
import os
import CF
import time
import argparse
import PDB
import urllib2
import sys

HOME = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
configfile = HOME+"/config/config.cfg"
defaults = {
            'File_Name' : "4eb0.fasta",
            'Email' : "None",
            'Maximum_Sequences' : "2000",
            'Blast_E_Value' : "1e-3",
            'Consensus_Threshold' : "0",#need to check if zero will give desired results. Might need to be 1
            'Consensus_Ratio' : "7",
            'Use_Complete_sequences' : "False",
            'Alignment_Iterations' : "1",
            'Maximum_Redundancy_Threshold' : "0.9",
            'Logging' : "False",
            'Keep_Temp_Files' : "False",
            'Chain' : 'A',
            'Angstrom' : "10",
            'Residue' : "56",
            'PDB_Name' : None,
            'Blast_binary' : HOME+"/binaries/blastp",
            'CDHIT_binary' : HOME+"/binaries/cd-hit",
            'ClustalO_binary' : HOME+"/binaries/clustalo-1.2.4-Ubuntu-x86_64",
            }


parser = argparse.ArgumentParser(description='Settings',prog='pyconsensus_finder.py',epilog='If no arguments specified on the command line, options will be read form the config file: '+configfile)
parser.add_argument('-q', '--query', metavar='FILENAME.FSTA',dest='FILENAME',type=str,default=defaults['File_Name'],help='query file to be analyzed')
parser.add_argument('-a', '--emailaddress',metavar='NAME@EXAMPLE.COM',dest='EMAIL',type=str,default=defaults['Email'],help='Entrez requires an email address to moitor usage')
parser.add_argument('-s', '--maxseqs',metavar='X',dest='MAXIMUMSEQUENCES',type=int,default=defaults['Maximum_Sequences'],help='Maximum sequences for BLAST search')
parser.add_argument('-e', '--evalue',metavar='1e-X',dest='BLASTEVALUE',type=float,default=defaults['Blast_E_Value'],help='Maximum e value for BLAST search')
parser.add_argument('-t', '--threshold',metavar='0.X',dest='CONSENSUSTHRESHOLD',type=float,default=defaults['Consensus_Threshold'],help='Minimum frequency for determining consensus')
parser.add_argument('--ratio',metavar='X',dest='RATIO',type=float,default=defaults['Consensus_Ratio'],help='Minimum ratio for determining consensus')
parser.add_argument('-p', '--partial',dest='USECOMPLETESEQUENCES',action='store_false',default=defaults['Use_Complete_sequences'],help='Use only matching partial sequences returned by BLAST, not complete sequences.')
parser.add_argument('-i', '--iter',metavar='X',dest='ALIGNMENTITERATIONS',type=int,default=defaults['Alignment_Iterations'],help='Number of clustal omega iterations')
parser.add_argument('-r', '--redundancy',metavar='0.X',dest='MAXIMUMREDUNDANCYTHRESHOLD',type=float,default=defaults['Maximum_Redundancy_Threshold'],help='Maximum identity for redundant sequence cutoff by CD-HIT')
parser.add_argument('-k', '--keeptemp',dest='KEEPTEMPFILES',action='store_true',default=defaults['Keep_Temp_Files'],help='Keep temporary files for troubleshooting')
parser.add_argument('-l', '--logging',dest='LOGGING',action='store_true',default=defaults['Logging'],help='Turn on logging for troubleshooting')
#options commented out for future use 
parser.add_argument('--chain',metavar="letter",dest="CHAIN",type=str,default=defaults['Chain'],help="Protein chain from PDB")
parser.add_argument('--residue',metavar="number",dest="RESIDUE",type=int,default=defaults['Residue'],help="Residue number from PDB")
parser.add_argument('--PDB',metavar="code",dest="PDB",type=str,default=defaults['PDB_Name'],help="Four letter PDB ID")
parser.add_argument('--angstroms',metavar="X",dest="ANG",type=float,default=defaults['Angstrom'],help="Distance in Angstroms")
parser.add_argument('--BLAST',type=str,default=defaults['Blast_binary'],help=argparse.SUPPRESS)
parser.add_argument('--CDHIT',type=str,default=defaults['CDHIT_binary'],help=argparse.SUPPRESS)
parser.add_argument('--CLUSTAL',type=str,default=defaults['ClustalO_binary'],help=argparse.SUPPRESS)
args = parser.parse_args()

if args.PDB and len(args.PDB) == 4:
    response = urllib2.urlopen('https://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList={}&compressionType=uncompressed'.format(args.PDB))
    html = response.read()
    with open(HOME+'/uploads/{}.fasta'.format(args.PDB), "w") as f:
        f.write(html)
    with open(HOME+'/uploads/{}.fasta'.format(args.PDB), "r") as f:
        for line in f:
            if line.strip('\n') in ("<head>"):
                CF.cleanexit('PDB code not found.')   
    args.FILENAME='{}.fasta'.format(args.PDB)

if args.PDB!="None" and args.FILENAME!="None.fasta":
    message = 'Entered both a PDB code and a default file name. Using the PDB code, {}'.format(args.PDB)
    print (message)


programstart = time.time()
#If command line variables are specified, use those.
if len(sys.argv) > 1:
    MainProgram=CF.CF(settings=args)
    if args.PDB!="None":
        PDBProgram=PDB.PDB(settings=args)
#If no command line variables are specified, use config file and defaults
else:
    print'reading settings from configfile ('+configfile+')'
    MainProgram=CF.CF(defaults=defaults,configfile=configfile)
    if args.PDB!="None":
        PDBProgram=PDB.PDB(defaults=defaults,configfile=configfile)
programend = time.time()
print '\nConsensus Finder Completed.'
os.rename(HOME+'/uploads/'+MainProgram.settings.FILENAME, HOME+'/completed/'+MainProgram.settings.FILENAME)
for i in MainProgram.output[:]:
    print(i)
print 'Your results are in the ./completed/ directory.'
print ''.join(MainProgram.warnings)
print 'Process took '+str(int(programend - programstart))+' seconds'
CF.cleanexit(0)
