#!/usr/bin/python
import _mypath
import os
import CF
import time
import argparse
import sys

HOME = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
configfile = HOME+"/config/config.cfg"
parser = argparse.ArgumentParser(description='Settings',prog='pyconsensus_finder.py',epilog='If no arguments specified on the command line, options will be read form the config file: '+configfile)
parser.add_argument('-q', '--query', metavar='FILENAME.FSTA',dest='FILENAME',type=str,help='query file to be analyzed')
parser.add_argument('-a', '-emailaddress',metavar='NAME@EXAMPLE.COM',dest='EMAIL',type=str,help='Entrez requires an email address to moitor usage')
parser.add_argument('-s', '-maxseqs',metavar='X',dest='MAXIMUMSEQUENCES',type=int,help='Maximum sequences for BLAST search')
parser.add_argument('-e', '--evalue',metavar='1e-X',dest='BLASTEVALUE',type=float,help='Maximum e value for BLAST search')
parser.add_argument('-t', '--threshold',metavar='0.X',dest='CONSENSUSTHRESHOLD',type=float,help='Minimum frequency for determining consensus')
parser.add_argument('--ratio',metavar='X',dest='RATIO',type=float,help='Minimum ratio for determining consensus')
parser.add_argument('-p', '--partial',dest='USECOMPLETESEQUENCES',action='store_false',help='Use only matching partial sequences returned by BLAST, not complete sequences.')
parser.add_argument('-i', '--iter',metavar='X',dest='ALIGNMENTITERATIONS',type=int,help='Number of clustal omega iterations')
parser.add_argument('-r', '--redundancy',metavar='0.X',dest='MAXIMUMREDUNDANCYTHRESHOLD',type=float,help='Maximum identity for redundant sequence cutoff by CD-HIT')
parser.add_argument('-k', '--keeptemp',dest='KEEPTEMPFILES',action='store_true',help='Keep temporary files for troubleshooting')
parser.add_argument('-l', '--logging',dest='LOGGING',action='store_true',help='Turn on logging for troubleshooting')
#parser.add_argument('--chain',metavar="letter",dest="CHAIN",type=str,help="Protein chain from PDB")
#parser.add_argument('--residue',metavar="number",dest="RESIDUE",type=int,help="Residue number from PDB")
#parser.add_argument('--PDB',metavar="code",dest="PDB",type=str,help="Four letter PDB ID")
#parser.add_argument('--angstroms',metavar="X",dest"ANG",type=float,help="Distance in Angstroms")
parser.add_argument('--BLAST',type=str,help=argparse.SUPPRESS)
parser.add_argument('--CDHIT',type=str,help=argparse.SUPPRESS)
parser.add_argument('--CLUSTAL',type=str,help=argparse.SUPPRESS)
args = parser.parse_args()
print args.File_Name

if len(sys.argv) > 1:
   configfile=None
else:
    print'reading settings from configfile ('+configfile+')'
