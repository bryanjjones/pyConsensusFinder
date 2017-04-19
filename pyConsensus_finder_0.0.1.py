#!/usr/bin/python

import numpy as np # need numpy for arrays and the like
#import subprocess # need to run binaries (blast, CD-hit, clustal)
import ConfigParser# need to read Config file
import time
import threading
import subprocess
import traceback
import shlex


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

def printf(x):
	print x
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
Config.read("./uploads/CONFIGFILE")

FILENAME=updatefromconfig('BasicSettings', 'filename')
if not FILENAME:
	print FILENAME
	print "No query file specified in CONFIGFILE. EXITING"
	sys.exit(0)
#overwrite defaults if settings specified in CONFIGFILE
MAXIMUMSEQUENCES=updatefromconfig('BlastSettings', 'MAXIMUMSEQUENCES')
BLASTEVALUE=updatefromconfig('BlastSettings', 'BLASTEVALUE')
CONSENSUSTHRESHOLD=updatefromconfig('AlignmentSettings', 'ConsensusThreshold')
USECOMPLETESEQUENCES=updatebooleanconfig('AlignmentSettings', 'UseCompleteSequences')
ALIGNMENTITERATIONS=updatefromconfig('AlignmentSettings', 'AlignmentIterations')
MAXIMUMREDUNDANCYTHRESHOLD=updatefromconfig('AlignmentSettings', 'MaximumRedundancyThreshold')
LOGGING=updatebooleanconfig('TroubleShooting', 'Logging')
KEEPTEMPFILES=updatebooleanconfig('TroubleShooting', 'KeepTempFiles')

#location of binaries to call
BLAST="./binaries/blastp" #change this if another version is installed locally
CDHIT="./binaries/cd-hit" #change this if another version is installed locally
CLUSTAL="./binaries/clustalo-1.2.0-Ubuntu-x86_64" #change this if another version is installed locally


#Run Blast

class Command(object):
	"""
	Enables to run subprocess commands in a different thread with TIMEOUT option.
	Adapted from kirpit's code http://gist.github.com/kirpit/1306188 and jcollado's solution:
	http://stackoverflow.com/questions/1191374/subprocess-with-timeout/4825933#4825933
	"""
	command = None
	process = None
	status = None
	output, error = '', ''

	def __init__(self, command):
		if isinstance(command, basestring):
			command = shlex.split(command)
		self.command = command

	def run(self, timeout=None, **kwargs):
		""" Run a command then return: (status, output, error). """
		def target(**kwargs):
			try:
				self.process = subprocess.Popen(self.command, **kwargs)
				self.output, self.error = self.process.communicate()
				self.status = self.process.returncode
			except:
				self.error = traceback.format_exc()
				self.status = -1
		# default stdout and stderr
		if 'stdout' not in kwargs:
			kwargs['stdout'] = subprocess.PIPE
		if 'stderr' not in kwargs:
			kwargs['stderr'] = subprocess.PIPE
		# thread
		thread = threading.Thread(target=target, kwargs=kwargs)
		thread.start()
		thread.join(timeout)
		if thread.is_alive():
			self.process.terminate()
			print 'Taking too long'
			print 'Terminating process'
			thread.join()
		return self.status, self.output, self.error

RUNBLAST = BLAST+' -db nr -query ./uploads/'+FILENAME+' -evalue '+BLASTEVALUE+' -max_target_seqs '+MAXIMUMSEQUENCES+' -outfmt "6 sacc sseq pident" -remote' 

start = time.time()
command = Command(RUNBLAST)
BLASTOUT=command.run(timeout=900)
end = time.time()
print 'Time elapsed (seconds)'
print(end - start)
print 'BLAST Results'
print BLASTOUT[1]

# BLASTOUT[1] needs to become an array
# when complete sequences =1, curl BLASTOUT[1] VERSIONS (first entry on line), ideally replace existing sequences with curled (complete sequences)



