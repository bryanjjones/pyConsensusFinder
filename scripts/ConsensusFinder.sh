#!/bin/bash
#bash wrapper for pyConsensusFinder. This will populate the required config.cfg file using the given command line arguments and run pyConsensusFinder.
#Default Values
MAXSEQ="2000"
EVALUE="1e-3"
THRESHOLD="0"
RATIO="7"
USE_COMPLETE_SEQ="1"
ALIGNMENTITER="1"
ID_REDUNDANCY="0.9"
LOGGING="0"
KEEPTEMP="0"
EMIAL="ConsensusFinderUser@xxx.com"
CFG="./config/config.cfg"
pyConsensusFinder="./scripts/pyconsensus_finder_0.0.2.py"

#get values form command line arguments:
while getopts ":hs:e:t:ci:r:klq:a:o:" opt; do
	case $opt in
		h)	
			usage
			exit 1
			;;
		s)
			MAXSEQ="$OPTARG"
			echo "Maximum BLAST sequences: $OPTARG" >&2
			;;
		e)
			EVALUE="$OPTARG"
			echo "Maximum BLAST e value: $OPTARG" >&2
			;;
		t)
			THRESHOLD="$OPTARG"
			RATIO=
			echo "Minimum threshold for consensus: $OPTARG" >&2
			;;
		c)
			USE_COMPLETE_SEQS="0"
			echo "Use complete sequences: off" >&2
			;;
		i)
			ALIGNMENTITER="$OPTARG"
			echo "Clustal Omega alingment iterations: $OPTARG" >&2
			;;
		r)
			ID_REDUNDANCY="$OPTARG"
			echo "Maximum identity for redundant sequence cutoff: $OPTARG" >&2
			;;			
		l)	
			LOGGING="1"
			echo "logging on" >&2
			;;
		k)	
			KEEPTEMP="1"
			echo "Temporary files will not be deleted" >&2
			;;
		q)
			SOURCE="$OPTARG"
			echo "query file: $OPTARG" >&2
			;;		
		a)
			EMIAL="$OPTARG"
			echo "Email Address: $OPTARG" >&2
			;;
		o) 
		    RATIO="$OPTARG"
		    echo "Minimum ratio for consensus: $OPTARG" >&2
		    ;;
		?)
			echo "Invalid option: -$OPTARG. Run Consensus Finder with -h flag for help." >&2
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument. Run Consensus Finder with -h flag for help." >&2
			exit 1
			;;
	esac
done
#build the config file
echo "[BasicSettings]" > $CFG
echo "FileName: "$SOURCE >> $CFG
echo "Email: "$EMIAL >> $CFG
echo "" >> $CFG
echo "[BlastSettings]" >> $CFG
echo "MaximumSequences: "$MAXSEQ >> $CFG
echo "BlastEValue: "$EVALUE >> $CFG
echo "" >> $CFG
echo "[AlignmentSettings]" >> $CFG
echo "UseCompleteSequences: "$USE_COMPLETE_SEQS >> $CFG 
echo "AlignmentIterations: "$ALIGNMENTITER >> $CFG
echo "MaximumRedundancyThreshold: "$ID_REDUNDANCY >> $CFG
echo "ConsensusRatio: "$RATIO >>$CFG
echo "ConsensusThreshold: "$THRESHOLD >> $CFG
echo "" >> $CFG
echo "[TroubleShooting]" >> $CFG
echo "Logging: "$LOGGING >> $CFG
echo "KeepTempFiles: "$KEEPTEMP >> $CFG
echo "" >> $CFG

python $pyConsensusFinder