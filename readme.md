# Consensus Finder

Consensus protein sequences are useful for numerous applications. Often, mutating a protein to be more like the consensus of homologs will often increase the stability of a protein, allowing it to function at higher temperatures, and have better soluble expression when expressed recombinatly in various hosts. "Consensus Finder" will help identify the consensus sequence and find potentially stabilizing mutations.

## What it does
Consensus Finder will take given your given protein sequence, find similar sequences from the NCBI database, align them, remove redundant/highly similar sequences, trim alignments to the size of the original query, and analyze consensus. Output is trimmed alignment, consensus sequence, frequency and count tables for amino acids at each position, as well as a list of suggested mutations to consensus that may be stabilizing.

## Web server
Here is a working web implementation:  http://kazlab.umn.edu/

## Setting up a local installation
pyConsensus Finder should run on most 64bit linux systems with python (2.7), numpy, and Biopython installeon](#nsensus Finder should run on most 64bit linux systems with python (2.7), numpy, and Biopython install-on)
pyConaries" folder) for: blastp (2.6.0+), CD-HIT (4.7), & Clustal Omega (1.2.4).

## How to run pyConsensus Finder
### Choose File
Choose a file containing your protein sequence. Needs to be a plain text file in FASTA format containing the sequence of your protein of interest (protein sequence, not DNA). 
If you don’t have that handy, you can download a FASTA file for proteins from NCBI. 
Go to https://www.ncbi.nlm.nih.gov/protein/
Search for your protein by name or GI, ACCESSION, UniProtKB/Swiss-Prot or any identifier you have.
In the list of results, check the box by the one you want.
Click the “Send To” drop down menu, check file, and select format “FASTA”, and click “Create File”.

### Run
Place your query file in the ./uploads directory.
Specify your query file name, email address (required by NCBI), and any optional settings (see below) in the config.cfg file in the ./config subdirectory. See 'EXAMPLEconfig.cfg" for formatting

To run, on the command line run

    python ./scripts/pyconsensus_finder_0.0.2.py

### Results
Note: Blast results can vary slightly between runs, even with the same query and settings. This can result in slightly different results from repeated identical Consensus Finder runs.
Results
The results shows you a list of mutations to make which are likely to increase the stability of your protein. The are listed in order of most conserved (and best bet for stabilization) to least. So if you get a long list, start at the top. If you have your gene and want to use site-directed mutagenesis to change the indicated residues, a site like primerX can aid you in designing primers to make these changes (http://www.bioinformatics.org/primerx/cgi-bin/DNA_1.cgi). 
If you don’t get any results, see the Options section below to change the parameters.

### Optional settings for config.cfg file

| Setting    | Default Value |  Explanation |
|:----------:|-------------|-------------|
| **`[BasicSettings]`** | | |
| `FileName:` | `FILENAME.fst` | Specifying filename for the query file. Mandantory. |
| `Email:` | `example@xxx.edu` | Entrez requires an email to monitor server usage. |
| **`[BlastSettings]`** | | |
| `MaximumSequences:` | `2000` | Set maximum sequences for BLAST search (Range: 10 - 10000) |
| `BlastEValue:` | `1e-3` | Set maximum e value for BLAST search (Range: 1e-30 - 1e-1) |
| **`[AlignmentSettings]`** | |
| `UseCompleteSequences:` | `yes` | Use only matched portions, not complete sequences |
| `AlignmentIterations:` | `1` | Iterations of ClastalW alignments (Range: 1 - 5) |
| `MaximumRedundancyThreshold:` | `.9` | CD-Hit redundancy (Range: .5 - 1.0) |
| `ConsensusRatio:` | `7` | Conservation threshold ratio of consensus residue to query residue for suggesting mutations (Range: 1 - 100) |
| `ConsensusThreshold:` | `.6` | Conservation threshold of consensus residue frequency for suggesting mutations (Range: .05 - .99)|

BLAST parameters can be changed to adjust maximum number of sequences and adjust the maxumum e value. One or the other will limit the number of sequences returned. For both, lower values will tend to return fewer more similar results. Default for “MaximumSequences” is 2000, with reasonable range would be 10-10,000, BLAST searches will occasionally time out if too many sequences are requested. Default for "BlastEValue"  is 1e-3 (scientific notation), reasonable range would be 1e-30 to 1e-1.

The threshold for defining consensus can by set by “MaximumRedundancyThreshold”, "ConsensusRatio", or both.  The default for "MaximumRedundancyThreshold" is 0.6, reasonable range 0.05-0.99, with lower numbers returning more suggestions.
The default for "ConsensusRatio" is 7, reasonable range is 1-100, with lower numbers returning more suggestions.

By default complete sequences will individually be downloaded from NCBI. This takes a lot of time, especially with a lot of BLAST hits.  Changing “UseCompleteSequences” to "no" uses only the partial sequences matching the target in the BLAST result, this can make the program run faster.

Increasing the number of “AlignmentIterations” will, in theory, will give better quality alignments, but it takes much longer and can give rise to other issues. Default is just 1 iteration, reasonable options: integers from 1-5.

The maximum threshold for eliminating redundant sequences with CD-HIT, "MaximumRedundancyThreshold", is defaulted to 0.9, which will remove any sequences with over 90% identity to prevent over sampling of over represented groups of proteins. Reasonable range would be 0.7-1.0, with a setting of 1.0 keeping all redundant sequences.

## Additional results files

There are other useful results in the 5 files in your results:

[your_query]_mutations.txt

This is a text file that contains the same information displayed on the results web page, i.e. a list of suggested mutations.

[your_query]_consensus.fst

This file contains the consensus sequence across the full length of your query protein. It is a plain text file in FASTA format. Each position is the most common amino acid at that position from all the homologous compared (regardless of how strong of a consensus there is).

[your_query]_trimmed_alignment.fst

This file is a text file in FASTA format containing an alignment of all the representative sequences used to calculate the consensus. This alignment has the redundant sequences removed, and is trimmed to the length of your query sequence. That means that any insertions relative to the query sequence are deleted out, and  any deletions are replaced with gaps (“-”).  The names are the VERSION numbers, so you can cross reference the individual sequences or look them up at NCBI.

[your_query]_counts.csv

This is a table (you can open it as a spread sheet in Excel) that shows the count of each amino acid at each position. The first column is a list of the 20 amino acids, plus gap (“-”), and other (this can be unknown amino acids or non-trad

[your_query]_frequencies.csv

## Authors
Copyright 2016 Bryan J. Jones (bry@njjon.es)
## License
Consensus Finder can be freely copied and distributed under the GNU General Public 
License version 2 (GPLv2) or later.

## Acknowledgments

Weizhong Li, Lukasz Jaroszewski & Adam Godzik. "Clustering of highly homologous sequences to reduce 
  thesize of large protein database",  Bioinformatics, (2001) 17:282-283

Weizhong Li, Lukasz Jaroszewski & Adam Godzik. "Tolerating some redundancy significantly speeds up 
  clustering of large protein databases", Bioinformatics, (2002) 18:77-82

Sievers F, Wilm A, Dineen D, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, 
  Thompson JD, Higgins DG. "Fast, scalable generation of high-quality protein multiple sequence 
  alignments using Clustal Omega." Mol Syst Biol. 2011 Oct 11;7:539. doi: 10.1038/msb.2011.75. 
  PMID: 21988835.

Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) 
  "BLAST+: architecture and applications." BMC Bioinformatics 10:421.

## F.A.Q.

I got no/too few suggested mutations what’s wrong?

Try changing some of the Advanced Options. You could degrease the “Conversation Threshold” so that a weaker consensus will still suggest a mutation. Changing the “maximum sequences” or “maximum e” value of the BLAST search can also affect the number of suggested mutations. Either increasing or decreasing the two BLAST options can result in more or fewer suggested mutations. More results can result in weaker consensus so fewer positions will be above your “Conversation Threshold”, but fewer results are more likely to have the same residues as your query, so no mutation can be made.


Why do I get so many suggested mutations?

If you got more suggestions than you want to actually make, just use the suggestions from the top of the list. Since the mutations are sorted with the “best” (i.e. most conserved) suggestions at the top. You can also increase the “Conservation Threshold” or change the BLAST settings to change the number of  suggested mutations.


Note: Blast results can vary slightly between runs, even with the same query and settings. This can result in slightly different results from repeated identical Consensus Finder runs.
