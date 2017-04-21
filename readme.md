Consensus Finder User Guide
Why: 
Consensus protein sequences are useful for numerous applications. Often, mutating a protein to be 
more like the consensus of homologs will often increase the stability of a protein, allowing it to  
function at higher temperatures, and have better soluble expression when expressed recombinatly in 
various hosts. "Consensus Finder" will help identify the consensus sequence and find potentially 
stabilizing mutations.
What it does:
Consensus Finder will take given your given protein sequence, find similar sequences 
from the NCBI database, align them, remove redundant/highly similar sequences, trim alignments to 
the size of the original query, and analyze consensus. Output is trimmed alignment, consensus 
sequence, frequency and count tables for amino acids at each position, as well as a list of 
suggested mutations to consensus that may be stabilizing.

How to:
Choose File
Choose a file containing your protein sequence. Needs to be a plain text file in FASTA format containing the sequence of your protein of interest (protein sequence, not DNA). 
If you don’t have that handy, you can download a FASTA file for proteins from NCBI. 
Go to https://www.ncbi.nlm.nih.gov/protein/
Search for your protein by name or GI, ACCESSION, UniProtKB/Swiss-Prot or any identifier you have.
In the list of results, check the box by the one you want.
Click the “Send To” drop down menu, check file, and select format “FASTA”, and click “Create File”.
Submit
Press Submit and wait. Click the link to see if it’s done or wait for an email if you provided an email address. Usually takes 10-30 minutes, but sometimes takes longer.
Results
The results page shows you a list of mutations to make which are likely to increase the stability of your protein. The are listed in order of most conserved (and best bet for stabilization) to least. So if you get a long list, start at the top. If you have your gene and want to use site-directed mutagenesis to change the indicated residues, a site like primerX can aid you in designing primers to make these changes (http://www.bioinformatics.org/primerx/cgi-bin/DNA_1.cgi). 
If you don’t get any results, see the Options section below to change the parameters.

Advanced Features
Optional operations
You can change extra options by clicking “show” near “Optional operations” before you submit your job.

2000 Set maximum sequences for BLAST search (Range: 10 - 10000)
1e-3  Set maximum e value for BLAST search (Range: 1e-30 - 1e-1)
.6      Conservation threshold for suggesting mutations (Range: .05 - .99)
[]      Use only matched portions, not complete sequences
1      Iterations of ClastalW alignments (Range: 1 - 5)
.9     CD-Hit redundancy (Range: .5 - 1.0)

BLAST parameters can be changed to adjust maximum number of sequences and adjust 
the maxumum e value. One or the other will limit the number of sequences returned. For both lower 
values will tend to return fewer more similar results. Default for “Set maximum sequences ...” is 2000, with reasonable range would be 10-10,000, BLAST searches will occasionally time out if too many sequences are requested. Default for "Set maximum e value ..."  is 1e-3 (scientific notation), reasonable range would be 1e-30 to 1e-1.

Changing the “Conservation threshold” will change how many mutation suggestions 
are returned. The default is 0.6, reasonable range 0.05-0.99, with lower numbers returning more suggestions.

By default complete sequences will individually be downloaded from NCBI. This takes a lot of time, 
especially with a lot of BLAST hits.  Checking the box for “Use only matched portions...” uses only the partial sequences matching the target in the BLAST result, this can make the program run faster.

Increasing the number of “Iterations of ClustalW alingments” will, in theory, will give better 
quality alignments, but it takes much longer and can give rise to other issues. Default is just 1 iteration, 
reasonable options: integers from 1-5.

The maximum threshold for eliminating redundant sequences with CD-HIT is defaulted to 0.9, which will remove any sequences with over 90% identity to prevent over sampling of over represented groups of proteins. Reasonable range would be 0.7-1.0, with a setting of 1.0 keeping all redundant sequences.

Copyright 2016 Bryan J. Jones (bryanjjones@gmail.com)
Consensus Finder can be freely copied and distributed under the GNU General Public 
License version 2 (GPLv2) or later.

Citations:  
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


Additional Outputs. i.e.  “Download full results”
There are other useful results in the .zip file you get! There are 6 files in your full results, one is your original query sequence. The others are:
[your_query]_mutations.txt
This is a text file that contains the same information displayed on the results web page, i.e. a list of suggested mutations.

[your_query]_consensus.fst
This file contains the consensus sequence across the full length of your query protein. It is a plain text file in FASTA format. Each position is the most common amino acid at that position from all the homologous compared (regardless of how strong of a consensus there is).

[your_query]_trimmed_alignment.fst
This file is a text file in FASTA format containing an alignment of all the representative sequences used to calculate the consensus. This alignment has the redundant sequences removed, and is trimmed to the length of your query sequence. That means that any insertions relative to the query sequence are deleted out, and  any deletions are replaced with gaps (“-”).  The names are the VERSION numbers, so you can cross reference the individual sequences or look them up at NCBI.

[your_query]_counts.csv
This is a table (you can open it as a spread sheet in Excel) that shows the count of each amino acid at each position. The first column is a list of the 20 amino acids, plus gap (“-”), and other (this can be unknown amino acids or non-trad
[your_query]_frequencies.csv
F.A.Q.
I got no/too few suggested mutations what’s wrong?
Try changing some of the Advanced Options. You could degrease the “Conversation Threshold” so that a weaker consensus will still suggest a mutation. Changing the “maximum sequences” or “maximum e” value of the BLAST search can also affect the number of suggested mutations. Either increasing or decreasing the two BLAST options can result in more or fewer suggested mutations. More results can result in weaker consensus so fewer positions will be above your “Conversation Threshold”, but fewer results are more likely to have the same residues as your query, so no mutation can be made.

Why do I get so many suggested mutations?
If you got more suggestions than you want to actually make, just use the suggestions from the top of the list. Since the mutations are sorted with the “best” (i.e. most conserved) suggestions at the top. You can also increase the “Conservation Threshold” or change the BLAST settings to change the number of  suggested mutations.

