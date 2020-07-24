#!/usr/bin/python
import sys #need sys to use system variables
import numpy as np # need numpy for arrays and the like
import Bio.SeqIO
import Bio.Seq
 
#Array of single letter amino acid cods for use in arrays. 
IDS = np.zeros([22,1],dtype=object)
IDS[:,0]=['G', 'P', 'A', 'V', 'L', 'I', 'M', 'C', 'F', 'Y', 'W', 'H', 'K', 'R', 'Q', 'N', 'E', 'D', 'S', 'T', '-', 'other']
#make a list of data types for suggested mutation tuples wt= wild type (given) residue, res=residue number, sug= suggested mutation (consensus residue), freq= frequency of consensus residue (0-1), wtfreq= frequency of wild type (given) residue (0-1)  
TYPES = [('wt', 'S1'), ('res', int), ('sug', 'S1'),('freq', float),('wtfreq', float)] 
class res_to_change(object):
    def __init__(self, wt, res, sug, freq, wtfreq):
        self.wt=wt
        self.res=res
        self.sug=sug
        self.freq=freq
        self.wtfreq=wtfreq

#Given an alingment as a Bio SeqRecord object,
#Returns an array of sequences with each amino acid as an element.
def aaaray(sequences, filename=None):
    entries = len(sequences) #how many sequences do we have
    LENGTH = len(sequences[0].seq)
    #make array from file with each AA as an element
    AAA = np.zeros((entries,LENGTH), dtype=np.str)
    for i in range(entries):
        x=str(sequences[i].seq)
        for position in range(len(x)):
            AAA[i,position] = (x[position])
    for index in range(len(AAA[0,:])) : #make sure everything is upper case
        AAA[0,index] = AAA[0,index].upper()
    return AAA

#Given an alignment of sequences as arrays with each amino acid as an element
#Returns array of amino acids with any gaps in the first sequence deleted from all sequences, i.e. all trimmed to the length of the first sequence. If filename and sequences (for sequence ids, given as list of SeqRecord objects) are both given, will write fasta formatted file with names from sequences.
def trimmer(AAA, sequenceids=None, filename=None):
    print('\nTrimming gaps from alignment')
    entries,LENGTH=(AAA.shape)
    # delete all positions that corrispond to gaps in the first sequence
    for index in range(LENGTH):#at each position
        if AAA[0,(LENGTH-1-index)] == "-":#LENGTH-1-index will start at the end, -1 to account for 0 based indexing, and find gaps
            AAA = np.delete(AAA, (LENGTH-1-index), 1)#delete the gaps

    # if filename and sequences given, format array of amino acids as list of fasta formatted sequences and save
    if filename and sequenceids:
        for index in range(len(AAA[:,0])): #alternate ">[GI number]" and sequences
            sequenceids[index].seq = Bio.Seq.Seq(''.join(AAA[(index),:]))
        Bio.SeqIO.write(sequenceids, filename, "fasta")
    return AAA

#Returns an array of amino acid counts given an array of aligned sequences with each element a single position.
#If given a filename, counts array is exported as a csv.
def aacounts(AAA, filename=None): 
    COUNTS = np.zeros([22, len(AAA[0,:])],dtype=object) #makes array the length of the alingment with 22 rows (20AAs + "-" + "other")
    for index in range(len(AAA[0,:])): # for each position along the alingment, count occourances of each AA
        COUNTS[0,index]=AAA[:,index].tolist().count("G")
        COUNTS[1,index]=AAA[:,index].tolist().count("P")
        COUNTS[2,index]=AAA[:,index].tolist().count("A")
        COUNTS[3,index]=AAA[:,index].tolist().count("V")
        COUNTS[4,index]=AAA[:,index].tolist().count("L")
        COUNTS[5,index]=AAA[:,index].tolist().count("I")
        COUNTS[6,index]=AAA[:,index].tolist().count("M")
        COUNTS[7,index]=AAA[:,index].tolist().count("C")
        COUNTS[8,index]=AAA[:,index].tolist().count("F")
        COUNTS[9,index]=AAA[:,index].tolist().count("Y")
        COUNTS[10,index]=AAA[:,index].tolist().count("W")
        COUNTS[11,index]=AAA[:,index].tolist().count("H")
        COUNTS[12,index]=AAA[:,index].tolist().count("K")
        COUNTS[13,index]=AAA[:,index].tolist().count("R")
        COUNTS[14,index]=AAA[:,index].tolist().count("Q")
        COUNTS[15,index]=AAA[:,index].tolist().count("N")
        COUNTS[16,index]=AAA[:,index].tolist().count("E")
        COUNTS[17,index]=AAA[:,index].tolist().count("D")
        COUNTS[18,index]=AAA[:,index].tolist().count("S")
        COUNTS[19,index]=AAA[:,index].tolist().count("T")
        COUNTS[20,index]=AAA[:,index].tolist().count("-") #empty spaces
        COUNTS[21,index]=(len(AAA[:,index]) - sum(COUNTS[:,index].tolist())) #other, not counted above
    IDCOUNTS = np.hstack((IDS,COUNTS)) #make list with AA counts and names of AAs
    header=list(1,range(len(AAA[0,:])))
    header=["AA number"]+header
    IDCOUNTS = np.vstack((header,IDCOUNTS))
    if filename:
        np.savetxt((filename),IDCOUNTS,delimiter=",",fmt="%s") #save file with AA names and counts
    return COUNTS

#Returns a frequency array from an array of amino acid counts.
#Frequencies represented as a decimal.
#If given a filename, frequency array is exported as a csv.
def aafrequencies(COUNTS, filename=None):
    print('\nCalculating amino acid frequencies')
    FREQS = np.zeros_like(COUNTS) #make an array for calculating frequencies of each AA
    FREQS = np.float64(FREQS) #it needs to be numbers
    for index in range(len(FREQS[0,:])): # calculate the frequencey of each AA as [occurrences]/[occurrences of all AAs], "-" and "other" not counted in total
        FREQS[0,index]=np.float64(COUNTS[0,index])/sum(np.float64(COUNTS[:20,index]))
        FREQS[1,index]=np.float64(COUNTS[1,index])/sum(np.float64(COUNTS[:20,index]))
        FREQS[2,index]=np.float64(COUNTS[2,index])/sum(np.float64(COUNTS[:20,index]))
        FREQS[3,index]=np.float64(COUNTS[3,index])/sum(np.float64(COUNTS[:20,index]))
        FREQS[4,index]=np.float64(COUNTS[4,index])/sum(np.float64(COUNTS[:20,index]))
        FREQS[5,index]=np.float64(COUNTS[5,index])/sum(np.float64(COUNTS[:20,index]))
        FREQS[6,index]=np.float64(COUNTS[6,index])/sum(np.float64(COUNTS[:20,index]))
        FREQS[7,index]=np.float64(COUNTS[7,index])/sum(np.float64(COUNTS[:20,index]))
        FREQS[8,index]=np.float64(COUNTS[8,index])/sum(np.float64(COUNTS[:20,index]))
        FREQS[9,index]=np.float64(COUNTS[9,index])/sum(np.float64(COUNTS[:20,index]))
        FREQS[10,index]=np.float64(COUNTS[10,index])/sum(np.float64(COUNTS[:20,index]))
        FREQS[11,index]=np.float64(COUNTS[11,index])/sum(np.float64(COUNTS[:20,index]))
        FREQS[12,index]=np.float64(COUNTS[12,index])/sum(np.float64(COUNTS[:20,index]))
        FREQS[13,index]=np.float64(COUNTS[13,index])/sum(np.float64(COUNTS[:20,index]))
        FREQS[14,index]=np.float64(COUNTS[14,index])/sum(np.float64(COUNTS[:20,index]))
        FREQS[15,index]=np.float64(COUNTS[15,index])/sum(np.float64(COUNTS[:20,index]))
        FREQS[16,index]=np.float64(COUNTS[16,index])/sum(np.float64(COUNTS[:20,index]))
        FREQS[17,index]=np.float64(COUNTS[17,index])/sum(np.float64(COUNTS[:20,index]))
        FREQS[18,index]=np.float64(COUNTS[18,index])/sum(np.float64(COUNTS[:20,index]))
        FREQS[19,index]=np.float64(COUNTS[19,index])/sum(np.float64(COUNTS[:20,index]))
        FREQS[20,index]=np.float64(COUNTS[20,index])/sum(np.float64(COUNTS[:,index])) #frequency of gaps "-" as fraction of all seqs
        FREQS[21,index]=np.float64(COUNTS[21,index])/sum(np.float64(COUNTS[:,index])) #frequency of gaps "other" as fraction of all seqs
 
    #IDS=aaletters()
    IDFREQS = np.hstack((IDS,FREQS)) #make list with names and AA frequencies
    header=list(range(1,len(IDFREQS[0,:])))
    header=["AA number"]+header
    IDFREQS = np.vstack((header,IDFREQS))
    if filename:
        np.savetxt((filename),IDFREQS,delimiter=",",fmt="%s") #save file with AA names and frequencies
    return FREQS

#Calculates the consensus sequence from given amino acid frequency array.
#Returns consensus sequence as an array of amino acid one letter codes.
#If given a filename, consensus sequence is saved in FASTA format to filename.
def consensus(FREQS, filename=None):
    CONSENSUS_SEQ = np.zeros([1, len(FREQS[0,:])],dtype=object) #make an array to store consensus sequence
    for index in range(len(FREQS[0,:])): #for each AA position
        CONSENSUS_SEQ[0,index] = IDS[np.argmax(FREQS[:20,index]),0] #find the largest value, and get the corrisponding AA from IDS, and add it to CONSENSUS_SEQ
    CONSENSUS=""
    for index in range(len(CONSENSUS_SEQ[0,:])):
        CONSENSUS=CONSENSUS+str(CONSENSUS_SEQ[0,index])
    CONSENSUS=">consensus_sequence",CONSENSUS # add header for FASTA format
    if filename:
        np.savetxt((filename),CONSENSUS,delimiter="",fmt="%s") #save file with AA sequence of consensus sequence
    return CONSENSUS_SEQ

#Returns a list of suggested amiono acid mutations when given a query sequence, 
#frequency array, and ratio. Will suggest mutations to consensus when query amino 
#acids that differ from the consensus (i.e. highest frequency) by at least the ratio.
#If given a filename, suggested mutations will be saved as txt file.
def ratioconsensus(query, FREQS, ratio):    
    mutations=[]
    aalist = IDS.flatten().tolist()
    for index in range(len(FREQS[0,:])): #for each AA position
        wtaa = query[index]
        consensus = IDS[np.argmax(FREQS[:20,index]),0]
        if wtaa != consensus: #check if the consens residue is different than the query sequence
            wtfreq = float(FREQS[(aalist.index(wtaa)),index])
            consensusfreq = float(FREQS[(aalist.index(consensus)),index])
            if (ratio * wtfreq) < consensusfreq: #if the consensus of a residue is greater than the threshold
                print "Residue number " + str(int(index) + 1)
                print str(int(100*consensusfreq)) + '% is at least ' + str(ratio) + ' times greater than ' + str(int(100*wtfreq)) + '%'
                thissuggestion=res_to_change(wtaa,(index+1), consensus, consensusfreq, wtfreq)
                mutations.append(thissuggestion)
    return mutations

#Returns a list of suggested amiono acid mutations when given a query sequence, 
#frequency array, and cutoff for the consensus threshold. Will suggest mutations to consensus when query amino 
#acids that differ from the consensus and the consensus is at least cutoff.
#If given a filename, suggested mutations will be saved as txt file.
def cutoffconsensus(query, FREQS, cutoff):    
    mutations=[]
    aalist = IDS.flatten().tolist()
    for index in range(len(FREQS[0,:])): #for each AA position
        wtaa = query[index]
        consensus = IDS[np.argmax(FREQS[:20,index]),0]
        if wtaa != consensus: #cehck if the consens residue is different than the query sequence
            wtfreq = float(FREQS[(aalist.index(wtaa)),index])
            consensusfreq = float(FREQS[(aalist.index(consensus)),index])
            if float(cutoff) < consensusfreq: #if the consensus of a residue is greater than the threshold
                print "Residue number " + str(int(index) + 1)
                print str(int(100*max(FREQS[:20,index]))) + "% is greater than or equal to " + str(int(100*float(cutoff))) + "%"
                thissuggestion=res_to_change(wtaa, (index+1), consensus, consensusfreq, wtfreq)
                mutations.append(thissuggestion)
    return mutations
    
#Takes a list of suggested mutations as 'res_to_change' objects, sorts by % conserved, removes duplicates.
#Returns modified mutations array in with TYPE data types, and human readable suggested mutations list.
#If given filename, will save human readable suggested mutation list as text file.
def formatmutations(mutations_with_dups):
    seen_mutations=set()
    mutations=[]
    for obj in mutations_with_dups:
        if obj.res not in seen_mutations:
            mutations.append(obj)
            seen_mutations.add(obj.res)
    mutations.sort(key=lambda x: x.freq, reverse=True)
    SUGGESTED_MUTATIONS=[]
    SUGGESTED_MUTATIONS.append("These mutations may stabilize your protein since they differ from the consensus residue")
    if not len(mutations):
            SUGGESTED_MUTATIONS.append("No mutations found. Try reducing the ConsensusRatio or ConsensusThreshold. You could also try changing the BLAST parameters to adjust the number of sequences being returned (MaximumSequences and BlastEValue).")
    else:
        for i in mutations: #for each suggested mutation
            SUGGESTED_MUTATIONS.append("Change " + i.wt + " " + str(i.res) + " to " + i.sug + " (" + str(int(round(100*i.freq))) + "% of similar proteins have " + i.sug +", only " + str(round(int(100*i.wtfreq))) + "% have "+ i.wt + ")" ) #add new suggestion on to any existing "SUGGESTED_MUTATIONS"
    return mutations, SUGGESTED_MUTATIONS

#define mutation list based on settings attributes of RATIO and/or CONSESUSTHRESHOLD and using trimmed alignment of sequences to identify query sequence (first sequence in alignment), and array of amino acid frequencies matching amino acid positions.
def mutations(settings, alignment, freqs):
    print('\nIdentifying suggested mutations')
    if settings.RATIO:
        ratiomutations = ratioconsensus(alignment[0], freqs, settings.RATIO)
    if settings.CONSENSUSTHRESHOLD:
        thresholdmutations = cutoffconsensus(alignment[0], freqs, settings.CONSENSUSTHRESHOLD)
        mutations = thresholdmutations
        if settings.RATIO:
            mutations = np.append(mutations, ratiomutations)
    else:
    	mutations = ratiomutations
    mutations, output = formatmutations(mutations)
    return mutations, output
    
#Save output suggestions file with any warnings that have been added
def saveoutput(settings, warnings, output, filename):
    file = open(filename,'wb')
    file.write(''.join(warnings) + '\n')# 's16\n')
    np.savetxt(file, output, delimiter=",", fmt='%s') 
