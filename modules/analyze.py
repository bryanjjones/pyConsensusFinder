#!/usr/bin/python
 
import sys #need sys to use system variables
import numpy as np # need numpy for arrays and the like
from Bio import SeqIO, Seq
 
IDS = np.zeros([22,1],dtype=object)
IDS[:,0]=['G', 'P', 'A', 'V', 'L', 'I', 'M', 'C', 'F', 'Y', 'W', 'H', 'K', 'R', 'Q', 'N', 'E', 'D', 'S', 'T', '-', 'other']   

#REDUCED_MAX_SEQS = os.environ['REDUCED_MAX_SEQS'] #pull REDUCED_MAX_SEQS to know if the max sequences had to be reduced

def trimmer(ALINGMENT, filename=None):
    
    LENGTH = len(next(SeqIO.parse(ALINGMENT, "fasta"))) #how many positions do we need to look for gaps and how big is the array?
    SEQUENCES = list(SeqIO.parse(ALINGMENT, "fasta"))
    entries = len(SEQUENCES)
    AA = np.zeros((entries,LENGTH), dtype=np.str) #make array from file with each AA as an element
    for i in range(entries):
        x=str(SEQUENCES[i].seq)
        for position in range(len(x)):
            AA[i,position] = (x[position])
    for index in range(len(AA[0,:])) : #make sure everything is upper case
        AA[0,index] = AA[0,index].upper()
    print (''.join(AA[0,:]))
    print
    for index in range(LENGTH):#at each position
        if AA[0,(LENGTH-1-index)] == "-":#LENGTH-1-index will start at the end, -1 to account for 0 based indexing, and find gaps
            AA = np.delete(AA, (LENGTH-1-index), 1)#delete the gaps
            print(str((''.join(AA[0,:])))) #print the target sequence to show progress
    for index in range(len(AA[:,0])): #alternate ">[GI number]" and sequences
        SEQUENCES[index].seq = Seq.Seq(''.join(AA[(index),:]))
    if filename:
        SeqIO.write(SEQUENCES, filename, "fasta")
    return AA

def aacounts(AA, filename=None): 
    COUNTS = np.zeros([22, len(AA[0,:])],dtype=object) #makes array the length of the alingment with 22 rows (20AAs + "-" + "other")
    for index in range(len(AA[0,:])): # for each position along the alingment, count occourances of each AA
        COUNTS[0,index]=AA[:,index].tolist().count("G")
        COUNTS[1,index]=AA[:,index].tolist().count("P")
        COUNTS[2,index]=AA[:,index].tolist().count("A")
        COUNTS[3,index]=AA[:,index].tolist().count("V")
        COUNTS[4,index]=AA[:,index].tolist().count("L")
        COUNTS[5,index]=AA[:,index].tolist().count("I")
        COUNTS[6,index]=AA[:,index].tolist().count("M")
        COUNTS[7,index]=AA[:,index].tolist().count("C")
        COUNTS[8,index]=AA[:,index].tolist().count("F")
        COUNTS[9,index]=AA[:,index].tolist().count("Y")
        COUNTS[10,index]=AA[:,index].tolist().count("W")
        COUNTS[11,index]=AA[:,index].tolist().count("H")
        COUNTS[12,index]=AA[:,index].tolist().count("K")
        COUNTS[13,index]=AA[:,index].tolist().count("R")
        COUNTS[14,index]=AA[:,index].tolist().count("Q")
        COUNTS[15,index]=AA[:,index].tolist().count("N")
        COUNTS[16,index]=AA[:,index].tolist().count("E")
        COUNTS[17,index]=AA[:,index].tolist().count("D")
        COUNTS[18,index]=AA[:,index].tolist().count("S")
        COUNTS[19,index]=AA[:,index].tolist().count("T")
        COUNTS[20,index]=AA[:,index].tolist().count("-") #empty spaces
        COUNTS[21,index]=(len(AA[:,index]) - sum(COUNTS[:,index].tolist())) #other, not counted above
    IDCOUNTS = np.hstack((IDS,COUNTS)) #make list with AA counts and names of AAs
    if filename:
        np.savetxt((filename),IDCOUNTS,delimiter=",",fmt="%s") #save file with AA names and counts
    return COUNTS
def aafrequencies(COUNTS, filename=None):
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
    
def ratioconsensus(FREQS, THRESHOLD, filename=None):    
    TYPES = [('wt', 'S1'), ('res', int), ('sug', 'S1'),('freq', float)] #make a list of data types for suggested mutation tuples
    MUTATIONS_ARRAY=np.empty([0,])
    MUTATIONS_ARRAY=np.array(MUTATIONS_ARRAY, dtype=TYPES)
    for index in range(len(FREQS[0,:])): #for each AA position
        if IDS[np.argmax(FREQS[:20,index]),0] != AA[0,index]: #and the consens residue is different than the first sequence
            if float(THRESHOLD) < max(FREQS[:20,index]): #if the consensus of a residue is greater than the threshold
                print "Residue number " + str(int(index) + 1)
                print str(int(100*max(FREQS[:20,index]))) + "% is greater than or equal to " + str(int(100*float(THRESHOLD))) + "%"
                SUGGESTION=np.array([(AA[0,index], int(index + 1), CONSENSUS_SEQ[0,index], (-1*max(FREQS[:20,index])))], dtype=TYPES) #entry with negative frequency to allow easy sorting.
                MUTATIONS_ARRAY = np.append(MUTATIONS_ARRAY,SUGGESTION, axis=0)#add new suggestion on to any existing "MUTATIONS_ARRAY"
    MUTATIONS_ARRAY=np.sort(MUTATIONS_ARRAY, order='freq')       
    SUGGESTED_MUTATIONS=np.zeros([1],dtype=object)

    SUGGESTED_MUTATIONS[0]="These mutations may stabilize your protein since they differ from a consensus of over " + str(100*float(THRESHOLD)) + "%"
#    if REDUCED_MAX_SEQS == '1':
#        SUGGESTED_MUTATIONS = np.vstack((SUGGESTED_MUTATIONS,("Results based upon up to 200 BLAST results. BLAST timed out when trying to retrieve the requested number of hits (default request is 2000)."))) #add note to let user know that results are based on fewer sequences
    for index in range(len(MUTATIONS_ARRAY[:,])): #for each suggested mutation
        SUGGESTED_MUTATIONS = np.vstack((SUGGESTED_MUTATIONS,("Change " + MUTATIONS_ARRAY[index,]['wt'] + " " + str(MUTATIONS_ARRAY[index,]['res']) + " to " + MUTATIONS_ARRAY[index,]['sug'] + " (found in " + str(int(-100*MUTATIONS_ARRAY[index,]['freq'])) + "% of similar proteins)" ))) #add new suggestion on to any existing "SUGGESTED_MUTATIONS"
    print SUGGESTED_MUTATIONS

    
    np.savetxt(('./completed/'+NAME+'_mutations.txt'),SUGGESTED_MUTATIONS,delimiter=",",fmt="%s") #save file with suggested stabilizing mutations

def ratioconsensus(FREQS):
    CONSENSUS_SEQ = np.zeros([1, len(FREQS[0,:])],dtype=object) #make an array to store consensus sequence
    TYPES = [('wt', 'S1'), ('res', int), ('sug', 'S1'),('freq', float)] #make a list of data types for suggested mutation tuples
    MUTATIONS_ARRAY=np.empty([0,])
    MUTATIONS_ARRAY=np.array(MUTATIONS_ARRAY, dtype=TYPES)
'''
    for index in range(len(FREQS[0,:])): #for each AA position
        CONSENSUS_SEQ[0,index] = IDS[np.argmax(COUNTS[:20,index]),0] #find the largest value, and get the corrisponding AA from IDS, and add it to CONSENSUS_SEQ
        if float(THRESHOLD) < max(FREQS[:20,index]): #if the consensus of a residue is greater than the threshold
            if IDS[np.argmax(FREQS[:20,index]),0] != AA[0,index]: #and the consens residue is different than the first sequence
                print "Residue number " + str(int(index) + 1)
                print str(int(100*max(FREQS[:20,index]))) + "% is greater than or equal to " + str(int(100*float(THRESHOLD))) + "%"
                SUGGESTION=np.array([(AA[0,index], int(index + 1), CONSENSUS_SEQ[0,index], (-1*max(FREQS[:20,index])))], dtype=TYPES) #entry with negative frequency to allow easy sorting.
                MUTATIONS_ARRAY = np.append(MUTATIONS_ARRAY,SUGGESTION, axis=0)#add new suggestion on to any existing "MUTATIONS_ARRAY"
    MUTATIONS_ARRAY=np.sort(MUTATIONS_ARRAY, order='freq')       
    SUGGESTED_MUTATIONS=np.zeros([1],dtype=object)

    SUGGESTED_MUTATIONS[0]="These mutations may stabilize your protein since they differ from a consensus of over " + str(100*float(THRESHOLD)) + "%"
#    if REDUCED_MAX_SEQS == '1':
#        SUGGESTED_MUTATIONS = np.vstack((SUGGESTED_MUTATIONS,("Results based upon up to 200 BLAST results. BLAST timed out when trying to retrieve the requested number of hits (default request is 2000)."))) #add note to let user know that results are based on fewer sequences
    for index in range(len(MUTATIONS_ARRAY[:,])): #for each suggested mutation
        SUGGESTED_MUTATIONS = np.vstack((SUGGESTED_MUTATIONS,("Change " + MUTATIONS_ARRAY[index,]['wt'] + " " + str(MUTATIONS_ARRAY[index,]['res']) + " to " + MUTATIONS_ARRAY[index,]['sug'] + " (found in " + str(int(-100*MUTATIONS_ARRAY[index,]['freq'])) + "% of similar proteins)" ))) #add new suggestion on to any existing "SUGGESTED_MUTATIONS"
    print SUGGESTED_MUTATIONS

    CONSENSUS=""
    for index in range(len(CONSENSUS_SEQ[0,:])):
        CONSENSUS=CONSENSUS+str(CONSENSUS_SEQ[0,index])
    CONSENSUS=">consensus_sequence",CONSENSUS # add header for FASTA format
    np.savetxt(('./completed/'+NAME+'_consensus.fst'),CONSENSUS,delimiter="",fmt="%s") #save file with AA sequence of consensus sequence
    np.savetxt(('./completed/'+NAME+'_mutations.txt'),SUGGESTED_MUTATIONS,delimiter=",",fmt="%s") #save file with suggested stabilizing mutations

    print 'Completed Trimming'
'''