import os
import sys
import csv
import numpy as np
import analyze
import CF
from Bio.PDB import *

class PDB(object):
	def __init__(self,defaults=None,configfile=None,settings=None,write=True,writeextended=False):
		#write will rewrite suggested mutations file (XXX_mutations.txt) writeextende will write extra files: detailed_output.csv, residue_list.csv, and simple_output.csv
		if settings is None:
			settings=CF.setsettings(defaults,configfile)
		##Make new directories
		filebase=os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
		outfilebase=filebase+'/completed'
		if not os.path.exists(outfilebase):
		    os.makedirs(outfilebase)

		##Defining lists
		sort1=[]
		sort2=[]
		sort3=[]
		sort4=[]
		sort5=[]
		info=[]
		pinfo=[]
		hitnum_list=[]
		hitting=[]
		atomlist=[]
		atomnamelist=[]
		resname=[]
		rescode=[]
		heteroname=[]
		heterocode=[]
		warnings = []

		##reading in config values
		try:
			if settings.RESIDUE <= 2:
				print('please enter a distance greater than 2 angstroms')
				sys.exit()
			else:
				pdb_name=settings.PDB.lower()
				achn=settings.CHAIN
				ares=settings.RESIDUE
				m=settings.ANG
				try:
					with open(filebase+'/completed/'+settings.FILENAME+'_mutations.txt', 'r') as consensus_file:
						for line in consensus_file:
							info.append(line.split())
				except:
					warnings.append('no output from consensus finder, cannot trim anything')
					print('no output from consensus finder, cannot trim anything')
		except:
			print('there is an error with the pdb or residue name provided')
			sys.exit()

		##creating new directories and initializing biopdb 
		plist = PDBList()
		plist.retrieve_pdb_file(pdb_name,pdir=filebase+'/uploads/pdb_data',file_format='pdb')
		pparse = PDBParser(QUIET=True)
		structure = pparse.get_structure('X', filebase+'/uploads/pdb_data/pdb'+pdb_name+'.ent')

		##generating residue list 
		rf = open(os.path.join(outfilebase, "residue_list"+'.csv'), 'wt')
		reswriter = csv.writer(rf, lineterminator='\n')
		reswriter.writerow(["-------------------------------------------------------------------"])
		model = structure[0]
		chain = model[achn]
		for res in chain.get_residues():
		    tags = res.get_full_id()
		    if res.get_resname() != 'HOH' and tags[3][0] == " ":
		        resname.append(res.get_resname())
		        resid = res.get_id()
		        rescode.append(resid[1])
		reswriter.writerow(["A total of"+" "+str(len(resname))+" "+"residues in chain "+achn])
		for i in range(0, len(resname), 10):
		    reswriter.writerow(["["+str(rescode[i])+"]"+" "+" ".join([str(v) for v in resname[i:i+10]])])
		reswriter.writerow(["-------------------------------------------------------------------"])
		reswriter.writerow(["List of Heteroatoms: "])
		for res in chain.get_residues():
		    tags = res.get_full_id()
		    if res.get_resname() != 'HOH' and tags[3][0] != " ":
		        heteroname.append(res.get_resname())
		        resid = res.get_id()
		        heterocode.append(resid[1])
		for i in range(0, len(heteroname), 10):
		        reswriter.writerow(["["+str(heterocode[i])+"]"+" "+ " ".join([str(v) for v in heteroname[i:i+10]])])
		if not rescode[0] <= ares <= rescode[len(resname)-1]:
			CF.cleanexit("Residue number not in range. Enter residue between "+str(rescode[0])+" and "+str(rescode[len(resname)-1]))
		##creating simple and detailed outputs
		sf = open(os.path.join(outfilebase, "simple_output" + '.csv'), 'wt')
		df = open(os.path.join(outfilebase, "detailed_output" + '.csv'), 'wt')
		try:
			writer = csv.writer(sf, lineterminator='\n')
			writer2 = csv.writer(df, lineterminator='\n')
			# Search with active site residue
			residue1 = chain[ares]
			res1id = residue1.get_id()
		
######## print the residue get id in case of error 
			for atom in residue1:
			    atomnamelist.append(atom.get_name())
			    atomlist.append(atom.get_vector())
			for i in range(0, len(atomlist)):
			    for model in structure:
			        for chain in model:
			            for residue in chain:
			                for atom in residue:
			                    dist = np.linalg.norm(atom.get_vector() - atomlist[i])
			                    if dist <= m:
			                        resid1 = residue.get_id()
			                        rdist = round(dist, 2)
			                        if residue.get_resname() != 'HOH':
			                            if resid1[1] != ares:
			                                sort1.append((residue.get_resname(), resid1[1], atom.get_name(), rdist,
			                                              residue1.get_resname(), res1id[1], atomnamelist[i]))
			                                sort1 = sorted(sort1, key=lambda e: (e[3]))
			writer2.writerow(['Residue Number: ' + str(ares)])
			writer2.writerow(['Residue', 'Distance', 'Origin'])
			for i in range(0, len(sort1)):
			    writer2.writerow([sort1[i][0] + " " + str(sort1[i][1]) + " " + sort1[i][2], sort1[i][3],
			                      sort1[i][4] + " " + str(sort1[i][5]) + " " + sort1[i][6]])
			    if sort1[i][1] not in sort2:
			        sort2.append(sort1[i][1])
			        sort3.append(sort1[i])
			    else:
			        sort4.append(sort1[i])
			        sort4 = sorted(sort4, key=lambda e: (e[1]))
			writer.writerow(['Residue', ' ', ' ', 'Distance', 'Origin', ' ', ' ', 'Other Atoms'])
			for i in range(0, len(sort3)):
			    num = sort3[i][1]
			    for r in range(0, len(sort4)):
			        if sort4[r][1] == num:
			            if sort4[r][2] not in sort5:
			                if sort4[r][2] != sort3[i][2]:
			                    sort5.append(sort4[r][2])
			    writer.writerow(
			        [sort3[i][0], str(sort3[i][1]), sort3[i][2], sort3[i][3], residue1.get_resname(), ares, sort3[i][6],
			         ','.join(sort5)])
			    del sort5[:]
		finally:
		    sf.close()
		    df.close()

		##trimming outputs
		with open(outfilebase+'//'+'simple_output.csv', 'r') as pdbresults:
		    for line in pdbresults:
		        pinfo.append(line.split(','))
		#Generating search hit
		for i in range(0, len(pinfo)):
		    hitnum_list.append(pinfo[i][1])
		hitnum_list.append(str(ares))

		#Writing output file
		if write:
			f = open(filebase+'/completed/'+settings.FILENAME+'_mutations.txt','w')
			for r in range(0, len(info)):
			    if len(info[r]) is 15:
			    	if str(info[r][2]) not in hitnum_list:
			        	f.write(' '.join(info[r])+'\n')
			    else:
			    	if len(info[r]) is 15:
			        	hitting.append(info[r][2])
			hitting.append(str(ares))
			if warnings:
				f.write(warnings[0])
			else:
				f.write('The following residues are removed because they are part of the active site: ')
				for i in range(0, len(hitting)):
					f.write(hitting[i]+' ')
