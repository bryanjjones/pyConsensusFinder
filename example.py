#!/usr/bin/python
import _mypath
import CF

de faults=
settingsfile=

settings=CF.setsettings(defaults,settingsfile)
MainProgram=CF.CF(settings,write=False)
Mutations=MainProgram.mutations

MutationsNotinActiveSite=[]

for i in Mutations:
    print(i.res) #residue number
    print(i.wt) #wildtype residue id
    print(i.sug) #suggested mutation id
    if i (test if i is in active site):
        MutationsNotinActiveSite.append(i)
settings.PDB