#!/usr/bin/python
import sys,os
import numpy as N
from collections import Counter
import re

# Input from .bash file which gives Interactions for a particular number of 
# Alpha Helix

# Function to populate dictionary that will contain PDB ids as Keys
def dictPDB_populate(dictPDB,tmp,helixgeom):
	print helixgeom
	inf = []
	inf.append(int(tmp[1]))      # Helix1
	inf.append(int(tmp[2]))      # Helix2
	inf.append(float(tmp[3]))      # Global Angle Helix1 & Helix2
	comb_intH = helixgeom[int(tmp[1])]+helixgeom[int(tmp[2])] ### Geometry Of Interacting helices.
	inf.append(comb_intH)
	dictPDB[tmp[0]].append(inf)    

# Function to extract Anti-parallel,Parallel and Orthogonal Orientations
# of Interacting Alpha Helices.Returns a Matrix hlx*hlx Conatinig Interaction
# patterns "iHelix".
def pattern(pdbid,detail):
	iHelix=N.zeros(shape=(hlx,hlx))
	PDBintHgeom[pdbid] = []	# list to save interacting helix geometry corresponding to PDBid
	for i in detail:
#	print i
		if abs(i[2]) >= float(135):
			iHelix[i[0],i[1]] = '1'	# Anti-Parallel
		elif abs(i[2]) <= float(45):
			iHelix[i[0],i[1]] = '2'	# Parallel
		else:
			iHelix[i[0],i[1]] = '3'	# Orthogonal
		PDBintHgeom[pdbid].append(i[3])
		
	return(iHelix)

def modiString(iHelix):
    lInt= []
    sepInt = []
    for step in range(1,N.shape(iHelix)[0]):
        i=0
        for j in range(i+step,N.shape(iHelix)[0]):
            if iHelix[i,j] == 1:
                lInt.append('a')
                sepInt.append('a')
            elif iHelix[i,j] == 2:
                lInt.append('p')
                sepInt.append('p')
            elif iHelix[i,j] == 3:
                lInt.append('r')
                sepInt.append('r')
            else:
                lInt.append('0')
                sepInt.append('0')
            i = i+1
        sepInt.append('-')
    return(lInt,sepInt)
# Interaction List By extracting upper-diagonal elements of the Interaction
# matrix iHelix. Converting -1 for Anti-parallel, +1 for parallel and 1 For
# Orthogonal. Also assigns type of Orientation 'A','P','M'. (Anti-paralle,Parallel & Mixed)

def changeto1(binary_pattern):
    for k in binary_pattern:
        if re.match(r'[0,-]',k):
            continue
        else:
            binary_pattern = binary_pattern.replace(k,'1')
    return binary_pattern
	
### Main Function Starts ##########################################
IntFile = open(sys.argv[1],'r').readlines()	# interaction File
hlx = int(sys.argv[2])	# Number of helix in 

dictPDB = {}

#Wr = open('Pat.txt','a')
Wr =open('IntPatterns_ContactNear_added.txt','a')
# Read and Extract Helix Shapes 
hlxShpFile = open('helSpaceTab','r').readlines()	# Tab delimited helix shape file
anaPDB = {}
PDBintHgeom = {} # 
for i in hlxShpFile:
	htmp = i.split()
	if htmp[0] in anaPDB.keys():
		anaPDB[htmp[0]].append(htmp[2])
	else:
		anaPDB[htmp[0]] = []
		anaPDB[htmp[0]].append(htmp[2])

#print anaPDB

for i in IntFile:
	i = i.rstrip()
	tmp = i.split()
#	0IntMatrix = N.zeros(shape=(totalInt,totalInt))
	print tmp,anaPDB[tmp[0]]
	if tmp[0] in dictPDB.keys():
		dictPDB_populate(dictPDB,tmp,anaPDB[tmp[0]])
	else:
		dictPDB[tmp[0]] = []
		dictPDB_populate(dictPDB,tmp,anaPDB[tmp[0]])

print len(dictPDB.keys())

c= 0
for k in dictPDB.keys():
    iH=pattern(k,dictPDB[k])
    #import ipdb; ipdb.set_trace() # BREAKPOINT
    print iH,k,dictPDB[k],PDBintHgeom[k]
    pat,sepPat = modiString(iH)
    sepPat = ''.join(sepPat)
    sepPat = sepPat.rstrip("-")
    #import ipdb; ipdb.set_trace() # BREAKPOINT
    #Wr.write("%s\t%d\t%d\t%d\t%s\t%s\n"%(k,hlx,len(iString),len(dictPDB[k]),iString,T))
    print k,hlx,len(dictPDB[k]),''.join(pat),sepPat
    pattern_binary = changeto1(sepPat)
    Wr.write("%s\t%d\t%d\t%s\t%s\n"%(k,hlx,len(dictPDB[k]),pattern_binary,sepPat))

Wr.close()
