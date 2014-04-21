#!/usr/bin/python
from collections import Counter
import sys

# this program processes the raw output from chain-helix.txt and used as
# ipput to extrPatDiff.py and extrModiPat.bash
# It maps the SSEs to helix Number

def process(Hs,Hpp):
	HsF = open(Hs,'r').readlines()
	pdbs = []
	PdbHlx = {}
	for i in HsF:
		tmp = i.split()
		pdbs.append(tmp[0])
		if tmp[0] in PdbHlx.keys():
			PdbHlx[tmp[0]].append(tmp[1])
		else:
			PdbHlx[tmp[0]] = []
			PdbHlx[tmp[0]].append(tmp[1])

	pdbsHlx = Counter(pdbs)
	HpF = open(Hpp,'r').readlines()
	
	for j in HpF:
		tmp = j.split()
		for k in PdbHlx.keys():
			if k == tmp[0]:
				print "%d\t%s\t%d\t%d\t%f"%(len(PdbHlx[k]),tmp[0],PdbHlx[k].index(tmp[1]),PdbHlx[k].index(tmp[2]),float(tmp[4]))
				#print len(PdbHlx[k]),tmp[0],PdbHlx[k].index(tmp[1]),PdbHlx[k].index(tmp[2]),tmp[4]
				break
			else:
				continue
		
	return(pdbsHlx)

Hs = sys.argv[1]
Hpp = sys.argv[2]
hlx = process(Hs,Hpp)
Chlx = Counter(hlx.values())
print Chlx,sum(Chlx.values())
