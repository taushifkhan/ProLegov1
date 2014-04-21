#!/usr/bin/python
import numpy as Np

class contact:
	def __init__(self):
		self.contString = ''
		self.hlxCount = 0
		self.mat = Np.array([])
		self.subPatterns = {}
	
	def feedMatrix(self,patternType):
		'''
		Build matrix from contact pattern. Build the object of contact and call feedMatrix 
		with contact pattern as argument.
		'''
		self.contString = patternType	
		tmp = patternType.split("-")	
		self.hlxCount = len(tmp)+1
#		mat = Np.zeros(shape = (len(tmp)+1,len(tmp)+1))	# Build matrix of helix * helix dimension
		mat = Np.chararray(shape = (len(tmp)+1,len(tmp)+1))
		# Enter Contact detail for matrix elements#######
		for d in range(len(tmp)):
			k = 0
			for m in range(Np.shape(mat)[0]-1):
				dist = d+1
				if m+dist <= (Np.shape(mat)[0]-1):
					mat[m,m+dist] = tmp[d][k]
					k += 1
				else:
					break
		self.mat = mat

	def getHlxCombPatt(self):
		cmbHlx = {}	# dictionary to store all possibilities of helix combination
		for h in range(3,self.hlxCount):
			cmbHlx[h] = []
			self.subPatterns[h] = []
			initialV = 0
			for i in range(self.hlxCount):
				lastVal = initialV+h
				if lastVal <= self.hlxCount:
					cmbHlx[h].append(range(initialV,lastVal)) # Possible sequencial helix combinations
					initialV += 1
				else:
					break
		
		subPatterns = {}
		# Build patterns of helix combinations
		for k1 in cmbHlx.values():
			for k in k1:
				contStrng = []
				for step in range(1,len(k)):
					i = 0
					if len(contStrng):
						contStrng.append("-")

					for j in range(i+step,len(k)):
						contStrng.append(str(self.mat[k[i],k[j]]))
						i += 1
				stng = ''.join(contStrng)
				
				self.subPatterns[len(k)].append(stng)	# append contact string to corresponding HGs
		print self.subPatterns
			
#------------------------------------------------------------------------------------------------------------####		

q = contact()
q.feedMatrix('0000-00r-00-a')
print q.getHlxCombPatt()
