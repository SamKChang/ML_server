#!/usr/bin/env python

import pickle
import cPickle
import time
import re
import json
import matplotlib.pyplot as plt
import parse_poscar
import DescriptorHandeler
from Constants import *
from numpy import *
from numpy import linalg, zeros, array, argsort
from scipy import special
from copy import deepcopy
from collections import Counter
from os import listdir,makedirs
from os.path import isfile, join,exists
from CustomOpperators import *
from pymatgen import Lattice, Structure
from pymatgen.io.vaspio import Vasprun, Poscar, Chgcar
from pymatgen.io.smartio import read_structure, write_structure, \
	read_mol, write_mol
from true_ewald_matrix import true_ewald_matrix

class ColumbPNorm(object):
	def __init__(self,size):
		self.size = size
		self.sumlen = 3
		self.p = 6
	
	def descibe(self,cell,coords,ocupationList):


		C = zeros(self.size**2).reshape(self.size,self.size)
		for i in range(0,len(ocupationList)):			  
				for k in range(i-1,len(ocupationList)):
						C[i,k] = self.lPNorm(cell*(coords[i]-coords[k]),cell, ocupationList[i],ocupationList[k])/linalg.det(cell)
						C[k,i] = C[i,k]


		return self.SortMatrix(C)


	def lPNorm(self,x,cell,z1,z2):

		c = 0
		for i in range(-self.sumlen,self.sumlen+1):
			for j in range(-self.sumlen,self.sumlen+1):
				for k in range(-self.sumlen,self.sumlen+1):
					
					dist =  linalg.norm(x + i*cell[0]+j*cell[1]+k*cell[2])
					if dist < 10**-6: 
						c += ((z1**(2.4))/(2.0))**self.p
					else:

						c += (z1*z2/dist)**self.p

		return c**(1.0/self.p)

		
class ExtendedColumbErfNorm(object):
	def __init__(self,size,sumlen,cutoff,debug=False):
		self.debug = debug
		self.size = size
		self.sumlen = sumlen
		self.cutoff = cutoff
		self.boxSize = 2*self.sumlen + 1
		set_printoptions(threshold= 'nan')
	
	def descibe(self,cell,coords,ocupationList):
		
		
		C = zeros(self.boxSize**3*self.size**2).reshape(self.size,self.boxSize**3*self.size)
		if self.debug == True:
			print 'Cell: ', cell
			print 'Coords: ', coords
		
		tmp = zeros(3*len(coords)*self.boxSize**3).reshape(len(coords),self.boxSize**3,3)
		tmp2 = zeros(len(coords)*self.boxSize**3)
		for i in range(len(coords)):
			tmp[i] = self.Extend(coords[i,:])
			for j in range(self.boxSize**3):
				tmp2[(self.boxSize**3)*(i):(self.boxSize**3)*(i+1)]\
				 = [ocupationList[i] for x in range(self.boxSize**3)]
		
		coordsExtended = asarray([val for sublist in tmp for val in sublist])
		ocupationListExtended = asarray(tmp2)
		if self.debug == True:
			print 'Coordlen: ',len(coords)
			print 'occLen: ', len(ocupationList)
		
		for i in range(0,len(ocupationList)):			   
			for k in range(0,len(ocupationListExtended)):
					atomDist = linalg.norm(cell*(coords[i]-coordsExtended[k]),2)
					if atomDist==0: 
						C[i,k] = (ocupationList[i]**2.4) /2
					else:		
						atomDist = linalg.norm(cell*(coords[i]-coordsExtended[k]),2)
						tmpDist = atomDist- self.cutoff
						C[i,k] = math.erfc(tmpDist)*ocupationList[i]*ocupationListExtended[k]/atomDist
		
		#C = sum(C, axis=0)
		if self.debug == True:
			print self.SortMatrix(C)
		return self.SortMatrix(C)


	def Extend(self,x):
		cellRange = range(-self.sumlen,self.sumlen+1)
		y = zeros(3*self.boxSize**3).reshape(self.boxSize**3,3)
		print x
		I = eye(3)
		for i in cellRange:
			for j in cellRange:
				for k in cellRange:
					index =  + (i+ self.sumlen)*(self.boxSize)**2 + (j+ self.sumlen)*(self.boxSize) + k + self.sumlen
					y[index,:] = x + i*I[0,:]+j*I[1,:]+k*I[2,:]
		return y
		
	def SortMatrix(self,X):
		dimX1 = len(X[:,1])
		dimX2 = len(X[1,:])
		NormList1 = zeros(dimX1)
		NormList2 = zeros(dimX2)
		for i in range(0,dimX1):
			NormList1[i] = linalg.norm(X[i,:])
		NormList1 = argsort(NormList1)[::-1]
		X[:,:] = X[NormList1,:]

		for i in range(0,dimX2):
			NormList2[i] = linalg.norm(X[:,i])
		NormList2 = argsort(NormList2)[::-1]
		X[:,:] = X[:,NormList2]		

		return X
		
		
class ExtendedColumbExpNorm(object):
	def __init__(self,size,sumlen,kernelWidth,debug=False):
		self.debug = debug
		self.size = size
		self.sumlen = sumlen
		self.kernelWidth = kernelWidth
		self.boxSize = 2*self.sumlen + 1
		set_printoptions(threshold= 'nan')
	
	def descibe(self,cell,coords,ocupationList):
		
		
		C = zeros(self.boxSize**3*self.size**2).reshape(self.size,self.boxSize**3*self.size)
		if self.debug == True:
			print 'Cell: ', cell
			print 'Coords: ', coords
		
		tmp = zeros(3*len(coords)*self.boxSize**3).reshape(len(coords),self.boxSize**3,3)
		tmp2 = zeros(len(coords)*self.boxSize**3)
		for i in range(len(coords)):
			tmp[i] = self.Extend(coords[i,:])
			for j in range(self.boxSize**3):
				tmp2[(self.boxSize**3)*(i):(self.boxSize**3)*(i+1)]\
				 = [ocupationList[i] for x in range(self.boxSize**3)]
		
		coordsExtended = [val for sublist in tmp for val in sublist]
		ocupationListExtended = tmp2
		if self.debug == True:
			print 'Coordlen: ',len(coords)
			print 'occLen: ', len(ocupationList)
		
		for i in range(0,len(ocupationList)):			   
			for k in range(0,len(ocupationListExtended)):
					C[i,k] = ocupationList[i]*ocupationListExtended[k]*\
					exp(-((linalg.norm(cell*(coords[i]-coordsExtended[k]),2)\
					/self.kernelWidth)**2)/2)
		
		#C = sum(C, axis=0)
		if self.debug == True:
			print self.SortMatrix(C)
		return self.SortMatrix(C)


	def Extend(self,x):
		cellRange = range(-self.sumlen,self.sumlen+1)
		y = zeros(3*self.boxSize**3).reshape(self.boxSize**3,3)
		print x
		I = eye(3)
		for i in cellRange:
			for j in cellRange:
				for k in cellRange:
					index =  + (i+ self.sumlen)*(self.boxSize)**2 + (j+ self.sumlen)*(self.boxSize) + k + self.sumlen
					y[index,:] = x + i*I[0,:]+j*I[1,:]+k*I[2,:]
		return y
		
	def SortMatrix(self,X):
		dimX1 = len(X[:,1])
		dimX2 = len(X[1,:])
		NormList1 = zeros(dimX1)
		NormList2 = zeros(dimX2)
		for i in range(0,dimX1):
			NormList1[i] = linalg.norm(X[i,:])
		NormList1 = argsort(NormList1)[::-1]
		X[:,:] = X[NormList1,:]

		for i in range(0,dimX2):
			NormList2[i] = linalg.norm(X[:,i])
		NormList2 = argsort(NormList2)[::-1]
		X[:,:] = X[:,NormList2]		

		return X

class RadialDistribution(object):
	def __init__(self,size,sumlen,kernelWidth,debug=False):
		self.debug = debug
		self.size = size
		self.sumlen = sumlen
		self.kernelWidth = kernelWidth
		self.boxSize = 2*self.sumlen + 1
		set_printoptions(threshold= 'nan')
	
	def descibe(self,cell,coords,ocupationList):
		
		
		C = zeros(self.boxSize**3*self.size**2).reshape(self.size,self.boxSize**3*self.size)
		if self.debug == True:
			print 'Cell: ', cell
			print 'Coords: ', coords
		
		tmp = zeros(3*len(coords)*self.boxSize**3).reshape(len(coords),self.boxSize**3,3)
		tmp2 = zeros(len(coords)*self.boxSize**3)
		for i in range(len(coords)):
			tmp[i] = self.Extend(coords[i,:])
			for j in range(self.boxSize**3):
				tmp2[(self.boxSize**3)*(i):(self.boxSize**3)*(i+1)]\
				 = [ocupationList[i] for x in range(self.boxSize**3)]
		
		coordsExtended = [val for sublist in tmp for val in sublist]
		ocupationListExtended = tmp2
		if self.debug == True:
			print 'Coordlen: ',len(coords)
			print 'occLen: ', len(ocupationList)
		
		for i in range(0,len(ocupationList)):			   
			for k in range(0,len(ocupationListExtended)):
					C[i,k] = ocupationList[i]*ocupationListExtended[k]*\
					exp(-((linalg.norm(cell*(coords[i]-coordsExtended[k]),2)\
					/self.kernelWidth)**2)/2)
		
		#C = sum(C, axis=0)
		if self.debug == True:
			print self.SortMatrix(C)
		return self.SortMatrix(C)


	def Extend(self,x):
		cellRange = range(-self.sumlen,self.sumlen+1)
		y = zeros(3*self.boxSize**3).reshape(self.boxSize**3,3)
		print x
		I = eye(3)
		for i in cellRange:
			for j in cellRange:
				for k in cellRange:
					index =  + (i+ self.sumlen)*(self.boxSize)**2 + (j+ self.sumlen)*(self.boxSize) + k + self.sumlen
					y[index,:] = x + i*I[0,:]+j*I[1,:]+k*I[2,:]
		return y

		
class Ewald(object):
	def __init__(self,size,eta = None):
		self.size = size
		self.eta  = eta
	
	def descibe(self,structure):
		C = zeros(self.size**2).reshape(self.size,self.size)
		tmp = true_ewald_matrix(structure,eta = self.eta)
		C[: tmp.shape[0],: tmp.shape[1]]  = tmp  
		return self.SortMatrix(C)

	def SortMatrix(self,X):
		dimX = len(X)
		NormList = zeros(dimX)
		for i in range(0,dimX):
			NormList[i] = linalg.norm(X[i])
		NormList = argsort(NormList)[::-1]
		X[:,:] = X[:,NormList]
		X[:,:] = X[NormList,:]
		return X

class EwaldWithSpaceGroup(object):
	def __init__(self,size):
		self.size = size

	
	def descibe(self,structure):
		C = zeros(self.size**2+ self.size).reshape(self.size,self.size + 1)
		tmp = EwaldSummation(structure).total_energy_matrix
		C[: tmp.shape[0],: tmp.shape[1]]  = tmp
		V = sym.SymmetryFinder(structure).get_spacegroup_number()
		C[:,0:-1] = self.SortMatrix(C[:,0:-1])
		C[0,-1] = V
		return	C

	def SortMatrix(self,X):
		dimX = len(X)
		NormList = zeros(dimX)
		for i in range(0,dimX):
			NormList[i] = linalg.norm(X[i])
		ArgList = argsort(NormList)[::-1]
		X[:,:] = X[:,ArgList]
		X[:,:] = X[ArgList,:]
		return X


class NonSpacial(object):
	def __init__(self,chargetype = 'PTP',debug = False):
		self.chargetype = chargetype
		self.debug = debug
		self.tag = 'wy'


		self.effeciveCharge = {5:[4.680,2.576,2.421,0,0,0,0,0,0,0,0], 6:[5.673,3.217,3.136,0,0,0,0,0,0,0,0],\
			7:[6.665,3.847,3.834,0,0,0,0,0,0,0,0], 8:[7.658,4.492,4.453,0,0,0,0,0,0,0,0],\
			9:[8.650,5.128,5.100,0,0,0,0,0,0,0,0], 13:[12.591,8.214,8.963,4.117,4.066,0,0,0,0,0,0],\
			14:[13.575,9.020,9.945,4.903,4.285,0,0,0,0,0,0],15:[14.558,9.825,10.961,5.642,4.886,0,0,0,0,0,0],\
			16:[15.541,10.629,11.977,6.367,5.482,0,0,0,0,0,0], 17:[16.524,11.430,12.993,7.068,6.116,0,0,0,0,0,0],\
			31:[30.309,22.599,27.091,16.996,16.204,7.067,15.093,6.222,0,0,0],\
			32:[31.294,23.365,28.082,17.790,17.014,8.044 ,16.251,6.780,0,0,0],\
			33:[32.278,24.127,29.074,18.596,17.850,8.944,17.378,7.449,0,0,0],\
			34:[33.262,24.888,30.065,19.403,18.705,9.758,18.477,8.287,0,0,0],\
			35:[34.247,25.643,31.056,20.219,19.571,10.553,19.559,9.028,0,0,0],\
			49:[48.010,36.124,44.898,31.631,31.521,21.761,34.678,20.369,9.512,16.942,8.470],\
			50:[48.992,36.859,45.885,32.420,32.353,22.658,35.742,21.265,10.629,17.970,9.102],\
			51:[49.974,37.595,46.873,33.209,33.184,23.544,36.800,22.181,11.617,18.974,9.995]}

		self.effeciveChargeOuter = {5:2.421, 6:3.136, 7:3.834, 8:4.453, 9:5.100, 13:4.066,\
			14:4.285,15:4.886, 16:5.482, 17:6.116, 31:6.222, 32:6.780, 33:7.449, 34:8.287,\
			31:9.028, 49:8.470, 50:9.102, 51:9.995}
		
		self.PeriodicTablePossition = {1:[1,1],
		2:[1,8], 3:[2,1], 4:[2,2], 5:[2,3],\
		6:[2,4], 7:[2,5], 8:[2,6], 9:[2,7],\
		10:[2,8],11:[3,1],12:[3,2],13:[3,3],\
		14:[3,4],15:[3,5],16:[3,6],17:[3,7],\
                18:[3,8],19:[4,1],20:[4,2],31:[4,3],32:[4,4],\
                33:[4,5],34:[4,6],35:[4,7],36:[4,8],37:[5,1],\
                38:[5,2],49:[5,3],50:[5,4],51:[5,5],52:[5,6],\
                53:[5,7],54:[5,8],55:[6,1],56:[6,2],81:[6,3],\
                82:[6,4],83:[6,5]}


		self.PTP = {\
		 1  :[1,1] ,2:  [1,8]#Row1
		
		,3  :[2,1] ,4:  [2,2]#Row2\
		,5  :[2,3] ,6:  [2,4] ,7  :[2,5] ,8  :[2,6] ,9  :[2,7] ,10 :[2,8]\
		
		,11 :[3,1] ,12: [3,2]#Row3\
		,13 :[3,3] ,14: [3,4] ,15 :[3,5] ,16 :[3,6] ,17 :[3,7] ,18 :[3,8]\
		
		,19 :[4,1] ,20: [4,2]#Row4\
		,31 :[4,3] ,32: [4,4] ,33 :[4,5] ,34 :[4,6] ,35 :[4,7] ,36 :[4,8]\
		,21 :[4,9] ,22: [4,10],23 :[4,11],24 :[4,12],25 :[4,13],26 :[4,14],27 :[4,15],28 :[4,16],29 :[4,17],30 :[4,18]\

		,37 :[5,1] ,38: [5,2]#Row5\
		,49 :[5,3] ,50: [5,4] ,51 :[5,5] ,52 :[5,6] ,53 :[5,7] ,54 :[5,8]\
		,39 :[5,9] ,40: [5,10],41 :[5,11],42 :[5,12],43 :[5,13],44 :[5,14],45 :[5,15],46 :[5,16],47 :[5,17],48 :[5,18]\

		,55 :[6,1] ,56: [6,2]#Row6\
		,81 :[6,3] ,82: [6,4] ,83 :[6,5] ,84 :[6,6] ,85 :[6,7] ,86 :[6,8]
			   ,72: [6,10],73 :[6,11],74 :[6,12],75 :[6,13],76 :[6,14],77 :[6,15],78 :[6,16],79 :[6,17],80 :[6,18]\
		,57 :[6,19],58: [6,20],59 :[6,21],60 :[6,22],61 :[6,23],62 :[6,24],63 :[6,25],64 :[6,26],65 :[6,27],66 :[6,28],67 :[6,29],68 :[6,30],69 :[6,31],70 :[6,32],71 :[6,33]\

		,87 :[7,1] ,88: [7,2]#Row7\
		,113:[7,3] ,114:[7,4] ,115:[7,5] ,116:[7,6] ,117:[7,7] ,118:[7,8]\
			   ,104:[7,10],105:[7,11],106:[7,12],107:[7,13],108:[7,14],109:[7,15],110:[7,16],111:[7,17],112:[7,18]\
		,89 :[7,19],90: [7,20],91 :[7,21],92 :[7,22],93 :[7,23],94 :[7,24],95 :[7,25],96 :[7,26],97 :[7,27],98 :[7,28],99 :[7,29],100:[7,30],101:[7,31],101:[7,32],102:[7,14],103:[7,33]}

		
	def descibe(self,name):
		#C = name.split('-')[0:-2]
		if self.debug:
			print 'Unrefined: ', name
		name = name[1::]
		C = re.split('0|1|2|3|4|5|6|7|8|9',name)
		if self.debug:		
			print 'Splited: ', C
		C = filter(None, C)
		if self.debug:	
			print 'Filtered: ',C	
		C = [parse_poscar.periodictable_numbers[c[1::]] for c in C]
		if self.debug:		
			print 'Refined: ', C		
		C = asarray([float(c) for c in C])
		if self.chargetype == 'effective':
			C = asarray([asarray(self.effeciveCharge[c]) for c in C])
		elif self.chargetype == 'outer':
			C = asarray([self.effeciveChargeOuter[c] for c in C])
		elif self.chargetype == 'possition':
			C = asarray([self.PeriodicTablePossition[c] for c in C])
		elif self.chargetype == 'PTP':
			C = asarray([self.PTP[c] for c in C])
		elif self.chargetype == 'normal':
			tmp = zeros(118)
			C = asarray([int(c) for c in C])
			for i in range(len(C)):		
				tmp[C[i]] =  i + 1
			tmp[0] = 1
			C = tmp
		if self.debug:		
			print 'Descriptor: ',C
		return	C

class SymGroupDescriptor(object):
	def __init__(self,wyLen = 10,debug = False):
		self.debug = debug
		self.wyLen = wyLen
		self.tag = 'wy'
	def descibe(self,wySec):
		#C = name.split('-')[0:-2]
		if self.debug:
			print 'Unrefined: ', wySec
		#wySec = wySec[1::]
		dOut = zeros([self.wyLen,118])
		wyStr = re.split('0|1|2|3|4|5|6|7|8|9',wySec)
		wyStr = filter(None, wyStr)
		tmp = ''.join(wyStr)[0] + ''.join(['|'+ i for i in ''.join(wyStr)])
		wyInt = re.split(tmp,wySec) 	
		wyInt = asarray(filter(None, wyInt)).astype(int)
		wyPos = [s[0] for s in wyStr]
		wyEle = [parse_poscar.periodictable_numbers[s[1::]] for s in wyStr]

		if self.debug:		
			print 'Num: ', wyInt
			print 'Pos: ', wyPos
			print 'Ele: ', wyEle
			#time.sleep(0.01)
		wyPos = [ord(s)- 97 for s in wyPos]
		for i in range(len(wyEle)):
			dOut[wyPos[i],wyEle[i]] = wyInt[i]
		if self.debug:		
			print 'descriptor: ', dOut

		return	dOut
class SymGroupDescriptorPTP(object):
	def __init__(self,wyLen = 10,nEle = 10,debug = False):
		self.debug = debug
		self.wyLen = wyLen
		self.nEle = nEle
		self.tag = 'wy'

	def descibe(self,wySec):
		#C = name.split('-')[0:-2]
		if self.debug:
			print 'Unrefined: ', wySec
		#wySec = wySec[1::]
		dOut = zeros([self.wyLen,self.nEle,7,33])
		wyStr = re.split('0|1|2|3|4|5|6|7|8|9',wySec)
		wyStr = filter(None, wyStr)
		tmp = ''.join(wyStr)[0] + ''.join(['|'+ i for i in ''.join(wyStr)])
		wyInt = re.split(tmp,wySec) 	
		wyInt = asarray(filter(None, wyInt)).astype(int)
		wyPos = [s[0] for s in wyStr]
		wyEle = asarray([PTP[parse_poscar.periodictable_numbers[s[1::]]] for s in wyStr])

		if self.debug:		
			print 'Num: ', wyInt
			print 'Pos: ', wyPos
			print 'Ele: ', wyEle
			#time.sleep(0.01)
		wyPos = [ord(s)- 97 for s in wyPos]
		for i in range(len(wyEle)):
			dOut[wyPos[i],wyInt[i],wyEle[i,0]-1,wyEle[i,1]-1] = 1.0
		if self.debug:		
			print 'descriptor: ', dOut

		return	dOut
class SymGroupDescriptorPTP2(object):
	def __init__(self,wyLen = 10,nEle = 10,debug = False):
		self.debug = debug
		self.wyLen = wyLen
		self.nEle = nEle
		self.tag = 'wy'
		
	def descibe(self,wySec):
		#C = name.split('-')[0:-2]
		if self.debug:
			print 'Unrefined: ', wySec
		#wySec = wySec[1::]
		dOut = zeros([self.wyLen,self.nEle,3])
		wyStr = re.split('0|1|2|3|4|5|6|7|8|9',wySec)
		wyStr = filter(None, wyStr)
		tmp = ''.join(wyStr)[0] + ''.join(['|'+ i for i in ''.join(wyStr)])
		wyInt = re.split(tmp,wySec) 	
		wyInt = asarray(filter(None, wyInt)).astype(int)
		wyPos = [s[0] for s in wyStr]
		wyEle = asarray([PTP[parse_poscar.periodictable_numbers[s[1::]]] for s in wyStr])

		if self.debug:		
			print 'Num: ', wyInt
			print 'Pos: ', wyPos
			print 'Ele: ', wyEle
			#time.sleep(0.01)
		wyPos = [ord(s)- 97 for s in wyPos]
		j = 0
		pWyPos = -1
		for i in range(len(wyEle)):
			if wyPos[i] == pWyPos:
				j = j+1
			else:
				j = 0
			dOut[wyPos[i],j] = [wyEle[i,0],wyEle[i,1],wyInt[i]]
			pWyPos = wyPos[i]
		if self.debug:		
			print 'descriptor: ', dOut

		return	dOut
class SymGroupDescriptorPTP3(object):
	def __init__(self,wyLen = 10,nEle = 10,debug = False):
		self.debug = debug
		self.wyLen = wyLen
		self.nEle = nEle
		self.tag = 'wy'
		
	def descibe(self,wySec):
		#C = name.split('-')[0:-2]
		if self.debug:
			print 'Unrefined: ', wySec
		#wySec = wySec[1::]
		dOut = zeros([self.wyLen,self.nEle,2])
		wyStr = re.split('0|1|2|3|4|5|6|7|8|9',wySec)
		wyStr = filter(None, wyStr)
		tmp = ''.join(wyStr)[0] + ''.join(['|'+ i for i in ''.join(wyStr)])
		wyInt = re.split(tmp,wySec) 	
		wyInt = asarray(filter(None, wyInt)).astype(int)
		wyPos = [s[0] for s in wyStr]
		wyEle = asarray([PTP[parse_poscar.periodictable_numbers[s[1::]]] for s in wyStr])

		if self.debug:		
			print 'Num: ', wyInt
			print 'Pos: ', wyPos
			print 'Ele: ', wyEle
			#time.sleep(0.01)
		wyPos = [ord(s)- 97 for s in wyPos]
		j = 0
		pWyPos = -1
		for i in range(len(wyEle)):
			if wyPos[i] == pWyPos:
				j = j+1
			else:
				j = 0
			for k in range(wyInt[i]):
				tmp = j + k				
				dOut[wyPos[i],tmp] = [wyEle[i,0],wyEle[i,1]]
			j = tmp
			pWyPos = wyPos[i]
		if self.debug:		
			print 'descriptor: ', dOut

		return	dOut
	
class SymGroupDescriptorPTP4(object):
	def __init__(self,desLen = 200,debug = False):
		self.debug = debug
		self.desLen = desLen
		self.tag = 'wy'
		inFi = open('./WyData.pkl', 'rb')
		self.wydata = cPickle.load(inFi)
		inFi.close()
		
	def descibe(self,wySec,spcGrp):
		wydata = self.wydata[spcGrp]
		
		if self.debug:
			print 'Unrefined: ', wySec
		#wySec = wySec[1::]
		dOut = zeros([self.desLen,14])
		wyStr = re.split('0|1|2|3|4|5|6|7|8|9',wySec)
		wyStr = filter(None, wyStr)
		tmp = ''.join(wyStr)[0] + ''.join(['|'+ i for i in ''.join(wyStr)])
		wyInt = re.split(tmp,wySec) 	
		wyInt = asarray(filter(None, wyInt)).astype(int)
		wyLet = [s[0] for s in wyStr]
		wyCoo = [wydata[s] for s in wyLet]
		wyEle = asarray([PTP[parse_poscar.periodictable_numbers[s[1::]]] for s in wyStr])
		tmpWyCoo = []
		tmpwyEle = []
		for i in range(len(wyInt)):
			for j in range(wyInt[i]):
				tmpWyCoo.append(wyCoo[i])
				tmpwyEle.append(wyEle[i])
		wyCoo = tmpWyCoo
		wyEle = tmpwyEle
		wyDes = []
		for i in range(len(wyEle)):
			for j in range(len(wyCoo[i])):
				wyDes.append(append(wyCoo[i][j].flatten(),wyEle[i]).flatten()) 
		wyDes = asarray(wyDes)	
		aSort = argsort([linalg.norm(d) for d in wyDes])
		wyDes = wyDes[aSort]
		
		dOut[0:len(wyDes)] = wyDes
		
		if self.debug:		
			print 'Num: ', wyInt
			#print 'Coo: ', wyCoo
			print 'Let: ', wyLet
			print 'Ele: ', wyEle
			print 'Des: ', wyDes
			#time.sleep(0.01)
		
		
		
		#if self.debug:		
		#	print 'descriptor: ', dOut
		return	dOut
	
	
	
class ElpasoDescriptor(object):
	def __init__(self,debug = False):
		self.debug = debug
		self.ElpasoOrder = [{9:1,8:2,2:3,3:4,17:5,16:6,34:7,1:8,20:9,56:10,35:11,38:12,52:13,31:14,13:15,15:16,12:17,53:18,50:19,49:20,
			33:21,51:22,11:23,14:24,32:25,81:26,83:27,82:28,4:29,19:30,55:31,37:32,10:33,5:34,7:35,18:36,36:37,54:38,6:39},
			{2:1,9:2,8:3,15:4,16:5,14:6,1:7,33:8,3:9,13:10,4:11,34:12,32:13,31:14,5:15,51:16,52:17,12:18,50:19,7:20,
			17:21,49:22,20:23,83:24,35:25,38:26,82:27,6:28,81:29,11:30,53:31,56:32,19:33,37:34,10:35,55:36,18:37,36:38,54:39},
			{56:1,38:2,20:3,37:4,19:5,55:6,11:7,81:8,83:9,49:10,82:11,3:12,12:13,52:14,50:15,51:16,53:17,35:18,31:19,17:20
			,10:21,18:22,36:23,34:24,54:25,13:26,33:27,32:28,9:29,16:30,14:31,15:32,4:33,1:34,8:35,2:36,5:37,7:38,6:39},
 			{9:1,17:2,35:3,8:4,56:5,38:6,20:7,3:8,55:9,53:10,11:11,19:12,37:13,81:14,12:15,49:16,31:17,34:18,82:19,52:20,
			2:21,50:22,1:23,16:24,83:25,10:26,51:27,13:28,54:29,36:30,18:31,32:32,33:33,15:34,14:35,4:36,7:37,5:38,6:39}]
		
	def descibe(self,name):
		#C = name.split('-')[0:-2]
		if self.debug:
			print 'Unrefined: ', name
		name = name[1::]
		C = re.split('0|1|2|3|4|5|6|7|8|9',name)
		if self.debug:		
			print 'Splited: ', C
		C = filter(None, C)
		if self.debug:	
			print 'Filtered: ',C	
		C = [parse_poscar.periodictable_numbers[c[1::]] for c in C]
		if self.debug:		
			print 'Refined: ', C		
		C = asarray([int(c) for c in C])
		
		C = asarray([self.ElpasoOrder[i][C[i]] for i in range(len(C))]).astype(float)

		if self.debug:		
			print 'Descriptor: ',C
		return	C
class QuantumNumberAWS(object):
	def __init__(self,column = 1,group = 1, row = 1):
		self.QuantumNumbersBCK = {\
		1:[1,1,1],2:[1,1,2],3:[2,1,1],4:[2,1,2]\
		,5:[2,2,1],6:[2,2,2],7:[2,2,3],8:[2,2,4],9:[2,2,5],10:[2,2,6]\
		,11:[3,1,1],12:[3,1,2]\
		,13:[3,2,1],14:[3,2,2],15:[3,2,3],16:[3,2,4],17:[3,2,5],18:[3,2,6]\
		,19:[4,1,1],20:[4,1,2]\
		,21:[4,3,1],22:[4,3,2],23:[4,3,3],24:[4,3,4],25:[4,3,5],26:[4,3,6],27:[4,3,7],28:[4,3,8],29:[4,3,9],30:[4,3,10]\
		,31:[4,2,1],32:[4,2,2],33:[4,2,3],34:[4,2,4],35:[4,2,5],36:[4,2,6]\
		,37:[5,1,1],38:[5,1,2]\
		,39:[5,3,1],40:[5,3,2],41:[5,3,3],42:[5,3,4],43:[5,3,5],44:[5,3,6],45:[5,3,7],46:[5,3,8],47:[5,3,9],48:[5,3,10]\
		,49:[5,2,1],50:[5,2,2],51:[5,2,3],52:[5,2,4],53:[5,2,5],54:[5,3,6]\
		,55:[6,1,1],56:[6,1,2]\
		,57:[6,4,1],58:[6,4,2],59:[6,4,3],60:[6,4,4],61:[6,4,5],62:[6,4,6],63:[6,4,7],64:[6,4,8],65:[6,4,9],66:[6,4,10],67:[6,4,11],68:[6,4,12],69:[6,4,13],70:[6,4,14],71:[6,4,15]\
			   ,72:[6,3,2],73:[6,3,3],74:[6,3,4],75:[6,3,5],76:[6,3,6],77:[6,3,7],78:[6,4,8],79:[6,3,9],80:[6,3,10]\
		,81:[6,2,1],82:[6,2,2],83:[6,2,3],84:[6,2,4],85:[6,3,5],86:[6,3,6]
		,87:[7,1,1],88:[7,1,2]\
		,89:[7,4,1],90:[7,4,2],91:[7,4,3],92:[7,4,4],93:[7,4,5],94:[7,4,6],95:[7,4,7],96:[7,4,8],97:[7,4,9],98:[7,4,10],99:[7,4,11],100:[7,4,12],101:[7,4,13],101:[7,4,14],102:[7,4,14],103:[7,4,15]\
			   ,104:[7,3,2],105:[7,3,3],106:[7,3,4],107:[7,3,5],108:[7,3,6],109:[7,3,7],110:[7,3,8],111:[7,3,9],112:[7,3,10]\
		,113:[7,2,1],114:[7,2,2],115:[7,2,3],116:[7,2,4],117:[7,2,5],118:[7,3,6]}
		
		self.QuantumNumbers = {\
		1:[1,1,1],2:[1,1,2]
		,3:[2,1,1],4:[2,1,2]\
		,5:[2,2,3],6:[2,2,4],7:[2,2,5],8:[2,2,6],9:[2,2,7],10:[2,2,8]\
		,11:[3,1,1],12:[3,1,2]\
		,13:[3,2,3],14:[3,2,4],15:[3,2,5],16:[3,2,6],17:[3,2,7],18:[3,2,8]\
		,19:[4,1,1],20:[4,1,2]\
		,31:[4,2,3],32:[4,2,4],33:[4,2,5],34:[4,2,6],35:[4,2,7],36:[4,2,8]\
		,21:[4,3,9],22:[4,3,10],23:[4,3,11],24:[4,3,12],25:[4,3,13],26:[4,3,14],27:[4,3,15],28:[4,3,16],29:[4,3,17],30:[4,3,18]\
		,37:[5,1,1],38:[5,1,2]\
		,49:[5,2,3],50:[5,2,4],51:[5,2,5],52:[5,2,6],53:[5,2,7],54:[5,3,8]\
		,39:[5,3,9],40:[5,3,10],41:[5,3,11],42:[5,3,12],43:[5,3,13],44:[5,3,14],45:[5,3,15],46:[5,3,16],47:[5,3,17],48:[5,3,18]\
		,55:[6,1,1],56:[6,1,2]\
		,81:[6,2,3],82:[6,2,4],83:[6,2,5],84:[6,2,6],85:[6,3,7],86:[6,3,8]
			   ,72:[6,3,10],73:[6,3,11],74:[6,3,12],75:[6,3,13],76:[6,3,14],77:[6,3,15],78:[6,4,16],79:[6,3,17],80:[6,3,18]\
		,57:[6,4,19],58:[6,4,20],59:[6,4,21],60:[6,4,22],61:[6,4,23],62:[6,4,24],63:[6,4,25],64:[6,4,26],65:[6,4,27],66:[6,4,28],67:[6,4,29],68:[6,4,30],69:[6,4,31],70:[6,4,32],71:[6,4,33]\
		,87:[7,1,1],88:[7,1,2]\
		,113:[7,2,3],114:[7,2,4],115:[7,2,5],116:[7,2,6],117:[7,2,7],118:[7,3,8]\
			   ,104:[7,3,10],105:[7,3,11],106:[7,3,12],107:[7,3,13],108:[7,3,14],109:[7,3,15],110:[7,3,16],111:[7,3,17],112:[7,3,18]\
		,89:[7,4,19],90:[7,4,20],91:[7,4,21],92:[7,4,22],93:[7,4,23],94:[7,4,24],95:[7,4,25],96:[7,4,26],97:[7,4,27],98:[7,4,28],99:[7,4,29],100:[7,4,30],101:[7,4,31],101:[7,4,32],102:[7,4,14],103:[7,4,33]}

		for j in range(len(self.QuantumNumbers)):
			i =j+1
			self.QuantumNumbers[i][0] = self.QuantumNumbers[i][0] * column
			self.QuantumNumbers[i][1] = self.QuantumNumbers[i][1] * group	
			self.QuantumNumbers[i][2] = self.QuantumNumbers[i][2] * row 
		self.tag = 'wy'

	def descibe(self,name):
		#C = name.split('-')[0:-2]
		print 'Unrefined: ', name
		name = name[1::]
		C = re.split('0|1|2|3|4|5|6|7|8|9',name)
		print 'Splited: ', C
		C = filter(None, C)	
		print 'Filtered: ',C	
		C = [parse_poscar.periodictable_numbers[c[1::]] for c in C]
		print 'Refined: ', C		
		C = asarray([float(c) for c in C])
		C = asarray([self.QuantumNumbers[c] for c in C])
		print 'Descriptor: ',C
		return	C
class EwaldForce(object):
	def __init__(self,size):
		self.size = size

	
	def descibe(self,structure):
		C = zeros(self.size**2+ 3*self.size).reshape(self.size,self.size + 3)
		tmp = EwaldSummation(structure).total_energy_matrix
		C[: tmp.shape[0],: tmp.shape[1]]  = tmp
		F = EwaldSummation(structure).forces
		C[:,0:-3] = self.SortMatrix(C[:,0:-3])
		C[: F.shape[0],-3:] =F
		C[:,-3:]= self.SortForce(C[:,-3:])
		return	C

	def SortMatrix(self,X):
		dimX = len(X)
		NormList = zeros(dimX)
		for i in range(0,dimX):
			NormList[i] = linalg.norm(X[i])
		ArgList = argsort(NormList)[::-1]
		X[:,:] = X[:,ArgList]
		X[:,:] = X[ArgList,:]
		self.arg=ArgList
		return X
	def SortForce(self,X):
		X[:,:] = X[self.arg,:]
		return X
class Coulumb(object):
	def __init__(self,size):
		self.size = size

	
	def descibe(self,coords,ocupationList,cell = eye(3)):
		
		C = zeros(self.size**2).reshape(self.size,self.size)
		for i in range(0,len(ocupationList)):			  
				for k in range(i,len(ocupationList)):
						if k == i: 
							C[i,k] = (ocupationList[i]**2.4) /2
						else:			
							C[i,k] = ocupationList[i] * ocupationList[k]\
							/linalg.norm(dot(cell,(coords[i]-coords[k])))
							C[k,i] = C[i,k]
			
		return self.SortMatrix(C)
		
		
	def SortMatrix(self,X):
		dimX = len(X)
		NormList = zeros(dimX)
		for i in range(0,dimX):
			NormList[i] = linalg.norm(X[i])
		NormList = argsort(NormList)[::-1]
		X[:,:] = X[:,NormList]
		X[:,:] = X[NormList,:]
		return X

class CoulumbFast(object):
	def __init__(self,size):
		self.size = size

	
	def descibe(self,coords,ocupationList):
		print coords.shape
		C = zeros([self.size,self.size])
		ocupationList = asarray(ocupationList)
		xCoords,yCoords = coords[:,newaxis],coords[newaxis,:]
		xOL,yOL = ocupationList[:,newaxis],ocupationList[newaxis,:]
		
		C_OL =xOL*yOL
		print C_OL.shape
		C_dist = xCoords-yCoords
		print C_dist.shape
		C_dist = linalg.norm(C_dist,axis = 2)
		print C_dist.shape
		
		r = range(len(ocupationList))
		C_dist[r,r] = 2*ones(len(ocupationList))
		C_OL[r,r] = C_OL[r,r]**1.2
		C_OL = C_OL/C_dist
		print C_OL.shape 
		C[0:len(C_OL),0:len(C_OL)] = C_OL
		print C
		return self.SortMatrix(C)
		
		
	def SortMatrix(self,X):
		dimX = len(X)
		NormList = zeros(dimX)
		
		NormList = linalg.norm(X,axis = 1)
		NormList = argsort(NormList)[::-1]
		X[:,:] = X[:,NormList]
		X[:,:] = X[NormList,:]
		return X
	
class CoulumbPTP(object):
	def __init__(self,size):
		self.size = size

	
	def descibe(self,coords,ocupationList):
		print coords.shape
		C = zeros([self.size,self.size])
		ocupationList = asarray([PTP[o] for o in ocupationList])
		xCoords,yCoords = coords[:,newaxis],coords[newaxis,:]
		xOL,yOL = ocupationList[:,newaxis],ocupationList[newaxis,:]
		
		C_OL = sum(xOL*yOL,axis = 2)
		
		C_dist = linalg.norm(xCoords-yCoords,axis = 2)
		
		fill_diagonal(C_dist, 2)
		C[0:len(C_OL),0:len(C_OL)] = C_OL/C_dist
		print C
		return self.SortMatrix(C)
		
		
	def SortMatrix(self,X):
		dimX = len(X)
		NormList = zeros(dimX)
		for i in range(0,dimX):
			NormList[i] = linalg.norm(X[i])
		NormList = argsort(NormList)[::-1]
		X[:,:] = X[:,NormList]
		X[:,:] = X[NormList,:]
		return X







class AbsSinColumbOld(object):
	def __init__(self,size):
		self.size = size

	
	def descibe(self,cell,coords,ocupationList):
		
		C = zeros(self.size**2).reshape(self.size,self.size)
		for i in range(0,len(ocupationList)):			  
				for k in range(i,len(ocupationList)):
						if k == i: 
							C[i,k] = (ocupationList[i]**2.4) /2
						else:
			
							C[i,k] = math.pi*ocupationList[i] * ocupationList[k]/linalg.norm(dot(cell,absolute(sin(math.pi*(coords[i]-coords[k])))))
							C[k,i] = C[i,k]
			
		return self.SortMatrix(C)

	def SortMatrix(self,X):
		dimX = len(X)
		NormList = zeros(dimX)
		for i in range(0,dimX):
			NormList[i] = linalg.norm(X[i])
		NormList = argsort(NormList)[::-1]
		X[:,:] = X[:,NormList]
		X[:,:] = X[NormList,:]
		return X

class Sin2ColumbOld(object):
	def __init__(self,size):
		self.size = size

	
	def descibe(self,cell,coords,ocupationList):
		
		C = zeros(self.size**2).reshape(self.size,self.size)
		for i in range(0,len(ocupationList)):			  
				for k in range(i,len(ocupationList)):
						if k == i: 
							C[i,k] = (ocupationList[i]**2.4) /2

						else:
			
							C[i,k] = ocupationList[i] * ocupationList[k]/linalg.norm(dot(cell,(sin(math.pi*(coords[i]-coords[k])))**2))

							C[k,i] = C[i,k]
		#print self.SortMatrix(C)
		return self.SortMatrix(C)

	def SortMatrix(self,X):
		dimX = len(X)
		NormList = zeros(dimX)
		for i in range(0,dimX):
			NormList[i] = linalg.norm(X[i])
		NormList = argsort(NormList)[::-1]
		X[:,:] = X[:,NormList]
		X[:,:] = X[NormList,:]
		return X



class Sin2ColumbTest(object):
	def __init__(self,size):
		self.size = size

	
	def descibe(self,cell,coords,ocupationList):
		
		C = zeros(self.size**2).reshape(self.size,self.size)
		for i in range(0,len(ocupationList)):			  
				for k in range(i,len(ocupationList)):
						if k == i: 
							C[i,k] = ocupationList[i]/2

						else:
			
							C[i,k] = sqrt(ocupationList[i] * ocupationList[k]) /linalg.norm(dot(cell,(sin(math.pi*(coords[i]-coords[k])))**2))

							C[k,i] = C[i,k]
		#print self.SortMatrix(C)
		return self.SortMatrix(C)

	def SortMatrix(self,X):
		dimX = len(X)
		NormList = zeros(dimX)
		for i in range(0,dimX):
			NormList[i] = linalg.norm(X[i])
		NormList = argsort(NormList)[::-1]
		X[:,:] = X[:,NormList]
		X[:,:] = X[NormList,:]
		return X
	

class ExtendedSin2ColumbTest(object):
	def __init__(self,size,sumlen = 1,debug=False):
		self.debug = debug
		self.size = size
		self.sumlen = sumlen
		self.boxSize = self.sumlen + 1
		set_printoptions(threshold= 'nan')
	
	def descibe(self,cell,coords,ocupationList):
		
		
		C = zeros(self.boxSize**6*self.size**2).reshape(self.size*self.boxSize**3,self.size*self.boxSize**3)
		if self.debug == True:
			print 'Cell: ', cell
			print 'Coords: ', coords
		
		tmp = zeros(3*len(coords)*self.boxSize**3).reshape(len(coords),self.boxSize**3,3)
		tmp2 = zeros(len(coords)*self.boxSize**3)
		for i in range(len(coords)):
			tmp[i] = self.Extend(coords[i,:])
			for j in range(self.boxSize**3):
				tmp2[(self.boxSize**3)*(i):(self.boxSize**3)*(i+1)]\
				 = [ocupationList[i] for x in range(self.boxSize**3)]
		
		coordsExtended = asarray([val for sublist in tmp for val in sublist])
		ocupationListExtended = asarray(tmp2)
		coordsExtended = coordsExtended/2
		cell = cell*2
		if self.debug == True:
			print 'Coordlen: ',len(coords)
			print 'occLen: ', len(ocupationList)
			print 'CoordsExtended: ',coordsExtended
			print 'Cell: ', cell
		
		for i in range(0,len(ocupationListExtended)):			  
			for k in range(0,len(ocupationListExtended)):
					if k == i: 
						C[i,k] = sqrt((ocupationListExtended[i]**2.4)/abs(linalg.det(cell)))
					else:
						C[i,k] = sqrt(ocupationListExtended[i] * ocupationListExtended[k]\
						/linalg.norm(dot(cell,(sin(math.pi*(coordsExtended[i]-coordsExtended[k])))**2)))
					

		#C = sum(C, axis=0)
		if self.debug == True:
			print self.SortMatrix(C)
		return self.SortMatrix(C)


	def Extend(self,x):
		cellRange = range(0,self.boxSize)
		y = zeros(3*self.boxSize**3).reshape(self.boxSize**3,3)
		print x
		I = eye(3)
		for i in cellRange:
			for j in cellRange:
				for k in cellRange:
					index =  + (i)*(self.boxSize)**2 + (j)*(self.boxSize) + k
					y[index,:] = x + i*I[0,:]+j*I[1,:]+k*I[2,:]
		return y
		
	def SortMatrix(self,X):
		dimX = len(X)
		NormList = zeros(dimX)
		for i in range(0,dimX):
			NormList[i] = linalg.norm(X[i])
		NormList = argsort(NormList)[::-1]
		X[:,:] = X[:,NormList]
		X[:,:] = X[NormList,:]
		return X

	



class MaxSin2Columb(object):
	def __init__(self,size):
		self.size = size

	
	def descibe(self,cell,coords,ocupationList):
		
		C = zeros(self.size**2).reshape(self.size,self.size)
		for i in range(0,len(ocupationList)):			  
				for k in range(i,len(ocupationList)):
						if k == i: 
							#C[i,k] = (ocupationList[i]**2.4) /2
							C[i,k] = (ocupationList[i]**1.0)/2
						else:			
							C[i,k] = ocupationList[k]/linalg.norm(dot(cell,(sin(math.pi*(coords[i]-coords[k])))**2))
							C[k,i] = C[i,k]
		#print self.SortMatrix(C)
		return self.SortMatrix(C)


	def SortMatrix(self,X):
		dimX = len(X)
		NormList = zeros(dimX)
		for i in range(0,dimX):
			NormList[i] = linalg.norm(X[i])
		NormList = argsort(NormList)[::-1]
		X[:,:] = X[:,NormList]
		X[:,:] = X[NormList,:]
		return X
