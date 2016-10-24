#!/usr/bin/env python

import pickle
import sys
import parse_poscar
import scipy
import time
import random
from List2Energies import *
from DiffEvolution import *
from MachineLearning import *
from scipy import optimize
from scipy.optimize import minimize
from numpy import *
from numpy import linalg, zeros, array, argsort, bincount
from copy import deepcopy
from collections import Counter
from os import listdir,makedirs
from os.path import isfile, join,exists,split

def genInitPop(SSpace,popSize):
	
	pop = []
	for i in range(popSize):
		ind = []
		for s in SSpace:
			ind.append(random.choice(s))
		pop.append(ind)
	return pop


def prepEve(target,NOpt,SSpace, m):
	#m 	= 	loadMachine(True)
	initPop = genInitPop(SSpace,NOpt + 10)
	d = DifferentialEvolution(initPop,lambda x: fitFun(x,target,m),
							  F = 1, CR = None,
							  metricFun = lambda a,b,c,F: metricFun(a,b,c,F,SSpace))
	popSet = set(tuple(ip) for ip in initPop)

	return d, popSet

def FindTargetE(target,NOpt,SSpace):

#	m 	= 	loadMachine(True)
#	initPop = genInitPop(SSpace,NOpt + 10)
#	d = DifferentialEvolution(initPop,lambda x: fitFun(x,target,m),
#							  F = 1, CR = None,
#							  metricFun = lambda a,b,c,F: metricFun(a,b,c,F,SSpace))
#	popSet = set(tuple(ip) for ip in initPop)
	for i in range(100):
		d.nextPop()
		for ip in d.pop:
			popSet.add(tuple(ip))
			
		E = [fitFun(a,target,m) for a in popSet]
		optIdxs = argsort(E)[0:NOpt]
		optInd = asarray(list(popSet))[optIdxs]
		popSet = set(tuple(ip) for ip in optInd)
		optE = [fitFun(ip,None,m) for ip in optInd]
		#print 'Gen:',i,' || Fit: ',optE, ' || '\
		#	, optInd
		minE = min(optE)
		print minE
def fitFun(ind,target,m):
	x,a = O2Desc(ind)
	e = m.Energy(x,a) 
	if target:
		return abs(target -  e)
	else:
		return e

def metricFun(a,b,c,F,SSpace):
	def S2D(S):
		N = asarray([PN[s] for s in S]).astype(int)
		return asarray([PTP[n] for n in N])
		
	a_d,b_d,c_d = S2D(a),S2D(b),S2D(c)
	S_d = [S2D(s) for s in SSpace]

	
	X = a_d + F * (b_d - c_d)
	newX = empty(len(X)).astype(str)
	for i in range(len(X)):
		tmp = X[i,newaxis,:] - asarray(S_d[i])
		newX[i] = SSpace[i][argmin(linalg.norm(tmp,1,axis = -1))]
	return newX
	
def main():
	FindTargetE(-1,50,[	['H'],
						['Al','Na','K','Cs','Rb'],
						['Al','Na','K','Cs','Rb','Ca','Mg','Be','Sr','Ba'],
						['O','F','Cl','Br','I']])

 
if __name__ == "__main__":
	main()
	
	print 'end'   
