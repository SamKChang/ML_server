#!/usr/bin/env python

import pickle
import sys
import scipy
from scipy import stats as st
from numpy import *
from copy import deepcopy
from collections import Counter
from os import listdir,makedirs
from os.path import isfile, join,exists,split

 
class DifferentialEvolution(object):
		

		def __init__(self,initPop,fitFun,CR = None,F = 1,metricFun = None):
			self.fitFun = fitFun
			self.metricFun = metricFun
			self.F = F
			self._pop = asarray(initPop)
			self.invShape = asarray(initPop[0]).shape
			self.invDim = int(prod(self.invShape))
			if CR:
				self.CR = CR
			else: 
				self.CR = 1/self.invDim
			self.popSize = len(initPop)
			self._pop = self._pop.reshape(self.popSize,self.invDim)
		def nextPop(self):
				
			newPop = copy(self._pop)
			for j in range(self.popSize):
				ranVal =  range(0, j) + range(j+1, self.popSize)
				random.shuffle(ranVal)
				ranVal = ranVal[0:3]

				tmp = self._pop[ranVal]
				a,b,c = tmp[0],tmp[1],tmp[2]
				if self.metricFun:
					newX = self.metricFun(a,b,c,self.F)
				else:
					newX = a + self.F * (b - c)  
				R = random.randint(0,prod(self.invDim))
				
				 
				for k in range(self.invDim):
					r = random.uniform(0,1)
					if r < self.CR or k == R:
						newPop[j,k] = newX[k]
			popFit = asarray([self.fitFun(x) for x in self._pop])
			newPopFit = asarray([self.fitFun(x) for x in newPop]) 
			comparison = popFit - newPopFit
			
			repIdx = where(comparison > 0)
			self._pop[repIdx] = newPop[repIdx]
			
		@property
		def pop(self):
			return self._pop.reshape(self.popSize,*self.invShape)	





