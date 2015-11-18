#!/usr/bin/env python

import pickle
import sys
import parse_poscar
import scipy
from scipy import optimize
from scipy.optimize import minimize
from numpy import *
from numpy import linalg, zeros, array, argsort, bincount
from copy import deepcopy
from collections import Counter
from os import listdir,makedirs
from os.path import isfile, join,exists,split
#from pymatgen.analysis.ewald import EwaldSummation
#from pymatgen import Lattice, Structure

class SphericalKer(object):

	def __init__(self):
		self.rad = True


	def SetSigma(self,sigma):
		self.const = sigma

	def computeKer(self,d):
		K = zeros(d.shape)
		idx = where(d < self.const)
		K[idx] = 1 - (3*d[idx])/(2*self.const) + (1/2) * ((d[idx]/self.const)**3)

		return K
		
	def computeDist(self,x1,x2):
		return linalg.norm(x1-x2,1,axis = len(x1.shape)-1)
	
	
	
class LogisticKer(object):

	def __init__(self):
		self.rad = True

	def SetSigma(self,sigma):
		self.lapconst = 1.0/(sigma)

	def computeKer(self,d):
		return 1/(exp(-d*self.lapconst) + 2 + exp(d*self.lapconst))

	def computeDist(self,x1,x2):
		return linalg.norm(x1-x2,1,axis = len(x1.shape)-1)	
	
class TestKer(object):

	def __init__(self):
		self.rad = False


	def SetSigma(self,sigma):
		self.const = sigma

	def computeKer(self,d):
		return d
		
	def computeDist(self,x1,x2):
		d =  abs(x1-x2)
		K = zeros(d.shape)
		idx = d != 0
		K[idx] = sin((pi*d[idx])/(2*self.const))/d[idx]
		
		K = prod(K,axis = len(d.shape)-1)
		
		return K
	
class PolyKer(object):

	def __init__(self):
		self.rad = False


	def SetSigma(self,sigma):
		self.const = sigma

	def computeKer(self,d):
		return d
		
	def computeDist(self,x1,x2):
		d =  (x1*x2)*self.const
		d = sum(d,axis = len(d.shape)-1)
		return d
	
class Test2Ker(object):

	def __init__(self):
		self.rad = False


	def SetSigma(self,sigma):
		self.const = sigma

	def computeKer(self,d):
		return d
		
	def computeDist(self,x1,x2):
		d = abs(x1-x2)
		K = zeros(d.shape)
		idx = d < self.const
		
		K[idx] = 1 - sqrt(d[idx]/self.const)
		K = prod(K,axis = len(d.shape)-1)
		
		return K
	
class QuaterCircleKer(object):

	def __init__(self):
		self.rad = False


	def SetSigma(self,sigma):
		self.const = sigma

	def computeKer(self,d):
		return d
		
	def computeDist(self,x1,x2):
		d = abs(x1-x2)
		K = zeros(d.shape)
		idx = d < self.const
		
		K[idx] = 1 - sqrt(1-(d[idx]/self.const - 1)**2)
		K = prod(K,axis = len(d.shape)-1)
		
		return K	
	
class SubTriangleKer(object):

	def __init__(self):
		self.rad = False


	def SetSigma(self,sigma):
		self.const = sigma

	def computeKer(self,d):
		return d
		
	def computeDist(self,x1,x2):
		d = abs(x1-x2)
		K = zeros(d.shape)
		idx = d < self.const
		
		K[idx] = 1 - sqrt(d[idx]/self.const)
		K = prod(K,axis = len(d.shape)-1)
		
		return K
	
class TriangleKer(object):

	def __init__(self):
		self.rad = False


	def SetSigma(self,sigma):
		self.const = sigma

	def computeKer(self,d):
		return d
		
	def computeDist(self,x1,x2):
		d = abs(x1-x2)
		K = zeros(d.shape)
		idx = where(d < self.const)
		
		K[idx] = 1 - d[idx]/self.const
		K = prod(K,axis = len(d.shape)-1)
		
		return K
	
	
class UniformKer(object):

	def __init__(self):
		self.rad = False


	def SetSigma(self,sigma):
		self.const = sigma

	def computeKer(self,d):
		return d
		
	def computeDist(self,x1,x2):
		d = abs(x1-x2)
		K = zeros(d.shape)
		idx = where(d < self.const)
		K[idx] = 1
		K = prod(K,axis = len(d.shape)-1)		
		return K
	
class CircularKer(object):

	def __init__(self):
		self.rad = True

	def SetSigma(self,sigma):
		self.const = sigma

	def computeKer(self,d):
		K = zeros(d.shape)
		idx = where(d < self.const)
		K[idx] = 2 -(2/pi) * arccos(-d[idx]/self.const) + (2/pi)*(d[idx]/ self.const)*sqrt(1-(d[idx]/self.const)**2)

		return K
		
	def computeDist(self,x1,x2):
		return linalg.norm(x1-x2,1,axis = len(x1.shape)-1)

class LapKer(object):

	def __init__(self):
		self.rad = True

	def SetSigma(self,sigma):
		self.lapconst = 1.0/(sigma)

	def computeKer(self,d):
		return exp(-d*self.lapconst)*(self.lapconst/2.)

	def computeDist(self,x1,x2):
		return linalg.norm(x1-x2,1,axis = len(x1.shape)-1)
	
class TriLapKer(object):

	def __init__(self):
		self.rad = True

	def SetSigma(self,sigma):
		self.lapconst = 1.0/(sigma)

	def computeKer(self,d):
		return exp(-d*self.lapconst)

	def computeDist(self,x1,x2,x3):
		d = (linalg.norm(x1-x2,1,axis = len(x1.shape)-1)\
			+linalg.norm(x1-x3,1,axis = len(x1.shape)-1)\
			+linalg.norm(x2-x3,1,axis = len(x1.shape)-1))
		return d
	
class CauchyKer(object):

	def __init__(self):
		self.rad = True

	def SetSigma(self,sigma):
		self.lapconst = 1/sigma

	def computeKer(self,d):
		return (1)/(1 + (self.lapconst*d)**2)

	def computeDist(self,x1,x2):
		return linalg.norm(x1-x2,2,axis = len(x1.shape)-1)

class LinKerL1(object):

	def __init__(self):
		self.rad = True

	def SetSigma(self,sigma):
		self.lapconst = sigma

	def computeKer(self,d):
		return d

	def computeDist(self,x1,x2):
		return linalg.norm(x1-x2,1,axis = len(x1.shape)-1)
	
class LinKerL2(object):

	def __init__(self):
		self.rad = True

	def SetSigma(self,sigma):
		self.lapconst = sigma

	def computeKer(self,d):
		return d

	def computeDist(self,x1,x2):
		return linalg.norm(x1-x2,2,axis = len(x1.shape)-1)
class GaussKer(object):

	def __init__(self):
		self.rad = True

	def SetSigma(self,sigma):
		self.gaussconst = 1.0/(2.0*sigma**2)

	def computeKer(self,d):
		return exp(-d**2*self.gaussconst)*sqrt(self.gaussconst/pi)

	def computeDist(self,x1,x2):
		return linalg.norm(x1-x2,2,axis = len(x1.shape)-1)


class Gaussl1Ker(object):

	def __init__(self):
		self.rad = True

	def SetSigma(self,sigma):
		self.gaussconst = 1.0/(2.0*sigma**2)

	def computeKer(self,d):
		return exp(-d**2*self.gaussconst)

	def computeDist(self,x1,x2):
		return linalg.norm(x1-x2,1,axis = len(x1.shape)-1)

class LapKerOld(object):

	def __init__(self):
		self.rad = True

	def SetSigma(self,sigma):
		self.lapconst = 1.0/(sigma)

	def computeKer(self,x1,x2):
		return math.exp(-linalg.norm(x1-x2,1)*self.lapconst)




class GaussKerOld(object):

	def __init__(self):
		self.rad = True

	def SetSigma(self,sigma):
		self.gaussconst = 1.0/(2.0*sigma**2)

	def computeKer(self,x1,x2):
		return math.exp(-linalg.norm(x1-x2,2)**2*self.gaussconst)
