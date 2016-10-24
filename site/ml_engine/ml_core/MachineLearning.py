#!/usr/bin/env python


import sys
import parse_poscar
import scipy
import multiprocessing
import parallelize as pa
from Kernels import *
from scipy import stats as st
from numpy import *
from numpy import linalg, zeros, array, argsort, bincount
from copy import deepcopy
from collections import Counter
from os import listdir,makedirs
from os.path import isfile, join,exists,split
from ObjectHandler import *


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ML(self,path,sigma = 1,oFPunish = 1E-6,linOfPunish = 10,
#	randomize = None,debug = None, ker = None,num_cores = None)
#
# The ML class is constructed to preform kernel ridge regression (KRR). The classs is both used to train on a 
# machine and to evaluate new data points. 
#
# A perliminary linear regression of the static contribution from all species in the crystal/molecule
# is done on the the data set before the KRR is applied. This can be turned of by setting linOfPunish to None.
#---------------------------------------------------------------------------------------------------------------
# Parameters,
# path: str or dictionary
# The path to a pickle file containing the data or a dictionary containing the data. 
#
# sigma: float or list of float (depending on the kernel)
# Specifies the parameters for the kernel 
#
# oFPunish: float
# Overfitting punish for the KRR
#
# linOfPunish: float or None
# Overfitting punish for the preliminary linear regression. Only mean shift will be preformed if set to None
#
# randomize: int or None
# Selects the random seed for randomizing the order in the data set. No randomization of the order will be done 
# if set to None
#
# debug: int or None
# Used to display debug mesage.
#
# ker: kernel class
# Specifies the type of kernel to be used.
#
# num_cores: int or None
# Specifies the number of cores to to be used (Multiprocessing not yet implimented).
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class ML(object):

	def __init__(self,path,sigma = 1,oFPunish = 1E-6,linOfPunish = 10,
				 randomize = None,debug = None, ker = None ,num_cores = None):
		
		self.debug = debug
		if type(path) is str:
			self.o = importData(path)		
			print 'path: ',path
		elif type(path) is dict:
			self.o = path
		

		self.alpha = None
	
		if num_cores == None:
			self.num_cores = multiprocessing.cpu_count() -2
		else:
			self.num_cores = num_cores

		self.sigma = sigma 
		self.oFPunish = oFPunish
		self.linOfPunish = linOfPunish
		if ker:
			self.ker = ker
		else:
			self.ker =  LapKer()
		
		self.DistMatrix = None 
		self.ker.SetSigma(self.sigma)

		if randomize != None:
			self.rng = random.RandomState(randomize)
			self.o = randOrdOfTraiSet(self.o,self.rng)
		print 'sigma: ',self.sigma,'oFPunish: ',self.oFPunish,'linOfPunish: ',self.linOfPunish
		print 'kertype: ', type(self.ker)
		set_printoptions(threshold= 'nan')

		
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# SetSigma(sigma)
# Sets a new set of parameters for the kernel
#---------------------------------------------------------------------------------------------------------------
# Parameters,
# sigma: float or list of float (depending on the kernel)
# New parameters for the kernel
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	def SetSigma(self,sigma):
		self.sigma = sigma
		self.ker.SetSigma(sigma)
 		

				
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# CreateDistMatrix(BlockSize = 500)
# Creats a matrix with the distance between all data points. This function is used to speed up the evaluation
# for kernels which uses a metric. The distance matrix will be stored in the class for use in evaluation of the 
# kernel. The distance calculation is split up into subblocks in order to save memmory.
#---------------------------------------------------------------------------------------------------------------
# Parameters,
# blockSize: int
# Determines the size of the subblocks for which the matrix is evaluated. Larger block sizes means more memmory 
# use.
#---------------------------------------------------------------------------------------------------------------
# Returns,
# ndarray
# A ndarray containing the distances of all data points in the class.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	def CreateDistMatrix(self,X1 = None,X2 = None, blockSize = 100):
		
		if X1 != None and X2!= None:
			d1,d2 =X1,X2
		else:
			self.o['X'] = asarray(list(self.o['X']))
			d1,d2 =self.o['X'],self.o['X']
		
 		DistMatrix = zeros(len(d1)*len(d2)).reshape(len(d1),len(d2))
		for i in range(0,len(d2),blockSize):
			for j in range(0,len(d1),blockSize):
				if i < len(d2) - blockSize:
					ii = i + blockSize
				else:
					ii = len(d2)
				if j < len(d1) - blockSize:
					jj = j + blockSize
				else:
					jj = len(d1)

				Xilen, Xjlen = asarray(d1[i:ii].shape), asarray(d2[j:jj].shape)

				DistMatrix[j:jj,i:ii] =  self.ker.computeDist(\
					d1[j:jj].reshape(jj - j,1,prod(Xjlen[1::]) )\
					,d2[i:ii].reshape(1,ii - i,prod(Xilen[1::])))
				
		if X1 == None or X2 == None:
			self.DistMatrix = DistMatrix
		return DistMatrix

	
	
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ChangeDistMatrix(D)
# Utility function for changing the distance matrix of the class, with a distance matrix outside of the class.
#---------------------------------------------------------------------------------------------------------------
# Parameters,
# D: ndarray
# A matrix containing all distances between all data points the class.
#---------------------------------------------------------------------------------------------------------------
# Returns,
# ndarray
# The old distance matrix. 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	def ChangeDistMatrix(self,D):
		tmp = self.DistMatrix
		self.DistMatrix = D
		return tmp

	
	
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# CreateKerMatrix()
# Evaluates kernel matrix. Note that this requires that there already exists a distance matrix in the class.
# A distance matrix can be created by calling CreateDistMatrix, or be added from outside of the class by 
# ChangeDistMatrix. This is a support function and should not be called.
#---------------------------------------------------------------------------------------------------------------
# Returns,
# ndarray
# Returns the kernel matrix.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	def CreateKerMatrix(self):
		return self.ker.computeKer(self.DistMatrix)

	
	
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# preAtomRegression()
# Preforms a preliminary linear regression on the data set. The regression is done with respect to the static
# contribution of all species in a data point. This is a support function and should not be called.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	def preAtomRegression(self):
		I = eye(len(self.o['atomCount'][0]))
		
		I[0,0] = 0 
		
		self.atomRegConst = linalg.solve(dot(transpose(self.o['atomCount']),self.o['atomCount'])\
		+ self.linOfPunish*I, dot(transpose(self.o['atomCount']), self.o['T']))

		

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FindAlpha()
# This is the main function for training the machine. Note it will generate a distance matrix if a previus 
# matrix does not exist (self.DistMatrix = None). Otherwise, it will use the distance matrix which already exits
# in the class. The distance matrix can be changed by using either ChangeDistMatrix or CreateDistMatrix.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	def FindAlpha(self):
		
		if self.DistMatrix == None:
			self.CreateDistMatrix()
	
		if self.linOfPunish != None:
			self.preAtomRegression()
			a = self.o['T'] - dot(self.o['atomCount'],transpose(self.atomRegConst))
		else:
			self.atomRegConst = mean(self.o['T'])
			a = self.o['T'] - self.atomRegConst

		print 'Standard diviation of dressed Atom: ', sqrt(mean(abs(a)**2))

		I = eye(len(self.o['X']))	
		kerMatrix = self.CreateKerMatrix()
		self.o['alpha'] = linalg.solve(kerMatrix + self.oFPunish*I,a) 
		

		
		
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Energy(x,atomCount = None)
# Evaluates a new data point with the machine. Note that FindAlpha must be called atleast one time before this 
# function.
#---------------------------------------------------------------------------------------------------------------
# Parameters,
# x: array_like
# The key of the new data point.
#
# atomCount: array_like or None
# A 118-by-1 array where each element is the amount of the species which has the same atomic number as the index
# of the site. atomCount[0] should always be set to equal to one. Note that this should be None if linOfPunish 
# is None.
#---------------------------------------------------------------------------------------------------------------
# Returns,
# float
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	def Energy(self,x,atomCount = None):

		xlen, Xlen = asarray(x.shape), asarray(self.o['X'].shape) 
		d = self.ker.computeDist(x.reshape(1,prod(xlen[0::])),\
			self.o['X'].reshape(Xlen[0],prod(Xlen[1::])))
		
		return self.Energy_KnownDist(d,atomCount)	
	
	
	
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Energy_KnownDist(d,atomCount = None)
# The same thing as Energy, but takes the distance between the key of a data point and all keys in the training
# set.
#---------------------------------------------------------------------------------------------------------------
# Parameters,
# d: array_like
# Array with the distance between the key of a data point and all keys in the training set.
#
# atomCount: array_like or None
# A 118-by-1 array where each element is the amount of the species which has the same atomic number as the index
# of the site. atomCount[0] should always be set to equal to one. Note that this should be None if linOfPunish 
# is None.
#---------------------------------------------------------------------------------------------------------------
# Returns,
# float
# The estimated value of that point.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	def Energy_KnownDist(self,d,atomCount = None):
		if self.linOfPunish != None:	
			linEnergy = dot(atomCount,self.atomRegConst)
		else:
			linEnergy = self.atomRegConst
		NonLinEnergy = dot(self.o['alpha'],self.ker.computeKer(d))
		return linEnergy + NonLinEnergy 


	
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Energy_KnownDist_list(d,atomCount = None)
# Same thing as Energy_KnownDist, but evaluates a list data points.
#---------------------------------------------------------------------------------------------------------------
# Parameters,
# d: array_like
# Array of arrays with the distances. See Energy_KnownDist
#
# atomCount: array_like or None
# A list of atom counts or None. See Energy_KnownDist
#---------------------------------------------------------------------------------------------------------------
# Returns,
# ndarray of floats
# List with estimated value of the data points.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	def Energy_KnownDist_list(self,d,atomCount = None):
		if atomCount == None:
			atomCount = [None for i in range(len(x))]		
		
		return asarray([self.Energy_KnownDist(d[i],atomCount[i]) for i in range(len(d))])
		
		


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Energy_list(x,atomCount = None)
# Same thing as Energy, but evaluates a list data points.
#---------------------------------------------------------------------------------------------------------------
# Parameters,
# x: array_like
# Array of keys with keys of the data points. See Energy
#
# atomCount: array_like or None
# A list of atom counts or None. See Energy
#---------------------------------------------------------------------------------------------------------------
# Returns,
# ndarray of floats
# List with estimated value of the data points.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	def Energy_list(self,x,atomCount = None):
	
		if atomCount == None:
			atomCount = [None for i in range(len(x))]			
		return asarray([self.Energy(x[i],atomCount[i]) for i in range(len(x))])
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Energy_list(x,atomCount = None)
# Same thing as Energy, but evaluates a list data points.
#---------------------------------------------------------------------------------------------------------------
# Parameters,
# x: array_like
# Array of keys with keys of the data points. See Energy
#
# atomCount: array_like or None
# A list of atom counts or None. See Energy
#---------------------------------------------------------------------------------------------------------------
# Returns,
# ndarray of floats
# List with estimated value of the data points.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	def Energy_list_parallel(self,x,atomCount = None,cores = None):
	
		if cores:
			c = cores
		else:
			c = self.num_cores
		if atomCount == None:
			atomCount = [None for i in range(len(x))]		
			
		return pa.parallelize(self.Energy,zip(*[x,atomCount]),c)

		
			#return asarray([self.Energy(x[i],None) for i in range(len(x))])
		
		
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ErrorVals(X,E,N = None ,atomCount = None)
# Gives the following error estimates of the preformence of the machine on a given data set in a dictionary:
# 'MAE': Mean absolute errror, 'RMSE': Root mean squared error, 'R^2': R squared and 'MaxErr': Maximal error. 
# the dictionary will also contain 'errIdx', which is the indices of all data points sorted by its error form 
# best to worst.
#---------------------------------------------------------------------------------------------------------------
# Parameters,
# X: array_like
# Array of keys with keys of the data points. See Energy
#
# E: array_like
# Array with the true values of the data points.
# 
# N: array_like or None
# List with the number of atoms for each data point or None. The error estimates will be given in unit/Atom if
# the list is provided and in unit if None.
#
# atomCount: array_like or None
# A list of atom counts or None. See Energy_KnownDist
#
#---------------------------------------------------------------------------------------------------------------
# Returns,
# dict
# The dictionary contains: float, float, folat float, list of ints
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	def ErrorVals(self,X,E,N = None ,atomCount = None,cores = None):
		if cores == -1:
			if N != None:
				E_atom = E / N	
				EEst_atom = self.Energy_list(X,atomCount)/N
			else: 
				E_atom = E
				EEst_atom = self.Energy_list(X,atomCount)
		else:
			if N != None:
				E_atom = E / N	
				EEst_atom = self.Energy_list_parallel(X,atomCount,cores)/N
			else: 
				E_atom = E
				EEst_atom = self.Energy_list_parallel(X,atomCount,cores)
			

		EDiff_atom = (EEst_atom - E_atom) 
		MAE = mean(abs(EDiff_atom)) 
		RMSE = sqrt(mean(abs(EDiff_atom)**2))
		maxerr = max(abs(EDiff_atom))
		RSqrd = st.pearsonr(EEst_atom,E_atom)[0]**2
		errIdx = argsort(abs(EDiff_atom)) 
		print 'RMSE:', RMSE, 'MAE:', MAE,'R^2: ',RSqrd , 'MaxErr: ', maxerr 
		return {'RMSE': RMSE, 'MAE': MAE,'R^2': RSqrd , 'MaxErr': maxerr,'errIdx':errIdx}
	



	

