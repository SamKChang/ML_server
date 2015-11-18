#!/usr/bin/env python

import pickle
import parse_poscar
import sys
from MachineLearning import *
from CustomOpperators import *
from numpy import *
from numpy import linalg, zeros, array, argsort
from copy import deepcopy
from collections import Counter
from os import listdir,makedirs
from os.path import isfile, join,exists

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The following functions used to conveniently handle data used for machine learning or other aplications.
# The data is stored as a dict of list. Each index in one of the lists corresponds to a data point (e.g. molecule
# or crystal). Each key corresponds to a propperty (feature type), eg. coordinates, representation, energy or symetry.
#
# The following properties are nessceary for machine learning are:
# 'X': Representation of the object, e.g the coulomb matrix 
# (http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.109.059802).
# 'T': Property you wish to learn, e.g. formation energy.
# 'N': Number of atoms in the molecule, or unit cell.
# 'Z': Number of each atom type in the molecule or unit cell.
#
# Note that this framework has been constructed mainly for preforming machinelearning on crystals and molecules.
# However, I see no reason why this should not work for other machine learning problems as well. Just set, 'N' as
# a list of ones and 'Z' as a list of list of ones.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# importData(path)
# Imports data from a pickle object (https://docs.python.org/2/library/pickle.html)
#---------------------------------------------------------------------------------------------------------------
# parameters.
# path: stirng
# Path to file
#---------------------------------------------------------------------------------------------------------------
# returns: dictionary
# Dictionary containing data, as discribed above.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def importData(path):
	f = open(path, 'rb')
	o = pickle.load(f)
	f.close()

	try:			
		o['X'] = asarray(o['X']).astype(float)
	except ValueError:

		print 'Warning! X is of type: ' ,type(o['X'])	
	try:
		tmp = asarray([bincount(z) for z in o['Z']])
	except TypeError:
		o['Z'] = asarray([[int(x) for x in z]  for z in o['Z']])
		tmp = asarray([bincount(z) for z in o['Z']])			
	a = zeros(len(o['Z'])*118).reshape(len(o['Z']),118)
	for i in range(len(tmp)):
		a[i,0:tmp[i].shape[0]]= tmp[i]            
		a[i,0] = 1	   
	o['atomCount'] = a
	return o	



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# exportData(path)
# Saves data from in pickle format (https://docs.python.org/2/library/pickle.html)
#---------------------------------------------------------------------------------------------------------------
# parameters.
# o: dictionary
# Dictionary to be saved.
# path: stirng
# Path to where the data is saved.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def exportData(o,path):
	pickle.dump(o, open( path +'.pkl', 'wb' ))

	
	
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# removeIndex(o,r)
# Removes data points specified by a list of indicices.
#---------------------------------------------------------------------------------------------------------------
# parameters.
# o: dictionary
# Data dictionary from where elements are to be removed.
# r: array_like (http://docs.scipy.org/doc/numpy/user/basics.creation.html)
# List containing indicices of the data points to be removed.
#---------------------------------------------------------------------------------------------------------------
# returns: dictionary
# Data dictionary without the elements contained in r.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
def removeIndex(o,r):
	o1 = o.copy()	
	for key, value in o1.iteritems():
		o1[key] = delete(value,r,0)
	return o1



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# pickIndex(o,r)
# Picks data points specified by a list of indicices.
#---------------------------------------------------------------------------------------------------------------
# parameters.
# o: dictionary
# Data dictionary from where elements are to be picked.
# r: array_like (http://docs.scipy.org/doc/numpy/user/basics.creation.html)
# List containing indicices of the data points to be picked.
#---------------------------------------------------------------------------------------------------------------
# returns: dictionary
# Data dictionary with ONLY the datapoints contained specified by r. Note that the dictionary will have the same
# order as r.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def pickIndex(o,r):
	o1 = o.copy()
	for key, value in o1.iteritems():
		o1[key] = value[r]
	return o1



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# splitIndex(o,r)
# Splits a data dictionaries into two dictionaries, specified by a list of indicices.
#---------------------------------------------------------------------------------------------------------------
# parameters.
# o: dictionary
# Data dictionary which is to be splited.
# r: array_like (http://docs.scipy.org/doc/numpy/user/basics.creation.html)
# List containing indicices specifing how the data dictionary is splited.
#---------------------------------------------------------------------------------------------------------------
# returns: dictionary, dictionary
# Returns two data dictionaries, first containg only data points form r. Second containing the remaining data
# points.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def splitIndex(o2,r):
	o1 = pickIndex(o2,r)
	o2 = removeIndex(o2,r)
	return o1,o2



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# appendObj(o,o2)
# Concancates two data dictionaries.
#---------------------------------------------------------------------------------------------------------------
# parameters.
# o: dictionary
# First Data dictionary.
# o2: dictionary
# Second Data dictionary.
#---------------------------------------------------------------------------------------------------------------
# returns: dictionary
# A concanced data dictionary from o and o2.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def appendObj(o,o2):

	for key, value in o.iteritems():
		o[key] =  concatenate((value,o2[key]),axis =0)
	return o



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# removeByFeature(o,feature,featureType)
# removes all data ponits which contains a feature from a feature type.
#---------------------------------------------------------------------------------------------------------------
# parameters.
# o: dictionary
# Data dictionary.
# feature: array_like (http://docs.scipy.org/doc/numpy/user/basics.creation.html)
# The feture to be removed.
# featureType: string
# The feture type from where to look for the feture, e.g. 'X'.
#---------------------------------------------------------------------------------------------------------------
# returns: dictionary
# Data dictionary containing no data points with the feature.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def removeByFeature(o,feature,featureType):  
	try:
		idx = where(o[featureType] == feature)
		o = removeIndex(o,idx)
	except AttributeError:
		pass
	return o

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# pickByFeature(o,feature,featureType)
# picks all data ponits which contains a feature from a feature type.
#---------------------------------------------------------------------------------------------------------------
# parameters.
# o: dictionary
# Data dictionary.
# feature: array_like (http://docs.scipy.org/doc/numpy/user/basics.creation.html)
# The feture to be picked.
# featureType: string
# The feture type from where to look for the feture, e.g. 'X'.
#---------------------------------------------------------------------------------------------------------------
# returns: dictionary
# Data dictionary containing only data points with the feature.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def pickByFeature(o,feature,featureType):  
	try:
		idx = where(o[featureType] == feature)
		o = pickIndex(o,idx)
	except AttributeError:
		pass
	return o

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# randOrdOfTraiSet(o)
# Randomizes the order of a data dictionary.
#---------------------------------------------------------------------------------------------------------------
# parameters.
# o: dictionary
# Data dictionary to be randomized.
#---------------------------------------------------------------------------------------------------------------
# returns: dictionary
# Data dictionary where the order of the data points have been randomized.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def randOrdOfTraiSet(o,rng):

	r = rng.random_sample([len(o['X'])])
	print len(r)
	r = argsort(r)
	o = pickIndex(o,r)
	return o

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# WARNING!
# All functions below are outdated. These should not be used and are therefore not commented.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


def SortMatrix(o):
	lenX = len(o['X'])

	for k in range(0,lenX):
		dimX = len(o['X'][k])
		NormList = zeros(dimX)
		for i in range(0,dimX):
			NormList[i] = linalg.norm(o['X'][k,i],1)
		NormList = argsort(NormList)[::-1]
		o['X'][k,:,:] = o['X'][k,:,NormList]
		o['X'][k,:,:] = o['X'][k,NormList,:]
	return o


def pickAllPolymorfs(o,pick = 'P'):

	Z_sorted = sort(o['Z'])
	tup = tuple(range(1,len(Z_sorted.shape)))
	Z_equal = [any(equal(Z_sorted[range(i) + range(i+1,len(Z_sorted))],Z_sorted[i]).all(tup)) for i in range(len(Z_sorted))]
	Z_where = where(Z_equal)


	if pick == 'R':
		o = removeIndex(o,Z_where)
		print 'Number of non-polymorfs found: ', len(o['Z'])
	else:
		o = pickIndex(o,Z_where)
		print 'Number of polymorfs found: ', len(o['Z'])
	return o
def pickSpaceGroup(o,spaceGroup,Type = 'S'):  
	try:
		tmp = []
		for i in range(len(o['S'])):
			if Type == 'S' or Type == 'R_S':
				if o['S'][i] in  spaceGroup:
					tmp.append(i)
			elif Type == 'Annon' or Type == 'R_Annon':
				if o['Annon'][i] in  spaceGroup:
					tmp.append(i)
		if Type == 'Annon' or Type == 'S':
			o = pickIndex(o,tmp)
		elif Type == 'R_S' or Type == 'R_Annon':
			o = removeIndex(o,tmp)

		if Type == 'S':
			print 'Crystals: ',len(o['T']),'in spacegroup: ', spaceGroup
		elif Type == 'Annon':
			print 'Crystals: ',len(o['T']),'in Annon: ', spaceGroup
		elif Type == 'R_S':
			print 'Crystals: ',len(o['T']),'left after removing spacegroup: ', spaceGroup
		elif Type == 'R_Annon':
			print 'Crystals: ',len(o['T']),'left after removing Annon: ', spaceGroup
	except AttributeError:
		pass
	return o

def pickCommonAtoms(o,n=None):

	A = argsort(sum(o['atomCount'][:,1:],axis=0))[::-1]
	B = sort(sum(o['atomCount'][:,1:],axis=0))[::-1]
	if n == None:
		print 'Order: ', A 
		print 'Amount: ', B
		return o
	tmp = []
	print 'Total amount of atoms removed: ', sum(B[n::])

	for a in A[n::]:
		for i in range(len(o['X'])):
			b = o['Z'][i]
			if a in b:
				tmp.append(i)
	tmp = asarray(unique(tmp))
	tmp2 = []
	for i in range(len(o['X'])):
		if i not in tmp:
			tmp2.append(i)
	o = pickIndex(o,tmp2)
	print 'Crystals left after picking atoms: ', len(o['X']) 
	return o


def removeSmallCrystals(o,r):

	tmp = []	
	for i in range(len(o['N'])):
		if o['N'][i] < r:
			tmp.append(i)
	tmp = asarray(tmp)
	o = removeIndex(o,tmp)
	print 'Crystals: ', len(o['N']), 'left after removing crystals with less than',r,'atoms'
	return o

def pickElements(o,r):

	tmp = []	
	for i in range(len(o['Z'])):
		ElementsInCrystal = True
		for z in o['Z'][i]:
			if z not in r:
				ElementsInCrystal = False
		if ElementsInCrystal == True:
			tmp.append(i)
	tmp = asarray(tmp)
	o = pickIndex(o,tmp)
	print 'Crystals: ', len(o['Z']), 'Containing only the following elements: ', r
	return o




