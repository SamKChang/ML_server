#!/usr/bin/env python

import pickle
import sys
import scipy
from collections import Counter
from multiprocessing import Process, Pipe
from itertools import izip
from scipy import optimize
from scipy.optimize import minimize
from numpy import *
from numpy import linalg, zeros, array, argsort, bincount
from collections import Counter
from os import listdir,makedirs
from os.path import isfile, join,exists,split

from functools import reduce

def meshgrid(x,y):
	x = asarray(x)
	y = asarray(y)
	numRows, numCols = asarray(y.shape), asarray(x.shape)  # yes, reversed
#	print x,"\n",y 
#	print 'shapes: ', numRows,"\n",numCols
	x = x.reshape(1,numCols[0],prod(numCols[1::]))
	x = x.repeat(numRows[0], axis=0)

	y = y.reshape(numRows[0],1,prod(numRows[1::]))
	y = y.repeat(numCols[0], axis=1)
	return x, y

def meshgrid_list(x,y):
	x = asarray(x)
	y = asarray(y)

#	y = y.repeat(len(x), axis=0)
	numRows, numCols = asarray(y.shape), asarray(x.shape)    # yes, reversed
	
	x = x.reshape(numCols[0],1,prod(numCols[1::]))
	X = x.repeat(numRows[0], axis=1)
	y = y.reshape(numRows[0],prod(numRows[1::]))
	Y = asarray([y]).repeat(numCols[0], axis= 0)
	return X, Y


def load(path):
	f = open(path,'rb')
	a = pickle.load(f)
	f.close()   
	return a

def spawn(f):
    def fun(pipe,x):
        pipe.send(f(x))
        pipe.close()
    return fun

def parmap(f,X):
    pipe=[Pipe() for x in X]
    proc=[Process(target=spawn(f),args=(c,x)) for x,(p,c) in izip(X,pipe)]
    [p.start() for p in proc]
    [p.join() for p in proc]
    return [p.recv() for (p,c) in pipe]

