#!/usr/bin/python
import re
import MachineLearning as ml
import cPickle as pkl
from Kernels import *
from numpy import *
from copy import deepcopy

PN={'H':1,'He':2,'Li':3,'Be':4,'B':5,'C':6,'N':7,'O':8,'F':9,'Ne':10,'Na':11,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,'Cl':17,'Ar':18,
							   'K':19,'Ca':20,'Sc':21,'Ti':22,'V':23,'Cr':24,'Mn':25,'Fe':26,'Co':27,'Ni':28,'Cu':29,'Zn':30,'Ga':31,'Ge':32,'As':33,'Se':34,
							   'Br':35,'Kr':36,'Rb':37,'Sr':38,'Y':39,'Zr':40,'Nb':41,'Mo':42,'Tc':43,'Ru':44,'Rh':45,'Pd':46,'Ag':47,'Cd':48,'In':49,'Sn':50,
							   'Sb':51,'Te':52,'I':53,'Xe':54,'Cs':55,'Ba':56,'La':57,'Ce':58,'Pr':59,'Nd':60,'Pm':61,'Sm':62,'Eu':63,'Gd':64,'Tb':65,'Dy':66,
							   'Ho':67,'Er':68,'Tm':69,'Yb':70,'Lu':71,'Hf':72,'Ta':73,'W':74,'Re':75,'Os':76,'Ir':77,'Pt':78,'Au':79,'Hg':80,'Tl':81,'Pb':82,
							   'Bi':83,'Po':84,'At':85,'Rn':86,'Fr':87,'Ra':88,'Ac':89,'Th':90,'Pa':91,'U':92,'Np':93,'Pu':94,'Am':95,'Cm':96,'Bk':97,'Cf':98,
							   'Es':99,'Fm':100,'Md':101,'No':102,'Lr':103,'Rf':104,'Db':105,'Sg':106,'Bh':107,'Hs':108,'Mt':109,'Ds':110,'Rg':111,'Cn':112,
							   'Uut':113,'Uuq':114,'Uup':115,'Uuh':116,'Uus':117,'Uuo':118}
	
PS=asarray([0,'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co',
                       'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',
                       'I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir',
                       'Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No',
                       'Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Uut','Uuq','Uup','Uuh','Uus','Uuo'])	
	
PTP = {\
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

PureStateE = {1:-.13582598E+02/4,2:-.23398159E-01/2,3:-.57148112E+01/3,4:-.74850320E+01/2,5:-.80162036E+02/12,6:-.36903874E+02/4,\
                7:-.66723762E+02/8,8:-.99042808E+01/2,9:-.76342895E+01/4,10:-.27979331E-01/1,11:-.39315371E+01/3,\
                12:-.32041556E+01/2,13:-.37605221E+01/1,14:-.10786680E+02/2,15:-.21512281E+02/4,16:-.13237729E+03/32,17:-.73953526E+01/4,\
                18:-.69673050E-01/1,19:-.10268679E+01/1,20:-.19806149E+01/1,31:-.12109324E+02/4,32:-.92361456E+01/2,\
                33:-.93147591E+01/2,34:-.10466293E+02/3,35:-.65344764E+01/4,36:-.53986466E-01/1,37:-.93246713E+00/1,\
                38:-.16870687E+01/1,49:-.27220994E+01/1,50:-.80051039E+01/2,51:-.82482387E+01/2,52:-.94307332E+01/3,\
                53:-.60693428E+01/4,54:-.35831113E-01/1,55:-.85401477E+00/1,56:-.19232545E+01,81:-.47273182E+01/2,\
                82:-.37408443E+01/1,83:-.24231858E+02/6   }

reversePTP = zeros(8*34).reshape(8,34)
for k, v in PTP.items():
	reversePTP[tuple(v)] = k 
	

	
	
def loadMachine(trained = False):
	def prepT(T,Z,N):		
		for i in range(len(T)):
			T[i] = ( T[i] - sum([PureStateE[z] for z in Z[i]]) )/N[i]
		return T
	path  =  'static/ml_data/trainingSet/NonSpacial-mp_AWG:225:1aA1bB1cC1eD_39Elements_Symmetrized.pkl'
	if not trained:
		m = ml.ML(path,sigma = 8,oFPunish = 1E-8,linOfPunish = 20, randomize = 1,ker = LapKer(),num_cores = 3)
		m.o['T'] = prepT(m.o['T'],m.o['Z'],m.o['N'])
		m.FindAlpha()
		print 'Machine Trained'
		print path.rstrip('pkl') + 'ml'
		f = open('static/ml_data/trainingSet/alphas', 'w' )
		pkl.dump([deepcopy(m.o),deepcopy(m.atomRegConst)],f)
		f.close()
	else:
		f = open('static/ml_data/trainingSet/alphas', 'r')
		tmp = pkl.load(f)
		f.close()
		m = ml.ML(tmp[0],sigma = 8,oFPunish = 1E-8,linOfPunish = 20, randomize = 1,ker = LapKer(),num_cores = 3)
		m.atomRegConst = tmp[1]
	return m


def O2Desc(O):
	C = asarray([PN[c] for c in O]).astype(int)
	X = asarray([PTP[c] for c in C])
	OL = asarray([C[0],C[1],C[2],C[2],C[3],C[3],C[3],C[3],C[3],C[3]]) # HARDCODED! Needs to be changed.
	tmp = bincount(OL)
	a = zeros(118)
	a[0:0+len(tmp)] = tmp            
	a[0] = 1	   
	#print a
	#print C
	#print X
	return X,a

def list2Desc(inpList):

	X = []
	A = []
	for o in inpList:
		x,a = O2Desc(o)
		X.append(x)
		A.append(a)
		
	return asarray(X),asarray(A)

def list2Energy(inpList):
	m 	= 	loadMachine(True)
	X,A =	list2Desc(inpList)
	E = m.Energy_list(X,A)
	return E
	
if __name__ == '__main__':
	path = '/home/felix/Dropbox/PHD/serverberk/trainingSet/NonSpacial-mp_AWG:225:1aA1bB1cC1eD_39Elements_Symmetrized.pkl'
	E = list2Energy([['H','Li','Ba','F'],['H','Xe','N','C'],['K','Al','Ca','F']])
	print E 
	
	
	
