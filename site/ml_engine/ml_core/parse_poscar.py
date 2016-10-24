#!/usr/bin/env python

import sys, os.path
from numpy import *

periodictable_numbers={'H':1,'He':2,'Li':3,'Be':4,'B':5,'C':6,'N':7,'O':8,'F':9,'Ne':10,'Na':11,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,'Cl':17,'Ar':18,
                       'K':19,'Ca':20,'Sc':21,'Ti':22,'V':23,'Cr':24,'Mn':25,'Fe':26,'Co':27,'Ni':28,'Cu':29,'Zn':30,'Ga':31,'Ge':32,'As':33,'Se':34,
                       'Br':35,'Kr':36,'Rb':37,'Sr':38,'Y':39,'Zr':40,'Nb':41,'Mo':42,'Tc':43,'Ru':44,'Rh':45,'Pd':46,'Ag':47,'Cd':48,'In':49,'Sn':50,
                       'Sb':51,'Te':52,'I':53,'Xe':54,'Cs':55,'Ba':56,'La':57,'Ce':58,'Pr':59,'Nd':60,'Pm':61,'Sm':62,'Eu':63,'Gd':64,'Tb':65,'Dy':66,
                       'Ho':67,'Er':68,'Tm':69,'Yb':70,'Lu':71,'Hf':72,'Ta':73,'W':74,'Re':75,'Os':76,'Ir':77,'Pt':78,'Au':79,'Hg':80,'Tl':81,'Pb':82,
                       'Bi':83,'Po':84,'At':85,'Rn':86,'Fr':87,'Ra':88,'Ac':89,'Th':90,'Pa':91,'U':92,'Np':93,'Pu':94,'Am':95,'Cm':96,'Bk':97,'Cf':98,
                       'Es':99,'Fm':100,'Md':101,'No':102,'Lr':103,'Rf':104,'Db':105,'Sg':106,'Bh':107,'Hs':108,'Mt':109,'Ds':110,'Rg':111,'Cn':112,
                       'Uut':113,'Uuq':114,'Uup':115,'Uuh':116,'Uus':117,'Uuo':118}

periodictable_symbols=asarray([0,'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co',
                       'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',
                       'I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir',
                       'Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No',
                       'Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Uut','Uuq','Uup','Uuh','Uus','Uuo'])

def poscar_parse(f):
    """
    Parses a file on VASPs POSCAR format. Returns 
      (cell, scale, vol, coords, coords_reduced, counts, occupations, comment)
    where
      cell: 3x3 nested list of *strings* designating the cell
      scale: *string* representing the overall scale of the cell
      vol: *string* representing the volume of the cell (only one of scale and vol will be set, the other one = None)
      coords: Nx3 nested list of *strings* designating the coordinates
      coords_reduced: bool, true = coords are given in reduced coordinate (in vasp D or Direct), false = coords are given in cartesian coordinates
      counts: how many atoms of each type
      occupations: which species of each atom type
      comment: the comment string given at the top of the file
    """
    fi = iter(f)

    comment = next(fi).strip()
    vol_or_scale = next(fi).strip()
    vol_or_scale_nbr = float(vol_or_scale)
    if vol_or_scale_nbr < 0:
        vol = vol_or_scale[1:]
        scale = None
    else:
        scale = vol_or_scale
        vol = None

    cell = [['','',''],['','',''],['','','']]        
    for i in [0,1,2]:
        cellline = next(fi).strip().split()
        for j,v in enumerate(cellline):
            cell[i][j] = v
   
    symbols_or_count = next(fi).strip().split()

    try:
        counts = map(int,symbols_or_count)
        symbols = None
        occupations = range(len(counts))
    except Exception:
        symbols = symbols_or_count
        counts = [int(s) for s in next(fi).strip().split()]
        occupations = [periodictable_numbers[symbol] for symbol in symbols]

    N = sum(counts)

    coordtype_or_selectivedynamics = next(fi).strip()
    if coordtype_or_selectivedynamics[0] in 'Ss':       
        # Skip row if selective dynamics specifier
        coordtype = next(fi).strip()        
    else:
        coordtype = coordtype_or_selectivedynamics

    if coordtype[0] in 'CcKk':
        coords_reduced = False
    else:
        coords_reduced = True

    coords = []
    for i in range(N):
        strcoord = next(fi).strip().split()[:3]
        coord = map(lambda x: x.strip(),strcoord)
        coords.append(coord)

    return (cell, scale, vol, coords, coords_reduced, counts, occupations, comment)

def poscar_parse_guss(f):
    """
    Parses a file on VASPs POSCAR format. Returns 
      (cell, scale, vol, coords, coords_reduced, counts, occupations, comment)
    where
      cell: 3x3 nested list of *strings* designating the cell
      scale: *string* representing the overall scale of the cell
      vol: *string* representing the volume of the cell (only one of scale and vol will be set, the other one = None)
      coords: Nx3 nested list of *strings* designating the coordinates
      coords_reduced: bool, true = coords are given in reduced coordinate (in vasp D or Direct), false = coords are given in cartesian coordinates
      counts: how many atoms of each type
      occupations: which species of each atom type
      comment: the comment string given at the top of the file

    """
    fi = iter(f)

    comment = next(fi).strip()
    vol_or_scale = next(fi).strip()
    vol_or_scale_nbr = float(vol_or_scale)
    if vol_or_scale_nbr < 0:
        vol = vol_or_scale[1:]
        scale = None
    else:
        scale = vol_or_scale
        vol = None

    cell = [['','',''],['','',''],['','','']]        
    for i in [0,1,2]:
        cellline = next(fi).strip().split()
        for j,v in enumerate(cellline):
            cell[i][j] = v
   
    symbols_or_count = next(fi).strip().split()

    try:
        counts = map(int,symbols_or_count)
        symbols = None
        occupations = []
    except Exception:
        symbols = symbols_or_count
        counts = [int(s) for s in next(fi).strip().split()]
        

    N = sum(counts)

    coordtype_or_selectivedynamics = next(fi).strip()
    if coordtype_or_selectivedynamics[0] in 'Ss':       
        # Skip row if selective dynamics specifier
        coordtype = next(fi).strip()        
    else:
        coordtype = coordtype_or_selectivedynamics

    if coordtype[0] in 'CcKk':
        coords_reduced = False
    else:
        coords_reduced = True

    coords = []
    old = None
    for i in range(N):
	tmp = next(fi).strip().split()
        strcoord = tmp[:3]
	if old != tmp[3]:
		occupations.append(periodictable_numbers[tmp[3]])
		old = tmp[3]

        coord = map(lambda x: x.strip(),strcoord)
        coords.append(coord)

    return (cell, scale, vol, coords, coords_reduced, counts, occupations, comment)


def poscar_output(f, cell, coords, coords_reduced, counts, occupations,comment="Comment",scale="1",vol=None):
    """
    Writes a file on VASPs POSCAR format. When string type are specified for input arguments below, any datatype that converts 
    cleanly to a string is also accepted, e.g., floats and numpy arrays.

    Input arguments  
      f: file stream to put output on  
      cell: 3x3 nested list of strings designating the cell
      coords: Nx3 nested list of strings designating the coordinates
      coords_reduced: bool, true = coords are given in reduced coordinate (in vasp D or Direct), false = coords are given in cartesian coordinates
      counts: how many atoms of each type
      occupations: which species of each atom type
      comment: (optional) the comment string given at the top of the file
      scale: (optional) string representing the overall scale of the cell
      vol: string representing the volume of the cell (only one of scale and vol can be set)
    """        
    f.write(str(comment)+"\n")
    if vol != None:
        f.write("-"+str(vol)+"\n")
    else:
        f.write(str(scale)+"\n")
    for c1, c2, c3 in cell:  
        f.write(str(c1)+" "+str(c2)+" "+str(c3)+"\n")

    for i in range(len(counts)):
        if occupations == None:
            f.write(periodictable_symbols[i] + " ")
        else:
            if isinstance(occupations[i], ( int, long ) ):
                f.write(periodictable_symbols[occupations[i]] + " ")
            else:
                f.write(str(occupations[i]) + " ")
    f.write("\n")
 
    for count in counts:
        f.write(str(count) + " ")
    f.write("\n")
    if coords_reduced:
        f.write("D\n")
    else:
        f.write("K\n")
    for c1, c2, c3 in coords:
        f.write(str(c1)+" "+str(c2)+" "+str(c3)+"\n")
    
def vol_to_scale(cell,vol):
    C = [[0,0,0],[0,0,0],[0,0,0]]        
    for i in [0,1,2]:
        for j in [0,1,2]:
            C[i][j] = float(cell[i][j])    
    # Done manually to avoid any dependency on numpy or similar. This works with *any* array-type datatype.
    celldet = C[0][0]*C[1][1]*C[2][2] + C[0][1]*C[1][2]*C[2][0] + C[0][2]*C[1][0]*C[2][1] - C[0][2]*C[1][1]*C[2][0] - C[0][1]*C[1][0]*C[2][2] - C[0][0]*C[1][2]*C[2][1]
    return "%.14f" % (float(vol)/abs(celldet))**(1.0/3.0) 

def cartesian_to_reduced(cell,scale,coords):
    C = [[0,0,0],[0,0,0],[0,0,0]]        
    for i in [0,1,2]:
        for j in [0,1,2]:
            C[i][j] = float(cell[i][j])        
        
    m = 1.0
    # Done manually to avoid any dependency on numpy or similar. This works with *any* array-type datatype.
    cellinv = [ [ m*(C[1][1]*C[2][2] - C[1][2]*C[2][1]), m*(C[0][2]*C[2][1]-C[0][1]*C[2][2]), m*(C[0][1]*C[1][2]-C[0][2]*C[1][1]) ],
                   [ m*(C[1][2]*C[2][0] - C[1][0]*C[2][2]), m*(C[0][0]*C[2][2]-C[0][2]*C[2][0]), m*(C[0][2]*C[1][0]-C[0][0]*C[1][2]) ],
                   [ m*(C[1][0]*C[2][1] - C[1][1]*C[2][0]), m*(C[0][1]*C[2][0]-C[0][0]*C[2][1]), m*(C[0][0]*C[1][1]-C[0][1]*C[1][0]) ] ]

    scale = float(scale)

    newcoords = []        

    for strcoord in coords:
        coord = map(float,strcoord)
        # Done manually to avoid any dependency on numpy or similar. This works with *any* array-type datatype.
        reduced_coord = [ scale * (coord[0] * cellinv[0][0] + coord[1] * cellinv[0][1] + coord[2] * cellinv[0][2]) , 
                        scale * (coord[0] * cellinv[1][0] + coord[1] * cellinv[1][1] + coord[2] * cellinv[1][2]) ,
                        scale * (coord[0] * cellinv[2][0] + coord[1] * cellinv[2][1] + coord[2] * cellinv[2][2]) ]
        newcoords.append(reduced_coord)

    return newcoords


def main():
    #print "TEST"
    #parser = argparse.ArgumentParser(description="Simple example to parse a poscar.")
    #parser.add_argument('file', metavar='run', help='a filename to a POSCAR file with the input crystal structure')
    #args = parser.parse_args()    
    filename = sys.argv[1]

    if not os.path.exists(filename):
        print "File not found."
        sys.exit(1)        

    f = open(filename)
    cell, scale, vol, coords, coords_reduced, counts, occupations, comment = poscar_parse(f)

    if vol != None and scale == None:
        scale = vol_to_scale(cell,vol)
        
    if coords_reduced == True:
        coords = cartesian_to_reduced(cell, scale, coords)

    # Convert to numpy objects
    cell = array(cell)
    scale = float(scale)
    coords = array(coords)

    # Just print out everything on POSCAR form again
    poscar_output(sys.stdout, cell, coords, True, counts, occupations, scale=scale, comment = comment)  

if __name__ == "__main__":
    main()

