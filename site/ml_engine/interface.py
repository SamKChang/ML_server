import itertools
import re
import numpy as np
from element_list import element_dict as element_dict
from ml_core.List2Energies import list2Energy
import operator
import dynamic_plot as pt

def outputString(cylList, m):
  """
  construct output string
  calling ML engine from here
  """
  supported = True
  cyl_count = 0
  atom_count = 0

  # crystal string
  cylStrList = []
  for cyl in cylList:
    cyl_count = cyl_count + 1
    atom_count = 0
    for a in range(len(cyl)):
      atom_count = atom_count + 1
      atom = cyl[a]
      if atom in element_dict:
        if a == 0:
          cylStr = atom
        else:
          cylStr = cylStr + '-' + atom
      else:
        supported = False
        cylStr = '\nElement: %s is not supported yet..' % atom
        break
    if supported and atom_count == 4:
      cylStrList.append(cylStr)
    elif atom_count == 4:
      cylStrList = cylStr
    else:
      cylStrList = '\n4 elements must be specified..'

  # calculate output
  if supported and atom_count == 4 :
    E = list2Energy(cylList, m)
    #E = List2Energies([atomList])
    #E = element_dict[atomList[0]][0]
  else:
    E = 0

  if supported and atom_count == 4: 
    result = [[s[0], s[1]] for s in zip(cylStrList, E)]
    return sorted(result, key=operator.itemgetter(1))
  else:
    return [[cylStrList]]

def optimizer(target, popsize, cylList, step, m):
  target = float(target)
  popsize = int(popsize)
  #cylList = [[str(s)] for t in s for s in cylList]
  cylList = [[str(s) for s in S] for S in cylList]
  for i in range(len(cylList)):
    if len(cylList[i]) == 0:
      cylList[i] = list(element_dict.iterkeys())

  step = int(step)
  out = pt.plotData(target, popsize, cylList, step, m)
  return out
