import itertools
import re
import numpy as np

element_dict = {
 'H': [1,1],
 'He': [2,2],
 'Li': [3,3],
 'Be': [4,4],
 'B': [5,5],
 'C': [6,5],
 'N': [7,5],
 'O': [8,5],
 'F': [9,5],
 'Na': [10,5],
}

def outputString(atomList):
  """
  construct output string
  calling ML engine from here
  """

  supported = True

  # crystal string
  for a in range(len(atomList)):
    atom = atomList[a]
    if atom in element_dict:
      if a == 0:
        cylStr = atom
      else:
        cylStr = cylStr + '-' + atom
    else:
      supported = False
      cylStr = 'Element: %s is not supported yet' % atom
      break

  # calculate output
  if supported:
    E = element_dict[atomList[0]][0]
  else:
    E = None

  return [cylStr, E]

def stringConverter(strInp):
  """
  text processing to convert single input field to unified format
  """
  spliter = re.compile(r'(, |-| |,)')
  text = re.sub(spliter, ' ', strInp)
  text = re.split(' ', text)
  text = filter(None, text)
  for t in range(len(text)):
    text[t] = text[t].title()
  #text = '-'.join(text)
  return text

def groupConstructor(atom_list):
  return [list(s)	 for s in list(itertools.product(*atom_list))]

import threading
def constraint(r_list, x_list):
  _list = list(element_dict.iterkeys())
  e_list = [s for s in _list and x not in x_list]
  args = [r_list, x_list]

  threading.Timer(0.5, constraint, *args).start()
  

