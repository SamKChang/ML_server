import itertools
import re
import numpy as np
from element_list import element_dict as element_dict

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

def constraint(r_list, x_list):
  _list = list(element_dict.iterkeys())
  e_list = []
  args = [r_list, x_list]

