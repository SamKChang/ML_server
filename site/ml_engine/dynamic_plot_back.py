import numpy as np

from datetime import timedelta
from functools import update_wrapper, wraps
from six import string_types

from bokeh.plotting import figure, show, output_file
from bokeh.models.sources import AjaxDataSource
from bokeh.models import HoverTool

#import sys
#import os.path
#sys.path.append(
#    os.path.abspath(os.path.join(os.path.dirname(__file__), 
#                                 os.path.pardir)))
#import application as app

source = AjaxDataSource(data_url='http://localhost:5000/data', 
                        mode="append",
                        if_modified=True, 
                        polling_interval=1000, 
                        max_size=125)
output_file('templates/ajax_plot.html')
hover = HoverTool(
        tooltips=[
            ("Crystal", "@cylStr"),
            ("Energy", "@y"),
        ]
    )
p = figure(tools=[hover])
p.plot_width=500
p.plot_height=400
p.line('x', 'y', 'cylStr', source=source)
p.xaxis.axis_label = 'Generations'
p.yaxis.axis_label = 'Best energy [eV/atom]'
#html = file_html(p, CDN, 'templates/ajax_plot.html')
#html_file = open('templates/ajax_plot.html', 'w')
#print >> html_file, html
#html_file.close()
#figJS, figDiv = autoload_static(p, CDN, 'static/js/plot.js')
#with open('static/js/plot.js', 'w') as f:
#  f.write(figJS)

import time
from threading import Thread
from collections import namedtuple, deque
import json
from ml_core.FindTargetE import prepEve, fitFun
from numpy import argsort, asarray


Entry = namedtuple('Entry', ['x', 'y', 'creation', 'cylStr'])

#entries = deque(maxlen=120)
entries = []

#def gen_entry(target,NOpt,SSpace):
def gen_entry(target=None, NOpt=None, SSpace=None, step=50):
#  target = -2
#  NOpt = 50
#  SSpace = [['H'],
#            ['Al','Na','K','Cs','Rb'],
#            ['Al','Na','K','Cs','Rb','Ca','Mg','Be','Sr','Ba'],
#            ['O','F','Cl','Br','I']]
#  print "yo yo yo yo yo yo "
#  print "target=", 
#  print target
#  print "NOpt=",
#  print NOpt
#  print "SSpace=",
#  print SSpace
#  print "step=",
#  print step

  NOpt = int(NOpt)
  d, popSet, m = prepEve(target,NOpt,SSpace)
  for i in range(int(step)):
    d.nextPop()
    for ip in d.pop:
      popSet.add(tuple(ip))
          
    E = [fitFun(a,target,m) for a in popSet]
    optIdxs = argsort(E)[0:NOpt]
    optInd = asarray(list(popSet))[optIdxs]
    popSet = set(tuple(ip) for ip in optInd)
    optE = [fitFun(ip,None,m) for ip in optInd]
    #print 'Gen:',i,' || Fit: ',optE, ' || '\
    # , optInd
    minE = optE[0]
    print minE,
    print optInd[0]
  
    last_entry = Entry(i, minE, time.time(), optInd[0])
    entries.append(last_entry)
    #entries = optE
    #time.sleep(1)
  #return optInd, optE
  del entries
  entries = []
  global entries


#import random
#def gen_entry():
##    global entries
#    x = 0
#    for e in entries:
#      del e
#    print entries
#    for _itr in range(10):
#        print _itr
#        rd = random.random()
#        last_entry = Entry(x, rd, time.time())
#        entries.append(last_entry)
#        x += 1
#        time.sleep(1)

from flask import Flask, jsonify, make_response, request,\
current_app, Response

def plotData(target, popsize, cylList, step):
  t = Thread(target=gen_entry, args=(target, popsize, cylList, step,))
  #t = Thread(target=gen_entry)
  t.daemon = True
  t.start()
  #print "start joining"
  #t.join()
  #print "end joining"
  
#t = Thread(target=gen_entry)
#t.daemon = True
#t.start()
show(p)

#########################################################
# Flask server related
#
# The following code has no relation to bokeh and it's only
# purpose is to serve data to the AjaxDataSource instantiated
# previously. Flask just happens to be one of the python
# web frameworks that makes it's easy and concise to do so
#########################################################

def crossdomain(origin=None, methods=None, headers=None,
                max_age=21600, attach_to_all=True,
                automatic_options=True):
    """
    Decorator to set crossdomain configuration on a Flask view
    For more details about it refer to:
    http://flask.pocoo.org/snippets/56/
    """
    if methods is not None:
        methods = ', '.join(sorted(x.upper() for x in methods))

    if headers is not None\
    and not isinstance(headers, string_types):
        headers = ', '.join(x.upper() for x in headers)

    if not isinstance(origin, string_types):
        origin = ', '.join(origin)

    if isinstance(max_age, timedelta):
        max_age = max_age.total_seconds()

    def get_methods():
        # called periodically, always return methods
        return methods
        print "  yoyo-get_methods" # not been called!
        options_resp = current_app.make_default_options_response()
        return options_resp.headers['allow']

    def decorator(f):
        print " yoyo-decorator" # only been called at begining
        @wraps(f)
        def wrapped_function(*args, **kwargs):
            if automatic_options and request.method == 'OPTIONS':
                resp = current_app.make_default_options_response()
            else:
                resp = make_response(f(*args, **kwargs))
            if not attach_to_all and request.method != 'OPTIONS':
                return resp

            h = resp.headers

            h['Access-Control-Allow-Origin'] = origin
            h['Access-Control-Allow-Methods'] = get_methods()
            h['Access-Control-Max-Age'] = str(max_age)
            requested_headers = request.headers.get(
                'Access-Control-Request-Headers'
            )
            if headers is not None:
                h['Access-Control-Allow-Headers'] = headers
            elif requested_headers :
                h['Access-Control-Allow-Headers'] =\
                  requested_headers
            return resp
        f.provide_automatic_options = False
        return update_wrapper(wrapped_function, f)

    return decorator
