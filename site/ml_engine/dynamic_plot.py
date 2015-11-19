import numpy as np
import os
from datetime import timedelta
from functools import update_wrapper, wraps
from six import string_types
from bokeh.plotting import figure, show, output_file
from bokeh.models.sources import AjaxDataSource
from bokeh.models import HoverTool
import time
from threading import Thread
from collections import namedtuple, deque
import json
from ml_core.FindTargetE import prepEve, fitFun
from numpy import argsort, asarray
from flask import Flask, jsonify, make_response, request,\
current_app, Response

Entry = namedtuple('Entry', ['x', 'y', 'creation', 'cylStr'])
max_length = 500
entries = deque(maxlen=max_length)

def gen_entry(target, NOpt, SSpace, step, m):
  NOpt = int(NOpt)
  d, popSet = prepEve(target,NOpt,SSpace, m)
  for i in range(int(step)):
    d.nextPop()
    for ip in d.pop:
      popSet.add(tuple(ip))
          
    E = [fitFun(a,target,m) for a in popSet]
    optIdxs = argsort(E)[0:NOpt]
    optInd = asarray(list(popSet))[optIdxs]
    popSet = set(tuple(ip) for ip in optInd)
    optE = [fitFun(ip,None,m) for ip in optInd]
    minE = optE[0]
    print minE,
    print optInd[0]
  
    last_entry = Entry(i, minE, time.time(), optInd[0])
    entries.append(last_entry)
  os.remove('templates/ajax_plot.html')



def plotData(target, popsize, cylList, step, m):
  for s in range(max_length):
    last_entry = Entry(None, None, None, None)
    entries.append(last_entry)
  source = AjaxDataSource(data_url='http://localhost:8080/data', 
                          mode="append",
                          if_modified=True, 
                          polling_interval=1000, 
                          max_size=max_length)
  output_file('templates/ajax_plot.html')
  hover = HoverTool(
          tooltips=[
              ("Crystal", "@cylStr"),
              ("Energy", "@y"),
              ("Generation", "@x"),
          ]
      )
  p = figure(tools=[hover], x_range=(0,step))
  p.plot_width=500
  p.plot_height=400
  p.line('x', 'y', 'cylStr', source=source)
  p.xaxis.axis_label = 'Generations'
  p.yaxis.axis_label = 'Best energy [eV/atom]'
  t = Thread(target=gen_entry, 
  args=(target, popsize, cylList, step, m,))
  t.daemon = True
  t.start()
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
