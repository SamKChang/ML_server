#!/usr/bin/python

from flask import Flask
from flask import request
from flask import render_template
app = Flask(__name__)

import ml_engine.converter as converter
import ml_engine.interface as interface
import ml_engine.dynamic_plot as pt
from ml_engine.ml_core.List2Energies import loadMachine

m = loadMachine(True)

@app.route('/')
def home():
  return render_template("home.html", mode=0)

@app.route('/', methods=['POST'])
def inpProcessor():
  if 'crystal' in request.form:
    mode = 0
    # input processing
    cyl = converter.stringConverter(request.form['crystal'])
    # output processing
    out = interface.outputString([cyl], m)
  elif 'target' not in request.form: 
    mode = 1
    # input processing
    atoms = []
    strList = ['a1', 'a2', 'a3', 'a4']
    for i in range(4):
      flag = strList[i]
      elist = list(set(
        converter.stringConverter(request.form[flag])))
      if len(elist) > 0:
        atoms.append(
          list(set(converter.stringConverter(request.form[flag])))
        )
        cylList = converter.groupConstructor(atoms)
        # output processing
        out = interface.outputString(cylList, m)
      else:
        out = [['\nEvery atom site must be specified.']]
  else:
    mode = 2
    # input processing
    atoms = []
    strList = ['a1', 'a2', 'a3', 'a4']
    for i in range(4):
      flag = strList[i]
      atoms.append(
        list(set(converter.stringConverter(request.form[flag])))
      )
    cylList = atoms
    target = request.form['target']
    popsize = request.form['ml_pop']
    step = request.form['ml_step']

    out = interface.optimizer(target, popsize, cylList, step, m)
    out = [['\nOptimizing crystals with\n' +\
            'target value: %.2f\n' % float(target) +\
            'total generations: %d\n' % int(step) +\
            'population size: %d\n' % int(popsize)]]


  return render_template("home.html", result=out, mode=mode)

from ml_engine.dynamic_plot import entries as entries
import json
from flask import Response
@app.route('/data', methods=['GET', 'OPTIONS', 'POST'])
@pt.crossdomain(origin="*", methods=['GET', 'POST'])
def hello_world():
    #entries
    #from ml_engine.dynamic_plot import entries as entries
    try:
        modified_since = float(
          request.headers.get('If-Modified-Since'))
    except TypeError:
        modified_since = 0

    new_entries = [e for e in entries\
                   if e.creation > modified_since]
    js = json.dumps({'x':[e.x for e in new_entries], 
                     'y':[e.y for e in new_entries],
                     'cylStr':[tuple(e.cylStr) for e in new_entries]})
    resp = Response(js, status=200, mimetype='application/json')
    print "resp:",
    print resp

    if new_entries:
        resp.headers['Last-Modified'] = new_entries[-1].creation
    elif modified_since:
        resp.headers['Last-Modified'] = modified_since

    return resp

@app.route('/done')
def done(out_cylStrLIst):
  print " yoyo from done"












@app.route('/date')
def date():
    return render_template('date.html')


@app.route('/add')
def test():
    return render_template('add.html')

@app.route('/ml_front')
def ml():
    return render_template('ml_server.html')

@app.route('/material_project')
def material():
    return render_template('materialsProject.html')

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=8080)
