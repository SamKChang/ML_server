#!/usr/bin/python

from flask import Flask
from flask import request
from flask import render_template
import ml_engine.converter as converter
import ml_engine.interface as interface
import ml_engine.dynamic_plot as pt

app = Flask(__name__)

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
    out = interface.outputString([cyl])
  elif 'target' not in request.form: 
    mode = 1
    # input processing
    atoms = []
    strList = ['a1', 'a2', 'a3', 'a4']
    for i in range(4):
      flag = strList[i]
      atoms.append(
        list(set(converter.stringConverter(request.form[flag])))
      )
    cylList = converter.groupConstructor(atoms)
    # output processing
    out = interface.outputString(cylList)
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

    out = interface.optimizer(target, popsize, cylList, step)
    out = [['Optimizing crystals with target value: %f'\
          % float(target)]]

    #pt.plotData()

  return render_template("home.html", result=out, mode=mode)

from ml_engine.dynamic_plot import entries as entries
import json
from flask import Response
@app.route('/data', methods=['GET', 'OPTIONS', 'POST'])
@pt.crossdomain(origin="*", methods=['GET', 'POST'])
def hello_world():
    #global entries
    entries
    try:
        modified_since = float(
          request.headers.get('If-Modified-Since'))
    except TypeError:
        modified_since = 0

    new_entries = [e for e in entries\
                   if e.creation > modified_since]
    js = json.dumps({'x':[e.x for e in new_entries], 
                     'y':[e.y for e in new_entries]})
    resp = Response(js, status=200, mimetype='application/json')
    print "resp:",
    print resp

    if new_entries:
        resp.headers['Last-Modified'] = new_entries[-1].creation
    elif modified_since:
        resp.headers['Last-Modified'] = modified_since

    return resp














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
    app.run(debug=True)
