#!/usr/bin/python

from flask import Flask
from flask import request
from flask import render_template
import ml_engine.converter as cv

app = Flask(__name__)

@app.route('/')
def home():
  return render_template("home.html", mode=0)

@app.route('/results', methods=['POST'])
def inpProcessor():
  if 'crystal' in request.form:
    mode = 0
    # input processing
    cyl = cv.stringConverter(request.form['crystal'])
    # output processing
    out = [cv.outputString(cyl)]
  elif 'a1' in request.form: 
    mode = 1
    atoms = []
    strList = ['a1', 'a2', 'a3', 'a4']
    for i in range(4):
      flag = strList[i]
      atoms.append(
        cv.stringConverter(request.form[flag])
      )
    # input processing
    cylList = cv.groupConstructor(atoms)
    out = []
    # output processing
    for c in cylList:
      entry = cv.outputString(c)
      if entry[1]:
        out.append(entry)
      else:
        out = [entry]
  elif 'include' in request.form:
    mode = 2
    # input processing
    include = stringConverter(request.form['include'])
    exclude = stringConverter(request.form['exclude'])
    element_range = cv.constraint(include, exclude)
    print element_range
  
  return render_template("home.html", result=out, mode=mode)








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
