#!/usr/bin/python

from flask import Flask
from flask import request
from flask import render_template
import ml_engine.converter as converter
import ml_engine.interface as interface
#import ml_engine.bokeh_plot as pt
import ml_engine.plot as pt

app = Flask(__name__)

@app.route('/')
def home():
  return render_template("test.html", 
    figJS=pt.figJS, figDiv=pt.figDiv)

if __name__ == '__main__':
    app.run(debug=True)
