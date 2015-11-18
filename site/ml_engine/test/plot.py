import pandas as pd
from bokeh.plotting import figure
from bokeh.io import output_notebook, show

AAPL = pd.read_csv(
    "http://ichart.yahoo.com/table.csv?s=AAPL&a=0&b=1&c=2000&d=0&e=1&f=2015",
    parse_dates=['Date'])
MSFT = pd.read_csv(
    "http://ichart.yahoo.com/table.csv?s=MSFT&a=0&b=1&c=2000&d=0&e=1&f=2015",
    parse_dates=['Date'])
IBM = pd.read_csv(
    "http://ichart.yahoo.com/table.csv?s=IBM&a=0&b=1&c=2000&d=0&e=1&f=2015",
    parse_dates=['Date'])

def make_figure(i):
    p = figure(x_axis_type="datetime", width=500, height=400)

    if i==0:
      p.line(AAPL['Date'], AAPL['Adj Close'], color='#A6CEE3', legend='AAPL')
    elif i==1:
      p.line(IBM['Date'], IBM['Adj Close'], color='#33A02C', legend='IBM')
    elif i==2:
      p.line(MSFT['Date'], MSFT['Adj Close'], color='#FB9A99', legend='MSFT')

    p.title = "Stock Closing Prices"
    p.grid.grid_line_alpha=0.3
    p.xaxis.axis_label = 'Date'
    p.yaxis.axis_label = 'Price'
    p.legend.orientation = "top_left"
    return p

import jinja2
from bokeh.embed import components, autoload_static
from bokeh.resources import CDN

#template = jinja2.Template("""
#<!DOCTYPE html>
#<html lang="en-US">
#
#<link
#    href="http://cdn.pydata.org/bokeh/release/bokeh-0.9.0.min.css"
#    rel="stylesheet" type="text/css"
#>
#<script 
#    src="http://cdn.pydata.org/bokeh/release/bokeh-0.9.0.min.js"
#></script>
#
#<body>
#
#    <h1>Hello Bokeh!</h1>
#    
#    <p> Below is a simple plot of stock closing prices </p>
#    
#    {{ script }}
#    
#    {{ div }}
#
#</body>
#
#</html>
#""")


p = make_figure(0)
#figJS,figDiv = components(p)
figJS, figDiv = autoload_static(p, CDN, '/static/js/plot.js')
with open('static/js/plot.js', 'w') as f:
  f.write(figJS)
#f.close()

#import random
#def plotData():
#  j = random.randint(0,1000)
#  print j%3
#  p = make_figure(j%3)
#  figJS, figDiv = autoload_static(p, CDN, '/static/js/plot.js')
#  with open('static/js/plot.js', 'w') as f:
#    f.write(figJS)
#  return figJS, figDiv
##  f.close()
