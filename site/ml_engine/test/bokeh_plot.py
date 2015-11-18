import random

# imports for Bokeh plotting
from bokeh.plotting import figure
from bokeh.resources import CDN
from bokeh.embed import file_html, components, autoload_static

# imports for matplotlib plotting
import tempfile
import matplotlib

matplotlib.use('Agg') # this allows PNG plotting

import matplotlib.pyplot as plt

# generate some random integers, sorted
exponent = .7+random.random()*.6
dta = []
for i in range(50):
  rnum = int((random.random()*10)**exponent)
  dta.append(rnum)
y = sorted(dta)
x = range(len(y))

# generate Bokeh HTML elements
# create a `figure` object
p = figure(title='crystal optimization',
    plot_width=500,plot_height=400)
# add the line
p.line(x,y)
# add axis labels
p.xaxis.axis_label = "step"
p.yaxis.axis_label = "penalty"
# create the HTML elements to pass to template
figJS,figDiv = components(p,CDN)
#figJS,figDiv = autoload_static(p,CDN,"static/temp/bokeh.js")
  
# generate matplotlib plot
fig = plt.figure(figsize=(5,4),dpi=100)
axes = fig.add_subplot(1,1,1)
# plot the data
axes.plot(x,y,'-')
# labels
axes.set_xlabel('step')
axes.set_ylabel('energy')
axes.set_title("crystal optimization")
# make the temporary file
f = tempfile.NamedTemporaryFile(
    dir='static/temp',
    suffix='.png',delete=False)
# save the figure to the temporary file
plt.savefig(f)
f.close() # close the file
# get the file's name (rather than the whole path)
# (the template will need that)
plotPng = f.name.split('/')[-1]
