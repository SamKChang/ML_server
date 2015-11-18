from bokeh.resources import CDN
from bokeh.plotting import figure
from bokeh.embed import autoload_static

plot = figure()
plot.circle([1,2], [3,4])

js, tag = autoload_static(plot, CDN, "some/path")

print js
print tag
