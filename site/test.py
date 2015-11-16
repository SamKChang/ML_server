#!/usr/bin/python

import random
import StringIO

from flask import Flask, make_response
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure


fig = Figure()
axis = fig.add_subplot(1, 1, 1)

xs = range(100)
ys = [random.randint(1, 50) for x in xs]

axis.plot(xs, ys)
canvas = FigureCanvas(fig)
output = StringIO.StringIO()
canvas.print_png(output)
#response = make_response(output.getvalue())
#response.mimetype = 'image/png'


